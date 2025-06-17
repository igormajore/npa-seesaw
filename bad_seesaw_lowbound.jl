using LinearAlgebra, Ket, JuMP, SCS 

#=

Ceci est une version non optimisée de seesaw_lowbound.jl. Dans la version initiale, je voulais utiliser des paramètres. Le problème est que si SCS (avec ParametricOptInterface) supporte les expressions de type Paramètre * Variable 
(ie ne les considère pas comme quadratiques), il ne supporte pas Paramètre * Paramètre * Variable (il considère que c'est une expression non linéaire). Je n'ai donc pas réussi à utiliser les paramètres. Dans cette version, je crée donc 
un nouveau modèle pour chaque problème que je veux résoudre. 




Un jeu G est un uplet G=(n,T,W,v0,v1,Pinit) où n est le nombre de joueurs, T est l'ensemble des types de départ autorisés, W est la liste des situations gagnantes (c'est une liste d'uplets (a,t)), v0 et v1 sont les valeurs
de gain, Pinit est la distribution de probabilité des types des joueurs (c'est une fonction de T -> [0,1]).

La variable InitialValues est un tableau de taille n+1 tel que : 
    - pour 1<=i<=n, InitialValues[i] est une matrice 2*2 telle que InitialValues[i][a_i,t_i] soit la mesure du joueur i, pour l'action a_i sachant qu'il a le type t_i. 
    - InitialValues[n+1] contient rho

Comme les actions et types sont à valeurs dans {0,1} mais les indices d'array sont dans {1,...}, il y a un décalage d'indices

=#






# Fonctions auxiliaires

function krons(lst) # calcule le produit tensoriel de toutes les matrices de lst, supposée non vide
    n = length(lst)
    @assert n>0
    prod = lst[1]
    for i in 2:n 
        prod = kron(prod,lst[i])
    end
    return prod 
end

povm_canonique = [ [[1,0] [0,0]] , [[0,0] [0,1]] ]
povm_canonique_renverse = [ [[0,0] [0,1]] , [[1,0] [0,0]] ]
povm_hadamard = [ (1/2)*[[1,1] [1,1]] , (1/2)*[[1,-1] [-1,1]] ]
povm_hadamard_renverse = [ (1/2)*[[1,-1] [-1,1]] , (1/2)*[[1,1] [1,1]] ]
povms1 = [povm_canonique,povm_hadamard]
povms2 = [povm_canonique,povm_canonique_renverse,povm_hadamard,povm_hadamard_renverse]
















# Fonctions de gain

function u(i,a,t,G) # Fonction de gain du joueur i 
    _,_,W,v0,v1,_ = G 
    if (a,t) in W 
        return (if a[i]==0 v0 else v1 end)
    else 
        return 0 
    end
end

function P(a,t,param,G) # corrélation P 
    n,_,_,_,_,_ = G
    to_kron(a,t) = [param[i][a[i]+1,t[i]+1] for i in 1:n] # liste des matrices dont on va calculer le produit tensoriel
    return LinearAlgebra.tr(param[n+1] * krons(to_kron(a,t)))
end

function SW(param,G)
    n,_,W,_,_,Pinit = G 


    return (1/n) * sum(u(i,a,t,G) * Pinit(t) * P(a,t,param,G) for (a,t) in W for i in 1:n)
end

function Gain_for_i(i,param,G)
    _,_,W,_,_,Pinit = G 

    return sum(u(i,a,t,G) * Pinit(t) * P(a,t,param,G) for (a,t) in W)
end













# Itérations du see-saw

function one_iteration_SW(InitialValues,G,k) # Applique une itération du see-saw, et modifie InitialValues avec les nouvelles valeurs. 
    n,_,_,_,_,_ = G

    # SDP sur rho 
    model = Model(SCS.Optimizer)
    set_silent(model)
    @variable(model,rho[1:(k^n),1:(k^n)],PSD) 
    @constraint(model,LinearAlgebra.tr(rho)==1)
    param = [if j==n+1 rho else InitialValues[j] end for j in 1:(n+1)]
    @objective(model,Max,SW(param,G))
    JuMP.optimize!(model)
    InitialValues[n+1] = JuMP.value(rho)

    # SDP sur les mesures de chaque joueur, dans l'ordre croissant 
    for i in 1:n
        model = Model(SCS.Optimizer)
        set_silent(model)
        @variable(model,N00[1:k,1:k],PSD)
        @variable(model,N10[1:k,1:k],PSD)
        @variable(model,N01[1:k,1:k],PSD)
        @variable(model,N11[1:k,1:k],PSD)
        @constraint(model,N00+N10==LinearAlgebra.I)
        @constraint(model,N01+N11==LinearAlgebra.I)
        N = [[N00,N10] [N01,N11]]
        param = [if j==i N else InitialValues[j] end for j in 1:(n+1)]
        @objective(model,Max,SW(param,G))
        JuMP.optimize!(model)
        InitialValues[i] = [[JuMP.value(N00),JuMP.value(N10)] [JuMP.value(N01),JuMP.value(N11)]]
    end
end

function one_iteration_Gain_for_i(InitialValues,G,k) # Applique une itération du see-saw où chaque joueur optimise son propre gain (pour arriver à un équilibre)
    n,_,_,_,_,_ = G 

    # SDP sur les mesures de chaque joueur, dans l'ordre croissant 
    for i in 1:n 
        model=Model(SCS.Optimizer)
        set_silent(model)
        @variable(model,N00[1:k,1:k],PSD)
        @variable(model,N10[1:k,1:k],PSD)
        @variable(model,N01[1:k,1:k],PSD)
        @variable(model,N11[1:k,1:k],PSD)
        @constraint(model,N00+N10==LinearAlgebra.I)
        @constraint(model,N01+N11==LinearAlgebra.I)
        N = [[N00,N10] [N01,N11]]
        param = [if j==i N else InitialValues[j] end for j in 1:(n+1)]
        @objective(model,Max,Gain_for_i(i,param,G))
        JuMP.optimize!(model)
        InitialValues[i] = [[JuMP.value(N00),JuMP.value(N10)] [JuMP.value(N01),JuMP.value(N11)]]   
    end
end 

function iterations_SW(InitialValues,G,k,eps) # fait des itérations jusqu'à obtenir un SW localement optimal, avec comme critère de convergence eps
    prevInitialValues = copy(InitialValues)

    while true 
        one_iteration_SW(InitialValues,G,k)
        sum(norm.(prevInitialValues-InitialValues)) > eps || break 
        prevInitialValues = copy(InitialValues)
    end
end

function iterations_Gain_for_i(InitialValues,G,k,eps) # fait des itérations jusqu'à obtenir un équilibre, avec comme critère de convergence eps
    prevInitialValues = copy(InitialValues)

    while true 
        one_iteration_Gain_for_i(InitialValues,G,k)
        # print("anotherone : delta = ",sum(norm.(prevInitialValues-InitialValues)),"\n")
        sum(norm.(prevInitialValues-InitialValues)) > eps || break 
        prevInitialValues = copy(InitialValues)
    end
end















# Valeurs initiales et tests

function given_InitialValues(G,k,povms) # povms est une liste de POVMs. Cette fonction renvoie toutes les InitialValues qu'on peut construire en partant des mesures dans la liste povms (et en prenant des états partagés aléatoires)
    n,_,_,_,_,_ = G
    
    un_joueur = [[m1 m2] for m1 in povms for m2 in povms] # toutes les combinaisons de mesures possibles pour un joueur : une mesure pour chaque type 

    function combinaisons(p) # renvoie la liste de toutes les InitialValues à tester (sans l'état partagé, qui sera choisi aléatoirement ensuite) s'il y a p joueurs (pour récursion)
        if p==0 
            return [[]]
        else
            res = []
            for InitialValues in combinaisons(p-1)
                for m in un_joueur 
                    newInitialValues = copy(InitialValues)
                    push!(newInitialValues,m)
                    push!(res,newInitialValues)
                end
            end
            return res 
        end
    end
    
    AllInitialValues = combinaisons(n)

    # Ajout des états partagés aléatoires 
    for InitialValues in AllInitialValues 
        push!(InitialValues,random_state(Float64,k^n))
    end 

    return AllInitialValues
end

function random_InitialValues(G,k)
    n,_,_,_,_,_ = G 

    InitialValues = []
    for i in 1:n 
        M0 = random_povm(Float64,k,2) # on prend des POVMs réels car les SDP ont pour variables des matrices symétriques positives
        M1 = random_povm(Float64,k,2)
        M = [M0 M1]
        push!(InitialValues,M)
    end
    rho = random_state(Float64,k^n)
    push!(InitialValues,rho) 
end

function many_tests(G,k,povms,eps,nb_tests,barre) # Fait plusieurs tests en partant de POVMs aléatoires et de POVMs donnés par la liste povms. Affiche ceux dont le social welfare dépasse barre
    # POVMs aléatoires
    print("POVMs aléatoires :\n")
    for i in 1:nb_tests 
        InitialValues = random_InitialValues(G,k)
        iterations_SW(InitialValues,G,k,eps)
        iterations_Gain_for_i(InitialValues,G,k,eps)
        sw = SW(InitialValues,G)
        if sw>barre
            print(SW(InitialValues,G))
            print("\n")
        end
    end

    # POVMs donnés 
    print("POVMs contrôlés :\n")
    for InitialValues in given_InitialValues(G,k,povms)
        iterations_SW(InitialValues,G,k,eps)
        iterations_Gain_for_i(InitialValues,G,k,eps)
        sw = SW(InitialValues,G)
        if sw>barre
            print(SW(InitialValues,G))
            print("\n")
        end
    end
end