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

function enum_tuples(n) # Énumère {0,1}^n sous forme d'une liste de n-uplets
    if n==0
       return [()]
    else 
       return vcat([(0,a...) for a in enum_tuples(n-1)],[(1,a...) for a in enum_tuples(n-1)])
    end
end

function change_tuple(x,i,yi,n) # remplace la i-ème coordonnée de l'uplet x (de taille n) par yi (renvoie un nouvel uplet)
    return ntuple(j -> if j==i yi else x[j] end,n)
end

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
Z = [[1,0] [0,-1]]
















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













# see-saw pour Q

function one_iteration_SW(InitialValues,G,k,delta) # Applique une itération du see-saw, et modifie InitialValues avec les nouvelles valeurs. 
    n,_,_,_,_,_ = G

    # SDP sur rho 
    current_SW = SW(InitialValues,G)
    current_rho = InitialValues[n+1]

    model = Model(SCS.Optimizer)
    set_silent(model) 

    @variable(model,rho[1:(k^n),1:(k^n)],PSD) 
    set_start_value.(rho,InitialValues[n+1])

    @constraint(model,LinearAlgebra.tr(rho)==1)

    param = [if j==n+1 rho else InitialValues[j] end for j in 1:(n+1)]
    @objective(model,Max,SW(param,G))
    JuMP.optimize!(model)
    InitialValues[n+1] = JuMP.value(rho)

    if SW(InitialValues,G) < current_SW + delta # si le sw n'est pas assez amélioré, on reste sur la solution qu'on avait de base 
        InitialValues[n+1] = current_rho 
    end


    # SDP sur les mesures de chaque joueur, dans l'ordre croissant 
    for i in 1:n
        current_SW = SW(InitialValues,G)
        current_Mi = InitialValues[i]

        model = Model(SCS.Optimizer)
        set_silent(model)

        @variable(model,N00[1:k,1:k],PSD)
        set_start_value.(N00,InitialValues[i][1,1])

        @variable(model,N10[1:k,1:k],PSD)
        set_start_value.(N10,InitialValues[i][2,1])

        @variable(model,N01[1:k,1:k],PSD)
        set_start_value.(N01,InitialValues[i][1,2])

        @variable(model,N11[1:k,1:k],PSD)
        set_start_value.(N11,InitialValues[i][2,2])

        @constraint(model,N00+N10==LinearAlgebra.I)
        @constraint(model,N01+N11==LinearAlgebra.I)

        N = [[N00,N10] [N01,N11]]
        param = [if j==i N else InitialValues[j] end for j in 1:(n+1)]
        @objective(model,Max,SW(param,G))
        JuMP.optimize!(model)
        InitialValues[i] = [[JuMP.value(N00),JuMP.value(N10)] [JuMP.value(N01),JuMP.value(N11)]]

        if SW(InitialValues,G) < current_SW + delta # si le sw n'est pas assez amélioré, on reste sur la solution qu'on avait de base 
            InitialValues[i] = current_Mi
        end
    end
end

function one_iteration_Gain_for_i(InitialValues,G,k,delta) # Applique une itération du see-saw où chaque joueur optimise son propre gain (pour arriver à un équilibre)
    n,_,_,_,_,_ = G 

    # SDP sur les mesures de chaque joueur, dans l'ordre croissant 
    for i in 1:n 
        current_Gain_for_i = Gain_for_i(i,InitialValues,G)
        current_Mi = InitialValues[i]

        model=Model(SCS.Optimizer)
        set_silent(model)

        @variable(model,N00[1:k,1:k],PSD)
        set_start_value.(N00,InitialValues[i][1,1])

        @variable(model,N10[1:k,1:k],PSD)
        set_start_value.(N10,InitialValues[i][2,1])

        @variable(model,N01[1:k,1:k],PSD)
        set_start_value.(N01,InitialValues[i][1,2])

        @variable(model,N11[1:k,1:k],PSD)
        set_start_value.(N11,InitialValues[i][2,2])

        @constraint(model,N00+N10==LinearAlgebra.I)
        @constraint(model,N01+N11==LinearAlgebra.I)
        
        N = [[N00,N10] [N01,N11]]
        param = [if j==i N else InitialValues[j] end for j in 1:(n+1)]
        @objective(model,Max,Gain_for_i(i,param,G))
        JuMP.optimize!(model)
        InitialValues[i] = [[JuMP.value(N00),JuMP.value(N10)] [JuMP.value(N01),JuMP.value(N11)]]   

        if Gain_for_i(i,InitialValues,G) < current_Gain_for_i + delta # si le joueur i n'améliore pas assez son gain, on reste sur la solution qu'on avait de base 
            InitialValues[i] = current_Mi
        end
    end
end 

function iterations_SW(InitialValues,G,k,delta,eps,seuil) # fait des itérations jusqu'à obtenir un SW localement optimal, avec comme critère de convergence eps
    n,_,_,_,_,_=G
    
    prevInitialValues = []
    nb_iter = 1

    while nb_iter <= seuil
        prevInitialValues = copy(InitialValues)
        one_iteration_SW(InitialValues,G,k,delta)
        sum(norm.(prevInitialValues-InitialValues)) > eps || break 
        nb_iter+=1
    end

    if nb_iter > seuil  # pas de convergence 
        print("Pas de convergence (iterations_SW). Erreur pour la dernière itération : ",sum(norm.(prevInitialValues-InitialValues)),"\n")
        # InitialValues[n+1]=0 # pour donner un social welfare de 0, pour pas interférer
    end
end

function iterations_Gain_for_i(InitialValues,G,k,delta,eps,seuil) # fait des itérations jusqu'à obtenir un équilibre, avec comme critère de convergence eps
    n,_,_,_,_,_=G
    
    prevInitialValues = []
    nb_iter = 1

    while nb_iter <= seuil
        prevInitialValues = copy(InitialValues)
        one_iteration_Gain_for_i(InitialValues,G,k,delta)
        sum(norm.(prevInitialValues-InitialValues)) > eps || break 
        nb_iter+=1
    end

    if nb_iter > seuil  # pas de convergence 
        print("Pas de convergence (iterations_Gain_for_i). Erreur pour la dernière itération : ",sum(norm.(prevInitialValues-InitialValues)),"\n")
        InitialValues[n+1]=0 # pour donner un social welfare de 0, pour pas interférer 
    end
end

function seesaw_Q(InitialValues,G,k,delta,eps,seuil)
    iterations_SW(InitialValues,G,k,delta,eps,seuil)
    iterations_Gain_for_i(InitialValues,G,k,delta,eps,seuil)
    return SW(InitialValues,G) # 0 si on a pas eu de convergence
end


























# see-saw pour Qcorr

function seesaw_NashConstraints_for_i_for_ti(i,ti,param,G)  # Contraintes d'équilibre de Nash pour le joueur i lorsqu'il reçoit le type ti
    n,T,W,v0,v1,Pinit = G

    T_such_that = [] # Ensemble des types t tels que t[i] = ti, c'est ce sur quoi on va sommer 
    for t in T 
        if t[i]==ti 
            push!(T_such_that,t)
        end
    end

    A = enum_tuples(n)

    Gain_for_i = sum(u(i,a,t,G) * Pinit(t) * P(a,t,param,G) for t in T_such_that for a in A) # gain actuel pour le joueur i : ce qu'on veut optimiser

    mus = [[0,0],[0,1],[1,0],[1,1]] # fonctions mu_i possibles, de A_i -> A_i, identifiées à des tableaux 

    Deviated_Gain_for_i(ri,mui) = sum(u(i,change_tuple(a,i,mui[a[i]+1],n),t,G) * P(a,change_tuple(t,i,ri,n),param,G) * Pinit(t) for t in T_such_that for a in A) # Gain pour i lorsqu'il dévie de sa stratégie

    return [Gain_for_i - Deviated_Gain_for_i(ri,mui) for mui in mus for ri in 0:1] # on veut que le gain actuel pour i soit toujours meilleur que lorsqu'il dévie de sa stratégie
end

function seesaw_NashConstraints_for_i(i,param,G) # Contraintes d'équilibre de Nash pour le joueur i 
    return [seesaw_NashConstraints_for_i_for_ti(i,0,param,G) ; seesaw_NashConstraints_for_i_for_ti(i,1,param,G)]
end

function seesaw_NashConstraints(param,G) # Liste de toutes les contraintes d'équilibre de Nash
    n,_,_,_,_,_ = G
    constraints = []
    
    for i in 1:n 
        constraints = [constraints ; seesaw_NashConstraints_for_i(i,param,G)]
    end

    return constraints # concaténation de toutes les contraintes pour chaque joueur i
end

function one_iteration_Qcorr(InitialValues,G,k,delta) # Applique une itération du see-saw, et modifie InitialValues avec les nouvelles valeurs. 
    n,_,_,_,_,_ = G

    # SDP sur rho 
    current_SW = SW(InitialValues,G)
    current_rho = InitialValues[n+1]

    model = Model(SCS.Optimizer)
    set_silent(model) 

    @variable(model,rho[1:(k^n),1:(k^n)],PSD) 
    set_start_value.(rho,InitialValues[n+1])

    @constraint(model,LinearAlgebra.tr(rho)==1)

    param = [if j==n+1 rho else InitialValues[j] end for j in 1:(n+1)]
    @objective(model,Max,SW(param,G))
    JuMP.optimize!(model)
    InitialValues[n+1] = JuMP.value(rho)

    if SW(InitialValues,G) < current_SW + delta # si le sw n'est pas assez amélioré, on reste sur la solution qu'on avait de base 
        InitialValues[n+1] = current_rho 
    end


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

        # Contraintes que la corrélation soit un équilibre pour le joueur i 
        for cons in seesaw_NashConstraints_for_i(i,param,G)
            @constraint(model,cons >= 0)
        end

        @objective(model,Max,SW(param,G))
        JuMP.optimize!(model)
        InitialValues[i] = [[JuMP.value(N00),JuMP.value(N10)] [JuMP.value(N01),JuMP.value(N11)]]
    end
end

function iterations_Qcorr(InitialValues,G,k,delta,eps,seuil) # fait des itérations jusqu'à obtenir un équilibre, avec comme critère de convergence eps
    n,_,_,_,_,_=G
    
    prevInitialValues = []
    nb_iter = 1

    while nb_iter <= seuil
        prevInitialValues = copy(InitialValues)
        one_iteration_Qcorr(InitialValues,G,k,delta)
        sum(norm.(prevInitialValues-InitialValues)) > eps || break 
        nb_iter+=1
    end

    if nb_iter > seuil  # pas de convergence 
        print("Pas de convergence (iterations_Qcorr). Erreur pour la dernière itération : ",sum(norm.(prevInitialValues-InitialValues)),"\n")
        InitialValues[n+1]=0 # pour donner un social welfare de 0, pour pas interférer 
    end
end

function seesaw_Qcorr(InitialValues,G,k,delta,eps,seuil)
    iterations_Qcorr(InitialValues,G,k,delta,eps,seuil)
    return SW(InitialValues,G) # 0 si on a pas eu de convergence
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

function many_tests(myseesaw,G,k,povms,nb_aleas,barre) # Fait plusieurs tests en partant de POVMs aléatoires et de POVMs donnés par la liste povms. Affiche ceux dont le social welfare dépasse barre
    best_sw = 0
    
    # POVMs aléatoires
    print("POVMs aléatoires :\n")
    for i in 1:nb_aleas
        InitialValues = random_InitialValues(G,k)
        sw = myseesaw(InitialValues,G)
        if sw>barre
            print(SW(InitialValues,G))
            print("\n")
        end
        if sw > best_sw 
            best_sw = sw 
        end
    end

    # POVMs donnés 
    print("POVMs contrôlés :\n")
    for InitialValues in given_InitialValues(G,k,povms)
        sw = myseesaw(InitialValues,G)
        if sw>barre
            print(SW(InitialValues,G))
            print("\n")
        end
        if sw > best_sw 
            best_sw = sw 
        end
    end

    return best_sw
end

function controlled_gate(n,i,j,U) # Matrice U contrôlée sur un système de n parties (l'espace total est donc de dimension 2^n). Le bit de contrôle est i, le bit cible est j. 
    @assert i!=j
    
    function to_kron_first_term(player) 
        if player == i 
            return povm_canonique[1]
        else 
            return Matrix(I,2,2)
        end
    end

    function to_kron_second_term(player)
        if player == i 
            return povm_canonique[2]
        elseif player == j 
            return U 
        else 
            return Matrix(I,2,2)
        end
    end

    return krons([to_kron_first_term(player) for player in 1:n]) + krons([to_kron_second_term(player) for player in 1:n])
end

function pseudo_telepathy(n)  # renvoie l'état et les mesures de la stratégie pseudo-télépathique pour le graphe cyclique à n joueurs 
    proj_plus = povm_hadamard[1] # |+><+|

    CZ(i,j) = controlled_gate(n,i,j,Z)

    rho = prod(CZ(i,i+1) for i in 1:(n-1))*CZ(n,1) * krons([proj_plus for i in 1:n]) * CZ(n,1)*prod(CZ(i,i+1) for i in (n-1):(-1):1)

    measures = [povm_canonique povm_hadamard]

    return [if i==n+1 rho else copy(measures) end for i in 1:(n+1)]
end

