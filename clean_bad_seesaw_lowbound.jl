using LinearAlgebra, Ket, JuMP, SCS 

#=

Ceci est une version non optimisée de seesaw_lowbound.jl. Dans la version initiale, je voulais utiliser des paramètres. Le problème est que si SCS (avec ParametricOptInterface) supporte les expressions de type Paramètre * Variable 
(ie ne les considère pas comme quadratiques), il ne supporte pas Paramètre * Paramètre * Variable (il considère que c'est une expression non linéaire). Je n'ai donc pas réussi à utiliser les paramètres. Dans cette version, je crée donc 
un nouveau modèle pour chaque problème que je veux résoudre. 

Ceci est une version plus propre et lisible de bad_seesaw_lowbound. Les noms des paramètres sont plus simples à comprendre et j'ai enlevé les fonctions inutiles. 





Un jeu G est un uplet G=(n,T,W,v0,v1,Pinit) où n est le nombre de joueurs, T est l'ensemble des types de départ autorisés, W est la liste des situations gagnantes (c'est une liste d'uplets (a,t)), v0 et v1 sont les valeurs
de gain, Pinit est la distribution de probabilité des types des joueurs (c'est une fonction de T -> [0,1]).

m (0<=m<=n) est le nombre de joueurs ayant un accès quantique aux ressources : les joueurs 1<=i<=m ont un accès quantique et les joueurs m+1<=i<=n ont un accès classique. m=n correspond à Q(G) et m=0 à Qcorr(G).

k est la dimension des espaces des joueurs (en général, k=2). L'espace total est donc de dimension k^n. 

Une solution quantique est représentée par un tableau QSol de taille n+1 tel que : 
    - pour 1<=i<=n, QSol[i] est une matrice 2*2 telle que QSol[i][a_i,t_i] soit la mesure du joueur i, pour l'action a_i sachant qu'il a le type t_i (matrices d'ordre k). 
    - QSol[n+1] contient rho

Approx = (ApproxStab,ApproxConv) contient les approximations numériques que le programme fait. ApproxStab est une approximation dont le but est de stabiliser le code, en l'empêchant de changer de solution si celle qu'il a est déjà optimale. 
Concrètement, le programme ne change de solution que si la nouvelle lui permet d'améliorer son gain d'au moins ApproxStab (c'est-à-dire qu'il ne change de solution que si cela améliore son gain de façon significative). ApproxConv est le critère 
de convergence : on considère qu'on a obtenu la convergence (et qu'on peut donc arrêter d'itérer) lorsque deux solutions successives diffèrent d'au plus ApproxConv. Dans la pratique, j'ai souvent pris ApproxStab = ApproxConv = 1e-6.

NbTentatives est le nombre d'itérations qu'on s'autorise avant de considérer qu'on n'a pas réussi à obtenir la convergence. J'ai généralement pris NbTentatives = 200. 

Les paramètres de la résolution numérique sont regroupés dans une variable Param = (Approx,NbTentatives).

J'ai aussi introduit une variable Version = (do_StabC : bool,do_OptRho : bool,to_OptC : Float dans [0,1]) qui permet de tester différentes versions de l'algorithme. do_StabC signifie qu'on active un stabilisateur pour les joueurs classiques. Lorsque 
cette variable vaut true, chaque joueur classique i commence par vérifier si QSol est un équilibre (corrélé) pour i, et ne fait rien si c'est le cas. do_OptRho est une variable qui indique si on fait l'optimisation sur rho dans iterations_equilibre. 
to_OptC permet de déterminer la fonction qu'on optimise pour les itérations sur les joueurs classiques : il s'agit de to_OptC*SW + (1-to_OptC)*Gain_for_i. 


Comme les actions et types sont à valeurs dans {0,1} mais les indices d'array sont dans {1,...}, il y a un décalage d'indices.







L'algorithme est le suivant : 
    Étape 1 : Optimiser SW sans contrainte (fonction iterations_SW).
    Étape 2 : Les joueurs 1<=i<=m optimisent leur propre gain (fonction Gain_for_i) sans contrainte, et les joueurs m+1<=i<=n optimisent SW (ou une combinaison convexe de SW et Gain_for_i, cf ci-dessus) sous la contrainte d'être un équilibre corrélé pour i (fonction iterations_equilibre). 

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

function krons(lst) # Calcule le produit tensoriel de toutes les matrices de lst, supposée non vide
    n = length(lst)
    @assert n>0
    prod = lst[1]
    for i in 2:n 
        prod = kron(prod,lst[i])
    end
    return prod 
end

function randfloat(a,b) # Nombre aléatoire flottant entre a et b
    return a + (b-a)*rand(Float64)
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

function P(a,t,QSol,G) # corrélation P(a|t) correspondant à la solution quantique QSol
    n,_,_,_,_,_ = G
    to_kron(a,t) = [QSol[i][a[i]+1,t[i]+1] for i in 1:n] # liste des matrices dont on va calculer le produit tensoriel
    return LinearAlgebra.tr(QSol[n+1] * krons(to_kron(a,t)))
end

function SW(QSol,G)
    n,_,W,_,_,Pinit = G 


    return (1/n) * sum(u(i,a,t,G) * Pinit(t) * P(a,t,QSol,G) for (a,t) in W for i in 1:n)
end

function Gain_for_i(i,QSol,G)
    _,_,W,_,_,Pinit = G 

    return sum(u(i,a,t,G) * Pinit(t) * P(a,t,QSol,G) for (a,t) in W)
end



















# Contraintes d'équilibres de Nash corrélé, c'est-à-dire, pour les joueurs ayant un accès classique aux ressources

function seesaw_NashConstraints_for_i_for_ti(i,ti,QSol,G)  # Contraintes d'équilibre de Nash pour le joueur i lorsqu'il reçoit le type ti (m+1<=i<=n)
    n,T,W,v0,v1,Pinit = G

    T_such_that = [] # Ensemble des types t tels que t[i] = ti, c'est ce sur quoi on va sommer 
    for t in T 
        if t[i]==ti 
            push!(T_such_that,t)
        end
    end

    A = enum_tuples(n)

    Gain_for_i = sum(u(i,a,t,G) * Pinit(t) * P(a,t,QSol,G) for t in T_such_that for a in A) # gain actuel pour le joueur i : ce qu'on veut optimiser

    mus = [[0,0],[0,1],[1,0],[1,1]] # fonctions mu_i possibles, de A_i -> A_i, identifiées à des tableaux 

    Deviated_Gain_for_i(ri,mui) = sum(u(i,change_tuple(a,i,mui[a[i]+1],n),t,G) * P(a,change_tuple(t,i,ri,n),QSol,G) * Pinit(t) for t in T_such_that for a in A) # Gain pour i lorsqu'il dévie de sa stratégie

    return [Gain_for_i - Deviated_Gain_for_i(ri,mui) for mui in mus for ri in 0:1] # on veut que le gain actuel pour i soit toujours meilleur que lorsqu'il dévie de sa stratégie
end

function seesaw_NashConstraints_for_i(i,QSol,G) # Contraintes d'équilibre de Nash pour le joueur i (m+1<=i<=n)
    return [seesaw_NashConstraints_for_i_for_ti(i,0,QSol,G) ; seesaw_NashConstraints_for_i_for_ti(i,1,QSol,G)]
end

function is_equilibrium_for_i(i,QSol,G,Param) # Teste si la solution vérifie les contraintes d'équilibre pour le joueur i, à ApproxStab près
    ((ApproxStab,_),_) = Param 

    for c in seesaw_NashConstraints_for_i(i,QSol,G)
        if c < -ApproxStab 
            return false 
        end
    end

    return true
end






























# Fonctions du see-saw 

function une_iteration_SW(QSol,G,Param) # Applique une itération du see-saw, en optimisant SW sans contrainte. La variable QSol est modifiée pour contenir la solution quantique obtenue à l'issue de l'itération.
    n,_,_,_,_,_ = G
    (ApproxStab,_),_ = Param 
    k,_ = LinearAlgebra.size(QSol[1][1,1])

    # SDP sur rho 
    current_SW = SW(QSol,G) 
    current_rho = QSol[n+1]

    model = Model(SCS.Optimizer)
    set_silent(model) 

    @variable(model,rho[1:(k^n),1:(k^n)],PSD) 
    set_start_value.(rho,QSol[n+1])

    @constraint(model,LinearAlgebra.tr(rho)==1)

    QSol_avec_variable = [if j==n+1 rho else QSol[j] end for j in 1:(n+1)]
    @objective(model,Max,SW(QSol_avec_variable,G))
    JuMP.optimize!(model)
    QSol[n+1] = JuMP.value(rho)

    if SW(QSol,G) < current_SW + ApproxStab # si le sw n'est pas assez amélioré, on reste sur la solution qu'on avait initialement 
        QSol[n+1] = current_rho 
    end


    # SDP sur les mesures de chaque joueur, dans l'ordre croissant 
    for i in 1:n
        current_SW = SW(QSol,G)
        current_Mi = QSol[i]

        model = Model(SCS.Optimizer)
        set_silent(model)

        @variable(model,N00[1:k,1:k],PSD)
        set_start_value.(N00,QSol[i][1,1])

        @variable(model,N10[1:k,1:k],PSD)
        set_start_value.(N10,QSol[i][2,1])

        @variable(model,N01[1:k,1:k],PSD)
        set_start_value.(N01,QSol[i][1,2])

        @variable(model,N11[1:k,1:k],PSD)
        set_start_value.(N11,QSol[i][2,2])

        @constraint(model,N00+N10==LinearAlgebra.I)
        @constraint(model,N01+N11==LinearAlgebra.I)

        N = [[N00,N10] [N01,N11]]
        QSol_avec_variable = [if j==i N else QSol[j] end for j in 1:(n+1)]
        @objective(model,Max,SW(QSol_avec_variable,G))
        JuMP.optimize!(model)
        QSol[i] = [[JuMP.value(N00),JuMP.value(N10)] [JuMP.value(N01),JuMP.value(N11)]]

        if SW(QSol,G) < current_SW + ApproxStab # si le sw n'est pas assez amélioré, on reste sur la solution qu'on avait initialement
            QSol[i] = current_Mi
        end
    end
end

function iterations_SW(QSol,G,Param) # Fait des itérations jusqu'à obtenir un SW localement optimal, avec comme critère de convergence ApproxConv
    n,_,_,_,_,_ = G
    (_,ApproxConv),NbTentatives = Param
    
    prevQSol = []
    nb_iter = 1

    while nb_iter <= NbTentatives
        prevQSol = copy(QSol)
        une_iteration_SW(QSol,G,Param)
        sum(norm.(prevQSol - QSol)) > ApproxConv || break 
        nb_iter+=1
    end
end

function une_iteration_equilibre(QSol,G,Param,Version,m) # Applique une itération du see-saw pour obtenir un équilibre, et modifie QSol avec la nouvelle solution obtenue.  
    n,_,_,_,_,_ = G
    (ApproxStab,_),_ = Param 
    (do_StabC,do_OptRho,to_OptC) = Version
    k,_ = LinearAlgebra.size(QSol[1][1,1])

    if do_OptRho      # SDP sur rho 
        current_SW = SW(QSol,G) 
        current_rho = QSol[n+1]

        model = Model(SCS.Optimizer)
        set_silent(model) 

        @variable(model,rho[1:(k^n),1:(k^n)],PSD) 
        set_start_value.(rho,QSol[n+1])

        @constraint(model,LinearAlgebra.tr(rho)==1)

        QSol_avec_variable = [if j==n+1 rho else QSol[j] end for j in 1:(n+1)]
        @objective(model,Max,SW(QSol_avec_variable,G))
        JuMP.optimize!(model)
        QSol[n+1] = JuMP.value(rho)

        if SW(QSol,G) < current_SW + ApproxStab # si le sw n'est pas assez amélioré, on reste sur la solution qu'on avait initialement 
            QSol[n+1] = current_rho 
        end
    end
    


    # SDP sur les mesures de chaque joueur ayant un accès quantique, dans l'ordre croissant 
    for i in 1:m
        current_Gain_for_i = Gain_for_i(i,QSol,G)
        current_Mi = QSol[i]

        model=Model(SCS.Optimizer)
        set_silent(model)

        @variable(model,N00[1:k,1:k],PSD)
        set_start_value.(N00,QSol[i][1,1])

        @variable(model,N10[1:k,1:k],PSD)
        set_start_value.(N10,QSol[i][2,1])

        @variable(model,N01[1:k,1:k],PSD)
        set_start_value.(N01,QSol[i][1,2])

        @variable(model,N11[1:k,1:k],PSD)
        set_start_value.(N11,QSol[i][2,2])

        @constraint(model,N00+N10==LinearAlgebra.I)
        @constraint(model,N01+N11==LinearAlgebra.I)
        
        N = [[N00,N10] [N01,N11]]
        QSol_avec_variable = [if j==i N else QSol[j] end for j in 1:(n+1)]
        @objective(model,Max,Gain_for_i(i,QSol_avec_variable,G))
        JuMP.optimize!(model)
        QSol[i] = [[JuMP.value(N00),JuMP.value(N10)] [JuMP.value(N01),JuMP.value(N11)]]   

        if Gain_for_i(i,QSol,G) < current_Gain_for_i + ApproxStab # si le joueur i n'améliore pas assez son gain, on reste sur la solution qu'on avait de base 
            QSol[i] = current_Mi
        end
    end



    # SDP sur les mesures de chaque joueur ayant un accès classique, dans l'ordre croissant 
    for i in (m+1):n

        if !do_StabC || !is_equilibrium_for_i(i,QSol,G,Param) # Stabilisateur pour les classiques : on ne modifie la mesure que si on n'est pas déjà sur un équilibre

            model = Model(SCS.Optimizer)
            set_silent(model)

            @variable(model,N00[1:k,1:k],PSD)

            @variable(model,N10[1:k,1:k],PSD)

            @variable(model,N01[1:k,1:k],PSD)

            @variable(model,N11[1:k,1:k],PSD)

            @constraint(model,N00+N10==LinearAlgebra.I)
            @constraint(model,N01+N11==LinearAlgebra.I)

            N = [[N00,N10] [N01,N11]]
            QSol_avec_variable = [if j==i N else QSol[j] end for j in 1:(n+1)]

            # Contraintes que la corrélation soit un équilibre pour le joueur i 
            for cons in seesaw_NashConstraints_for_i(i,QSol_avec_variable,G)
                @constraint(model,cons >= 0)
            end

            @objective(model,Max,to_OptC*SW(QSol_avec_variable,G)+(1-to_OptC)*Gain_for_i(i,QSol_avec_variable,G))
            JuMP.optimize!(model)
            QSol[i] = [[JuMP.value(N00),JuMP.value(N10)] [JuMP.value(N01),JuMP.value(N11)]]
        end
    end
end

function iterations_equilibre(QSol,G,Param,Version,m) # Fait des itérations jusqu'à obtenir un équilibre, avec comme critère de convergence ApproxConv
    n,_,_,_,_,_ = G
    (_,ApproxConv),NbTentatives = Param
    
    prevQSol = []
    nb_iter = 1

    while nb_iter <= NbTentatives
        prevQSol = copy(QSol)
        une_iteration_equilibre(QSol,G,Param,Version,m)
        sum(norm.(prevQSol - QSol)) > ApproxConv || break 
        nb_iter+=1
    end

    if nb_iter > NbTentatives  # pas de convergence 
        print("Pas de convergence (iterations_equilibre). Erreur pour la dernière itération : ",sum(norm.(prevQSol - QSol)),"\n")
        QSol[n+1]=0 # pour donner un social welfare de 0, pour pas interférer 
    end
end

function seesaw_mixte(QSol,G,Param,Version,m)
    iterations_SW(QSol,G,Param)
    iterations_equilibre(QSol,G,Param,Version,m)
    return SW(QSol,G) # =0 si on n'a pas eu de convergence
end
































# Valeurs initiales 

function given_QSol(G,povms) # povms est une liste de POVMs. Cette fonction renvoie toutes les solutions quantiques qu'on peut construire en partant des mesures dans la liste povms (et en prenant des états partagés aléatoires)
    if povms==[]
        return []
    end
    
    n,_,_,_,_,_ = G
    k,_ = LinearAlgebra.size(povms[1][1])
    
    un_joueur = [[m1 m2] for m1 in povms for m2 in povms] # toutes les combinaisons de mesures possibles pour un joueur : une mesure pour chaque type 

    function combinaisons(p) # renvoie la liste de toutes les QSol à tester (sans l'état partagé, qui sera choisi aléatoirement ensuite) s'il y a p joueurs (pour récursion)
        if p==0 
            return [[]]
        else
            res = []
            for QSol in combinaisons(p-1)
                for m in un_joueur 
                    newQSol = copy(QSol)
                    push!(newQSol,m)
                    push!(res,newQSol)
                end
            end
            return res 
        end
    end
    
    AllQSol = combinaisons(n)

    # Ajout des états partagés aléatoires 
    for QSol in AllQSol
        push!(QSol,random_state(Float64,k^n))
    end 

    return AllQSol
end

function random_QSol(G,k)
    n,_,_,_,_,_ = G 

    QSol = []
    for i in 1:n 
        M0 = random_povm(Float64,k,2) # on prend des POVMs réels car les SDP ont pour variables des matrices symétriques positives
        M1 = random_povm(Float64,k,2)
        M = [M0 M1]
        push!(QSol,M)
    end
    rho = random_state(Float64,k^n)
    push!(QSol,rho) 
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


















# Fonctions de test

function many_tests(G,Param,m,k,NbAlea,povms) # Fait plusieurs tests en partant de POVMs aléatoires et de POVMs donnés par la liste povms. Renvoie le meilleur SW obtenu. 
    best_sw = 0
    
    # POVMs aléatoires
    print("\nPOVMs aléatoires :\n")
    for _ in 1:NbAlea
        QSol = random_QSol(G,k)
        sw = seesaw_mixte(QSol,G,Param,m)
        print(sw,"\n")
        if sw > best_sw 
            best_sw = sw 
        end
    end

    # POVMs donnés 
    print("\nPOVMs contrôlés :\n")
    for QSol in given_QSol(G,povms)
        sw = seesaw_mixte(QSol,G,Param,m)
        print(sw,"\n")
        if sw > best_sw 
            best_sw = sw 
        end
    end

    return best_sw
end

