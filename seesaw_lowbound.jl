using LinearAlgebra, JuMP, SCS 

#=

Un jeu G est un uplet G=(n,T,W,v0,v1,Pinit) où n est le nombre de joueurs, T est l'ensemble des types de départ autorisés, W est la liste des situations gagnantes (c'est une liste d'uplets (a,t)), v0 et v1 sont les valeurs
de gain, Pinit est la distribution de probabilité des types des joueurs (c'est une fonction de T -> [0,1]).

La variable MODEL fait référence à un couple (model,param) où model est un modèle de JuMP et param est un tableau de taille n+1 qui contient les variables et paramètres du problème. Plus précisément, param[n+1] contient rho 
et pour 1<=i<=n, param[i] contient une matrice 2*2 M telle que M[a_i,t_i] soit la matrice de mesure du joueur i pour a_i|t_i. Pour chaque modèle, exactement une case de param contient de "vraies" variables et les autres ont 
une valeur de paramètre.

La variable MODELS contient une liste de tous les modèles MODEL du problème. Plus précisément, MODELS[n+1] est le problème où la variable est rho, et MODELS[i] pour 1<=i<=n est le modèle où les variables sont les mesures du 
joueur i.

La variable InitialValues a la même structure que param mais contient de vraies matrices

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









# Fonctions de gain

function u(i,a,t,G) # Fonction de gain du joueur i 
    _,_,W,v0,v1,_ = G 
    if (a,t) in W 
        return (if a[i]==0 v0 else v1 end)
    else 
        return 0 
    end
end

function SW(param,G)
    n,_,W,_,_,Pinit = G 

    to_kron(a,t) = [param[i][a[i]+1,t[i]+1] for i in 1:n] # liste des matrices dont on va calculer le produit tensoriel

    return (1/n) * sum(u(i,a,t,G) * Pinit(t) * LinearAlgebra.tr(krons(to_kron(a,t))) for (a,t) in W for i in 1:n)
end
















# Initialisation de modèles

function init_MODELS_SWopt(G,k) # G : un jeu, k : taille des matrices de mesures. La fonction construit MODELS lorsqu'on veut optimiser le social welfare
    n,_,_,_,_,_ = G
    
    MODELS = []

    for i in 1:n  # Création du MODEL où les variables sont les mesures pour le joueur i
        model = Model(SCS.Optimizer)

        param = [] # Initialisation des variables et des paramètres
        for j in 1:n
            if j==i # ce sont les variables
                N00 = @variable(model,[1:k,1:k],PSD)
                N10 = @variable(model,[1:k,1:k],PSD)
                N01 = @variable(model,[1:k,1:k],PSD)
                N11 = @variable(model,[1:k,1:k],PSD)
                N = [[N00, N10] [N01,N11]]
                push!(param,N)

                @constraint(model,N00 + N10 == LinearAlgebra.I)
                @constraint(model,N01 + N11 == LinearAlgebra.I)
            else # Ce sont des paramètres, initialisés à 0
                M00 = @variable(model,[1:k,1:k] in Parameter(0))
                M10 = @variable(model,[1:k,1:k] in Parameter(0))
                M01 = @variable(model,[1:k,1:k] in Parameter(0))
                M11 = @variable(model,[1:k,1:k] in Parameter(0))
                M = [[M00, M10] [M01,M11]]
                push!(param,M)
            end
        end
        rho = @variable(model,[1:(k*n),1:(k*n)] in Parameter(0))
        push!(param,rho)

        @objective(model,Max,SW(param,G))

        MODEL = (model,param)
        push!(MODELS,MODEL)
    end

    # Création du MODEL où la variable est rho  
    model = Model(SCS.Optimizer)
    param = []
    for j in 1:n   # Toutes les mesures sont des paramètres
        M00 = @variable(model,[1:k,1:k] in Parameter(0))
        M10 = @variable(model,[1:k,1:k] in Parameter(0))
        M01 = @variable(model,[1:k,1:k] in Parameter(0))
        M11 = @variable(model,[1:k,1:k] in Parameter(0))
        M = [[M00, M10] [M01,M11]]
        push!(param,M)
    end
    rho = @variable(model,[1:(k*n),1:(k*n)],PSD) # La variable est rho
    push!(param,rho)
    @constraint(model,LinearAlgebra.tr(rho)==1)

    @objective(model,Max,SW(param,G))

    MODEL = (model,param)
    push!(MODELS,MODEL)

    return MODELS
end























# Calcul effectif des SDP

function one_iteration(MODELS,InitialValues) # Applique une itération du see-saw, et modifie InitialValues avec les nouvelles valeurs

    # see-saw sur rho
    model,param = MODELS[n+1]
    for i in 1:n # initialisation des paramètres 
        set_parameter_value.(param[i][1,1],InitialValues[i][1,1])
        set_parameter_value.(param[i][2,1],InitialValues[i][2,1])
        set_parameter_value.(param[i][1,2],InitialValues[i][1,2])
        set_parameter_value.(param[i][2,2],InitialValues[i][2,2])
    end

    JuMP.optimize!(model)
    InitialValues[n+1] = JuMP.value(rho)


    # see-saw sur les mesures de chaque joueur, dans l'ordre croissant 
    for i in 1:n 
        model,param = MODELS[i]

        for j in 1:n # initialisation des paramètres 
            if j!=i 
                set_parameter_value.(param[j][1,1],InitialValues[j][1,1])
                set_parameter_value.(param[j][2,1],InitialValues[j][2,1])
                set_parameter_value.(param[j][1,2],InitialValues[j][1,2])
                set_parameter_value.(param[j][2,2],InitialValues[j][2,2])
            end
        end
        set_parameter_value.(param[n+1],InitialValues[n+1])

        JuMP.optimize!(model)
        InitialValues[i] = JuMP.value.(param[i])
    end
end












function bad_initial_values(G,k)
    n,_,_,_,_,_ = G

    InitialValues=[]
    for i in 1:n 
        M00 = I 
        M10 = zeros(k,k)
        M01 = I
        M11 = zeros(k,k)
        push!(InitialValues,[[M00, M10] [M01, M11]])
    end
    push!(InitialValues,(1/n)*I)
    
    return InitialValues
end




