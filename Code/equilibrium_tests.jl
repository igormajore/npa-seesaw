#=

Le but de ce code est de vérifier qu'une solution quantique donnée QSol (au sens de seesaw.jl) est un équilibre (classique, quantique, ou mixte), à une approximation ApproxStab près.

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


















# Tests d'équilibre 

function is_ClassicalEquilibrium_for_i(i,QSol,G,ApproxStab)
    n,T,W,_,_,Pinit = G

    complet_types(ti) = filter(t -> (t[i]==ti),T) # Ensemble des types t tels que t[i] = ti, c'est ce sur quoi on va sommer 
    A = enum_tuples(n)

    Current_Gain_for_i(ti) = sum(u(i,a,t,G) * Pinit(t) * P(a,t,QSol,G) for t in complet_types(ti) for a in A) # gain actuel pour le joueur i : ce qu'on veut optimiser

    Deviated_Gain_for_i(ti,ri,mui) = sum(u(i,change_tuple(a,i,mui[a[i]+1],n),t,G) * P(a,change_tuple(t,i,ri,n),QSol,G) * Pinit(t) for t in complet_types(ti) for a in A) # Gain pour i lorsqu'il dévie de sa stratégie

    for ti in 0:1
        for ri in 0:1
            for mui in [[0,0],[0,1],[1,0],[1,1]] # fonctions mu_i possibles, de A_i -> A_i, identifiées à des tableaux 
                if Deviated_Gain_for_i(ti,ri,mui) > Current_Gain_for_i(ti) + ApproxStab 
                    #print("La solution proposée n'est pas un équilibre : le joueur ",i,", qui est classique, peut améliorer son gain lorsqu'il a le type ",ti,", en effectuant la déviation ri = ",ri," et mui = ",mui,", faisant passer son gain de ",Current_Gain_for_i(ti)," à ",Deviated_Gain_for_i(ti,ri,mui),".\n")
                    return false
                end
            end
        end
    end
    
    return true
end

function is_QuantumEquilibrium_for_i(i,QSol,G,ApproxStab)
    n,_,_,_,_,_ = G
    k,_ = LinearAlgebra.size(QSol[1][1,1])

    model=Model(SCS.Optimizer)
    set_silent(model)

    @variable(model,N00[1:k,1:k],PSD)

    @variable(model,N10[1:k,1:k],PSD)

    @variable(model,N01[1:k,1:k],PSD)

    @variable(model,N11[1:k,1:k],PSD)

    @constraint(model,N00+N10==LinearAlgebra.I)
    @constraint(model,N01+N11==LinearAlgebra.I)
        
    N = [[N00,N10] [N01,N11]]

    QSol_avec_variable = [if j==i N else QSol[j] end for j in 1:(n+1)]
    @objective(model,Max,Gain_for_i(i,QSol_avec_variable,G))
    JuMP.optimize!(model)

    if objective_value(model) > Gain_for_i(i,QSol,G) + ApproxStab
        #print("La solution proposée n'est pas un équilibre : le joueur ",i,", qui est quantique, peut améliorer son gain en utilisant les mesures ",JuMP.value.(N)," faisant passer son gain de ",Gain_for_i(i,QSol,G)," à ",objective_value(model),".\n")
        return false
    end

    return true
end

function is_Equilibrium(QSol,G,m,ApproxStab)
    n,_,_,_,_,_ = G 

    is_eq = true 

    for i in 1:m 
        is_eq = is_QuantumEquilibrium_for_i(i,QSol,G,ApproxStab) && is_eq
    end

    for i in (m+1):n 
        is_eq = is_ClassicalEquilibrium_for_i(i,QSol,G,ApproxStab) && is_eq
    end

    return is_eq
end
























