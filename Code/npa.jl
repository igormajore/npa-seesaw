using QuantumNPA
import Base: +  # pour ajouter des méthodes à l'addition

#= 

Un jeu G est un uplet G=(n,T,W,v0,v1,Pinit) où n est le nombre de joueurs, T est l'ensemble des types de départ autorisés, W est la liste des situations gagnantes (c'est une liste d'uplets (a,t)), v0 et v1 sont les valeurs
de gain, Pinit est la distribution de probabilité des types des joueurs (c'est une fonction de T -> [0,1]).

Note importante : dans W, les types t_i et a_i sont à valeurs dans 0,1, alors que les projecteurs de NPA sont indexés sur 1,2 à chaque fois, d'où le '+1' qui 
apparaît dans les expressions. 

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

function +(X::Polynomial,n::Number) # définit X+0 où X est un polynôme NPA et 0 est le nombre
    @assert n==0 
    return X
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

function SWnpa(G,P) # Expression du social welfare d'un jeu G sous la forme d'un polynôme NPA, en connaissant la corrélation P
    n,_,W,_,_,Pinit = G 

    return (1/n) * sum(u(i,a,t,G) * Pinit(t) * P(a,t) for (a,t) in W for i in 1:n)
end














# Contraintes d'équilibre de Nash

function npa_NashConstraints_for_i_for_ti(i,ti,G,P)  # Contraintes d'équilibre de Nash pour le joueur i lorsqu'il reçoit le type ti, en connaissant la corrélation P
    n,T,W,v0,v1,Pinit = G

    T_such_that = filter(t -> (t[i]==ti),T) # Ensemble des types t tels que t[i] = ti, c'est ce sur quoi on va sommer

    A = enum_tuples(n)

    Gain_for_i = sum(u(i,a,t,G) * Pinit(t) * P(a,t) for t in T_such_that for a in A) # gain actuel pour le joueur i : ce qu'on veut optimiser

    mus = [[0,0],[0,1],[1,0],[1,1]] # fonctions mu_i possibles, de A_i -> A_i, identifiées à des tableaux 

    Deviated_Gain_for_i(ri,mui) = sum(u(i,change_tuple(a,i,mui[a[i]+1],n),t,G) * P(a,change_tuple(t,i,ri,n)) * Pinit(t) for t in T_such_that for a in A) # Gain pour i lorsqu'il dévie de sa stratégie

    return [Gain_for_i - Deviated_Gain_for_i(ri,mui) for mui in mus for ri in 0:1] # on veut que le gain actuel pour i soit toujours meilleur que lorsqu'il dévie de sa stratégie
end

function npa_NashConstraints_for_i(i,G,P)
    return [npa_NashConstraints_for_i_for_ti(i,0,G,P) ; npa_NashConstraints_for_i_for_ti(i,1,G,P)]
end

function npa_NashConstraints(G,P) # Liste de toutes les contraintes d'équilibre de Nash, sous forme d'une liste de polynômes NPA qui doivent être positifs, comme attendu par la fonction npa_max 
    n,_,_,_,_,_ = G
    constraints = []
    
    for i in 1:n 
        constraints = [constraints ; npa_NashConstraints_for_i(i,G,P)]
    end

    return constraints # concaténation de toutes les contraintes pour chaque joueur i, chaque type ti
end 
















# Optimisation NPA 

function NPAopt(G,lv)
    n,_,_,_,_ = G

    # Définition des projecteurs pour chaque joueur
    M = [projector(i,1:2,1:2,full=true) for i in 1:n]

    # Définition de la corrélation
    P(a,t) = prod(M[j][a[j]+1,t[j]+1] for j in 1:n)

    return npa_max(SWnpa(G,P),lv,ge=npa_NashConstraints(G,P)) # On optimise le social welfare sous la contrainte "équilibre de Nash"
end












































