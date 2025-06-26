

#=

Un jeu G est un uplet G=(n,T,W,v0,v1,Pinit) où n est le nombre de joueurs, T est l'ensemble des types de départ autorisés, W est la liste des situations gagnantes (c'est une liste d'uplets (a,t)), v0 et v1 sont les valeurs
de gain, Pinit est la distribution de probabilité des types des joueurs (c'est une fonction de T -> [0,1]).

=#




# Fonctions auxiliaires

function cyclic_idx(k,n) # modulo décalé : pour n=5, k=0,1,..., on a la séquence 5,1,2,3,4,5,1,2,3,...
    return mod(k-1,n) + 1
end

function enum_tuples(n) # Énumère {0,1}^n sous forme d'une liste de n-uplets
    if n==0
       return [()]
    else 
       return vcat([(0,a...) for a in enum_tuples(n-1)],[(1,a...) for a in enum_tuples(n-1)])
    end
end













# Initialisation d'un jeu 

function NC_00_cycle(n,rat) # Initialisation NC_00(C_n). rat=v0/(v0+v1), et on met v0+v1=2 
    T = vcat([ntuple(j -> if j==i 1 else 0 end,n) for i in 1:n],[ntuple(_ -> 1,n)]) # Ensemble des types (inputs), tout autre type est de probabilité 0
    A = enum_tuples(n) # Ensemble des actions (outputs) = {0,1}^n 

    # Construction de W : on parcourt A et si a dans A satisfait une winning condition pour un type t, on ajoute (a,t) à W
    win1(i,a) = ((a[cyclic_idx(i-1,n)] + a[i] + a[cyclic_idx(i+1,n)]) % 2 == 0)
    win2(a) = (sum(a[i] for i in 1:n) % 2 == 1)

    W = []
    for a in A 
        for i in 1:n 
            if win1(i,a)
                t = ntuple(j -> if j==i 1 else 0 end,n)
                push!(W,(a,t))
            end
        end
        if win2(a)
            t = ntuple(_ -> 1,n)
            push!(W,(a,t))
        end
    end

    v0 = 2rat 
    v1 = 2-v0 


    Pinit(t) =              # Distribution initiale : loi uniforme sur T
        if t in T 
            return 1/(n+1)
        else 
            return 0 
        end

    return (n,T,W,v0,v1,Pinit)
end









# Modification des gains

function change_game(G,rat)
    n,T,W,_,_,Pinit = G
    v0 = 2rat
    v1=2-v0
    return (n,T,W,v0,v1,Pinit)
end