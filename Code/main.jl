using Plots

#= 

Les différents fichiers sont indépendants (on peut donc ouvrir l'un d'eux sans avoir à ouvrir les autres). Le fichier main se contente de tous les inclure.

Le fichier main sert aussi à faire les plots

=#

include("game_generator.jl")
include("npa.jl")
include("seesaw.jl")
include("equilibrium_tests.jl")




function plot_NPAopt(G,lv,ratios) 
    upbounds = [NPAopt(change_game(G,rat),lv) for rat in ratios]
    #plot(ratios,upbounds,title="Majorant NPA",xlabel="v0/(v0+v1)")
    return upbounds
end





function plot_seesaw(G,Param,m,ratios;k=2,NbAlea=50,povms_demandés=[],QSols_demandés=[],state=[]) 
    n,_,_,_,_,_ = G
    # Pour les valeurs initiales : 
    #       - une est la meilleure pour le ratio v0/(v0+v1) précédent (pour la première valeur de ratio, une valeur initiale est choisie aléatoirement) ; 
    #       - QSols_demandés est un tableau de taille length(ratios). Pour chaque valeur rat de ratios, fait beaucoup de tests pour QSols_demandés[rat] avec des Version différentes
    #       - NbAlea sont choisies aléatoirement ; 
    #       - on teste aussi tous les éléments de given_QSol(G,povms_demandés) 

    # Si state!=[], tous les see-saw seront calculés en partant de l'état fixé state.

    constant_state = (state!=[])

    G = change_game(G,ratios[1])
    QSols_povms_demandés = given_QSol(G,povms_demandés)
    if constant_state
        for QSol in QSols_povms_demandés
            QSol[n+1]=state 
        end
    end

    best_QSols = [] # Les tableaux contenant tous les records pour chaque ratio
    best_Versions = []
    best_SWs = []


    best_QSol = random_QSol(G,k) # Initialisation des variables records
    if constant_state 
        best_QSol[n+1]=state 
    end
    best_Version = ("",(rand((false,true)),rand((false,true)),randfloat(0,1))) 
    best_SW = SW(best_QSol,G)

    compteur=1


    for rat in ratios 
        # best_QSol, best_Version, best_SW contiennent les records pour l'itération précédente

        print("\n\nNouveau ratio : v0/(v0+v1) = ",rat,"\n\n")
        G = change_game(G,rat)

        # Meilleure pour le ratio précédent 
        print("\nMeilleure pour le ratio précédent\n")

        QSol = copy(best_QSol)
        _,Version = best_Version

        sw = seesaw_mixte(QSol,G,Param,Version,m,constant_state=constant_state)
        print(sw,"\n")
        best_QSol = QSol # Mise à jour de best_QSol, best_Version, best_SW, qui contiennent à présent le record pour l'itération actuelle
        best_Version = ("PRÉCÉDENTE",Version)
        best_SW = sw

        
        # QSols_demandés 
        if QSols_demandés != []
            print("\nValeurs initiales demandées\n")
            for do_StabC in [false,true]
                for do_OptRho in [false,true]
                    for to_OptC in range(0,1,length=div(NbAlea,8))
                        QSol = copy(QSols_demandés[compteur])
                        if constant_state
                            QSol[n+1]=state 
                        end
                        Version = (do_StabC,do_OptRho,to_OptC)

                        sw = seesaw_mixte(QSol,G,Param,Version,m,constant_state=constant_state)
                        print(sw,"\n")
                        if sw>best_SW
                            best_QSol = QSol
                            best_Version = ("QSOL DEMANDÉE",Version)
                            best_SW = sw
                            print("Amélioration : nouveau sw est ",best_SW,", obtenu par la version ",Version,"\n")
                        end
                    end
                end
            end
        end

        # aléatoires 
        print("\nMesures aléatoires\n")
        for _ in 1:NbAlea
            QSol=random_QSol(G,k)
            if constant_state
                QSol[n+1]=state 
            end
            Version = (rand((false,true)),rand((false,true)),randfloat(0,1))

            sw = seesaw_mixte(QSol,G,Param,Version,m,constant_state=constant_state)
            print(sw,"\n")
            if sw>best_SW
                best_QSol = QSol
                best_Version = ("ALÉATOIRE",Version)
                best_SW = sw
                print("Amélioration : nouveau sw est ",best_SW,", obtenu par la version ",Version,"\n")
            end
        end


        # given_QSol(G,povms)
        print("\nMesures avec POVMs demandés\n")
        to_do = copy.(QSols_povms_demandés)
        for QSol in to_do 
            Version = (rand((false,true)),rand((false,true)),randfloat(0,1))

            sw = seesaw_mixte(QSol,G,Param,Version,m,constant_state=constant_state)
            print(sw,"\n")
            if sw>best_SW
                best_QSol = QSol
                best_Version = ("PRESCRITE",Version)
                best_SW = sw
                print("Amélioration : nouveau sw est ",best_SW,", obtenu par la version ",Version,"\n")
            end
        end

        init,Version = best_Version
        print("\nConclusion : le meilleur SW a été obtenu pour une valeur initiale dans ",init,", avec la version ",Version,", donnant un record de ",best_SW,"\n")

        push!(best_QSols,best_QSol)
        push!(best_Versions,best_Version)
        push!(best_SWs,best_SW)

        compteur+=1
    end

    #plot(ratios,best_sws,title="Minorant see-saw",xlabel="v0/(v0+v1)")
    return best_QSols,best_Versions,best_SWs
end






function plot_equilibria(ratios,thetas,ApproxStab)
    G=NC_00_cycle(3,0)
    best_ms = []

    for rat in ratios 
        print("\n\nv0/(v0+v1) = ",rat,"\n\n")
        G = change_game(G,rat)

        for theta in thetas 
            for u in 1:3 
                print("\ntheta = ",theta,", u = ",u,"\n")
                best_m = -1 
                for m in 0:3 
                    is_eq = is_Equilibrium(tilted_solution(3,u,theta),G,m,ApproxStab)

                    print("m = ",m," : ",is_eq,"\n")
                    if is_eq 
                        best_m=m 
                    end
                end

                push!(best_ms,best_m)
            end
        end
    end

    return best_ms
end







