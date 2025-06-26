using Plots

#= 

Les différents fichiers sont indépendants (on peut donc ouvrir l'un d'eux sans avoir à ouvrir les autres). Le fichier main se contente de tous les inclure.

Le fichier main sert aussi à faire les plots

=#

include("game_gen.jl")
include("npa_upbound_Qcorr.jl")
include("clean_bad_seesaw_lowbound.jl")




function plot_NPAopt(G,lv,ratios) 
    upbounds = [NPAopt(change_game(G,rat),lv) for rat in ratios]
    plot(ratios,upbounds,title="Majorant NPA",xlabel="v0/(v0+v1)")
end





function plot_seesaw(G,Param,m,k,NbAlea,povms,ratios) 
    # Pour les valeurs initiales : 
    #       - une est la meilleure pour le ratio v0/(v0+v1) précédent (pour la première valeur de ratio, une valeur initiale est choisie aléatoirement) ; 
    #       - NbAlea sont choisies aléatoirement ; 
    #       - on teste aussi tous les éléments de given_QSol(G,povms) 

    G = change_game(G,ratios[1])
    QSols = given_QSol(G,povms)
    lowbounds= []

    best_QSol = random_QSol(G,k)
    best_sw = SW(best_QSol,G)
    from_where = ("",(rand((false,true)),rand((false,true)),randfloat(0.5,1))) # Dit d'où vient le record qu'on a obtenu (record de l'itération précédente, aléatoire, given_QSol) et la Version correspondante


    for rat in ratios 
        # best_sw,best_InitialValues contiennent les records pour l'itération précédente

        print("\n\nNouveau ratio : v0/(v0+v1) = ",rat,"\n\n")
        G = change_game(G,rat)

        # meilleure pour le ratio précédent 
        print("\nMeilleure pour le ratio précédent\n")

        QSol = copy(best_QSol)
        _,Version = from_where

        sw = seesaw_mixte(QSol,G,Param,Version,m)
        print(sw,"\n")
        best_QSol = QSol # Mise à jour de best_QSol et best_sw, qui contiennent à présent le record pour l'itération actuelle
        best_sw = sw
        from_where = ("PRÉCÉDENTE",Version)
        


        # aléatoires 
        print("\nMesures aléatoires\n")
        for _ in 1:NbAlea
            QSol=random_QSol(G,k)
            Version = (rand((false,true)),rand((false,true)),randfloat(0.5,1))

            sw = seesaw_mixte(QSol,G,Param,Version,m)
            print(sw,"\n")
            if sw>best_sw
                best_QSol = QSol
                best_sw = sw
                from_where = ("ALÉATOIRE",Version)
                print("Amélioration : nouveau sw est ",best_sw,", obtenu par la version ",Version,"\n")
            end
        end


        # given_QSol(G,povms)
        print("\nMesures prescrites\n")
        to_do = copy.(QSols)
        for QSol in to_do 
            Version = (rand((false,true)),rand((false,true)),randfloat(0.5,1))

            sw = seesaw_mixte(QSol,G,Param,Version,m)
            print(sw,"\n")
            if sw>best_sw
                best_QSol = QSol
                best_sw = sw
                from_where = ("PRESCRITE",Version)
                print("Amélioration : nouveau sw est ",best_sw,", obtenu par la version ",Version,"\n")
            end
        end

        init,Version = from_where
        print("\nConclusion : le meilleur SW a été obtenu pour une valeur initiale dans ",init,", avec la version",Version,", donnant un record de ",best_sw,"\n")

        push!(lowbounds,best_sw)
    end

    #plot(ratios,lowbounds,title="Minorant see-saw",xlabel="v0/(v0+v1)")
    return lowbounds
end







