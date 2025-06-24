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
    best_sw = SW(QSol,G)

    for rat in ratios 
        # best_sw,best_InitialValues contiennent les records pour l'itération précédente

        print("\n\nNouveau ratio : v0/(v0+v1) = ",rat,"\n\n")
        G = change_game(G,rat)

        # meilleure pour le ratio précédent 
        print("\nMeilleure pour le ratio précédent\n")
        QSol = copy(best_QSol)
        sw = seesaw_mixte(QSol,G,Param,m)
        print(sw,"\n")

        best_QSol = QSol # Mise à jour de best_QSol et best_sw, qui contiennent à présent le record pour l'itération actuelle
        best_sw = sw


        # aléatoires 
        print("\nMesures aléatoires\n")
        for _ in 1:NbAlea
            QSol=random_QSol(G,k)
            sw = seesaw_mixte(QSol,G,Param,m)
            print(sw,"\n")
            if sw>best_sw
                best_QSol = QSol
                best_sw = sw
                print("Amélioration : nouveau sw est ",best_sw,"\n")
            end
        end


        # given_QSol(G,povms)
        print("\nMesures prescrites\n")
        to_do = copy.(QSols)
        for QSol in to_do 
            sw = seesaw_mixte(QSol,G,Param,m)
            print(sw,"\n")
            if sw>best_sw
                best_QSol = QSol
                best_sw = sw
                print("Amélioration : nouveau sw est ",best_sw,"\n")
            end
        end

        push!(lowbounds,best_sw)
    end

    plot(ratios,lowbounds,title="Minorant see-saw",xlabel="v0/(v0+v1)")
end





