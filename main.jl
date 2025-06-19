using Plots

#= 

Les différents fichiers sont indépendants (on peut donc ouvrir l'un d'eux sans avoir à ouvrir les autres). Le fichier main se contente de tous les inclure.

Le fichier main sert aussi à faire les plots

=#

include("game_gen.jl")
include("npa_upbound_Qcorr.jl")
include("bad_seesaw_lowbound.jl")




function plot_NPAopt(G,lv,ratios) 
    upbounds = [NPAopt(change_game(G,rat),lv) for rat in ratios]
    plot(ratios,upbounds,title="Majorant NPA",xlabel="v0/(v0+v1)")
end




function plot_seesaw(G,k,povms,eps,seuil,nb_alea,ratios)
    # Pour les valeurs initiales : 
    #       - une est la meilleure pour le ratio v0/(v0+v1) précédent (pour la première valeur de ratio, une valeur initiale est choisie aléatoirement) ; 
    #       - nb_alea sont choisies aléatoirement ; 
    #       - on teste aussi tous les éléments de given_InitialValues(G,k,povms) (tableau AllInitialValues)

    G = change_game(G,ratios[1])
    AllInitialValues = given_InitialValues(G,k,povms)
    lowbounds= []

    best_InitialValues = random_InitialValues(G,k)
    best_sw = SW(best_InitialValues,G)

    for rat in ratios 
        # best_sw,best_InitialValues contiennent les records pour l'itération précédente

        print("\n\nNouveau ratio : v0/(v0+v1) = ",rat,"\n\n")
        G = change_game(G,rat)

        # meilleure pour le ratio précédent 
        print("\nMeilleure pour le ratio précédent\n")
        InitialValues = copy(best_InitialValues)
        sw = seesaw(InitialValues,G,k,delta,eps,seuil)
        print(sw,"\n")

        best_InitialValues = InitialValues # Mise à jour de best_sw = best_InitialValues, qui contiennent le record pour l'itération actuelle
        best_sw = sw


        # aléatoires 
        print("\nMesures aléatoires\n")
        for i in 1:nb_alea
            InitialValues=random_InitialValues(G,k)
            sw = seesaw(InitialValues,G,k,delta,eps,seuil)
            print(sw,"\n")
            if sw>best_sw
                best_InitialValues = InitialValues
                best_sw = sw
                print("Amélioration : nouveau sw est ",best_sw,"\n")
            end
        end


        # given_InitialValues(G,k,povms)
        print("\nMesures prescrites\n")
        to_do = copy.(AllInitialValues)
        for InitialValues in to_do 
            sw = seesaw(InitialValues,G,k,delta,eps,seuil)
            print(sw,"\n")
            if sw>best_sw
                best_InitialValues = InitialValues
                best_sw = sw 
                print("Amélioration : nouveau sw est ",best_sw,"\n")
            end
        end

        push!(lowbounds,best_sw)
    end

    plot(ratios,lowbounds,title="Minorant see-saw",xlabel="v0/(v0+v1)")
end

# Pour la convergence, seuil=30 est trop peu (essayer seuil=200)

#= pour plot seesaw 
G=NC_00_cycle(3,0);k=2;eps=1e-5;seuil=200;
plot_seesaw(G,k,povms1,eps,seuil,30,range(0,0.8,length=20))
=#

