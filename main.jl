
#= 

Les différents fichiers sont indépendants (on peut donc ouvrir l'un d'eux sans avoir à ouvrir les autres). Le fichier main se contente de tous les inclure.

=#

include("game_gen.jl")
include("npa_upbound_Qcorr.jl")
include("seesaw_lowbound.jl")