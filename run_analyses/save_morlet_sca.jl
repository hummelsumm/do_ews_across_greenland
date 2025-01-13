using Pkg
Pkg.activate(".")

include("../data_and_functions.jl")

#1-6
begin
    which_params = parse(Int64,ARGS[1])
    global param_counter = 0  
    sr = Tuple{Int64, Int64}[(10,50)]
    folders = ["NGRIP5", "NGRIP5", "10y", "10y", "10y", "10y"]
    for (i,ice) in enumerate([ice_cores[1], ice_cores[3], ice_cores_10...])
        global param_counter +=1
        if param_counter == which_params
            save_sca_slopes_and_p(n_surs, [ice], sr, mother = "MORLET", param = -1, 
                sur_method = "TFTS",
                dat_type = "orig", cold_type = "gs",
                filter_indicator = false, smooth_w = false,  
                exclude_coi = true, normalise = false, only_full = true, folder = folders[i])
        end
    end
    
end