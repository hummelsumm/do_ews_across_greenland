using Pkg
Pkg.activate(".")

include("../data_and_functions.jl")

#1-4
begin
    n_surs = 10_000
   
    which_params = parse(Int64,ARGS[1])
    global param_counter = 0  

    # sranges_10 #55 elements
    # sranges_5 #55 elements
    # sranges_20 #10 elements
    lls = 1:5:length(sranges_20) #2 elements
    uls = 5:5:length(sranges_20) #2 elements
    for (l,u) in zip(lls,uls)
        sr = sranges_20[l:u]
        global param_counter +=1
        if param_counter == which_params
            save_sca_slopes_and_p(n_surs, ice_cores_20, sr, mother = "PAUL", param = -1, 
                    sur_method = "TFTS", dat_type = "orig", cold_type = "gs",
                    filter_indicator = false, smooth_w = false,  
                    exclude_coi = true, normalise = false, only_full = true, folder = "20y")
        end

        global param_counter +=1
        if param_counter == which_params
            save_hurst_slopes_and_p(n_surs, ice_cores_20, sr, mother = "PAUL", param = -1, 
                    sur_method = "TFTS", dat_type = "orig", cold_type = "gs",
                    filter_indicator = false, smooth_w = false,  
                    exclude_coi = true, normalise = false, only_full = true, folder = "20y")
        end
    end
end


