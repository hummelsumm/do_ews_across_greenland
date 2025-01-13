using Pkg
Pkg.activate(".")

include("../data_and_functions.jl")

#1-11
begin
    n_surs = 10_000
   
    which_params = parse(Int64,ARGS[1])
    global param_counter = 0  

    lls = 1:5:length(sranges_5) #11 elements
    uls = 5:5:length(sranges_5) #11 elements
    for (l,u) in zip(lls,uls)
        sr = sranges_5[l:u]
        global param_counter +=1
        if param_counter == which_params
            save_sca_slopes_and_p(n_surs, [ice_cores[1], ice_cores[3]], sr, mother = "PAUL", param = -1, 
                sur_method = "TFTS", dat_type = "orig", cold_type = "gs",
                filter_indicator = false, smooth_w = false,  
                exclude_coi = true, normalise = false, only_full = true, folder = "NGRIP5")
        end
    end
end
