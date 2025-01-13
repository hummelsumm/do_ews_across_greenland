using Pkg
Pkg.activate(".")

include("../data_and_functions.jl")

begin
    n_surs = 10_000
    save_var_slopes_and_p(n_surs, ice_cores_10,
                    windows = [200], sur_method = "TFTS", 
                    filt_type = ["filt","normed_filt"][2], 
                    cold_type = ["gs", "gs_short"][1],
                    filter_indicator = false, 
                    only_full = true, 
                    folder = "10y")

    save_ac_slopes_and_p(n_surs, ice_cores_10,
                    windows = [200], sur_method = "TFTS", 
                    filt_type = ["filt","normed_filt"][2], 
                    cold_type = ["gs", "gs_short"][1],
                    filter_indicator = false, 
                    only_full = true, 
                    folder = "10y")
end