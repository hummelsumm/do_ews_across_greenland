using Pkg
Pkg.activate(".")

include("../data_and_functions.jl")


begin
    n_surrs = 2_000
    n_surrs_per_gs = 1_000
    p_thresholds = [0.05,0.1,0.2,0.3]
    smooth_w = false
    
    new_sign_number_of_sca_numeric([ice_cores[3]], 
                    n_surrs,
                    n_surrs_per_gs,  
                    p_thresholds, 
                    [(10,50)],
                    mother = "PAUL", param = -1, 
                    sur_method= "TFTS", 
                    dat_type = "orig", 
                    cold_type = "gs",
                    filter_indicator = false, smooth_w = smooth_w,  exclude_coi = true, normalise = false, only_full = true,
                    folder = "NGRIP5", gs_lengths = gs_lengths, gs_lengths_n = gs_lengths_n, batch_no = -42)
end
