using Pkg
Pkg.activate(".")

include("../data_and_functions.jl")

begin
    which_params = parse(Int64,ARGS[1])
    global param_counter = 0

    for c in ["lp", "no_lp"]
        for i = 1:2
            global param_counter +=1
            
            if param_counter == which_params
                n_surrs = 2_000
                n_surrs_per_gs = 1_000
                p_thresholds = [0.05,0.1,0.2,0.3]
                smooth_w = false

                if c== "lp"
                    ice_core = ice_cores[1]
                else
                    ice_core = ice_cores[3]
                end
                
                if param_counter%2 == 0
                    sign_number_of_ac_numeric([ice_core],n_surrs, n_surrs_per_gs,  p_thresholds,
                        windows= [200], sur_method = "TFTS", filt_type="normed_filt", cold_type = "gs",
                        filter_indicator = false, only_full= true, folder = "NGRIP5", gs_lengths = gs_lengths, gs_lengths_n=gs_lengths_n)
                else 
                    sign_number_of_var_numeric([ice_core],n_surrs, n_surrs_per_gs,  p_thresholds,
                    windows= [200], sur_method = "TFTS", filt_type="normed_filt", cold_type = "gs",
                    filter_indicator = false, only_full= true, folder = "NGRIP5", gs_lengths = gs_lengths, gs_lengths_n=gs_lengths_n)
                end
            end
        end
    end
end