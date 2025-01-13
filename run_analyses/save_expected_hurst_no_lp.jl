using Pkg
Pkg.activate(".")

include("../data_and_functions.jl")

#only lowpass=false
begin
    which_params = parse(Int64,ARGS[1])
    #@show which_params

    #n_surrs = 1_060 #overall
    if which_params <=4
        n_surrs = 25
    elseif 5 <= which_params <=12
        n_surrs = 100
    elseif which_params == 13
        n_surrs = 160
    end

    n_surrs_per_gs = 1_000
    p_thresholds = [0.05,0.1,0.2,0.3]

    
    smooth_w = false
    sign_number_of_hurst_numeric([ice_cores[3]],  
                    n_surrs,
                    n_surrs_per_gs,  
                    p_thresholds, 
                    [(10,50)],
                    mother = "PAUL", param = -1, 
                    sur_method= "TFTS", 
                    dat_type = "orig", 
                    cold_type = "gs",
                    filter_indicator = false, smooth_w = false,  exclude_coi = true, normalise = false, only_full = true,
                    folder = "NGRIP5", gs_lengths = gs_lengths, gs_lengths_n = gs_lengths_n, batch_no = which_params)

end




