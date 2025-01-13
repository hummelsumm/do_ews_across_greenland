using Pkg
Pkg.activate(".")

include("../data_and_functions.jl")

#1-7
begin
    i = parse(Int64,ARGS[1])
    #ice cores: c_l, n_l, c_no, n_no
    n_surs = 10_000
    py_core = load("new_surrogate_files/ice_cores_py/N_lowpass_py2.jld2")["ice"]
    
    cores = [py_core, py_core, py_core, ice_cores[2], ice_cores[4],ice_cores[3], ice_cores[1]]
    which_periods = ["gs_short","gs_short","gs","gs","gs","gs","gs"]
    filts = [true, false, false, false, false, false, false]
    ofs = [false, false, true, true, true, true, true]

    
    save_var_slopes_and_p(n_surs, [cores[i]],
                    windows = [200], sur_method = "TFTS", 
                    filt_type = "normed_filt", 
                    cold_type = which_periods[i],
                    filter_indicator =filts[i], 
                    only_full = ofs[i], 
                    folder = "NGRIP5")

    save_ac_slopes_and_p(n_surs, [cores[i]],
                    windows = [200], sur_method = "TFTS", 
                    filt_type = ["filt","normed_filt"][2], 
                    cold_type = which_periods[i],
                    filter_indicator = filts[i], 
                    only_full = ofs[i], 
                    folder = "NGRIP5")

end

