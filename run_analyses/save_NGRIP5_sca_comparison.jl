using Pkg
Pkg.activate(".")

include("../data_and_functions.jl")

#1-9
begin
    i = parse(Int64,ARGS[1])
    #ice cores: c_l, n_l, c_no, n_no
    n_surs = 10_000
    sr = Tuple{Int64, Int64}[(10,50)]
    
    py_core = load("new_surrogate_files/ice_cores_py/N_lowpass_py2.jld2")["ice"]
    
    
    cores = [py_core, py_core, py_core, py_core, py_core, ice_cores[2], ice_cores[4],ice_cores[3], ice_cores[1]]
    dtype = ["normed","normed","normed","normed","orig","orig","orig","orig","orig"]
    which_periods = ["gs_short","gs_short","gs_short","gs","gs","gs","gs","gs","gs"]
    filts = [true, false, false, false, false, false, false, false, false]
    nocois = [false, false, false, true, true, true, true, true, true]
    ofs = [false, false, false, false, false, false, false, true, true]
    sws = [true, true, false, false, false, false, false, false, false]

    
    save_sca_slopes_and_p(n_surs, cores[i], sr, mother = "PAUL", param = -1, sur_method = "TFTS",
            dat_type = dtype[i], 
            cold_type = which_periods[i],
            filter_indicator = filts[i], 
            smooth_w = sws[i],  
            exclude_coi = nocois[i], 
            normalise = false, only_full = ofs[i], folder = "NGRIP5")

end

