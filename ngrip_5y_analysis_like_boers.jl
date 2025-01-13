# applying Niklas' method to GS only (using a corrected and py2-to-jl-translated version of the code used to obtain the results in 
# Boers, Niklas. “Early-Warning Signals for Dansgaard-Oeschger Events in a High-Resolution Ice Core Record.” Nature Communications 9, no. 1 (July 2, 2018): 2556. https://doi.org/10.1038/s41467-018-04881-7.
using Pkg
Pkg.activate(".")

using Plots
using LinearAlgebra
using LaTeXStrings

include("data_and_functions.jl")

plotlyjs()


function runstd(x, w; only_full = false)
    n = length(x)
    xs = zeros(Union{Missing, Float64},n)
    half_w = div(w, 2)
    #@show half_w
    
    for i in 1:half_w
        # @show i
        # @show x[1:i + half_w]
        if only_full
            xs[i] = missing
        else
            xs[i] = std(x[1:min(i + half_w, n)])
        end
    end
    for i in (n - half_w + 1):n
        if only_full
            xs[i] = missing
        else
            xs[i] = std(x[max(1,i - half_w):end])
        end
    end
    for i in (half_w + 1):(n - half_w)
        xs[i] = std(x[(i - half_w):min(i + half_w, n)])
    end
    return xs
end


###OBS: N uses std and not VAR!!!

function run_niklas_analysis(dat, cold_indices, time; mother = "PAUL", ew_type = ["CSD", "WAVE"][1], s1 = 10, s2 = 50,
                            dt = 5, only_full = false, w = 200, vac = 1, 
                            filter_indicator = true, smooth_w = true, 
                            n_surrs = 200, detrend = ["GLSAR", "DIRECT", "NONE"][1], sur_substract_mean = true,
                            method = ["FOURIER", "TFTS"][1], p_tail = [:both, :right][1], 
                            plotting = true, saving = false, folder = "NGRIP5", ice_name = "N_lowpass",
                            save_entire=false)
    n = length(dat)
    scavf = Int(w/dt)
    samp = Int(dt)

    nodo = length(cold_indices)


    ew_type = uppercase(ew_type)
    
    if ew_type == "CSD"
        ind_names = ["std", "ac"]
        hf_var =  runstd(dat, scavf, only_full=only_full) #std
        hurst = runac(dat, scavf, k=1, v= vac, only_full = only_full) #ac

        #glsar not implemented in jl -> values taken from runs of Niklas Boers' python (2) code
        slopes_glsar_n = reverse([0.0008022767846264092,0.0011695365075209481,0.00022117216034612687,-1.0465302859653137e-05,0.0004409870140254746,-0.00035389913076197966,0.0006044250120980157,-0.00035391568271526025,0.00024788629541157444,0.00076486321393256,0.000702973590445941,0.00017978294475871453,0.00016039834192000104,-0.00014022199444269128,-3.3668688373372495e-05,4.171618965366219e-05,0.0001682085091579384])
        slopes_glsar_n_h = reverse([0.0010087018372028503,0.0004541064623133223,0.0015139924520043253,-0.00016988855564256612,0.0002764944931529208,4.461365113879774e-05,0.00135654136824095,-0.0018099203710023298,0.0005839370020318598,0.00030833932571630834,-8.377983128589686e-05,0.0003336270530099929,-1.7395533978551404e-07,0.00020161860801111566,0.0001122394988878023,0.00023511444131098916,-0.00010712702100554362])

        if saving
            saveto_hf_var = "new_surrogate_files/$(folder)/var/NIKLAS_w_$(w)_$(ice_name)_FiltInd_$(filter_indicator)_onlyfull_$(only_full)_detrend_$(detrend)_DemeanSur_$(sur_substract_mean)_$(n_surrs)_$(method).jld2"
            saveto_hurst = "new_surrogate_files/$(folder)/ac/NIKLAS_w_$(w)_$(ice_name)_FiltInd_$(filter_indicator)_onlyfull_$(only_full)_detrend_$(detrend)_DemeanSur_$(sur_substract_mean)_$(n_surrs)_$(method).jld2"
            save_names = ["Var", "AC"]
            ss1 = 0
            ss2 = w
            println("calculating",saveto_hf_var)
            println("calculating", saveto_hurst)
        end

    elseif ew_type == "WAVE"
        ind_names = ["w", "H"]
        pad = 1  # pad the time series with zeroes (recommended)
        # dj = 0.25  # this will do 4 sub-octaves per octave
        dj = 0.1  
        s0 = 2 * dt  # this says start at a scale of 6 months
        j1 = Int(11 / dj)  # this says do 7 powers-of-two with dj sub-octaves each
        mother = uppercase(mother)

        if mother == "MORLET"
            Cdelta = 0.776  # this is for the MORLET wavelet
        elseif mother == "PAUL"
            #@show mother
            Cdelta = 1.132  # this is for the PAUL wavelet
        end

        wave, period, scale, coi = wavelet(dat,dt,pad=pad,dj=dj,s0=s0,j1=j1,mother=mother, param=-1)
        pow = abs.(wave).^2

        gr = scale .>=s1 
        sm = scale .<=s2 
        avg = gr.*sm
        scale_avg = scale * ones(n)'
        scale_avg = pow ./ scale_avg  # [Eqn(24) in Torrence and Campo, 1998] 
        if smooth_w
            scav = zeros(size(scale_avg[avg,:]))
            for i in 1:size(scav)[1]
                scav[i,:] = runmean(scale_avg[avg, :][i,:], scavf, only_full = only_full)
            end
        else
            scav = scale_avg[avg,:]
        end 

        pa = [Polynomials.fit(log.(scale[avg]), log.(scav)[:,i], 1)[1] for i = 1:size(log.(scav))[2]]
        hurst = (pa .+1)./2

        scale_avg = dj * dt / Cdelta * sum(scale_avg[avg, :], dims = 1)[:]  # [Eqn(24)]

        if smooth_w == true
            scale_avg = runmean(scale_avg, scavf, only_full = only_full)
        end
        hf_var = scale_avg
        
        #this one for N's unfiltered (normalised) data
        slopes_glsar_n = reverse([0.0006925771242266657, 0.0007279313213808535,0.0002132725450006559,-3.323563256125345e-05,0.0002638017646690965,-0.00020863916203724857,0.0004050129395217431,-0.00045816732481727816,0.00021447775821652643,0.0006901404446660931,0.0007619668472772513,0.0002818774420150065,7.582922815768085e-05,-0.00018343748325594958,-4.583862507070955e-05,6.570747523197839e-05,0.00020655194839774133])
        slopes_glsar_n_h = reverse([0.002463966894775794,0.003676498953079716,0.0024449065470704637,-0.0005341508827329744,0.00037253616210878064,-0.0007427402574957651,0.0019248055387631703,-0.002273552344244769,0.0012257557461855377,-0.0012318763313910268,-0.00014815817825908643,0.0006853278436977103,-0.00015062562631962286,0.001534022614636041,0.0008798766096766266,0.000125997915326764,-0.0004341193829324191])
        #this one for N's filtered data:
        #slopes_glsar_n = reverse([0.0006007364657524586,0.0006999215575096019,0.00023563422564215093,-3.300887004856792e-05,0.0002576178104488837,-0.0001946444775865609,0.0005423247174559536,-0.0004458848041202882,0.0002609598360193786,0.0006961190935120936,0.0007983879851254848,0.0002677523880112797,8.936302218275072e-05,-0.00021079224512353858,-4.964670367802907e-05,6.0490315800441826e-05,0.00020655116894138033])
        #slopes_glsar_n_h = reverse([0.0018459833977997322, 0.0034713043334136245,0.0032811205164832605,-0.0005184310596709715,0.0004957625725742515,-0.0008808999893788312,0.0015191414769403417,-0.00227939590068507,0.001218540320680838,-0.0010063940700524264,-0.00011124162629498446,0.0004760500816885633,-9.90655743303126e-05,0.0015322341011990927,0.000876635750366599,0.00013948170804936214,-0.0004514390399015883])
                        #reverse([0.0018430861653613182,0.00340523692064618,0.003262880361225577,-0.0005178011368494553,0.0004920650320694975,-0.0008925118312985695,0.0015035607803211525,-0.0022720454495644324,0.0012179556398722016,-0.001002592286119909,-0.00010887485431048177,0.0004827191015710159,-9.360701260236224e-05,0.0015385454874835805,0.0008795906154912298,0.00013939051091372342,-0.0004571740718922533])
        
        
        if saving
            saveto_hf_var = "new_surrogate_files/$(folder)/sca/NIKLAS_s1_$(s1)_s2_$(s2)_$(ice_name)_$(mother)_FiltInd_$(filter_indicator)_onlyfull_$(only_full)_detrend_$(detrend)_DemeanSur_$(sur_substract_mean)_smoothw_$(smooth_w)_$(n_surrs)_$(method).jld2"
            saveto_hurst = "new_surrogate_files/$(folder)/hurst/NIKLAS_s1_$(s1)_s2_$(s2)_$(ice_name)_$(mother)_FiltInd_$(filter_indicator)_onlyfull_$(only_full)_detrend_$(detrend)_DemeanSur_$(sur_substract_mean)_smoothw_$(smooth_w)_$(n_surrs)_$(method).jld2"
            save_names = ["sca", "hurst"]
            ss1 = s1
            ss2 = s2
            println("calculating",saveto_hf_var)
            println("calculating", saveto_hurst)
        end
    end

    if filter_indicator == true
        hf_var_hf = copy(hf_var)
        hf_var2 = convert.(Float64, hf_var)
        hf_var = cheby_lowpass_filter(hf_var2, .95 * 1. /800, 1. / samp, 8, .05)

        hurst_hf = copy(hurst)
        hurst2 = convert.(Float64, hurst)
        hurst = cheby_lowpass_filter(hurst2, .95 * 1. / 800, 1. / samp, 8, .05)
    else
        hf_var = convert.(Float64, hf_var)
        hurst = convert.(Float64, hurst)
    end

    p= zeros(Float64, nodo)
    kt = Array{Union{Nothing,Float64}}(nothing,nodo) #zeros(nodo)
    ph =  zeros(Float64, nodo)
    kth = Array{Union{Nothing,Float64}}(nothing,nodo) #zeros(nodo)

    pone = Array{Union{Nothing,Float64}}(nothing,nodo) #zeros(Float64, nodo)
    ptwo = Array{Union{Nothing,Float64}}(nothing,nodo) #zeros(Float64, nodo)
    pone_h =Array{Union{Nothing,Float64}}(nothing,nodo) #zeros(Float64, nodo)
    ptwo_h =Array{Union{Nothing,Float64}}(nothing,nodo) #zeros(Float64, nodo)

    surr_slopes = Array{Union{Nothing,Float64}}(nothing,n_surrs, nodo)
    surr_slopes_h = Array{Union{Nothing,Float64}}(nothing,n_surrs, nodo)

    detrend = uppercase(detrend)
    method = uppercase(method)

    function q2(x)
        # NB would subtract linear trend here ...
        sur_slope = Polynomials.fit(1:length(x), x, 1)[1]
        return sur_slope
    end

    for (i,cold) in enumerate(cold_indices)
        ww = hf_var[cold]
        ww_h = hurst[cold] #convert.(Float64, hurst[cold_idx])
        tw = 1:length(ww)

        kt[i] = Polynomials.fit(tw, ww, 1)[1]
        kth[i] = Polynomials.fit(tw,ww_h,1)[1]

        if sur_substract_mean == true
            x = ww.-mean(ww) 
            x_h = ww_h .- mean(ww_h)
        else
            x = ww
            x_h = ww_h
        end 

        if detrend == "GLSAR"
            slope = slopes_glsar_n[i]
            slope_h = slopes_glsar_n_h[i]
        elseif detrend == "DIRECT"
            slope = Polynomials.fit(1:length(x), x, 1)[1]
            slope_h = Polynomials.fit(1:length(x_h), x_h, 1)[1]
        elseif detrend == "NONE"
            slope = zeros(size(x))
            slope_h = zeros(size(x_h))
        end

        x = x .- (slope.*collect(0:1:length(x)-1))
        x_h = x_h .- (slope_h.*collect(0:1:length(x_h)-1))

        if method == "FOURIER"
            test = TimeseriesSurrogates.SurrogateTest(q2,x,RandomFourier(true), n=n_surrs)
            pp = pvalue(test, tail = p_tail)
            
            pone[i] = count(v -> v ≥ kt[i], test.vals)/n_surrs
            pr = count(v -> v ≥ kt[i], test.vals)
            pl = count(v -> v ≤ kt[i], test.vals)
            ptwo[i] =2min(pr, pl)/n_surrs
            if p_tail == :right
                p[i] = pone[i]
            elseif  p_tail == :both 
                p[i] = ptwo[i]
            end
            
            test_h = TimeseriesSurrogates.SurrogateTest(q2,x_h,RandomFourier(true), n=n_surrs)
            pph = pvalue(test_h, tail = p_tail)

            pone_h[i] =count(v -> v ≥ kth[i], test_h.vals)/n_surrs
            pr_h = count(v -> v ≥ kth[i], test_h.vals)
            pl_h = count(v -> v ≤ kth[i], test_h.vals)
            ptwo_h[i] =2min(pr_h, pl_h)/n_surrs
            if p_tail == :right
                ph[i] = pone_h[i]
            elseif  p_tail == :both 
                ph[i] = ptwo_h[i]
            end

            surr_slopes[:,i] = test.vals
            surr_slopes_h[:,i] = test_h.vals

        elseif method == "TFTS"
            test = TimeseriesSurrogates.SurrogateTest(q2,ww,TFTS(0.05), n=n_surrs) # or on ww and ww_h???
            test_h = TimeseriesSurrogates.SurrogateTest(q2,ww_h,TFTS(0.05), n=n_surrs)
            p[i] = TimeseriesSurrogates.pvalue(test, tail=p_tail)
            ph[i] = TimeseriesSurrogates.pvalue(test_h, tail=p_tail)

            pone[i] = count(v -> v ≥ kt[i], test.vals)/n_surrs
            pr = count(v -> v ≥ kt[i], test.vals)
            pl = count(v -> v ≤ kt[i], test.vals)
            ptwo[i] =2min(pr, pl)/n_surrs
            

            pone_h[i] =count(v -> v ≥ kth[i], test_h.vals)/n_surrs
            pr_h = count(v -> v ≥ kth[i], test_h.vals)
            pl_h = count(v -> v ≤ kth[i], test_h.vals)
            ptwo_h[i] =2min(pr_h, pl_h)/n_surrs

            surr_slopes[:,i] = test.vals
            surr_slopes_h[:,i] = test_h.vals
            
        end
    end

    if saving
        if !isfile(saveto_hf_var)
            obj_hf_var = indicator_and_significance("NIKLAS_"*ice_name, [time[cold] for cold in cold_indices], 
                [hf_var[cold] for cold in cold_indices], kt, surr_slopes, save_names[1], ss1, ss2, pone, ptwo, dt)
            save(saveto_hf_var, "slopes", obj_hf_var)
        else
            println(saveto_hf_var, "already exists :)")
        end
        if !isfile(saveto_hurst)
            obj_hurst = indicator_and_significance("NIKLAS_"*ice_name, [time[cold] for cold in cold_indices], 
                [hurst[cold] for cold in cold_indices], kth, surr_slopes_h, save_names[2], ss1, ss2, pone_h, ptwo_h, dt)
            save(saveto_hurst, "slopes", obj_hurst)
        else
            println(saveto_hurst, "already exists :)")
        end

        if save_entire
            save("new_surrogate_files/ice_cores_py/indicators/jl_$(save_names[1])_N_lowpass_py2.jld2", "indicator",hf_var)
            save("new_surrogate_files/ice_cores_py/indicators/jl_$(save_names[2])_N_lowpass_py2.jld2", "indicator",hurst)
            save("new_surrogate_files/ice_cores_py/indicators/jl_$(save_names[1])_hf_N_lowpass_py2.jld2", "indicator",hf_var_hf)
            save("new_surrogate_files/ice_cores_py/indicators/jl_$(save_names[2])_hf_N_lowpass_py2.jld2", "indicator",hurst_hf)
        end
    end
    
    ews_at_all =  findall(x->x>0, kt)
    ews_at_all_h = findall(x->x>0, kth)
    
    ews = intersect(findall(x->x>0, kt), findall(x-> x < 0.05, p))
    ews2 = intersect(findall(x->x>0, kt), findall(x-> x < 0.1, p))
    
    ews_h = intersect(findall(x->x>0, kth), findall(x-> x < 0.05, ph))
    ews2_h = intersect(findall(x->x>0, kth), findall(x-> x < 0.1, ph))
    
    ews_both = intersect(ews, ews_h)
    ews_both2 = intersect(ews2, ews2_h)

    # ew_lengths = [(findlast(c) - findfirst(c)) * samp for c in cold_idx]
    # println("mean ew duration: ", mean(ew_lengths))
    # println("min ew duration: ", minimum(ew_lengths))
    # println("max ew duration: ", maximum(ew_lengths))
    # return

    println("method: $(method), detrended: $(detrend), filter_indicator: $(filter_indicator), n_surrs: $(n_surrs)")
    if ew_type == "WAVE"
        println("mother: $(mother), s1: $(s1), s2: $(s2)")
    end
    #println(ind_names[1])
    #println("# EWS (no sign): ", length(ews_at_all), " out of ", nodo)
    println(ind_names[1], " # sign. EWS (0.05): ", length(ews), " out of ", nodo)
    #@show kt #p
    #println(ind_names[1], "# sign. EWS (0.1): ", length(ews2), " out of ", nodo)
    #println(" ")
    @show reverse(p)
    
    #println(ind_names[2])
    #println("# EWS (no sign): ", length(ews_at_all_h), " out of ", nodo)
    println(ind_names[2], " # sign. EWS (0.05): ", length(ews_h), " out of ", nodo)
    #println(ind_names[2], "# sign. EWS (0.1): ", length(ews2_h), " out of ", nodo)
    #@show kth #ph
    @show reverse(ph)
   # println("both")
    println("# both (0.05): ", length(ews_both), " out of ", nodo)
    #println("# both(0.1): ", length(ews_both2), " out of ", nodo)
    


    if plotting        
        if ew_type == "CSD"
            title = "EWS 100y highpass, method=$method, detrend = $detrend, p=$p_tail, σ: $(length(ews)), α_1: $(length(ews_h)), both: $(length(ews_both))"
            l1 = "σ"
            l2 = "α_1"
        elseif  ew_type == "WAVE"
            title = "EWS s1=$s1, s2=$s2, method=$method, wave = $mother, detrend = $detrend, p=$p_tail, w: $(length(ews)), H: $(length(ews_h)), both: $(length(ews_both))"
            l1 = "w"
            l2 = "H"
        end
        tsp = plot(time, dat, xflip = true, c=:darkred, label = nothing, alpha = 0.6, ylabel = "data")
       
        gs_end = [11705, 14690, 23375, 27790, 28910, 32520, 33735, 35505, 38235, 40165, 41480, 43365, 46860, 49315, 54235, 55815, 58280]
        gs_start = [12900, 23105, 27460, 28510, 32025, 33390, 34725, 36590, 39935, 40790, 42100, 44285, 49105, 51660, 54745, 56440, 58515]

        i1 = plot(time, hf_var, alpha = 0.6, c=:gray, label = nothing, xflip = true, ylabel = l1)
        i2 = plot(time, hurst, alpha = 0.6, c=:gray, label = nothing, xflip = true, ylabel = l2)
        for (i,c) in enumerate(cold_indices)
            co1 = co2 = :darkblue
            lw1 = lw2 = 1
            if i in ews_at_all
                co1=:darkred
            end 
            if i in ews_at_all_h
                co2=:darkred
            end
            if i in ews
                Plots.vspan!(i1, [time[c][1], time[c][end]], alpha = 0.3, c=:darkred, label = nothing)
                lw1 = 3
            else
                Plots.vspan!(i1, [time[c][1], time[c][end]], alpha = 0.3, c=:lightgrey, label = nothing)
            end
            
            if i in ews_h
                Plots.vspan!(i2, [time[c][1], time[c][end]], alpha = 0.3, c=:darkred, label = nothing)
                lw2 = 3
            else
                Plots.vspan!(i2, [time[c][1], time[c][end]], alpha = 0.3, c=:lightgrey, label = nothing)
            end
           
            plot!(tsp, time[c], dat[c], c=:darkblue, label = nothing)
            vline!(tsp, gs_start, label = nothing, ls = :dot, alpha = 0.2)
            vline!(tsp, gs_end, label = nothing, ls = :dash, alpha = 0.2)
            
            plot!(i1, time[c], hf_var[c], c=:black, label = nothing)
            pol1 = Polynomials.fit(time[c], hf_var[c],1)
            plot!(i1, x-> pol1(x), time[c], label = nothing, c=co1, lw = lw1)
            vline!(i1, gs_start, label = nothing, ls = :dot, alpha = 0.2)
            vline!(i1, gs_end, label = nothing, ls = :dash, alpha = 0.2)
            
            plot!(i2, time[c], hurst[c], c=:black, label = nothing)
            pol2 = Polynomials.fit(time[c], hurst[c],1)
            plot!(i2, x-> pol2(x), time[c], label = nothing, c=co2, lw = lw2)
            vline!(i2, gs_start, label = nothing, ls = :dot, alpha = 0.2)
            vline!(i2, gs_end, label = nothing, ls = :dash, alpha = 0.2)
        end

        display(plot(tsp, i1,i2, layout = (3,1), plot_title = title, plot_titlefontsize=10, size = (800,600)))
    end
end


N_core_py = load("new_surrogate_files/ice_cores_py/N_lowpass_py2.jld2")["ice"]
N_core = load("new_surrogate_files/ice_cores/NGRIP5/N_lowpass.jld2")["ice"]

#### what N does:
#ew_type = ["CSD", "WAVE"][1]
#detrend = ["GLSAR", "DIRECT", "NONE"][1]
#method = ["FOURIER", "TFTS"][1]
#p_tail = [:both, :right][2]

for core in [N_core_py, N_core]
    #CSD based:
    run_niklas_analysis(core.δ_normed_filt_100, core.cold_idx_n, core.age, mother = "PAUL",
                ew_type = "CSD", s1 = 10, s2 = 50,
                dt = core.resolution, only_full = false, w = 200, vac = 0, 
                filter_indicator = true, smooth_w = true, 
                n_surrs = 10_000, detrend = "GLSAR", sur_substract_mean = true,
                method = "FOURIER", p_tail = :right, 
                plotting = false, saving = true, folder = "NGRIP5", ice_name = N_core.name,
                save_entire = true)


    #wavelet-based
    run_niklas_analysis(core.δ_normed, core.cold_idx_n, core.age, mother = "PAUL",
                ew_type = "WAVE", s1 = 10, s2 = 50,
                dt = N_core.resolution, only_full = false, w = 200, vac = 0, 
                filter_indicator = true, smooth_w = true, 
                n_surrs = 10_000, detrend = "GLSAR", sur_substract_mean = true,
                method = "FOURIER", p_tail = :right, 
                plotting = false, saving = true, folder = "NGRIP5", ice_name = N_core.name,
                save_entire=true)
end


#Niklas' (corrected) py2 code gives for p<0.05:
#w: 12
#H: 9
#w and H: 8
#σ: 11
#α: 7
#σ and α: 5


