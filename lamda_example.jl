# using Pkg
# Pkg.activate(".")

using JLD2, FileIO
using NPZ
using CairoMakie
using Polynomials
using LaTeXStrings
using Statistics


mutable struct ice_core
    name::String
    age::Array{Float64}
    δ::Array{Float64}
    δ_normed::Array{Float64}
    δ_filt_100::Array{Float64}
    δ_normed_filt_100::Array{Float64} # 100y filtered of normed data
    resolution::Int64
    cold_idx::Array
    warm_idx::Array
    cold_idx_n::Array #NB searches for EWS in this interval (and not the entire one!!)
end

struct indicator_and_significance
    name::String
    times::Array{Vector{Float64}}  
    vals::Array{Vector{Union{Missing,Float64}}} 
    slopes::Array{Union{Nothing,Float64}}
    surr_slopes::Array{Union{Nothing, Float64}}
    type::String
    s1::Int32
    s2::Int32
    p_one::Vector{Union{Nothing,Float64}}
    p_two::Vector{Union{Nothing,Float64}}
    resolution::Int32
end


#save ice core as npz file to calculate lambda in python
if !isfile("do_ews_across_greenland_ice_cores/ngrip_5.npz")  
    xx = load("new_surrogate_files/ice_cores/NGRIP5/C_lowpass.jld2")["ice"]

    #save each GS
    d = Dict{String, Any}()
    for (i,c) in enumerate(xx.cold_idx)
        t = xx.age[c]
        x= xx.δ_normed_filt_100[c]
        d["$i"] = [t x]
    end
    npzwrite("do_ews_across_greenland_ice_cores/ngrip_5.npz", d)
end

#convert to the same data objects we use here
saveto_lambda = "do_ews_across_greenland_ice_cores/new_surrogate_files/NGRIP5/lambda_w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2"
if !isfile(saveto_lambda) 
    lambda_file = load("do_ews_across_greenland_ice_cores/lambda_ngrip5_w_40_10000_TFTS.npz")

    l_times = []
    l_vals = []
    l_p_ones = []
    for i in 1:17
        t = lambda_file["times"][i,:]
        x = lambda_file["vals"][i,:]
        new_t = t[isfinite.(t)]
        @show typeof(new_t)
        push!(l_times, new_t)

        new_x = x[isfinite.(t)]
        new_x = convert(Vector{Union{Float64, Missing}}, new_x)
        new_x[findall(isnan, new_x)] .= missing
        push!(l_vals, new_x)

        #mistake in the python code...
        new_p = count(v -> v ≥ lambda_file["slopes"][i], lambda_file["surr_slopes"][i,:])/10_000
        push!(l_p_ones, new_p)
    end

    lambda_obj = indicator_and_significance(
        "NGRIP5",
        l_times,
        l_vals,
        lambda_file["slopes"],
        lambda_file["surr_slopes"],
        "lambda",
        0,
        100,
        l_p_ones,
        [],
        5)
    save(saveto_lambda, "slopes", lambda_obj)
end

#load the right datasets to plot
v = load("do_ews_across_greenland_ice_cores/new_surrogate_files/NGRIP5/var/w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2")["slopes"]
ac = load("do_ews_across_greenland_ice_cores/new_surrogate_files/NGRIP5/ac/w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2")["slopes"]
lambda = load(saveto_lambda)["slopes"]
w = load("do_ews_across_greenland_ice_cores/new_surrogate_files/NGRIP5/sca/s1_10_s2_50_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2")["slopes"]
H =load("do_ews_across_greenland_ice_cores/new_surrogate_files/NGRIP5/hurst/s1_10_s2_50_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2")["slopes"]

# for i in corrs
#     @show i
# end

# eachin
# plot correlations

corrs = zeros(Float64, 17,4)
for i = 1:17
    mask_csd = findall(!ismissing, lambda.vals[i])
    ca = cor(lambda.vals[i][mask_csd], ac.vals[i][mask_csd])
    cv = cor(lambda.vals[i][mask_csd], v.vals[i][mask_csd])
    mask_h = intersect(findall(!ismissing, lambda.vals[i]), findall(!ismissing, H.vals[i]))
    ch = cor(lambda.vals[i][mask_h], H.vals[i][mask_h])
    mask_w = intersect(findall(!ismissing, lambda.vals[i]), findall(!ismissing, w.vals[i]))
    cw = cor(lambda.vals[i][mask_w], w.vals[i][mask_w])
    corrs[i,1] = cv
    corrs[i,2] = ca
    corrs[i,3] = cw
    corrs[i,4] = ch
    #cor(lambda.vals[1][findall(!ismissing, lambda.vals[2])], H.vals[1][findall(!ismissing, lambda.vals[2])])
end

corrs
mean(corrs, dims = 1)
#That properly!!
begin
    f = Figure(size=(1000,300))

    event_labels = ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]
    ylabs = [L"$\mathrm{𝖵}$",L"$\mathrm{\alpha_𝟣}$",  L"$\mathrm{\hat{𝗐}^𝟤_{𝟣𝟢,𝟧𝟢}}$", L"$\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}_{𝟣𝟢,𝟧𝟢}}$"]
    ax = Axis(f[1,1],
        #title = "Significance of the observed number of EWS",
        xlabel = "Event",
        xticks = (1:17, event_labels),
        xminorgridwidth = 1.0, 
        xminorgridcolor = :grey30 ,
        xminorgridvisible = true, xgridvisible = false,
        xticklabelrotation = pi/3,
        #yminorticks = 1.5:1:8.5,
        ylabel = "Indicator",
        yticks = (1:4, ylabs),
       yticklabelsize=18,
        yminorgridwidth = 1, 
        yminorgridcolor = :grey30, 
        yminorgridvisible = true, ygridvisible = false)
    ax.yreversed=true
    ax.xreversed=true


    hm = heatmap!(ax, corrs, colormap = :balance, colorrange=(-1,1))
    translate!(hm, 0,0,-100)

    for i in 1:17
        for j in 1:4
            text!(i,j,text="$(round(corrs[i,j], digits=2))", align = (:center, :center),
                color=:grey45)
        end
    end

    ga = f[1,2] = GridLayout()
    Colorbar(ga[1,1], hm,
        label = "Cor(λ,indicator) prior to event" )
        # halign=:left)
    #display(f)
    #save("do_ews_across_greenland_ice_cores/lambda_corrs_NGRIP5_text2.png", f)
    save("do_ews_across_greenland_ice_cores/figures/lambda_corrs_NGRIP5_text2.pdf", f)
    f
end

#plot:
#Make nice
begin
    plim = 0.05
    GS_onsets = [12_900,23_105,27_460,28_510,32_025,33_390,34_725,36_590,39_935,40_790,42_100,44_285,49_105,51_650,54_745,56_440,58_515]
    GI_onsets = [11_705,14_690,23_375,27_790,28_910,32_520,33_735,35_505,38_235,40_165,41_480,43_365,46_860,49_315,54_235,55_815,58_280]
    
    event_labels = ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]
    
    letters = ["(a)", "(b)", "(c)", "(d)", "(e)"]
    ylabs= [L"Restoring rate $\mathrm{λ}$",  L"Variance $\mathrm{V}$",L"Autocorrelation $\mathrm{\alpha_1}$", L"Wavelet coefficient $\mathrm{\hat{𝗐}^𝟤_{10,50}}$", L"Hurst exponent $\mathrm{\hat{H}^\text{𝗅𝗈𝖼}_{10,50}}$"]
    ylabs_short = [L"$\mathrm{λ}$",  L"$\mathrm{𝖵}$", L"$\mathrm{\alpha_𝟣}$", L"$\mathrm{\hat{𝗐}^𝟤_{𝟣𝟢,𝟧𝟢}}$", L"$\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}_{𝟣𝟢,𝟧𝟢}}$"]
    f = Figure(size=(1000,1000))
    for (i,ind) in enumerate([lambda, v, ac, w, H])#enumerate([var, ac, lambda, w, H])
        gl = f[i,:] = GridLayout()
        ax = Axis(gl[1,1], 
            xlabel = "Age (kyr b2k)", 
            ylabel = ylabs_short[i], 
            ylabelsize = 18,
            #xlabelsize= 20,
            ylabelfont = :bold,
            xticks = (10_000:5_000:60_000, string.(10:5:60)), 
            )
        ax.xreversed = true

        ax_event = Axis(gl[1,1],xaxisposition = :top, xticklabelrotation = pi*0.4, ygridvisible=false)
        ax_event.xticks = (GI_onsets, event_labels)
        ax_event.xreversed = true
        hidespines!(ax_event)
        hideydecorations!(ax_event, grid = false)
        
        n_v=0
        for ev in 1:17
            va = 0.4
            vl = 2.0
            vsc = (:white,1.0)
            vc = :black
            if typeof(ind.slopes[ev]) != Nothing
                if ind.slopes[ev] >0 
                    vc = :red
                    vsc = (:darkred,0.1)
                    if ind.p_one[ev] < plim
                        va = 0.9
                        vl = 3.5
                        vsc = (:darkred,0.8)
                        n_v +=1
                    end      
                else
                    vc = :blue
                    vsc = (:steelblue,0.1)
                end
            end
            vspan!(ax, [GI_onsets[ev]], [GS_onsets[ev]], color = vsc)
            pv = Polynomials.fit(ind.times[ev][findall(!ismissing, ind.vals[ev])], convert.(Float64,ind.vals[ev][findall(!ismissing, ind.vals[ev])]),1)
            lines!(ax,ind.times[ev], ind.vals[ev], color = :black, alpha = va)#, linewidth = vl)
            lines!(ax, ind.times[ev], pv.(ind.times[ev]), color = vc, alpha = va, linewidth = vl)
            vlines!(ax, [GI_onsets[ev]], color = :darkred, alpha = va, linewidth = vl)
        end
        xlims!(ax,nothing,nothing)
        ax.xreversed = true
        xlims!(ax_event, ax.xaxis.attributes.limits[]...)
        ax_event.xreversed = true
        if i>1
            hidespines!(ax,:t)
            hidexdecorations!(ax_event, grid=false)
        end
        if i<5
            hidespines!(ax,:b)
            hidexdecorations!(ax, grid=false)
        end
        
        Label(f[i,1,TopLeft()], letters[i], fontsize = 20,
            font = :bold, padding = (0,20,-10,15),
            halign = :left,valign =:bottom)
        rowgap!(f.layout, 10)
    end
    elems = [
                #[LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(0.5, 0), (0.5, 1)])],
                [LineElement(color= (:black, 0.4), points = Point2f[(0, 0.5), (0.5, 0.5)]), LineElement(color= (:black, 0.9), points = Point2f[(0.5, 0.5), (1, 0.5)])],
                [PolyElement(color = (:darkred, 0.1), strokewidth = 0), LineElement(color= (:red,0.4), linewidth=2.0)],# LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
                [PolyElement(color = (:steelblue,0.1), strokewidth = 0), LineElement(color= (:blue, 0.4), linewidth=2.0)],#,LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
                [PolyElement(color = (:darkred, 0.8), strokewidth = 0), LineElement(color= (:red,0.9), linewidth=3.5)],#,LineElement(color = (:darkred,0.9), linewidth= 3.5,points = Point2f[(1, 0), (1, 1)])]]
                ]
    labels = [#"GI onset", 
            "EWS indicator", "increasing", "decreasing", "significantly increasing (𝗉<0.05)"]
    Legend(f[end+1,1:end], elems, labels, orientation = :horizontal, valign = :bottom, margin = (10,10,10,10)) 
    # save("do_ews_across_greenland_ice_cores/lambda_EWS_comparison_NGRIP5_v3.png", f)
    save("do_ews_across_greenland_ice_cores/figures/lambda_EWS_comparison_NGRIP5_v3.pdf", f)
    f
end


# NO LP!!

begin
    v_n = load("do_ews_across_greenland_ice_cores/new_surrogate_files/NGRIP5/var/w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2")["slopes"]
    ac_n = load("do_ews_across_greenland_ice_cores/new_surrogate_files/NGRIP5/ac/w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2")["slopes"]
    #lambda = load(saveto_lambda)["slopes"]
    w_n = load("do_ews_across_greenland_ice_cores/new_surrogate_files/NGRIP5/sca/s1_10_s2_50_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2")["slopes"]
    H_n =load("do_ews_across_greenland_ice_cores/new_surrogate_files/NGRIP5/hurst/s1_10_s2_50_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2")["slopes"]
    
    plim = 0.05
    GS_onsets = [12_900,23_105,27_460,28_510,32_025,33_390,34_725,36_590,39_935,40_790,42_100,44_285,49_105,51_650,54_745,56_440,58_515]
    GI_onsets = [11_705,14_690,23_375,27_790,28_910,32_520,33_735,35_505,38_235,40_165,41_480,43_365,46_860,49_315,54_235,55_815,58_280]
    
    event_labels = ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]
    
    letters = ["(a)", "(b)", "(c)", "(d)"]
    ylabs= [L"Variance $\mathrm{V}$",L"Autocorrelation $\mathrm{\alpha_1}$", L"Wavelet coefficient $\mathrm{\hat{𝗐}^𝟤_{10,50}}$", L"Hurst exponent $\mathrm{\hat{H}^\text{𝗅𝗈𝖼}_{10,50}}$"]
    ylabs_short = [L"$\mathrm{𝖵}$", L"$\mathrm{\alpha_𝟣}$", L"$\mathrm{\hat{𝗐}^𝟤_{𝟣𝟢,𝟧𝟢}}$", L"$\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}_{𝟣𝟢,𝟧𝟢}}$"]
    f = Figure(size=(1000,800)) #1000
    for (i,ind) in enumerate([v_n, ac_n, w_n, H_n])#enumerate([var, ac, lambda, w, H])
        gl = f[i,:] = GridLayout()
        ax = Axis(gl[1,1], 
            xlabel = "Age (kyr b2k)", 
            ylabel = ylabs_short[i], 
            ylabelsize = 18,
            #xlabelsize= 20,
            ylabelfont = :bold,
            xticks = (10_000:5_000:60_000, string.(10:5:60)), 
            )
        ax.xreversed = true

        ax_event = Axis(gl[1,1],xaxisposition = :top, xticklabelrotation = pi*0.4, ygridvisible=false)
        ax_event.xticks = (GI_onsets, event_labels)
        ax_event.xreversed = true
        hidespines!(ax_event)
        hideydecorations!(ax_event, grid = false)
        
        n_v=0
        for ev in 1:17
            va = 0.4
            vl = 2.0
            vsc = (:white,1.0)
            vc = :black
            if typeof(ind.slopes[ev]) != Nothing
                if ind.slopes[ev] >0 
                    vc = :red
                    vsc = (:darkred,0.1)
                    if ind.p_one[ev] < plim
                        va = 0.9
                        vl = 3.5
                        vsc = (:darkred,0.8)
                        n_v +=1
                    end      
                else
                    vc = :blue
                    vsc = (:steelblue,0.1)
                end
            end
            vspan!(ax, [GI_onsets[ev]], [GS_onsets[ev]], color = vsc)
            pv = Polynomials.fit(ind.times[ev][findall(!ismissing, ind.vals[ev])], convert.(Float64,ind.vals[ev][findall(!ismissing, ind.vals[ev])]),1)
            lines!(ax,ind.times[ev], ind.vals[ev], color = :black, alpha = va)#, linewidth = vl)
            lines!(ax, ind.times[ev], pv.(ind.times[ev]), color = vc, alpha = va, linewidth = vl)
            vlines!(ax, [GI_onsets[ev]], color = :darkred, alpha = va, linewidth = vl)
        end
        xlims!(ax,nothing,nothing)
        ax.xreversed = true
        xlims!(ax_event, ax.xaxis.attributes.limits[]...)
        ax_event.xreversed = true
        if i>1
            hidespines!(ax,:t)
            hidexdecorations!(ax_event, grid=false)
        end
        if i<5
            hidespines!(ax,:b)
            hidexdecorations!(ax, grid=false)
        end
        
        Label(f[i,1,TopLeft()], letters[i], fontsize = 20,
            font = :bold, padding = (0,20,-10,15),
            halign = :left,valign =:bottom)
        rowgap!(f.layout, 10)
    end
    elems = [
                #[LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(0.5, 0), (0.5, 1)])],
                [LineElement(color= (:black, 0.4), points = Point2f[(0, 0.5), (0.5, 0.5)]), LineElement(color= (:black, 0.9), points = Point2f[(0.5, 0.5), (1, 0.5)])],
                [PolyElement(color = (:darkred, 0.1), strokewidth = 0), LineElement(color= (:red,0.4), linewidth=2.0)],# LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
                [PolyElement(color = (:steelblue,0.1), strokewidth = 0), LineElement(color= (:blue, 0.4), linewidth=2.0)],#,LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
                [PolyElement(color = (:darkred, 0.8), strokewidth = 0), LineElement(color= (:red,0.9), linewidth=3.5)],#,LineElement(color = (:darkred,0.9), linewidth= 3.5,points = Point2f[(1, 0), (1, 1)])]]
                ]
    labels = [#"GI onset", 
            "EWS indicator", "increasing", "decreasing", "significantly increasing (𝗉<0.05)"]
    Legend(f[end+1,1:end], elems, labels, orientation = :horizontal, valign = :bottom, margin = (10,10,10,10)) 
    save("do_ews_across_greenland_ice_cores/figures/EWS_comparison_NGRIP5_no_lp.pdf", f)
    f
end