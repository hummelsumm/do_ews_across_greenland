using Pkg
Pkg.activate(".")

using CairoMakie
using LaTeXStrings
using JLD2, FileIO
using Distributions
using StatsBase, Statistics
using Polynomials
using Interpolations
using NPZ
using Roots
using DelimitedFiles, XLSX
using TimeseriesSurrogates
using GeoMakie
using NaturalEarth

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

mutable struct distribution_significant_increases
    name::String
    num_inc_one::Array
    num_inc_two::Array 
    n::Int32
    n_GS::Int32
    pvals::Any 
    type::String
    s1::Int32
    s2::Int32
    resolution::Int32
end


GS_onsets = [12_900,23_105,27_460,28_510,32_025,33_390,34_725,36_590,39_935,40_790,42_100,44_285,49_105,51_650,54_745,56_440,58_515]
GI_onsets = [11_705,14_690,23_375,27_790,28_910,32_520,33_735,35_505,38_235,40_165,41_480,43_365,46_860,49_315,54_235,55_815,58_280]

missing_idx_GISP2 = [1126, 1127, 1128, 1129, 1130, 1197, 1265, 1364, 1415, 1416, 1526, 1618, 1619, 1620, 1621, 1784, 1785, 1892, 1893, 1894, 1935, 2419, 2472, 2482]
journal = "esd" #nothing #"esd","cd"
plim = 0.05
lowpass = true
smoothw = false
saving = true
showing = true

letters = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"]

letters_esd = "(".*lowercase.(letters).*")"
letters_cd = lowercase.(letters)

if journal == "esd"
    letters = letters_esd
elseif journal == "cd"
    letters = letters_cd
end

event_labels = ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]

ylabs2 = [["𝖭𝖦𝖱𝖨𝖯"],["𝖭𝖤𝖤𝖬","𝖭𝖦𝖱𝖨𝖯"],["𝖦𝖨𝖲𝖯𝟤","𝖦𝖱𝖨𝖯","𝖭𝖦𝖱𝖨𝖯"]]
ylabs2b = [["𝗡𝗚𝗥𝗜𝗣"],["𝗡𝗘𝗘𝗠","𝗡𝗚𝗥𝗜𝗣"],["𝗚𝗜𝗦𝗣𝟮","𝗚𝗥𝗜𝗣","𝗡𝗚𝗥𝗜𝗣"]]


#get sans serif version of integers for text in plots
function sstring(i::Int)
    if i==0
        return "𝟢"
    elseif i<=9
        return ["𝟣","𝟤","𝟥","𝟦","𝟧","𝟨","𝟩","𝟪","𝟫"][i]
    else
        st = string(i)
        ret = ""
        for d in st
            ret *= ["𝟣","𝟤","𝟥","𝟦","𝟧","𝟨","𝟩","𝟪","𝟫"][parse(Int64, d)]
        end
        return ret
    end
end

#get sans serif version of strings for text in plots
function sstring(s::String)
    Letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    letters = lowercase(Letters)
    numbers = "0123456789"
    ks = [string(i) for i in Letters*letters*numbers]

    sans_Letters = "𝖠𝖡𝖢𝖣𝖤𝖥𝖦𝖧𝖨𝖩𝖪𝖫𝖬𝖭𝖮𝖯𝖰𝖱𝖲𝖳𝖴𝖵𝖶𝖷𝖸𝖹"
    sans_letters = "𝖺𝖻𝖼𝖽𝖾𝖿𝗀𝗁𝗂𝗃𝗄𝗅𝗆𝗇𝗈𝗉𝗊𝗋𝗌𝗍𝗎𝗏𝗐𝗑𝗒𝗓"
    sans_numbers = "𝟢𝟣𝟤𝟥𝟦𝟧𝟨𝟩𝟪𝟫"
    vs = [string(i) for i in sans_Letters*sans_letters*sans_numbers]
    
    trad_dict = Dict(ks .=> vs)
    #@show trad_dict
    new_str = ""
    for i in s
        #@show i
        try
            new_str*= trad_dict[string(i)]
        catch
            new_str*=string(i)
        end

    end
    return new_str
end

my_cdf(d,n) = cdf(d,n-1)
my_ecdf(arr, n) = ecdf(arr)(n-1)

function significance_threshold(dist, p; test_vals = [1,2,3,4,5,6,7])
    if typeof(dist) == distribution_significant_increases #numeric
        pv = round(1-p , digits=2)
        if pv ∉ dist.pvals
            error("Please choose p s.t. 1-p ∈ $(dist.pvals)")
        end
        nums = dist.num_inc_one[dist.pvals .== pv,:][:]
        threshold = findfirst(>(p),[my_ecdf(nums,i) for i in test_vals])

    elseif typeof(dist) == Binomial{Float64} #analytical
        threshold = test_vals[findfirst(>(p),my_cdf.(dist,test_vals))]
    else 
        error("dist must be of type distribution_significant_increases or Binomial{Float64}")
    end
    threshold
end


#Fig. 10, S14
function get_numbers_base_cases(;lowpass= lowpass, showing = true, saving = false, saveto = "paper/result_overview_n_inc_lowpass_$(lowpass)_smoothw_false_p_0.05.pdf")
    function get_n(x,plim)
        n=0
        for ev in 1:17
            if typeof(x.slopes[ev]) != Nothing
                if x.slopes[ev] >0  && x.p_one[ev] < plim
                    n+=1
                end
            end
        end
        return n
    end

    b = Distributions.Binomial(17,0.05)
    thr95 = significance_threshold(b, 0.95)
    thr90 = significance_threshold(b, 0.9)

    cases = [L"\textbf{\mathrm{\sigma^𝟤}}",L"\mathrm{\alpha_𝟣}",L"$\textbf{\mathrm{\hat{𝗐}^𝟤}}$, 𝟣𝟢-𝟧𝟢 𝗒, 𝖯𝖺𝗎𝗅", L"$\mathrm{\hat{𝖧}}$, 𝟣𝟢-𝟧𝟢 𝗒, 𝖯𝖺𝗎𝗅", L"$\mathrm{\hat{𝗐}^𝟤}$, 𝟣𝟢-𝟧𝟢 𝗒, 𝖬𝗈𝗋𝗅𝖾𝗍", L"$\mathrm{\hat{𝖧}}$, 𝟣𝟢-𝟧𝟢 𝗒, 𝖬𝗈𝗋𝗅𝖾𝗍", L"$\mathrm{\hat{𝗐}^𝟤}$, 𝟤𝟢-𝟨𝟢 𝗒, 𝖯𝖺𝗎𝗅", L"$\mathrm{\hat{𝖧}_{𝟤𝟢-𝟨𝟢 𝗒, 𝖯𝖺𝗎𝗅}}$"]

    ngrip_irreg_paths = ["new_surrogate_files/NGRIP_irreg/var_filt_1000_tfts.jld2",
                        "new_surrogate_files/NGRIP_irreg/ac_filt_1000_tfts.jld2",
                        missing,
                        missing,
                        "new_surrogate_files/NGRIP_irreg/sca_10_50_1000_tfts.jld2",
                        "new_surrogate_files/NGRIP_irreg/hurst_10_50_1000_tfts.jld2",
                        missing,
                        missing]
    
    ngrip5_paths = ["new_surrogate_files/NGRIP5/var/w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/NGRIP5/ac/w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/NGRIP5/sca/s1_10_s2_50_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/NGRIP5/hurst/s1_10_s2_50_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/NGRIP5/sca/s1_10_s2_50_C_lowpass_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/NGRIP5/hurst/s1_10_s2_50_C_lowpass_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/NGRIP5/sca/s1_20_s2_60_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/NGRIP5/hurst/s1_20_s2_60_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]
    
    ngrip10_paths = ["new_surrogate_files/10y/var/w_200_normed_filt_NGRIP_lp_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/10y/ac/w_200_normed_filt_NGRIP_lp_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/10y/sca/s1_10_s2_50_NGRIP_lp_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/10y/hurst/s1_10_s2_50_NGRIP_lp_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/10y/sca/s1_10_s2_50_NGRIP_lp_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/10y/hurst/s1_10_s2_50_NGRIP_lp_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/10y/sca/s1_20_s2_60_NGRIP_lp_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/10y/hurst/s1_20_s2_60_NGRIP_lp_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]
    
    neem_paths = ["new_surrogate_files/10y/var/w_200_normed_filt_NEEM_lp_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/10y/ac/w_200_normed_filt_NEEM_lp_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/10y/sca/s1_10_s2_50_NEEM_lp_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/10y/hurst/s1_10_s2_50_NEEM_lp_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/10y/sca/s1_10_s2_50_NEEM_lp_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/10y/hurst/s1_10_s2_50_NEEM_lp_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/10y/sca/s1_20_s2_60_NEEM_lp_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/10y/hurst/s1_20_s2_60_NEEM_lp_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]
    
    ngrip20_paths = ["new_surrogate_files/20y/var/w_200_normed_filt_NGRIP_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/20y/ac/w_200_normed_filt_NGRIP_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    missing,
                    missing,
                    missing,
                    missing,
                    "new_surrogate_files/20y/sca/s1_20_s2_60_NGRIP_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/20y/hurst/s1_20_s2_60_NGRIP_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]
    
    grip_paths = ["new_surrogate_files/20y/var/w_200_normed_filt_GRIP_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/20y/ac/w_200_normed_filt_GRIP_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    missing,
                    missing,
                    missing,
                    missing,
                    "new_surrogate_files/20y/sca/s1_20_s2_60_GRIP_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/20y/hurst/s1_20_s2_60_GRIP_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]

    
    gisp2_paths = ["new_surrogate_files/20y/var/w_200_normed_filt_GISP2_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/20y/ac/w_200_normed_filt_GISP2_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    missing,
                    missing,
                    missing,
                    missing,
                    "new_surrogate_files/20y/sca/s1_20_s2_60_GISP2_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/20y/hurst/s1_20_s2_60_GISP2_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]

    if lowpass == false
        ngrip5_paths = ["new_surrogate_files/NGRIP5/var/w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                        "new_surrogate_files/NGRIP5/ac/w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                        "new_surrogate_files/NGRIP5/sca/s1_10_s2_50_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                        "new_surrogate_files/NGRIP5/hurst/s1_10_s2_50_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                        "new_surrogate_files/NGRIP5/sca/s1_10_s2_50_C_no_lowpass_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                        "new_surrogate_files/NGRIP5/hurst/s1_10_s2_50_C_no_lowpass_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                        "new_surrogate_files/NGRIP5/sca/s1_20_s2_60_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                        "new_surrogate_files/NGRIP5/hurst/s1_20_s2_60_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]

        ngrip10_paths = ["new_surrogate_files/10y/var/w_200_normed_filt_NGRIP_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                        "new_surrogate_files/10y/ac/w_200_normed_filt_NGRIP_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                        "new_surrogate_files/10y/sca/s1_10_s2_50_NGRIP_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                        "new_surrogate_files/10y/hurst/s1_10_s2_50_NGRIP_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                        "new_surrogate_files/10y/sca/s1_10_s2_50_NGRIP_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                        "new_surrogate_files/10y/hurst/s1_10_s2_50_NGRIP_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                        #"new_surrogate_files/10y/hurst/s1_10_s2_50_NGRIP_lp_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS",
                        "new_surrogate_files/10y/sca/s1_20_s2_60_NGRIP_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                        "new_surrogate_files/10y/hurst/s1_20_s2_60_NGRIP_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]
        
        neem_paths = ["new_surrogate_files/10y/var/w_200_normed_filt_NEEM_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                        "new_surrogate_files/10y/ac/w_200_normed_filt_NEEM_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                        "new_surrogate_files/10y/sca/s1_10_s2_50_NEEM_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                        "new_surrogate_files/10y/hurst/s1_10_s2_50_NEEM_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                        "new_surrogate_files/10y/sca/s1_10_s2_50_NEEM_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                        "new_surrogate_files/10y/hurst/s1_10_s2_50_NEEM_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                        "new_surrogate_files/10y/sca/s1_20_s2_60_NEEM_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                        "new_surrogate_files/10y/hurst/s1_20_s2_60_NEEM_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]
    end
    
    issig = Array{Union{Missing,Int}}(missing,7,length(cases)) 
    justn = Array{Union{Missing,Int}}(missing,7,length(cases))  
    for (i,core) in enumerate([ngrip_irreg_paths, ngrip5_paths, ngrip10_paths, neem_paths, ngrip20_paths, grip_paths, gisp2_paths])
        for (j,path) in enumerate(core)
            if !ismissing(path)
                ind = load(path)["slopes"]
                n = get_n(ind,0.05)
                justn[i,j]=n
                if n>=thr95
                    #sign at 95%
                    issig[i,j] = 95
                elseif n>=thr90
                    #sign at 90
                    issig[i,j] = 90
                else
                    #not significant
                    issig[i,j] = 0
                end
            end
        end
    end
    #plot that:
    f= Figure(size=(1200,1000))
    cases2 = [L"\mathrm{𝖵}",
              L"\mathrm{\alpha_𝟣}",
              L"$\mathrm{\hat{𝗐}^𝟤_{𝟣𝟢,𝟧𝟢}}$ (𝖯𝖺𝗎𝗅)", 
              L"$\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}_{𝟣𝟢,𝟧𝟢}}$ (𝖯𝖺𝗎𝗅)",
              L"$\mathrm{\hat{𝗐}^𝟤_{𝟣𝟢,𝟧𝟢}}$ (𝖬𝗈𝗋𝗅𝖾𝗍)", 
              L"$\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}_{𝟣𝟢,𝟧𝟢}}$ (𝖬𝗈𝗋𝗅𝖾𝗍)", 
              L"$\mathrm{\hat{𝗐}^𝟤_{𝟤𝟢,𝟨𝟢}}$ (𝖯𝖺𝗎𝗅)", 
              L"$\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}_{𝟤𝟢,𝟨𝟢}}$ (𝖯𝖺𝗎𝗅)"]

    ax = Axis(f[1,1],
        title = "Significance of the observed number of EWS",
        xlabel = "Ice core",
        xticks = (1:7, ["NGRIP irregular", "NGRIP 5 y", "NGRIP 10 y", "NEEM 10 y", "NGRIP 20 y", "GRIP 20 y", "GISP2 20 y"]),
        xminorgridwidth = 1.0, 
        xminorgridcolor = :grey30 ,
        xminorgridvisible = true, xgridvisible = false,
        yminorticks = 1.5:1:8.5,
        ylabel = "Indicator",
        yticks = (1:8, cases2),
        yminorgridwidth = 1, 
        yminorgridcolor = :grey30, 
        yminorgridvisible = true, ygridvisible = false)
    
    n_max = maximum(skipmissing(justn))
    co = cgrad(:amp, n_max+1, categorical = true, alpha=0.9)

    hm = heatmap!(ax, justn, colormap = co, nan_color=:grey60)#, colorrange= (0,200))
    translate!(hm, 0,0,-100)
    ga = f[1,2] = GridLayout()
    co_bar = cgrad([:grey60,co...], n_max+2, categorical=true)
    Colorbar(ga[2,1], colormap = co_bar, colorrange = (0,n_max+2),
            ticks = (0.5:1:Int(n_max+1)+0.5, ["undefined", string.(0:1:Int(n_max))...]),
                label = "No. of significant EWS (p<0.05)",
                halign=:left)
    for i in 1:size(issig)[1]
        for j in 1:size(issig)[2]
            if issig[i,j] === 90
                scatter!(Point2f(i,j), strokecolor = :black, strokewidth=3, marker= :circle, markersize=40, color=:transparent)
            elseif issig[i,j] ===95
                scatter!(Point2f(i,j), strokecolor = :black, strokewidth=3, marker= :circle, markersize=40, color=:black)
            end
        end
    end
    elems = [[MarkerElement(strokecolor = :black, strokewidth=1, marker= :circle, markersize=20, color=:black)],
        [MarkerElement(strokecolor = :black, strokewidth=1, marker= :circle, markersize=20, color=:transparent)]]
        
    Legend(ga[1,1], elems, 
                ["significant at 95%", "significant at 90%"],
                rowgap = 20, 
                framevisible = false, #false,
                framecolor = :grey70,
                tellwidth = true,
                tellheight=true,
                orientation = :vertical,
                )
    colgap!(f.layout,40)
    if showing
        display(f)
    end
    if saving
        save(saveto, f)
    end
    
end


get_numbers_base_cases(showing=showing,saving=true, lowpass = true, saveto = "do_ews_across_greenland_ice_cores/figures/new_fig10.pdf")
get_numbers_base_cases(showing=showing,saving=true, lowpass = false, saveto = "do_ews_across_greenland_ice_cores/figures/new_figS14.pdf")



####### NEW aggregated plot #################

function aggregated_plot(;lowpass= lowpass, showing = true, saving = false, saveto = "paper/result_overview_n_inc_lowpass_$(lowpass)_smoothw_false_p_0.05.pdf")
    

    #cases = [L"\textbf{\mathrm{\sigma^𝟤}}",L"\mathrm{\alpha_𝟣}", L"$\mathrm{\hat{𝗐}^𝟤}$, 𝟤𝟢-𝟨𝟢 𝗒, 𝖯𝖺𝗎𝗅", L"$\mathrm{\hat{𝖧}_{𝟤𝟢-𝟨𝟢 𝗒, 𝖯𝖺𝗎𝗅}}$"]
    

    
    ngrip5_paths = ["new_surrogate_files/NGRIP5/var/w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/NGRIP5/ac/w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/NGRIP5/sca/s1_20_s2_60_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/NGRIP5/hurst/s1_20_s2_60_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]
    
    ngrip10_paths = ["new_surrogate_files/10y/var/w_200_normed_filt_NGRIP_lp_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/10y/ac/w_200_normed_filt_NGRIP_lp_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/10y/sca/s1_20_s2_60_NGRIP_lp_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/10y/hurst/s1_20_s2_60_NGRIP_lp_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]
    
    neem_paths = ["new_surrogate_files/10y/var/w_200_normed_filt_NEEM_lp_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/10y/ac/w_200_normed_filt_NEEM_lp_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/10y/sca/s1_20_s2_60_NEEM_lp_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/10y/hurst/s1_20_s2_60_NEEM_lp_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]
    
    ngrip20_paths = ["new_surrogate_files/20y/var/w_200_normed_filt_NGRIP_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/20y/ac/w_200_normed_filt_NGRIP_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/20y/sca/s1_20_s2_60_NGRIP_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/20y/hurst/s1_20_s2_60_NGRIP_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]
    
    grip_paths = ["new_surrogate_files/20y/var/w_200_normed_filt_GRIP_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/20y/ac/w_200_normed_filt_GRIP_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/20y/sca/s1_20_s2_60_GRIP_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/20y/hurst/s1_20_s2_60_GRIP_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]

    
    gisp2_paths = ["new_surrogate_files/20y/var/w_200_normed_filt_GISP2_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/20y/ac/w_200_normed_filt_GISP2_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/20y/sca/s1_20_s2_60_GISP2_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/20y/hurst/s1_20_s2_60_GISP2_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]

    if lowpass == false
        ngrip5_paths = ["new_surrogate_files/NGRIP5/var/w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                        "new_surrogate_files/NGRIP5/ac/w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                        "new_surrogate_files/NGRIP5/sca/s1_20_s2_60_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                        "new_surrogate_files/NGRIP5/hurst/s1_20_s2_60_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]

        ngrip10_paths = ["new_surrogate_files/10y/var/w_200_normed_filt_NGRIP_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                        "new_surrogate_files/10y/ac/w_200_normed_filt_NGRIP_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                        "new_surrogate_files/10y/sca/s1_20_s2_60_NGRIP_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                        "new_surrogate_files/10y/hurst/s1_20_s2_60_NGRIP_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]
        
        neem_paths = ["new_surrogate_files/10y/var/w_200_normed_filt_NEEM_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                        "new_surrogate_files/10y/ac/w_200_normed_filt_NEEM_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                        "new_surrogate_files/10y/sca/s1_20_s2_60_NEEM_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                        "new_surrogate_files/10y/hurst/s1_20_s2_60_NEEM_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]
    end


    f= Figure(size=(1200,500))


    cases2 = [L"\mathrm{𝖵𝖺𝗋}",
              L"\mathrm{\alpha_𝟣}",
              L"$\mathrm{\hat{𝗐}^𝟤_{𝟤𝟢,𝟨𝟢}}$",# (𝖯𝖺𝗎𝗅)", 
              L"$\mathrm{\hat{𝖧}_{𝟤𝟢,𝟨𝟢}}$"]# (𝖯𝖺𝗎𝗅)"]

    ax = Axis(f[1,1],
        xticks = (1:17, ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]) ,
        xreversed = true,
        yreversed=true,
        xlabel = "Event",
        xminorgridwidth = 3.5,
        xgridwidth=0.7,
        xminorgridcolor = :grey10 ,
        xgridcolor=:grey30,
        xminorgridvisible = true, 
        xgridvisible = true,
        yticks = (1:6, ["NGRIP 5 y", "NGRIP 10 y", "NEEM 10 y", "NGRIP 20 y", "GRIP 20 y", "GISP2 20 y"]),
        ylabel = "Ice core",
        yminorgridwidth = 4.5,
        ygridwidth=0.7,
        yminorgridcolor = :grey10 ,
        ygridcolor=:grey30,
        yminorgridvisible = true, 
        ygridvisible = true,
        )

    
    new_mat = zeros(2*17,2*6)
    for (i,core) in enumerate([ngrip5_paths, ngrip10_paths, neem_paths, ngrip20_paths, grip_paths, gisp2_paths])
        for (j,path) in enumerate(core)     
            ind = load(path)["slopes"]
            for ev in 1:17
                if j==1
                    new_mat[2*ev,2*i-1]=1
                elseif j==2
                    new_mat[2*ev-1,2*i-1]=2
                elseif j==3
                    new_mat[2*ev,2*i]=3
                elseif j==4
                    new_mat[2*ev-1,2*i]=4
                end
                if typeof(ind.slopes[ev]) != Nothing
                    if ind.slopes[ev]>0 && ind.p_one[ev] < plim
                        if j==1
                            text!(ax, ev+0.25,i-0.25, text = L"𝖵",align = (:center, :center))
                            new_mat[2*ev,2*i-1]+=4
                        elseif j==2
                            text!(ax, ev-0.25,i-0.25, text = L"\alpha_𝟣",align = (:center, :center))
                            new_mat[2*ev-1,2*i-1]+=4
                        elseif j==3
                            text!(ax, ev+0.25,i+0.25, text = L"\hat{𝗐}^𝟤",align = (:center, :center))
                            new_mat[2*ev,2*i]+=4
                        elseif j==4
                            text!(ax, ev-0.25,i+0.25, text = L"\hat{𝖧}^{\text{𝗅𝗈𝖼}}",align = (:center, :center))
                            new_mat[2*ev-1,2*i]+=4
                        end
                    end
                end
            end
            
        end
    end
    #t1 = cgrad(:managua, 4, categorical=true )
    #t2 = cgrad(:roma, 4, categorical=true )
    t3 = cgrad(:Spectral_4, 4, categorical=true )
    co2 = cgrad([[i for i in t3]...,[i for i in t3]...],8, categorical=true, alpha=[0.08,0.08,0.08,0.08,0.8,0.8,0.8,0.8])   
    
    #co2 = cgrad([:seagreen3, :dodgerblue3, :goldenrod1, :firebrick, :seagreen3, :dodgerblue3, :goldenrod1, :firebrick], 8, categorical = true, alpha=[0.1,0.1,0.1,0.1,0.6,0.6,0.6,0.6])#[0.1,0.1,0.1,0.1,1,1,1,1])    
    ##co2 = cgrad([[i for i in t1]...,[i for i in t1]...],8, categorical=true, alpha=[0.1,0.1,0.1,0.1,0.7,0.7,0.7,0.7])
    #co2 = cgrad([t2[3],t2[4],t2[1],t2[1],t2[3],t2[4],t2[1],t2[2]],8, categorical=true, alpha=[0.08,0.08,0.08,0.08,0.75,0.75,0.75,0.75])    
    #co2 = cgrad([[i for i in t3]...,[i for i in t3]...],8, categorical=true, alpha=[0.1,0.1,0.1,0.1,0.75,0.75,0.75,0.75])     
    #co2 = cgrad([:steelblue1, :steelblue2, :steelblue3, :steelblue, :firebrick1, :firebrick2, :firebrick3, :firebrick], 8, categorical = true, alpha=[0.2,0.2,0.2,0.2,0.8,0.8,0.8,0.8])
    #co2 = cgrad([:grey70, :grey70, :grey70, :grey70, :firebrick3, :firebrick3, :firebrick3, :firebrick3], 8, categorical = true, alpha=[0.1,0.1,0.1,0.1,0.8,0.8,0.8,0.8])
    hm2 = heatmap!(ax, 0.5:0.5:17.5,0.5:0.5:6.5, new_mat, colormap = co2, colorrange=extrema(new_mat))
    #Colorbar(f[1,2],hm2)
    translate!(hm2, 0,0,-100)

    elems= [
        [PolyElement(color = co2[3], strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 0.5), (0.5, 0.5), (0.5, 0),(0,0)]),
        PolyElement(color = co2[4], strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 0.5), (1, 0.5), (1, 0),(0.5,0)]),
        PolyElement(color = co2[1], strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 1), (0.5, 1), (0.5, 0.5),(0,0.5)]),
        PolyElement(color = co2[2], strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 1), (1, 1), (1, 0.5),(0.5,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey10, strokewidth = 2, points = Point2f[(0, 1), (1, 1), (1, 0),(0,0)])
        ],
        [PolyElement(color = co2[3], strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 0.5), (0.5, 0.5), (0.5, 0),(0,0)]),
        PolyElement(color = co2[4], strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 0.5), (1, 0.5), (1, 0),(0.5,0)]),
        PolyElement(color = co2[5], strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 1), (0.5, 1), (0.5, 0.5),(0,0.5)]),
        PolyElement(color = co2[2], strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 1), (1, 1), (1, 0.5),(0.5,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey10, strokewidth = 2, points = Point2f[(0, 1), (1, 1), (1, 0),(0,0)])
        ],
        [PolyElement(color = co2[3], strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 0.5), (0.5, 0.5), (0.5, 0),(0,0)]),
        PolyElement(color = co2[4], strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 0.5), (1, 0.5), (1, 0),(0.5,0)]),
        PolyElement(color = co2[1], strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 1), (0.5, 1), (0.5, 0.5),(0,0.5)]),
        PolyElement(color = co2[6], strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 1), (1, 1), (1, 0.5),(0.5,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey10, strokewidth = 2, points = Point2f[(0, 1), (1, 1), (1, 0),(0,0)])
        ],
        [PolyElement(color = co2[7], strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 0.5), (0.5, 0.5), (0.5, 0),(0,0)]),
        PolyElement(color = co2[4], strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 0.5), (1, 0.5), (1, 0),(0.5,0)]),
        PolyElement(color = co2[1], strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 1), (0.5, 1), (0.5, 0.5),(0,0.5)]),
        PolyElement(color = co2[2], strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 1), (1, 1), (1, 0.5),(0.5,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey10, strokewidth = 2, points = Point2f[(0, 1), (1, 1), (1, 0),(0,0)])
        ],
        [PolyElement(color = co2[3], strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 0.5), (0.5, 0.5), (0.5, 0),(0,0)]),
        PolyElement(color = co2[8], strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 0.5), (1, 0.5), (1, 0),(0.5,0)]),
        PolyElement(color = co2[1], strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 1), (0.5, 1), (0.5, 0.5),(0,0.5)]),
        PolyElement(color = co2[2], strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 1), (1, 1), (1, 0.5),(0.5,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey10, strokewidth = 2, points = Point2f[(0, 1), (1, 1), (1, 0),(0,0)])]
    ]

     Legend(f[1,2], elems, [L"𝗇𝗈𝗇𝖾",  L"𝖵", L"\alpha_𝟣", L"\hat{𝗐}^𝟤", L"\hat{𝖧}^{\text{𝗅𝗈𝖼}}"],#[L"𝗇𝗈𝗇𝖾 $ $", cases2...], 
                "significant increase in",
                titlefont=:regular,
                titlegap=15,
                rowgap = 20, 
                framevisible = false, #false,
                framecolor = :grey70,
                tellwidth = true,
                #tellheight=true,
                orientation = :vertical,
                patchsize=(30,40),
                patchlabelgap=10,
                )
    colgap!(f.layout,40)
    if showing
        display(f)
    end
    if saving
        save(saveto, f)
    end
    
end
aggregated_plot(lowpass=true, saving=false, saveto="do_ews_across_greenland_ice_cores/figures/aggregated_overview_lp_colourful2.pdf")
aggregated_plot(lowpass=false, saving=false, saveto="do_ews_across_greenland_ice_cores/figures/aggregated_overview_no_lp_colourful2.pdf")


function aggregated_plot2(;lowpass= lowpass, showing = true, saving = false, saveto = "paper/result_overview_n_inc_lowpass_$(lowpass)_smoothw_false_p_0.05.pdf")
    #cases = [L"\textbf{\mathrm{\sigma^𝟤}}",L"\mathrm{\alpha_𝟣}", L"$\mathrm{\hat{𝗐}^𝟤}$, 𝟤𝟢-𝟨𝟢 𝗒, 𝖯𝖺𝗎𝗅", L"$\mathrm{\hat{𝖧}_{𝟤𝟢-𝟨𝟢 𝗒, 𝖯𝖺𝗎𝗅}}$"]
    
    ngrip5_paths = ["new_surrogate_files/NGRIP5/var/w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/NGRIP5/ac/w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/NGRIP5/sca/s1_20_s2_60_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/NGRIP5/hurst/s1_20_s2_60_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]
    
    ngrip10_paths = ["new_surrogate_files/10y/var/w_200_normed_filt_NGRIP_lp_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/10y/ac/w_200_normed_filt_NGRIP_lp_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/10y/sca/s1_20_s2_60_NGRIP_lp_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/10y/hurst/s1_20_s2_60_NGRIP_lp_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]
    
    neem_paths = ["new_surrogate_files/10y/var/w_200_normed_filt_NEEM_lp_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/10y/ac/w_200_normed_filt_NEEM_lp_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/10y/sca/s1_20_s2_60_NEEM_lp_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/10y/hurst/s1_20_s2_60_NEEM_lp_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]
    
    ngrip20_paths = ["new_surrogate_files/20y/var/w_200_normed_filt_NGRIP_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/20y/ac/w_200_normed_filt_NGRIP_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/20y/sca/s1_20_s2_60_NGRIP_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/20y/hurst/s1_20_s2_60_NGRIP_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]
    
    grip_paths = ["new_surrogate_files/20y/var/w_200_normed_filt_GRIP_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/20y/ac/w_200_normed_filt_GRIP_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/20y/sca/s1_20_s2_60_GRIP_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/20y/hurst/s1_20_s2_60_GRIP_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]

    
    gisp2_paths = ["new_surrogate_files/20y/var/w_200_normed_filt_GISP2_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/20y/ac/w_200_normed_filt_GISP2_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/20y/sca/s1_20_s2_60_GISP2_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/20y/hurst/s1_20_s2_60_GISP2_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]

    if lowpass == false
        ngrip5_paths = ["new_surrogate_files/NGRIP5/var/w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                        "new_surrogate_files/NGRIP5/ac/w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                        "new_surrogate_files/NGRIP5/sca/s1_20_s2_60_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                        "new_surrogate_files/NGRIP5/hurst/s1_20_s2_60_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]

        ngrip10_paths = ["new_surrogate_files/10y/var/w_200_normed_filt_NGRIP_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                        "new_surrogate_files/10y/ac/w_200_normed_filt_NGRIP_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                        "new_surrogate_files/10y/sca/s1_20_s2_60_NGRIP_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                        "new_surrogate_files/10y/hurst/s1_20_s2_60_NGRIP_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]
        
        neem_paths = ["new_surrogate_files/10y/var/w_200_normed_filt_NEEM_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                        "new_surrogate_files/10y/ac/w_200_normed_filt_NEEM_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                        "new_surrogate_files/10y/sca/s1_20_s2_60_NEEM_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                        "new_surrogate_files/10y/hurst/s1_20_s2_60_NEEM_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]
    end


    f= Figure(size=(1200,600))

    function get_n(x,plim)
        n=0
        for ev in 1:17
            if typeof(x.slopes[ev]) != Nothing
                if x.slopes[ev] >0  && x.p_one[ev] < plim
                    n+=1
                end
            end
        end
        return n
    end


    cases2 = [L"\mathrm{𝖵𝖺𝗋}",
              L"\mathrm{\alpha_𝟣}",
              L"$\mathrm{\hat{𝗐}^𝟤_{𝟤𝟢,𝟨𝟢}}$",# (𝖯𝖺𝗎𝗅)", 
              L"$\mathrm{\hat{𝖧}_{𝟤𝟢,𝟨𝟢}}$"]# (𝖯𝖺𝗎𝗅)"]

    gl = f[1,1] = GridLayout()
    ax = Axis(gl[1,1],
        xticks = (1:17, ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]) ,
        xreversed = true,
        yreversed=true,
        xlabel = "Event",
        xminorgridwidth = 3.5,
        xgridwidth=0.7,
        xminorgridcolor = :grey10 ,
        xgridcolor=:grey30,
        xminorgridvisible = true, 
        xgridvisible = true,
        yticks = (1:6, ["NGRIP 5 y", "NGRIP 10 y", "NEEM 10 y", "NGRIP 20 y", "GRIP 20 y", "GISP2 20 y"]),
        ylabel = "Ice core",
        yminorgridwidth = 4.5,
        ygridwidth=0.7,
        yminorgridcolor = :grey10 ,
        ygridcolor=:grey30,
        yminorgridvisible = true, 
        ygridvisible = true,
        )
    
    gr = f[1,2] = GridLayout(width=80)
    ax2 = Axis(gr[1,1],
        #xticks = (1:17, ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]) ,
        xticks=[0.5],
        xreversed = true,
        yreversed=true,
        #xlabel = "Transition",
        xminorgridwidth = 3.5,
        xgridwidth=0.7,
        xminorgridcolor = :grey10 ,
        xgridcolor=:grey30,
        xminorgridvisible = true, 
        xgridvisible = true,
        yticks = (1:6, ["NGRIP 5 y", "NGRIP 10 y", "NEEM 10 y", "NGRIP 20 y", "GRIP 20 y", "GISP2 20 y"]),
        #ylabel = "Method",
        yminorgridwidth = 4.5,
        ygridwidth=0.7,
        yminorgridcolor = :grey10 ,
        ygridcolor=:grey30,
        yminorgridvisible = true, 
        ygridvisible = true,
        xticksvisible = false,
        xticklabelsvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false
        )

    
    new_mat = zeros(2*17,2*6)
    num_mat = zeros(2,2*6)
    for (i,core) in enumerate([ngrip5_paths, ngrip10_paths, neem_paths, ngrip20_paths, grip_paths, gisp2_paths])
        for (j,path) in enumerate(core)     
            ind = load(path)["slopes"]
            n_h = get_n(ind,0.05)
            t_c =  :black #:grey60
            if n_h >= 4
                t_f = :bold
                f_s = 16
                #t_c =  :grey60
            else
                t_f = :regular
                f_s = 15
                #t_c =  :black
            end

            if j==1
                text!(ax2, 0.75,i-0.25, text = "$n_h",align = (:center, :center), color=t_c, font = t_f, fontsize = f_s)
                num_mat[2,2*i-1]=n_h   
            elseif j==2
                text!(ax2, 0.25,i-0.25, text = "$n_h",align = (:center, :center), color=t_c, font = t_f, fontsize = f_s)
                num_mat[1,2*i-1]=n_h
            elseif j==3
                text!(ax2, 0.75,i+0.25, text = "$n_h",align = (:center, :center), color=t_c, font = t_f, fontsize = f_s)
                num_mat[2,2*i]=n_h
            elseif j==4
                text!(ax2, 0.25,i+0.25, text = "$n_h",align = (:center, :center), color=t_c, font = t_f, fontsize = f_s)
                num_mat[1,2*i]=n_h
            end

            for ev in 1:17
                # if j==1
                #     new_mat[2*ev,2*i-1]=1
                # elseif j==2
                #     new_mat[2*ev-1,2*i-1]=2
                # elseif j==3
                #     new_mat[2*ev,2*i]=3
                # elseif j==4
                #     new_mat[2*ev-1,2*i]=4
                # end
                if typeof(ind.slopes[ev]) != Nothing
                    if ind.slopes[ev]>0
                        if j==1
                            #text!(ax, ev+0.25,i-0.25, text = L"𝖵",align = (:center, :center))
                            new_mat[2*ev,2*i-1]+=1
                        elseif j==2
                            #text!(ax, ev-0.25,i-0.25, text = L"\alpha_𝟣",align = (:center, :center))
                            new_mat[2*ev-1,2*i-1]+=1
                        elseif j==3
                            #text!(ax, ev+0.25,i+0.25, text = L"\hat{𝗐}^𝟤",align = (:center, :center))
                            new_mat[2*ev,2*i]+=1
                        elseif j==4
                            #text!(ax, ev-0.25,i+0.25, text = L"\hat{𝖧}^{\text{𝗅𝗈𝖼}}",align = (:center, :center))
                            new_mat[2*ev-1,2*i]+=1
                        end
                        
                        tc2 = :black
                        fs2 = 16
                        if ind.p_one[ev] < plim
                            if j==1
                                text!(ax, ev+0.25,i-0.25, text = L"\mathbf{𝖵}",align = (:center, :center), color=tc2, fontsize = fs2)
                                new_mat[2*ev,2*i-1]+=1
                            elseif j==2
                                text!(ax, ev-0.25,i-0.25, text = L"\mathbf{\alpha_𝟣}",align = (:center, :center), color=tc2, fontsize = fs2)
                                new_mat[2*ev-1,2*i-1]+=1
                            elseif j==3
                                text!(ax, ev+0.25,i+0.25, text = L"\mathbf{\hat{𝗐}^𝟤}",align = (:center, :center), color=tc2, fontsize = fs2)
                                new_mat[2*ev,2*i]+=1
                            elseif j==4
                                text!(ax, ev-0.25,i+0.25, text = L"\mathbf{\hat{𝖧}^{\text{𝗅𝗈𝖼}}}",align = (:center, :center), color=tc2, fontsize = fs2)
                                new_mat[2*ev-1,2*i]+=1
                            end
                        end
                    end
                end
            end
            
        end
    end
    #t1 = cgrad(:managua, 4, categorical=true )
    #t2 = cgrad(:roma, 4, categorical=true )
    #t3 = cgrad(:Spectral_4, 4, categorical=true )
    #co2 = cgrad([[i for i in t3]...,[i for i in t3]...],8, categorical=true, alpha=[0.08,0.08,0.08,0.08,0.8,0.8,0.8,0.8])   
    
    #co2 = cgrad([:seagreen3, :dodgerblue3, :goldenrod1, :firebrick, :seagreen3, :dodgerblue3, :goldenrod1, :firebrick], 8, categorical = true, alpha=[0.1,0.1,0.1,0.1,0.6,0.6,0.6,0.6])#[0.1,0.1,0.1,0.1,1,1,1,1])    
    ##co2 = cgrad([[i for i in t1]...,[i for i in t1]...],8, categorical=true, alpha=[0.1,0.1,0.1,0.1,0.7,0.7,0.7,0.7])
    #co2 = cgrad([t2[3],t2[4],t2[1],t2[1],t2[3],t2[4],t2[1],t2[2]],8, categorical=true, alpha=[0.08,0.08,0.08,0.08,0.75,0.75,0.75,0.75])    
    #co2 = cgrad([[i for i in t3]...,[i for i in t3]...],8, categorical=true, alpha=[0.1,0.1,0.1,0.1,0.75,0.75,0.75,0.75])     
    #co2 = cgrad([:steelblue1, :steelblue2, :steelblue3, :steelblue, :firebrick1, :firebrick2, :firebrick3, :firebrick], 8, categorical = true, alpha=[0.2,0.2,0.2,0.2,0.8,0.8,0.8,0.8])
    # co2 = cgrad([:grey70, :grey70, :grey70, :grey70, :firebrick3, :firebrick3, :firebrick3, :firebrick3], 8, categorical = true, alpha=[0.1,0.1,0.1,0.1,0.8,0.8,0.8,0.8])
    
    co2 = cgrad([:steelblue, :darkred, :firebrick], 3, categorical = true, alpha=[0.12,0.12,0.8])
    #co2 = cgrad([:steelblue, :darkred, :darkred], 3, categorical = true, alpha=[0.15,0.15,1])
    hm2 = heatmap!(ax, 0.5:0.5:17.5,0.5:0.5:6.5, new_mat, colormap = co2, colorrange=extrema(new_mat))
    #Colorbar(f[1,2],hm2)
    translate!(hm2, 0,0,-100)

    n_max = Int(maximum(num_mat))
    #@show n_max
    co = cgrad(:amp, n_max+1, categorical = true, alpha=0.55)
    # if lowpass
    #     co = cgrad(:amp, n_max+1, categorical = true, alpha=[0.7,0.7,0.7,0.7,0.9,0.9])
    # else
    #     co = cgrad(:amp, n_max+1, categorical = true, alpha=[0.7,0.7,0.7,0.7,0.9])
    # end
    # hm_n = heatmap!(ax2, 0.5:0.5:1,0.5:0.5:4.5, num_mat, colormap = co)
    hm_n = heatmap!(ax2, 0:0.5:1,0.5:0.5:6.5, num_mat, colormap = co, colorrange=(0,n_max))
    #hm_n = heatmap!(ax2, 0:0.5:1,0.5:0.5:4.5, num_mat, colormap =Makie.Categorical(:amp),colorrange=(0,n_max))
    Colorbar(gr[1,end+1], #hm_n,
        colormap = co, 
        colorrange=(0,n_max+1),
        ticks = (0.5:1:n_max+0.5, string.(0:1:n_max)), 
        #ticks = 0.5:1:n_max+1.5,
        label = "No. of significant EWS (p<0.05)",
        halign=:left)

    translate!(hm_n, 0,0,-100)



    Label(gl[1,1,TopLeft()], "(a)", fontsize = 20,
                        font = :bold, padding = (0,20,-10,15),
                        halign = :left,valign =:bottom)
    Label(gr[1,1,TopLeft()], "(b)", fontsize = 20,
                        font = :bold, padding = (0,20,-10,15),
                        halign = :left,valign =:bottom)
    elems2 = [
        [PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 0.5), (0.5, 0.5), (0.5, 0),(0,0)]),
        #MarkerElement(color = :black, marker = 'x', markersize = 15, points = Point2f[(0.25, 0.75)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 0.5), (1, 0.5), (1, 0),(0.5,0)]),
        PolyElement(color = :grey20, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 1), (0.5, 1), (0.5, 0.5),(0,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 1), (1, 1), (1, 0.5),(0.5,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey10, strokewidth = 2, points = Point2f[(0, 1), (1, 1), (1, 0),(0,0)])
        ],
        [PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 0.5), (0.5, 0.5), (0.5, 0),(0,0)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 0.5), (1, 0.5), (1, 0),(0.5,0)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 1), (0.5, 1), (0.5, 0.5),(0,0.5)]),
        PolyElement(color = :grey20, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 1), (1, 1), (1, 0.5),(0.5,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey10, strokewidth = 2, points = Point2f[(0, 1), (1, 1), (1, 0),(0,0)])
        ],
        [PolyElement(color = :grey20, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 0.5), (0.5, 0.5), (0.5, 0),(0,0)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 0.5), (1, 0.5), (1, 0),(0.5,0)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 1), (0.5, 1), (0.5, 0.5),(0,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 1), (1, 1), (1, 0.5),(0.5,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey10, strokewidth = 2, points = Point2f[(0, 1), (1, 1), (1, 0),(0,0)])
        ],
        [PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 0.5), (0.5, 0.5), (0.5, 0),(0,0)]),
        PolyElement(color = :grey20, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 0.5), (1, 0.5), (1, 0),(0.5,0)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 1), (0.5, 1), (0.5, 0.5),(0,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 1), (1, 1), (1, 0.5),(0.5,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey10, strokewidth = 2, points = Point2f[(0, 1), (1, 1), (1, 0),(0,0)])],
        [MarkerElement(color = co2[1], marker = :rect, markersize = 25,
        strokecolor = :grey, strokewidth = 0.7)],
        [MarkerElement(color =co2[2], marker = :rect, markersize = 25,
        strokecolor = :grey, strokewidth = 0.7)],
        [MarkerElement(color = co2[3], marker = :rect, markersize = 25,
        strokecolor = :grey, strokewidth = 0.7)],
        ]
    Legend(f[end+1,1], elems2, [ L"𝖵", L"\alpha_𝟣", L"\hat{𝗐}_{𝟤𝟢,𝟨𝟢}^𝟤", L"\hat{𝖧}_{𝟤𝟢,𝟨𝟢}^{\text{𝗅𝗈𝖼}}", "decreasing", "increasing","significantly increasing (p<$(plim))"],
                "EWS indicators",
                titleposition=:top,
                #titlefont=:regular,
                #titlegap=15,
                rowgap = 15, 
                framevisible = false, #false,
                framecolor = :grey70,
                tellwidth = false,
                tellheight=true,
                orientation = :horizontal,
                patchsize=(30,40),
                patchlabelgap=10,
                )
    
    colgap!(f.layout,20)
    if showing
        display(f)
    end
    if saving
        save(saveto, f)
    end
    
end
aggregated_plot2(lowpass = true, saving=true, saveto="do_ews_across_greenland_ice_cores/figures/aggregated_overview2_incl_n_lp.pdf")
aggregated_plot2(lowpass = false, saving=true, saveto="do_ews_across_greenland_ice_cores/figures/aggregated_overview2_incl_n_no_lp.pdf")


function aggregated_plot2_method(;lowpass= lowpass, showing = true, saving = false, saveto = "paper/result_overview_method_lowpass_$(lowpass)_smoothw_false_p_0.05.pdf")
    #cases = [L"\textbf{\mathrm{\sigma^𝟤}}",L"\mathrm{\alpha_𝟣}", L"$\mathrm{\hat{𝗐}^𝟤}$, 𝟤𝟢-𝟨𝟢 𝗒, 𝖯𝖺𝗎𝗅", L"$\mathrm{\hat{𝖧}_{𝟤𝟢-𝟨𝟢 𝗒, 𝖯𝖺𝗎𝗅}}$"]
    
    boers_paths = ["new_surrogate_files/NGRIP5/var/NIKLAS_w_200_N_lowpass_py2_FiltInd_true_onlyfull_false_detrend_GLSAR_DemeanSur_true_10000_FOURIER.jld2",
                    "new_surrogate_files/NGRIP5/ac/NIKLAS_w_200_N_lowpass_py2_FiltInd_true_onlyfull_false_detrend_GLSAR_DemeanSur_true_10000_FOURIER.jld2",
                    "new_surrogate_files/NGRIP5/sca/NIKLAS_s1_10_s2_50_N_lowpass_py2_PAUL_FiltInd_true_onlyfull_false_detrend_GLSAR_DemeanSur_true_smoothw_true_10000_FOURIER.jld2",
                    "new_surrogate_files/NGRIP5/hurst/NIKLAS_s1_10_s2_50_N_lowpass_py2_PAUL_FiltInd_true_onlyfull_false_detrend_GLSAR_DemeanSur_true_smoothw_true_10000_FOURIER.jld2"]
    
    step1_paths = ["new_surrogate_files/NGRIP5/var/w_200_normed_filt_N_lowpass_py2_gs_short_FiltInd_true_onlyfull_false_10000_TFTS.jld2",
                    "new_surrogate_files/NGRIP5/ac/w_200_normed_filt_N_lowpass_py2_gs_short_FiltInd_true_onlyfull_false_10000_TFTS.jld2",
                    "new_surrogate_files/NGRIP5/sca/s1_10_s2_50_N_lowpass_py2_PAUL_normed_gs_short_FiltInd_true_normalise_false_nocoi_false_onlyfull_false_smoothw_true_10000_TFTS.jld2",
                    "new_surrogate_files/NGRIP5/hurst/s1_10_s2_50_N_lowpass_py2_PAUL_normed_gs_short_FiltInd_true_normalise_false_nocoi_false_onlyfull_false_smoothw_true_10000_TFTS.jld2"]
    
    step2_paths = ["new_surrogate_files/NGRIP5/var/w_200_normed_filt_N_lowpass_py2_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/NGRIP5/ac/w_200_normed_filt_N_lowpass_py2_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/NGRIP5/sca/s1_10_s2_50_N_lowpass_py2_PAUL_normed_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_false_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/NGRIP5/hurst/s1_10_s2_50_N_lowpass_py2_PAUL_normed_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_false_smoothw_false_10000_TFTS.jld2"]
                    
    step3_paths = ["new_surrogate_files/NGRIP5/var/w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/NGRIP5/ac/w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/NGRIP5/sca/s1_10_s2_50_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                    "new_surrogate_files/NGRIP5/hurst/s1_10_s2_50_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]
    
   #irregular_paths = []

    
    #Morlet_paths = []

    if lowpass == false
       
        step3_paths = ["new_surrogate_files/NGRIP5/var/w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                        "new_surrogate_files/NGRIP5/ac/w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                        "new_surrogate_files/NGRIP5/sca/s1_10_s2_50_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                        "new_surrogate_files/NGRIP5/hurst/s1_10_s2_50_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]

        #irregular_paths = []


        #Morlet_paths = []
    end

    
    f= Figure(size=(1200,500))

    function get_n(x,plim)
        n=0
        for ev in 1:17
            if typeof(x.slopes[ev]) != Nothing
                if x.slopes[ev] >0  && x.p_one[ev] < plim
                    n+=1
                end
            end
        end
        return n
    end

    cases2 = [L"\mathrm{𝖵𝖺𝗋}",
              L"\mathrm{\alpha_𝟣}",
              L"$\mathrm{\hat{𝗐}^𝟤_{𝟤𝟢,𝟨𝟢}}$",# (𝖯𝖺𝗎𝗅)", 
              L"$\mathrm{\hat{𝖧}_{𝟤𝟢,𝟨𝟢}}$"]# (𝖯𝖺𝗎𝗅)"]
    
    gl = f[1,1] = GridLayout()
    ax = Axis(gl[1,1],
        xticks = (1:17, ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]) ,
        xreversed = true,
        yreversed=true,
        xlabel = "Event",
        xminorgridwidth = 3.5,
        xgridwidth=0.7,
        xminorgridcolor = :grey10 ,
        xgridcolor=:grey30,
        xminorgridvisible = true, 
        xgridvisible = true,
        yticks = (1:4, ["Boers, 2018", "Modified significance testing\n (Step 1)", "Modified EWS calculation\n (Step 2)", "Modified data preprocessing\n(Step 3)"]),
        #ylabel = "Method",
        yminorgridwidth = 4.5,
        ygridwidth=0.7,
        yminorgridcolor = :grey10 ,
        ygridcolor=:grey30,
        yminorgridvisible = true, 
        ygridvisible = true,
        )
    
    gr = f[1,2] = GridLayout(width=80)
    ax2 = Axis(gr[1,1],
        #xticks = (1:17, ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]) ,
        xticks=[0.5],
        xreversed = true,
        yreversed=true,
        #xlabel = "Transition",
        xminorgridwidth = 3.5,
        xgridwidth=0.7,
        xminorgridcolor = :grey10 ,
        xgridcolor=:grey30,
        xminorgridvisible = true, 
        xgridvisible = true,
        yticks = (1:4, ["Boers, 2018", "Modified significance testing\n (Step 1)", "Modified EWS calculation\n (Step 2)", "Modified data preprocessing\n(Step 3)"]),
        #ylabel = "Method",
        yminorgridwidth = 4.5,
        ygridwidth=0.7,
        yminorgridcolor = :grey10 ,
        ygridcolor=:grey30,
        yminorgridvisible = true, 
        ygridvisible = true,
        xticksvisible = false,
        xticklabelsvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false
        )

    
    new_mat = zeros(2*17,2*4)
    num_mat = zeros(2,2*4)
    for (i,core) in enumerate([boers_paths, step1_paths, step2_paths, step3_paths])
        for (j,path) in enumerate(core)     
            ind = load(path)["slopes"]
            n_h = get_n(ind,0.05)
            t_c =  :black #:grey60
            if n_h >= 4
                t_f = :bold
                f_s = 16
                #t_c =  :grey60
            else
                t_f = :regular
                f_s = 15
                #t_c =  :black
            end

            if j==1
                text!(ax2, 0.75,i-0.25, text = "$n_h",align = (:center, :center), color=t_c, font = t_f, fontsize = f_s)
                num_mat[2,2*i-1]=n_h   
            elseif j==2
                text!(ax2, 0.25,i-0.25, text = "$n_h",align = (:center, :center), color=t_c, font = t_f, fontsize = f_s)
                num_mat[1,2*i-1]=n_h
            elseif j==3
                text!(ax2, 0.75,i+0.25, text = "$n_h",align = (:center, :center), color=t_c, font = t_f, fontsize = f_s)
                num_mat[2,2*i]=n_h
            elseif j==4
                text!(ax2, 0.25,i+0.25, text = "$n_h",align = (:center, :center), color=t_c, font = t_f, fontsize = f_s)
                num_mat[1,2*i]=n_h
            end
            # if n_h >=4
            #     if j==1
            #         scatter!(ax2, [0.75],[i-0.25], strokecolor = :black, strokewidth=3, marker= :circle, markersize=20, color=:black)
                  
            #     elseif j==2
            #         text!(ax2, 0.25,i-0.25, text =string(n_h),align = (:center, :center), color=:grey70)
                 
            #     elseif j==3
            #         text!(ax2, 0.75,i+0.25, text = string(n_h),align = (:center, :center), color=:grey70)
            #     elseif j==4
            #         text!(ax2, 0.25,i+0.25, text = string(n_h),align = (:center, :center), color=:grey70)
    
            #     end
            # end
            for ev in 1:17
      
                if typeof(ind.slopes[ev]) != Nothing
                    if ind.slopes[ev]>0
                        if j==1
                            #text!(ax, ev+0.25,i-0.25, text = L"𝖵",align = (:center, :center))
                            new_mat[2*ev,2*i-1]+=1
                        elseif j==2
                            #text!(ax, ev-0.25,i-0.25, text = L"\alpha_𝟣",align = (:center, :center))
                            new_mat[2*ev-1,2*i-1]+=1
                        elseif j==3
                            #text!(ax, ev+0.25,i+0.25, text = L"\hat{𝗐}^𝟤",align = (:center, :center))
                            new_mat[2*ev,2*i]+=1
                        elseif j==4
                            #text!(ax, ev-0.25,i+0.25, text = L"\hat{𝖧}^{\text{𝗅𝗈𝖼}}",align = (:center, :center))
                            new_mat[2*ev-1,2*i]+=1
                        end
                        tc2 = :black
                        fs2 = 16
                        if ind.p_one[ev] < plim
                            if j==1
                                text!(ax, ev+0.25,i-0.25, text = L"\mathbf{𝖵}",align = (:center, :center), color=tc2, fontsize = fs2)
                                new_mat[2*ev,2*i-1]+=1
                            elseif j==2
                                text!(ax, ev-0.25,i-0.25, text = L"\mathbf{\alpha_𝟣}",align = (:center, :center), color=tc2, fontsize = fs2)
                                new_mat[2*ev-1,2*i-1]+=1
                            elseif j==3
                                text!(ax, ev+0.25,i+0.25, text = L"\mathbf{\hat{𝗐}^𝟤}",align = (:center, :center), color=tc2, fontsize = fs2)
                                new_mat[2*ev,2*i]+=1
                            elseif j==4
                                text!(ax, ev-0.25,i+0.25, text = L"\mathbf{\hat{𝖧}^{\text{𝗅𝗈𝖼}}}",align = (:center, :center), color=tc2, fontsize = fs2)
                                new_mat[2*ev-1,2*i]+=1
                            end
                        end
                    end
                end
            end
            
        end
    end
    #t1 = cgrad(:managua, 4, categorical=true )
    #t2 = cgrad(:roma, 4, categorical=true )
    #t3 = cgrad(:Spectral_4, 4, categorical=true )
    #co2 = cgrad([[i for i in t3]...,[i for i in t3]...],8, categorical=true, alpha=[0.08,0.08,0.08,0.08,0.8,0.8,0.8,0.8])   
    
    #co2 = cgrad([:seagreen3, :dodgerblue3, :goldenrod1, :firebrick, :seagreen3, :dodgerblue3, :goldenrod1, :firebrick], 8, categorical = true, alpha=[0.1,0.1,0.1,0.1,0.6,0.6,0.6,0.6])#[0.1,0.1,0.1,0.1,1,1,1,1])    
    ##co2 = cgrad([[i for i in t1]...,[i for i in t1]...],8, categorical=true, alpha=[0.1,0.1,0.1,0.1,0.7,0.7,0.7,0.7])
    #co2 = cgrad([t2[3],t2[4],t2[1],t2[1],t2[3],t2[4],t2[1],t2[2]],8, categorical=true, alpha=[0.08,0.08,0.08,0.08,0.75,0.75,0.75,0.75])    
    #co2 = cgrad([[i for i in t3]...,[i for i in t3]...],8, categorical=true, alpha=[0.1,0.1,0.1,0.1,0.75,0.75,0.75,0.75])     
    #co2 = cgrad([:steelblue1, :steelblue2, :steelblue3, :steelblue, :firebrick1, :firebrick2, :firebrick3, :firebrick], 8, categorical = true, alpha=[0.2,0.2,0.2,0.2,0.8,0.8,0.8,0.8])
    # co2 = cgrad([:grey70, :grey70, :grey70, :grey70, :firebrick3, :firebrick3, :firebrick3, :firebrick3], 8, categorical = true, alpha=[0.1,0.1,0.1,0.1,0.8,0.8,0.8,0.8])
    co2 = cgrad([:steelblue, :darkred, :firebrick], 3, categorical = true, alpha=[0.12,0.12,0.8])
    hm2 = heatmap!(ax, 0.5:0.5:17.5,0.5:0.5:4.5, new_mat, colormap = co2, colorrange=extrema(new_mat))
    #Colorbar(f[1,2],hm2)
    translate!(hm2, 0,0,-100)
    
    #@show num_mat
    n_max = Int(maximum(num_mat))
    #@show n_max
    co = cgrad(:amp, n_max+1, categorical = true, alpha=0.55)
    #co = cgrad(:amp, n_max+1, categorical = true, alpha=[0.7,0.7,0.7,0.7,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9])
    # hm_n = heatmap!(ax2, 0.5:0.5:1,0.5:0.5:4.5, num_mat, colormap = co)
    hm_n = heatmap!(ax2, 0:0.5:1,0.5:0.5:4.5, num_mat, colormap = co, colorrange=(0,n_max))
    #hm_n = heatmap!(ax2, 0:0.5:1,0.5:0.5:4.5, num_mat, colormap =Makie.Categorical(:amp),colorrange=(0,n_max))
    Colorbar(gr[1,end+1], #hm_n,
        colormap = co, 
        colorrange=(0,n_max+1),
        ticks = (0.5:1:n_max+0.5, string.(0:1:n_max)), 
        #ticks = 0.5:1:n_max+1.5,
        label = "No. of significant EWS (p<0.05)",
        halign=:left)

    translate!(hm_n, 0,0,-100)


    Label(gl[1,1,TopLeft()], "(a)", fontsize = 20,
                        font = :bold, padding = (0,20,-10,15),
                        halign = :left,valign =:bottom)
    Label(gr[1,1,TopLeft()], "(b)", fontsize = 20,
                        font = :bold, padding = (0,20,-10,15),
                        halign = :left,valign =:bottom)
    elems2 = [
        [PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 0.5), (0.5, 0.5), (0.5, 0),(0,0)]),
        #MarkerElement(color = :black, marker = 'x', markersize = 15, points = Point2f[(0.25, 0.75)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 0.5), (1, 0.5), (1, 0),(0.5,0)]),
        PolyElement(color = :grey20, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 1), (0.5, 1), (0.5, 0.5),(0,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 1), (1, 1), (1, 0.5),(0.5,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey10, strokewidth = 2, points = Point2f[(0, 1), (1, 1), (1, 0),(0,0)])
        ],
        [PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 0.5), (0.5, 0.5), (0.5, 0),(0,0)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 0.5), (1, 0.5), (1, 0),(0.5,0)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 1), (0.5, 1), (0.5, 0.5),(0,0.5)]),
        PolyElement(color = :grey20, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 1), (1, 1), (1, 0.5),(0.5,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey10, strokewidth = 2, points = Point2f[(0, 1), (1, 1), (1, 0),(0,0)])
        ],
        [PolyElement(color = :grey20, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 0.5), (0.5, 0.5), (0.5, 0),(0,0)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 0.5), (1, 0.5), (1, 0),(0.5,0)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 1), (0.5, 1), (0.5, 0.5),(0,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 1), (1, 1), (1, 0.5),(0.5,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey10, strokewidth = 2, points = Point2f[(0, 1), (1, 1), (1, 0),(0,0)])
        ],
        [PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 0.5), (0.5, 0.5), (0.5, 0),(0,0)]),
        PolyElement(color = :grey20, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 0.5), (1, 0.5), (1, 0),(0.5,0)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 1), (0.5, 1), (0.5, 0.5),(0,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 1), (1, 1), (1, 0.5),(0.5,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey10, strokewidth = 2, points = Point2f[(0, 1), (1, 1), (1, 0),(0,0)])],
        [MarkerElement(color = co2[1], marker = :rect, markersize = 25,
        strokecolor = :grey, strokewidth = 0.7)],
        [MarkerElement(color =co2[2], marker = :rect, markersize = 25,
        strokecolor = :grey, strokewidth = 0.7)],
        [MarkerElement(color = co2[3], marker = :rect, markersize = 25,
        strokecolor = :grey, strokewidth = 0.7)],
        ]
    Legend(f[end+1,1], elems2, [ L"𝖵", L"\alpha_𝟣", L"\hat{𝗐}^𝟤_{𝟣𝟢,𝟧𝟢}", L"\hat{𝖧}^{\text{𝗅𝗈𝖼}}_{𝟣𝟢,𝟧𝟢}", "decreasing", "increasing","significantly increasing (p<$(plim))"],
                "EWS indicators",
                titleposition=:top,
                #titlefont=:regular,
                #titlegap=15,
                rowgap = 15, 
                framevisible = false, #false,
                framecolor = :grey70,
                tellwidth = false,
                tellheight=true,
                orientation = :horizontal,
                patchsize=(30,40),
                patchlabelgap=10,
                )
    
    colgap!(f.layout,20)
    if showing
        display(f)
    end
    if saving
        save(saveto, f)
    end
    
end


aggregated_plot2_method(lowpass = true, saving=true, saveto="do_ews_across_greenland_ice_cores/figures/aggregated_overview2_incl_n_method_lp.pdf")
aggregated_plot2_method(lowpass = false, saving=true, saveto="do_ews_across_greenland_ice_cores/figures/aggregated_overview2_incl_n_method_no_lp.pdf")


function aggregated_plot2_irreg(;lowpass= lowpass, showing = true, saving = false, saveto = "paper/result_overview_method_lowpass_$(lowpass)_smoothw_false_p_0.05.pdf")
    #cases = [L"\textbf{\mathrm{\sigma^𝟤}}",L"\mathrm{\alpha_𝟣}", L"$\mathrm{\hat{𝗐}^𝟤}$, 𝟤𝟢-𝟨𝟢 𝗒, 𝖯𝖺𝗎𝗅", L"$\mathrm{\hat{𝖧}_{𝟤𝟢-𝟨𝟢 𝗒, 𝖯𝖺𝗎𝗅}}$"]
    
    irreg_paths = ["new_surrogate_files/NGRIP_irreg/var_filt_1000_tfts.jld2",
                    "new_surrogate_files/NGRIP_irreg/ac_filt_1000_tfts.jld2",
                    "new_surrogate_files/NGRIP_irreg/sca_10_50_1000_tfts.jld2",
                    "new_surrogate_files/NGRIP_irreg/hurst_10_50_1000_tfts.jld2"]
    
    morlet_paths = ["new_surrogate_files/NGRIP5/var/w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/NGRIP5/ac/w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                     "new_surrogate_files/NGRIP5/sca/s1_10_s2_50_C_lowpass_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                     "new_surrogate_files/NGRIP5/hurst/s1_10_s2_50_C_lowpass_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]
                    
    paul_paths = ["new_surrogate_files/NGRIP5/var/w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                    "new_surrogate_files/NGRIP5/ac/w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                     "new_surrogate_files/NGRIP5/sca/s1_10_s2_50_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                     "new_surrogate_files/NGRIP5/hurst/s1_10_s2_50_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]
   #irregular_paths = []

    
    #Morlet_paths = []

    if lowpass == false
       
        paul_paths = ["new_surrogate_files/NGRIP5/var/w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                        "new_surrogate_files/NGRIP5/ac/w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                        "new_surrogate_files/NGRIP5/sca/s1_10_s2_50_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                        "new_surrogate_files/NGRIP5/hurst/s1_10_s2_50_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]

        morlet_paths = ["new_surrogate_files/NGRIP5/var/w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                        "new_surrogate_files/NGRIP5/ac/w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",
                        "new_surrogate_files/NGRIP5/sca/s1_10_s2_50_C_no_lowpass_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2",
                        "new_surrogate_files/NGRIP5/hurst/s1_10_s2_50_C_no_lowpass_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"]

        #irregular_paths = []


        #Morlet_paths = []
    end

    
    f= Figure(size=(1200,450))

    function get_n(x,plim)
        n=0
        for ev in 1:17
            if typeof(x.slopes[ev]) != Nothing
                if x.slopes[ev] >0  && x.p_one[ev] < plim
                    n+=1
                end
            end
        end
        return n
    end

    cases2 = [L"\mathrm{𝖵𝖺𝗋}",
              L"\mathrm{\alpha_𝟣}",
              L"$\mathrm{\hat{𝗐}^𝟤_{𝟤𝟢,𝟨𝟢}}$",# (𝖯𝖺𝗎𝗅)", 
              L"$\mathrm{\hat{𝖧}_{𝟤𝟢,𝟨𝟢}}$"]# (𝖯𝖺𝗎𝗅)"]
    
    gl = f[1,1] = GridLayout()
    ax = Axis(gl[1,1],
        xticks = (1:17, ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]) ,
        xreversed = true,
        yreversed=true,
        xlabel = "Event",
        xminorgridwidth = 3.5,
        xgridwidth=0.7,
        xminorgridcolor = :grey10 ,
        xgridcolor=:grey30,
        xminorgridvisible = true, 
        xgridvisible = true,
        yticks = (1:3, ["irregular\n (Morlet)", "5 years\n (Morlet)", "5 years\n (Paul)"]),
        ylabel = "Temporal reolution\n (wavelet basis)",
        yminorgridwidth = 4.5,
        ygridwidth=0.7,
        yminorgridcolor = :grey10 ,
        ygridcolor=:grey30,
        yminorgridvisible = true, 
        ygridvisible = true,
        )
    
    gr = f[1,2] = GridLayout(width=80)
    ax2 = Axis(gr[1,1],
        #xticks = (1:17, ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]) ,
        xticks=[0.5],
        xreversed = true,
        yreversed=true,
        #xlabel = "Transition",
        xminorgridwidth = 3.5,
        xgridwidth=0.7,
        xminorgridcolor = :grey10 ,
        xgridcolor=:grey30,
        xminorgridvisible = true, 
        xgridvisible = true,
        yticks = (1:3, ["irregular\n (Morlet)", "5 years\n (Morlet)", "5 years\n (Paul)"]),
        #ylabel = "Method",
        yminorgridwidth = 4.5,
        ygridwidth=0.7,
        yminorgridcolor = :grey10 ,
        ygridcolor=:grey30,
        yminorgridvisible = true, 
        ygridvisible = true,
        xticksvisible = false,
        xticklabelsvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false
        )

    
    new_mat = zeros(2*17,2*3)
    num_mat = zeros(2,2*3)
    for (i,core) in enumerate([irreg_paths, morlet_paths, paul_paths])
        for (j,path) in enumerate(core)     
            ind = load(path)["slopes"]
            n_h = get_n(ind,0.05)
            t_c =  :black #:grey60
            if n_h >= 4
                t_f = :bold
                f_s = 16
                #t_c =  :grey60
            else
                t_f = :regular
                f_s = 15
                #t_c =  :black
            end

            if j==1
                text!(ax2, 0.75,i-0.25, text = "$n_h",align = (:center, :center), color=t_c, font = t_f, fontsize = f_s)
                num_mat[2,2*i-1]=n_h   
            elseif j==2
                text!(ax2, 0.25,i-0.25, text = "$n_h",align = (:center, :center), color=t_c, font = t_f, fontsize = f_s)
                num_mat[1,2*i-1]=n_h
            elseif j==3
                text!(ax2, 0.75,i+0.25, text = "$n_h",align = (:center, :center), color=t_c, font = t_f, fontsize = f_s)
                num_mat[2,2*i]=n_h
            elseif j==4
                text!(ax2, 0.25,i+0.25, text = "$n_h",align = (:center, :center), color=t_c, font = t_f, fontsize = f_s)
                num_mat[1,2*i]=n_h
            end
            # if n_h >=4
            #     if j==1
            #         scatter!(ax2, [0.75],[i-0.25], strokecolor = :black, strokewidth=3, marker= :circle, markersize=20, color=:black)
                  
            #     elseif j==2
            #         text!(ax2, 0.25,i-0.25, text =string(n_h),align = (:center, :center), color=:grey70)
                 
            #     elseif j==3
            #         text!(ax2, 0.75,i+0.25, text = string(n_h),align = (:center, :center), color=:grey70)
            #     elseif j==4
            #         text!(ax2, 0.25,i+0.25, text = string(n_h),align = (:center, :center), color=:grey70)
    
            #     end
            # end
            for ev in 1:17
                # if j==1
                #     new_mat[2*ev,2*i-1]=1
                # elseif j==2
                #     new_mat[2*ev-1,2*i-1]=2
                # elseif j==3
                #     new_mat[2*ev,2*i]=3
                # elseif j==4
                #     new_mat[2*ev-1,2*i]=4
                # end
                if typeof(ind.slopes[ev]) != Nothing
                    if ind.slopes[ev]>0
                        if j==1
                            #text!(ax, ev+0.25,i-0.25, text = L"𝖵",align = (:center, :center))
                            new_mat[2*ev,2*i-1]+=1
                        elseif j==2
                            #text!(ax, ev-0.25,i-0.25, text = L"\alpha_𝟣",align = (:center, :center))
                            new_mat[2*ev-1,2*i-1]+=1
                        elseif j==3
                            #text!(ax, ev+0.25,i+0.25, text = L"\hat{𝗐}^𝟤",align = (:center, :center))
                            new_mat[2*ev,2*i]+=1
                        elseif j==4
                            #text!(ax, ev-0.25,i+0.25, text = L"\hat{𝖧}^{\text{𝗅𝗈𝖼}}",align = (:center, :center))
                            new_mat[2*ev-1,2*i]+=1
                        end
                        tc2 = :black
                        fs2 = 16
                        if ind.p_one[ev] < plim
                            if j==1
                                text!(ax, ev+0.25,i-0.25, text = L"\mathbf{𝖵}",align = (:center, :center), color=tc2, fontsize = fs2)
                                new_mat[2*ev,2*i-1]+=1
                            elseif j==2
                                if i == 1 
                                    text!(ax, ev-0.25,i-0.25, text = L"\mathbf{\hat{\alpha}_𝟣}",align = (:center, :center), color=tc2, fontsize = fs2)
                                else
                                    text!(ax, ev-0.25,i-0.25, text = L"\mathbf{\alpha_𝟣}",align = (:center, :center), color=tc2, fontsize = fs2)
                                end
                                new_mat[2*ev-1,2*i-1]+=1
                            elseif j==3
                                text!(ax, ev+0.25,i+0.25, text = L"\mathbf{\hat{𝗐}^𝟤}",align = (:center, :center), color=tc2, fontsize = fs2)
                                new_mat[2*ev,2*i]+=1
                            elseif j==4
                                text!(ax, ev-0.25,i+0.25, text = L"\mathbf{\hat{𝖧}^{\text{𝗅𝗈𝖼}}}",align = (:center, :center), color=tc2, fontsize = fs2)
                                new_mat[2*ev-1,2*i]+=1
                            end
                        end
                    end
                end
            end
            
        end
    end
    #t1 = cgrad(:managua, 4, categorical=true )
    #t2 = cgrad(:roma, 4, categorical=true )
    #t3 = cgrad(:Spectral_4, 4, categorical=true )
    #co2 = cgrad([[i for i in t3]...,[i for i in t3]...],8, categorical=true, alpha=[0.08,0.08,0.08,0.08,0.8,0.8,0.8,0.8])   
    
    #co2 = cgrad([:seagreen3, :dodgerblue3, :goldenrod1, :firebrick, :seagreen3, :dodgerblue3, :goldenrod1, :firebrick], 8, categorical = true, alpha=[0.1,0.1,0.1,0.1,0.6,0.6,0.6,0.6])#[0.1,0.1,0.1,0.1,1,1,1,1])    
    ##co2 = cgrad([[i for i in t1]...,[i for i in t1]...],8, categorical=true, alpha=[0.1,0.1,0.1,0.1,0.7,0.7,0.7,0.7])
    #co2 = cgrad([t2[3],t2[4],t2[1],t2[1],t2[3],t2[4],t2[1],t2[2]],8, categorical=true, alpha=[0.08,0.08,0.08,0.08,0.75,0.75,0.75,0.75])    
    #co2 = cgrad([[i for i in t3]...,[i for i in t3]...],8, categorical=true, alpha=[0.1,0.1,0.1,0.1,0.75,0.75,0.75,0.75])     
    #co2 = cgrad([:steelblue1, :steelblue2, :steelblue3, :steelblue, :firebrick1, :firebrick2, :firebrick3, :firebrick], 8, categorical = true, alpha=[0.2,0.2,0.2,0.2,0.8,0.8,0.8,0.8])
    # co2 = cgrad([:grey70, :grey70, :grey70, :grey70, :firebrick3, :firebrick3, :firebrick3, :firebrick3], 8, categorical = true, alpha=[0.1,0.1,0.1,0.1,0.8,0.8,0.8,0.8])
    co2 = cgrad([:steelblue, :darkred, :firebrick], 3, categorical = true, alpha=[0.12,0.12,0.8])
    hm2 = heatmap!(ax, 0.5:0.5:17.5,0.5:0.5:3.5, new_mat, colormap = co2, colorrange=extrema(new_mat))
    #Colorbar(f[1,2],hm2)
    translate!(hm2, 0,0,-100)
    
    #@show num_mat
    n_max = Int(maximum(num_mat))
    #@show n_max
    co = cgrad(:amp, n_max+1, categorical = true, alpha=0.55)
    # if lowpass
    #     co = cgrad(:amp, n_max+1, categorical = true, alpha=[0.7,0.7,0.7,0.7,0.9,0.9])
    # else
    #     co = cgrad(:amp, n_max+1, categorical = true, alpha=[0.7,0.7,0.7,0.7,0.9])
    # end
    # hm_n = heatmap!(ax2, 0.5:0.5:1,0.5:0.5:4.5, num_mat, colormap = co)
    hm_n = heatmap!(ax2, 0:0.5:1,0.5:0.5:3.5, num_mat, colormap = co, colorrange=(0,n_max))
    #hm_n = heatmap!(ax2, 0:0.5:1,0.5:0.5:4.5, num_mat, colormap =Makie.Categorical(:amp),colorrange=(0,n_max))
    Colorbar(gr[1,end+1], #hm_n,
        colormap = co, 
        colorrange=(0,n_max+1),
        ticks = (0.5:1:n_max+0.5, string.(0:1:n_max)), 
        #ticks = 0.5:1:n_max+1.5,
        label = "No. of significant EWS (p<0.05)",
        halign=:left)

    translate!(hm_n, 0,0,-100)


    Label(gl[1,1,TopLeft()], "(a)", fontsize = 20,
                        font = :bold, padding = (0,20,-10,15),
                        halign = :left,valign =:bottom)
    Label(gr[1,1,TopLeft()], "(b)", fontsize = 20,
                        font = :bold, padding = (0,20,-10,15),
                        halign = :left,valign =:bottom)
    elems2 = [
        [PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 0.5), (0.5, 0.5), (0.5, 0),(0,0)]),
        #MarkerElement(color = :black, marker = 'x', markersize = 15, points = Point2f[(0.25, 0.75)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 0.5), (1, 0.5), (1, 0),(0.5,0)]),
        PolyElement(color = :grey20, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 1), (0.5, 1), (0.5, 0.5),(0,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 1), (1, 1), (1, 0.5),(0.5,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey10, strokewidth = 2, points = Point2f[(0, 1), (1, 1), (1, 0),(0,0)])
        ],
        [PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 0.5), (0.5, 0.5), (0.5, 0),(0,0)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 0.5), (1, 0.5), (1, 0),(0.5,0)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 1), (0.5, 1), (0.5, 0.5),(0,0.5)]),
        PolyElement(color = :grey20, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 1), (1, 1), (1, 0.5),(0.5,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey10, strokewidth = 2, points = Point2f[(0, 1), (1, 1), (1, 0),(0,0)])
        ],
        [PolyElement(color = :grey20, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 0.5), (0.5, 0.5), (0.5, 0),(0,0)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 0.5), (1, 0.5), (1, 0),(0.5,0)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 1), (0.5, 1), (0.5, 0.5),(0,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 1), (1, 1), (1, 0.5),(0.5,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey10, strokewidth = 2, points = Point2f[(0, 1), (1, 1), (1, 0),(0,0)])
        ],
        [PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 0.5), (0.5, 0.5), (0.5, 0),(0,0)]),
        PolyElement(color = :grey20, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 0.5), (1, 0.5), (1, 0),(0.5,0)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0, 1), (0.5, 1), (0.5, 0.5),(0,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey30, strokewidth = 1, points = Point2f[(0.5, 1), (1, 1), (1, 0.5),(0.5,0.5)]),
        PolyElement(color = :transparent, strokecolor = :grey10, strokewidth = 2, points = Point2f[(0, 1), (1, 1), (1, 0),(0,0)])],
        [MarkerElement(color = co2[1], marker = :rect, markersize = 25,
        strokecolor = :grey, strokewidth = 0.7)],
        [MarkerElement(color =co2[2], marker = :rect, markersize = 25,
        strokecolor = :grey, strokewidth = 0.7)],
        [MarkerElement(color = co2[3], marker = :rect, markersize = 25,
        strokecolor = :grey, strokewidth = 0.7)],
        ]
    Legend(f[end+1,1], elems2, [ L"𝖵", L"\alpha_𝟣", L"\hat{𝗐}^𝟤_{𝟣𝟢,𝟧𝟢}", L"\hat{𝖧}^{\text{𝗅𝗈𝖼}}_{𝟣𝟢,𝟧𝟢}", "decreasing", "increasing","significantly increasing (p<$(plim))"],
                "EWS indicators",
                titleposition=:top,
                #titlefont=:regular,
                #titlegap=15,
                rowgap = 15, 
                framevisible = false, #false,
                framecolor = :grey70,
                tellwidth = false,
                tellheight=true,
                orientation = :horizontal,
                patchsize=(30,40),
                patchlabelgap=10,
                )
    
    colgap!(f.layout,20)
    if showing
        display(f)
    end
    if saving
        save(saveto, f)
    end
    
end

aggregated_plot2_irreg(lowpass = true, saving=true, saveto="do_ews_across_greenland_ice_cores/figures/aggregated_overview2_incl_n_irreg_lp.pdf")
aggregated_plot2_irreg(lowpass = false, saving=true, saveto="do_ews_across_greenland_ice_cores/figures/aggregated_overview2_incl_n_irreg_no_lp.pdf")

# elems2 = [[MarkerElement(color = t1[ci], marker = :rect, markersize = 25,
#                 strokecolor = :grey, strokewidth = 0.7)] for ci = 1:2]
# cgrad([:seagreen3, :darkcyan, :goldenrod1, :firebrick, :seagreen3, :darkcyan, :goldenrod1, :firebrick], 8, categorical = true, alpha=[0.1,0.1,0.1,0.1,1,1,1,1])
# t1 = cgrad(:managua, 4, categorical=true )
# t2 = cgrad(:roma, 4, categorical=true )
# t3 = cgrad(:Spectral_4, 4, categorical=true )
# t3
# cgrad(:balance, 4, categorical=true )
# cgrad([[i for i in t1]...,[i for i in t1]...],8, categorical=true, alpha=[0.1,0.1,0.1,0.1,1,1,1,1])
# cgrad(:davos, 4, categorical=true )

# cgrad([:steelblue1, :steelblue2, :steelblue3, :steelblue, :firebrick1, :firebrick2, :firebrick3, :firebrick], 8, categorical = true, alpha=[0.2,0.2,0.2,0.2,1,1,1,1])

1049-976
#################################
###########################

function load_data(type)
    ice_core_paths = readdir("new_surrogate_files/ice_cores/$(type)/", join = true)
    ice_core_names = readdir("new_surrogate_files/ice_cores/$(type)/", join = false)
    ice_core_names = [k for k in ice_core_names if occursin("jld2", k)]

    var_paths = readdir("new_surrogate_files/$(type)/var/", join = true)
    var_names = readdir("new_surrogate_files/$(type)/var/", join = false)
    var_names = [k for k in var_names if occursin("jld2", k)]

    ac_paths = readdir("new_surrogate_files/$(type)/ac/", join = true)
    ac_names = readdir("new_surrogate_files/$(type)/ac/", join = false)
    ac_names = [k for k in ac_names if occursin("jld2", k)]

    sca_paths = readdir("new_surrogate_files/$(type)/sca/", join = true)
    sca_names = readdir("new_surrogate_files/$(type)/sca/", join = false)
    sca_names = [k for k in sca_names if occursin("jld2", k)]

    hurst_paths = readdir("new_surrogate_files/$(type)/hurst/", join = true)
    hurst_names = readdir("new_surrogate_files/$(type)/hurst/", join = false)
    hurst_names = [k for k in hurst_names if occursin("jld2", k)]

    ice_cores = [load(k)["ice"] for k in ice_core_paths if occursin("jld2", k)]
    vars = [load(k)["slopes"] for k in var_paths if occursin("jld2", k)]
    acs = [load(k)["slopes"] for k in ac_paths if occursin("jld2", k)]
    scas = [load(k)["slopes"] for k in sca_paths if occursin("jld2", k)]
    hursts = [load(k)["slopes"] for k in hurst_paths if occursin("jld2", k)]

    all_sranges = Tuple{Int64, Int64}[]
    for i in 10:10:100
        for j in 20:10:110
            if i<j
                push!(all_sranges,(i,j))
            end
        end
    end
    all_sranges

    sr10_id = [s1>=10 && s2-s1 >=10 && s1% 10 == 0 && s2% 10 == 0  for (s1,s2) in all_sranges]
    sranges_10 = all_sranges[sr10_id]

    sr20_id = [s1>=20 && s2-s1 >=20 && s1% 20 == 0 && s2% 20 == 0 for (s1,s2) in all_sranges]
    sranges_20 = all_sranges[sr20_id]

    if type == "NGRIP5"
        sranges = all_sranges
    elseif type == "10y"
        sranges = sranges_10
    elseif type == "20y"
        sranges = sranges_20
    end
    return ice_cores, var_names, vars, ac_names, acs, sca_names, scas, hurst_names, hursts, sranges
end


function plot_ts(ice_cores; showing = true, saving = false, saveto= "test.pdf")
    f = Figure(size=(900,170*length(ice_cores)))
    ga = f[1,1] = GridLayout()
    for (i,ice) in enumerate(ice_cores)
        if i == length(ice_cores)
            ax = Axis(ga[i, 1], xlabel = "Age (kyr b2k)", ylabel = L"$\delta^{18}$O (‰)")
        else
            ax = Axis(ga[i, 1], ylabel = L"$\delta^{18}$O (‰)")
            hidexdecorations!(ax, grid=false)
        end
        ax.xticks = (10_000:5_000:60_000, string.(10:5:60))
        ax.xreversed = true

        for c in ice.cold_idx
            lines!(ax, ice.age[c], ice.δ[c], color = :darkblue)
        end
        for w in ice.warm_idx
            lines!(ax, ice.age[w], ice.δ[w], color = :darkred)
        end
        vlines!(ax, GS_onsets, color = :darkblue, alpha = 0.4)
        vlines!(ax, GI_onsets, color = :darkred, alpha = 0.4)

        Box(ga[i,1,Right()], color = :gray90)
        Label(ga[i,1,Right()], ice.name, rotation = pi/2)

    end
    linkyaxes!(f.content[1:3:end]...)
    linkxaxes!(f.content[1:3:end]...)
    if saving
        save(saveto, f)
    end
    if showing
        display(f)
    end
end


#Fig. 2 and S14
function plot_all_ts_with_labels(;do_labels = true, lowpass = lowpass, showing = true, saving = false, saveto= "paper/all_ts.pdf")
    dim2 = 1200
    f_all_ts = Figure(size=(1200,dim2))
    function good_files(v,type,lowpass)
        if type == "NGRIP5"
            if lowpass
                return occursin("C_lowpass",v)
            else
                return occursin("C_no_lowpass",v)
            end
        elseif type == "10y"
            if lowpass
                return occursin("_lp",v)
            else
                return !occursin("_lp",v)
            end
        elseif type == "20y"
            return true
        end
    end
    
    whichletter = 1
    for (it,type) in enumerate(["NGRIP5", "10y", "20y"])
        ice_core_paths = readdir("new_surrogate_files/ice_cores/$(type)/", join = true)
        ice_core_names = readdir("new_surrogate_files/ice_cores/$(type)/", join = false)
        ice_core_names = [k for k in ice_core_names if occursin("jld2", k)]
        ice_core_paths = [k for k in ice_core_paths if occursin("jld2", k)]

        good = @. good_files(ice_core_names,type,lowpass)
        ices = [load(k)["ice"] for k in ice_core_paths[good]]
        
        for (i,ice) in enumerate(ices[end:-1:1])
            gl = f_all_ts[whichletter,1] = GridLayout()
            ax = Axis(gl[1,1], 
                    ylabel = L"$\delta^{𝟣𝟪}$𝖮 (‰)",
                    ylabelsize = 15,
                    ylabelfont = :bold,)
            
            ax.xticks = (10_000:5_000:60_000, string.(10:5:60))
            ax.xreversed = true
            xlims!(ax,ice.age[1], ice.age[end])
            
            #extra axis for the DO labels
            if do_labels
                ax_event = Axis(gl[1,1],xaxisposition = :top, xticklabelrotation = pi/5)
                ax_event.xticks = (GI_onsets, event_labels)
                ax_event.xreversed = true
                xlims!(ax_event,ice.age[1], ice.age[end])
                hidespines!(ax_event)
                hideydecorations!(ax_event)
            end
            
            if whichletter ∈ [1,3,6]
                ax.xlabel = "Age (kyr b2k)"
            else
                hidexdecorations!(ax, grid=false)
                hidespines!(ax, :b)

            end

            if i == 1
                if do_labels
                    ax_event.title = "Ice core records with $([5,10,20][it])-year resolution"
                    ax_event.titlesize = 20
                else
                    ax.title = "Ice core records with $([5,10,20][it])-year resolution"
                    ax.titlesize = 20
                end
            elseif i>1
                hidespines!(ax, :t)
                if do_labels
                    hidexdecorations!(ax_event, grid=false)#, ticks = false)
                end
            end
            for c in ice.cold_idx
                lines!(ax, ice.age[c], ice.δ[c], color = :darkblue)
            end
            for w in ice.warm_idx
                lines!(ax, ice.age[w], ice.δ[w], color = :darkred)
            end
            vlines!(ax, GS_onsets, color = :darkblue, alpha = 0.4)
            vlines!(ax, GI_onsets, color = :darkred, alpha = 0.4)
            nam = ice.name
            if it == 1
                nam = "NGRIP"
            elseif it == 2 && lowpass
                nam = ice.name[1:end-3]
            end
            Label(gl[1,1,Left()], nam, padding = (0.0,70.0,0.0,0.0), valign = :center, font = :bold, fontsize = 15, rotation = pi/2)
            Label(gl[1,1,TopLeft()], letters[whichletter], fontsize = 20,
                            font = :bold, padding = (0,5,5,0),
                            halign = :left,valign=:bottom)
            whichletter +=1
        end
    end
    if do_labels
        linkyaxes!(f_all_ts.content[1:4:end][2:3]...)
        linkyaxes!(f_all_ts.content[1:4:end][4:6]...)
        
        linkxaxes!(f_all_ts.content[1:4:end]...)
    else
        linkyaxes!(f_all_ts.content[1:3:end][2:3]...)
        linkyaxes!(f_all_ts.content[1:3:end][4:6]...)
        
        linkxaxes!(f_all_ts.content[1:3:end]...)
    end

    if do_labels
        Label(f_all_ts[2,1:end,Top()], " ", padding = (0.0,0.0,30.0,30.0), valign = :bottom, font = :bold, fontsize = 45)
        Label(f_all_ts[4,1:end,Top()], " ", padding = (0.0,0.0,30.0,30.0), valign = :bottom, font = :bold, fontsize = 45)
    else
        Label(f_all_ts[2,1:end,Top()], " ", padding = (0.0,0.0,30.0,30.0), valign = :bottom, font = :bold, fontsize = 10)
        Label(f_all_ts[4,1:end,Top()], " ", padding = (0.0,0.0,30.0,30.0), valign = :bottom, font = :bold, fontsize = 10)
    end
    
    rowgap!(f_all_ts.layout, -10)
    if saving
        save(saveto,f_all_ts)
    end
    if showing
        display(f_all_ts)
    end

end

plot_all_ts_with_labels(lowpass = true, showing = showing, saving = saving, saveto= "figures/fig2.pdf")
plot_all_ts_with_labels(lowpass = false, showing = showing, saving = saving, saveto= "figures/figS14.pdf")


function get_num_incr_csd(vars, acs, var_names, type; plim = 0.05)
    incr_matrix_var = zeros(Union{Missing, Float64},length(vars), 17)
    incr_matrix_ac = zeros(Union{Missing, Float64},length(vars), 17)
    sig_incr_matr_var_one = zeros(Union{Missing, Float64},length(vars), 17)
    sig_incr_matr_var_two = zeros(Union{Missing, Float64},length(vars), 17)
    sig_incr_matr_ac_one = zeros(Union{Missing, Float64},length(vars), 17)
    sig_incr_matr_ac_two = zeros(Union{Missing, Float64},length(vars), 17)
    sig_incr_matr_both_one = zeros(Union{Missing, Float64},length(vars), 17)
    sig_incr_matr_both_two = zeros(Union{Missing, Float64},length(vars), 17)
    
    ylabs = []
    for (i,var) in enumerate(vars)
        if type != "NGRIP5"
            push!(ylabs, split(var_names[i],"_")[5])
        else
            push!(ylabs, "NGRIP")
        end
        exist_v = findall(x -> !isnothing(x), var.slopes)
        notexists_v = findall(x -> isnothing(x), var.slopes)

        incr_matrix_var[i,exist_v[findall(x->x>0, var.slopes[exist_v])]] .=1
        incr_matrix_var[i, notexists_v] .= missing
        ews_v1 = intersect(findall(x->x>0, var.slopes[exist_v]), findall(x-> x < plim, var.p_one[exist_v]))
        ews_v2 = intersect(findall(x->x>0, var.slopes[exist_v]), findall(x-> x < plim, var.p_two[exist_v]))
        sig_incr_matr_var_one[i,exist_v[ews_v1]] .= 1
        sig_incr_matr_var_two[i,exist_v[ews_v2]] .=1
        sig_incr_matr_var_one[i, notexists_v] .= missing
        sig_incr_matr_var_two[i, notexists_v] .= missing

        exists_a = findall(x -> !isnothing(x), acs[i].slopes)
        notexists_a = findall(x -> isnothing(x), acs[i].slopes)
        incr_matrix_ac[i,exists_a[findall(x->x>0, acs[i].slopes[exists_a])]] .= 1
        ews_a1 = intersect(findall(x->x>0, acs[i].slopes[exists_a]), findall(x-> x < plim, acs[i].p_one[exists_a]))
        ews_a2 = intersect(findall(x->x>0, acs[i].slopes[exists_a]), findall(x-> x < plim, acs[i].p_two[exists_a]))
        sig_incr_matr_ac_one[i,exists_a[ews_a1]] .= 1
        sig_incr_matr_ac_two[i,exists_a[ews_a2]] .= 1
        incr_matrix_ac[i, notexists_a] .= missing
        sig_incr_matr_ac_one[i, notexists_a] .= missing
        sig_incr_matr_ac_two[i, notexists_a] .= missing

        exist_both = intersect(exists_a, exist_v)
        notexists_both = union(notexists_a, notexists_v)
        
        ews_both1 = intersect(exist_both, ews_v1, ews_a1)
        ews_both2 = intersect(exist_both, ews_v2, ews_a2)

        sig_incr_matr_both_one[i,exist_both[ews_both1]] .= 1
        sig_incr_matr_both_two[i,exist_both[ews_both2]] .= 1
        sig_incr_matr_both_one[i,notexists_both] .= missing
        sig_incr_matr_both_two[i,notexists_both] .= missing

    end
    incr_matrix_both = incr_matrix_ac .* incr_matrix_var
    
    return incr_matrix_var ,
        incr_matrix_ac,
        incr_matrix_both,
        sig_incr_matr_var_one,
        sig_incr_matr_var_two,
        sig_incr_matr_ac_one,
        sig_incr_matr_ac_two,
        sig_incr_matr_both_one,
        sig_incr_matr_both_two, 
        ylabs 
end


function get_num_incr_wavelet(scana, ws, hurstna, hs, type, sranges; plim = 0.05, given_cnames = true, cn = ["C"])
    if type == "NGRIP5"
        ncores = 1
        cnames = ["C"]
    elseif type == "10y"
        ncores = 2
        cnames = ["NEEM", "NGRIP"]
    elseif  type == "20y"
        ncores = 3
        cnames = ["GISP2", "GRIP", "NGRIP"]
    end

    if given_cnames
        cnames = cnames
    else
        cnames = cn
    end

    incr_matrix_sca = Array{Union{Missing, Float64}}(missing,ncores,length(sranges), 17)
    sig_incr_matr_sca_one = Array{Union{Missing, Float64}}(missing,ncores,length(sranges), 17)
    sig_incr_matr_sca_two = Array{Union{Missing, Float64}}(missing,ncores,length(sranges), 17)
    
    incr_matrix_hurst = Array{Union{Missing, Float64}}(missing,ncores,length(sranges), 17)
    sig_incr_matr_hurst_one = Array{Union{Missing, Float64}}(missing,ncores,length(sranges), 17)
    sig_incr_matr_hurst_two = Array{Union{Missing, Float64}}(missing,ncores,length(sranges), 17)
    
    for (ii,sn) in enumerate(scana)
        splits = split(sn,"_")
        sh1 = parse(Int64, splits[2])
        sh2 = parse(Int64, splits[4])
        
        name = splits[5]

        for (kk, cname) in enumerate(cnames)
            for (jj,(s1,s2)) in enumerate(sranges)
                if name == cname && sh1 == s1 && sh2 == s2
                    sca = ws[ii]
                    exists = findall(x -> !isnothing(x), sca.slopes)
                
                    incr_matrix_sca[kk,jj,exists] .=0
                    incr_matrix_sca[kk,jj,exists[findall(x->x>0, sca.slopes[exists])]] .=1

                    sig_incr_matr_sca_one[kk,jj,exists] .=0
                    sig_incr_matr_sca_two[kk,jj,exists] .=0
        
                    ews_v1 = intersect(findall(x->x>0, sca.slopes[exists]), findall(x-> x < plim, sca.p_one[exists]))
                    ews_v2 = intersect(findall(x->x>0, sca.slopes[exists]), findall(x-> x < plim, sca.p_two[exists]))
                    sig_incr_matr_sca_one[kk,jj,exists[ews_v1]] .= 1
                    sig_incr_matr_sca_two[kk,jj,exists[ews_v2]] .=1
                end
            end
        end
        
    end


    for (ii,sn) in enumerate(hurstna)
        splits = split(sn,"_")
        sh1 = parse(Int64, splits[2])
        sh2 = parse(Int64, splits[4])
        
        name = splits[5]
       
        for (kk, cname) in enumerate(cnames)
            for (jj,(s1,s2)) in enumerate(sranges)
                if name == cname && sh1 == s1 && sh2 == s2
                    hurst = hs[ii]
                    exists = findall(x -> !isnothing(x), hurst.slopes)

                    incr_matrix_hurst[kk,jj,exists] .=0
                    incr_matrix_hurst[kk,jj,exists[findall(x->x>0, hurst.slopes[exists])]] .=1

                    sig_incr_matr_hurst_one[kk,jj,exists] .=0
                    sig_incr_matr_hurst_two[kk,jj,exists] .=0
                    
                    ews_v1 = intersect(findall(x->x>0, hurst.slopes[exists]), findall(x-> x < plim, hurst.p_one[exists]))
                    ews_v2 = intersect(findall(x->x>0, hurst.slopes[exists]), findall(x-> x < plim, hurst.p_two[exists]))
                    sig_incr_matr_hurst_one[kk,jj,exists[ews_v1]] .= 1
                    sig_incr_matr_hurst_two[kk,jj,exists[ews_v2]] .=1


                    # incr_matrix_hurst[kk,jj,exists[findall(x->x>0, hurst.slopes[exists])]] .=1
                    # incr_matrix_hurst[kk,jj,exists[findall(x->x≤ 0, hurst.slopes[exists])]] .=0
                    
                    # ews_v1 = intersect(findall(x->x>0, hurst.slopes[exists]), findall(x-> x < plim, hurst.p_one[exists]))
                    # ews_v2 = intersect(findall(x->x>0, hurst.slopes[exists]), findall(x-> x < plim, hurst.p_two[exists]))
                    # sig_incr_matr_hurst_one[kk,jj,exists[ews_v1]] .= 1
                    # sig_incr_matr_hurst_two[kk,jj,exists[ews_v2]] .=1
                    # sig_incr_matr_hurst_one[kk,jj, exists[findall(x->x ≤ 0, hurst.slopes[exists])]] .= 0
                    # sig_incr_matr_hurst_two[kk,jj, exists[findall(x->x ≤ 0, hurst.slopes[exists])]] .= 0
                end
            end
        end
        
    end

    # incr_matrix_both_wave = Array{Union{Missing, Float64}}(missing,size(incr_matrix_sca))
    # sig_incr_matr_both_wave_one = Array{Union{Missing, Float64}}(missing,size(incr_matrix_sca))
    # sig_incr_matr_both_wave_two = Array{Union{Missing, Float64}}(missing,size(incr_matrix_sca))



    # @show exists_w = findall(x -> !ismissing(x), incr_matrix_sca)
    # @show exists_h = findall(x -> !ismissing(x), incr_matrix_hurst)
    # @show ex_both = intersect(exists_w, exists_h)
    

    # both_incr = intersect(findall(x->x==(1), incr_matrix_sca[ex_both]), findall(x->x==(1), incr_matrix_hurst[ex_both]))
    # both_sig_incr_one = intersect(findall(x->x==(1), sig_incr_matr_sca_one[ex_both]), findall(x->x==(1), sig_incr_matr_hurst_one[ex_both]))
    # both_sig_incr_two = intersect(findall(x->x==(1), sig_incr_matr_sca_two[ex_both]), findall(x->x==(1), sig_incr_matr_hurst_two[ex_both]))

    # not_both_incr = intersect(findall(x->x==(0), incr_matrix_sca[ex_both]), findall(x->x==(0), incr_matrix_hurst[ex_both]))
    # not_both_sig_incr_one = intersect(findall(x->x==(0), sig_incr_matr_sca_one[ex_both]), findall(x->x==(0), sig_incr_matr_hurst_one[ex_both]))
    # not_both_sig_incr_two = intersect(findall(x->x==(0), sig_incr_matr_sca_two[ex_both]), findall(x->x==(0), sig_incr_matr_hurst_two[ex_both]))

    # incr_matrix_both_wave[ex_both[both_incr]] .=1
    # sig_incr_matr_both_wave_one[ex_both[both_sig_incr_one]] .=1
    # sig_incr_matr_both_wave_two[ex_both[both_sig_incr_two]] .=1

    # incr_matrix_wave_both[ex_both[not_both_incr]] .=0
    # sig_incr_matr_both_wave_one[ex_both[not_both_sig_incr_one]] .=0
    # sig_incr_matr_both_wave_two[ex_both[not_both_sig_incr_two]] .=0

    incr_matrix_both_wave = incr_matrix_sca .* incr_matrix_hurst
    sig_incr_matr_both_wave_one = sig_incr_matr_sca_one .* sig_incr_matr_hurst_one
    sig_incr_matr_both_wave_two = sig_incr_matr_sca_two .* sig_incr_matr_hurst_two

    return incr_matrix_sca, 
        sig_incr_matr_sca_one, 
        sig_incr_matr_sca_two, 
        incr_matrix_hurst, 
        sig_incr_matr_hurst_one, 
        sig_incr_matr_hurst_two, 
        incr_matrix_both_wave,
        sig_incr_matr_both_wave_one,
        sig_incr_matr_both_wave_two
end



function get_incr_matrices(;lowpass = true, smoothw =false, plim = 0.05)
    incr_matrix_var = Array{Matrix{Union{Missing, Float64}}}(undef,3)
    incr_matrix_ac = Array{Matrix{Union{Missing, Float64}}}(undef,3)
    incr_matrix_both = Array{Matrix{Union{Missing, Float64}}}(undef,3)
    sig_incr_matr_var_one = Array{Matrix{Union{Missing, Float64}}}(undef,3)
    sig_incr_matr_var_two = Array{Matrix{Union{Missing, Float64}}}(undef,3)
    sig_incr_matr_ac_one = Array{Matrix{Union{Missing, Float64}}}(undef,3)
    sig_incr_matr_ac_two = Array{Matrix{Union{Missing, Float64}}}(undef,3)
    sig_incr_matr_both_one = Array{Matrix{Union{Missing, Float64}}}(undef,3)
    sig_incr_matr_both_two =Array{Matrix{Union{Missing, Float64}}}(undef,3)
    ylabs = Array{Array{String}}(undef,3)

    incr_matrix_sca = Array{Array{Union{Missing, Float64}}}(undef,3)
    sig_incr_matr_sca_one = Array{Array{Union{Missing, Float64}}}(undef,3)
    sig_incr_matr_sca_two = Array{Array{Union{Missing, Float64}}}(undef,3)
    incr_matrix_hurst = Array{Array{Union{Missing, Float64}}}(undef,3)
    sig_incr_matr_hurst_one = Array{Array{Union{Missing, Float64}}}(undef,3)
    sig_incr_matr_hurst_two = Array{Array{Union{Missing, Float64}}}(undef,3)
    incr_matrix_wave_both = Array{Array{Union{Missing, Float64}}}(undef,3)
    sig_incr_matr_wave_both_one = Array{Array{Union{Missing, Float64}}}(undef,3)
    sig_incr_matr_wave_both_two = Array{Array{Union{Missing, Float64}}}(undef,3)

    scaleranges = Array{Array{Tuple{Int64, Int64}}}(undef, 3)


    function good_csd_files(v,type,lowpass)
        if type == "NGRIP5"
            if lowpass
                return occursin("w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",v) && !occursin("NIKLAS_", v)
            else
                return occursin("w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",v) && !occursin("NIKLAS_", v)
            end
        elseif type == "10y"
            if lowpass
                return occursin("w_200_normed_filt", v) && occursin("lp_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2", v)
            else
                return occursin("w_200_normed_filt", v) && occursin("gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2", v) && !occursin("_lp_", v)
            end
        elseif type == "20y"
            return occursin("w_200_normed_filt", v) && occursin("gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2", v)
        end
    end

    function good_wavelet_files(v,type,lowpass, smoothw)
        if type == "NGRIP5"
            if lowpass
                return occursin("C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_$(smoothw)_10000_TFTS.jld2",v) && !occursin("NIKLAS_", v)
            else
                return occursin("C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_$(smoothw)_10000_TFTS.jld2",v) && !occursin("NIKLAS_", v)
            end
        elseif type == "10y"
            if lowpass
                return occursin("lp_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_$(smoothw)_10000_TFTS.jld2", v)
            else
                return occursin("PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_$(smoothw)_10000_TFTS.jld2", v) && !occursin("_lp_", v)
            end
        elseif type == "20y"
            return occursin("PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_$(smoothw)_10000_TFTS.jld2", v)
        end
    end

    for (it,type) in enumerate(["NGRIP5", "10y", "20y"])

        ice_cores, var_names, vars, ac_names, acs, sca_names, scas, hurst_names, hursts, sranges = load_data(type)
        scaleranges[it] = sranges
        
        good = @. good_csd_files(var_names,type,lowpass)
        if type == "NGRIP5"
            vnames = ["NGRIP"]
        else
            vnames = var_names[good]
        end
        v = vars[good]
        a = acs[good]
        
        incr_matrix_var[it] ,
        incr_matrix_ac[it],
        incr_matrix_both[it],
        sig_incr_matr_var_one[it],
        sig_incr_matr_var_two[it],
        sig_incr_matr_ac_one[it],
        sig_incr_matr_ac_two[it],
        sig_incr_matr_both_one[it],
        sig_incr_matr_both_two[it], 
        ylabs[it] = get_num_incr_csd(v, a, vnames, type, plim = plim)


        good_sca = @. good_wavelet_files(sca_names,type,lowpass, smoothw)
        good_hurst = @. good_wavelet_files(hurst_names,type,lowpass, smoothw)


        scana = sca_names[good_sca] 
        hurstna = hurst_names[good_hurst] 

        ws = scas[good_sca] 
        hs = hursts[good_hurst] 
        
        incr_matrix_sca[it], 
        sig_incr_matr_sca_one[it], 
        sig_incr_matr_sca_two[it], 
        incr_matrix_hurst[it], 
        sig_incr_matr_hurst_one[it], 
        sig_incr_matr_hurst_two[it], 
        incr_matrix_wave_both[it],
        sig_incr_matr_wave_both_one[it],
        sig_incr_matr_wave_both_two[it] = get_num_incr_wavelet(scana, ws, hurstna, hs, type, sranges; plim = plim)
        
    end
    return  incr_matrix_var,
    incr_matrix_ac,
    incr_matrix_both,
    sig_incr_matr_var_one,
    sig_incr_matr_var_two,
    sig_incr_matr_ac_one,
    sig_incr_matr_ac_two,
    sig_incr_matr_both_one,
    sig_incr_matr_both_two,
    ylabs,
    incr_matrix_sca,
    sig_incr_matr_sca_one,
    sig_incr_matr_sca_two,
    incr_matrix_hurst,
    sig_incr_matr_hurst_one,
    sig_incr_matr_hurst_two,
    incr_matrix_wave_both,
    sig_incr_matr_wave_both_one,
    sig_incr_matr_wave_both_two,
    scaleranges

end


# Fig A5 and S23
function plot_num_incr_wavelet_own_scales2(sig_incr_matr_sca_one,sig_incr_matr_hurst_one,sig_incr_matr_wave_both_one,
    scaleranges, ylabs;showing = true, saving = false, saveto = "paper/num_incr_wavelet_own_scales_lowpass_$(lowpass)_smoothw_$(smoothw)_p_$(plim).pdf")

    letters_h = reshape(letters[1:18],3,6)
    #all num of sig increases in own scales
    ff = Figure(size=(1500,1600))
    siggi = @. replace(sig_incr_matr_sca_one, missing => 0)
    allma_sca = @. maximum(sum(siggi, dims=3),dims=2)
    max_sca = maximum(@. maximum(allma_sca))
    n_max = Int(max_sca)
    co = cgrad(:amp,  n_max+1, rev = false, categorical = true)

    for j = 1:3  
        if j == 1
            sig_incr_matr = sig_incr_matr_sca_one
            ind_n = L"\textbf{\mathrm{\hat{𝗐}^𝟤}}"
        elseif j == 2
            sig_incr_matr = sig_incr_matr_hurst_one
            ind_n =L"\textbf{\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}}}"
        elseif j == 3
            sig_incr_matr = sig_incr_matr_wave_both_one
            ind_n = L"\textbf{\mathrm{\hat{𝗐}^𝟤} \text{𝖺𝗇𝖽 }\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}}}"
        end
        whichletter = 1
        for i = 1:3
            s1s = sort(unique([s[1] for s in scaleranges[i]])) 
            s2s = sort(unique([s[2] for s in scaleranges[i]]))
            replace!(sig_incr_matr[i], missing => 0)
            n_incr = sum(sig_incr_matr[i], dims = 3)
            for k =  size(n_incr)[1]:-1:1
                ga = ff[whichletter,j] = GridLayout()
                ax =  Axis(ga[1,1], 
                            title = ind_n,
                            titlesize=20,
                            titlefont = :bold,
                            xticks = (1:length(s1s), string.(s1s)),
                            #xreversed = true,
                            xlabel = "s₁ (years)",
                            yticksvisible = true,
                            xminorgridwidth = 1.0, 
                            xminorgridcolor = :grey30 ,
                            xminorgridvisible = true, xgridvisible = false,
                            yticks = (1:length(s2s), string.(s2s)), 
                            #ylabel = "(s_1, s_2)",
                            ylabel = "s₂ (years)",
                            #yticklabelrotation = pi/2,
                            #yminorticks = 1.5:1:length(srange),
                            yminorgridwidth = 1.0, 
                            yminorgridcolor = :grey30, #:grey5,
                            yminorgridvisible = true, ygridvisible = false)
                n_incr_matrix = Array{Union{Missing, Float64}}(missing, length(s1s), length(s2s))
                for (isr, sr) in enumerate(scaleranges[i])
                    n_incr_matrix[findfirst(==(sr[1]), s1s), findfirst(==(sr[2]), s2s)] =  n_incr[k,isr,1]
                end

                rowsize!(ga,1, Auto(1))

                hm = heatmap!(ax, n_incr_matrix, colormap = co, colorrange = (0,n_max))
                translate!(hm, 0,0,-100)
                Colorbar(ga[1,2], colormap = co, colorrange = (0,n_max+1),
                ticks = (0.5:1:Int(n_max)+0.5, string.(0:1:Int(n_max))),
                label = "No. of significant EWS")

                Label(ga[1,1,TopLeft()], letters_h[j,whichletter], fontsize = 20,
                    font = :bold,  padding = (0,5,5,0),
                    halign = :left)

                if j == 1
                    Label(ga[1,1,Left()], ylabs[i][k], padding = (0.0,70.0,0.0,0.0), valign = :center, font = :bold, fontsize = 20, rotation = pi/2)
                end

                whichletter +=1
            end
        end      
    end


    Label(ff[1,1:end,Top()], "Number of significant EWS in (s₁-s₂) year bands of records with 5-year resolution", padding = (0.0,0.0,30.0,0.0), valign = :bottom, font = :bold, fontsize = 20)
    Label(ff[2,1:end,Top()], "Number of significant EWS in (s₁-s₂) year bands of records with 10-year resolution", padding = (0.0,0.0,30.0,3.0), valign = :bottom, font = :bold, fontsize = 20)
    Label(ff[4,1:end,Top()], "Number of significant EWS in (s₁-s₂) year bands of records with 20-year resolution", padding = (0.0,0.0,30.0,3.0), valign = :bottom, font = :bold, fontsize = 20)

    colgap!(ff.layout, 80)
    if saving
        save(saveto, ff)
    end
    if showing
        display(ff)
    end
end


# Fig S10
function plot_num_incr_wavelet_20_scales(sig_incr_matr_sca_one,sig_incr_matr_hurst_one,sig_incr_matr_wave_both_one,
    scaleranges, ylabs;showing = true, saving = false, saveto = "paper/num_incr_wavelet_20_scales_lowpass_$(lowpass)_smoothw_$(smoothw)_p_$(plim).pdf")
    letters_h = reshape(letters[1:18],3,6)
    #all num of sig increases in 20y scales
    ff = Figure(size=(1500,2000))
    siggi = @. replace(sig_incr_matr_sca_one, missing => 0)
    allma_sca = @. maximum(sum(siggi, dims=3),dims=2)
    max_sca = maximum(@. maximum(allma_sca))
    n_max = Int(max_sca)
    co = cgrad(:amp,  n_max+1, rev = false, categorical = true)

    indin5s = []
    indin10s = []
    for el10 in scaleranges[3]
        push!(indin5s, findfirst(==(el10), scaleranges[1]))
        push!(indin10s, findfirst(==(el10), scaleranges[2]))
    end

    for j = 1:3  
        if j == 1
            sig_incr_matr = sig_incr_matr_sca_one
            ind_n = L"\textbf{\mathrm{\hat{𝗐}^𝟤}}"
        elseif j == 2
            sig_incr_matr = sig_incr_matr_hurst_one
            ind_n = L"\textbf{\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}}}"
        elseif j == 3
            sig_incr_matr = sig_incr_matr_wave_both_one
            ind_n = L"\textbf{\mathrm{\hat{𝗐}^𝟤} \text{𝖺𝗇𝖽 }\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}}}"
        end
        whichletter = 1
        for i = 1:3
            if i == 1
                s1s = sort(unique([s[1] for s in scaleranges[i][indin5s]])) #reshape(sranges, -1,2)
                s2s = sort(unique([s[2] for s in scaleranges[i][indin5s]]))
            elseif i ==2
                s1s = sort(unique([s[1] for s in scaleranges[i][indin10s]])) #reshape(sranges, -1,2)
                s2s = sort(unique([s[2] for s in scaleranges[i][indin10s]]))
            else
                s1s = sort(unique([s[1] for s in scaleranges[i]])) #reshape(sranges, -1,2)
                s2s = sort(unique([s[2] for s in scaleranges[i]]))
            end
            replace!(sig_incr_matr[i], missing => 0)
            n_incr = sum(sig_incr_matr[i], dims = 3)

            for k =  size(n_incr)[1]:-1:1
                ga = ff[whichletter,j] = GridLayout()
                ax =  Axis(ga[1,1], 
                            title = ind_n,
                            titlesize=20,
                            titlefont = :bold,
                            xticks = (1:length(s1s), string.(s1s)),
                            xlabel = "s₁ (years)",
                            yticksvisible = true,
                            xminorgridwidth = 1.0, 
                            xminorgridcolor = :grey30 ,
                            xminorgridvisible = true, xgridvisible = false,
                            yticks = (1:length(s2s), string.(s2s)), 
                            ylabel = "s₂ (years)",
                            yminorgridwidth = 1.0, 
                            yminorgridcolor = :grey30, 
                            yminorgridvisible = true, ygridvisible = false)
                n_incr_matrix = Array{Union{Missing, Float64}}(missing, length(s1s), length(s2s))
                if i == 1
                    for (isr, sr) in enumerate(scaleranges[i])
                        if sr in scaleranges[i][indin5s]
                            n_incr_matrix[findfirst(==(sr[1]), s1s), findfirst(==(sr[2]), s2s)] =  n_incr[k,isr,1]
                        end
                    end
                elseif i == 2
                    for (isr, sr) in enumerate(scaleranges[i])
                        if sr in scaleranges[i][indin10s]
                            n_incr_matrix[findfirst(==(sr[1]), s1s), findfirst(==(sr[2]), s2s)] =  n_incr[k,isr,1]
                        end
                    end
                elseif i == 3
                    for (isr, sr) in enumerate(scaleranges[i])
                        n_incr_matrix[findfirst(==(sr[1]), s1s), findfirst(==(sr[2]), s2s)] =  n_incr[k,isr,1]
                    end
                end
               
                rowsize!(ga,1, Auto(1))

                hm = heatmap!(ax, n_incr_matrix, colormap = co, colorrange = (0,n_max))
                translate!(hm, 0,0,-100)
                Colorbar(ga[1,2], colormap = co, colorrange = (0,n_max+1),
                ticks = (0.5:1:Int(n_max)+0.5, string.(0:1:Int(n_max))),
                label = "No. of significant EWS")

                if j == 1
                    Label(ga[1,1,Left()], ylabs[i][k], padding = (0.0,70.0,0.0,0.0), valign = :center, font = :bold, fontsize = 20, rotation = pi/2)
                end

                Label(ga[1,1,TopLeft()], letters_h[j,whichletter], fontsize = 20,
                        font = :bold,  padding = (0,5,5,0),
                        halign = :left)
                whichletter +=1
                
            end
        end      

    end
    Label(ff[1,1:end,Top()], "Number of significant EWS in (s₁-s₂) year bands of records with 5-year resolution", padding = (0.0,0.0,40.0,40.0), valign = :bottom, font = :bold, fontsize = 20)
    Label(ff[2,1:end,Top()], "Number of significant EWS in (s₁-s₂) year bands of records with 10-year resolution", padding = (0.0,0.0,40.0,40.0), valign = :bottom, font = :bold, fontsize = 20)
    Label(ff[4,1:end,Top()], "Number of significant EWS in (s₁-s₂) year bands of records with 20-year resolution", padding = (0.0,0.0,40.0,40.0), valign = :bottom, font = :bold, fontsize = 20)

    colgap!(ff.layout, 40)
    if saving
        save(saveto,ff)
    end
    if showing
        display(ff)
    end
end


#Fig S8 and S9
function plot_increases_wavelet_own_scales_split(incr_matrix_sca, sig_incr_matr_sca_one,incr_matrix_hurst,sig_incr_matr_hurst_one,incr_matrix_wave_both,sig_incr_matr_wave_both_one,
    scaleranges, ylabs;showing = true, saving = false, saveto = "paper/incr_wavelet_own_scales_lowpass_$(lowpass)_smoothw_$(smoothw)_p_$(plim).pdf")
    letters_h = reshape(letters[1:18],3,6)
    # all wavelet csd for all cores - in own ranges
    co = cgrad([:steelblue, :darkred, :darkred],[0,0.49, 0.9], alpha = [0.15,0.15,1.0])
    for i = 1:3
        if i == 1
            s = 55*24+65
        elseif i == 2
            s=55*2*24+2*65 -100
        else
            s = 10*3*24+3*65 +120
        end
        ff = Figure(size=(1500,s))
        for j = 1:3  
            if j == 1
                incr_ma = incr_matrix_sca
                sig_incr_matr = sig_incr_matr_sca_one
                ind_n = L"\textbf{\mathrm{\hat{𝗐}^𝟤}}" 
            elseif j == 2
                incr_ma = incr_matrix_hurst
                sig_incr_matr = sig_incr_matr_hurst_one
                ind_n = L"\textbf{\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}}}" 
            elseif j == 3
                incr_ma = incr_matrix_wave_both
                sig_incr_matr = sig_incr_matr_wave_both_one
                ind_n = L"\textbf{\mathrm{\hat{𝗐}^𝟤} \text{𝖺𝗇𝖽 }\mathrm{\hat{𝖧}}^\text{𝗅𝗈𝖼}}"
            end

            whichletter = 1

            for k = size(incr_ma[i])[1]:-1:1
                ga = ff[whichletter,j] = GridLayout()
                ax =  Axis(ga[1,1], 
                            title = ind_n,
                            titlesize=20,
                            titlefont = :bold,
                            xticks = (1:17, ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]) ,
                            xticklabelrotation=pi/3,
                            xreversed = true,
                            xlabel = "Event",
                            yticksvisible = false,
                            xminorgridwidth = 1.0, 
                            xminorgridcolor = :grey30 ,
                            xminorgridvisible = true, xgridvisible = false,
                            yticks = (1:length(scaleranges[i]), string.(scaleranges[i])),
                            yminorgridwidth = 1, 
                            yminorgridcolor = :grey30, #:grey5,
                            yminorgridvisible = true, ygridvisible = false)
                hm = heatmap!(ax,(incr_ma[i][k,:,:] + sig_incr_matr[i][k,:,:])', colormap = co, colorrange = (0,2))
                translate!(hm, 0,0,-100)

                ax.ylabel = "(s₁,s₂) (years)"

                Label(ga[1,1,TopLeft()], letters_h[j,whichletter], fontsize = 20,
                font = :bold,  padding = (0,5,5,0),
                halign = :left)
                whichletter +=1

                if j == 1
                    Label(ga[1,1,Left()], ylabs[i][k], padding = (0.0,110.0,0.0,0.0), valign = :center, font = :bold, fontsize = 20, rotation = pi/2)
                end
                
                rowsize!(ga,1,Auto(1))
            end
        end
        elems = [
            [MarkerElement(color = :white, marker = :rect, markersize = 25,
            strokecolor = :grey, strokewidth = 0.7)],
            [MarkerElement(color = co[1], marker = :rect, markersize = 25,
            strokecolor = :grey, strokewidth = 0.7)],
            [MarkerElement(color = (:darkred, 0.15), marker = :rect, markersize = 25,
            strokecolor = :grey, strokewidth = 0.7)],
            [MarkerElement(color = (:darkred, 1.0), marker = :rect, markersize = 25,
            strokecolor = :grey, strokewidth = 0.7)],
            ]
        Legend(ff[end+1,:], 
            elems,
            ["undefined", "decreasing", "increasing","significantly increasing (p<$(plim))"],
            rowgap = 20, 
            framevisible = true, 
            framecolor = :grey70,
            tellwidth = false,
            orientation = :horizontal)

        lbls = ["Linear trends of indicators in (s₁-s₂) year bands of records with 5-year resolution",
                "Linear trends of indicators in (s₁-s₂) year bands of records with 10-year resolution",
                "Linear trends of indicators in (s₁-s₂) year bands of records with 20-year resolution"]
        
        Label(ff[1,1:end,Top()], lbls[i], padding = (0.0,0.0,40.0,40.0), valign = :bottom, font = :bold, fontsize = 20)

        colgap!(ff.layout, 40)
        vs = [5,10,20]
        saveto_local = join([split(saveto,".pdf")[1]*"_$(vs[i])_yrs",".pdf"])
  
        if saving
            save(saveto_local,ff)
        end
        if showing
            display(ff)
        end
    end
end


#Fig A4 and S22    
function plot_increases_wavelet_20_scales2(incr_matrix_sca, sig_incr_matr_sca_one,incr_matrix_hurst,sig_incr_matr_hurst_one,incr_matrix_wave_both,sig_incr_matr_wave_both_one,
    scaleranges, ylabs;showing = true, saving = false, saveto = "paper/incr_wavelet_20_scales_lowpass_$(lowpass)_smoothw_$(smoothw).pdf")
    # all wavelet csd for all cores - in 20y res
    letters_h = reshape(letters[1:18],3,6)
    ff = Figure(size=(1600,1800))
    co = cgrad([:steelblue, :darkred, :darkred],[0,0.49, 0.9], alpha = [0.15,0.15,1.0])
    
    
    indin5s = []
    indin10s = []
    for el20 in scaleranges[3]
        push!(indin5s, findfirst(==(el20), scaleranges[1]))
        push!(indin10s, findfirst(==(el20), scaleranges[2]))
    end
    for j = 1:3  
        if j == 1
            incr_ma = incr_matrix_sca
            sig_incr_matr = sig_incr_matr_sca_one
            ind_n = L"\textbf{\mathrm{\hat{𝗐}^𝟤}}"
        elseif j == 2
            incr_ma = incr_matrix_hurst
            sig_incr_matr = sig_incr_matr_hurst_one
            ind_n = L"\textbf{\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}}}"
        elseif j == 3
            incr_ma = incr_matrix_wave_both
            sig_incr_matr = sig_incr_matr_wave_both_one
            ind_n = L"\textbf{\mathrm{\hat{𝗐}^𝟤} \text{𝖺𝗇𝖽 }\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}}}"
        end

        whichletter = 1

        for i = 1:3
            for k = size(incr_ma[i])[1]:-1:1
                ga = ff[whichletter,j] = GridLayout()
                ax =  Axis(ga[1,1], 
                            title = ind_n,
                            titlesize=20,
                            titlefont = :bold,
                            xticks = (1:17, ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]) ,
                            xticklabelrotation=pi/3,
                            xreversed = true,
                            xlabel = "Event",
                            yticksvisible = false,
                            xminorgridwidth = 1.0, 
                            xminorgridcolor = :grey30 ,
                            xminorgridvisible = true, xgridvisible = false,
                            yticks = (1:length(scaleranges[3]), string.(scaleranges[3])),
                            yminorgridwidth = 1, 
                            yminorgridcolor = :grey30, 
                            yminorgridvisible = true, ygridvisible = false)
                if i == 1
                    hm = heatmap!(ax,(incr_ma[i][k,indin5s,:] + sig_incr_matr[i][k,indin5s,:])', colormap = co, colorrange = (0,2))    
                elseif i == 2
                    hm = heatmap!(ax,(incr_ma[i][k,indin10s,:] + sig_incr_matr[i][k,indin10s,:])', colormap = co, colorrange = (0,2))
                else
                    hm = heatmap!(ax,(incr_ma[i][k,:,:] + sig_incr_matr[i][k,:,:])', colormap = co, colorrange = (0,2))
                end
                translate!(hm, 0,0,-100)

                ax.ylabel = "(s₁,s₂) (years)"

                Label(ga[1,1,TopLeft()], letters_h[j,whichletter], fontsize = 20,
                    font = :bold, padding = (0,5,5,0),
                    halign = :left)
                whichletter +=1

                if j == 1
                    Label(ga[1,1,Left()], ylabs[i][k], padding = (0.0,110.0,0.0,0.0), valign = :center, font = :bold, fontsize = 20, rotation = pi/2)
                end
                
                rowsize!(ga,1,Auto(1))
            end
        end
    end
    elems = [
        [MarkerElement(color = :white, marker = :rect, markersize = 25,
        strokecolor = :grey, strokewidth = 0.7)],
        [MarkerElement(color = co[1], marker = :rect, markersize = 25,
        strokecolor = :grey, strokewidth = 0.7)],
        [MarkerElement(color = (:darkred, 0.15), marker = :rect, markersize = 25,
        strokecolor = :grey, strokewidth = 0.7)],
        [MarkerElement(color = (:darkred, 1.0), marker = :rect, markersize = 25,
        strokecolor = :grey, strokewidth = 0.7)],
        ]
    Legend(ff[end+1,:],
        elems,
        ["undefined", "decreasing", "increasing","significantly increasing (p<$(plim))"],
        rowgap = 20, 
        framevisible = true,
        framecolor = :grey70,
        tellwidth = false,
        orientation = :horizontal)
    
    Label(ff[1,1:end,Top()], "Linear trends of indicators in (s₁-s₂) year bands of records with 5-year resolution", padding = (0.0,0.0,30.0,0.0), valign = :bottom, font = :bold, fontsize = 20)
    Label(ff[2,1:end,Top()], "Linear trends of indicators in (s₁-s₂) year bands of records with 10-year resolution", padding = (0.0,0.0,30.0,3.0), valign = :bottom, font = :bold, fontsize = 20)
    Label(ff[4,1:end,Top()], "Linear trends of indicators in (s₁-s₂) year bands of records with 20-year resolution", padding = (0.0,0.0,30.0,3.0), valign = :bottom, font = :bold, fontsize = 20)
    
    colgap!(ff.layout, 80)
    if saving
        save(saveto,ff)
    end
    if showing
        display(ff)
    end
end


#Fig S11
function plot_common_increases_wavelet(sig_incr_matr_sca_one,sig_incr_matr_hurst_one,sig_incr_matr_wave_both_one,
    scaleranges;showing = true, saving = false, saveto = "paper/common_incr_wavelet_lowpass_$(lowpass)_smoothw_$(smoothw)_p_$(plim).pdf")
    # common signidicant wavelet increases
    letters_h = reshape(letters[1:18],3,6)
    ff = Figure(size=(1500,1800))
    co5 = cgrad(:gist_earth, 5, categorical = true, rev=true, alpha = [1.0,0.5,0.5,0.5,1.0])
    co2 = cgrad(:gist_earth, 2, categorical = true, rev=true, alpha = [1.0,1.0])
    for j = 1:3
        if j == 1
            sig_incr_matr = sig_incr_matr_sca_one
            ind_n =L"\textbf{\mathrm{\hat{𝗐}^𝟤}}" 
        elseif j == 2
            sig_incr_matr = sig_incr_matr_hurst_one
            ind_n = L"\textbf{\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}}}" 
        elseif j == 3
            sig_incr_matr = sig_incr_matr_wave_both_one
            ind_n = L"\textbf{\mathrm{\hat{𝗐}^𝟤} \text{𝖺𝗇𝖽 }\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}}}"
        end

        for l = 1:3
            replace!(sig_incr_matr[l], missing =>0)
        end

        whichletter = 1

        for i = 1:3
            ga = ff[whichletter,j] = GridLayout()

            ax =  Axis(ga[1,1], 
                            title = ind_n,
                            titlesize=20,
                            titlefont = :bold,
                            xticks = (1:17, ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]) ,
                            xticklabelrotation=pi/3,
                            xreversed = true,
                            xlabel = "Event",
                            yticksvisible = false,
                            xminorgridwidth = 1.0, 
                            xminorgridcolor = :grey30 ,
                            xminorgridvisible = true, xgridvisible = false,
                            yticks = (1:length(scaleranges[i]), string.(scaleranges[i])),
                            yminorgridwidth = 1, 
                            yminorgridcolor = :grey30, 
                            yminorgridvisible = true, ygridvisible = false)
            
            if i == 1
                # compare NGRIP
                commons = zeros(size(sig_incr_matr[3][end,:,:]))
                indin10s = []
                indin5s = []
                for el20 in scaleranges[3]
                    push!(indin10s, findfirst(==(el20), scaleranges[2]))
                    push!(indin5s, findfirst(==(el20), scaleranges[1]))
                end
                inc12 = sig_incr_matr[1][end,indin5s,:] .== sig_incr_matr[2][end,indin10s,:] .==1
                inc13 = sig_incr_matr[1][end,indin5s,:] .== sig_incr_matr[3][end,:,:] .==1
                inc23 = sig_incr_matr[2][end,indin10s,:] .== sig_incr_matr[3][end,:,:] .==1
                incall = sig_incr_matr[1][end,indin5s,:] .== sig_incr_matr[2][end,indin10s,:] .== sig_incr_matr[3][end,:,:] .== 1
                commons[inc12] .= 1
                commons[inc13] .= 2
                commons[inc23] .= 3
                commons[incall] .= 4
                hm = heatmap!(ax, commons', colormap = co5, colorrange=(0,4))
                translate!(hm, 0,0,-100)
                ax.yticks = (1:length(scaleranges[3]), string.(scaleranges[3]))

            elseif i == 2
                commons = zeros(size(sig_incr_matr[2][1,:,:]))
                inc12 = sig_incr_matr[2][1,:,:] .== sig_incr_matr[2][2,:,:] .== 1
                commons[inc12] .= 1
                hm = heatmap!(ax, commons', colormap = co2, colorrange=(0,1))
                translate!(hm, 0,0,-100)
                
                ax.yticks = (1:length(scaleranges[2]), string.(scaleranges[2]))

            elseif i == 3
                #20y
                commons = zeros(size(sig_incr_matr[3][1,:,:]))
                inc12 = sig_incr_matr[3][1,:,:] .== sig_incr_matr[3][2,:,:] .== 1
                inc13 = sig_incr_matr[3][1,:,:] .== sig_incr_matr[3][3,:,:] .== 1
                inc23 = sig_incr_matr[3][2,:,:] .== sig_incr_matr[3][3,:,:] .== 1
                incall = sig_incr_matr[3][1,:,:] .== sig_incr_matr[3][2,:,:] .== sig_incr_matr[3][3,:,:] .== 1
                commons[inc12] .= 1
                commons[inc13] .= 2
                commons[inc23] .= 3
                commons[incall] .= 4
                
                ax.yticks = (1:length(scaleranges[3]), string.(scaleranges[3]))
               

                hm = heatmap!(ax, commons', colormap = co5, colorrange=(0,4))
                translate!(hm, 0,0,-100)
            end

            ax.ylabel = "(s₁,s₂) (years)"
            Label(ga[1,1,TopLeft()], letters_h[j,whichletter], fontsize = 20,
                font = :bold, padding = (0,5,5,0),
                halign = :left)
            whichletter +=1
        end
    end
    elems5 = [[MarkerElement(color = co5[ci], marker = :rect, markersize = 25,
    strokecolor = :grey, strokewidth = 0.7)] for ci = 1:5]
    elems2 = [[MarkerElement(color = co2[ci], marker = :rect, markersize = 25,
        strokecolor = :grey, strokewidth = 0.7)] for ci = 1:2]
    labels1 = ["none", "5 & 10 y","5 & 20 y", "10 & 20 y","5, 10 & 20 y"]
    labels2 = ["none", "NGRIP & NEEM"]
    labels3 = ["none", "GRIP & GISP2 ","NGRIP & GISP2", "NGRIP & GRIP","NGRIP, GRIP & GISP2"]


    Legend(ff[1,4], 
        elems5,
        labels1,
        colgap = 10, 
        framevisible = false,
        framecolor =:grey70,
        tellwidth = true,
        tellheight = false,
        halign = :left
        )
    Legend(ff[2,4],
        elems2,
        labels2,
        framevisible = false, 
        framecolor =:grey70,
        tellwidth = true,
        halign = :left,
        gridshalign = true,
        tellheight = false,
        )
    Legend(ff[3,4], 
        elems5,
        labels3,
        colgap = 10, 
        framevisible = false, 
        framecolor =:grey70,
        tellwidth = true,
        tellheight =false,
        halign = :left,
        gridshalign = true,
        )

    Label(ff[1,1:end,Top()], "Common significant EWS in in (s₁-s₂) year bands of the NGRIP record with 5-, 10-, and 20-year resolutions", padding = (0.0,0.0,40.0,40.0), valign = :bottom, font = :bold, fontsize = 20)
    Label(ff[2,1:end,Top()], "Common significant EWS in in (s₁-s₂) year bands of records with 10-year resolution", padding = (0.0,0.0,40.0,40.0), valign = :bottom, font = :bold, fontsize = 20)
    Label(ff[3,1:end,Top()], "Common significant EWS in in (s₁-s₂) year bands of records with 20-year resolution", padding = (0.0,0.0,40.0,40.0), valign = :bottom, font = :bold, fontsize = 20)
    colgap!(ff.layout, 40)
    rowsize!(ff.layout, 2, Auto(6))
    if saving
        save(saveto,ff)
    end
    if showing
        display(ff)
    end

end


#compare sign. CSD increases in NGRIP, 10y, and 20y
#Fig S7
function plot_common_increases_csd(sig_incr_matr_var_one, sig_incr_matr_ac_one,sig_incr_matr_both_one ;showing = true, saving = false, saveto = "paper/common_increases_csd_lowpass_$(lowpass)_p_$(plim).pdf")
    ff = Figure(size=(1200,750))
    letters_h = reshape(letters[1:18],3,6)
    co5 = cgrad(:gist_earth, 5, categorical = true, rev=true, alpha = [1.0,0.5,0.5,0.5,1.0])
    co2 = cgrad(:gist_earth, 2, categorical = true, rev=true, alpha = [1.0,1.0])
    for j = 1:3
        if j ==1
            sig_incr_matr = sig_incr_matr_var_one
            ind_n = L"\textbf{\mathrm{𝖵}}"
        elseif j == 2
            sig_incr_matr = sig_incr_matr_ac_one
            ind_n = L"\textbf{\mathrm{\alpha_𝟣}}"
        elseif j == 3
            sig_incr_matr = sig_incr_matr_both_one
            ind_n = L"\textbf{\mathrm{𝖵} \text{𝖺𝗇𝖽 } \mathrm{\alpha_𝟣}}"
        end
        for l in 1:3
            replace!(sig_incr_matr[l], missing =>0)
        end

        whichletter = 1

        for i = 1:3
            commons = zeros(17)
            ga = ff[whichletter,j] = GridLayout()

            ax =  Axis(ga[1,1], 
                        title = ind_n,
                        titlesize=20,
                        titlefont = :bold,
                        xticks = (1:17, ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]) ,
                        xticklabelrotation=pi/3,
                        xreversed = true,
                        yticksvisible = false,
                        xminorgridwidth = 1.0, 
                        xminorgridcolor = :grey30 ,
                        xminorgridvisible = true, xgridvisible = false,
                        yminorgridwidth = 3.2, 
                        yminorgridcolor = :grey30, 
                        yminorgridvisible = false, ygridvisible = false)
            hideydecorations!(ax, label = true)

            if i == 1
                #NGRIP
                inc12 = sig_incr_matr[1][end,:] .== sig_incr_matr[2][end,:] .==1
                inc13 = sig_incr_matr[1][end,:] .== sig_incr_matr[3][end,:] .==1
                inc23 = sig_incr_matr[2][end,:] .== sig_incr_matr[3][end,:] .==1
                incall = sig_incr_matr[1][end,:] .== sig_incr_matr[2][end,:] .== sig_incr_matr[3][end,:] .== 1
                commons[inc12] .= 1
                commons[inc13] .= 2
                commons[inc23] .= 3
                commons[incall] .= 4

                commons = reshape(commons, 17,1)
                hm = heatmap!(ax, commons, colormap = co5, colorrange=(0,4))
                translate!(hm, 0,0,-100)
                                    
            elseif i == 2
                #10y
                inc12 = sig_incr_matr[2][1,:] .== sig_incr_matr[2][2,:] .== 1
                commons[inc12] .= 1
                commons = reshape(commons, 17,1)
                hm = heatmap!(ax, commons, colormap = co2, colorrange=(0,1))
                translate!(hm, 0,0,-100)
            elseif i == 3
                #20y
                inc12 = sig_incr_matr[3][1,:] .== sig_incr_matr[3][2,:] .== 1
                inc13 = sig_incr_matr[3][1,:] .== sig_incr_matr[3][3,:] .== 1
                inc23 = sig_incr_matr[3][2,:] .== sig_incr_matr[3][3,:] .== 1
                incall = sig_incr_matr[3][1,:] .== sig_incr_matr[3][2,:] .== sig_incr_matr[3][3,:] .== 1
                commons[inc12] .= 1
                commons[inc13] .= 2
                commons[inc23] .= 3
                commons[incall] .= 4
                commons = reshape(commons, 17,1)
                hm = heatmap!(ax, commons, colormap = co5, colorrange=(0,4))
                translate!(hm, 0,0,-100)
            end
            ax.xlabel = "Event"
            rowsize!(ff.layout, i, Relative(1/3))
            Label(ga[1,1,TopLeft()], letters_h[j,whichletter], fontsize = 20,
                font = :bold, padding = (0,5,5,0),
                halign = :left, valign =:bottom)
            whichletter +=1
        end
    end
    elems5 = [[MarkerElement(color = co5[ci], marker = :rect, markersize = 25,
                strokecolor = :grey, strokewidth = 0.7)] for ci = 1:5]
    elems2 = [[MarkerElement(color = co2[ci], marker = :rect, markersize = 25,
                strokecolor = :grey, strokewidth = 0.7)] for ci = 1:2]
    labels1 = ["none", "5 & 10 y","5 & 20 y", "10 & 20 y","5, 10 & 20 y"]
    labels2 = ["none", "NGRIP & NEEM"]
    labels3 = ["none", "GRIP & GISP2 ","NGRIP & GISP2", "NGRIP & GRIP","NGRIP, GRIP & GISP2"]
    
    
    Legend(ff[1,4],
        elems5,
        labels1,
        colgap = 10, 
        framevisible = false, 
        framecolor =:grey70,
        tellwidth = true,
        tellheight = false,
        halign = :left
        )
    Legend(ff[2,4],
        elems2,
        labels2,
        framevisible = false,
        framecolor =:grey70,
        tellwidth = true,
        halign = :left,
        gridshalign = true,
        tellheight = false,
        )
    Legend(ff[3,4], 
        elems5,
        labels3,
        colgap = 10, 
        framevisible = false, 
        framecolor =:grey70,
        tellwidth = true,
        tellheight =false,
        halign = :left,
        gridshalign = true
        )
        Label(ff[1,1:end,Top()], "Common significant EWS in 100-year high-pass filtered NGRIP record with 5-, 10-, and 20-year resolutions", padding = (0.0,0.0,40.0,40.0), valign = :bottom, font = :bold, fontsize = 20)
        Label(ff[2,1:end,Top()], "Common significant EWS in 100-year high-pass filtered records with 10-year resolution", padding = (0.0,0.0,40.0,40.0), valign = :bottom, font = :bold, fontsize = 20)
        Label(ff[3,1:end,Top()], "Common significant EWS in 100-year high-pass filtered records with 20-year resolution", padding = (0.0,0.0,40.0,40.0), valign = :bottom, font = :bold, fontsize = 20)
    colsize!(ff.layout, 4, Relative(1/5))
    rowgap!(ff.layout, 20)
    if saving
        save(saveto,ff)
    end
    if showing
        display(ff)
    end
end


saving = true
for lowpass in [true, false]
    smoothw=false
    #@show lowpass, smoothw
    local incr_matrix_var,
    incr_matrix_ac,
    incr_matrix_both,
    sig_incr_matr_var_one,
    sig_incr_matr_var_two,
    sig_incr_matr_ac_one,
    sig_incr_matr_ac_two,
    sig_incr_matr_both_one,
    sig_incr_matr_both_two,
    ylabs,
    incr_matrix_sca,
    sig_incr_matr_sca_one,
    sig_incr_matr_sca_two,
    incr_matrix_hurst,
    sig_incr_matr_hurst_one,
    sig_incr_matr_hurst_two,
    incr_matrix_wave_both,
    sig_incr_matr_wave_both_one,
    sig_incr_matr_wave_both_two,
    scaleranges = get_incr_matrices(lowpass = lowpass, smoothw =smoothw, plim = plim);

    if lowpass
        plot_num_incr_wavelet_own_scales2(sig_incr_matr_sca_one,sig_incr_matr_hurst_one,sig_incr_matr_wave_both_one,
            scaleranges, ylabs;showing = showing, saving = saving, saveto = "do_ews_across_greenland_ice_cores/figures/figA5.pdf")
        plot_num_incr_wavelet_20_scales(sig_incr_matr_sca_one,sig_incr_matr_hurst_one,sig_incr_matr_wave_both_one,
            scaleranges, ylabs;showing = showing, saving = saving, saveto = "do_ews_across_greenland_ice_cores/figures/figS10.pdf")
        plot_increases_wavelet_own_scales_split(incr_matrix_sca, sig_incr_matr_sca_one,incr_matrix_hurst,sig_incr_matr_hurst_one,incr_matrix_wave_both,sig_incr_matr_wave_both_one,
            scaleranges, ylabs;showing = showing, saving = saving, saveto = "do_ews_across_greenland_ice_cores/figures/figS8_S9.pdf")
        plot_increases_wavelet_20_scales2(incr_matrix_sca, sig_incr_matr_sca_one,incr_matrix_hurst,sig_incr_matr_hurst_one,incr_matrix_wave_both,sig_incr_matr_wave_both_one,
            scaleranges, ylabs;showing = showing, saving = saving, saveto = "do_ews_across_greenland_ice_cores/figures/figA4.pdf")
        plot_common_increases_wavelet(sig_incr_matr_sca_one,sig_incr_matr_hurst_one,sig_incr_matr_wave_both_one,
            scaleranges;showing = showing, saving = saving, saveto = "do_ews_across_greenland_ice_cores/figures/figS11.pdf")
        plot_common_increases_csd(sig_incr_matr_var_one, sig_incr_matr_ac_one,sig_incr_matr_both_one; 
            showing = showing, saving = saving, saveto = "do_ews_across_greenland_ice_cores/figures/figS7.pdf")
    else
        plot_num_incr_wavelet_own_scales2(sig_incr_matr_sca_one,sig_incr_matr_hurst_one,sig_incr_matr_wave_both_one,
            scaleranges, ylabs;showing = showing, saving = saving, saveto = "do_ews_across_greenland_ice_cores/figures/figS23.pdf")
        plot_increases_wavelet_20_scales2(incr_matrix_sca, sig_incr_matr_sca_one,incr_matrix_hurst,sig_incr_matr_hurst_one,incr_matrix_wave_both,sig_incr_matr_wave_both_one,
            scaleranges, ylabs;showing = showing, saving = saving, saveto = "do_ews_across_greenland_ice_cores/figures/figS22.pdf")
    end
end


## all CSD for all the cores with indicator TS plotted
#Fig 6 and S17
function plot_all_csd_ts2_label(ylabs;lowpass = lowpass, plim = plim, showing = true, saving = false, legend = false, saveto = "paper/all_csd_ts2_label_lowpass_$(lowpass)_p_$(plim).pdf")
    f_csd_ts = Figure(size=(1200,1250))
    letters_h = reshape(letters[1:18],2,9)
    
    whichletter = 0

    function good_csd_files(v,type,lowpass)
        if type == "NGRIP5"
            if lowpass
                return occursin("w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",v) && !occursin("NIKLAS_", v)
            else
                return  occursin("w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",v) && !occursin("NIKLAS_", v)
            end
        elseif type == "10y"
            if lowpass
                return occursin("w_200_normed_filt", v) && occursin("lp_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2", v)
            else
                return occursin("w_200_normed_filt", v) && occursin("gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2", v) && !occursin("_lp_", v)
            end
        elseif type == "20y"
            return occursin("w_200_normed_filt", v) && occursin("gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2", v)
        end
    end

    for (it,type) in enumerate(["NGRIP5", "10y", "20y"])
        ice_cores, var_names, vars, ac_names, acs, sca_names, scas, hurst_names, hursts, sranges = load_data(type)

        good = @. good_csd_files(var_names,type,lowpass)
        
        for (i,(v,a)) in enumerate(zip(reverse(vars[good]), reverse(acs[good])))
            whichletter +=1
            n_v = 0 
            n_a = 0
            n_both = 0
            which_v = []
            for k = 1:2
                gl = f_csd_ts[whichletter,k] = GridLayout()
                
                

                if k == 1
                    axv = Axis(gl[1,1], 
                            xlabel = "Age (kyr b2k)", 
                            ylabel = L"\mathrm{𝖵}", 
                            ylabelsize = 18,
                            ylabelfont = :bold,
                            xticks = (10_000:5_000:60_000, string.(10:5:60)), 
                            )
                    axv.xreversed = true
                    
                    ax_event = Axis(gl[1,1],xaxisposition = :top, xticklabelrotation = pi*0.4)
                    ax_event.xticks = (GI_onsets, event_labels)
                    ax_event.xreversed = true
                    hidespines!(ax_event)
                    hideydecorations!(ax_event)

                    if whichletter ∈ [1,3,6]
                        axv.xlabel = "Age (kyr b2k)"
                    else
                        hidexdecorations!(axv, grid=false)
                        hidespines!(axv, :b) 
                    end

                    if i>1
                        hidespines!(axv,:t)
                        hidexdecorations!(ax_event, grid=false)
                    end

                    for ev in 1:17
                        va = 0.4
                        vl = 2.0
                        vsc = (:white,1.0)
                        vc = :black
                        if typeof(v.slopes[ev]) != Nothing
                            if v.slopes[ev] >0 
                                vc = :red
                                vsc = (:darkred,0.1)
                                if v.p_one[ev] < plim
                                    va = 0.9
                                    vl = 3.5
                                    vsc = (:darkred,0.8)
                                    n_v +=1
                                    push!(which_v, ev)
                                end
                            
                            else
                                vc = :blue
                                vsc = (:steelblue,0.1)
                            end
                        end
                        vspan!(axv, [GI_onsets[ev]], [GS_onsets[ev]], color = vsc)
                        pv = Polynomials.fit(v.times[ev][findall(!ismissing, v.vals[ev])], convert.(Float64,v.vals[ev][findall(!ismissing, v.vals[ev])]),1)
                        lines!(axv,v.times[ev], v.vals[ev], color = :black, alpha = va)#, linewidth = vl)
                        lines!(axv, v.times[ev], pv.(v.times[ev]), color = vc, alpha = va, linewidth = vl)
                        vlines!(axv, [GI_onsets[ev]], color = :darkred, alpha = va, linewidth = vl)
                    end
                    xlims!(axv,nothing,nothing)
                    axv.xreversed = true
                    xlims!(ax_event, axv.xaxis.attributes.limits[]...)
                    ax_event.xreversed = true
                elseif k == 2
                    axa = Axis(gl[1,1], 
                        ylabel =  L"\mathrm{\alpha_𝟣}",
                        ylabelsize = 18,
                        ylabelfont = :bold,
                        xticks = (10_000:5_000:60_000, string.(10:5:60)), 
                        )
                    
                    axa.xreversed = true
                   

                    ax_event = Axis(gl[1,1],xaxisposition = :top, xticklabelrotation = pi*0.4)
                    ax_event.xticks = (GI_onsets, event_labels)
                    
                    ax_event.xreversed = true
                    hidespines!(ax_event)
                    hideydecorations!(ax_event)
                    
                    if whichletter ∈ [1,3,6]
                        axa.xlabel = "Age (kyr b2k)"
                    else
                        hidexdecorations!(axa, grid=false)
                        hidespines!(axa, :b) 
                    end

                    if i>1
                        hidespines!(axa,:t)
                        hidexdecorations!(ax_event, grid=false)
                    end
                    
                    
                    for ev in 1:17
                        aa = 0.4
                        al = 2.0
                        asc = (:white,1.0)
                        ac = :black
                        if typeof(a.slopes[ev]) != Nothing
                            if a.slopes[ev] >0 
                                ac = :red
                                asc = (:darkred,0.1)
                                if a.p_one[ev] < plim
                                    aa = 0.9
                                    al = 3.5
                                    asc = (:darkred,0.8)
                                    n_a +=1
                                    if ev in which_v
                                        n_both +=1
                                    end
                                end
                            
                            else
                                ac = :blue
                                asc = (:steelblue,0.1)
                            end
                        end
                        vspan!(axa, [GI_onsets[ev]], [GS_onsets[ev]], color = asc)
                        pa = Polynomials.fit(a.times[ev][findall(!ismissing, a.vals[ev])], convert.(Float64,a.vals[ev][findall(!ismissing, a.vals[ev])]),1)
                        lines!(axa,a.times[ev], a.vals[ev], color = :black, alpha = aa)#, linewidth = vl)
                        lines!(axa, a.times[ev], pa.(a.times[ev]), color = ac, alpha = aa, linewidth = al)
                        vlines!(axa, [GI_onsets[ev]], color = :darkred, alpha = aa, linewidth = al)

                    end

                    xlims!(axa,nothing,nothing)
                    axa.xreversed = true
                    xlims!(ax_event, axa.xaxis.attributes.limits[]...)
                    ax_event.xreversed = true
                end


                Label(gl[1,1,TopLeft()], letters_h[k,whichletter], fontsize = 20,
                            font = :bold, padding = (0,20,-10,15),
                            halign = :left, valign =:bottom)
            end
            if i == 1
                Label(f_csd_ts[whichletter,1:end,Top()], L"\textbf{\textrm{%$(reverse(ylabs2b[it])[i]): }\mathrm{𝗇_{𝖵} = %$(sstring(n_v)),} \mathrm{𝗇_{\alpha_𝟣} = %$(sstring(n_a)),} \mathrm{𝗇_{\text{𝖻𝗈𝗍𝗁}} = %$(sstring(n_both))}}", padding = (0.0,0.0,20.0,20.0), valign = :center, font=:bold,fontsize = 20)
            else
                Label(f_csd_ts[whichletter,1:end,Top()], L"\textbf{\textrm{%$(reverse(ylabs2b[it])[i]): }\mathrm{𝗇_{𝖵} = %$(sstring(n_v)),} \mathrm{𝗇_{\alpha_𝟣} = %$(sstring(n_a)),} \mathrm{𝗇_{\text{𝖻𝗈𝗍𝗁}} = %$(sstring(n_both))}}", padding = (0.0,0.0,14.0,12.0), valign = :center, font=:bold,fontsize = 20)
            end
        end
    end

    Label(f_csd_ts[1,1:end,Top()], "EWS in 100-year high-pass filtered records with 5-year resolution", padding = (0.0,0.0,85.0,25.0), valign = :bottom, font = :bold, fontsize = 20)
    Label(f_csd_ts[2,1:end,Top()], "EWS in 100-year high-pass filtered records with 10-year resolution", padding = (0.0,0.0,85.0,30.0), valign = :bottom, font = :bold, fontsize = 20)
    Label(f_csd_ts[4,1:end,Top()], "EWS in 100-year high-pass filtered records with 20-year resolution", padding = (0.0,0.0,85.0,30.0), valign = :bottom, font = :bold, fontsize = 20)
    
    if legend
        elems = [
            [LineElement(color= (:black, 0.4), points = Point2f[(0, 0.5), (0.5, 0.5)]), LineElement(color= (:black, 0.9), points = Point2f[(0.5, 0.5), (1, 0.5)])],
            [PolyElement(color = (:darkred, 0.1), strokewidth = 0), LineElement(color= (:red,0.4), linewidth=2.0)],# LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:steelblue,0.1), strokewidth = 0), LineElement(color= (:blue, 0.4), linewidth=2.0)],#,LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:darkred, 0.8), strokewidth = 0), LineElement(color= (:red,0.9), linewidth=3.5)],#,LineElement(color = (:darkred,0.9), linewidth= 3.5,points = Point2f[(1, 0), (1, 1)])]]
            ]
        labels = ["EWS indicator", "increasing", "decreasing", "significantly increasing (𝗉<0.05)"]
        Legend(f_csd_ts[end+1,1:end], elems, labels, orientation = :horizontal, valign = :bottom, margin = (10,10,10,10)) 
    end

    colgap!(f_csd_ts.layout, 20)
    rowgap!(f_csd_ts.layout, 0)

    if saving
        save(saveto,f_csd_ts)
    end
    if showing
        display(f_csd_ts)
    end
end

plot_all_csd_ts2_label(ylabs2,lowpass = true, plim = plim, showing = showing, saving = saving, legend = true, saveto = "do_ews_across_greenland_ice_cores/figures/fig6.pdf")
plot_all_csd_ts2_label(ylabs2,lowpass = false, plim = plim, showing = showing, saving = saving, legend = true,  saveto = "do_ews_across_greenland_ice_cores/figures/figS17.pdf")


## again, but aggregated 
function plot_all_csd_ts2_label_aggr(ylabs;lowpass = lowpass, plim = plim, showing = true, saving = false, legend = false, saveto = "paper/all_csd_ts2_label_lowpass_$(lowpass)_p_$(plim).pdf")
    f_csd_ts = Figure(size=(1200,1250))
    letters_h = reshape(letters[1:18],2,9)
    
    whichletter = 0

    function good_csd_files(v,type,lowpass)
        if type == "NGRIP5"
            if lowpass
                return occursin("w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",v) && !occursin("NIKLAS_", v)
            else
                return  occursin("w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",v) && !occursin("NIKLAS_", v)
            end
        elseif type == "10y"
            if lowpass
                return occursin("w_200_normed_filt", v) && occursin("lp_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2", v)
            else
                return occursin("w_200_normed_filt", v) && occursin("gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2", v) && !occursin("_lp_", v)
            end
        elseif type == "20y"
            return occursin("w_200_normed_filt", v) && occursin("gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2", v)
        end
    end

    for (it,type) in enumerate(["NGRIP5", "10y", "20y"])
        ice_cores, var_names, vars, ac_names, acs, sca_names, scas, hurst_names, hursts, sranges = load_data(type)

        good = @. good_csd_files(var_names,type,lowpass)
        
        for (i,(v,a)) in enumerate(zip(reverse(vars[good]), reverse(acs[good])))
            whichletter +=1
            n_v = 0 
            n_a = 0
            n_both = 0
            which_v = []
            for k = 1:2
                gl = f_csd_ts[whichletter,k] = GridLayout()
                
                

                if k == 1
                    ax = Axis(gl[1,1],
                        xticks = (1:17, ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]) ,
                        xreversed = true,
                        xticklabelrotation = pi*0.4,
                        #yreversed=true,
                        ylabel =  L"\mathrm{𝖵}",
                        #xlabel = "Transition",
                        xminorgridwidth = 1,
                        xgridwidth=0,
                        xminorgridcolor = :grey10 ,
                        xgridcolor=:transparent,
                        xminorgridvisible = true, 
                        xgridvisible = false,
                        #yticks = (1:6, ["NGRIP 5 y", "NGRIP 10 y", "NEEM 10 y", "NGRIP 20 y", "GRIP 20 y", "GISP2 20 y"]),
                        #ylabel = "Ice core",
                        #yminorgridwidth = 1,
                        #ygridwidth=0.7,
                        #yminorgridcolor = :grey10 ,
                        #ygridcolor=:grey30,
                        yminorgridvisible = false, 
                        ygridvisible = false,
                        yticksvisible=false,
                        yticklabelsvisible=false,
                        ylabelsize = 18,
                        ylabelfont = :bold,
                        )

                    if whichletter ∈ [1,3,6]
                        ax.xlabel = "Event"
                    else
                        hidexdecorations!(ax, grid=false, minorgrid=false)
                        hidespines!(ax, :b) 
                    end

                    if i>1
                        hidespines!(ax,:t)
                        #hidexdecorations!(ax_event, grid=false)
                    end

                    ma = zeros(17,1)
                    for ev in 1:17
                        #va = 0.4
                        #vl = 2.0
                        #vsc = (:white,1.0)
                        #vc = :black
                        if typeof(v.slopes[ev]) != Nothing
                            if v.slopes[ev] >0 
                                ma[ev,1]+=1
                        #        vc = :red
                        #        vsc = (:darkred,0.1)
                                if v.p_one[ev] < plim
                                    ma[ev,1]+=1
                        #            va = 0.9
                        #            vl = 3.5
                        #            vsc = (:darkred,0.8)
                                    n_v +=1
                                    push!(which_v, ev)
                                end
                            
                            else
                        #        vc = :blue
                        #        vsc = (:steelblue,0.1)
                            end
                        end
                        # vspan!(axv, [GI_onsets[ev]], [GS_onsets[ev]], color = vsc)
                        # pv = Polynomials.fit(v.times[ev][findall(!ismissing, v.vals[ev])], convert.(Float64,v.vals[ev][findall(!ismissing, v.vals[ev])]),1)
                        # lines!(axv,v.times[ev], v.vals[ev], color = :black, alpha = va)#, linewidth = vl)
                        # lines!(axv, v.times[ev], pv.(v.times[ev]), color = vc, alpha = va, linewidth = vl)
                        # vlines!(axv, [GI_onsets[ev]], color = :darkred, alpha = va, linewidth = vl)
                    end
                    co2 = cgrad([:steelblue, :darkred, :darkred], 3, categorical = true, alpha=[0.1,0.1,0.8])
                    hm = heatmap!(ax, ma, colormap=co2, colorrange=(0,2))
                    # xlims!(axv,nothing,nothing)
                    # axv.xreversed = true
                    # xlims!(ax_event, axv.xaxis.attributes.limits[]...)
                elseif k == 2
                    ax = Axis(gl[1,1],
                        xticks = (1:17, ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]) ,
                        xreversed = true,
                        xticklabelrotation = pi*0.4,
                        #yreversed=true,
                        ylabel =  L"\mathrm{\alpha_𝟣}",
                        #xlabel = "Transition",
                        xminorgridwidth = 1,
                        xgridwidth=0,
                        xminorgridcolor = :grey10 ,
                        xgridcolor=:transparent,
                        xminorgridvisible = true, 
                        xgridvisible = false,
                        #yticks = (1:6, ["NGRIP 5 y", "NGRIP 10 y", "NEEM 10 y", "NGRIP 20 y", "GRIP 20 y", "GISP2 20 y"]),
                        #ylabel = "Ice core",
                        yminorgridwidth = 1,
                        ygridwidth=0.7,
                        yminorgridcolor = :grey10 ,
                        ygridcolor=:grey30,
                        yminorgridvisible = false, 
                        ygridvisible = false,
                        yticksvisible=false,
                        yticklabelsvisible=false,
                        ylabelsize = 18,
                        ylabelfont = :bold,
                        )

                    if whichletter ∈ [1,3,6]
                        ax.xlabel = "Event"
                    else
                        hidexdecorations!(ax, grid=false, minorgrid=false)
                        hidespines!(ax, :b) 
                    end
                    #hidedecorations!(ax, label=true)

                    if i>1
                        hidespines!(ax,:t)
                        #hidexdecorations!(ax_event, grid=false)
                    end

                    ma = zeros(17,1)
                    
                    for ev in 1:17
                        #aa = 0.4
                        #al = 2.0
                        #asc = (:white,1.0)
                        #ac = :black
                        if typeof(a.slopes[ev]) != Nothing
                            if a.slopes[ev] >0 
                                #ac = :red
                                #asc = (:darkred,0.1)
                                ma[ev,1]+=1
                                if a.p_one[ev] < plim
                                    ma[ev,1]+=1
                                    #aa = 0.9
                                    #al = 3.5
                                    #asc = (:darkred,0.8)
                                    n_a +=1
                                    if ev in which_v
                                        n_both +=1
                                    end
                                end
                            
                            #else
                            #    ac = :blue
                            #    asc = (:steelblue,0.1)
                            end
                        end
                        # vspan!(axa, [GI_onsets[ev]], [GS_onsets[ev]], color = asc)
                        # pa = Polynomials.fit(a.times[ev][findall(!ismissing, a.vals[ev])], convert.(Float64,a.vals[ev][findall(!ismissing, a.vals[ev])]),1)
                        # lines!(axa,a.times[ev], a.vals[ev], color = :black, alpha = aa)#, linewidth = vl)
                        # lines!(axa, a.times[ev], pa.(a.times[ev]), color = ac, alpha = aa, linewidth = al)
                        # vlines!(axa, [GI_onsets[ev]], color = :darkred, alpha = aa, linewidth = al)

                    end
                    co2 = cgrad([:steelblue, :darkred, :darkred], 3, categorical = true, alpha=[0.1,0.1,0.8])
                    hm = heatmap!(ax, ma, colormap=co2, colorrange=(0,2))
                    #ax.xlabel = "Transition"
                    #xlims!(axa,nothing,nothing)
                    #axa.xreversed = true
                    #xlims!(ax_event, axa.xaxis.attributes.limits[]...)
                    #ax_event.xreversed = true
                end


                Label(gl[1,1,TopLeft()], letters_h[k,whichletter], fontsize = 20,
                            font = :bold, padding = (0,20,-10,15),
                            halign = :left, valign =:bottom)
            end
            if i == 1
                Label(f_csd_ts[whichletter,1:end,Top()], L"\textbf{\textrm{%$(reverse(ylabs2b[it])[i]): }\mathrm{𝗇_{𝖵} = %$(sstring(n_v)),} \mathrm{𝗇_{\alpha_𝟣} = %$(sstring(n_a)),} \mathrm{𝗇_{\text{𝖻𝗈𝗍𝗁}} = %$(sstring(n_both))}}", padding = (0.0,0.0,0.0,30.0), valign = :center, font=:bold,fontsize = 20)
            else
                Label(f_csd_ts[whichletter,1:end,Top()], L"\textbf{\textrm{%$(reverse(ylabs2b[it])[i]): }\mathrm{𝗇_{𝖵} = %$(sstring(n_v)),} \mathrm{𝗇_{\alpha_𝟣} = %$(sstring(n_a)),} \mathrm{𝗇_{\text{𝖻𝗈𝗍𝗁}} = %$(sstring(n_both))}}", padding = (0.0,0.0,14.0,12.0), valign = :center, font=:bold,fontsize = 20)
            end
        end
    end

    Label(f_csd_ts[1,1:end,Top()], "EWS in 100-year high-pass filtered records with 5-year resolution", padding = (0.0,0.0,60.0,25.0), valign = :bottom, font = :bold, fontsize = 20)
    Label(f_csd_ts[2,1:end,Top()], "EWS in 100-year high-pass filtered records with 10-year resolution", padding = (0.0,0.0,60.0,30.0), valign = :bottom, font = :bold, fontsize = 20)
    Label(f_csd_ts[4,1:end,Top()], "EWS in 100-year high-pass filtered records with 20-year resolution", padding = (0.0,0.0,60.0,30.0), valign = :bottom, font = :bold, fontsize = 20)
    
    if legend
        # elems = [
        #     [LineElement(color= (:black, 0.4), points = Point2f[(0, 0.5), (0.5, 0.5)]), LineElement(color= (:black, 0.9), points = Point2f[(0.5, 0.5), (1, 0.5)])],
        #     [PolyElement(color = (:darkred, 0.1), strokewidth = 0), LineElement(color= (:red,0.4), linewidth=2.0)],# LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
        #     [PolyElement(color = (:steelblue,0.1), strokewidth = 0), LineElement(color= (:blue, 0.4), linewidth=2.0)],#,LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
        #     [PolyElement(color = (:darkred, 0.8), strokewidth = 0), LineElement(color= (:red,0.9), linewidth=3.5)],#,LineElement(color = (:darkred,0.9), linewidth= 3.5,points = Point2f[(1, 0), (1, 1)])]]
        #     ]
        # labels = ["EWS indicator", "increasing", "decreasing", "significantly increasing (𝗉<0.05)"]
        elems = [
            [PolyElement(color = (:steelblue,0.1), strokewidth = 0)],#,LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:darkred, 0.1), strokewidth = 0)],# LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:darkred, 0.8), strokewidth = 0)],#,LineElement(color = (:darkred,0.9), linewidth= 3.5,points = Point2f[(1, 0), (1, 1)])]]
            ]
        labels = ["increasing", "decreasing", "significantly increasing (𝗉<0.05)"]
        Legend(f_csd_ts[end+1,1:end], elems, labels, orientation = :horizontal, valign = :bottom, margin = (10,10,10,10)) 
    end

    colgap!(f_csd_ts.layout, 20)
    rowgap!(f_csd_ts.layout, 0)

    if saving
        save(saveto,f_csd_ts)
    end
    if showing
        display(f_csd_ts)
    end
end

plot_all_csd_ts2_label_aggr(ylabs2,lowpass = true, plim = plim, showing = showing, saving = saving, legend = true, saveto = "do_ews_across_greenland_ice_cores/figures/fig6_aggr.pdf")
plot_all_csd_ts2_label_aggr(ylabs2,lowpass = false, plim = plim, showing = showing, saving = saving, legend = true,  saveto = "do_ews_across_greenland_ice_cores/figures/figS17_aggr.pdf")

#Fig 9 and S21
function plot_all_wave_ts2_label(;lowpass = lowpass, 
    smoothw = smoothw, 
    plim = plim, 
    scales = [(20,60),(20,100)],
    showing = true,
    saving = false,
    legend = false,
    savetos = ["paper/all_wave_ts2_label_20_60_lowpass_$(lowpass)_smoothw_$(smoothw)_p_$(plim).pdf",
    "paper/all_wave_ts2_label_20_100_lowpass_$(lowpass)_smoothw_$(smoothw)_p_$(plim).pdf"] )

    function good_wavelet_files(v,type,lowpass, smoothw,s1,s2)
        if type == "NGRIP5"
            if lowpass
                return occursin("s1_$(s1)_s2_$(s2)_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_$(smoothw)_10000_TFTS.jld2",v) && !occursin("NIKLAS_", v)
            else
                return occursin("s1_$(s1)_s2_$(s2)_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_$(smoothw)_10000_TFTS.jld2",v) && !occursin("NIKLAS_", v)
            end
        elseif type == "10y"
            if lowpass
                return occursin("s1_$(s1)_s2_$(s2)", v) && occursin("_lp_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_$(smoothw)_10000_TFTS.jld2", v)
            else
                return occursin("s1_$(s1)_s2_$(s2)", v) && occursin("PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_$(smoothw)_10000_TFTS.jld2", v) && !occursin("_lp_", v)
            end
        elseif type == "20y"
            return occursin("s1_$(s1)_s2_$(s2)", v) && occursin("PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_$(smoothw)_10000_TFTS.jld2", v)
        end
    end

    for (snum,(s1,s2)) in enumerate(scales)
        f_wave_ts = Figure(size=(1200,1250))
        letters_h = reshape(letters[1:18],2,9)

        whichletter = 0

        for (it,type) in enumerate(["NGRIP5", "10y", "20y"])
            ice_cores, var_names, vars, ac_names, acs, sca_names, scas, hurst_names, hursts, sranges = load_data(type)

            goods = @. good_wavelet_files(sca_names,type,lowpass, smoothw,s1,s2)
            goodh = @. good_wavelet_files(hurst_names,type,lowpass, smoothw,s1,s2)


            for (i,(v,a)) in enumerate(zip(reverse(scas[goods]), reverse(hursts[goodh])))
                whichletter +=1
                n_v = 0 
                n_a = 0
                n_both = 0
                which_v = []
                for k = 1:2
                    gl = f_wave_ts[whichletter,k] = GridLayout()
                    if k == 1
                        axv = Axis(gl[1,1], 
                            xlabel = "Age (kyr b2k)", 
                            ylabel = L"\mathrm{\hat{𝗐}^𝟤}", 
                            ylabelsize = 18,
                            ylabelfont = :bold,
                            xticks = (10_000:5_000:60_000, string.(10:5:60)), 
                            )
                        axv.xreversed = true

                        ax_event = Axis(gl[1,1],xaxisposition = :top, xticklabelrotation = pi*0.4)
                        ax_event.xticks = (GI_onsets, event_labels)
                        ax_event.xreversed = true
                        hidespines!(ax_event)
                        hideydecorations!(ax_event)


                        if whichletter ∈ [1,3,6]
                            axv.xlabel = "Age (kyr b2k)"
                        else
                            hidexdecorations!(axv, grid=false)
                            hidespines!(axv, :b) 
                        end

                        if i>1
                            hidespines!(axv,:t)
                            hidexdecorations!(ax_event, grid=false)
                        end
                        for ev in 1:17
                            va = 0.4
                            vl = 2.0
                            vsc = (:white,1.0)
                            vc = :black
                            if typeof(v.slopes[ev]) != Nothing
                                if v.slopes[ev] >0 
                                    vc = :red
                                    vsc = (:darkred,0.1)
                                    if v.p_one[ev] < plim
                                        va = 0.9
                                        vl = 3.5
                                        vsc = (:darkred,0.8)
                                        n_v +=1
                                        push!(which_v, ev)
                                    end
                                
                                else
                                    vc = :blue
                                    vsc = (:steelblue,0.1)
                                end
                            end
                            vspan!(axv, [GI_onsets[ev]], [GS_onsets[ev]], color = vsc)
                            pv = Polynomials.fit(v.times[ev][findall(!ismissing, v.vals[ev])], convert.(Float64,v.vals[ev][findall(!ismissing, v.vals[ev])]),1)
                            lines!(axv,v.times[ev], v.vals[ev], color = :black, alpha = va)#, linewidth = vl)
                            lines!(axv, v.times[ev], pv.(v.times[ev]), color = vc, alpha = va, linewidth = vl)
                            vlines!(axv, [GI_onsets[ev]], color = :darkred, alpha = va, linewidth = vl)
                        end
                        xlims!(axv,nothing,nothing)
                        axv.xreversed = true
                        xlims!(ax_event, axv.xaxis.attributes.limits[]...)
                        ax_event.xreversed = true
                    elseif k == 2
                        axa = Axis(gl[1,1], 
                            ylabel =  L"\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}}",
                            ylabelsize = 18,
                            ylabelfont = :bold,
                            xticks = (10_000:5_000:60_000, string.(10:5:60)),
                            )

                        axa.xreversed = true

                        ax_event = Axis(gl[1,1],xaxisposition = :top, xticklabelrotation = pi*0.4)
                        ax_event.xticks = (GI_onsets, event_labels)
                        ax_event.xreversed = true
                        hidespines!(ax_event)
                        hideydecorations!(ax_event)

                        if whichletter ∈ [1,3,6]
                            axa.xlabel = "Age (kyr b2k)"
                        else
                            hidexdecorations!(axa, grid=false)
                            hidespines!(axa, :b) 
                        end

                        if i>1
                            hidespines!(axa,:t)
                            hidexdecorations!(ax_event, grid=false)
                        end

                        for ev in 1:17      
                            aa = 0.4
                            al = 2.0
                            asc = (:white,1.0)
                            ac = :black
                            if typeof(a.slopes[ev]) != Nothing
                                if a.slopes[ev] >0 
                                    ac = :red
                                    asc = (:darkred,0.1)
                                    if a.p_one[ev] < plim
                                        aa = 0.9
                                        al = 3.5
                                        asc = (:darkred,0.8)
                                        n_a +=1
                                        if ev in which_v
                                            n_both +=1
                                        end
                                    end
                                
                                else
                                    ac = :blue
                                    asc = (:steelblue,0.1)
                                end
                            end
                            vspan!(axa, [GI_onsets[ev]], [GS_onsets[ev]], color = asc)
                            pa = Polynomials.fit(a.times[ev][findall(!ismissing, a.vals[ev])], convert.(Float64,a.vals[ev][findall(!ismissing, a.vals[ev])]),1)
                            lines!(axa,a.times[ev], a.vals[ev], color = :black, alpha = aa)#, linewidth = vl)
                            lines!(axa, a.times[ev], pa.(a.times[ev]), color = ac, alpha = aa, linewidth = al)
                            vlines!(axa, [GI_onsets[ev]], color = :darkred, alpha = aa, linewidth = al)
                        end
                        xlims!(axa,nothing,nothing)
                        axa.xreversed = true
                        xlims!(ax_event, axa.xaxis.attributes.limits[]...)
                        ax_event.xreversed = true
                    end
                    Label(gl[1,1,TopLeft()], letters_h[k,whichletter], fontsize = 20,
                        font = :bold, padding = (0,20,-10,15),
                        halign = :left,valign =:bottom)
                end
                if i == 1
                    Label(f_wave_ts[whichletter,1:end,Top()], L"\textbf{\textrm{%$(reverse(ylabs2b[it])[i]): }\mathrm{𝗇_{\hat{𝗐}^𝟤} = %$(sstring(n_v)),} \mathrm{𝗇_{\hat{𝖧}^\text{𝗅𝗈𝖼}} = %$(sstring(n_a)),} \mathrm{𝗇_{\text{𝖻𝗈𝗍𝗁}} = %$(sstring(n_both))}}", padding = (0.0,0.0,20.0,20.0), valign = :center, font=:bold,fontsize = 20)
                else
                    Label(f_wave_ts[whichletter,1:end,Top()], L"\textbf{\textrm{%$(reverse(ylabs2b[it])[i]): }\mathrm{𝗇_{\hat{𝗐}^𝟤} = %$(sstring(n_v)),} \mathrm{𝗇_{\hat{𝖧}^\text{𝗅𝗈𝖼}} = %$(sstring(n_a)),} \mathrm{𝗇_{\text{𝖻𝗈𝗍𝗁}} = %$(sstring(n_both))}}", padding = (0.0,0.0,14.0,10.0), valign = :center, font=:bold,fontsize = 20)
                end
            end
        end
        Label(f_wave_ts[1,1:end,Top()], "EWS in ($(s1)-$(s2)) year band of records with 5-year resolution", padding = (0.0,0.0,85.0,25.0), valign = :bottom, font = :bold, fontsize = 20)
        Label(f_wave_ts[2,1:end,Top()], "EWS in ($(s1)-$(s2)) year band of records with 10-year resolution", padding = (0.0,0.0,85.0,30.0), valign = :bottom, font = :bold, fontsize = 20)
        Label(f_wave_ts[4,1:end,Top()], "EWS in ($(s1)-$(s2)) year band of records with 20-year resolution", padding = (0.0,0.0,85.0,30.0), valign = :bottom, font = :bold, fontsize = 20)

        if legend
            elems = [
                #[LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(0.5, 0), (0.5, 1)])],
                [LineElement(color= (:black, 0.4), points = Point2f[(0, 0.5), (0.5, 0.5)]), LineElement(color= (:black, 0.9), points = Point2f[(0.5, 0.5), (1, 0.5)])],
                [PolyElement(color = (:darkred, 0.1), strokewidth = 0), LineElement(color= (:red,0.4), linewidth=2.0)],# LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
                [PolyElement(color = (:steelblue,0.1), strokewidth = 0), LineElement(color= (:blue, 0.4), linewidth=2.0)],#,LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
                [PolyElement(color = (:darkred, 0.8), strokewidth = 0), LineElement(color= (:red,0.9), linewidth=3.5)],#,LineElement(color = (:darkred,0.9), linewidth= 3.5,points = Point2f[(1, 0), (1, 1)])]]
                ]
            labels = [#"GI onset", 
                    "EWS indicator", "increasing", "decreasing", "significantly increasing (𝗉<0.05)"]
            Legend(f_wave_ts[end+1,1:end], elems, labels, orientation = :horizontal, valign = :bottom, margin = (10,10,10,10)) 
        end

        colgap!(f_wave_ts.layout, 20)
        rowgap!(f_wave_ts.layout, 0)

        if saving
            save(savetos[snum],f_wave_ts)
        end
        if showing
            display(f_wave_ts)
        end
    end
end


plot_all_wave_ts2_label(lowpass = true, smoothw = false, plim = plim, scales = [(20,60)],showing = showing,saving = saving,legend = true,
                savetos = ["do_ews_across_greenland_ice_cores/figures/fig9.pdf"])
plot_all_wave_ts2_label(lowpass = false, smoothw = false, plim = plim, scales = [(20,60)],showing = showing,saving = saving,legend = true,
                savetos = ["do_ews_across_greenland_ice_cores/figures/figS21.pdf"])

function plot_all_wave_ts2_label_aggr(;lowpass = lowpass, smoothw = smoothw, plim = plim, scales = [(20,60),(20,100)], showing = true, saving = false, legend = false, savetos = ["paper/all_wave_ts2_label_20_60_lowpass_$(lowpass)_smoothw_$(smoothw)_p_$(plim).pdf",
    "paper/all_wave_ts2_label_20_100_lowpass_$(lowpass)_smoothw_$(smoothw)_p_$(plim).pdf"] )
    
    function good_wavelet_files(v,type,lowpass, smoothw,s1,s2)
        if type == "NGRIP5"
            if lowpass
                return occursin("s1_$(s1)_s2_$(s2)_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_$(smoothw)_10000_TFTS.jld2",v) && !occursin("NIKLAS_", v)
            else
                return occursin("s1_$(s1)_s2_$(s2)_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_$(smoothw)_10000_TFTS.jld2",v) && !occursin("NIKLAS_", v)
            end
        elseif type == "10y"
            if lowpass
                return occursin("s1_$(s1)_s2_$(s2)", v) && occursin("_lp_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_$(smoothw)_10000_TFTS.jld2", v)
            else
                return occursin("s1_$(s1)_s2_$(s2)", v) && occursin("PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_$(smoothw)_10000_TFTS.jld2", v) && !occursin("_lp_", v)
            end
        elseif type == "20y"
            return occursin("s1_$(s1)_s2_$(s2)", v) && occursin("PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_$(smoothw)_10000_TFTS.jld2", v)
        end
    end

    for (snum,(s1,s2)) in enumerate(scales)
        f_wave_ts = Figure(size=(1200,1250))
        letters_h = reshape(letters[1:18],2,9)
        
        whichletter = 0


        for (it,type) in enumerate(["NGRIP5", "10y", "20y"])
            ice_cores, var_names, vars, ac_names, acs, sca_names, scas, hurst_names, hursts, sranges = load_data(type)

            goods = @. good_wavelet_files(sca_names,type,lowpass, smoothw,s1,s2)
            goodh = @. good_wavelet_files(hurst_names,type,lowpass, smoothw,s1,s2)
            
            for (i,(v,a)) in enumerate(zip(reverse(scas[goods]), reverse(hursts[goodh])))
                whichletter +=1
                n_v = 0 
                n_a = 0
                n_both = 0
                which_v = []
                for k = 1:2
                    gl = f_wave_ts[whichletter,k] = GridLayout()
        
                    if k == 1
                        ax = Axis(gl[1,1],
                            xticks = (1:17, ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]) ,
                            xreversed = true,
                            xticklabelrotation = pi*0.4,
                            #yreversed=true,
                            ylabel =  L"\mathrm{\hat{𝗐}^𝟤}",
                            #xlabel = "Transition",
                            xminorgridwidth = 1,
                            xgridwidth=0,
                            xminorgridcolor = :grey10 ,
                            xgridcolor=:transparent,
                            xminorgridvisible = true, 
                            xgridvisible = false,
                            #yticks = (1:6, ["NGRIP 5 y", "NGRIP 10 y", "NEEM 10 y", "NGRIP 20 y", "GRIP 20 y", "GISP2 20 y"]),
                            #ylabel = "Ice core",
                            #yminorgridwidth = 1,
                            #ygridwidth=0.7,
                            #yminorgridcolor = :grey10 ,
                            #ygridcolor=:grey30,
                            yminorgridvisible = false, 
                            ygridvisible = false,
                            yticksvisible=false,
                            yticklabelsvisible=false,
                            ylabelsize = 18,
                            ylabelfont = :bold,
                            )

                        if whichletter ∈ [1,3,6]
                            ax.xlabel = "Event"
                        else
                            hidexdecorations!(ax, grid=false, minorgrid=false)
                            hidespines!(ax, :b) 
                        end

                        if i>1
                            hidespines!(ax,:t)
                            #hidexdecorations!(ax_event, grid=false)
                        end

                        ma = zeros(17,1)
                        for ev in 1:17
                            #va = 0.4
                            #vl = 2.0
                            #vsc = (:white,1.0)
                            #vc = :black
                            if typeof(v.slopes[ev]) != Nothing
                                if v.slopes[ev] >0 
                                    ma[ev,1]+=1
                            #        vc = :red
                            #        vsc = (:darkred,0.1)
                                    if v.p_one[ev] < plim
                                        ma[ev,1]+=1
                            #            va = 0.9
                            #            vl = 3.5
                            #            vsc = (:darkred,0.8)
                                        n_v +=1
                                        push!(which_v, ev)
                                    end
                                
                                else
                            #        vc = :blue
                            #        vsc = (:steelblue,0.1)
                                end
                            end
                            # vspan!(axv, [GI_onsets[ev]], [GS_onsets[ev]], color = vsc)
                            # pv = Polynomials.fit(v.times[ev][findall(!ismissing, v.vals[ev])], convert.(Float64,v.vals[ev][findall(!ismissing, v.vals[ev])]),1)
                            # lines!(axv,v.times[ev], v.vals[ev], color = :black, alpha = va)#, linewidth = vl)
                            # lines!(axv, v.times[ev], pv.(v.times[ev]), color = vc, alpha = va, linewidth = vl)
                            # vlines!(axv, [GI_onsets[ev]], color = :darkred, alpha = va, linewidth = vl)
                        end
                        co2 = cgrad([:steelblue, :darkred, :darkred], 3, categorical = true, alpha=[0.1,0.1,0.8])
                        hm = heatmap!(ax, ma, colormap=co2, colorrange=(0,2))
                        # xlims!(axv,nothing,nothing)
                        # axv.xreversed = true
                        # xlims!(ax_event, axv.xaxis.attributes.limits[]...)
                    elseif k == 2
                        ax = Axis(gl[1,1],
                            xticks = (1:17, ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]) ,
                            xreversed = true,
                            xticklabelrotation = pi*0.4,
                            #yreversed=true,
                            ylabel =  L"\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}}",
                            #xlabel = "Transition",
                            xminorgridwidth = 1,
                            xgridwidth=0,
                            xminorgridcolor = :grey10 ,
                            xgridcolor=:transparent,
                            xminorgridvisible = true, 
                            xgridvisible = false,
                            #yticks = (1:6, ["NGRIP 5 y", "NGRIP 10 y", "NEEM 10 y", "NGRIP 20 y", "GRIP 20 y", "GISP2 20 y"]),
                            #ylabel = "Ice core",
                            yminorgridwidth = 1,
                            ygridwidth=0.7,
                            yminorgridcolor = :grey10 ,
                            ygridcolor=:grey30,
                            yminorgridvisible = false, 
                            ygridvisible = false,
                            yticksvisible=false,
                            yticklabelsvisible=false,
                            ylabelsize = 18,
                            ylabelfont = :bold,
                            )

                        if whichletter ∈ [1,3,6]
                            ax.xlabel = "Event"
                        else
                            hidexdecorations!(ax, grid=false, minorgrid=false)
                            hidespines!(ax, :b) 
                        end
                        #hidedecorations!(ax, label=true)

                        if i>1
                            hidespines!(ax,:t)
                            #hidexdecorations!(ax_event, grid=false)
                        end

                        ma = zeros(17,1)
                        
                        for ev in 1:17
                            #aa = 0.4
                            #al = 2.0
                            #asc = (:white,1.0)
                            #ac = :black
                            if typeof(a.slopes[ev]) != Nothing
                                if a.slopes[ev] >0 
                                    #ac = :red
                                    #asc = (:darkred,0.1)
                                    ma[ev,1]+=1
                                    if a.p_one[ev] < plim
                                        ma[ev,1]+=1
                                        #aa = 0.9
                                        #al = 3.5
                                        #asc = (:darkred,0.8)
                                        n_a +=1
                                        if ev in which_v
                                            n_both +=1
                                        end
                                    end
                                
                                #else
                                #    ac = :blue
                                #    asc = (:steelblue,0.1)
                                end
                            end
                            # vspan!(axa, [GI_onsets[ev]], [GS_onsets[ev]], color = asc)
                            # pa = Polynomials.fit(a.times[ev][findall(!ismissing, a.vals[ev])], convert.(Float64,a.vals[ev][findall(!ismissing, a.vals[ev])]),1)
                            # lines!(axa,a.times[ev], a.vals[ev], color = :black, alpha = aa)#, linewidth = vl)
                            # lines!(axa, a.times[ev], pa.(a.times[ev]), color = ac, alpha = aa, linewidth = al)
                            # vlines!(axa, [GI_onsets[ev]], color = :darkred, alpha = aa, linewidth = al)

                        end
                        co2 = cgrad([:steelblue, :darkred, :darkred], 3, categorical = true, alpha=[0.1,0.1,0.8])
                        hm = heatmap!(ax, ma, colormap=co2, colorrange=(0,2))
                        #ax.xlabel = "Transition"
                        #xlims!(axa,nothing,nothing)
                        #axa.xreversed = true
                        #xlims!(ax_event, axa.xaxis.attributes.limits[]...)
                        #ax_event.xreversed = true
                    end


                    Label(gl[1,1,TopLeft()], letters_h[k,whichletter], fontsize = 20,
                                font = :bold, padding = (0,20,-10,15),
                                halign = :left, valign =:bottom)
                end
                if i == 1
                    Label(f_wave_ts[whichletter,1:end,Top()], L"\textbf{\textrm{%$(reverse(ylabs2b[it])[i]): }\mathrm{𝗇_{\hat{𝗐}^𝟤} = %$(sstring(n_v)),} \mathrm{𝗇_{\hat{𝖧}^\text{𝗅𝗈𝖼}} = %$(sstring(n_a)),} \mathrm{𝗇_{\text{𝖻𝗈𝗍𝗁}} = %$(sstring(n_both))}}", padding = (0.0,0.0,0.0,30.0), valign = :center, font=:bold,fontsize = 20)
                else
                    Label(f_wave_ts[whichletter,1:end,Top()], L"\textbf{\textrm{%$(reverse(ylabs2b[it])[i]): }\mathrm{𝗇_{\hat{𝗐}^𝟤} = %$(sstring(n_v)),} \mathrm{𝗇_{\hat{𝖧}^\text{𝗅𝗈𝖼}} = %$(sstring(n_a)),} \mathrm{𝗇_{\text{𝖻𝗈𝗍𝗁}} = %$(sstring(n_both))}}", padding = (0.0,0.0,14.0,12.0), valign = :center, font=:bold,fontsize = 20)
                end
            end
        end

        Label(f_wave_ts[1,1:end,Top()], "EWS in ($(s1)-$(s2)) year band of records with 5-year resolution", padding = (0.0,0.0,60.0,25.0), valign = :bottom, font = :bold, fontsize = 20)
        Label(f_wave_ts[2,1:end,Top()], "EWS in ($(s1)-$(s2)) year band of records with 10-year resolution", padding = (0.0,0.0,60.0,30.0), valign = :bottom, font = :bold, fontsize = 20)
        Label(f_wave_ts[4,1:end,Top()], "EWS in ($(s1)-$(s2)) year band of records with 20-year resolution", padding = (0.0,0.0,60.0,30.0), valign = :bottom, font = :bold, fontsize = 20)
        
        if legend
            # elems = [
            #     [LineElement(color= (:black, 0.4), points = Point2f[(0, 0.5), (0.5, 0.5)]), LineElement(color= (:black, 0.9), points = Point2f[(0.5, 0.5), (1, 0.5)])],
            #     [PolyElement(color = (:darkred, 0.1), strokewidth = 0), LineElement(color= (:red,0.4), linewidth=2.0)],# LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            #     [PolyElement(color = (:steelblue,0.1), strokewidth = 0), LineElement(color= (:blue, 0.4), linewidth=2.0)],#,LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            #     [PolyElement(color = (:darkred, 0.8), strokewidth = 0), LineElement(color= (:red,0.9), linewidth=3.5)],#,LineElement(color = (:darkred,0.9), linewidth= 3.5,points = Point2f[(1, 0), (1, 1)])]]
            #     ]
            # labels = ["EWS indicator", "increasing", "decreasing", "significantly increasing (𝗉<0.05)"]
            elems = [
                [PolyElement(color = (:steelblue,0.1), strokewidth = 0)],#,LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
                [PolyElement(color = (:darkred, 0.1), strokewidth = 0)],# LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
                [PolyElement(color = (:darkred, 0.8), strokewidth = 0)],#,LineElement(color = (:darkred,0.9), linewidth= 3.5,points = Point2f[(1, 0), (1, 1)])]]
                ]
            labels = ["increasing", "decreasing", "significantly increasing (𝗉<0.05)"]
            Legend(f_wave_ts[end+1,1:end], elems, labels, orientation = :horizontal, valign = :bottom, margin = (10,10,10,10)) 
        end

        colgap!(f_wave_ts.layout, 20)
        rowgap!(f_wave_ts.layout, 0)

        if saving
            save(savetos[snum],f_wave_ts)
        end
        if showing
            display(f_wave_ts)
        end
    end
end

plot_all_wave_ts2_label_aggr(lowpass = true, smoothw = false, plim = plim, scales = [(20,60)],showing = showing,saving = saving,legend = true,
                savetos = ["do_ews_across_greenland_ice_cores/figures/fig9_aggr.pdf"])
plot_all_wave_ts2_label_aggr(lowpass = false, smoothw = false, plim = plim, scales = [(20,60)],showing = showing,saving = saving,legend = true,
                savetos = ["do_ews_across_greenland_ice_cores/figures/figS21_aggr.pdf"])


#Fig S12
function plot_compare_methods_ngrip_csd2(; method = "TFTS", include_filtering = false,
                                        lowpass = lowpass, 
                                        plim = plim, 
                                        py2=true,
                                        showing = true,
                                        saving = false,
                                        legend = false,
                                        saveto ="paper/compare_methods_csd2_ngrip5_$(method)_inclFilt_$(include_filtering)_lowpass_$(lowpass)_p_$(plim).pdf")
    #CSD comparison with all the steps
    type = "NGRIP5"
    _, var_names, vars, _, acs, _ = load_data(type)

    file1(v) = findfirst(==("NIKLAS_w_200_N_lowpass_FiltInd_true_onlyfull_false_detrend_GLSAR_DemeanSur_true_10000_FOURIER.jld2"),v)
    file2(v) = findfirst(==("w_200_normed_filt_N_lowpass_gs_short_FiltInd_true_onlyfull_false_10000_$method.jld2"),v)
    
    file3(v) = findfirst(==("w_200_normed_filt_N_lowpass_gs_short_FiltInd_false_onlyfull_false_10000_$method.jld2"),v)
    
    file4(v) = findfirst(==("w_200_normed_filt_N_lowpass_gs_FiltInd_false_onlyfull_true_10000_$method.jld2"),v)
    
    file5(v) = findfirst(==("w_200_normed_filt_N_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_$method.jld2"),v)
    file6(v) = findfirst(==("w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_$method.jld2"),v)

    file7(v) = findfirst(==("w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_$method.jld2"),v)

    file1p(v) = findfirst(==("NIKLAS_w_200_N_lowpass_py2_FiltInd_true_onlyfull_false_detrend_GLSAR_DemeanSur_true_10000_FOURIER.jld2"),v)
    file2p(v) = findfirst(==("w_200_normed_filt_N_lowpass_py2_gs_short_FiltInd_true_onlyfull_false_10000_TFTS.jld2"),v)
    file3p(v) = findfirst(==("w_200_normed_filt_N_lowpass_py2_gs_short_FiltInd_false_onlyfull_false_10000_TFTS.jld2"),v)
    file4p(v) = findfirst(==("w_200_normed_filt_N_lowpass_py2_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2"),v)
    
    
    
    if include_filtering
        if lowpass
            file_list = [file1, file2, file3, file4, file7]
            version_names = ["Boers, 2018", "TFTS surrogates on data (Step 1)", "No indicator filtering (Step 2.1a)", "Entire GS and only windows with 200 y (Step 2.2a)", "Raw ages not rounded (Step 3.2)"]
        else
            file_list = [file1, file2, file3, file4, file5, file6]
            version_names = ["Boers, 2018", "TFTS surrogates on data (Step 1)", "No indicator filtering (Step 2.1a)", "Entire GS and only windows with 200 y (Step 2.2a)","No low-pass filtering after interpolation","Raw ages not rounded (Step 3.2)"]
        end
        if py2
            if lowpass
                file_list = [file1p,file2p,file3p,file4p,file4,file7]
                version_names = ["Boers, 2018", "TFTS surrogates on data (Step 1)", "No indicator filtering (Step 2.1a)", "Entire GS and only windows with 200 y (Step 2.2a)" , "Preprocessing in Julia (Step 3.1)", "Raw ages not rounded (Step 3.2)"]
            else
                file_list = [file1p,file2p,file3p,file4p,file4, file5, file6]
                version_names = ["Boers, 2018", "TFTS surrogates on data (Step 1)", "No indicator filtering (Step 2.1a)", "Entire GS and only windows with 200 y (Step 2.2a)", "Preprocessing in Julia (Step 3.1)", "No low-pass filtering after interpolation", "Raw ages not rounded (Step 3.2)"]
            end
        end
    else
        if lowpass
            file_list = [file1, file2, file4, file7]
            version_names = ["Boers, 2018", "TFTS surrogates on data (Step 1)", "Entire GS and only windows with 200 y (Step 2a)","Raw ages not rounded (Step 3.2)"]
        else
            file_list = [file1, file2, file4, file5, file6]
            version_names = ["Boers, 2018", "TFTS surrogates on data (Step 1)", "Entire GS and only windows with 200  y (Step 2a)","No low-pass filtering after interpolation","Raw ages not rounded (Step 3.2)"]
        end
        if py2
            if lowpass
                file_list = [file1p, file2p, file4p, file4, file7]
                version_names = ["Boers, 2018", "TFTS surrogates on data (Step 1)", "Entire GS and only windows with 200 y (Step 2a)", "Preprocessing in Julia (Step 3.1)", "Raw ages not rounded (Step 3.2)"]
            else
                file_list = [file1p, file2p, file4p, file4, file5, file6]
                version_names = ["Boers, 2018", "TFTS surrogates on data (Step 1)", "Entire GS and only windows with 200 y (Step 2a)", "Preprocessing in Julia (Step 3.1)", "No low-pass filtering after interpolation", "Raw ages not rounded (Step 3.2)"]
            end
        end
    end

    fts = Figure(size= (1600,250*length(file_list)))
    letters_h = reshape(letters[1:18],2,9)
    
    whichletter = 0

    for (i,ff) in enumerate(file_list)
        whichletter +=1
        
        good = ff(var_names)
        v = vars[good]
        a = acs[good]

        n_v = 0 
        n_a = 0
        n_both = 0
        which_v = []

        for k = 1:2
            gl = fts[whichletter,k] = GridLayout()
            if k == 1
                axv = Axis(gl[1,1], 
                    ylabel = L"\mathrm{𝖵}", 
                    ylabelsize = 18,
                    ylabelfont = :bold,
                    xticks = (10_000:5_000:60_000, string.(10:5:60)), 
                    )
                axv.xreversed = true
                if whichletter == length(file_list)
                    axv.xlabel = "Age (kyr b2k)"
                else
                    hidexdecorations!(axv, grid=false)
                    hidespines!(axv, :b) 
                end

                if i>1
                    hidespines!(axv,:t)
                end

                for ev in 1:17
                    va = 0.4
                    vl = 2.0
                    vsc = (:white,1.0)
                    vc = :black
                    if typeof(v.slopes[ev]) != Nothing
                        if v.slopes[ev] >0 
                            vc = :red
                            vsc = (:darkred,0.1)
                            if v.p_one[ev] < plim
                                va = 0.9
                                #vas[ev] = 1.0
                                vl = 3.5
                                vsc = (:darkred,0.8)
                                n_v +=1
                                push!(which_v, ev)
                            end
                        
                        else
                            vc = :blue
                            vsc = (:steelblue,0.1)
                        end
                    end
                
                    vspan!(axv, [GI_onsets[ev]], [GS_onsets[ev]], color = vsc)
                    pv = Polynomials.fit(v.times[ev][findall(!ismissing, v.vals[ev])], convert.(Float64,v.vals[ev][findall(!ismissing, v.vals[ev])]),1)
                    lines!(axv,v.times[ev], v.vals[ev], color = :black, alpha = va)#, linewidth = vl)
                    lines!(axv, v.times[ev], pv.(v.times[ev]), color = vc, alpha = va, linewidth = vl)
                    vlines!(axv, [GI_onsets[ev]], color = :darkred, alpha = va, linewidth = vl)
                end
            elseif k == 2
                axa = Axis(gl[1,1], 
                        ylabel =  L"\mathrm{\alpha_𝟣}",
                        ylabelsize = 18,
                        ylabelfont = :bold,
                        xticks = (10_000:5_000:60_000, string.(10:5:60)), 
                        )
                    
                axa.xreversed = true

                if whichletter == length(file_list)
                    axa.xlabel = "Age (kyr b2k)"
                else
                    hidexdecorations!(axa, grid=false)
                    hidespines!(axa, :b) 
                end

                if i>1
                    hidespines!(axa,:t)
                end

                for ev in 1:17
                    aa = 0.4
                    al = 2.0
                    asc = (:white,1.0)
                    ac = :black
                    if typeof(a.slopes[ev]) != Nothing
                        if a.slopes[ev] >0 
                            ac = :red
                            asc = (:darkred,0.1)
                            if a.p_one[ev] < plim
                                aa = 0.9
                                al = 3.5
                                asc = (:darkred,0.8)
                                n_a +=1
                                if ev in which_v
                                    n_both +=1
                                end
                            end
                        
                        else
                            ac = :blue
                            asc = (:steelblue,0.1)
                        end
                    end
                    vspan!(axa, [GI_onsets[ev]], [GS_onsets[ev]], color = asc)
                    pa = Polynomials.fit(a.times[ev][findall(!ismissing, a.vals[ev])], convert.(Float64,a.vals[ev][findall(!ismissing, a.vals[ev])]),1)
                    lines!(axa,a.times[ev], a.vals[ev], color = :black, alpha = aa)#, linewidth = vl)
                    lines!(axa, a.times[ev], pa.(a.times[ev]), color = ac, alpha = aa, linewidth = al)
                    vlines!(axa, [GI_onsets[ev]], color = :darkred, alpha = aa, linewidth = al)
                end
            end
            Label(gl[1,1,TopLeft()], letters_h[k,whichletter], fontsize = 20,
                            font = :bold, padding = (0,20,-10,15),#(0,5,5,0),
                            halign = :left,valign =:bottom)
        end
        Label(fts[i,1:end,Top()], L"\textbf{\textrm{%$(sstring(version_names[i])): } \mathrm{𝗇_{𝖵} = %$(sstring(n_v)), } \mathrm{𝗇_{\alpha_𝟣} = %$(sstring(n_a)), } \mathrm{𝗇_{𝖻𝗈𝗍𝗁} = %$(sstring(n_both))}}" ,
                     padding = (0.0,0.0,20.0,20.0), valign = :bottom, font = :bold, fontsize = 20)
    end
    if legend
        elems = [
            #[LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(0.5, 0), (0.5, 1)])],
            [LineElement(color= (:black, 0.4), points = Point2f[(0, 0.5), (0.5, 0.5)]), LineElement(color= (:black, 0.9), points = Point2f[(0.5, 0.5), (1, 0.5)])],
            [PolyElement(color = (:darkred, 0.1), strokewidth = 0), LineElement(color= (:red,0.4), linewidth=2.0)],# LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:steelblue,0.1), strokewidth = 0), LineElement(color= (:blue, 0.4), linewidth=2.0)],#,LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:darkred, 0.8), strokewidth = 0), LineElement(color= (:red,0.9), linewidth=3.5)],#,LineElement(color = (:darkred,0.9), linewidth= 3.5,points = Point2f[(1, 0), (1, 1)])]]
            ]
        labels = [#"GI onset", 
                "EWS indicator", "increasing", "decreasing", "significantly increasing (𝗉<0.05)"]
        Legend(fts[end+1,1:end], elems, labels, orientation = :horizontal, valign = :bottom, margin = (10,10,10,10)) 
    end

    colgap!(fts.layout, 20)
    rowgap!(fts.layout, 0)
    Label(fts[1,1:end,Top()], "EWS in 100-year high-pass filtered NGRIP record with 5-year resolution", 
    padding = (0.0,0.0,60.0,20.0), valign = :top, font = :bold, fontsize = 20)
    if saving
        save(saveto,fts)
    end
    if showing
        display(fts)
    end
end

plot_compare_methods_ngrip_csd2(method = "TFTS", include_filtering = true,lowpass = true, plim = plim, py2=true,showing = showing,saving = saving,
                        legend=true, saveto ="do_ews_across_greenland_ice_cores/figures/figS12.pdf")



# Fig 4 and S15
function plot_compare_methods_shorter_ngrip_csd2_label(; method = "TFTS",
        lowpass = lowpass, 
        plim = plim, 
        py2=true,
        showing = true,
        saving = false,
        label_w_steps = true,
        legend = true,
        saveto ="paper/compare_methods_csd_shorter2_ngrip5_$(method)_lowpass_$(lowpass)_p_$(plim)_stepnumber_$(label_w_steps).pdf")
    #CSD comparison with all the steps
    type = "NGRIP5"
    _, var_names, vars, _, acs, _ = load_data(type)

    file1(v) = findfirst(==("NIKLAS_w_200_N_lowpass_FiltInd_true_onlyfull_false_detrend_GLSAR_DemeanSur_true_10000_FOURIER.jld2"),v)
    file2(v) = findfirst(==("w_200_normed_filt_N_lowpass_gs_short_FiltInd_true_onlyfull_false_10000_$method.jld2"),v)

    file4(v) = findfirst(==("w_200_normed_filt_N_lowpass_gs_FiltInd_false_onlyfull_true_10000_$method.jld2"),v)

    file6(v) = findfirst(==("w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_$method.jld2"),v)

    file7(v) = findfirst(==("w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_$method.jld2"),v)

    file1p(v) = findfirst(==("NIKLAS_w_200_N_lowpass_py2_FiltInd_true_onlyfull_false_detrend_GLSAR_DemeanSur_true_10000_FOURIER.jld2"),v)
    file2p(v) = findfirst(==("w_200_normed_filt_N_lowpass_py2_gs_short_FiltInd_true_onlyfull_false_10000_TFTS.jld2"),v)
    file3p(v) = findfirst(==("w_200_normed_filt_N_lowpass_py2_gs_short_FiltInd_false_onlyfull_false_10000_TFTS.jld2"),v)
    file4p(v) = findfirst(==("w_200_normed_filt_N_lowpass_py2_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2"),v)


    if lowpass
        file_list = [file1, file2, file4, file7]  
    else
        file_list = [file1, file2, file4, file6]
    end

    if py2
        if lowpass
            file_list = [file1p, file2p, file4p, file7]
        else
            file_list = [file1p, file2p, file4p, file6]
        end
    end
    version_names = ["Boers, 2018", "Modified significance testing", "Modified EWS calculation", "Modified data preprocessing"]
    if label_w_steps == true
        version_names = ["Boers, 2018", "Modified significance testing (Step 1)", "Modified EWS calculation (Step 2a)", "Modified data preprocessing (Step 3)"]
    end

    fts = Figure(size=(1200, 300*length(file_list)))
    letters_h = reshape(letters[1:18],2,9)
    whichletter = 0

    for (i,ff) in enumerate(file_list)
        whichletter +=1

        good = ff(var_names)
        v = vars[good]
        a = acs[good]

        n_v = 0 
        n_a = 0
        n_both = 0
        which_v = []

        for k = 1:2
            gl = fts[whichletter,k] = GridLayout()
            if k == 1
                axv = Axis(gl[1,1], 
                    ylabel = L"\mathrm{𝖵}", 
                    ylabelsize = 18,
                    ylabelfont = :bold,
                    xticks = (10_000:5_000:60_000, string.(10:5:60)), 
                    )
                axv.xreversed = true

                ax_event = Axis(gl[1,1],xaxisposition = :top, xticklabelrotation = pi*0.4)
                ax_event.xticks = (GI_onsets, event_labels)
                ax_event.xreversed = true
                hidespines!(ax_event)
                hideydecorations!(ax_event)

                if whichletter == length(file_list)
                    axv.xlabel = "Age (kyr b2k)"
                else
                    hidexdecorations!(axv, grid=false)
                    hidespines!(axv, :b) 
                end

                if i>1
                    hidespines!(axv,:t)
                    hidexdecorations!(ax_event, grid=false)
                end
                for ev in 1:17
                    va = 0.4
                    vl = 2.0
                    vsc = (:white,1.0)
                    vc = :black
                    if typeof(v.slopes[ev]) != Nothing
                        if v.slopes[ev] >0 
                            vc = :red
                            vsc = (:darkred,0.1)
                            if v.p_one[ev] < plim
                                va = 0.9
                                #vas[ev] = 1.0
                                vl = 3.5
                                vsc = (:darkred,0.8)
                                n_v +=1
                                push!(which_v, ev)
                            end
                        else
                            vc = :blue
                            vsc = (:steelblue,0.1)
                        end
                    end

                    vspan!(axv, [GI_onsets[ev]], [GS_onsets[ev]], color = vsc, label = "significantly increasing")
                    pv = Polynomials.fit(v.times[ev][findall(!ismissing, v.vals[ev])], convert.(Float64,v.vals[ev][findall(!ismissing, v.vals[ev])]),1)
                    lines!(axv,v.times[ev], v.vals[ev], color = :black, alpha = va, label = "variance")#, linewidth = vl)
                    lines!(axv, v.times[ev], pv.(v.times[ev]), color = vc, alpha = va, linewidth = vl, label = "slope")
                    vlines!(axv, [GI_onsets[ev]], color = :darkred, alpha = va, linewidth = vl, label = "GI onset")
                end
                xlims!(axv,nothing,nothing)
                axv.xreversed = true
                xlims!(ax_event, axv.xaxis.attributes.limits[]...)
                ax_event.xreversed = true
            elseif k == 2
                axa = Axis(gl[1,1], 
                        ylabel =  L"\mathrm{\alpha_𝟣}",
                        ylabelsize = 18,
                        ylabelfont = :bold,
                        xticks = (10_000:5_000:60_000, string.(10:5:60)), 
                        )

                axa.xreversed = true

                ax_event = Axis(gl[1,1],xaxisposition = :top, xticklabelrotation = pi*0.4)
                ax_event.xticks = (GI_onsets, event_labels)
                ax_event.xreversed = true
                hidespines!(ax_event)
                hideydecorations!(ax_event)

                if whichletter == length(file_list)
                    axa.xlabel = "Age (kyr b2k)"
                else
                    hidexdecorations!(axa, grid=false)
                    hidespines!(axa, :b) 
                end

                if i>1
                    hidespines!(axa,:t)
                    hidexdecorations!(ax_event, grid=false)
                end

                for ev in 1:17
                    aa = 0.4
                    al = 2.0
                    asc = (:white,1.0)
                    ac = :black
                    if typeof(a.slopes[ev]) != Nothing
                        if a.slopes[ev] >0 
                            ac = :red
                            asc = (:darkred,0.1)
                            if a.p_one[ev] < plim
                                aa = 0.9
                                al = 3.5
                                asc = (:darkred,0.8)
                                n_a +=1
                                if ev in which_v
                                    n_both +=1
                                end
                            end

                        else
                            ac = :blue
                            asc = (:steelblue,0.1)
                            
                        end
                    end
                    vspan!(axa, [GI_onsets[ev]], [GS_onsets[ev]], color = asc)
                    pa = Polynomials.fit(a.times[ev][findall(!ismissing, a.vals[ev])], convert.(Float64,a.vals[ev][findall(!ismissing, a.vals[ev])]),1)
                    lines!(axa,a.times[ev], a.vals[ev], color = :black, alpha = aa)#, linewidth = vl)
                    lines!(axa, a.times[ev], pa.(a.times[ev]), color = ac, alpha = aa, linewidth = al)
                    vlines!(axa, [GI_onsets[ev]], color = :darkred, alpha = aa, linewidth = al)
                end
                xlims!(axa,nothing,nothing)
                axa.xreversed = true
                xlims!(ax_event, axa.xaxis.attributes.limits[]...)
                ax_event.xreversed = true
            end
            Label(gl[1,1,TopLeft()], letters_h[k,whichletter], fontsize = 20,
                font = :bold, padding = (0,20,-10,15),
                halign = :left,valign =:bottom)
        end
        Label(fts[i,1:end,Top()], L"\textbf{\textrm{%$(sstring(version_names[i])): } \mathrm{𝗇_{𝖵} = %$(sstring(n_v)), } \mathrm{𝗇_{\alpha_𝟣} = %$(sstring(n_a)), } \mathrm{𝗇_{𝖻𝗈𝗍𝗁} = %$(sstring(n_both))}}" ,
                padding = (0.0,0.0,20.0,20.0), valign = :center, font = :bold, fontsize = 20)

    end
    if legend
        elems = [
            #[LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(0.5, 0), (0.5, 1)])],
            [LineElement(color= (:black, 0.4), points = Point2f[(0, 0.5), (0.5, 0.5)]), LineElement(color= (:black, 0.9), points = Point2f[(0.5, 0.5), (1, 0.5)])],
            [PolyElement(color = (:darkred, 0.1), strokewidth = 0), LineElement(color= (:red,0.4), linewidth=2.0)],# LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:steelblue,0.1), strokewidth = 0), LineElement(color= (:blue, 0.4), linewidth=2.0)],#,LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:darkred, 0.8), strokewidth = 0), LineElement(color= (:red,0.9), linewidth=3.5)],#,LineElement(color = (:darkred,0.9), linewidth= 3.5,points = Point2f[(1, 0), (1, 1)])]]
            ]
        labels = [#"GI onset", 
                "EWS indicator", "increasing", "decreasing", "significantly increasing (𝗉<0.05)"]
        Legend(fts[end+1,1:end], elems, labels, orientation = :horizontal, valign = :bottom, margin = (10,10,10,10)) 
    end

    colgap!(fts.layout, 20)
    rowgap!(fts.layout, 0)
    Label(fts[1,1:end,Top()], "EWS in 100-year high-pass filtered NGRIP record with 5-year resolution", 
        padding = (0.0,0.0,110.0,3.0), valign = :bottom, font = :bold, fontsize = 20)
    if saving
        save(saveto,fts)
    end
    if showing
        display(fts)
    end
end

plot_compare_methods_shorter_ngrip_csd2_label(method = "TFTS",lowpass = true, plim = plim, label_w_steps=true,py2=true,legend = true,
            showing = showing, saving = saving, saveto = "do_ews_across_greenland_ice_cores/figures/fig4.pdf")
plot_compare_methods_shorter_ngrip_csd2_label(method = "TFTS",lowpass = false, plim = plim, label_w_steps=true,py2=true,legend = true,
            showing = showing, saving = saving, saveto = "do_ews_across_greenland_ice_cores/figures/figS15.pdf")


function plot_compare_methods_shorter_ngrip_csd2_label_aggr(; method = "TFTS",
            lowpass = lowpass, 
            plim = plim, 
            py2=true,
            showing = true,
            saving = false,
            label_w_steps = true,
            legend = true,
            saveto ="paper/compare_methods_csd_shorter2_aggr_ngrip5_$(method)_lowpass_$(lowpass)_p_$(plim)_stepnumber_$(label_w_steps).pdf")
    
    type = "NGRIP5"
    _, var_names, vars, _, acs, _ = load_data(type)

    file1(v) = findfirst(==("NIKLAS_w_200_N_lowpass_FiltInd_true_onlyfull_false_detrend_GLSAR_DemeanSur_true_10000_FOURIER.jld2"),v)
    file2(v) = findfirst(==("w_200_normed_filt_N_lowpass_gs_short_FiltInd_true_onlyfull_false_10000_$method.jld2"),v)

    file4(v) = findfirst(==("w_200_normed_filt_N_lowpass_gs_FiltInd_false_onlyfull_true_10000_$method.jld2"),v)

    file6(v) = findfirst(==("w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_$method.jld2"),v)

    file7(v) = findfirst(==("w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_$method.jld2"),v)

    file1p(v) = findfirst(==("NIKLAS_w_200_N_lowpass_py2_FiltInd_true_onlyfull_false_detrend_GLSAR_DemeanSur_true_10000_FOURIER.jld2"),v)
    file2p(v) = findfirst(==("w_200_normed_filt_N_lowpass_py2_gs_short_FiltInd_true_onlyfull_false_10000_TFTS.jld2"),v)
    file3p(v) = findfirst(==("w_200_normed_filt_N_lowpass_py2_gs_short_FiltInd_false_onlyfull_false_10000_TFTS.jld2"),v)
    file4p(v) = findfirst(==("w_200_normed_filt_N_lowpass_py2_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2"),v)


    if lowpass
        file_list = [file1, file2, file4, file7]  
    else
        file_list = [file1, file2, file4, file6]
    end

    if py2
        if lowpass
            file_list = [file1p, file2p, file4p, file7]
        else
            file_list = [file1p, file2p, file4p, file6]
        end
    end
    version_names = ["Boers, 2018", "Modified significance testing", "Modified EWS calculation", "Modified data preprocessing"]
    if label_w_steps == true
        version_names = ["Boers, 2018", "Modified significance testing (Step 1)", "Modified EWS calculation (Step 2a)", "Modified data preprocessing (Step 3)"]
    end

    fts = Figure(size=(1200, 200*length(file_list)))
    letters_h = reshape(letters[1:18],2,9)
    
    whichletter = 0



        
    for (i,ff) in enumerate(file_list)
        whichletter +=1

        good = ff(var_names)
        v = vars[good]
        a = acs[good]

        n_v = 0 
        n_a = 0
        n_both = 0
        which_v = []
        
        for k = 1:2
            gl = fts[whichletter,k] = GridLayout()
            
            

            if k == 1
                ax = Axis(gl[1,1],
                    xticks = (1:17, ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]) ,
                    xreversed = true,
                    xticklabelrotation = pi*0.4,
                    #yreversed=true,
                    ylabel =  L"\mathrm{𝖵}",
                    #xlabel = "Transition",
                    xminorgridwidth = 1,
                    xgridwidth=0,
                    xminorgridcolor = :grey10 ,
                    xgridcolor=:transparent,
                    xminorgridvisible = true, 
                    xgridvisible = false,
                    #yticks = (1:6, ["NGRIP 5 y", "NGRIP 10 y", "NEEM 10 y", "NGRIP 20 y", "GRIP 20 y", "GISP2 20 y"]),
                    #ylabel = "Ice core",
                    #yminorgridwidth = 1,
                    #ygridwidth=0.7,
                    #yminorgridcolor = :grey10 ,
                    #ygridcolor=:grey30,
                    yminorgridvisible = false, 
                    ygridvisible = false,
                    yticksvisible=false,
                    yticklabelsvisible=false,
                    ylabelsize = 18,
                    ylabelfont = :bold,
                    )

                if whichletter == length(file_list)
                    ax.xlabel = "Event"
                else
                    hidexdecorations!(ax, grid=false, minorgrid=false)
                    hidespines!(ax, :b) 
                end

                if i>1
                    hidespines!(ax,:t)
                    #hidexdecorations!(ax_event, grid=false)
                end

                ma = zeros(17,1)
                for ev in 1:17
                    #va = 0.4
                    #vl = 2.0
                    #vsc = (:white,1.0)
                    #vc = :black
                    if typeof(v.slopes[ev]) != Nothing
                        if v.slopes[ev] >0 
                            ma[ev,1]+=1
                    #        vc = :red
                    #        vsc = (:darkred,0.1)
                            if v.p_one[ev] < plim
                                ma[ev,1]+=1
                    #            va = 0.9
                    #            vl = 3.5
                    #            vsc = (:darkred,0.8)
                                n_v +=1
                                push!(which_v, ev)
                            end
                        
                        else
                    #        vc = :blue
                    #        vsc = (:steelblue,0.1)
                        end
                    end
                    # vspan!(axv, [GI_onsets[ev]], [GS_onsets[ev]], color = vsc)
                    # pv = Polynomials.fit(v.times[ev][findall(!ismissing, v.vals[ev])], convert.(Float64,v.vals[ev][findall(!ismissing, v.vals[ev])]),1)
                    # lines!(axv,v.times[ev], v.vals[ev], color = :black, alpha = va)#, linewidth = vl)
                    # lines!(axv, v.times[ev], pv.(v.times[ev]), color = vc, alpha = va, linewidth = vl)
                    # vlines!(axv, [GI_onsets[ev]], color = :darkred, alpha = va, linewidth = vl)
                end
                co2 = cgrad([:steelblue, :darkred, :darkred], 3, categorical = true, alpha=[0.1,0.1,0.8])
                hm = heatmap!(ax, ma, colormap=co2, colorrange=(0,2))
                # xlims!(axv,nothing,nothing)
                # axv.xreversed = true
                # xlims!(ax_event, axv.xaxis.attributes.limits[]...)
            elseif k == 2
                ax = Axis(gl[1,1],
                    xticks = (1:17, ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]) ,
                    xreversed = true,
                    xticklabelrotation = pi*0.4,
                    #yreversed=true,
                    ylabel =  L"\mathrm{\alpha_𝟣}",
                    #xlabel = "Transition",
                    xminorgridwidth = 1,
                    xgridwidth=0,
                    xminorgridcolor = :grey10 ,
                    xgridcolor=:transparent,
                    xminorgridvisible = true, 
                    xgridvisible = false,
                    #yticks = (1:6, ["NGRIP 5 y", "NGRIP 10 y", "NEEM 10 y", "NGRIP 20 y", "GRIP 20 y", "GISP2 20 y"]),
                    #ylabel = "Ice core",
                    yminorgridwidth = 1,
                    ygridwidth=0.7,
                    yminorgridcolor = :grey10 ,
                    ygridcolor=:grey30,
                    yminorgridvisible = false, 
                    ygridvisible = false,
                    yticksvisible=false,
                    yticklabelsvisible=false,
                    ylabelsize = 18,
                    ylabelfont = :bold,
                    )

                if whichletter == length(file_list)
                    ax.xlabel = "Event"
                else
                    hidexdecorations!(ax, grid=false, minorgrid=false)
                    hidespines!(ax, :b) 
                end
                #hidedecorations!(ax, label=true)

                if i>1
                    hidespines!(ax,:t)
                    #hidexdecorations!(ax_event, grid=false)
                end

                ma = zeros(17,1)
                
                for ev in 1:17
                    #aa = 0.4
                    #al = 2.0
                    #asc = (:white,1.0)
                    #ac = :black
                    if typeof(a.slopes[ev]) != Nothing
                        if a.slopes[ev] >0 
                            #ac = :red
                            #asc = (:darkred,0.1)
                            ma[ev,1]+=1
                            if a.p_one[ev] < plim
                                ma[ev,1]+=1
                                #aa = 0.9
                                #al = 3.5
                                #asc = (:darkred,0.8)
                                n_a +=1
                                if ev in which_v
                                    n_both +=1
                                end
                            end
                        
                        #else
                        #    ac = :blue
                        #    asc = (:steelblue,0.1)
                        end
                    end
                    # vspan!(axa, [GI_onsets[ev]], [GS_onsets[ev]], color = asc)
                    # pa = Polynomials.fit(a.times[ev][findall(!ismissing, a.vals[ev])], convert.(Float64,a.vals[ev][findall(!ismissing, a.vals[ev])]),1)
                    # lines!(axa,a.times[ev], a.vals[ev], color = :black, alpha = aa)#, linewidth = vl)
                    # lines!(axa, a.times[ev], pa.(a.times[ev]), color = ac, alpha = aa, linewidth = al)
                    # vlines!(axa, [GI_onsets[ev]], color = :darkred, alpha = aa, linewidth = al)

                end
                co2 = cgrad([:steelblue, :darkred, :darkred], 3, categorical = true, alpha=[0.1,0.1,0.8])
                hm = heatmap!(ax, ma, colormap=co2, colorrange=(0,2))
                #ax.xlabel = "Transition"
                #xlims!(axa,nothing,nothing)
                #axa.xreversed = true
                #xlims!(ax_event, axa.xaxis.attributes.limits[]...)
                #ax_event.xreversed = true
            end


            Label(gl[1,1,TopLeft()], letters_h[k,whichletter], fontsize = 20,
                        font = :bold, padding = (0,20,-10,15),
                        halign = :left, valign =:bottom)
        end
            Label(fts[whichletter,1:end,Top()], L"\textbf{\textrm{%$(sstring(version_names[i])): } \mathrm{𝗇_{𝖵} = %$(sstring(n_v)), } \mathrm{𝗇_{\alpha_𝟣} = %$(sstring(n_a)), } \mathrm{𝗇_{𝖻𝗈𝗍𝗁} = %$(sstring(n_both))}}", padding = (0.0,0.0,20.0,20.0), valign = :center, font=:bold,fontsize = 20)

    end
    
    if legend
        # elems = [
        #     [LineElement(color= (:black, 0.4), points = Point2f[(0, 0.5), (0.5, 0.5)]), LineElement(color= (:black, 0.9), points = Point2f[(0.5, 0.5), (1, 0.5)])],
        #     [PolyElement(color = (:darkred, 0.1), strokewidth = 0), LineElement(color= (:red,0.4), linewidth=2.0)],# LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
        #     [PolyElement(color = (:steelblue,0.1), strokewidth = 0), LineElement(color= (:blue, 0.4), linewidth=2.0)],#,LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
        #     [PolyElement(color = (:darkred, 0.8), strokewidth = 0), LineElement(color= (:red,0.9), linewidth=3.5)],#,LineElement(color = (:darkred,0.9), linewidth= 3.5,points = Point2f[(1, 0), (1, 1)])]]
        #     ]
        # labels = ["EWS indicator", "increasing", "decreasing", "significantly increasing (𝗉<0.05)"]
        elems = [
            [PolyElement(color = (:steelblue,0.1), strokewidth = 0)],#,LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:darkred, 0.1), strokewidth = 0)],# LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:darkred, 0.8), strokewidth = 0)],#,LineElement(color = (:darkred,0.9), linewidth= 3.5,points = Point2f[(1, 0), (1, 1)])]]
            ]
        labels = ["increasing", "decreasing", "significantly increasing (𝗉<0.05)"]
        Legend(fts[end+1,1:end], elems, labels, orientation = :horizontal, valign = :bottom, margin = (10,10,10,10)) 
    end

    colgap!(fts.layout, 20)
    rowgap!(fts.layout, 0)
    Label(fts[1,1:end,Top()], "EWS in 100-year high-pass filtered NGRIP record with 5-year resolution", 
        padding = (0.0,0.0,80.0,3.0), valign = :bottom, font = :bold, fontsize = 20)

    if saving
        save(saveto,fts)
    end
    if showing
        display(fts)
    end
end

plot_compare_methods_shorter_ngrip_csd2_label_aggr(method = "TFTS",lowpass = true, plim = plim, label_w_steps=true,py2=true,legend = true,
            showing = showing, saving = saving, saveto = "do_ews_across_greenland_ice_cores/figures/fig4_aggr.pdf")
plot_compare_methods_shorter_ngrip_csd2_label_aggr(method = "TFTS",lowpass = false, plim = plim, label_w_steps=true,py2=true,legend = true,
            showing = showing, saving = saving, saveto = "do_ews_across_greenland_ice_cores/figures/figS15_aggr.pdf")

# Fig S13
function plot_compare_methods_ngrip_wave2(;method = "TFTS", smooth_w = smoothw, 
        s1 = 10, s2=50, 
        include_filtering = false,
        lowpass = lowpass, 
        plim = plim, 
        py2=true,
        showing = true,
        saving = false,
        legend = false,
        saveto ="paper/compare_methods_wave2_ngrip5_$(s1)_$(s2)_$(method)_inclFilt_$(include_filtering)_lowpass_$(lowpass)_smoothw_$(smooth_w)_p_$(plim).pdf")

    #wavelet comparison with all the steps
    type = "NGRIP5"
    _, _, _, _, _, sca_names, scas, hurst_names, hursts, _ = load_data(type)

    file1(v) = findfirst(==("NIKLAS_s1_$(s1)_s2_$(s2)_N_lowpass_PAUL_FiltInd_true_onlyfull_false_detrend_GLSAR_DemeanSur_true_smoothw_true_10000_FOURIER.jld2"), v)   

    #if smoothw
    file2w(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_PAUL_normed_filt_gs_short_FiltInd_true_normalise_false_nocoi_false_onlyfull_false_smoothw_true_10000_$method.jld2"), v) 
    file3w(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_PAUL_normed_filt_gs_short_FiltInd_false_normalise_false_nocoi_false_onlyfull_false_smoothw_true_10000_$method.jld2"), v)
    file4w(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_PAUL_normed_filt_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_true_10000_$method.jld2"), v)
    file5w(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_no_lowpass_PAUL_normed_filt_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_true_10000_$method.jld2"), v)  
    file6w(v) = findfirst(==("s1_$(s1)_s2_$(s2)_C_no_lowpass_PAUL_normed_filt_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_true_10000_$method.jld2"), v)
    file7w(v) = findfirst(==("s1_$(s1)_s2_$(s2)_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_true_10000_$method.jld2"), v)
    file8w(v) = findfirst(==("s1_$(s1)_s2_$(s2)_C_lowpass_PAUL_normed_filt_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_true_10000_$method.jld2"), v)
    file9w(v) = findfirst(==("s1_$(s1)_s2_$(s2)_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_true_10000_$method.jld2"), v)

    #else
    file2(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_PAUL_normed_filt_gs_short_FiltInd_true_normalise_false_nocoi_false_onlyfull_false_smoothw_true_10000_$method.jld2"), v)
    file3(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_PAUL_normed_filt_gs_short_FiltInd_false_normalise_false_nocoi_false_onlyfull_false_smoothw_true_10000_$method.jld2"), v)
    file4(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_PAUL_normed_filt_gs_short_FiltInd_false_normalise_false_nocoi_false_onlyfull_false_smoothw_false_10000_$method.jld2"), v)
    file5(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_PAUL_normed_filt_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_false_smoothw_false_10000_$method.jld2"), v)
    file6(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_no_lowpass_PAUL_normed_filt_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_false_smoothw_false_10000_$method.jld2"), v)
    file7(v) = findfirst(==("s1_$(s1)_s2_$(s2)_C_no_lowpass_PAUL_normed_filt_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_false_smoothw_false_10000_$method.jld2"), v)
    file8(v) = findfirst(==("s1_$(s1)_s2_$(s2)_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_$method.jld2"), v)
    file9(v) = findfirst(==("s1_$(s1)_s2_$(s2)_C_lowpass_PAUL_normed_filt_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_false_smoothw_false_10000_$method.jld2"), v)
    file10(v) = findfirst(==("s1_$(s1)_s2_$(s2)_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_$method.jld2"), v)

    file1p(v) = findfirst(==("NIKLAS_s1_$(s1)_s2_$(s2)_N_lowpass_py2_PAUL_FiltInd_true_onlyfull_false_detrend_GLSAR_DemeanSur_true_smoothw_true_10000_FOURIER.jld2"),v)
    file2p(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_py2_PAUL_normed_gs_short_FiltInd_true_normalise_false_nocoi_false_onlyfull_false_smoothw_true_10000_TFTS.jld2"),v)
    file3p(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_py2_PAUL_normed_gs_short_FiltInd_false_normalise_false_nocoi_false_onlyfull_false_smoothw_true_10000_TFTS.jld2"),v)
    file4p(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_py2_PAUL_normed_gs_short_FiltInd_false_normalise_false_nocoi_false_onlyfull_false_smoothw_false_10000_TFTS.jld2"),v)
    file5p(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_py2_PAUL_normed_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_false_smoothw_false_10000_TFTS.jld2"),v)
    file6p(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_py2_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_false_smoothw_false_10000_TFTS.jld2"),v)
    file7p(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_false_smoothw_false_10000_TFTS.jld2"),v)
    file8p(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_false_smoothw_false_10000_TFTS.jld2"),v)

    if smooth_w
        if include_filtering
            if lowpass
                file_list = [file1, file2w, file3w, file4w, file8w, file9w]
                version_names = ["Boers, 2018","TFTS surrogates on data","No indicator filtering","EWS in entire GS, excluding COI" ,"Raw ages not rounded","CWT on unfiltered data"]
            else
                file_list = [file1, file2w, file3w, file4w, file5w, file6w, file7w]
                version_names = ["Boers, 2018","TFTS surrogates on data","No indicator filtering","EWS in entire GS, excluding COI", "no lowpass after interpolation","Raw ages not rounded","CWT on unfiltered data"]
            end
        else
            if lowpass
                file_list = [file1, file2w,  file4w, file8w, file9w]
                version_names = ["Boers, 2018","TFTS surrogates on data","EWS in entire GS, excluding COI", "Raw ages not rounded","CWT on unfiltered data"]
            else
                file_list = [file1, file2w,  file4w, file5w, file6w, file7w]
                version_names = ["Boers, 2018","TFTS surrogates on data","EWS in entire GS, excluding COI", "No lowpass after interpolation","Raw ages not rounded","CWT on unfiltered data"]
            end
        end
    else
        if include_filtering
            if lowpass
                file_list =[file1, file2, file3, file4, file5, file9, file10]
                version_names = ["Boers, 2018","TFTS surrogates on data","No indicator filtering","no smoothing","EWS in entire GS, excluding COI","Raw ages not rounded","CWT on unfiltered data"]
                #version_names = ["Boers, 2018","TFTS surrogates on data (Step 1)","No indicator filtering (Step 2.1b)","No smoothing ov indicators (Step 2.2b)","EWS in entire GS, excluding COI (Step 2.3b)","Raw ages not rounded","CWT on unfiltered data"]
            else
                file_list =[file1, file2, file3, file4, file5, file6, file7, file8]
                version_names = ["Boers, 2018","TFTS surrogates on data","No indicator filtering","no smoothing","EWS in entire GS, excluding COI","No lowpass after interpolation","Raw ages not rounded","CWT on unfiltered data"]
            end
        else
            if lowpass
                file_list = [file1, file2, file4, file5, file9, file10]
                version_names = ["Boers, 2018","TFTS surrogates on data","No smoothing","EWS in entire GS, excluding COI","Raw ages not rounded","CWT on unfiltered data"]
            else
                file_list = [file1, file2, file4, file5, file6, file7, file8]
                version_names = ["Boers, 2018","TFTS surrogates on data","No smoothing","EWS in entire GS, excluding COI","No lowpass after interpolation","Raw ages not rounded","CWT on unfiltered data"]
            end
        end
    end

    if py2
        if include_filtering
            if lowpass
                file_list = [file1p, file2p, file3p, file4p, file5p, 
                                file6p, 
                                file7p, 
                                file10]
                version_names = ["Boers, 2018","TFTS surrogates on data (Step 1)","No indicator filtering (Step 2.1b)","No indicator smoothing (Step 2.2b)","EWS in entire GS, excluding COI (Step 2.3b)", 
                            "CWT on original data (not normed)", 
                            "Preprocessing in Julia (Step 3.1)", 
                            "Raw ages not rounded (Step 3.2)"]
            else
                file_list = [file1p, file2p, file3p, file4p, file5p, 
                                file6p, 
                                file7p, 
                                file8p,
                                file8]
                version_names = ["Boers, 2018","TFTS surrogates on data (Step 1)","No indicator filtering (Step 2.1b)","No indicator smoothing (Step 2.2b)","EWS in entire GS, excluding COI (Step 2.3b)", 
                            "CWT on original data (not normed)", 
                            "Preprocessing in Julia (Step 3.1)", 
                            "No low-pass filtering after interpolation",
                            "Raw ages not rounded (Step 3.2)"]
            end
        else
            if lowpass
                file_list = [file1p, file2p, file4p, file5p, 
                                file6p, 
                                file7p, 
                                file10]
                version_names = ["Boers, 2018","TFTS surrogates on data (Step 1)","No indicator smoothing (Step 2.2b)","EWS in entire GS, excluding COI (Step 2.3b)", 
                            "CWT on original data (not normed)", 
                            "Preprocessing in Julia (Step 3.1)", 
                            "Raw ages not rounded (Step 3.2)"]
            else
                file_list = [file1p, file2p, file4p, file5p, 
                                file6p, 
                                file7p, 
                                file8p,
                                file8]
                version_names = ["Boers, 2018","TFTS surrogates on data (Step 1)","No indicator smoothing (Step 2.2b)","EWS in entire GS, excluding COI (Step 2.3b)", 
                            "CWT on original data (not normed)", 
                            "Preprocessing in Julia (Step 3.1)", 
                            "No low-pass filtering after interpolation",
                           "Raw ages not rounded (Step 3.2)"]
            end
        end
    end

    letters_h = reshape(letters[1:18],2,9)
    whichletter = 0

    if method == "FOURIER"
        file_list = file_list[1:end-1]
    end

    fts = Figure(size=(1600, 250*length(file_list)))

    for (i,ff) in enumerate(file_list)
        whichletter +=1
        good = ff(sca_names)
        v = scas[good]
        a = hursts[good]

        n_v = 0 
        n_a = 0
        n_both = 0
        which_v = []

        for k = 1:2
            gl = fts[whichletter,k] = GridLayout()
            if k == 1
                axv = Axis(gl[1,1], 
                    ylabel = L"\mathrm{\hat{𝗐}^𝟤}", 
                    ylabelsize = 18,
                    ylabelfont = :bold,
                    xticks = (10_000:5_000:60_000, string.(10:5:60)), 
                    )
                axv.xreversed = true
                if whichletter == length(file_list)
                    axv.xlabel = "Age (kyr b2k)"
                else
                    hidexdecorations!(axv, grid=false)
                    hidespines!(axv, :b) 
                end

                if i>1
                    hidespines!(axv,:t)
                end

                for ev in 1:17
                    va = 0.4
                    vl = 2.0
                    vsc = (:white,1.0)
                    vc = :black
                    if typeof(v.slopes[ev]) != Nothing
                        if v.slopes[ev] >0 
                            vc = :red
                            vsc = (:darkred,0.1)
                            if v.p_one[ev] < plim
                                va = 0.9
                                vl = 3.5
                                vsc = (:darkred,0.8)
                                n_v +=1
                                push!(which_v, ev)
                            end

                        else
                            vc = :blue
                            vsc = (:steelblue,0.1)
                        end
                    end
                    vspan!(axv, [GI_onsets[ev]], [GS_onsets[ev]], color = vsc)
                    pv = Polynomials.fit(v.times[ev][findall(!ismissing, v.vals[ev])], convert.(Float64,v.vals[ev][findall(!ismissing, v.vals[ev])]),1)
                    lines!(axv,v.times[ev], v.vals[ev], color = :black, alpha = va)
                    lines!(axv, v.times[ev], pv.(v.times[ev]), color = vc, alpha = va, linewidth = vl)
                    vlines!(axv, [GI_onsets[ev]], color = :darkred, alpha = va, linewidth = vl)
                end
            elseif k == 2
                axa = Axis(gl[1,1], 
                    ylabel =  L"\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}}",
                    ylabelsize = 18,
                    ylabelfont = :bold,
                    xticks = (10_000:5_000:60_000, string.(10:5:60)), 
                    )

                axa.xreversed = true
                if whichletter == length(file_list)
                    axa.xlabel = "Age (kyr b2k)"
                else
                    hidexdecorations!(axa, grid=false)
                    hidespines!(axa, :b) 
                end

                if i>1
                    hidespines!(axa,:t)
                end
                for ev in 1:17
                    aa = 0.4
                    al = 2.0
                    asc = (:white,1.0)
                    ac = :black
                    if typeof(a.slopes[ev]) != Nothing
                        if a.slopes[ev] >0 
                            ac = :red
                            asc = (:darkred,0.1)
                            if a.p_one[ev] < plim
                                aa = 0.9
                                al = 3.5
                                asc = (:darkred,0.8)
                                n_a +=1
                                if ev in which_v
                                    n_both +=1
                                end
                            end

                        else
                            ac = :blue
                            asc = (:steelblue,0.1)
                        end
                    end
                    vspan!(axa, [GI_onsets[ev]], [GS_onsets[ev]], color = asc)
                    pa = Polynomials.fit(a.times[ev][findall(!ismissing, a.vals[ev])], convert.(Float64,a.vals[ev][findall(!ismissing, a.vals[ev])]),1)
                    lines!(axa,a.times[ev], a.vals[ev], color = :black, alpha = aa)#, linewidth = vl)
                    lines!(axa, a.times[ev], pa.(a.times[ev]), color = ac, alpha = aa, linewidth = al)
                    vlines!(axa, [GI_onsets[ev]], color = :darkred, alpha = aa, linewidth = al)
                end
            end
            Label(gl[1,1,TopLeft()], letters_h[k,whichletter], fontsize = 20,
                font = :bold, padding = (0,20,-10,15),#(0,5,5,0),
                halign = :left,valign =:bottom)
        end

        Label(fts[i,1:end,Top()], L"\textbf{\textrm{%$(sstring(version_names[i])): } \mathrm{𝗇_{\hat{𝗐}^𝟤} = %$(sstring(n_v)), } \mathrm{𝗇_{\hat{𝖧}^\text{𝗅𝗈𝖼}} = %$(sstring(n_a)), } \mathrm{𝗇_{𝖻𝗈𝗍𝗁} = %$(sstring(n_both))}}" ,
                            padding = (0.0,0.0,20.0,20.0), valign = :bottom, font = :bold, fontsize = 20)
    end
    if legend
        elems = [
            #[LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(0.5, 0), (0.5, 1)])],
            [LineElement(color= (:black, 0.4), points = Point2f[(0, 0.5), (0.5, 0.5)]), LineElement(color= (:black, 0.9), points = Point2f[(0.5, 0.5), (1, 0.5)])],
            [PolyElement(color = (:darkred, 0.1), strokewidth = 0), LineElement(color= (:red,0.4), linewidth=2.0)],# LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:steelblue,0.1), strokewidth = 0), LineElement(color= (:blue, 0.4), linewidth=2.0)],#,LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:darkred, 0.8), strokewidth = 0), LineElement(color= (:red,0.9), linewidth=3.5)],#,LineElement(color = (:darkred,0.9), linewidth= 3.5,points = Point2f[(1, 0), (1, 1)])]]
            ]
        labels = [#"GI onset", 
                "EWS indicator", "increasing", "decreasing", "significantly increasing (𝗉<0.05)"]
        Legend(fts[end+1,1:end], elems, labels, orientation = :horizontal, valign = :bottom, margin = (10,10,10,10)) 
    end

    colgap!(fts.layout, 20)
    rowgap!(fts.layout, 0)
    Label(fts[1,1:end,Top()], "EWS in ($(s1)-$(s2)) year band of NGRIP record with 5-year resolution", 
    padding = (0.0,0.0,60.0,20.0), valign = :top, font = :bold, fontsize = 20)
    if saving
        save(saveto,fts)
    end
    if showing
        display(fts)
    end
end

plot_compare_methods_ngrip_wave2(method = "TFTS", smooth_w = false, s1 = 10, s2=50, include_filtering = true, lowpass = true, 
    plim = plim, py2=true,showing = showing,saving = saving,legend=true,saveto ="do_ews_across_greenland_ice_cores/figures/figS13.pdf")


# Fig 7 and S18
function plot_compare_methods_shorter_ngrip_wave2_label(;method = "TFTS", smooth_w = smoothw, s1 = 10, s2=50,
        lowpass = lowpass, 
        plim = plim, 
        py2=true,
        label_w_steps = true,
        showing = true,
        saving = false,
        legend = false,
        saveto ="paper/compare_methods_wave_shorter2_ngrip5_$(s1)_$(s2)_$(method)_lowpass_$(lowpass)_smoothw_$(smoothw)_p_$(plim).pdf")
    #wavelet comparison with all the steps
    type = "NGRIP5"
    _, _, _, _, _, sca_names, scas, hurst_names, hursts, _ = load_data(type)

    file1(v) = findfirst(==("NIKLAS_s1_$(s1)_s2_$(s2)_N_lowpass_PAUL_FiltInd_true_onlyfull_false_detrend_GLSAR_DemeanSur_true_smoothw_true_10000_FOURIER.jld2"), v)   

    #if smoothw
    file2w(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_PAUL_normed_filt_gs_short_FiltInd_true_normalise_false_nocoi_false_onlyfull_false_smoothw_true_10000_$method.jld2"), v) 
    file4w(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_PAUL_normed_filt_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_true_10000_$method.jld2"), v)
    file7w(v) = findfirst(==("s1_$(s1)_s2_$(s2)_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_true_10000_$method.jld2"), v)
    file9w(v) = findfirst(==("s1_$(s1)_s2_$(s2)_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_true_10000_$method.jld2"), v)

    #else
    file2(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_PAUL_normed_filt_gs_short_FiltInd_true_normalise_false_nocoi_false_onlyfull_false_smoothw_true_10000_$method.jld2"), v)
    file5(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_PAUL_normed_filt_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_false_smoothw_false_10000_$method.jld2"), v)
    file8(v) = findfirst(==("s1_$(s1)_s2_$(s2)_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_$method.jld2"), v)
    file10(v) = findfirst(==("s1_$(s1)_s2_$(s2)_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_$method.jld2"), v)

    #py files
    file1p(v) = findfirst(==("NIKLAS_s1_$(s1)_s2_$(s2)_N_lowpass_py2_PAUL_FiltInd_true_onlyfull_false_detrend_GLSAR_DemeanSur_true_smoothw_true_10000_FOURIER.jld2"),v)
    file2p(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_py2_PAUL_normed_gs_short_FiltInd_true_normalise_false_nocoi_false_onlyfull_false_smoothw_true_10000_TFTS.jld2"),v)
    file5p(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_py2_PAUL_normed_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_false_smoothw_false_10000_TFTS.jld2"),v)

    if smooth_w
        if lowpass
            file_list = [file1, file2w, file4w, file9w]
        else
            file_list = [file1, file2w, file4w, file7w]
        end

        else
        if lowpass
            file_list = [file1, file2, file5, file10]
        else
            file_list = [file1, file2, file5, file8]
        end 
    end

    if py2
        if lowpass
            file_list = [file1p, file2p, file5p, file10]
        else
            file_list = [file1p, file2p, file5p, file8]
        end
    end

    version_names = ["Boers, 2018", "Modified significance testing", "Modified EWS calculation", "Modified data preprocessing"]
    if label_w_steps == true
        version_names = ["Boers, 2018", "Modified significance testing (Step 1)", "Modified EWS calculation (Step 2b)", "Modified data preprocessing (Step 3)"]
    end

    letters_h = reshape(letters[1:18],2,9)
    whichletter = 0

    if method == "FOURIER"
        file_list = file_list[1:end-1]
    end

    fts = Figure(size=(1200, 300*length(file_list)))

    for (i,ff) in enumerate(file_list)
        whichletter +=1
        good = ff(sca_names)
        v = scas[good]
        a = hursts[good]
        n_v = 0 
        n_a = 0
        n_both = 0
        which_v = []

        for k = 1:2
            gl = fts[whichletter,k] = GridLayout()
            if k == 1
                axv = Axis(gl[1,1], 
                        ylabel = L"\mathrm{\hat{𝗐}^𝟤}", 
                        ylabelsize = 18,
                        ylabelfont = :bold,
                        xticks = (10_000:5_000:60_000, string.(10:5:60)), 
                        )
                axv.xreversed = true

                ax_event = Axis(gl[1,1],xaxisposition = :top, xticklabelrotation = pi*0.4)
                ax_event.xticks = (GI_onsets, event_labels)
                ax_event.xreversed = true
                hidespines!(ax_event)
                hideydecorations!(ax_event)

                if whichletter == length(file_list)
                    axv.xlabel = "Age (kyr b2k)"
                else
                    hidexdecorations!(axv, grid=false)
                    hidespines!(axv, :b) 
                end
                
                if i>1
                    hidespines!(axv,:t)
                    hidexdecorations!(ax_event, grid=false)
                end

                for ev in 1:17
                    va = 0.4
                    vl = 2.0
                    vsc = (:white,1.0)
                    vc = :black
                    if typeof(v.slopes[ev]) != Nothing
                        if v.slopes[ev] >0 
                            vc = :red
                            vsc = (:darkred,0.1)
                            if v.p_one[ev] < plim
                                va = 0.9
                                vl = 3.5
                                vsc = (:darkred,0.8)
                                n_v +=1
                                push!(which_v, ev)
                            end
                        else
                            vc = :blue
                            vsc = (:steelblue,0.1)
                        end
                    end
                    vspan!(axv, [GI_onsets[ev]], [GS_onsets[ev]], color = vsc)
                    pv = Polynomials.fit(v.times[ev][findall(!ismissing, v.vals[ev])], convert.(Float64,v.vals[ev][findall(!ismissing, v.vals[ev])]),1)
                    lines!(axv,v.times[ev], v.vals[ev], color = :black, alpha = va)#, linewidth = vl)
                    lines!(axv, v.times[ev], pv.(v.times[ev]), color = vc, alpha = va, linewidth = vl)
                    vlines!(axv, [GI_onsets[ev]], color = :darkred, alpha = va, linewidth = vl)
                end
                xlims!(axv,nothing,nothing)
                axv.xreversed = true
                xlims!(ax_event, axv.xaxis.attributes.limits[]...)
                ax_event.xreversed = true
            elseif k == 2
                axa = Axis(gl[1,1], 
                        ylabel =  L"\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}}",
                        ylabelsize = 18,
                        ylabelfont = :bold,
                        xticks = (10_000:5_000:60_000, string.(10:5:60)), 
                        )

                axa.xreversed = true

                ax_event = Axis(gl[1,1],xaxisposition = :top, xticklabelrotation = pi*0.4)
                ax_event.xticks = (GI_onsets, event_labels)
                ax_event.xreversed = true
                hidespines!(ax_event)
                hideydecorations!(ax_event)

                if whichletter == length(file_list)
                    axa.xlabel = "Age (kyr b2k)"
                else
                    hidexdecorations!(axa, grid=false)
                    hidespines!(axa, :b) 
                end

                if i>1
                    hidespines!(axa,:t)
                    hidexdecorations!(ax_event, grid=false)
                end

                for ev in 1:17
                    aa = 0.4
                    al = 2.0
                    asc = (:white,1.0)
                    ac = :black
                    if typeof(a.slopes[ev]) != Nothing
                        if a.slopes[ev] >0 
                            ac = :red
                            asc = (:darkred,0.1)
                            if a.p_one[ev] < plim
                                aa = 0.9
                                al = 3.5
                                asc = (:darkred,0.8)
                                n_a +=1
                                if ev in which_v
                                    n_both +=1
                                end
                            end

                        else
                            ac = :blue
                            asc = (:steelblue,0.1)
                        end
                    end
                    vspan!(axa, [GI_onsets[ev]], [GS_onsets[ev]], color = asc)
                    pa = Polynomials.fit(a.times[ev][findall(!ismissing, a.vals[ev])], convert.(Float64,a.vals[ev][findall(!ismissing, a.vals[ev])]),1)
                    lines!(axa,a.times[ev], a.vals[ev], color = :black, alpha = aa)#, linewidth = vl)
                    lines!(axa, a.times[ev], pa.(a.times[ev]), color = ac, alpha = aa, linewidth = al)
                    vlines!(axa, [GI_onsets[ev]], color = :darkred, alpha = aa, linewidth = al)
                end
                xlims!(axa,nothing,nothing)
                axa.xreversed = true
                xlims!(ax_event, axa.xaxis.attributes.limits[]...)
                ax_event.xreversed = true
            end
            Label(gl[1,1,TopLeft()], letters_h[k,whichletter], fontsize = 20,
                font = :bold, padding = (0,20,-10,15),#(0,5,5,0),
                halign = :left,valign =:bottom)
        end
        Label(fts[i,1:end,Top()], L"\textbf{\textrm{%$(sstring(version_names[i])): } \mathrm{𝗇_{\hat{𝗐}^𝟤} = %$(sstring(n_v)), } \mathrm{𝗇_{\hat{𝖧}^\text{𝗅𝗈𝖼}} = %$(sstring(n_a)), } \mathrm{𝗇_{𝖻𝗈𝗍𝗁} = %$(sstring(n_both))}}" ,
                    padding = (0.0,0.0,20.0,20.0), valign = :center, font = :bold, fontsize = 20)
    end
    if legend
        elems = [
            #[LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(0.5, 0), (0.5, 1)])],
            [LineElement(color= (:black, 0.4), points = Point2f[(0, 0.5), (0.5, 0.5)]), LineElement(color= (:black, 0.9), points = Point2f[(0.5, 0.5), (1, 0.5)])],
            [PolyElement(color = (:darkred, 0.1), strokewidth = 0), LineElement(color= (:red,0.4), linewidth=2.0)],# LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:steelblue,0.1), strokewidth = 0), LineElement(color= (:blue, 0.4), linewidth=2.0)],#,LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:darkred, 0.8), strokewidth = 0), LineElement(color= (:red,0.9), linewidth=3.5)],#,LineElement(color = (:darkred,0.9), linewidth= 3.5,points = Point2f[(1, 0), (1, 1)])]]
            ]
        labels = [#"GI onset", 
                "EWS indicator", "increasing", "decreasing", "significantly increasing (𝗉<0.05)"]
        Legend(fts[end+1,1:end], elems, labels, orientation = :horizontal, valign = :bottom, margin = (10,10,10,10)) 
    end
    colgap!(fts.layout, 20)
    rowgap!(fts.layout, 0)
    Label(fts[1,1:end,Top()], "EWS in ($(s1)-$(s2)) year band of NGRIP record with 5-year resolution", 
    padding = (0.0,0.0,110.0,30.0), valign = :bottom, font = :bold, fontsize = 20)
    if saving
        save(saveto,fts)
    end
    if showing
        display(fts)
    end
end

plot_compare_methods_shorter_ngrip_wave2_label(method = "TFTS", smooth_w =false, s1 = 10, s2=50,lowpass = true, label_w_steps = true,
            plim = plim, py2=true,showing = showing,saving = saving,legend = true, saveto ="do_ews_across_greenland_ice_cores/figures/fig7.pdf")
plot_compare_methods_shorter_ngrip_wave2_label(method = "TFTS", smooth_w =false, s1 = 10, s2=50,lowpass = false, label_w_steps = true,
            plim = plim, py2=true,showing = showing,saving = saving,legend = true, saveto ="do_ews_across_greenland_ice_cores/figures/figS18.pdf")


function plot_compare_methods_shorter_ngrip_wave2_label_aggr(;method = "TFTS", smooth_w = smoothw, s1 = 10, s2=50,
    lowpass = lowpass, 
    plim = plim, 
    py2=true,
    label_w_steps = true,
    showing = true,
    saving = false,
    legend = false,
    saveto ="paper/compare_methods_wave_shorter2_aggr_ngrip5_$(s1)_$(s2)_$(method)_lowpass_$(lowpass)_smoothw_$(smoothw)_p_$(plim).pdf")
    
    type = "NGRIP5"
    _, _, _, _, _, sca_names, scas, hurst_names, hursts, _ = load_data(type)

    file1(v) = findfirst(==("NIKLAS_s1_$(s1)_s2_$(s2)_N_lowpass_PAUL_FiltInd_true_onlyfull_false_detrend_GLSAR_DemeanSur_true_smoothw_true_10000_FOURIER.jld2"), v)   

    #if smoothw
    file2w(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_PAUL_normed_filt_gs_short_FiltInd_true_normalise_false_nocoi_false_onlyfull_false_smoothw_true_10000_$method.jld2"), v) 
    file4w(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_PAUL_normed_filt_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_true_10000_$method.jld2"), v)
    file7w(v) = findfirst(==("s1_$(s1)_s2_$(s2)_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_true_10000_$method.jld2"), v)
    file9w(v) = findfirst(==("s1_$(s1)_s2_$(s2)_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_true_10000_$method.jld2"), v)

    #else
    file2(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_PAUL_normed_filt_gs_short_FiltInd_true_normalise_false_nocoi_false_onlyfull_false_smoothw_true_10000_$method.jld2"), v)
    file5(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_PAUL_normed_filt_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_false_smoothw_false_10000_$method.jld2"), v)
    file8(v) = findfirst(==("s1_$(s1)_s2_$(s2)_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_$method.jld2"), v)
    file10(v) = findfirst(==("s1_$(s1)_s2_$(s2)_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_$method.jld2"), v)

    #py files
    file1p(v) = findfirst(==("NIKLAS_s1_$(s1)_s2_$(s2)_N_lowpass_py2_PAUL_FiltInd_true_onlyfull_false_detrend_GLSAR_DemeanSur_true_smoothw_true_10000_FOURIER.jld2"),v)
    file2p(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_py2_PAUL_normed_gs_short_FiltInd_true_normalise_false_nocoi_false_onlyfull_false_smoothw_true_10000_TFTS.jld2"),v)
    file5p(v) = findfirst(==("s1_$(s1)_s2_$(s2)_N_lowpass_py2_PAUL_normed_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_false_smoothw_false_10000_TFTS.jld2"),v)

    if smooth_w
        if lowpass
            file_list = [file1, file2w, file4w, file9w]
        else
            file_list = [file1, file2w, file4w, file7w]
        end

        else
        if lowpass
            file_list = [file1, file2, file5, file10]
        else
            file_list = [file1, file2, file5, file8]
        end 
    end

    if py2
        if lowpass
            file_list = [file1p, file2p, file5p, file10]
        else
            file_list = [file1p, file2p, file5p, file8]
        end
    end

    version_names = ["Boers, 2018", "Modified significance testing", "Modified EWS calculation", "Modified data preprocessing"]
    if label_w_steps == true
        version_names = ["Boers, 2018", "Modified significance testing (Step 1)", "Modified EWS calculation (Step 2b)", "Modified data preprocessing (Step 3)"]
    end


    fts = Figure(size=(1200, 200*length(file_list)))
    letters_h = reshape(letters[1:18],2,9)
    
    whichletter = 0

    if method == "FOURIER"
        file_list = file_list[1:end-1]
    end


        
    for (i,ff) in enumerate(file_list)
        whichletter +=1

        good = ff(sca_names)
        v = scas[good]
        a = hursts[good]

        n_v = 0 
        n_a = 0
        n_both = 0
        which_v = []
        
        for k = 1:2
            gl = fts[whichletter,k] = GridLayout()
            

            if k == 1
                ax = Axis(gl[1,1],
                    xticks = (1:17, ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]) ,
                    xreversed = true,
                    xticklabelrotation = pi*0.4,
                    #yreversed=true,
                    ylabel = L"\mathrm{\hat{𝗐}^𝟤}", 
                    #xlabel = "Transition",
                    xminorgridwidth = 1,
                    xgridwidth=0,
                    xminorgridcolor = :grey10 ,
                    xgridcolor=:transparent,
                    xminorgridvisible = true, 
                    xgridvisible = false,
                    #yticks = (1:6, ["NGRIP 5 y", "NGRIP 10 y", "NEEM 10 y", "NGRIP 20 y", "GRIP 20 y", "GISP2 20 y"]),
                    #ylabel = "Ice core",
                    #yminorgridwidth = 1,
                    #ygridwidth=0.7,
                    #yminorgridcolor = :grey10 ,
                    #ygridcolor=:grey30,
                    yminorgridvisible = false, 
                    ygridvisible = false,
                    yticksvisible=false,
                    yticklabelsvisible=false,
                    ylabelsize = 18,
                    ylabelfont = :bold,
                    )

                if whichletter == length(file_list)
                    ax.xlabel = "Event"
                else
                    hidexdecorations!(ax, grid=false, minorgrid=false)
                    hidespines!(ax, :b) 
                end

                if i>1
                    hidespines!(ax,:t)
                    #hidexdecorations!(ax_event, grid=false)
                end

                ma = zeros(17,1)
                for ev in 1:17
                    #va = 0.4
                    #vl = 2.0
                    #vsc = (:white,1.0)
                    #vc = :black
                    if typeof(v.slopes[ev]) != Nothing
                        if v.slopes[ev] >0 
                            ma[ev,1]+=1
                    #        vc = :red
                    #        vsc = (:darkred,0.1)
                            if v.p_one[ev] < plim
                                ma[ev,1]+=1
                    #            va = 0.9
                    #            vl = 3.5
                    #            vsc = (:darkred,0.8)
                                n_v +=1
                                push!(which_v, ev)
                            end
                        
                        else
                    #        vc = :blue
                    #        vsc = (:steelblue,0.1)
                        end
                    end
                    # vspan!(axv, [GI_onsets[ev]], [GS_onsets[ev]], color = vsc)
                    # pv = Polynomials.fit(v.times[ev][findall(!ismissing, v.vals[ev])], convert.(Float64,v.vals[ev][findall(!ismissing, v.vals[ev])]),1)
                    # lines!(axv,v.times[ev], v.vals[ev], color = :black, alpha = va)#, linewidth = vl)
                    # lines!(axv, v.times[ev], pv.(v.times[ev]), color = vc, alpha = va, linewidth = vl)
                    # vlines!(axv, [GI_onsets[ev]], color = :darkred, alpha = va, linewidth = vl)
                end
                co2 = cgrad([:steelblue, :darkred, :darkred], 3, categorical = true, alpha=[0.1,0.1,0.8])
                hm = heatmap!(ax, ma, colormap=co2, colorrange=(0,2))
                # xlims!(axv,nothing,nothing)
                # axv.xreversed = true
                # xlims!(ax_event, axv.xaxis.attributes.limits[]...)
            elseif k == 2
                ax = Axis(gl[1,1],
                    xticks = (1:17, ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]) ,
                    xreversed = true,
                    xticklabelrotation = pi*0.4,
                    #yreversed=true,
                    ylabel =  L"\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}}",
                    #xlabel = "Transition",
                    xminorgridwidth = 1,
                    xgridwidth=0,
                    xminorgridcolor = :grey10 ,
                    xgridcolor=:transparent,
                    xminorgridvisible = true, 
                    xgridvisible = false,
                    #yticks = (1:6, ["NGRIP 5 y", "NGRIP 10 y", "NEEM 10 y", "NGRIP 20 y", "GRIP 20 y", "GISP2 20 y"]),
                    #ylabel = "Ice core",
                    yminorgridwidth = 1,
                    ygridwidth=0.7,
                    yminorgridcolor = :grey10 ,
                    ygridcolor=:grey30,
                    yminorgridvisible = false, 
                    ygridvisible = false,
                    yticksvisible=false,
                    yticklabelsvisible=false,
                    ylabelsize = 18,
                    ylabelfont = :bold,
                    )

                if whichletter == length(file_list)
                    ax.xlabel = "Event"
                else
                    hidexdecorations!(ax, grid=false, minorgrid=false)
                    hidespines!(ax, :b) 
                end
                #hidedecorations!(ax, label=true)

                if i>1
                    hidespines!(ax,:t)
                    #hidexdecorations!(ax_event, grid=false)
                end

                ma = zeros(17,1)
                
                for ev in 1:17
                    #aa = 0.4
                    #al = 2.0
                    #asc = (:white,1.0)
                    #ac = :black
                    if typeof(a.slopes[ev]) != Nothing
                        if a.slopes[ev] >0 
                            #ac = :red
                            #asc = (:darkred,0.1)
                            ma[ev,1]+=1
                            if a.p_one[ev] < plim
                                ma[ev,1]+=1
                                #aa = 0.9
                                #al = 3.5
                                #asc = (:darkred,0.8)
                                n_a +=1
                                if ev in which_v
                                    n_both +=1
                                end
                            end
                        
                        #else
                        #    ac = :blue
                        #    asc = (:steelblue,0.1)
                        end
                    end
                    # vspan!(axa, [GI_onsets[ev]], [GS_onsets[ev]], color = asc)
                    # pa = Polynomials.fit(a.times[ev][findall(!ismissing, a.vals[ev])], convert.(Float64,a.vals[ev][findall(!ismissing, a.vals[ev])]),1)
                    # lines!(axa,a.times[ev], a.vals[ev], color = :black, alpha = aa)#, linewidth = vl)
                    # lines!(axa, a.times[ev], pa.(a.times[ev]), color = ac, alpha = aa, linewidth = al)
                    # vlines!(axa, [GI_onsets[ev]], color = :darkred, alpha = aa, linewidth = al)

                end
                co2 = cgrad([:steelblue, :darkred, :darkred], 3, categorical = true, alpha=[0.1,0.1,0.8])
                hm = heatmap!(ax, ma, colormap=co2, colorrange=(0,2))
                #ax.xlabel = "Transition"
                #xlims!(axa,nothing,nothing)
                #axa.xreversed = true
                #xlims!(ax_event, axa.xaxis.attributes.limits[]...)
                #ax_event.xreversed = true
            end


            Label(gl[1,1,TopLeft()], letters_h[k,whichletter], fontsize = 20,
                        font = :bold, padding = (0,20,-10,15),
                        halign = :left, valign =:bottom)
        end
            Label(fts[whichletter,1:end,Top()], L"\textbf{\textrm{%$(sstring(version_names[i])): } \mathrm{𝗇_{\hat{𝗐}^𝟤} = %$(sstring(n_v)), } \mathrm{𝗇_{\hat{𝖧}^\text{𝗅𝗈𝖼}} = %$(sstring(n_a)), } \mathrm{𝗇_{𝖻𝗈𝗍𝗁} = %$(sstring(n_both))}}", 
            padding = (0.0,0.0,20.0,20.0), valign = :center, font=:bold,fontsize = 20)

    end
    
    if legend
        # elems = [
        #     [LineElement(color= (:black, 0.4), points = Point2f[(0, 0.5), (0.5, 0.5)]), LineElement(color= (:black, 0.9), points = Point2f[(0.5, 0.5), (1, 0.5)])],
        #     [PolyElement(color = (:darkred, 0.1), strokewidth = 0), LineElement(color= (:red,0.4), linewidth=2.0)],# LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
        #     [PolyElement(color = (:steelblue,0.1), strokewidth = 0), LineElement(color= (:blue, 0.4), linewidth=2.0)],#,LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
        #     [PolyElement(color = (:darkred, 0.8), strokewidth = 0), LineElement(color= (:red,0.9), linewidth=3.5)],#,LineElement(color = (:darkred,0.9), linewidth= 3.5,points = Point2f[(1, 0), (1, 1)])]]
        #     ]
        # labels = ["EWS indicator", "increasing", "decreasing", "significantly increasing (𝗉<0.05)"]
        elems = [
            [PolyElement(color = (:steelblue,0.1), strokewidth = 0)],#,LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:darkred, 0.1), strokewidth = 0)],# LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:darkred, 0.8), strokewidth = 0)],#,LineElement(color = (:darkred,0.9), linewidth= 3.5,points = Point2f[(1, 0), (1, 1)])]]
            ]
        labels = ["increasing", "decreasing", "significantly increasing (𝗉<0.05)"]
        Legend(fts[end+1,1:end], elems, labels, orientation = :horizontal, valign = :bottom, margin = (10,10,10,10)) 
    end

    colgap!(fts.layout, 20)
    rowgap!(fts.layout, 0)
    Label(fts[1,1:end,Top()], "EWS in($(s1)-$(s2)) year band of NGRIP record with 5-year resolution", 
        padding = (0.0,0.0,80.0,3.0), valign = :bottom, font = :bold, fontsize = 20)

    if saving
        save(saveto,fts)
    end
    if showing
        display(fts)
    end
end

plot_compare_methods_shorter_ngrip_wave2_label_aggr(method = "TFTS", smooth_w =false, s1 = 10, s2=50,lowpass = true, label_w_steps = true,
            plim = plim, py2=true,showing = showing,saving = true,legend = true, saveto ="do_ews_across_greenland_ice_cores/figures/fig7_aggr.pdf")
plot_compare_methods_shorter_ngrip_wave2_label_aggr(method = "TFTS", smooth_w =false, s1 = 10, s2=50,lowpass = false, label_w_steps = true,
            plim = plim, py2=true,showing = showing,saving = true,legend = true, saveto ="do_ews_across_greenland_ice_cores/figures/figS18_aggr.pdf")

# Fig3 and S16
function plot_overall_significance_example_csd(;lowpass = true, plim = plim, showing = true, saving = false, saveto = "paper/ngrip5_example_overall_signif_csd_with_both_lowpass_$(lowpass)_p_$(plim).pdf")
    # plot distribution of number of "expected" false ews
    function good_files(v,lowpass)
        if lowpass
            return occursin("w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_n_2000_nGS_1000_p_0.05_0.1_0.2_0.3_TFTS.jld2", v)
        else
            return occursin("w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_n_2000_nGS_1000_p_0.05_0.1_0.2_0.3_TFTS.jld2", v)
        end
    end

    function get_ns(x,y,plim)
        n_x=0
        n_y=0
        n_both = 0
        which_v = []
        for ev in 1:17
            if typeof(x.slopes[ev]) != Nothing
                if x.slopes[ev] >0  && x.p_one[ev] < plim
                    n_x+=1
                    push!(which_v, ev)
                end
            end
            if typeof(y.slopes[ev]) != Nothing
                if y.slopes[ev] >0  && y.p_one[ev] < plim
                    n_y+=1
                    if ev in which_v
                        n_both+=1
                    end
                end
            end

        end
        return n_x, n_y, n_both
    end

    if lowpass
        v = load("new_surrogate_files/NGRIP5/var/w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2")["slopes"]
        a= load("new_surrogate_files/NGRIP5/ac/w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2")["slopes"]    
    else   
        v = load("new_surrogate_files/NGRIP5/var/w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2")["slopes"]
        a= load("new_surrogate_files/NGRIP5/ac/w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2")["slopes"]
    end       
    
    n_x, n_y, n_both = get_ns(v,a,0.05)
    
    whichletter = 0
    
    fhist = Figure(size = (1200,1000))
    cm = cgrad(:managua,7, categorical=true, rev=true)[[2,6,4]]
    for (i,ind) in enumerate(["var","ac"])
        
        dist_path = readdir("new_surrogate_files/NGRIP5/n_sig_dists/$ind/", join=true)
        dist_name = readdir("new_surrogate_files/NGRIP5/n_sig_dists/$ind", join=false)
        dists = [load(k)["distribution"] for k in dist_path]
        
        good = @. good_files(dist_name,lowpass)
        ga = fhist[1,i] = GridLayout()
        for (k,test_dist) in enumerate(dists[good])
            for (j,p) in enumerate(test_dist.pvals)
                if p == plim
                    jj = 1
                    whichletter +=1
                    ax = Axis(ga[jj,k], 
                            titlesize=25,
                            titlefont = :bold,
                            ylabel = "Relative frequency", 
                            yticks = 0:0.2:1, 
                            xlabel = "No. of false significant EWS (out of 17)",
                            xticks = 0:1:17,
                            limits= (-0,6.6,nothing,nothing)
                            )
                    if ind == "var"
                        ax.title = L"\textbf{\mathrm{𝖵}}"
                        vlines!(ax, [n_x], color=cm[2], linewidth=4, label = "Number of significant EWS",alpha=0.7)
                    elseif ind == "ac"
                        ax.title = L"\textbf{\mathrm{\alpha_𝟣}}"
                        vlines!(ax, [n_y], color=cm[2], linewidth=4, label = "Number of significant EWS",alpha=0.7)
                    end
                    ma = Int(maximum(test_dist.num_inc_one[jj,:]))
                    h= hist!(ax,test_dist.num_inc_one[jj,:], bins = -0.5:ma+0.5,#ma+1, 
                                    normalization = :pdf,
                                    color = (cm[1],0.3),
                                    label = "Histogram numeric"
                                    )
                    
                    density!(ax,test_dist.num_inc_one[jj,:], npoints = 17+2,#17+2,
                            color = (:steelblue,0.0), 
                            strokecolor = cm[1],#:steelblue, 
                            strokewidth = 2, #strokearound = true,
                            boundary = (-1,17),
                            )
                    lines!(ax, 1:5,1:5, color = cm[1], linewidth=2, label = "Distribution numeric", visible = false)


                    # vlines!(ax,[quantile(test_dist.num_inc_one[jj,:],1-p)], 
                    #     color = cm[1],
                    #     linewidth = 2, linestyle=(:dot,:dense),
                    #     label = "0.95 significance level numeric")

                    thr_num = significance_threshold(test_dist,round(1-p, digits=2))
                    vlines!(ax,[thr_num-0.5], 
                        color = cm[1],
                        linewidth = 2, linestyle=(:dot,:dense),
                        label = "0.95 significance level numeric")

                    n=17
                    b = Distributions.Binomial(n,p)

                    lines!(ax,0:17, pdf.(b,0:17), color = cm[3],
                        linewidth = 2, 
                        label ="Distribution analytic")
                    # vlines!(ax,[quantile(b,1-p)], color = cm[3],
                    #     linestyle = :dash,
                    #     linewidth = 2, 
                    #     label = "$(1-plim) significance level analytic")
                    thr_ana = significance_threshold(b,round(1-p, digits=2))
                    vlines!(ax,[thr_ana-0.5], color = cm[3],
                        linestyle = :dash,
                        linewidth = 2, 
                        label = "$(1-plim) significance level analytic")

                    Label(ga[1,1,TopLeft()], letters[whichletter], fontsize = 20,
                        font = :bold, padding = (0,5,5,0),
                        halign = :left, valign=:bottom)

                    if k==1
                        Legend(fhist[2,2], ax,
                            framevisible=false,
                            tellwidth = false,
                            halign = :center,
                            valign = :center,
                            rowgap = 20,
                            framecolor = :grey50)
                    end
                end
            end
        end
    end
    
    gboth = fhist[2,1] = GridLayout()
    whichletter +=1
    ax3 = Axis(gboth[1,1], title = L"\textbf{\mathrm{𝖵} 𝖺𝗇𝖽 \mathrm{\alpha_𝟣}}",
                            titlesize=25,
                            titlefont = :bold,
                            ylabel = "Relative frequency",
                            yticks = 0:0.2:1, 
                            xlabel = "No. of false significant EWS (out of 17)",
                            xticks = 0:1:17,
                            limits= (-0.1,6.6,nothing,nothing)
                            )
    pboth = round(plim*plim, digits=5)
    bboth = Distributions.Binomial(17,pboth)
    
    lines!(ax3,
            0:17, 
            pdf.(bboth,0:17), 
            color=cm[3],
            linewidth = 2, 
            label ="Distribution analytic (B(17,$(round(pboth, digits=4))))")
    vlines!(ax3, [n_both], color=cm[2], linewidth=4, alpha = 0.7)
    # vlines!(ax3,[quantile(bboth,1-plim)], color = cm[3],
    #             linestyle = :dash,
    #             linewidth = 2, 
    #             label = "$(1-plim) significance level analytic")
    thr_ana2 = significance_threshold(bboth,round(1-plim, digits=2))
    vlines!(ax3,[thr_ana2-0.5], color = cm[3],
                linestyle = :dash,
                linewidth = 2, 
                label = "$(1-plim) significance level analytic")

    
    Label(gboth[1,1,TopLeft()], letters[whichletter], fontsize = 20,
                    font = :bold, padding = (0,5,5,0),
                    halign = :left, valign=:bottom)

    linkyaxes!(fhist.content[1:3:end-1]...)
    colgap!(fhist.layout, 40)
    rowgap!(fhist.layout, 40)
    if saving
        save(saveto,fhist)
    end
    if showing
        display(fhist)
    end
end


plot_overall_significance_example_csd(lowpass=true, showing=showing, saving = saving, saveto = "do_ews_across_greenland_ice_cores/figures/new_fig3.pdf")
plot_overall_significance_example_csd(lowpass=false,showing=showing,  saving = saving, saveto = "do_ews_across_greenland_ice_cores/figures/new_figS16.pdf")


#create struct object for hurst (no lowpass, no smoothw) of the combined overall significance files
begin
    disth_path = readdir("new_surrogate_files/NGRIP5/n_sig_dists/hurst/", join=true)
    disth_name = readdir("new_surrogate_files/NGRIP5/n_sig_dists/hurst", join=false)
    
    #### Hurst no lowpass
    good_dists = [load(k)["distribution"] for k in disth_path[findall(occursin.("nGS_1000", disth_name) .* occursin.("s1_10_s2_50_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false", disth_name))] ]
    new_ninc_one = zeros(4,2000)
    new_ninc_two = zeros(4,2000)
    curr_ind = 1
    for (i,d) in enumerate(good_dists)
        global curr_ind
        gotsofar = findlast(>(0), sum(d.num_inc_one, dims = 1))[2]
        new_ninc_one[:,curr_ind:gotsofar+curr_ind-1] = d.num_inc_one[:,1:gotsofar]
        new_ninc_two[:,curr_ind:gotsofar+curr_ind-1] = d.num_inc_two[:,1:gotsofar]
        curr_ind += gotsofar
    end
    gotsofar_overall = findlast(>(0), sum(new_ninc_one, dims = 1))[2]
    @show "Hurst no lowpass:", gotsofar_overall
    
    #make object:
    combined_hurst_no_lp = distribution_significant_increases(
        "C_no_lowpass", 
        new_ninc_one[:,1:gotsofar_overall],
        new_ninc_two[:,1:gotsofar_overall], 
        gotsofar_overall,
        1000,
        [0.05,0.1,0.2,0.3],
        "hurst",
        10,
        50,
        5)

    #### Hurst lowpass
    good_dists = [load(k)["distribution"] for k in disth_path[findall(occursin.("nGS_1000", disth_name) .* occursin.("s1_10_s2_50_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false", disth_name))] ]
    new_ninc_one = zeros(4,2000)
    new_ninc_two = zeros(4,2000)
    curr_ind = 1
    for (i,d) in enumerate(good_dists)
        global curr_ind
        gotsofar = findlast(>(0), sum(d.num_inc_one, dims = 1))[2]
        new_ninc_one[:,curr_ind:gotsofar+curr_ind-1] = d.num_inc_one[:,1:gotsofar]
        new_ninc_two[:,curr_ind:gotsofar+curr_ind-1] = d.num_inc_two[:,1:gotsofar]
        curr_ind += gotsofar
    end
    gotsofar_overall = findlast(>(0), sum(new_ninc_one, dims = 1))[2]
    @show "Hurst lowpass:", gotsofar_overall
    
    #make object:
    combined_hurst_lp = distribution_significant_increases(
        "C_lowpass", 
        new_ninc_one[:,1:gotsofar_overall],
        new_ninc_two[:,1:gotsofar_overall], 
        gotsofar_overall,
        1000,
        [0.05,0.1,0.2,0.3],
        "hurst",
        10,
        50,
        5)


    #### Sca lowpass
    dists_path = readdir("new_surrogate_files/NGRIP5/n_sig_dists/sca/", join=true)
    dists_name = readdir("new_surrogate_files/NGRIP5/n_sig_dists/sca", join=false)

    good_dists = [load(k)["distribution"] for k in dists_path[findall(occursin.("nGS_1000", dists_name) .* occursin.("s1_10_s2_50_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false", dists_name))] ]
    new_ninc_one = zeros(4,2000)
    new_ninc_two = zeros(4,2000)
    curr_ind = 1
    for (i,d) in enumerate(good_dists)
        global curr_ind
        gotsofar = findlast(>(0), sum(d.num_inc_one, dims = 1))[2]
        new_ninc_one[:,curr_ind:gotsofar+curr_ind-1] = d.num_inc_one[:,1:gotsofar]
        new_ninc_two[:,curr_ind:gotsofar+curr_ind-1] = d.num_inc_two[:,1:gotsofar]
        curr_ind += gotsofar
    end
    gotsofar_overall = findlast(>(0), sum(new_ninc_one, dims = 1))[2]
    @show "Sca lowpass:", gotsofar_overall
    
    #make object:
    combined_sca_lp = distribution_significant_increases(
        "C_lowpass", 
        new_ninc_one[:,1:gotsofar_overall],
        new_ninc_two[:,1:gotsofar_overall], 
        gotsofar_overall,
        1000,
        [0.05,0.1,0.2,0.3],
        "sca",
        10,
        50,
        5)
end;


#Fig A3 and S19
function plot_overall_significance_example_wave(;lowpass = true, plim = plim, showing = true, saving = false, saveto = "paper/ngrip5_example_overall_signif_wave_with_both_lowpass_$(lowpass)_p_$(plim).pdf")
    # plot distribution of number of "expected" false ews

    function get_ns(x,y,plim)
        n_x=0
        n_y=0
        n_both = 0
        which_v = []
        for ev in 1:17
            if typeof(x.slopes[ev]) != Nothing
                if x.slopes[ev] >0  && x.p_one[ev] < plim
                    n_x+=1
                    push!(which_v, ev)
                end
            end
            if typeof(y.slopes[ev]) != Nothing
                if y.slopes[ev] >0  && y.p_one[ev] < plim
                    n_y+=1
                    if ev in which_v
                        n_both+=1
                    end
                end
            end

        end
        return n_x, n_y, n_both
    end

    if lowpass
        v = load("new_surrogate_files/NGRIP5/sca/s1_10_s2_50_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2")["slopes"]
        a= load("new_surrogate_files/NGRIP5/hurst/s1_10_s2_50_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2")["slopes"]  
    else
        v = load("new_surrogate_files/NGRIP5/sca/s1_10_s2_50_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2")["slopes"]
        a= load("new_surrogate_files/NGRIP5/hurst/s1_10_s2_50_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2")["slopes"]  
    end
    
    n_x, n_y, n_both = get_ns(v,a,0.05)
    
    whichletter = 0
    
    fhist = Figure(size = (1200,1000))
    cm = cgrad(:managua,7, categorical=true, rev=true)[[2,6,4]]
    for (i,ind) in enumerate(["sca","hurst"])
        
        dist_path = readdir("new_surrogate_files/NGRIP5/n_sig_dists/$ind/", join=true)
        dist_name = readdir("new_surrogate_files/NGRIP5/n_sig_dists/$ind", join=false)
        if ind == "sca"
            if lowpass
                global combined_sca_lp 
                dists = [combined_sca_lp]
            else
                dists = [load(k)["distribution"] for k in dist_path[findall(occursin.("s1_10_s2_50_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_n_2000_nGS_1000_p_0.05_0.1_0.2_0.3_TFTS.jld2", dist_name))]]
            end
        elseif ind == "hurst"
            if lowpass
                global combined_hurst_lp
                dists = [combined_hurst_lp]
            else
                global combined_hurst_no_lp
                dists = [combined_hurst_no_lp]
            end
        end
        ga = fhist[1,i] = GridLayout()
        for (k,test_dist) in enumerate(dists)
            gotsofar = findlast(>(0), sum(test_dist.num_inc_one, dims = 1))[2]
            if gotsofar < 2000
                println("WARNING: only got $(gotsofar)/2000 surrogates for $(ind)")
            end
            for (j,p) in enumerate(test_dist.pvals)
                if p == plim
                    jj = 1
                    whichletter +=1
                    ax = Axis(ga[jj,k], 
                            titlesize=25,
                            titlefont = :bold,
                            ylabel = "Relative frequency", 
                            yticks = 0:0.2:1, 
                            xlabel = "No. of false significant EWS (out of 17)",
                            xticks = 0:1:17,
                            #limits= (-0,6.7,nothing,nothing)
                            limits= (-0.1,6.7,nothing,nothing)
                            )
                    if ind == "sca"
                        ax.title = L"\textbf{\mathrm{\hat{𝗐}^𝟤}}" 
                        vlines!(ax, [n_x], color=cm[2], linewidth=4, label = "Number of significant EWS",alpha=0.7)
                    elseif ind == "hurst"
                        ax.title = L"\textbf{\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}}}" 
                        vlines!(ax, [n_y], color=cm[2], linewidth=4, label = "Number of significant EWS",alpha=0.7)
                    end
                    ma = Int(maximum(test_dist.num_inc_one[jj,:]))
                    hist!(ax,test_dist.num_inc_one[j,:], bins = -0.5:ma+0.5,
                                    normalization = :pdf,
                                    color = (cm[1],0.3),
                                    label = "Histogram numeric"
                                    )
                    density!(ax,test_dist.num_inc_one[jj,:], npoints = 17+2,
                            color = (:steelblue,0.0), 
                            strokecolor = cm[1],
                            strokewidth = 2, 
                            boundary = (-1,17),
                            )
                    lines!(ax, 1:5,1:5, color = cm[1], linewidth=2, label = "Distribution numeric", visible = false)

                    # vlines!(ax,[quantile(test_dist.num_inc_one[jj,:],1-p)], 
                    #     color = cm[1],
                    #     linewidth = 2, linestyle=(:dot,:dense),
                    #     label = "0.95 significance level numeric")
                    thr_num = significance_threshold(test_dist,round(1-p, digits=2))
                    vlines!(ax,[thr_num-0.5], 
                        color = cm[1],
                        linewidth = 2, linestyle=(:dot,:dense),
                        label = "0.95 significance level numeric")
                    

                    n=17
                    b = Distributions.Binomial(n,p)

                    lines!(ax,0:17, pdf.(b,0:17), color = cm[3], 
                        linewidth = 2, 
                        label ="Distribution analytic")# (B(17,$(plim)))")
                    # vlines!(ax,[quantile(b,1-p)], color = cm[3],
                    #     linestyle = :dash,
                    #     linewidth = 2, 
                    #     label = "$(1-plim) significance level analytic")
                    thr_ana = significance_threshold(b,round(1-p, digits=2))
                    vlines!(ax,[thr_ana-0.5], color = cm[3],
                        linestyle = :dash,
                        linewidth = 2, 
                        label = "$(1-plim) significance level analytic")

                    Label(ga[1,1,TopLeft()], letters[whichletter], fontsize = 20,
                        font = :bold, padding = (0,5,5,0),
                        halign = :left, valign=:bottom)

                    if k==1
                        Legend(fhist[2,2], ax,
                            framevisible=false,
                            tellwidth = false,
                            halign = :center,
                            valign = :center,
                            rowgap = 20,
                            framecolor = :grey50)
                    end
                end
            end
        end
    end
    
    gboth = fhist[2,1] = GridLayout()
    whichletter +=1
    ax3 = Axis(gboth[1,1], title = L"\textbf{\mathrm{\hat{𝗐}^𝟤} 𝖺𝗇𝖽 \mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}}}",
                            titlesize=25,
                            titlefont = :bold,
                            ylabel = "Relative frequency",
                            yticks = 0:0.2:1, 
                            xlabel = "No. of false significant EWS (out of 17)",
                            xticks = 0:1:17,
                            limits= (-0.1,6.7,nothing,nothing)
                            )
    pboth = round(plim*plim, digits=5)
    bboth = Distributions.Binomial(17,pboth)

    lines!(ax3,0:17, pdf.(bboth,0:17), color=cm[3],
                        linewidth = 2, 
                        label ="Distribution analytic (B(17,$(round(pboth, digits=4))))")
    vlines!(ax3, [n_both], color=cm[2], linewidth=4, alpha = 0.7)
    # vlines!(ax3,[quantile(bboth,1-plim)], color = cm[3],
    #             linestyle = :dash,
    #             linewidth = 2, 
    #             label = "$(1-plim) significance level analytic")
    thr_ana2 = significance_threshold(bboth,round(1-plim, digits=2))
    vlines!(ax3,[thr_ana2-0.5], color = cm[3],
                linestyle = :dash,
                linewidth = 2, 
                label = "$(1-plim) significance level analytic")

    
    Label(gboth[1,1,TopLeft()], letters[whichletter], fontsize = 20,
                    font = :bold, padding = (0,5,5,0),
                    halign = :left, valign=:bottom)
    
    linkyaxes!(fhist.content[1:3:end-1]...)
    colgap!(fhist.layout, 40)
    rowgap!(fhist.layout, 40)
    if saving
        save(saveto,fhist)
    end
    if showing
        display(fhist)
    end
end

#plot_overall_significance_example_wave(lowpass = true, plim = plim, showing = showing, saving = saving, saveto = "do_ews_across_greenland_ice_cores/figures/new_figA3.pdf")
plot_overall_significance_example_wave(lowpass = false, plim = plim, showing = showing, saving = true, saveto = "do_ews_across_greenland_ice_cores/figures/new_figS19.pdf")




#histogram of the slopes for csd (event 2 / DO-1 as example) for all cores
# Fig S3
function plot_surrogates_distribution_csd2(ylabs2;event_no = [2], lowpass = lowpass, plim = plim,
            showing = true, saving = false, 
            saveto = "paper/surrogate_distribution_csd_event_2_lowpass_$(lowpass)_p_$(plim).pdf")

    fhist = Figure(size = (1200,1800))
    letters_h = reshape(letters[1:18],2,9)

    cm = cgrad(:managua,7, categorical=true, rev=true)[[2,6,4]]

    whichletter = 0

    function good_csd_files(v,type,lowpass)
        if type == "NGRIP5"
            if lowpass
                return occursin("w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",v) && !occursin("NIKLAS_", v)
            else
                return  occursin("w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",v) && !occursin("NIKLAS_", v)
            end
        elseif type == "10y"
            if lowpass
                return occursin("w_200_normed_filt", v) && occursin("lp_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2", v)
            else
                return occursin("w_200_normed_filt", v) && occursin("gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2", v) && !occursin("_lp_", v)
            end
        elseif type == "20y"
            return occursin("w_200_normed_filt", v) && occursin("gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2", v)
        end
    end


    for (it,type) in enumerate(["NGRIP5", "10y", "20y"])
        ice_cores, var_names, vars, ac_names, acs, sca_names, scas, hurst_names, hursts, sranges = load_data(type)


        good = @. good_csd_files(var_names,type,lowpass)


        for (i,(v,a)) in enumerate(zip(reverse(vars[good]), reverse(acs[good])))
            for (ik,k) in enumerate(event_no)
                whichletter +=1
                if k>1
                    txt = " prior to DO-$(k-1) in "
                else
                    txt = " prior to YD/PB in "
                end
                for kk = 1:2
                    gl = fhist[whichletter,kk] = GridLayout()
                    if kk == 1
                        axv = Axis(gl[1,1],
                            title = L"\textbf{\mathrm{𝖵}} \textbf{\textrm{%$(sstring(txt))%$(reverse(ylabs2[it])[i]) }}", 
                            titlesize=22,
                            titlefont = :bold,
                            ylabel = "Frequency",
                            xlabel = "Linear trend (‰ yr⁻¹)")

                        if typeof(v.slopes[k]) != Nothing
                            hist!(axv,v.surr_slopes[:,k], 
                                    bins = range(0.9*minimum(v.surr_slopes[:,k]),1.1*maximum(v.surr_slopes[:,k]),50),
                                    normalization = :pdf,
                                    label = "Histogram",
                                    color = (cm[1],0.5))

                            density!(axv,convert.(Float64,v.surr_slopes[:,k]), npoints = 50,
                                    color = (:blue,0.0), 
                                    strokecolor = cm[1], 
                                    strokewidth = 2, 
                                    boundary = (0.9*minimum(v.surr_slopes[:,k]),1.1*maximum(v.surr_slopes[:,k]))
                                    )
                            lines!(axv, 1:5,1:5, color = cm[1], linewidth=2, label = "KDE", visible = false)

                            vlines!(axv,[v.slopes[k]], 
                                    color = cm[2], 
                                    linewidth = 2,
                                    label = "observed trend")
                            vlines!(axv,[quantile(v.surr_slopes[:,k],1-plim)], 
                                    color = cm[3], 
                                    linewidth = 2,
                                    linestyle = :dash,
                                    label = "$(1-plim) confidence level")

                            
                        end
                        Legend(gl[1,1], axv,
                            tellwidth = false,
                            halign = :left,
                            valign = :top,
                            margin = (10,10,10,10),
                            backgroundcolor = (:white,0.5),
                            framecolor = :grey50)
                    elseif kk ==2
                        axa = Axis(gl[1,1], 
                                title = L"\textbf{\mathrm{\alpha_𝟣}} \textbf{\textrm{%$(sstring(txt))%$(reverse(ylabs2[it])[i]) }}", 
                                titlesize=22,
                                titlefont = :bold,
                                ylabel = "Frequency",
                                xlabel = "Linear trend (‰ yr⁻¹)")#, title = reverse(ylabs[it])[i] )
                        if typeof(a.slopes[k]) != Nothing
                            hist!(axa,a.surr_slopes[:,k], 
                                    bins = range(0.9*minimum(a.surr_slopes[:,k]),1.1*maximum(a.surr_slopes[:,k]),50),
                                    normalization = :pdf,
                                    color = (cm[1],0.5),
                                    label = "Histogram")

                            density!(axa,convert.(Float64,a.surr_slopes[:,k]), npoints = 50,
                                    color = (:blue,0.0), 
                                    strokecolor = cm[1],# :orange, 
                                    strokewidth = 2, #strokearound = true,
                                    boundary = (0.9*minimum(a.surr_slopes[:,k]),1.1*maximum(a.surr_slopes[:,k]))
                                    )

                            lines!(axa, 1:5,1:5, color = cm[1], linewidth=2, label = "KDE", visible = false)
                            
                            vlines!(axa,[a.slopes[k]], color = cm[2],#:red, 
                                    linewidth = 2,
                                    label = "observed trend")
                            vlines!(axa,[quantile(a.surr_slopes[:,k],1-plim)], color = cm[3],
                                    linewidth = 2, linestyle = :dash,label = "$(1-plim) confidence level" )
                            
                        end

                        Legend(gl[1,1], axa,
                            tellwidth = false,
                            halign = :left,
                            valign = :top,
                            margin = (10,10,10,10),
                            backgroundcolor = (:white,0.5),
                            framecolor = :grey50)

                    end
                    Label(gl[1,1,TopLeft()], letters_h[kk,whichletter], fontsize = 20,
                            font = :bold, padding = (0,30,-10,15),
                            halign = :left, valign =:top)
                end
            end
        end
    end
    Label(fhist[1,1:end,Top()], "Null-model distibution of linear trends in EWS of 100-year high-pass filtered records with 5-year resolution", padding = (0.0,0.0,60.0,20.0), valign = :bottom, font = :bold, fontsize = 20)
    Label(fhist[2,1:end,Top()], "Null-model distibution of linear trends in EWS of 100-year high-pass filtered records with 10-year resolution", padding = (0.0,0.0,60.0,40.0), valign = :bottom, font = :bold, fontsize = 20)
    Label(fhist[4,1:end,Top()], "Null-model distibution of linear trends in EWS of 100-year high-pass filtered records with 20-year resolution", padding = (0.0,0.0,60.0,40.0), valign = :bottom, font = :bold, fontsize = 20)

    colgap!(fhist.layout, 20)
    rowgap!(fhist.layout, 8)
    if saving
        save(saveto,fhist)
    end
    if showing
        display(fhist)
    end
end

plot_surrogates_distribution_csd2(ylabs2,event_no = [2], lowpass = true, plim = plim, showing = showing, saving = saving, 
    saveto = "do_ews_across_greenland_ice_cores/figures/new_figS3.pdf")



#histogram of the slopes for wavelet (event 2 / DO-1 as example)
#Fig S4
function plot_surrogates_distribution_wave2(ylabs2;event_no = [2], srs = [(20,60),(20,100)], 
        lowpass = lowpass, plim = plim,
        smoothw = smoothw,
        showing = true, saving = false, 
        savetos = ["paper/surrogate_distribution_wave2_20_60_event_2_lowpass_$(lowpass)_smoothw_$(smoothw)_p_$(plim).pdf",
                "paper/surrogate_distribution_wave2_20_100_event_2_lowpass_$(lowpass)_smoothw_$(smoothw)_p_$(plim).pdf"])

    function good_wavelet_files(v,type,lowpass, smoothw,s1,s2)
        if type == "NGRIP5"
            if lowpass
                return occursin("s1_$(s1)_s2_$(s2)_C_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_$(smoothw)_10000_TFTS.jld2",v) && !occursin("NIKLAS_", v)
            else
                return occursin("s1_$(s1)_s2_$(s2)_C_no_lowpass_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_$(smoothw)_10000_TFTS.jld2",v) && !occursin("NIKLAS_", v)
            end
        elseif type == "10y"
            if lowpass
                return occursin("s1_$(s1)_s2_$(s2)", v) && occursin("_lp_PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_$(smoothw)_10000_TFTS.jld2", v)
            else
                return occursin("s1_$(s1)_s2_$(s2)", v) && occursin("PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_$(smoothw)_10000_TFTS.jld2", v) && !occursin("_lp_", v)
            end
        elseif type == "20y"
            return occursin("s1_$(s1)_s2_$(s2)", v) && occursin("PAUL_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_$(smoothw)_10000_TFTS.jld2", v)
        end
    end


    for (snr,(s1,s2)) in enumerate(srs)
        fhist = Figure(size = (1200,1800)) #(2900,1800))
        letters_h = reshape(letters[1:18],2,9)

        cm = cgrad(:managua,7, categorical=true, rev=true)[[2,6,4]]

        whichletter = 0

        for (it,type) in enumerate(["NGRIP5", "10y", "20y"])
            ice_cores, var_names, vars, ac_names, acs, sca_names, scas, hurst_names, hursts, sranges = load_data(type)

            goods = @. good_wavelet_files(sca_names,type,lowpass, smoothw,s1,s2)
            goodh = @. good_wavelet_files(hurst_names,type,lowpass, smoothw,s1,s2)

            for (i,(v,a)) in enumerate(zip(reverse(scas[goods]), reverse(hursts[goodh])))
                for (ik,k) in enumerate(event_no)#1:17
                    whichletter +=1
                    if k>1
                        txt = " prior to DO-$(k-1) in "
                    else
                        txt = " prior to YD/PB in "
                    end
                    for kk = 1:2
                        gl = fhist[whichletter,kk] = GridLayout()
                        if kk == 1
                            axv = Axis(gl[1,1], 
                                    title = L"\textbf{\mathrm{\hat{𝗐}^𝟤}} \textbf{\textrm{%$(sstring(txt))%$(reverse(ylabs2[it])[i]) }}", 
                                    titlesize=22,
                                    titlefont = :bold,
                                    ylabel = "Frequency",
                                    xlabel = "Linear trend (‰ yr⁻¹)")
                            if typeof(v.slopes[k]) != Nothing
                                hist!(axv,v.surr_slopes[:,k], 
                                        bins = range(0.9*minimum(v.surr_slopes[:,k]),1.1*maximum(v.surr_slopes[:,k]),50),
                                        normalization = :pdf,
                                        label = "Histogram",
                                        color = (cm[1],0.5))

                                density!(axv,convert.(Float64,v.surr_slopes[:,k]), npoints = 50,
                                        color = (:blue,0.0), 
                                        strokecolor = cm[1], strokewidth = 2, #strokearound = true,
                                        boundary = (0.9*minimum(v.surr_slopes[:,k]),1.1*maximum(v.surr_slopes[:,k]))
                                        )
                                lines!(axv, 1:5,1:5, color = cm[1], linewidth=2, label = "KDE", visible = false)
                                
                                vlines!(axv,[v.slopes[k]], color = cm[2], 
                                        linewidth = 2,
                                        label = "observed trend")
                                vlines!(axv,[quantile(v.surr_slopes[:,k],1-plim)], color = cm[3], 
                                        linewidth = 2,
                                        linestyle = :dash,
                                        label = "$(1-plim) confidence level")
                            end
                            Legend(gl[1,1], axv,
                                        tellwidth = false,
                                        halign = :left,
                                        valign = :top,
                                        margin = (10,10,10,10),
                                        backgroundcolor = (:white,0.5),
                                        framecolor = :grey50)
                        elseif kk == 2
                            if typeof(a.slopes[k]) != Nothing
                                axa = Axis(gl[1,1], 
                                title = L"\textbf{\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}}} \textbf{\textrm{%$(sstring(txt))%$(reverse(ylabs2[it])[i]) }}", 
                                titlesize=22,
                                titlefont = :bold,
                                ylabel = "Frequency",
                                xlabel = "Linear trend (‰ yr⁻¹)")
                                hist!(axa,a.surr_slopes[:,k], 
                                        bins = range(0.9*minimum(a.surr_slopes[:,k]),1.1*maximum(a.surr_slopes[:,k]),50),
                                        normalization = :pdf,
                                        color = (cm[1],0.5),
                                        label = "Histogram")

                                density!(axa,convert.(Float64,a.surr_slopes[:,k]), npoints = 50,
                                        color = (:blue,0.0), 
                                        strokecolor = cm[1], strokewidth = 2,
                                        boundary = (0.9*minimum(a.surr_slopes[:,k]),1.1*maximum(a.surr_slopes[:,k]))
                                        )
 
                                vlines!(axa,[a.slopes[k]], color = cm[2],#:red, 
                                        linewidth = 2,
                                        label = "observed trend")
                                vlines!(axa,[quantile(a.surr_slopes[:,k],1-plim)], color = cm[3], 
                                        linewidth = 2, linestyle=:dash, label = "$(1-plim) confidence level" )
                            end
                            Legend(gl[1,1], axa,
                                        tellwidth = false,
                                        halign = :left,
                                        valign = :top,
                                        margin = (10,10,10,10),
                                        backgroundcolor = (:white,0.5),
                                        framecolor = :grey50)

                        end
                        Label(gl[1,1,TopLeft()], letters_h[kk,whichletter], fontsize = 20,
                                font = :bold, padding = (0,5,-10,15),
                                halign = :left, valign =:top)
                    end
                end
            end
        end
        Label(fhist[1,1:end,Top()], "Null-model distibution of linear trends in EWS of ($(s1)-$(s2)) year band of records with 5-year resolution", padding = (0.0,0.0,60.0,20.0), valign = :bottom, font = :bold, fontsize = 20)
        Label(fhist[2,1:end,Top()], "Null-model distibution of linear trends in EWS of ($(s1)-$(s2)) year band of records with 10-year resolution", padding = (0.0,0.0,60.0,40.0), valign = :bottom, font = :bold, fontsize = 20)
        Label(fhist[4,1:end,Top()], "Null-model distibution of linear trends in EWS of ($(s1)-$(s2)) year band of records with 20-year resolution", padding = (0.0,0.0,60.0,40.0), valign = :bottom, font = :bold, fontsize = 20)

        colgap!(fhist.layout, 20)
        rowgap!(fhist.layout, 8)
        
        if saving
            save(savetos[snr],fhist)
        end
        if showing
            display(fhist)
        end
    end
end

plot_surrogates_distribution_wave2(ylabs2;event_no = [2], srs = [(20,60)], lowpass = true, plim = plim,smoothw = false,
        showing = showing, saving = saving ,savetos = ["do_ews_across_greenland_ice_cores/figures/new_figS4.pdf"])




begin
    wavepal_path = "do_ews_across_greenland_ice_cores/new_surrogate_files/NGRIP_irreg/indicators_NGRIP_5y_c_wavepal_theta_time_sranges_10_0_50_tfts_ns_1000.npz"
    wavepal_file = npzread(wavepal_path)
end



#plot EWS for irregular time stepping (from wavepal)
#Fig 5 
function plot_wavepal_csd_ts_label(;plim = plim, legend = false,showing = true, saving = false, saveto = "paper/wavepal_csd_ts2_ngrip_p_$(plim).pdf")
    global wavepal_file
    times = wavepal_file["thetas"]
    ff = Figure(size=(1200,400))
    n_v = 0 
    n_a = 0
    n_both = 0
    which_v = []
    for k = 1:2
        gl = ff[1,k] = GridLayout()
        if k == 1
            ind = wavepal_file["stds"][:,1,:]
            true_slopes = wavepal_file["slopes_sigma"][:,1]
            surr_slopes = wavepal_file["slopes_sigma"][:,2:end]
            ind_name = L"\mathrm{𝖵}"
            ind_title = L"\textbf{\mathrm{𝖵}}"
        elseif k == 2
            ind = wavepal_file["acs"][:,1,:]
            true_slopes = wavepal_file["slopes_alpha"][:,1]
            surr_slopes = wavepal_file["slopes_alpha"][:,2:end]
            ind_name = L"\mathrm{\hat{\alpha}_𝟣}"
            ind_title = L"\textbf{\mathrm{\hat{\alpha}_𝟣}}"
        end
        axv = Axis(gl[1,1], 
                xlabel = "Age (kyr b2k)", 
                ylabel = ind_name, 
                ylabelsize = 18,
                ylabelfont = :bold,
                xticks = (10_000:5_000:60_000, string.(10:5:60)))
        axv.xreversed = true
        ax_event = Axis(gl[1,1],xaxisposition = :top, xticklabelrotation = pi*0.4)
        ax_event.xticks = (GI_onsets, event_labels)
        ax_event.xreversed = true
        hidespines!(ax_event)
        hideydecorations!(ax_event)
        for ev in 1:17
            va = 0.4
            vl = 2.0
            vsc = (:white,1.0)
            vc = :black
            if typeof(ind[ev,:]) != Nothing
                if true_slopes[ev] >0 
                    vc = :red
                    vsc = (:darkred,0.1)
                    pone = count(v -> v ≥ true_slopes[ev], surr_slopes[ev,:])/length(surr_slopes[ev,:])
                    if pone < plim
                        va = 0.9
                        vl = 3.5
                        vsc = (:darkred,0.8)
                        if k==1
                            n_v +=1
                            push!(which_v, ev)
                        elseif k==2
                            n_a +=1
                            if ev in which_v
                                n_both +=1
                            end
                        end
                    end
                else
                    vc = :blue
                    vsc = (:steelblue,0.1)
                end
            end
            vspan!(axv, [GI_onsets[ev]], [GS_onsets[ev]], color = vsc)
            lines!(axv, times[ev,:], ind[ev,:], color = :black, alpha = va)
            pv = Polynomials.fit(times[ev,:][findall(!isnan, ind[ev,:])], convert.(Float64,ind[ev,:][findall(!isnan, ind[ev,:])]),1)
            lines!(axv, times[ev,:], pv.(times[ev,:]), color = vc, alpha = va, linewidth = vl)
            vlines!(axv, [GI_onsets[ev]], color = :darkred, alpha = va, linewidth = vl)

        end
        xlims!(axv,nothing,nothing)
        axv.xreversed = true
        xlims!(ax_event, axv.xaxis.attributes.limits[]...)
        ax_event.xreversed = true

        Label(gl[1,1,TopLeft()], letters[k], fontsize = 20,
                font = :bold, padding = (0,20,-10,15),
                halign = :left, valign =:bottom)
    end
    if legend
        elems = [
            #[LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(0.5, 0), (0.5, 1)])],
            [LineElement(color= (:black, 0.4), points = Point2f[(0, 0.5), (0.5, 0.5)]), LineElement(color= (:black, 0.9), points = Point2f[(0.5, 0.5), (1, 0.5)])],
            [PolyElement(color = (:darkred, 0.1), strokewidth = 0), LineElement(color= (:red,0.4), linewidth=2.0)],# LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:steelblue,0.1), strokewidth = 0), LineElement(color= (:blue, 0.4), linewidth=2.0)],#,LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:darkred, 0.8), strokewidth = 0), LineElement(color= (:red,0.9), linewidth=3.5)],#,LineElement(color = (:darkred,0.9), linewidth= 3.5,points = Point2f[(1, 0), (1, 1)])]]
            ]
        labels = [#"GI onset", 
                "EWS indicator", "increasing", "decreasing", "significantly increasing (𝗉<0.05)"]
        Legend(ff[end+1,1:end], elems, labels, orientation = :horizontal, valign = :bottom, margin = (10,10,10,10)) 
    end
    colgap!(ff.layout, 20)
    rowgap!(ff.layout, 0)
    txt = "Irregular temporal resolution"
    Label(ff[1,1:end,Top()], "EWS in 100-year high-pass filtered NGRIP record", padding = (0.0,0.0,110.0,30.0), valign = :bottom, font = :bold, fontsize = 20)
    Label(ff[1,1:end,Top()], L"\textbf{\textrm{%$(sstring(txt)):} \mathrm{𝗇_{𝖵} = %$(sstring(n_v)),} \mathrm{𝗇_{\hat{\alpha}_𝟣} = %$(sstring(n_a)),} \mathrm{𝗇_{\text{𝖻𝗈𝗍𝗁}} = %$(sstring(n_both))}}", padding = (0.0,0.0,20.0,20.0), valign = :center, font = :bold, fontsize = 20)
    if saving
        save(saveto,ff)
    end
    if showing
        display(ff)
    end
end


plot_wavepal_csd_ts_label(plim=0.05,legend = true,showing = showing, saving = saving, saveto = "do_ews_across_greenland_ice_cores/figures/fig5.pdf")


function plot_wavepal_csd_ts_label_aggr(;plim = plim, legend = false,showing = true, saving = false, saveto = "paper/wavepal_csd_ts2_ngrip_p_$(plim).pdf")
    
    global wavepal_file

    fts = Figure(size=(1200, 370))
    
        n_v = 0 
        n_a = 0
        n_both = 0
        which_v = []
        
        for k = 1:2
            gl = fts[1,k] = GridLayout()
            
            if k == 1
                ind = wavepal_file["stds"][:,1,:]
                true_slopes = wavepal_file["slopes_sigma"][:,1]
                surr_slopes = wavepal_file["slopes_sigma"][:,2:end]
                ind_name = L"\mathrm{𝖵}"
                ind_title = L"\textbf{\mathrm{𝖵}}"
            elseif k == 2
                ind = wavepal_file["acs"][:,1,:]
                true_slopes = wavepal_file["slopes_alpha"][:,1]
                surr_slopes = wavepal_file["slopes_alpha"][:,2:end]
                ind_name = L"\mathrm{\hat{\alpha}_𝟣}"
                ind_title = L"\textbf{\mathrm{\hat{\alpha}_𝟣}}"
            end
            ax = Axis(gl[1,1],
                xticks = (1:17, ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]) ,
                xreversed = true,
                xticklabelrotation = pi*0.4,
                #yreversed=true,
                ylabel =  ind_name,

                xlabel = "Event",
                xminorgridwidth = 1,
                xgridwidth=0,
                xminorgridcolor = :grey10 ,
                xgridcolor=:transparent,
                xminorgridvisible = true, 
                xgridvisible = false,
                #yticks = (1:6, ["NGRIP 5 y", "NGRIP 10 y", "NEEM 10 y", "NGRIP 20 y", "GRIP 20 y", "GISP2 20 y"]),
                #ylabel = "Ice core",
                #yminorgridwidth = 1,
                #ygridwidth=0.7,
                #yminorgridcolor = :grey10 ,
                #ygridcolor=:grey30,
                yminorgridvisible = false, 
                ygridvisible = false,
                yticksvisible=false,
                yticklabelsvisible=false,
                ylabelsize = 18,
                ylabelfont = :bold,
                )


            ma = zeros(17,1)
            for ev in 1:17
                #va = 0.4
                #vl = 2.0
                #vsc = (:white,1.0)
                #vc = :black
                if typeof(ind[ev,:]) != Nothing
                    if true_slopes[ev] >0 
                        ma[ev,1]+=1
                        pone = count(v -> v ≥ true_slopes[ev], surr_slopes[ev,:])/length(surr_slopes[ev,:])
                #        vc = :red
                #        vsc = (:darkred,0.1)
                        if pone < plim
                            ma[ev,1]+=1
                #            va = 0.9
                #            vl = 3.5
                #            vsc = (:darkred,0.8)
                            if k==1
                                n_v +=1
                                push!(which_v, ev)
                            elseif k==2
                                n_a +=1
                                if ev in which_v
                                    n_both +=1
                                end
                            end
                        end
                    end
                end
                # vspan!(axv, [GI_onsets[ev]], [GS_onsets[ev]], color = vsc)
                # pv = Polynomials.fit(v.times[ev][findall(!ismissing, v.vals[ev])], convert.(Float64,v.vals[ev][findall(!ismissing, v.vals[ev])]),1)
                # lines!(axv,v.times[ev], v.vals[ev], color = :black, alpha = va)#, linewidth = vl)
                # lines!(axv, v.times[ev], pv.(v.times[ev]), color = vc, alpha = va, linewidth = vl)
                # vlines!(axv, [GI_onsets[ev]], color = :darkred, alpha = va, linewidth = vl)
            end
            co2 = cgrad([:steelblue, :darkred, :darkred], 3, categorical = true, alpha=[0.1,0.1,0.8])
            hm = heatmap!(ax, ma, colormap=co2, colorrange=(0,2))
            # xlims!(axv,nothing,nothing)
            # axv.xreversed = true
            # xlims!(ax_event, axv.xaxis.attributes.limits[]...)
    


        Label(gl[1,1,TopLeft()], letters[k], fontsize = 20,
                    font = :bold, padding = (0,20,-10,15),
                    halign = :left, valign =:bottom)
    end
    

    
    
    if legend
        # elems = [
        #     [LineElement(color= (:black, 0.4), points = Point2f[(0, 0.5), (0.5, 0.5)]), LineElement(color= (:black, 0.9), points = Point2f[(0.5, 0.5), (1, 0.5)])],
        #     [PolyElement(color = (:darkred, 0.1), strokewidth = 0), LineElement(color= (:red,0.4), linewidth=2.0)],# LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
        #     [PolyElement(color = (:steelblue,0.1), strokewidth = 0), LineElement(color= (:blue, 0.4), linewidth=2.0)],#,LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
        #     [PolyElement(color = (:darkred, 0.8), strokewidth = 0), LineElement(color= (:red,0.9), linewidth=3.5)],#,LineElement(color = (:darkred,0.9), linewidth= 3.5,points = Point2f[(1, 0), (1, 1)])]]
        #     ]
        # labels = ["EWS indicator", "increasing", "decreasing", "significantly increasing (𝗉<0.05)"]
        elems = [
            [PolyElement(color = (:steelblue,0.1), strokewidth = 0)],#,LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:darkred, 0.1), strokewidth = 0)],# LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:darkred, 0.8), strokewidth = 0)],#,LineElement(color = (:darkred,0.9), linewidth= 3.5,points = Point2f[(1, 0), (1, 1)])]]
            ]
        labels = ["increasing", "decreasing", "significantly increasing (𝗉<0.05)"]
        Legend(fts[end+1,1:end], elems, labels, orientation = :horizontal, valign = :bottom, margin = (10,10,10,10)) 
    end

    colgap!(fts.layout, 20)
    rowgap!(fts.layout, 0)

    txt = "Irregular temporal resolution"
    Label(fts[1,1:end,Top()], "EWS in 100-year high-pass filtered NGRIP record", padding = (0.0,0.0,80.0,3.0), valign = :bottom, font = :bold, fontsize = 20)
    Label(fts[1,1:end,Top()], L"\textbf{\textrm{%$(sstring(txt)):} \mathrm{𝗇_{𝖵} = %$(sstring(n_v)),} \mathrm{𝗇_{\hat{\alpha}_𝟣} = %$(sstring(n_a)),} \mathrm{𝗇_{\text{𝖻𝗈𝗍𝗁}} = %$(sstring(n_both))}}", padding = (0.0,0.0,20.0,20.0), valign = :center, font = :bold, fontsize = 20)


    if saving
        save(saveto,fts)
    end
    if showing
        display(fts)
    end
end

plot_wavepal_csd_ts_label_aggr(plim=0.05,legend = true,showing = showing, saving = saving, saveto = "do_ews_across_greenland_ice_cores/figures/fig5_aggr.pdf")

#not used:
function plot_wavepal_wave_ts_label(;plim = plim, legend = false, showing = true, saving = true, saveto = "paper/wavepal_wave_ts2_10_50_ngrip_p_$(plim).pdf")
    global wavepal_file
    times = wavepal_file["thetas"]
    ff = Figure(size=(1200,400))
    n_v = 0 
    n_a = 0
    n_both = 0
    which_v = []
    for k = 1:2
        gl = ff[1,k] = GridLayout()
        if k == 1
            ind = wavepal_file["scavs"][:,1,1,:]
            true_slopes = wavepal_file["slopes_w"][:,1,1]
            surr_slopes = wavepal_file["slopes_w"][:,2:end,1]
            ind_name = L"\mathrm{\hat{𝗐}^𝟤}"
            ind_title = L"\textbf{\mathrm{\hat{𝗐}^𝟤}}"
        elseif k == 2
            ind = wavepal_file["hursts"][:,1,1,:]
            true_slopes = wavepal_file["slopes_h"][:,1,1]
            surr_slopes = wavepal_file["slopes_h"][:,2:end,1]
            ind_name = L"\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}}"
            ind_title = L"\textbf{\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}}}"
        end
        axv = Axis(gl[1,1], 
                xlabel = "Age (kyr b2k)", 
                ylabel = ind_name, 
                ylabelsize = 18,
                ylabelfont = :bold,
                xticks = (10_000:5_000:60_000, string.(10:5:60)))
        axv.xreversed = true
        ax_event = Axis(gl[1,1],xaxisposition = :top, xticklabelrotation = pi*0.4)
        ax_event.xticks = (GI_onsets, event_labels)
        ax_event.xreversed = true
        hidespines!(ax_event)
        hideydecorations!(ax_event)
        for ev in 1:17
            va = 0.4
            vl = 2.0
            vsc = (:white,1.0)
            vc = :black
            if typeof(ind[ev,:]) != Nothing
                if true_slopes[ev] >0 
                    vc = :red
                    vsc = (:darkred,0.1)
                    pone = count(v -> v ≥ true_slopes[ev], surr_slopes[ev,:])/length(surr_slopes[ev,:])
                    if pone < plim
                        va = 0.9
                        vl = 3.5
                        vsc = (:darkred,0.8)
                        if k==1
                            n_v +=1
                            push!(which_v, ev)
                        elseif k==2
                            n_a +=1
                            if ev in which_v
                                n_both +=1
                            end
                        end 
                    end
                else
                    vc = :blue
                    vsc = (:steelblue,0.1)
                end
            end
            vspan!(axv, [GI_onsets[ev]], [GS_onsets[ev]], color = vsc)
            lines!(axv, times[ev,:], ind[ev,:], color = :black, alpha = va)
            pv = Polynomials.fit(times[ev,:][findall(!isnan, ind[ev,:])], convert.(Float64,ind[ev,:][findall(!isnan, ind[ev,:])]),1)
            lines!(axv, times[ev,:], pv.(times[ev,:]), color = vc, alpha = va, linewidth = vl)
            vlines!(axv, [GI_onsets[ev]], color = :darkred, alpha = va, linewidth = vl)
        end
        xlims!(axv,nothing,nothing)
        axv.xreversed = true
        xlims!(ax_event, axv.xaxis.attributes.limits[]...)
        ax_event.xreversed = true

        Label(gl[1,1,TopLeft()], letters[k], fontsize = 20,
                font = :bold, padding = (0,5,-10,15),
                halign = :left, valign =:bottom)
    end
    if legend
        elems = [
            #[LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(0.5, 0), (0.5, 1)])],
            [LineElement(color= (:black, 0.4), points = Point2f[(0, 0.5), (0.5, 0.5)]), LineElement(color= (:black, 0.9), points = Point2f[(0.5, 0.5), (1, 0.5)])],
            #[LineElement(color= (:black, 0.4), points = Point2f[(0, 0), (1, 0)]), LineElement(color= (:black, 0.9), points = Point2f[(0, 1), (1, 1)])],
            [PolyElement(color = (:darkred, 0.1), strokewidth = 0), LineElement(color= (:red,0.4), linewidth=2.0)],# LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:steelblue,0.1), strokewidth = 0), LineElement(color= (:blue, 0.4), linewidth=2.0)],#,LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:darkred, 0.8), strokewidth = 0), LineElement(color= (:red,0.9), linewidth=3.5)],#,LineElement(color = (:darkred,0.9), linewidth= 3.5,points = Point2f[(1, 0), (1, 1)])]]
            ]
        labels = [#"GI onset", 
                "EWS indicator", "increasing", "decreasing", "significantly increasing (𝗉<0.05)"]
        Legend(ff[end+1,1:end], elems, labels, orientation = :horizontal, valign = :bottom, margin = (10,10,10,10)) 
    end
    colgap!(ff.layout, 20)
    rowgap!(ff.layout, 0)
    txt = "Irregular temporal resolution"
    Label(ff[1,1:end,Top()], "EWS in (10-50) year band of NGRIP record", padding = (0.0,0.0,110.0,30.0), valign = :bottom, font = :bold, fontsize = 20)
    Label(ff[1,1:end,Top()], L"\textbf{\textrm{%$(sstring(txt)):} \mathrm{𝗇_{\hat{𝗐}^𝟤} = %$(sstring(n_v)),} \mathrm{𝗇_{\hat{𝖧}^\text{𝗅𝗈𝖼}} = %$(sstring(n_a)),} \mathrm{𝗇_{\text{𝖻𝗈𝗍𝗁}} = %$(sstring(n_both))}}", padding = (0.0,0.0,20.0,20.0), valign = :center, font = :bold, fontsize = 20)
    if saving
        save(saveto,ff)
    end
    if showing
        display(ff)
    end
    
end

# plot_wavepal_wave_ts_label(legend = true,showing=showing, saving = saving, saveto = "paper/figures/new/wavepal_wave_ts2_label_10_50_ngrip_p_$(plim).pdf")



# Fig S5
function plot_wavepal_surrogates_distribution_csd(;event_no = 2, plim = plim,
                                showing = true, saving = false, 
                                saveto = "paper/wavepal_surrogate_distribution_csd_event_2_p_$(plim).pdf")
    fhist = Figure(size = (1200,400))
    cm = cgrad(:managua,7, categorical=true, rev=true)[[2,6,4]]
    global wavepal_file
    for k = 1:2
        gl = fhist[1,k] = GridLayout()
        if k == 1
            true_slopes = wavepal_file["slopes_sigma"][:,1]
            surr_slopes = wavepal_file["slopes_sigma"][:,2:end]
            ind_title = L"\textbf{\mathrm{𝖵}}"
        elseif k == 2
            true_slopes = wavepal_file["slopes_alpha"][:,1]
            surr_slopes = wavepal_file["slopes_alpha"][:,2:end]
            ind_title = L"\textbf{\mathrm{\hat{\alpha}_𝟣}}"
        end
        axv = Axis(gl[1,1], 
            title = ind_title, 
            titlesize=25,
            titlefont = :bold,
            ylabel = "Frequency",
            xlabel = "Linear trend (‰/yr)")
        
        hist!(axv,surr_slopes[event_no,:], 
                bins = range(0.9*minimum(surr_slopes[event_no,:]),1.1*maximum(surr_slopes[event_no,:]),50),
                normalization = :pdf,
                label = "Histogram",
                color = (cm[1],0.5))
        
        density!(axv,convert.(Float64,surr_slopes[event_no,:]), npoints = 50,
                color = (:blue,0.0), 
                strokecolor = cm[1], strokewidth = 2, #strokearound = true,
                boundary = (0.9*minimum(surr_slopes[event_no,:]),1.1*maximum(surr_slopes[event_no,:]))
                )
        lines!(axv, 1:5,1:5, color = cm[1], linewidth=2, label = "KDE", visible = false)
        
        vlines!(axv,[true_slopes[event_no]], color = cm[2], 
                linewidth = 2,
                label = "observed trend")
        vlines!(axv,[quantile(surr_slopes[event_no,:],1-plim)], color = cm[3], linewidth = 2, linestyle=:dash, label = "$(1-plim) confidence level" )

        Label(gl[1,1,TopLeft()], letters[k], fontsize = 20,
                font = :bold, padding = (0,5,5,0),
                halign = :left)

        Legend(gl[1,1], axv,
                tellwidth = false,
                halign = :left,
                valign = :top,
                margin = (10,10,10,10),
                backgroundcolor = (:white,0.5),
                framecolor = :grey50)

    end
    colgap!(fhist.layout, 40)
    if saving
        save(saveto,fhist)
    end
    if showing
        display(fhist)
    end
end

plot_wavepal_surrogates_distribution_csd(event_no = 2, plim = plim, showing = showing, saving = saving, 
                saveto = "do_ews_across_greenland_ice_cores/figures/new_figS5.pdf")



#Fig S6
function plot_wavepal_surrogates_distribution_wave(;event_no = 2,plim = plim,
                        showing = true, saving = false, 
                        saveto = "paper/wavepal_surrogate_distribution_wave_10_50_event_2_p_$(plim).pdf")
    fhist = Figure(size = (1200,400))
    cm = cgrad(:managua,7, categorical=true, rev=true)[[2,6,4]]
    global wavepal_file
    for k = 1:2
        gl = fhist[1,k] = GridLayout()
        if k == 1
            true_slopes = wavepal_file["slopes_w"][:,1,1]
            surr_slopes = wavepal_file["slopes_w"][:,2:end,1]
            ind_title = L"\textbf{\mathrm{\hat{𝗐}^𝟤}}"
        elseif k == 2
            true_slopes = wavepal_file["slopes_h"][:,1,1]
            surr_slopes = wavepal_file["slopes_h"][:,2:end,1]
            ind_title = L"\textbf{\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}}}"
        end
        axv = Axis(gl[1,1], 
            title = ind_title, 
            titlesize=25,
            titlefont = :bold,
            ylabel = "Frequency",
            xlabel = "Linear trend (‰/yr)")
        
        hist!(axv,surr_slopes[event_no,:],
                bins = range(0.9*minimum(surr_slopes[event_no,:]),1.1*maximum(surr_slopes[event_no,:]),50),
                normalization = :pdf,
                label = "Histogram",
                color = (cm[1],0.5))
        density!(axv,convert.(Float64,surr_slopes[event_no,:]), npoints = 50,
                color = (:blue,0.0), 
                strokecolor = cm[1], strokewidth = 2, 
                boundary = (0.9*minimum(surr_slopes[event_no,:]),1.1*maximum(surr_slopes[event_no,:]))
                )
        lines!(axv, 1:5,1:5, color = cm[1], linewidth=2, label = "KDE", visible = false)
        
        
        vlines!(axv,[true_slopes[event_no]], color = cm[2], 
                linewidth = 2,
                label = "observed trend")
        vlines!(axv,[quantile(surr_slopes[event_no,:],1-plim)], color = cm[3], 
                linewidth = 2,
                linestyle=:dash,
                label = "$(1-plim) confidence level")

        Label(gl[1,1,TopLeft()], letters[k], fontsize = 20,
                font = :bold, padding = (0,5,5,0),
                halign = :left)

        Legend(gl[1,1], axv,
                tellwidth = false,
                halign = :left,
                valign = :top,
                margin = (10,10,10,10),
                backgroundcolor = (:white,0.5),
                framecolor = :grey50)

    end
    colgap!(fhist.layout, 40)
    if saving
        save(saveto,fhist)
    end
    if showing
        display(fhist)
    end
end

plot_wavepal_surrogates_distribution_wave(event_no = 2,plim = plim, showing = showing, saving = saving, 
                        saveto = "do_ews_across_greenland_ice_cores/figures/new_figS6.pdf")



 
## make wavepal data into object like the others
begin 
    wavepal_indicators = indicator_and_significance[]
    wtimes = wavepal_file["thetas"]
    resolution = 1
    name = "NGRIP wavepal"
   for (i,na) in enumerate(["stds", "acs", "scavs", "hursts"])
        times = Array{Vector{Float64}}(undef, 17)
        vals = Array{Vector{Union{Missing,Float64}}}(undef, 17)
        slopes = Array{Union{Nothing,Float64}}(nothing, 17)
        surr_slopes = Array{Union{Nothing,Float64}}(nothing, 1_000, 17)
        p_one = Array{Union{Nothing,Float64}}(nothing, 17)
        p_two= Array{Union{Nothing,Float64}}(nothing, 17)
       
        if na == "stds"
            wvals = wavepal_file[na][:,1,:]
            wslopes = wavepal_file["slopes_sigma"][:,1]
            wsurr_slopes = wavepal_file["slopes_sigma"][:,2:end]
            wtype = "Var"
            ws1=0
            ws2 = 100
            saveto = "new_surrogate_files/NGRIP_irreg/var_filt_1000_tfts.jld2"
        elseif na == "acs"
            wvals = wavepal_file[na][:,1,:]
            wslopes = wavepal_file["slopes_alpha"][:,1]
            wsurr_slopes = wavepal_file["slopes_alpha"][:,2:end]
            wtype = "AC"
            ws1=0
            ws2 = 100
            saveto = "new_surrogate_files/NGRIP_irreg/ac_filt_1000_tfts.jld2"
        elseif na == "scavs"
            wvals = wavepal_file[na][:,1,1,:]
            wslopes = wavepal_file["slopes_w"][:,1,1]
            wsurr_slopes = wavepal_file["slopes_w"][:,2:end,1]
            wtype = "sca"
            ws1 = 10
            ws2 = 50
            saveto = "new_surrogate_files/NGRIP_irreg/sca_10_50_1000_tfts.jld2"
        elseif na == "hursts"
            wvals = wavepal_file[na][:,1,1,:]
            wslopes = wavepal_file["slopes_h"][:,1,1]
            wsurr_slopes = wavepal_file["slopes_h"][:,2:end,1]
            wtype = "hurst"
            ws1 = 10
            ws2 = 50
            saveto = "new_surrogate_files/NGRIP_irreg/hurst_10_50_1000_tfts.jld2"
        end
        for ev in 1:17
            times[ev] = wtimes[ev,:]
            vals[ev] = wvals[ev,:]
            replace!(x->isnan(x) ? missing : x, vals[ev])
            slopes[ev] = wslopes[ev]
            surr_slopes[:,ev] = wsurr_slopes[ev,:]
            p_one[ev] = count(v -> v ≥ wslopes[ev], wsurr_slopes[ev,:])/length(wsurr_slopes[ev,:])
            pr = count(v -> v ≥ wslopes[ev], wsurr_slopes[ev,:])
            pl = count(v -> v ≤ wslopes[ev], wsurr_slopes[ev,:])
            p_two[ev]= 2min(pr, pl)/length(wsurr_slopes[ev,:])
        end

        ind_obj = indicator_and_significance(
            name,
            times,    
            vals, 
            slopes,
            surr_slopes,
            wtype,
            ws1,
            ws2,
            p_one,
            p_two,
            resolution)
        push!(wavepal_indicators, ind_obj)

        
        if !isfile(saveto)    
            save(saveto, "slopes", ind_obj)
        end
   end
end



#not used:
function plot_wavepal_compare_csd_ts_label(;lowpass = false, plim = plim, legend = false, showing = true, saving = true, saveto = "paper/wavepal_compare_csd_ts_label_ngrip_lowpass_$(lowpass)_p_$(plim).pdf")
    f_csd_ts = Figure(size=(1200,600))
    
    type = "NGRIP5"
    letters_h = reshape(letters[1:18],2,9)

    whichletter = 0
    
    _, var_names, vars, ac_names, acs, _ = load_data(type)

    function good_csd_file_ngrip(v, lowpass) 
        if lowpass
            return occursin("w_200_normed_filt_C_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",v) && !occursin("NIKLAS_", v)
        else
            return occursin("w_200_normed_filt_C_no_lowpass_gs_FiltInd_false_onlyfull_true_10000_TFTS.jld2",v) && !occursin("NIKLAS_", v)
        end
    end


    good = @. good_csd_file_ngrip(var_names, lowpass)

    version_names = ["Regular 5-year resolution","Irregular temporal resolution"]
    for (i,(v,a)) in enumerate(zip([vars[good]...,wavepal_indicators[1]], [acs[good]..., wavepal_indicators[2]]))
        glv = f_csd_ts[i,1] = GridLayout()
        gla = f_csd_ts[i,2] = GridLayout()
        axv = Axis(glv[1,1], 
                ylabel = L"\mathrm{𝖵}", 
                ylabelsize = 18,
                ylabelfont = :bold,
                xticks = (10_000:5_000:60_000, string.(10:5:60)), 
                )
        axa = Axis(gla[1,1], 
                ylabel = L"\mathrm{\hat{\alpha}_𝟣}", 
                ylabelsize = 18,
                ylabelfont = :bold,
                xticks = (10_000:5_000:60_000, string.(10:5:60)), 
                )
        axv.xreversed = true
        axa.xreversed = true

        axv_event = Axis(glv[1,1],xaxisposition = :top, xticklabelrotation = pi*0.4)
        axv_event.xticks = (GI_onsets, event_labels)
        axv_event.xreversed = true
        hidespines!(axv_event)
        hideydecorations!(axv_event)

        axa_event = Axis(gla[1,1],xaxisposition = :top, xticklabelrotation = pi*0.4)
        axa_event.xticks = (GI_onsets, event_labels)
        axa_event.xreversed = true
        hidespines!(axa_event)
        hideydecorations!(axa_event)

        whichletter +=1

        n_v = 0 
        n_a = 0
        n_both = 0
        which_v = []

        if whichletter == 2
            axv.xlabel = "Age (kyr b2k)"
            axa.xlabel = "Age (kyr b2k)"
        else
            hidexdecorations!(axv, grid=false)
            hidespines!(axv, :b) 
            hidexdecorations!(axa, grid=false)
            hidespines!(axa, :b) 
        end

        if i>1
            hidespines!(axv,:t)
            hidespines!(axa,:t)
            hidexdecorations!(axv_event, grid=false)
            hidexdecorations!(axa_event, grid=false)
        end

        for ev in 1:17
            va = 0.4
            vl = 2.0
            vsc = (:white,1.0)
            vc = :black
            if typeof(v.slopes[ev]) != Nothing
                if v.slopes[ev] >0 
                    vc = :red
                    vsc = (:darkred,0.1)
                    if v.p_one[ev] < plim
                        va = 0.9
                        vl = 3.5
                        vsc = (:darkred,0.8)
                        n_v +=1
                        push!(which_v, ev)
                    end
                
                else
                    vc = :blue
                    vsc = (:steelblue,0.1)
                end
            end
            vspan!(axv, [GI_onsets[ev]], [GS_onsets[ev]], color = vsc)
            pv = Polynomials.fit(v.times[ev][findall(!ismissing, v.vals[ev])], convert.(Float64,v.vals[ev][findall(!ismissing, v.vals[ev])]),1)
            lines!(axv,v.times[ev], v.vals[ev], color = :black, alpha = va)#, linewidth = vl)
            lines!(axv, v.times[ev], pv.(v.times[ev]), color = vc, alpha = va, linewidth = vl)
            vlines!(axv, [GI_onsets[ev]], color = :darkred, alpha = va, linewidth = vl)

            aa = 0.4
            al = 2.0
            asc = (:white,1.0)
            ac = :black
            if typeof(a.slopes[ev]) != Nothing
                if a.slopes[ev] >0 
                    ac = :red
                    asc = (:darkred,0.1)
                    if a.p_one[ev] < plim
                        aa = 0.9
                        al = 3.5
                        asc = (:darkred,0.8)
                        n_a +=1
                        if ev in which_v
                            n_both +=1
                        end
                    end
                else
                    ac = :blue
                    asc = (:steelblue,0.1)
                end
            end
            vspan!(axa, [GI_onsets[ev]], [GS_onsets[ev]], color = asc)
            pa = Polynomials.fit(a.times[ev][findall(!ismissing, a.vals[ev])], convert.(Float64,a.vals[ev][findall(!ismissing, a.vals[ev])]),1)
            lines!(axa,a.times[ev], a.vals[ev], color = :black, alpha = aa)#, linewidth = vl)
            lines!(axa, a.times[ev], pa.(a.times[ev]), color = ac, alpha = aa, linewidth = al)
            vlines!(axa, [GI_onsets[ev]], color = :darkred, alpha = aa, linewidth = al)
            
            xlims!(axv,nothing,nothing)
            axv.xreversed = true
            xlims!(axv_event, axv.xaxis.attributes.limits[]...)
            axv_event.xreversed = true
            xlims!(axa_event, axv.xaxis.attributes.limits[]...)
            axa_event.xreversed = true
        end
        Label(glv[1,1,TopLeft()], letters_h[1,whichletter], fontsize = 20,
                font = :bold, padding = (0,30,-10,15),
                halign = :left,valign =:bottom)
        Label(gla[1,1,TopLeft()], letters_h[2,whichletter], fontsize = 20,
                font = :bold, padding = (0,30,-10,15),
                halign = :left,valign =:bottom)

        Label(f_csd_ts[i,1:end,Top()], L"\textbf{\textrm{%$(sstring(version_names[i])): } \mathrm{𝗇_{𝖵} = %$(sstring(n_v)), } \mathrm{𝗇_{\hat{\alpha}_𝟣} = %$(sstring(n_a)), } \mathrm{𝗇_{𝖻𝗈𝗍𝗁} = %$(sstring(n_both))}}" ,
                padding = (0.0,0.0,20.0,20.0), valign = :center, font = :bold, fontsize = 20)
    end
    if legend
        elems = [
            #[LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(0.5, 0), (0.5, 1)])],
            [LineElement(color= (:black, 0.4), points = Point2f[(0, 0.5), (0.5, 0.5)]), LineElement(color= (:black, 0.9), points = Point2f[(0.5, 0.5), (1, 0.5)])],
            #[LineElement(color= (:black, 0.4), points = Point2f[(0, 0), (1, 0)]), LineElement(color= (:black, 0.9), points = Point2f[(0, 1), (1, 1)])],
            [PolyElement(color = (:darkred, 0.1), strokewidth = 0), LineElement(color= (:red,0.4), linewidth=2.0)],# LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:steelblue,0.1), strokewidth = 0), LineElement(color= (:blue, 0.4), linewidth=2.0)],#,LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:darkred, 0.8), strokewidth = 0), LineElement(color= (:red,0.9), linewidth=3.5)],#,LineElement(color = (:darkred,0.9), linewidth= 3.5,points = Point2f[(1, 0), (1, 1)])]]
            ]
        labels = [#"GI onset", 
                "EWS indicator", "increasing", "decreasing", "significantly increasing (𝗉<0.05)"]
        Legend(f_csd_ts[end+1,1:end], elems, labels, orientation = :horizontal, valign = :bottom, margin = (10,10,10,10)) 
    end
    colgap!(f_csd_ts.layout, 20)
    rowgap!(f_csd_ts.layout, 0)

    
    Label(f_csd_ts[1,1:end,Top()], "EWS in 100-year high-pass filtered NGRIP record", padding = (0.0,0.0,110.0,30.0), valign = :bottom, font = :bold, fontsize = 20)
    if saving
        save(saveto, f_csd_ts)
    end
    if showing
        display(f_csd_ts)
    end
end

# for lp in [true, false]
#     plot_wavepal_compare_csd_ts_label(;lowpass = lp, legend = true, plim = plim, showing = showing, saving = saving, saveto = "paper/figures/new/wavepal_compare_csd_ts_label_ngrip_lowpass_$(lp)_p_$(plim).pdf")
# end



# wavelet incicators (10,50) for "good" NGRIP 5y AND WAVEPAL with indicator TS plotted (using Morlet)
#Fig 8 and S20
function plot_wavepal_compare_wave_ts_label(;lowpass=lowpass, plim = plim, legend = false, showing = true, saving = false, saveto = "paper/wavepal_compare_wave_ts_label_10_50_ngrip_lowpass_$(lowpass)_p_$(plim).pdf")
    f_csd_ts = Figure(size=(1200,600))
    
    type = "NGRIP5"
    letters_h = reshape(letters[1:18],2,9)

    whichletter = 0

    if lowpass
        sca_path = "new_surrogate_files/NGRIP5/sca/s1_10_s2_50_C_lowpass_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"
        hurst_path = "new_surrogate_files/NGRIP5/hurst/s1_10_s2_50_C_lowpass_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"
    else
        sca_path = "new_surrogate_files/NGRIP5/sca/s1_10_s2_50_C_no_lowpass_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"
        hurst_path = "new_surrogate_files/NGRIP5/hurst/s1_10_s2_50_C_no_lowpass_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"
    end
    scas = load(sca_path)["slopes"]
    hursts = load(hurst_path)["slopes"]
    version_names = ["Regular 5-year resolution","Irregular temporal resolution"]
    for (i,(v,a)) in enumerate(zip([scas,wavepal_indicators[3]], [hursts, wavepal_indicators[4]]))
        glv = f_csd_ts[i,1] = GridLayout()
        gla = f_csd_ts[i,2] = GridLayout()
        axv = Axis(glv[1,1], 
                ylabel = L"\mathrm{\hat{𝗐}^𝟤}", 
                ylabelsize = 18,
                ylabelfont = :bold,
                xticks = (10_000:5_000:60_000, string.(10:5:60)), 
                )
        axa = Axis(gla[1,1], 
                ylabel = L"\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}}", 
                ylabelsize = 18,
                ylabelfont = :bold,
                xticks = (10_000:5_000:60_000, string.(10:5:60)), 
                )
        axv.xreversed = true
        axa.xreversed = true

        axv_event = Axis(glv[1,1],xaxisposition = :top, xticklabelrotation = pi*0.4)
        axv_event.xticks = (GI_onsets, event_labels)
        axv_event.xreversed = true
        hidespines!(axv_event)
        hideydecorations!(axv_event)

        axa_event = Axis(gla[1,1],xaxisposition = :top, xticklabelrotation = pi*0.4)
        axa_event.xticks = (GI_onsets, event_labels)
        axa_event.xreversed = true
        hidespines!(axa_event)
        hideydecorations!(axa_event)

        whichletter +=1
        n_v = 0 
        n_a = 0
        n_both = 0
        which_v = []

        if whichletter == 2
            axv.xlabel = "Age (kyr b2k)"
            axa.xlabel = "Age (kyr b2k)"
        else
            hidexdecorations!(axv, grid=false)
            hidespines!(axv, :b) 
            hidexdecorations!(axa, grid=false)
            hidespines!(axa, :b) 
        end

        if i>1
            hidespines!(axv,:t)
            hidespines!(axa,:t)
            hidexdecorations!(axv_event, grid=false)
            hidexdecorations!(axa_event, grid=false)
        end

        for ev in 1:17
            va = 0.4
            vl = 2.0
            vsc = (:white,1.0)
            if typeof(v.slopes[ev]) != Nothing
                if v.slopes[ev] >0 
                    vc = :red
                    vsc = (:darkred,0.1)
                    if v.p_one[ev] < plim
                        va = 0.9
                        vl = 3.5
                        vsc = (:darkred,0.8)
                        n_v +=1
                        push!(which_v, ev)
                    end       
                else
                    vc = :blue
                    vsc = (:steelblue,0.1)
                end
            end
            vspan!(axv, [GI_onsets[ev]], [GS_onsets[ev]], color = vsc)
            pv = Polynomials.fit(v.times[ev][findall(!ismissing, v.vals[ev])], convert.(Float64,v.vals[ev][findall(!ismissing, v.vals[ev])]),1)
            lines!(axv,v.times[ev], v.vals[ev], color = :black, alpha = va)#, linewidth = vl)
            lines!(axv, v.times[ev], pv.(v.times[ev]), color = vc, alpha = va, linewidth = vl)
            vlines!(axv, [GI_onsets[ev]], color = :darkred, alpha = va, linewidth = vl)

            aa = 0.4
            al = 2.0
            asc = (:white,1.0)
            if typeof(a.slopes[ev]) != Nothing
                if a.slopes[ev] >0 
                    ac = :red
                    asc = (:darkred,0.1)
                    if a.p_one[ev] < plim
                        aa = 0.9
                        al = 3.5
                        asc = (:darkred,0.8)
                        n_a +=1
                        if ev in which_v
                            n_both +=1
                        end
                    end
                
                else
                    ac = :blue
                    asc = (:steelblue,0.1)
                end
            end
            vspan!(axa, [GI_onsets[ev]], [GS_onsets[ev]], color = asc)
            pa = Polynomials.fit(a.times[ev][findall(!ismissing, a.vals[ev])], convert.(Float64,a.vals[ev][findall(!ismissing, a.vals[ev])]),1)
            lines!(axa,a.times[ev], a.vals[ev], color = :black, alpha = aa)#, linewidth = vl)
            lines!(axa, a.times[ev], pa.(a.times[ev]), color = ac, alpha = aa, linewidth = al)
            vlines!(axa, [GI_onsets[ev]], color = :darkred, alpha = aa, linewidth = al)

            xlims!(axv,nothing,nothing)
            axv.xreversed = true
            xlims!(axv_event, axv.xaxis.attributes.limits[]...)
            axv_event.xreversed = true
            xlims!(axa_event, axv.xaxis.attributes.limits[]...)
            axa_event.xreversed = true
        end
        Label(glv[1,1,TopLeft()], letters_h[1,whichletter], fontsize = 20,
                font = :bold, padding = (0,30,-10,15),
                halign = :left,valign =:bottom)
        Label(gla[1,1,TopLeft()], letters_h[2,whichletter], fontsize = 20,
                font = :bold, padding = (0,30,-10,15),
                halign = :left,valign =:bottom)
        Label(f_csd_ts[i,1:end,Top()], L"\textbf{\textrm{%$(sstring(version_names[i])): } \mathrm{𝗇_{\hat{𝗐}^𝟤} = %$(sstring(n_v)), } \mathrm{𝗇_{\hat{𝖧}^\text{𝗅𝗈𝖼}} = %$(sstring(n_a)), } \mathrm{𝗇_{𝖻𝗈𝗍𝗁} = %$(sstring(n_both))}}" ,
                padding = (0.0,0.0,20.0,20.0), valign = :center, font = :bold, fontsize = 20)
    end
    if legend
        elems = [
            #[LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(0.5, 0), (0.5, 1)])],
            [LineElement(color= (:black, 0.4), points = Point2f[(0, 0.5), (0.5, 0.5)]), LineElement(color= (:black, 0.9), points = Point2f[(0.5, 0.5), (1, 0.5)])],
            #[LineElement(color= (:black, 0.4), points = Point2f[(0, 0), (1, 0)]), LineElement(color= (:black, 0.9), points = Point2f[(0, 1), (1, 1)])],
            [PolyElement(color = (:darkred, 0.1), strokewidth = 0), LineElement(color= (:red,0.4), linewidth=2.0)],# LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:steelblue,0.1), strokewidth = 0), LineElement(color= (:blue, 0.4), linewidth=2.0)],#,LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:darkred, 0.8), strokewidth = 0), LineElement(color= (:red,0.9), linewidth=3.5)],#,LineElement(color = (:darkred,0.9), linewidth= 3.5,points = Point2f[(1, 0), (1, 1)])]]
            ]
        labels = [#"GI onset", 
                "EWS indicator", "increasing", "decreasing", "significantly increasing (𝗉<0.05)"]
        Legend(f_csd_ts[end+1,1:end], elems, labels, orientation = :horizontal, valign = :bottom, margin = (10,10,10,10)) 
    end

    colgap!(f_csd_ts.layout, 20)
    rowgap!(f_csd_ts.layout, 0)

    
    Label(f_csd_ts[1,1:end,Top()], "EWS in (10-50) year band of NGRIP record", padding = (0.0,0.0,110.0,30.0), valign = :bottom, font = :bold, fontsize = 20)
    if saving
        save(saveto, f_csd_ts)
    end
    if showing
        display(f_csd_ts)
    end
end

plot_wavepal_compare_wave_ts_label(lowpass=true, plim = plim, legend = true, showing = showing, saving = saving, saveto = "do_ews_across_greenland_ice_cores/figures/fig8.pdf")
plot_wavepal_compare_wave_ts_label(lowpass=false, plim = plim, legend = true, showing = showing, saving = saving, saveto = "do_ews_across_greenland_ice_cores/figures/figS20.pdf")


function plot_wavepal_compare_wave_ts_label_aggr(;lowpass=lowpass, plim = plim, legend = false, showing = true, saving = false, saveto = "paper/wavepal_compare_wave_ts_label_aggr_10_50_ngrip_lowpass_$(lowpass)_p_$(plim).pdf")
    

    f_csd_ts = Figure(size=(1200, 500))
    type = "NGRIP5"
    letters_h = reshape(letters[1:18],2,9)

    whichletter = 0

    if lowpass
        sca_path = "new_surrogate_files/NGRIP5/sca/s1_10_s2_50_C_lowpass_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"
        hurst_path = "new_surrogate_files/NGRIP5/hurst/s1_10_s2_50_C_lowpass_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"
    else
        sca_path = "new_surrogate_files/NGRIP5/sca/s1_10_s2_50_C_no_lowpass_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"
        hurst_path = "new_surrogate_files/NGRIP5/hurst/s1_10_s2_50_C_no_lowpass_MORLET_orig_gs_FiltInd_false_normalise_false_nocoi_true_onlyfull_true_smoothw_false_10000_TFTS.jld2"
    end
    scas = load(sca_path)["slopes"]
    hursts = load(hurst_path)["slopes"]
    version_names = ["Regular 5-year resolution","Irregular temporal resolution"]
    for (i,(v,a)) in enumerate(zip([scas,wavepal_indicators[3]], [hursts, wavepal_indicators[4]]))
        glv = f_csd_ts[i,1] = GridLayout()
        gla = f_csd_ts[i,2] = GridLayout()
        
        
        axv = Axis(glv[1,1],
            xticks = (1:17, ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]) ,
            xreversed = true,
            xticklabelrotation = pi*0.4,
            ylabel =  L"\mathrm{\hat{𝗐}^𝟤}", 
            #xlabel = "Transition",
            xminorgridwidth = 1,
            xminorgridcolor = :grey10 ,
            xminorgridvisible = true, 
            xgridvisible = false,
            yminorgridvisible = false, 
            ygridvisible = false,
            yticksvisible=false,
            yticklabelsvisible=false,
            ylabelsize = 18,
            ylabelfont = :bold,
            )
        axa = Axis(gla[1,1],
            xticks = (1:17, ["YD/PB","DO-1","DO-2","DO-3","DO-4","DO-5","DO-6","DO-7","DO-8","DO-9","DO-10","DO-11","DO-12","DO-13","DO-14","DO-15","DO-16"]) ,
            xreversed = true,
            xticklabelrotation = pi*0.4,
            ylabel =  L"\mathrm{\hat{𝖧}^\text{𝗅𝗈𝖼}}", 
            #xlabel = "Transition",
            xminorgridwidth = 1,
            xminorgridcolor = :grey10 ,
            xminorgridvisible = true, 
            xgridvisible = false,
            yminorgridvisible = false, 
            ygridvisible = false,
            yticksvisible=false,
            yticklabelsvisible=false,
            ylabelsize = 18,
            ylabelfont = :bold,
            )

        whichletter +=1
        n_v = 0 
        n_a = 0
        n_both = 0
        which_v = []

        if whichletter == 2
            axv.xlabel = "Event"
            axa.xlabel = "Event"
        else
            hidexdecorations!(axv, grid=false, minorgrid=false)
            hidespines!(axv, :b) 
            hidexdecorations!(axa, grid=false, minorgrid=false)
            hidespines!(axa, :b) 
        end

        if i>1
            hidespines!(axv,:t)
            hidespines!(axa,:t)
        end


        mv = zeros(17,1)
        ma = zeros(17,1)
        for ev in 1:17
            #va = 0.4
            #vl = 2.0
            #vsc = (:white,1.0)
            #vc = :black
            if typeof(v.slopes[ev]) != Nothing
                if v.slopes[ev] >0 
                    mv[ev,1]+=1
                    #pone = count(v -> v ≥ true_slopes[ev], surr_slopes[ev,:])/length(surr_slopes[ev,:])
            #        vc = :red
            #        vsc = (:darkred,0.1)
                    if v.p_one[ev] < plim
                        mv[ev,1]+=1
            #            va = 0.9
            #            vl = 3.5
            #            vsc = (:darkred,0.8)
                        
                        n_v +=1
                        push!(which_v, ev)
                        
                    end
                end
            end

            if typeof(a.slopes[ev]) != Nothing
                if a.slopes[ev] >0 
                    ma[ev,1]+=1
                    #pone = count(v -> v ≥ true_slopes[ev], surr_slopes[ev,:])/length(surr_slopes[ev,:])
            #        vc = :red
            #        vsc = (:darkred,0.1)
                    if a.p_one[ev] < plim
                        ma[ev,1]+=1
            #            va = 0.9
            #            vl = 3.5
            #            vsc = (:darkred,0.8)
                        n_a +=1
                        if ev in which_v
                            n_both +=1
                        end
                    end
                end
            end


        end
        co2 = cgrad([:steelblue, :darkred, :darkred], 3, categorical = true, alpha=[0.1,0.1,0.8])
        hmv = heatmap!(axv, mv, colormap=co2, colorrange=(0,2))
        hma = heatmap!(axa, ma, colormap=co2, colorrange=(0,2))
        # xlims!(axv,nothing,nothing)
        # axv.xreversed = true
        # xlims!(ax_event, axv.xaxis.attributes.limits[]...)
    


        Label(glv[1,1,TopLeft()], letters_h[1,whichletter], fontsize = 20,
                    font = :bold, padding = (0,20,-10,15),
                    halign = :left, valign =:bottom)
        Label(gla[1,1,TopLeft()], letters_h[2,whichletter], fontsize = 20,
                    font = :bold, padding = (0,20,-10,15),
                    halign = :left, valign =:bottom)
        Label(f_csd_ts[i,1:end,Top()], L"\textbf{\textrm{%$(sstring(version_names[i])): } \mathrm{𝗇_{\hat{𝗐}^𝟤} = %$(sstring(n_v)), } \mathrm{𝗇_{\hat{𝖧}^\text{𝗅𝗈𝖼}} = %$(sstring(n_a)), } \mathrm{𝗇_{𝖻𝗈𝗍𝗁} = %$(sstring(n_both))}}" ,
                    padding = (0.0,0.0,20.0,20.0), valign = :center, font = :bold, fontsize = 20)
    end
    

    
    
    if legend
        # elems = [
        #     [LineElement(color= (:black, 0.4), points = Point2f[(0, 0.5), (0.5, 0.5)]), LineElement(color= (:black, 0.9), points = Point2f[(0.5, 0.5), (1, 0.5)])],
        #     [PolyElement(color = (:darkred, 0.1), strokewidth = 0), LineElement(color= (:red,0.4), linewidth=2.0)],# LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
        #     [PolyElement(color = (:steelblue,0.1), strokewidth = 0), LineElement(color= (:blue, 0.4), linewidth=2.0)],#,LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
        #     [PolyElement(color = (:darkred, 0.8), strokewidth = 0), LineElement(color= (:red,0.9), linewidth=3.5)],#,LineElement(color = (:darkred,0.9), linewidth= 3.5,points = Point2f[(1, 0), (1, 1)])]]
        #     ]
        # labels = ["EWS indicator", "increasing", "decreasing", "significantly increasing (𝗉<0.05)"]
        elems = [
            [PolyElement(color = (:steelblue,0.1), strokewidth = 0)],#,LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:darkred, 0.1), strokewidth = 0)],# LineElement(color = (:darkred,0.4), linewidth= 2.0,points = Point2f[(1, 0), (1, 1)])],
            [PolyElement(color = (:darkred, 0.8), strokewidth = 0)],#,LineElement(color = (:darkred,0.9), linewidth= 3.5,points = Point2f[(1, 0), (1, 1)])]]
            ]
        labels = ["increasing", "decreasing", "significantly increasing (𝗉<0.05)"]
        Legend(f_csd_ts[end+1,1:end], elems, labels, orientation = :horizontal, valign = :bottom, margin = (10,10,10,10)) 
    end

    colgap!(f_csd_ts.layout, 20)
    rowgap!(f_csd_ts.layout, 0)

    
    Label(f_csd_ts[1,1:end,Top()], "EWS in (10-50) year band of NGRIP record", padding = (0.0,0.0,80.0,3.0), valign = :bottom, font = :bold, fontsize = 20)

    if saving
        save(saveto,f_csd_ts)
    end
    if showing
        display(f_csd_ts)
    end
end

plot_wavepal_compare_wave_ts_label_aggr(lowpass=true, plim = plim, legend = saving, showing = showing, saving = saving, saveto = "do_ews_across_greenland_ice_cores/figures/fig8_aggr.pdf")
plot_wavepal_compare_wave_ts_label_aggr(lowpass=false, plim = plim, legend = true, showing = showing, saving = saving, saveto = "do_ews_across_greenland_ice_cores/figures/figS20_aggr.pdf")

#Fig A1
function plot_sampling_steps_ngrip(;showing=true, saving=false, saveto = "paper/raw_NGRIP_sampling.pdf")
    NGRIP_data = XLSX.readxlsx("ice_core_data/NGRIP_d18O_and_dust_5cm.xlsx")
    NGRIP2_data = NGRIP_data["NGRIP-2 d18O and Dust"]
    NGRIP2_age = NGRIP2_data[:][2:end,4]
    NGRIP2_age_short = NGRIP2_age[10295 .<= NGRIP2_age .<= 59920];
    diff_NGRIP2_short = diff(NGRIP2_age_short)  

    fngrip = Figure(size = (1000,500))
    g1 = fngrip[1,1] = GridLayout()
    g2 = fngrip[1,2] = GridLayout()
    ax1 = Axis(g1[1,1],
           xlabel = "Age (kyr b2k)",
            ylabel = "Sampling steps (year)",
            yticks = 0:7,
            xticks = (10_000:10_000:60_000, string.(10:10:60)),)
    ax1.xreversed = true
    
    lines!(ax1, convert.(Float64, NGRIP2_age_short[1:end-1]),diff_NGRIP2_short, label = "Sampling steps", color = :black)
    Label(g1[1,1,TopLeft()], "(a)", fontsize = 20,
                        font = :bold, padding = (0,5,5,0),
                        halign = :left)
    
    ax2 = Axis(g2[1,1], xlabel = "Relative freuquency",
                    yticks = 0:7,
                    ylabel = "nothing",
                    ylabelcolor = :transparent
                    )
    CairoMakie.density!(ax2,diff_NGRIP2_short, direction =:y,
                    color = (:black,0.6), 
                    strokecolor = :black, 
                    strokewidth = 2, 
                    boundary = extrema(diff_NGRIP2_short),
                    label = "Distribution"
                    )
    hlines!(ax2, [5], color = :red)
    Label(g2[1,1,TopLeft()], "(b)", fontsize = 20,
                        font = :bold, padding = (0,5,5,0),
                        halign = :left)
    axislegend(ax1)
    axislegend(ax2)
    CairoMakie.ylims!(ax1, (0,7))
    CairoMakie.ylims!(ax2, (0,7))
    linkyaxes!(ax1,ax2)

    if saving
        save(saveto, fngrip)
    end
    if showing
        display(fngrip)
    end
end

plot_sampling_steps_ngrip(showing=showing, saving=saving, saveto = "figures/figA1.pdf")


#Fig A2
function plot_sampling_steps_neem(;showing=true, saving=false, saveto = "paper/raw_NEEM_sampling.pdf")
    neem_data = readdlm("ice_core_data/NEEM_d18O.tab",'\t', skipstart = 46)
    neem_data[neem_data .== ""] .= missing
    NEEM_age_GICC05 = neem_data[:,2]

    NEEM_age = NEEM_age_GICC05 .*1_000 # in years instead of ky
    NEEM_age_nomiss = NEEM_age[1:findfirst(ismissing, NEEM_age)-1]
    convert.(Float64,NEEM_age_nomiss)
    NEEM_age_short = NEEM_age_nomiss[10295 .<= NEEM_age_nomiss .<= 59920]
    diff_NEEM_short = diff(NEEM_age_short);

    fneem = Figure(size = (1000,500))
    g1 = fneem[1,1] = GridLayout()
    g2 = fneem[1,2] = GridLayout()
    ax1 = Axis(g1[1,1],
            xlabel = "Age (kyr b2k)",
            ylabel = "Sampling steps (year)",
            yticks = 0:12,
            xticks = (10_000:10_000:60_000, string.(10:10:60)),
            )
    ax1.xreversed = true

    lines!(ax1, NEEM_age_short[1:end-1],diff_NEEM_short, label = "Sampling steps", color = :black)
    Label(g1[1,1,TopLeft()], "(a)", fontsize = 20,
                        font = :bold, padding = (0,5,5,0),
                        halign = :left)

    ax2 = Axis(g2[1,1], xlabel = "Relative frequency",
            yticks = 0:12,
            ylabel = "nothing",
            ylabelcolor = :transparent)
    CairoMakie.density!(ax2,diff_NEEM_short, direction =:y,
                    color = (:black,0.6), 
                    strokecolor = :black, 
                    strokewidth = 2,
                    boundary = extrema(diff_NEEM_short),
                    label = "Distribution"
                    )
    hlines!(ax2, [10], color = :red)
    Label(g2[1,1,TopLeft()], "(b)", fontsize = 20,
                        font = :bold, padding = (0,5,5,0),
                        halign = :left)
    axislegend(ax1)
    axislegend(ax2)
    linkyaxes!(ax1,ax2)

    if saving
        save(saveto, fneem)
    end
    if showing
        display(fneem)
    end
end

plot_sampling_steps_neem(showing=showing, saving=saving, saveto = "figures/figA2.pdf")



#Fig S1 and S2
function plot_surrogate_acf_and_var(;gs=2, n_surrs=10_000, showing=true, saving = false, saveto = "paper/surrogate_acf_var_NGRIP5_lowpass_true_gs_2_tfts.pdf")
    ngrip_core = load("new_surrogate_files/ice_cores/NGRIP5/C_lowpass.jld2")["ice"]
    x = ngrip_core.δ[ngrip_core.cold_idx[gs]]
    s = surrogate(x, TFTS(0.05))
    n = n_surrs
    lags = 0:1:Int(round(min(size(x,1)-1, 10*log10(size(x,1)))))
    
    f= Figure(size = (1000,400))
    g1 = f[1,1] = GridLayout()
    g3 = f[1,2] = GridLayout()
    ax1 = Axis(g1[1,1],
                xlabel = "lag",
                ylabel = "ACF",
                xticks = 0:2:100,
                yticks=-1.2:0.2:1.2
                )

    sc_s=CairoMakie.scatter!(ax1, [1],[var(s)], color = :cadetblue3, label = "surrogate", alpha = 0.1, visible=false)
    var_s =Float64[]
    for i in 1:n
            s= surrogate(x, TFTS(0.05))
            lines!(ax1, lags, autocor(s,lags), color = :cadetblue3, alpha = 0.2, linewidth = 2, label = "surrogate")
            push!(var_s, var(s))
    end

    ax3 = Axis(g3[2,1], 
                ylabel = "Variance",
                xticks = 1:1:1,
                )
    CairoMakie.xlims!(ax3, (0.9,1.1))
    CairoMakie.ylims!(ax3, extrema(var_s).+(-1e-4,1e-4))

    rb= rangebars!(ax3, [1], [minimum(var_s)], [maximum(var_s)],
        label = "TFTS",
        color = :cadetblue3,
        linecap=true,
        whiskerwidth = 30)
    sc = CairoMakie.scatter!(ax3,[1],[var(x)],  label = "data",color = :black)

    lines!(ax1, lags,autocor(x,lags), label = "data", color = :black, linewidth = 2)

    hidexdecorations!(ax3)


    Label(g1[1,1,TopLeft()], "(a)", fontsize = 20,
                        font = :bold, padding = (0,5,5,0),
                        halign = :left)
    Label(g3[1,1,TopLeft()], "(b)", fontsize = 20,
                        font = :bold, padding = (0,5,5,0),
                        halign = :left)

    colsize!(f.layout,1,Relative(18/20))
    colsize!(f.layout,2,Relative(2/20))
    colgap!(f.layout,40)

    rowsize!(g3,1, Relative(1/5))
    rowsize!(g3,2, Relative(4/5))
    rowgap!(g3,4)

    axislegend(ax1, merge=true,unique=true)

    Legend(g3[1,1], [[
            LineElement(linepoints = [Point2f(0.5, 0), Point2f(0.5, 1)], color = :cadetblue3),
            MarkerElement(points = [Point2f(0.5, 0), Point2f(0.5, 1)], marker = :hline, markersize = 10, color = :cadetblue3), 
            sc_s],
        sc],
        ["surrogate","data"], merge=true, unique=true, tellheight=true)

    if saving
        save(saveto, f)
    end
    if showing
        display(f)
    end
end

plot_surrogate_acf_and_var(gs=2, n_surrs=10_000, showing=showing, saving = saving, saveto = "figures/figS1.pdf")
plot_surrogate_acf_and_var(gs=17, n_surrs=10_000, showing=showing, saving = saving, saveto = "figures/figS2.pdf")


#Fig 1
function plot_map(;showing=true, saving = false, saveto = "paper/map2.pdf")
    fig = Figure(size = (760, 1000))
    ga = GeoAxis(fig[1, 1]; 
            limits = ((-60, 20), (56, 77)), 
            source="+proj=latlong", 
            dest = "+proj=stere +lat_0=90 +lon_0=-45", 
            xticks = -105:15:15,
            yticks = 55:5:85,
            tellwidth=true
            )
    
    # plot coastlines from Natural Earth as a reference
    li = lines!(ga, naturalearth("ne_10m_coastline"), color = "black")


    scatter!(ga,[-42.32],[75.10], color=:red, markersize=10) #NGRIP
    scatter!(ga,[-37.64],[72.58], color=:red, markersize=10) #GRIP
    scatter!(ga,[-38.48],[72.58], color=:red, markersize=10) #GISP2
    scatter!(ga,[-51.06],[77.45], color=:red, markersize=10) #NEEM

    text!(ga,[-42.32], [75.10], text=" NGRIP") #NGRIP
    text!(ga,[-37.64],[72.58], text=" GRIP") #GRIP
    text!(ga,[-38.48],[72.58], text="GISP2 ", align=(:right,:top)) #GISP2
    text!(ga,[-51.06],[77.45], text=" NEEM") #NEEM

    ga.elements[:yticklabels].text = ["", "",  "60° W", "45° W", "30° W", "", "", ""]
    ga.elements[:xticklabels].text = ["60° N  ", "65° N  ", "70° N  ", "75° N  ", "80° N               "]
    if saving
        save(saveto, fig)
    end
    if showing
        display(fig)
    end
end


plot_map(showing=showing, saving = saving, saveto = "figures/fig1.pdf")







