using Pkg
Pkg.activate(".")

using XLSX
using Interpolations 
using DSP #for filtering
using TimeseriesSurrogates
using Random
using NPZ


include("base_functions.jl")

#Times of GS and GIs as used in Boers, 2018
global GS_onsets_niklas = [12_900,23_105,27_460,28_510,32_025,33_390,34_725,36_590,39_935,40_790,42_100,44_285,49_105,51_650,54_745,56_440,58_515]
global GI_onsets_niklas = [11_705,14_690,23_375,27_790,28_910,32_520,33_735,35_505,38_235,40_165,41_480,43_365,46_860,49_315,54_235,55_815,58_280]

#Alternatively:
# gi_onsets = XLSX.readxlsx("ice_core_data/GI_onsets.xlsx")
# XLSX.sheetnames(gi_onsets)
# GI_data = gi_onsets["GI_onsets"]
# GI_onsets_luc = GI_data[:][2:end,1]

# gs_onsets = XLSX.readxlsx("ice_core_data/GS_onsets.xlsx")
# XLSX.sheetnames(gs_onsets)
# GS_data = gs_onsets["GS_onsets"]
# GS_onsets_luc = GS_data[:][2:end,1]


####################### raw NGRIP data #########################

# what Niklas uses (raw ages are rounded in this dataset):
NGRIP_data_n0 = XLSX.readxlsx("ice_core_data/ngrip_age_depth_dust_d18o_unc_10k_60k.xlsx")
NGRIP_data_n = NGRIP_data_n0["ngrip_age_depth_dust_d18o_unc_1"]
age_n = NGRIP_data_n[:][:,2][1:2:end]
δ_n = NGRIP_data_n[:][:,3][1:2:end]
dust_n = NGRIP_data_n[:][:,4][1:2:end]
convert.(Float64, age_n)
end_age = 10277
start = argmin(abs.(age_n .- end_age))
age_int_n = Int.(floor.(age_n))

age_n = age_n[start:end-1]
age_int_n = age_int_n[start:end-1]
age_ipd_n = collect(minimum(age_int_n):maximum(age_int_n))
δ_n = δ_n[start:end-1]
dust_n = dust_n[start:end-1]

samp = 5
cutoff = 0.5 *1 /samp
fs = 1.0
rp = 0.05
order = 8
sec_factor = 0.95 
filtb = 100
offset = 19 #offset at beginning and end
w_order = 4

f_δ_n = Spline1D(age_n, δ_n, k=3, s= 0.0, bc = "extrapolate")
d18o_ipd_cub_temp = f_δ_n(age_ipd_n)
d18o_ipd_cub = cheby_lowpass_filter(d18o_ipd_cub_temp,sec_factor * cutoff, fs, order, rp)

dat_n_low = d18o_ipd_cub[offset+1 : end-offset][1:samp:end][end:-1:1]
time_n_low = age_ipd_n[offset+1 : end-offset][1:samp:end][end:-1:1]

dat_n = d18o_ipd_cub_temp[offset+1 : end-offset][1:samp:end][end:-1:1]
time_n = time_n_low

mdat_n = mean(dat_n)
sdat_n = std(dat_n)
dat_n_norm = (dat_n .- mdat_n) / sdat_n

mdat_n_low = mean(dat_n_low)
sdat_n_low = std(dat_n_low)
dat_n_low_norm = (dat_n_low .- mdat_n_low) / sdat_n_low


##########################
#NGRIP data from the KU PICE homepage (newer version)
NGRIP_data = XLSX.readxlsx("ice_core_data/NGRIP_d18O_and_dust_5cm.xlsx")
NGRIP2_data = NGRIP_data["NGRIP-2 d18O and Dust"]
#NGRIP2_depth = NGRIP2_data[:][2:end,1]
NGRIP2_δ = NGRIP2_data[:][2:end,2]
#NGRIP2_dust = NGRIP2_data[:][2:end,3]
NGRIP2_age = NGRIP2_data[:][2:end,4]
#NGRIP2_mce = NGRIP2_data[:][2:end,5]

NGRIP2_age = convert.(Float64, NGRIP2_age)
NGRIP2_δ = convert.(Float64, NGRIP2_δ)

end_age = 10277
start = argmin(abs.(NGRIP2_age .- end_age))
stopp = argmin(abs.(NGRIP2_age .- age_n[end]))
age_c = NGRIP2_age[start:stopp]
age_int_c = Int.(floor.(age_c))
δ_c = NGRIP2_δ[start:stopp]
age_ipd_c = collect(minimum(age_int_c):maximum(age_int_c))

spl_NGRIP = Spline1D(age_c, δ_c, k=3, s= 0.0)
δ_NGRIP_5_temp = spl_NGRIP(age_ipd_c)
δ_NGRIP_5 = cheby_lowpass_filter(δ_NGRIP_5_temp,sec_factor * cutoff, fs, order, rp)

dat_c_low = δ_NGRIP_5[offset+1 : end-offset][1:samp:end][end:-1:1]
time_c_low = age_ipd_c[offset+1 : end-offset][1:samp:end][end:-1:1]

dat_c = δ_NGRIP_5_temp[offset+1 : end-offset][1:samp:end][end:-1:1]
time_c = time_c_low

mdat_c = mean(dat_c)
sdat_c = std(dat_c)
dat_c_norm = (dat_c .- mdat_c) / sdat_c

mdat_c_low = mean(dat_c_low)
sdat_c_low = std(dat_c_low)
dat_c_low_norm = (dat_c_low .- mdat_c_low) / sdat_c_low

####################
# Filtering of the different NGRIP version
offset = div(200,samp)

times = [time_c_low, time_n_low, time_c, time_n]
dats = [dat_c_low, dat_n_low, dat_c, dat_n]
norm_dats = [dat_c_low_norm, dat_n_low_norm, dat_c_norm, dat_n_norm]
labels = ["C (not rounded), lowpass",
            "N (rounded), lowpass",
            "C (not rounded), no lowpass",
            "N (rounded), no lowpass"]
names = ["C_lowpass", "N_lowpass", "C_no_lowpass", "N_no_lowpass"]

filt_dats_100 = []
for dat in dats
    filt_dat = cheby_highpass_filter(dat, .95 * 1. / filtb, 1. / samp, 8, .05)
    push!(filt_dats_100, filt_dat)
end

normed_filt_dats_100 = []
for dat in norm_dats
    filt_dat = cheby_highpass_filter(dat, .95 * 1. / filtb, 1. / samp, 8, .05)
    push!(normed_filt_dats_100, filt_dat)
end

#data structure for ice core data
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
    cold_idx_n::Array #N searches for EWS in this interval (and not the entire one!!)
end


ice_cores = ice_core[]
for (i,d) in enumerate(dats)
    ice = ice_core(names[i], times[i], 
        d, norm_dats[i], filt_dats_100[i],normed_filt_dats_100[i],
        5, [],[],[])
    push!(ice_cores, ice)
end

#saving indices for the GS (cold periods)
for ice in ice_cores
    cold_indices = []
    for (k,event) in enumerate(GI_onsets_niklas)
        cold_idx = @. event <= ice.age <=GS_onsets_niklas[k]
        push!(cold_indices, cold_idx)
    end
    ice.cold_idx =cold_indices
end

#saving indices for the GI (warm periods)
for ice in ice_cores 
    warm_idx = []
    end_idx = length(ice.age)
    first = 1:findfirst(ice.cold_idx[end])
    push!(warm_idx, first)
    for ic in length(ice.cold_idx)-1:-1:1
        warm = findlast(ice.cold_idx[ic+1]):findfirst(ice.cold_idx[ic])
        push!(warm_idx, warm)
    end
    last = findlast(ice.cold_idx[1]):end_idx
    push!(warm_idx, last)
    ice.warm_idx = warm_idx
end

#saving indices for the truncated GS (cold periods) as done by Niklas in Boers, 2018
for ice in ice_cores
    niklas_colds = []
    for (i,ci) in enumerate(ice.cold_idx)
        niklas_cold = falses(length(ci))
        niklas_cold[findall(ci)[1:end-offset-1]] .= true
        push!(niklas_colds, niklas_cold)
    # @show niklas_colds
    end
    ice.cold_idx_n = niklas_colds
end


global ice_cores #make list of the different NGRIP core versions a global variable for easier use

#data structure for indicators and their significance
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

gs_lengths = []
for ice in ice_cores
    gs_lengths_per_core = Int64[]
    for gs in 1:17
        #@show ice.name
        gs_length =  sum(ice.cold_idx[gs])
        push!(gs_lengths_per_core, gs_length)
    end
    push!(gs_lengths, gs_lengths_per_core)
end


gs_lengths_n = []
for ice in ice_cores
    gs_lengths_per_core = Int64[]
    for gs in 1:17
        #@show ice.name
        gs_length =  sum(ice.cold_idx_n[gs])
        push!(gs_lengths_per_core, gs_length)
    end
    push!(gs_lengths_n, gs_lengths_per_core)
end



#data structure for the "overall" significant test (how many EWS are expected to occur by chance)
mutable struct distribution_significant_increases
    name::String
    num_inc_one::Array #{Int64}
    num_inc_two::Array #{Int64}
    n::Int32
    n_GS::Int32
    pvals::Any #Array{Float64}
    type::String
    s1::Int32
    s2::Int32
    resolution::Int32
end

##################################
function save_ice_cores(ice_cores) #function to save NGRIP at 5y resolution as ice_core object
    for ice in ice_cores
        saveto = "new_surrogate_files/ice_cores/NGRIP5/$(ice.name).jld2"
        if !isfile(saveto)
            println("saving $saveto")
            save(saveto, "ice", ice)
        else
            println("$saveto already exists")
        end
    end
    println("Done :)")
end


#println("SAVE different version of NGRIP (5y) as ICE CORES")
save_ice_cores(ice_cores)

### NGRIP (5y) preprocessed in python 2 --> transform to ice core objects as well!
py_core_paths = readdir("new_surrogate_files/ice_cores_py/", join = true)
py_core_paths = [p for p in py_core_paths if occursin(".npz",p)]
py_core_names = readdir("new_surrogate_files/ice_cores_py/", join = false)
py_core_names = [p for p in py_core_names if occursin(".npz",p)]

##########
py2_c_file = npzread(py_core_paths[1])
py2_n_file = npzread(py_core_paths[2])

new_times = [py2_n_file["time"],py2_c_file["time"]]
new_dats = [py2_n_file["dat"], py2_c_file["dat"]]
new_dats_normed = [py2_n_file["dat_normed"],py2_c_file["dat_normed"]]
new_filt_dats= [py2_n_file["filt_dat"], py2_c_file["filt_dat"]]
new_names = ["N_lowpass_py2", "C_lowpass_py2"]

for (i,dat) in enumerate(new_dats)
    sampl = 5
    offset5 = div(200,sampl)
    cold_indices = []
    for (k,event) in enumerate(GI_onsets_niklas)
        cold_idx = @. event <= new_times[i] <=GS_onsets_niklas[k]
        push!(cold_indices, cold_idx)
    end

    warm_idx = []
    end_idx = length(new_times[i])
    first = 1:findfirst(cold_indices[end])
    push!(warm_idx, first)
    for ic in length(cold_indices)-1:-1:1
        warm = findlast(cold_indices[ic+1]):findfirst(cold_indices[ic])
        push!(warm_idx, warm)
    end
    last = findlast(cold_indices[1]):end_idx
    push!(warm_idx, last)

    niklas_colds = []
    for (j,ci) in enumerate(cold_indices)
        niklas_cold = falses(length(ci))
        niklas_cold[findall(ci)[1:end-offset5]] .= true
        push!(niklas_colds, niklas_cold)
    end
    ice = ice_core(new_names[i], new_times[i], new_dats[i], new_dats_normed[i], [9999.999], new_filt_dats[i], sampl, cold_indices,warm_idx,niklas_colds)
    saveto = "new_surrogate_files/ice_cores_py/$(ice.name).jld2"
    if !isfile(saveto)    
        save(saveto, "ice", ice)
    end
end



#### 10y ice core data ####
## NEEM ##
neem_data = readdlm("ice_core_data/NEEM_d18O.tab",'\t', skipstart = 46)#, Float64)
neem_data[neem_data .== ""] .= missing
# #NEEM_depth = neem_data[:,1]
NEEM_age_GICC05 = neem_data[:,2]
NEEM_age_GICC05_modelext = neem_data[:,3]
# #NEEM_age_GICC05_AICC2012 = neem_data[:,4]
# #NEEM_mce_GICC05 = neem_data[:,5]
NEEM_δ = neem_data[:,6]
# # NEEM_δD = neem_data[:,7]
# # NEEM_δ_std =  neem_data[:,8]
NEEM_age = NEEM_age_GICC05 .*1_000 # in years instead of kyrs

## NGRIP ##
NGRIP_data = XLSX.readxlsx("ice_core_data/NGRIP_d18O_and_dust_5cm.xlsx")
NGRIP2_data = NGRIP_data["NGRIP-2 d18O and Dust"]
#NGRIP2_depth = NGRIP2_data[:][2:end,1]
NGRIP2_δ = NGRIP2_data[:][2:end,2]
#NGRIP2_dust = NGRIP2_data[:][2:end,3]
NGRIP2_age = NGRIP2_data[:][2:end,4]


#### 20y ice core data ####
## NGRIP, GRIP and GISP2 ##
xf = XLSX.readxlsx("ice_core_data/GICC05modelext_GRIP_and_GISP2_and_resampled_data_series_Seierstad_et_al._2014_version_10Dec2014-2.xlsx")
d18O_20_year_ds = xf["3) d18O and Ca 20 yrs mean"]
age_20_orig = d18O_20_year_ds[:][53:end,1]

NGRIP2_20_δ = d18O_20_year_ds[:][53:end,5]
NGRIP2_20_δ[NGRIP2_20_δ.=="NaN"] .=missing

GRIP_20_δ = d18O_20_year_ds[:][53:end,8]
GRIP_20_δ[GRIP_20_δ.=="NaN"] .=missing

GISP2_20_δ = d18O_20_year_ds[:][53:end,11]
GISP2_20_δ[GISP2_20_δ.=="NaN"] .=missing

age_20_new =  convert(Array{Int64},age_20_orig)[1:2:end].+20
age_20 = age_20_new[10300 .<= age_20_new .<= 59920]

NGRIP2_20_δ_new = NGRIP2_20_δ[1:2:end][10300 .<= age_20_new .<= 59920]
GRIP_20_δ_new = GRIP_20_δ[1:2:end][10300 .<= age_20_new .<= 59920]
GISP2_20_δ_new = GISP2_20_δ[1:2:end][10300 .<= age_20_new .<= 59920]

δ_20s = [NGRIP2_20_δ_new[end:-1:1], GRIP_20_δ_new[end:-1:1],GISP2_20_δ_new[end:-1:1]]
time20 = age_20[end:-1:1]


time10 = collect(10290:10:59920)
δ_10s = []
for (t,d) in zip([NEEM_age,NGRIP2_age], [NEEM_δ, NGRIP2_δ])
    #end_age = 10277
    start10 = argmin(abs.(skipmissing(t .- end_age)))
    stopp10 = argmin(abs.(skipmissing(t .- age_n[end])))
    age = t[start10:stopp10]
    δ = d[start10:stopp10]
    age = age[findall(!ismissing, age)]
    δ = δ[findall(!ismissing, age)]
    spl = Spline1D(age, δ, k=3, s= 0.0)
    
    δ_10 = spl(time10)
    push!(δ_10s, δ_10[end:-1:1])

    time_10_yearly = time10[1]:1:time10[end]
    δ_10_yearly = spl(time_10_yearly)
    samp_10 = 10
    cutoff_10 = (0.5 *1)/samp_10
    lp_filtered_10  = cheby_lowpass_filter(δ_10_yearly,sec_factor * cutoff_10, fs, order, rp)
    δ_10_lp = lp_filtered_10[1:10:end]
    push!(δ_10s, δ_10_lp[end:-1:1])
end
time10 = time10[end:-1:1]


names_10 = ["NEEM", "NGRIP", "NEEM_lp", "NGRIP_lp"]
ice_cores_10 = ice_core[]

for (i,dat) in enumerate(δ_10s)
        #filtb = 100
        sampl = 10
        offset10 = div(200,sampl)
        filt_dat = cheby_highpass_filter(dat, .95 * 1. / filtb, 1. / sampl, 8, .05)
        dat_norm = (dat .- mean(dat))/std(dat)
        normed_filt_dat = cheby_highpass_filter(dat_norm, .95 * 1. / filtb, 1. / sampl, 8, .05)
        cold_indices = []
        for (k,event) in enumerate(GI_onsets_niklas)
            cold_idx = @. event <= time10 <=GS_onsets_niklas[k]
            push!(cold_indices, cold_idx)
        end

        warm_idx = []
        end_idx = length(time10)
        first = 1:findfirst(cold_indices[end])
        push!(warm_idx, first)
        for ic in length(cold_indices)-1:-1:1
            warm = findlast(cold_indices[ic+1]):findfirst(cold_indices[ic])
            push!(warm_idx, warm)
        end
        last = findlast(cold_indices[1]):end_idx
        push!(warm_idx, last)

        niklas_colds = []
        for (j,ci) in enumerate(cold_indices)
            niklas_cold = falses(length(ci))
            niklas_cold[findall(ci)[1:end-offset10]] .= true
            push!(niklas_colds, niklas_cold)
        end

        ice = ice_core(names_10[i], time10, dat, dat_norm, filt_dat, normed_filt_dat, sampl, cold_indices,warm_idx,niklas_colds)
        saveto10 = "new_surrogate_files/ice_cores/10y/$(ice.name).jld2"
        if !isfile(saveto10)    
            save(saveto10, "ice", ice)
        end
        push!(ice_cores_10, ice)
end

names_20 = ["NGRIP", "GRIP", "GISP2"]
ice_cores_20 = ice_core[]

missings_gisp2 = []

for (i,dat) in enumerate(δ_20s)
    #filtb = 100
    sampl = 20
    offset20 = div(200,sampl)

    cold_indices = []
    for (k,event) in enumerate(GI_onsets_niklas)
        cold_idx = @. event <= time20 <=GS_onsets_niklas[k]
        #@show time20[findfirst(cold_idx)], time20[findlast(cold_idx)]
        push!(cold_indices, cold_idx)
    end

    warm_idx = []
    end_idx = length(time20)
    first = 1:findfirst(cold_indices[end])
    push!(warm_idx, first)
    for ic in length(cold_indices)-1:-1:1
        warm = findlast(cold_indices[ic+1]):findfirst(cold_indices[ic])
        push!(warm_idx, warm)
    end
    last = findlast(cold_indices[1]):end_idx
    push!(warm_idx, last)

    niklas_colds = []
    for (j,ci) in enumerate(cold_indices)
        niklas_cold = falses(length(ci))
        niklas_cold[findall(ci)[1:end-offset20]] .= true
        push!(niklas_colds, niklas_cold)
    end
    #handle missing data in GISP2 by replacing with random value from the same period
    if i==3
        push!(missings_gisp2, findall(ismissing, dat)...)
        for mv in findall(ismissing, dat)
            id_cold = findfirst([cold[mv] for cold in cold_indices])
            if typeof(id_cold) == Nothing
                id= findfirst([mv in warm for warm in warm_idx])
                close_indices = intersect(mv-1:mv+1, warm_idx[id])
                non_mis = findall(!ismissing, dat[close_indices])
                dat[mv] = mean(dat[non_mis])#rand(dat[non_mis])
            else
                id = id_cold
                close_indices = intersect(mv-div(60,sampl):mv+div(60,sampl), findall(cold_indices[id]))
                non_mis = findall(!ismissing, dat[close_indices])

                ndist_loc = Normal(mean(dat[non_mis]), std(dat[non_mis])/2)
                dat[mv] = rand(ndist_loc)
            end
        end
        #@show findall(ismissing, dat)
    end

    filt_dat = cheby_highpass_filter(dat, .95 * 1. / filtb, 1. / sampl, 8, .05)
    dat_norm = (dat .- mean(dat))/std(dat)
    normed_filt_dat = cheby_highpass_filter(dat_norm, .95 * 1. / filtb, 1. / sampl, 8, .05)

    ice = ice_core(names_20[i], time20, dat, dat_norm, filt_dat, normed_filt_dat, sampl, cold_indices,warm_idx,niklas_colds)
    saveto20 = "new_surrogate_files/ice_cores/20y/$(ice.name).jld2"
    if !isfile(saveto20)    
        save(saveto20, "ice", ice)
    end
    push!(ice_cores_20, ice)
end

#@show missings_gisp2

gs_lengths_10 = []
for ice in ice_cores_10
    gs_lengths_per_core = Int64[]
    for gs in 1:17
        #@show ice.name
        gs_length =  sum(ice.cold_idx[gs])
        push!(gs_lengths_per_core, gs_length)
    end
    push!(gs_lengths_10, gs_lengths_per_core)
end


gs_lengths_n_10 = []
for ice in ice_cores_10
    gs_lengths_per_core = Int64[]
    for gs in 1:17
        #@show ice.name
        gs_length =  sum(ice.cold_idx_n[gs])
        push!(gs_lengths_per_core, gs_length)
    end
    push!(gs_lengths_n_10, gs_lengths_per_core)
end

gs_lengths_20 = []
for ice in ice_cores_20
    gs_lengths_per_core = Int64[]
    for gs in 1:17
        #@show ice.name
        gs_length =  sum(ice.cold_idx[gs])
        push!(gs_lengths_per_core, gs_length)
    end
    push!(gs_lengths_20, gs_lengths_per_core)
end


gs_lengths_n_20 = []
for ice in ice_cores_20
    gs_lengths_per_core = Int64[]
    for gs in 1:17
        #@show ice.name
        gs_length =  sum(ice.cold_idx_n[gs])
        push!(gs_lengths_per_core, gs_length)
    end
    push!(gs_lengths_n_20, gs_lengths_per_core)
end

# sranges = Tuple{Int64, Int64}[]
# for i in 10:10:100
#     for j in 20:10:110
#         if i<j
#             push!(sranges,(i,j))
#         end
#     end
# end
# sranges

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

sranges_5 = all_sranges

n_surs = 10_000
cold_type = "gs" #["gs", "gs_short"][1]
filter_indicator = false
only_full = true


function save_var_slopes_and_p(n_surrs::Int64, ice_core_list::Array{ice_core}; 
            windows::Array{Int64} = [200], sur_method::String = "TFTS", filt_type::String = ["filt","normed_filt"][2], cold_type::String = ["gs", "gs_short"][1],
            filter_indicator::Bool = false, only_full::Bool = true, folder::String = "NGRIP5")
    sur_method == uppercase(sur_method)
    if sur_method == "TFTS"
        method = TFTS(0.05)
    elseif sur_method == "TFTD"
        method = TFTDRandomFourier(true, 0.02)
    else
        phases = true
        method = RandomFourier(phases)
    end

    filt_type = lowercase(filt_type)
    cold_type = lowercase(cold_type)
    for ice in ice_core_list #[NGRIP_10_spl, NEEM_10_spl]# #ice_cores
        if filt_type == "filt"
            filt_data = ice.δ_filt_100
        elseif filt_type == "normed_filt"
            filt_data = ice.δ_normed_filt_100
        else
            error("filt_type must be \"filt\" or \"normed_filt\".")
        end

        if cold_type == "gs"
            cold = ice.cold_idx
            n_cold = length(cold)
        elseif cold_type == "gs_short"
            cold = ice.cold_idx_n
            n_cold = length(cold)
        else
            error("cold_type must be \"gs\" or \"gs_short\".")
        end

        for w in windows
            saveto_var = "new_surrogate_files/$(folder)/var/w_$(w)_$(filt_type)_$(ice.name)_$(cold_type)_FiltInd_$(filter_indicator)_onlyfull_$(only_full)_$(n_surrs)_$(sur_method).jld2"
            if !isfile(saveto_var)
                println("Calculating ", saveto_var)
                slopes_var = Array{Union{Nothing,Float64}}(nothing,n_cold )
                surr_slopes_var = Array{Union{Nothing,Float64}}(nothing,n_surrs, n_cold )
                p_one_var = Array{Union{Nothing,Float64}}(nothing,n_cold )
                p_two_var = Array{Union{Nothing,Float64}}(nothing,n_cold )
                vars = Array{Vector{Union{Missing,Float64}}}(undef, n_cold)
                times = Array{Vector{Float64}}(undef, n_cold)

                for (ic, c) in enumerate(cold)
                    times[ic] = @view ice.age[c]
                    x = filt_data[c]
                    x = convert.(Float64,x)
                    
                    dt = ice.resolution
                    scavf = Int(w / dt)

                    true_var = runvar(x,scavf, only_full=only_full)
                    if filter_indicator
                        v_nomiss = true_var[findall(!ismissing,true_var)]
                        if length(v_nomiss)>=25
                            v_nomiss = convert.(Float64, v_nomiss)
                            v_nomiss =  cheby_lowpass_filter(v_nomiss, .95 * 1. /800, 1. / dt, 8, .05)
                        else
                            true_var = missings(size(true_var))
                        end
                    end
                    true_slope_var = get_slopes(1:length(true_var), true_var)
                    
                    vars[ic] = true_var
                    slopes_var[ic] =  true_slope_var
                    
                    function q_var_no_filt(x)
                        v= runvar(x,scavf, only_full=only_full)
                        slope = get_slopes(1:length(v), v)#, return_pred = false)
                        return slope
                    end

                    function q_var_filt(x) #function for significance testing
                        v= runvar(x,scavf, only_full=only_full)
                        v_nomiss = @view v[findall(!ismissing,v)]
                        v_nomiss = convert.(Float64, v_nomiss)
                        v_nomiss = cheby_lowpass_filter(v_nomiss, 0.95 * 1. /800, 1. / dt, 8, .05)
                        slope = get_slopes(1:length(v), v)#, return_pred = false)
                        return slope
                    end


                    if filter_indicator == false
                       q_var = q_var_no_filt
                    elseif filter_indicator == true
                        q_var = q_var_filt
                    end

                    if (typeof(true_slope_var) != Nothing) #&& (true_slope > 0)
                        #only for increases, if they exist
                        var_test = TimeseriesSurrogates.SurrogateTest(q_var,x,method, n=n_surrs)
                        p_two_var[ic] = TimeseriesSurrogates.pvalue(var_test, tail=:both)
                        p_one_var[ic] = count(v -> v ≥ true_slope_var, var_test.vals)/length(var_test.vals)
                        
                        surr_slopes_var[:,ic] =  var_test.vals
                    end
                end
                obj_var = indicator_and_significance(ice.name,times,
                    vars, slopes_var, surr_slopes_var, "Var",0,100,p_one_var,p_two_var,ice.resolution)
                save(saveto_var, "slopes", obj_var)
            else
                println("$(saveto_var) already exists")
            end

        end
    end
    println("Done :)")
end


function save_ac_slopes_and_p(n_surrs::Int64, ice_core_list::Array{ice_core}; 
    windows::Array{Int64} = [200], sur_method::String = "TFTS", filt_type::String = ["filt","normed_filt"][2], cold_type::String = ["gs", "gs_short"][1],
    filter_indicator::Bool = false, only_full::Bool = true, folder::String = "NGRIP5")
    sur_method == uppercase(sur_method)
    if sur_method == "TFTS"
        method = TFTS(0.05)
    elseif sur_method == "TFTD"
        method = TFTDRandomFourier(true, 0.02)
    else
        phases = true
        method = RandomFourier(phases)
    end

    filt_type = lowercase(filt_type)
    cold_type = lowercase(cold_type)

    for ice in ice_core_list 
        if filt_type == "filt"
            filt_data = ice.δ_filt_100
        elseif filt_type == "normed_filt"
            filt_data = ice.δ_normed_filt_100
        else
            error("ERROR: filt_type must be \"filt\" or \"normed_filt\".")
        end

        if cold_type == "gs"
            cold = ice.cold_idx
            n_cold = length(cold)
        elseif cold_type == "gs_short"
            cold = ice.cold_idx_n
            n_cold = length(cold)
        else
            error("ERROR: cold_type must be \"gs\" or \"gs_short\".")
        end

        for w in windows
            w = Int(w)
            saveto = "new_surrogate_files/$(folder)/ac/w_$(w)_$(filt_type)_$(ice.name)_$(cold_type)_FiltInd_$(filter_indicator)_onlyfull_$(only_full)_$(n_surrs)_$(sur_method).jld2"
            if !isfile(saveto)
                println("Calculating ", saveto)
    

                slopes = Array{Union{Nothing,Float64}}(nothing,n_cold)
                surr_slopes = Array{Union{Nothing,Float64}}(nothing,n_surrs, n_cold )
                p_one_ac = Array{Union{Nothing,Float64}}(nothing,n_cold )
                p_two_ac = Array{Union{Nothing,Float64}}(nothing,n_cold )
                acs = Array{Vector{Union{Missing,Float64}}}(undef, n_cold)
                times = Array{Vector{Float64}}(undef, n_cold)

                for (ic, c) in enumerate(cold)
                    times[ic] = @view ice.age[c]
                    x = filt_data[c]
                    x = convert.(Float64,x)
                    
                    dt = ice.resolution
                    samp = Int(dt)

                    scavf = Int(div(w,dt))

                    true_ac = runac(x,scavf, only_full=only_full)
                    if filter_indicator
                        v_nomiss = @view true_ac[findall(!ismissing,true_ac)]
                        if length(v_nomiss)>=25
                            v_nomiss = convert.(Float64, v_nomiss)
                            v_nomiss =  cheby_lowpass_filter(v_nomiss, .95 * 1. /800, 1. / samp, 8, .05)
                        else
                            true_ac = missings(size(true_ac))
                        end


                    end
                    true_slope_ac = get_slopes(1:length(true_ac), true_ac)

                    acs[ic] = true_ac
                    slopes[ic] =  true_slope_ac
                    
                    function q_ac_nofilt(x)
                        v= runac(x,scavf, only_full=only_full)
                        slope = get_slopes(1:length(v), v)#, return_pred = false)
                        return slope
                    end
                    function q_ac_filt(x)
                        v= runac(x,scavf, only_full=only_full)
                        v_nomiss = @view v[findall(!ismissing,v)]
                        v_nomiss = convert.(Float64, v_nomiss)
                        v_nomiss = cheby_lowpass_filter(v_nomiss, 0.95 * 1. /800, 1. / samp, 8, .05)
                        slope = get_slopes(1:length(v), v)#, return_pred = false)
                        return slope
                    end

                    if filter_indicator == false
                        q_ac = q_ac_nofilt

                    elseif filter_indicator == true
                        q_ac = q_ac_filt
                    end

                    if (typeof(true_slope_ac) != Nothing) #&& (true_slope > 0)
                        #only for increases, if they exist
                        ac_test = TimeseriesSurrogates.SurrogateTest(q_ac,x,method, n=n_surrs)
                        p_two_ac[ic] = TimeseriesSurrogates.pvalue(ac_test, tail=:both)
                        p_one_ac[ic] = count(v -> v ≥ true_slope_ac, ac_test.vals)/n_surrs

                        surr_slopes[:,ic] = ac_test.vals
                    end
                end
                obj = indicator_and_significance(ice.name,times,
                    acs, slopes, surr_slopes, "AC",0,100,p_one_ac,p_two_ac,ice.resolution)
                save(saveto, "slopes", obj)
            else
                println("$(saveto) already exists")
            end

        end
    end
    println("Done :)")
end


function get_fake_ac_slope_pvals(x::Array{Float64},len_gs::Array{Int64}, max_len_gs::Int64, scavf::Int64, filter_indicator::Bool, method::Surrogate, n_surrs_per_gs::Int64, only_full::Bool)
    n_gs = length(len_gs)
    fake_gs_end_locs = rand(max_len_gs+1:length(x), 17)
    fake_gs_start_locs = fake_gs_end_locs .- len_gs .+1
    fake_gs_locs = [fake_gs_start_locs[i]: fake_gs_end_locs[i] for i in 1:n_gs]
    slopes = Array{Union{Nothing,Float64}}(nothing,n_gs)
    pvals_one = Array{Union{Nothing,Float64}}(nothing,n_gs)
    pvals_two = Array{Union{Nothing,Float64}}(nothing,n_gs)


    function q_ac_no_filt(xx)
        v= runac(xx,scavf, only_full=only_full)
        slope = get_slopes(1:length(v), v)#, return_pred = false)
        return slope
    end
    function q_ac_filt(xx)
        v= runac(xx,scavf, only_full=only_full)
        v_nomiss = @view v[findall(!ismissing,v)]
        if length(v_nomiss)>=25
            v_nomiss = convert.(Float64, v_nomiss)
            v_nomiss = cheby_lowpass_filter(v_nomiss, .95 * 1. /800, 1. / samp, 8, .05)
            slope = get_slopes(1:length(v), v)#, return_pred = false)
        else
            slope = nothing
        end
        return slope
    end

    if filter_indicator == false
        q_ac = q_ac_no_filt

    elseif filter_indicator == true
        q_ac = q_ac_filt
    end
    
    for (j,gs) in enumerate(fake_gs_locs)
        #calc ews (and significance with how many??? 
        # and at different pvals?!) here
        #@show gs
        xx = x[gs]
        #t = ice.age[gs]
        true_slope = q_ac(xx)
        slopes[j] = true_slope
        if (typeof(true_slope) != Nothing) #&& (true_slope > 0)
            #only for increases, if they exist
            ac_test = TimeseriesSurrogates.SurrogateTest(q_ac,xx,method, n=n_surrs_per_gs)
            pvals_two[j] = TimeseriesSurrogates.pvalue(ac_test, tail=:both)
            pvals_one[j] = count(v -> v ≥ true_slope, ac_test.vals)/n_surrs_per_gs
        end
    end
    return slopes, pvals_one, pvals_two
end


function sign_number_of_ac_numeric(ice_core_list::Array{ice_core}, n_surrs::Int64, n_surrs_per_gs::Int64,  p_thresholds::Array{Float64};
                    windows::Array{Int64} = [200], sur_method::String = "TFTS", filt_type::String = ["filt","normed_filt"][2], cold_type::String = ["gs", "gs_short"][1],
                    filter_indicator::Bool = false, only_full::Bool = true, folder::String = "NGRIP5", gs_lengths = gs_lengths, gs_lengths_n=gs_lengths_n)
    #1. make n_surrs surrogates of δ
    #2. place 17 GS randomly
    #3. calc ews (and significance at different pvals) there
    #4. get number of sign. positive trends per time series
    #5. get distribution of that
    sur_method == uppercase(sur_method)
    if sur_method == "TFTS"
        method = TFTS(0.05)
    elseif sur_method == "TFTD"
        method = TFTDRandomFourier(true, 0.02)
    else
        phases = true
        method = RandomFourier(phases)
    end

    filt_type = lowercase(filt_type)
    cold_type = lowercase(cold_type)

    
    for (i,ice) in enumerate(ice_core_list)
        if filt_type == "filt"
            filt_data = ice.δ_filt_100
        elseif filt_type == "normed_filt"
            filt_data = ice.δ_normed_filt_100
        else
            error("filt_type must be \"filt\" or \"normed_filt\".")
        end
        
        if cold_type == "gs"
            len_gs = gs_lengths[i]
            max_len_gs = maximum(len_gs)
        elseif cold_type == "gs_short"
            len_gs = gs_lengths_n[i]
            max_len_gs = maximum(len_gs)
        else
            error("cold_type must be \"gs\" or \"gs_short\".")
        end

        dt = ice.resolution
        sgen = surrogenerator(filt_data, method)
        siter = (sgen() for _ in 1:n_surrs)

        for w in windows
            scavf = Int(w / dt)
            qs_one = zeros(Int32,length(p_thresholds),n_surrs)
            qs_two = zeros(Int32,length(p_thresholds),n_surrs)
            pstring = string.(p_thresholds)
            pstring = join(pstring,"_")

            
            saveto = "new_surrogate_files/$(folder)/n_sig_dists/ac/w_$(w)_$(filt_type)_$(ice.name)_$(cold_type)_FiltInd_$(filter_indicator)_onlyfull_$(only_full)_n_$(n_surrs)_nGS_$(n_surrs_per_gs)_p_$(pstring)_$(sur_method).jld2"
            saveto_temp = "new_surrogate_files/$(folder)/n_sig_dists/ac/TEMP_w_$(w)_$(filt_type)_$(ice.name)_$(cold_type)_FiltInd_$(filter_indicator)_onlyfull_$(only_full)_n_$(n_surrs)_nGS_$(n_surrs_per_gs)_p_$(pstring)_$(sur_method).jld2"
            if !isfile(saveto)
                println("Calculating $saveto")
                if !isfile(saveto_temp)
                    obj = distribution_significant_increases(ice.name, qs_one, qs_two, 
                            0, n_surrs_per_gs, p_thresholds, "AC", 0, 100, ice.resolution)
                else 
                    obj = load(saveto_temp)["distribution"]
                end
                for (isu,xsur) in enumerate(siter)
                    #@show isu
                    if isu > obj.n
                        slopes, pvals_one, pvals_two = get_fake_ac_slope_pvals(xsur,len_gs, max_len_gs, scavf, filter_indicator,method, n_surrs_per_gs, only_full)
                        for (ip,p) in enumerate(p_thresholds)
                            function q_nsig(slopes, pvals,p)
                                real_slopes = findall(!isnothing, slopes) 
                                pos_slopes = findall(>(0),slopes[real_slopes])
                                good_pval = pvals[real_slopes][pos_slopes] .< p
                                sum(good_pval)
                            end
                            obj.num_inc_one[ip,isu] = q_nsig(slopes,pvals_one,p)
                            obj.num_inc_two[ip,isu] = q_nsig(slopes,pvals_two,p)
                        end
                        #save every 20th iteration
                        if isu %20 == 0
                            println("calculated $(isu) surrogates so far.")
                            obj.n = isu
                            save(saveto_temp, "distribution", obj)
                        end
                    end
                end
                obj.n = n_surrs
                save(saveto, "distribution", obj)
                if isfile(saveto_temp)
                    rm(saveto_temp)
                end
            else
                println("$saveto already exists")
            end
        end
    end
    println("Done :)")
end


function get_fake_var_slope_pvals(x::Array{Float64},len_gs::Array{Int64}, max_len_gs::Int64, scavf::Int64, filter_indicator::Bool, method::Surrogate, n_surrs_per_gs::Int64,only_full::Bool)
    n_gs = length(len_gs)
    fake_gs_end_locs = rand(max_len_gs+1:length(x), 17)
    fake_gs_start_locs = fake_gs_end_locs .- len_gs .+1
    fake_gs_locs = [fake_gs_start_locs[i]: fake_gs_end_locs[i] for i in 1:n_gs]
    slopes = Array{Union{Nothing,Float64}}(nothing,n_gs)
    pvals_one = Array{Union{Nothing,Float64}}(nothing,n_gs)
    pvals_two = Array{Union{Nothing,Float64}}(nothing,n_gs)


    function q_var_no_filt(xx)
        v= runvar(xx,scavf, only_full=only_full)
        slope = get_slopes(1:length(v), v)#, return_pred = false)
        return slope
    end
    function q_var_filt(xx)
        v= runvar(xx,scavf, only_full=only_full)
        v_nomiss = @view v[findall(!ismissing,v)]
        if length(v_nomiss)>=25
            v_nomiss = convert.(Float64, v_nomiss)
            v_nomiss = cheby_lowpass_filter(v_nomiss, .95 * 1. /800, 1. / samp, 8, .05)
            slope = get_slopes(1:length(v), v)#, return_pred = false)
        else
            slope = nothing
        end
        return slope
    end

    if filter_indicator == false
        q_var = q_var_no_filt

    elseif filter_indicator == true
        q_var = q_var_filt
    end

    for (j,gs) in enumerate(fake_gs_locs)
        xx = x[gs]
        #t = ice.age[gs]
        true_slope = q_var(xx)
        slopes[j] = true_slope
        if (typeof(true_slope) != Nothing) #&& (true_slope > 0)
            #only for increases, if they exist
            test = TimeseriesSurrogates.SurrogateTest(q_var,xx,method, n=n_surrs_per_gs)
            pvals_two[j] = TimeseriesSurrogates.pvalue(test, tail=:both)
            pvals_one[j] = count(v -> v ≥ true_slope, test.vals)/n_surrs_per_gs
        end
    end
    return slopes, pvals_one, pvals_two
end

function sign_number_of_var_numeric(ice_core_list, n_surrs, n_surrs_per_gs,  p_thresholds;
            windows = [200], sur_method = "TFTS", filt_type = ["filt","normed_filt"][2], cold_type = ["gs", "gs_short"][1],
            filter_indicator = false, only_full::Bool = true, folder::String = "NGRIP5", gs_lengths = gs_lengths, gs_lengths_n = gs_lengths_n)
    #1. make n_surrs surrogates of δ
    #2. place 17 GS randomly
    #3. calc ews (and significance at different pvals) there
    #4. get number of positive trends per time series
    #5. get distribution of that
    sur_method == uppercase(sur_method)
    if sur_method == "TFTS"
        method = TFTS(0.05)
    elseif sur_method == "TFTD"
        method = TFTDRandomFourier(true, 0.02)
    else
        phases = true
        method = RandomFourier(phases)
    end

    filt_type = lowercase(filt_type)
    cold_type = lowercase(cold_type)


    for (i,ice) in enumerate(ice_core_list)
        if filt_type == "filt"
            filt_data = ice.δ_filt_100
        elseif filt_type == "normed_filt"
            filt_data = ice.δ_normed_filt_100
        else
            error("filt_type must be \"filt\" or \"normed_filt\".")
        end
        
        if cold_type == "gs"
            len_gs = gs_lengths[i]
            max_len_gs = maximum(len_gs)
        elseif cold_type == "gs_short"
            len_gs = gs_lengths_n[i]
            max_len_gs = maximum(len_gs)
        else
            error("cold_type must be \"gs\" or \"gs_short\".")
        end

        dt = ice.resolution
        sgen = surrogenerator(filt_data, method)
        siter = (sgen() for _ in 1:n_surrs)


        for w in windows
            scavf = Int(w / dt)
            qs_one = zeros(Int32,length(p_thresholds),n_surrs)
            qs_two = zeros(Int32,length(p_thresholds),n_surrs)
            pstring = string.(p_thresholds)
            pstring = join(pstring,"_")
            saveto = "new_surrogate_files/$(folder)/n_sig_dists/var/w_$(w)_$(filt_type)_$(ice.name)_$(cold_type)_FiltInd_$(filter_indicator)_onlyfull_$(only_full)_n_$(n_surrs)_nGS_$(n_surrs_per_gs)_p_$(pstring)_$(sur_method).jld2"
            saveto_temp = "new_surrogate_files/$(folder)/n_sig_dists/var/TEMP_w_$(w)_$(filt_type)_$(ice.name)_$(cold_type)_FiltInd_$(filter_indicator)_onlyfull_$(only_full)_n_$(n_surrs)_nGS_$(n_surrs_per_gs)_p_$(pstring)_$(sur_method).jld2"
            if !isfile(saveto)
                println("Calculating $saveto")
                if !isfile(saveto_temp)
                    obj = distribution_significant_increases(ice.name, qs_one, qs_two, 
                            0, n_surrs_per_gs, p_thresholds, "Var", 0, 100, ice.resolution)
                else 
                    obj = load(saveto_temp)["distribution"]
                end

                for (isu,xsur) in enumerate(siter)
                    if isu > obj.n
                        slopes, pvals_one, pvals_two = get_fake_var_slope_pvals(xsur,len_gs, max_len_gs, scavf, filter_indicator,method, n_surrs_per_gs, only_full)
                        for (ip,p) in enumerate(p_thresholds)
                            function q_nsig(slopes, pvals,p)
                                real_slopes = findall(!isnothing, slopes) 
                                pos_slopes = findall(>(0),slopes[real_slopes])
                                good_pval = pvals[real_slopes][pos_slopes] .< p
                                sum(good_pval)
                            end
                            obj.num_inc_one[ip,isu] = q_nsig(slopes,pvals_one,p)
                            obj.num_inc_two[ip,isu] = q_nsig(slopes,pvals_two,p)
                        end
                        #save every 20th iteration
                        if isu %20 == 0
                            println("calculated $(isu) surrogates so far.")
                            obj.n = isu
                            save(saveto_temp, "distribution", obj)
                        end
                    end
                end

                obj.n = n_surrs
                save(saveto, "distribution", obj)
                if isfile(saveto_temp)
                    rm(saveto_temp)
                end
            else
                println("$saveto already exists")
            end

        end
    end
    println("Done :)")
end


function get_w(x::Array{Float64},dt::Int64, s1::Int64, s2::Int64,
    mother::String, param::Int,
    exclude_coi::Bool = true, smooth_w::Bool = false, filter_indicator::Bool = false; 
    pad = 1, s0=-1, j1 = -1, dj::Float64 = 0.1,
    normalise::Bool = false, only_full::Bool = true)

    scavf = Int(200/dt)
    
    mother = uppercase(mother)
    #get the reconstrucion factor C_δ
    if mother == "MORLET"
        if param == -1
            param = 6
        end
        k0 = param
        if k0 == 6; C_δ = 0.776; end
    elseif mother == "PAUL"
        if param == -1
            param = 4
        end
        m = param
        if m == 4; C_δ = 1.132; end
    elseif mother == "DOG"
        if param == -1
            param = 2
        end
        m = param
        if m == 2; C_δ = 3.541; end
        if m == 6; C_δ = 1.966; end
    else
        error("Mother must be one of MORLET,PAUL,DOG")
    end
    #@show C_δ

    if normalise
        x = (x .- mean(x)) ./ std(x)
    end
    
    wave, period, scale, coi = wavelet(x,dt,pad=pad,dj=dj,s0=s0,j1=j1,mother=mother, param=param)
    pow = abs.(wave).^2

    if exclude_coi
        avg = @. s1 <= scale <=s2
        scale2 = scale[avg]

        pow2 = Matrix{Union{Float64, Missing}}(copy(pow[avg, :]./ scale2))
        coi2 = repeat(coi, inner = (1,size(pow2)[1]))'
        period2 = repeat(period[avg], inner = (1,size(pow2)[2]))
        pow2[findall(period2 .> coi2)] .= missing 
        
        ŵ² = (dj*dt / C_δ) .* sum(pow2, dims = 1)[:]
    else
        avg = @. s1 <= scale <=s2
        scale2 = scale[avg]
        pow2 = Matrix{Float64}(copy(pow[avg, :]./ scale2))
        ŵ² = (dj*dt / C_δ) .* sum(pow2, dims = 1)[:]
    end

    if smooth_w
        scavf = Int(200/dt) 
        ŵ² = runmean(ŵ², scavf, only_full = only_full) #??
    end

    if filter_indicator
        lpf = 800.0
        w_nomiss_ind = findall(!ismissing, ŵ²)
        w_nomiss = @view ŵ²[w_nomiss_ind]
        if length(w_nomiss)>=25#>0
            if maximum(diff(w_nomiss_ind)) >1
                error("there are missings in the middle of ŵ²")
            end
            w_nomiss = cheby_lowpass_filter(w_nomiss, .95 * 1. / lpf, 1. / dt, 8, .05)

        end
    end
    return ŵ²
end


function get_h(x::Array{Float64},dt::Int64, s1::Int64, s2::Int64,
    mother::String, param,
    exclude_coi::Bool = true, smooth_w::Bool = false, filter_indicator::Bool = false; 
    pad = 1, s0=-1, j1 = -1, dj::Float64 = 0.1,
    normalise::Bool = false, only_full::Bool = true)

    scavf = Int(200/dt)
    
    mother = uppercase(mother)
    #get the reconstrucion factor C_δ
    if mother == "MORLET"
        if param == -1
            param = 6
        end
        k0 = param
        if k0 == 6; C_δ = 0.776; end
    elseif mother == "PAUL"
        if param == -1
            param = 4
        end
        m = param
        if m == 4; C_δ = 1.132; end
    elseif mother == "DOG"
        if param == -1
            param = 2
        end
        m = param
        if m == 2; C_δ = 3.541; end
        if m == 6; C_δ = 1.966; end
    else
        error("Mother must be one of MORLET,PAUL,DOG")
    end

    if normalise
        x = (x .- mean(x)) ./ std(x)
    end
    
    wave, period, scale, coi = wavelet(x,dt,pad=pad,dj=dj,s0=s0,j1=j1,mother=mother, param=param)
    pow = abs.(wave).^2

    if exclude_coi
        avg = @. s1 <= scale <=s2
        scale2 = scale[avg]

        pow2 = Matrix{Union{Float64, Missing}}(copy(pow[avg, :]./ scale2))
        coi2 = repeat(coi, inner = (1,size(pow2)[1]))'
        period2 = repeat(period[avg], inner = (1,size(pow2)[2]))
        pow2[findall(period2 .> coi2)] .= missing 

    else
        avg = @. s1 <= scale <=s2
        scale2 = scale[avg]
        pow2 = copy(pow[avg, :]./ scale2)
    end

    if smooth_w
        scavf = Int(200/dt) 
        scav = zeros(Union{Missing, Float64},size(pow2))#zeros()
        for i in 1:size(scav)[1]
            scav[i,:] = runmean(pow2[i,:], scavf, only_full = only_full)
        end
        pow2 = scav
    end

    if exclude_coi
        pa = Union{Missing,Float64,     Nothing}[all(!ismissing, pow2[:,i]) ? get_slopes(log.(scale2), log.(pow2)[:,i]) : missing for i = 1:size(pow2)[2]]
            pa[findall(isnothing,pa)] .= missing
        hurst = @. (pa +1)/2
            hurst = convert(Array{Union{Missing, Float64}},hurst)
    else
        pa = Union{Missing, Nothing, Float64}[get_slopes(log.(scale2), log.(pow2)[:,i]) for i = 1:size(pow2)[2]]
        pa[findall(isnothing,pa)] .= missing
        hurst = @. (pa +1)/2
        hurst = convert(Array{Union{Missing, Float64}},hurst)
    end


    if filter_indicator
        lpf = 800.
        h_nomiss_ind = findall(!ismissing, hurst)
        h_nomiss = @view hurst[h_nomiss_ind]
        if length(h_nomiss)>=25 #>0
            if maximum(diff(h_nomiss_ind)) >1
                error("there are missings in the middle of hurst")
            end
            h_nomiss= cheby_lowpass_filter(h_nomiss, .95 * 1. / lpf, 1. / samp, 8, .05)
        end
    end
    return hurst
end


function save_sca_slopes_and_p(n_surrs::Int64, ice_core_list::Array{ice_core}, sranges::Vector{Tuple{Int64, Int64}}; mother::String = ["PAUL","MORLET"][1], param::Int64 = -1, 
    sur_method::String = "TFTS",
    dat_type::String = ["orig", "normed", "filt", "normed_filt"][1], cold_type = ["gs", "gs_short"][1],
    filter_indicator::Bool = false, smooth_w::Bool = false,  exclude_coi::Bool = true, normalise::Bool=false, only_full::Bool = true, folder::String = "NGRIP5")

    sur_method == uppercase(sur_method)
    if sur_method == "TFTS"
        method = TFTS(0.05)
    elseif sur_method == "TFTD"
        method = TFTDRandomFourier(true, 0.02)
    else
        phases = true
        method = RandomFourier(phases)
    end

    dj=0.1
    filt_type = lowercase(dat_type)
    cold_type = lowercase(cold_type)
    #for (is,splitt) in enumerate(["Niklas"])#,"UiT"])
        for ice in ice_core_list 
            if filt_type == "orig"
                data = ice.δ
            elseif filt_type == "normed"
                data = ice.δ_normed
            elseif filt_type == "filt"
                data = ice.δ_filt_100
            elseif filt_type == "normed_filt"
                data = ice.δ_normed_filt_100
            else
                error("ERROR: dat_type must be \"orig\", \"normed\", \"filt\" or \"normed_filt\".")
            end
    
            if cold_type == "gs"
                cold = ice.cold_idx
                n_cold = length(cold)
            elseif cold_type == "gs_short"
                cold = ice.cold_idx_n
                n_cold = length(cold)
            else
                error("ERROR: cold_type must be \"gs\" or \"gs_short\".")
            end
    
            for (i,(s1,s2)) in enumerate(sranges)
                    if s1 < s2
                        saveto = "new_surrogate_files/$(folder)/sca/s1_$(s1)_s2_$(s2)_$(ice.name)_$(mother)_$(filt_type)_$(cold_type)_FiltInd_$(filter_indicator)_normalise_$(normalise)_nocoi_$(exclude_coi)_onlyfull_$(only_full)_smoothw_$(smooth_w)_$(n_surrs)_$(sur_method).jld2"
                        if !isfile(saveto)
                            println("Calculating $saveto")
                            slopes =Array{Union{Nothing,Float64}}(nothing, n_cold)
                            surr_slopes = Array{Union{Nothing,Float64}}(nothing,n_surrs, n_cold )
                            pvals_one = Array{Union{Nothing,Float64}}(nothing, n_cold)
                            pvals_two = Array{Union{Nothing,Float64}}(nothing, n_cold)
                            scas = Array{Vector{Union{Missing,Float64}}}(undef, n_cold)
                            times = Array{Vector{Float64}}(undef, n_cold)

                            for (ic, c) in enumerate(cold)
                                times[ic] = ice.age[c]
                                
                                x = data[c]
                                x = convert.(Float64,x)
                                dt = ice.resolution
                                
                                true_sca = get_w(x,dt, s1,s2,mother, param, exclude_coi, smooth_w,
                                             filter_indicator, pad = 1, s0=-1, j1=-1, dj=0.1, normalise = normalise, only_full=only_full)
 
                                scas[ic] = true_sca
                                
                                function q_sca(x)
                                    sca = get_w(x,dt, s1,s2,mother, param, exclude_coi, smooth_w,
                                        filter_indicator, pad = 1, s0=-1, j1=-1, dj=0.1,normalise=normalise, only_full=only_full)
                                    slope = get_slopes(1:length(sca), sca)#, return_pred = false)
                                    slope
                                end
                                true_slope = q_sca(x)
                                slopes[ic] = true_slope
                                if (typeof(true_slope) != Nothing) #&& (true_slope > 0)
                                    #only for increases, if they exist
                                    sca_test_tfts = TimeseriesSurrogates.SurrogateTest(q_sca,x,method, n=n_surrs)
                                    pvals_two[ic] = TimeseriesSurrogates.pvalue(sca_test_tfts, tail=:both)
                                    pvals_one[ic] = count(v -> v ≥ true_slope, sca_test_tfts.vals)/n_surrs
                                    surr_slopes[:,ic] = sca_test_tfts.vals
                                end
                            end
                            obj = indicator_and_significance(ice.name,times,
                                scas, slopes, surr_slopes, "sca",s1,s2,pvals_one, pvals_two,ice.resolution)
                            save(saveto, "slopes", obj)
                        else
                            println("$saveto already exists")
                        end
                    end 
            end
        end
    println("Done :)")
end


function save_hurst_slopes_and_p(n_surrs::Int64, ice_core_list::Array{ice_core}, sranges::Vector{Tuple{Int64, Int64}}; mother::String = ["PAUL","MORLET"][1], param::Int64 = -1, 
    sur_method::String = "TFTS",
    dat_type::String = ["orig", "normed", "filt", "normed_filt"][1], cold_type = ["gs", "gs_short"][1],
    filter_indicator::Bool = false, smooth_w::Bool = false,  exclude_coi::Bool = true, normalise::Bool = false, only_full::Bool = true, folder::String = "NGRIP5")

    sur_method == uppercase(sur_method)
    if sur_method == "TFTS"
        method = TFTS(0.05)
    elseif sur_method == "TFTD"
        method = TFTDRandomFourier(true, 0.02)
    else
        phases = true
        method = RandomFourier(phases)
    end

    dj=0.1
    filt_type = lowercase(dat_type)
    cold_type = lowercase(cold_type)
    #for (is,splitt) in enumerate(["Niklas"])#,"UiT"])
        for ice in ice_core_list
            if filt_type == "orig"
                data = ice.δ
            elseif filt_type == "normed"
                data = ice.δ_normed
            elseif filt_type == "filt"
                data = ice.δ_filt_100
            elseif filt_type == "normed_filt"
                data = ice.δ_normed_filt_100
            else
                error("ERROR: dat_type must be \"orig\", \"normed\", \"filt\" or \"normed_filt\".")
            end
    
            if cold_type == "gs"
                cold = ice.cold_idx
                n_cold = length(cold)
            elseif cold_type == "gs_short"
                cold = ice.cold_idx_n
                n_cold = length(cold)
            else
                error("ERROR: cold_type must be \"gs\" or \"gs_short\".")
            end
    
            for (i,(s1,s2)) in enumerate(sranges)
                    if s1 < s2
                        saveto = "new_surrogate_files/$(folder)/hurst/s1_$(s1)_s2_$(s2)_$(ice.name)_$(mother)_$(filt_type)_$(cold_type)_FiltInd_$(filter_indicator)_normalise_$(normalise)_nocoi_$(exclude_coi)_onlyfull_$(only_full)_smoothw_$(smooth_w)_$(n_surrs)_$(sur_method).jld2"
                        if !isfile(saveto)
                            println("Calculating $saveto")
                            slopes =Array{Union{Nothing,Float64}}(nothing, n_cold)
                            surr_slopes = Array{Union{Nothing,Float64}}(nothing,n_surrs, n_cold )
                            pvals_one = Array{Union{Nothing,Float64}}(nothing, n_cold)
                            pvals_two = Array{Union{Nothing,Float64}}(nothing, n_cold)
                            hursts = Array{Vector{Union{Missing,Float64}}}(undef, n_cold)
                            times = Array{Vector{Float64}}(undef, n_cold)

                            for (ic, c) in enumerate(cold)
                                times[ic] = ice.age[c]
                                
                                x = data[c]
                                x = convert.(Float64,x)
                                dt = ice.resolution
                                

                                true_hurst = get_h(x,dt, s1,s2,mother, param, exclude_coi, smooth_w,
                                             filter_indicator, pad = 1, s0=-1, j1=-1, dj=0.1,normalise = normalise, only_full=only_full)

                                hursts[ic] = true_hurst
                                
                                function q_hurst(x)
                                    hurst = get_h(x,dt, s1,s2,mother, param, exclude_coi, smooth_w,
                                        filter_indicator, pad = 1, s0=-1, j1=-1, dj=0.1,normalise = normalise, only_full=only_full)
                                    slope = get_slopes(1:length(hurst), hurst)#, return_pred = false)
                                    slope
                                end
                                true_slope = q_hurst(x)
                                slopes[ic] = true_slope
                                if (typeof(true_slope) != Nothing) #&& (true_slope > 0)
                                    #only for increases, if they exist
                                    test_tfts = TimeseriesSurrogates.SurrogateTest(q_hurst,x,method, n=n_surrs)
                                    pvals_two[ic] = TimeseriesSurrogates.pvalue(test_tfts, tail=:both)
                                    pvals_one[ic] = count(v -> v ≥ true_slope, test_tfts.vals)/n_surrs
                                    surr_slopes[:,ic] = test_tfts.vals
                                end
                            end
                            obj = indicator_and_significance(ice.name,times,
                                hursts, slopes, surr_slopes, "hurst",s1,s2,pvals_one, pvals_two,ice.resolution)
                            save(saveto, "slopes", obj)
                        else
                            println("$saveto already exists")
                        end
                    end 
            end
        end
    println("Done :)")
end

# Overall significance test: 
function get_sca_slope_pvals(x::Array{Float64},s1::Int64,s2::Int64,
        len_gs::Array{Int64}, max_len_gs::Int64, n_surrs_per_gs::Int64,
        dj::Float64, dt::Int64,
        mother::String, param::Int64,
        exclude_coi::Bool,
        method::Surrogate,
        filter_indicator::Bool, smooth_w::Bool, normalise::Bool, only_full::Bool)
    
        n_gs = length(len_gs)
    fake_gs_end_locs = rand(max_len_gs+1:length(x), 17)
    fake_gs_start_locs = fake_gs_end_locs .- len_gs .+1
    fake_gs_locs = [fake_gs_start_locs[i]: fake_gs_end_locs[i] for i in 1:17]

    slopes = Array{Union{Nothing,Float64}}(nothing, n_gs)
    pvals_one =  Array{Union{Nothing,Float64}}(nothing, n_gs)
    pvals_two =  Array{Union{Nothing,Float64}}(nothing, n_gs)
    for (ic,gs) in enumerate(fake_gs_locs)
        #calc ews (and significance at different pvals) here
        xx = x[gs]
        xx = convert.(Float64,xx)

        function q_sca(xx)
            sca = get_w(xx,dt, s1,s2,mother, param, exclude_coi, smooth_w,
                filter_indicator, pad = 1, s0=-1, j1=-1, dj=0.1,normalise = normalise, only_full=only_full)
            slope = get_slopes(1:length(sca), sca)#, return_pred = false)
            slope
        end

        true_slope = q_sca(xx)
        slopes[ic]= true_slope
        if (typeof(true_slope) != Nothing) #&& (true_slope > 0)
            #only for increases, if they exist
            sca_test_tfts = TimeseriesSurrogates.SurrogateTest(q_sca,xx,method, n=n_surrs_per_gs)
            pvals_two[ic] = TimeseriesSurrogates.pvalue(sca_test_tfts, tail=:both)
            pvals_one[ic] = count(v -> v ≥ true_slope, sca_test_tfts.vals)/n_surrs_per_gs
                            
        end
    end
    return slopes, pvals_one, pvals_two
end

function sign_number_of_sca_numeric(ice_core_list::Array{ice_core}, n_surrs::Int64, n_surrs_per_gs::Int64,  p_thresholds::Array{Float64}, sranges::Vector{Tuple{Int64, Int64}};
                mother::String = ["PAUL","MORLET"][1], param::Int64 = -1, 
                sur_method::String = "TFTS", 
                dat_type::String = ["orig", "normed", "filt", "normed_filt"][1], cold_type::String = ["gs", "gs_short"][1],
                filter_indicator::Bool = false, smooth_w::Bool = false,  exclude_coi::Bool = true, normalise::Bool = false, only_full::Bool = true,
                folder::String = "NGRIP5", gs_lengths = gs_lengths, gs_lengths_n = gs_lengths_n, batch_no::Int64 = 1)
    
    #1. make n_surrs surrogates of δ
    #2. place 17 GS randomly
    #3. calc ews (and significance at different pvals) there
    #4. get number of positive trends per time series
    #5. get distribution of that
    sur_method == uppercase(sur_method)
    if sur_method == "TFTS"
        method = TFTS(0.05)
    elseif sur_method == "TFTD"
        method = TFTDRandomFourier(true, 0.02)
    else
        method = RandomFourier(true)
    end

    filt_type = lowercase(dat_type)
    cold_type = lowercase(cold_type)

    for (i,ice) in enumerate(ice_core_list)
       if filt_type == "orig"
            data = ice.δ
        elseif filt_type == "normed"
            data = ice.δ_normed
        elseif filt_type == "filt"
            data = ice.δ_filt_100
        elseif filt_type == "normed_filt"
            data = ice.δ_normed_filt_100
        else
         error("ERROR: dat_type must be \"orig\", \"normed\", \"filt\" or \"normed_filt\".")
        end
        
        if cold_type == "gs"
            len_gs = gs_lengths[i]
            max_len_gs = maximum(len_gs)
            #cold = ice.cold_idx
        elseif cold_type == "gs_short"
            len_gs = gs_lengths_n[i]
            max_len_gs = maximum(len_gs)
            #cold = ice.cold_idx_n
        else
            error("cold_type must be \"gs\" or \"gs_short\".")
        end

        dj = 0.1
        dt = ice.resolution

        sgen = surrogenerator(data, method)
        siter = (sgen() for _ in 1:n_surrs)

        for (i,(s1,s2)) in enumerate(sranges)
            
            if s1 < s2
                qs_one = zeros(Int32,length(p_thresholds),n_surrs)
                qs_two = zeros(Int32,length(p_thresholds),n_surrs)
                pstring = string.(p_thresholds)
                pstring = join(pstring,"_")
                if batch_no>0
                    saveto = "new_surrogate_files/$(folder)/n_sig_dists/sca/batch_$(batch_no)_s1_$(s1)_s2_$(s2)_$(ice.name)_$(mother)_$(filt_type)_$(cold_type)_FiltInd_$(filter_indicator)_normalise_$(normalise)_nocoi_$(exclude_coi)_onlyfull_$(only_full)_smoothw_$(smooth_w)_n_$(n_surrs)_nGS_$(n_surrs_per_gs)_p_$(pstring)_$(sur_method).jld2"
                    saveto_temp = "new_surrogate_files/$(folder)/n_sig_dists/sca/TEMP_batch_$(batch_no)__s1_$(s1)_s2_$(s2)_$(ice.name)_$(mother)_$(filt_type)_$(cold_type)_FiltInd_$(filter_indicator)_normalise_$(normalise)_nocoi_$(exclude_coi)_onlyfull_$(only_full)_smoothw_$(smooth_w)_n_$(n_surrs)_nGS_$(n_surrs_per_gs)_p_$(pstring)_$(sur_method).jld2"
                else
                    saveto = "new_surrogate_files/$(folder)/n_sig_dists/sca/s1_$(s1)_s2_$(s2)_$(ice.name)_$(mother)_$(filt_type)_$(cold_type)_FiltInd_$(filter_indicator)_normalise_$(normalise)_nocoi_$(exclude_coi)_onlyfull_$(only_full)_smoothw_$(smooth_w)_n_$(n_surrs)_nGS_$(n_surrs_per_gs)_p_$(pstring)_$(sur_method).jld2"
                    saveto_temp = "new_surrogate_files/$(folder)/n_sig_dists/sca/TEMP_s1_$(s1)_s2_$(s2)_$(ice.name)_$(mother)_$(filt_type)_$(cold_type)_FiltInd_$(filter_indicator)_normalise_$(normalise)_nocoi_$(exclude_coi)_onlyfull_$(only_full)_smoothw_$(smooth_w)_n_$(n_surrs)_nGS_$(n_surrs_per_gs)_p_$(pstring)_$(sur_method).jld2"
                end
                if !isfile(saveto)
                    println("Calculating $saveto")
                    if !isfile(saveto_temp)
                        obj= distribution_significant_increases(ice.name, qs_one, qs_two, 
                                    0, n_surrs_per_gs, p_thresholds, "sca", s1, s2, ice.resolution)
                    else 
                        obj = load(saveto_temp)["distribution"]
                    end
                    for (isu,xsur) in enumerate(siter)
                        if isu > obj.n
                            slopes, pvals_one, pvals_two = new_get_sca_slope_pvals(xsur,s1,s2,
                                    len_gs, max_len_gs, n_surrs_per_gs,
                                    dj, dt,
                                    mother, param,
                                    exclude_coi,
                                    method,
                                    filter_indicator, smooth_w, normalise, only_full)
                            for (ip,p) in enumerate(p_thresholds)
                                function q_nsig(slopes, pvals,p)
                                    real_slopes = findall(!isnothing, slopes) 
                                    pos_slopes = findall(>(0),slopes[real_slopes])
                                    good_pval = pvals[real_slopes][pos_slopes] .< p
                                    sum(good_pval)
                                end
                                obj.num_inc_one[ip,isu] = q_nsig(slopes,pvals_one,p)
                                obj.num_inc_two[ip,isu] = q_nsig(slopes,pvals_two,p)
                            end
                            #save every 20th iteration
                            if isu %20 == 0
                                println("calculated $(isu) surrogates so far.")
                                obj.n = isu
                                save(saveto_temp, "distribution", obj)
                            end
                        end
                    end

                    obj.n = n_surrs
                    save(saveto, "distribution", obj)
                    if isfile(saveto_temp)
                        rm(saveto_temp)
                    end
                else
                    println("$saveto already exists")
                end
            end
        end
    end
    println("Done :)")
end

function get_hurst_slope_pvals(x::Array{Float64},s1::Int64,s2::Int64,
    len_gs::Array{Int64}, max_len_gs::Int64, n_surrs_per_gs::Int64,
    dj::Float64, dt::Int64,
    mother::String, param::Int64,
    exclude_coi::Bool,
    method::Surrogate,
    filter_indicator::Bool, smooth_w::Bool, normalise::Bool, only_full::Bool)

    n_gs = length(len_gs)
    fake_gs_end_locs = rand(max_len_gs+1:length(x), 17)
    fake_gs_start_locs = fake_gs_end_locs .- len_gs .+1
    fake_gs_locs = [fake_gs_start_locs[i]: fake_gs_end_locs[i] for i in 1:17]

    slopes = Array{Union{Nothing,Float64}}(nothing, n_gs)
    pvals_one =  Array{Union{Nothing,Float64}}(nothing, n_gs)
    pvals_two =  Array{Union{Nothing,Float64}}(nothing, n_gs)
    for (ic,gs) in enumerate(fake_gs_locs)
        #calc ews (and significance at different pvals) here
        xx = x[gs]
        xx = convert.(Float64,xx)

        function q_hurst(xx)
            sca = get_h(xx,dt, s1,s2,mother, param, exclude_coi, smooth_w,
                filter_indicator, pad = 1, s0=-1, j1=-1, dj=0.1,normalise = normalise, only_full=only_full)
            slope = get_slopes(1:length(sca), sca)#, return_pred = false)
            slope
        end

        true_slope = q_hurst(xx)
        slopes[ic]= true_slope
        if (typeof(true_slope) != Nothing) #&& (true_slope > 0)
            #only for increases, if they exist
            test_tfts = TimeseriesSurrogates.SurrogateTest(q_hurst,xx,method, n=n_surrs_per_gs)
            pvals_two[ic] = TimeseriesSurrogates.pvalue(test_tfts, tail=:both)
            pvals_one[ic] = count(v -> v ≥ true_slope, test_tfts.vals)/n_surrs_per_gs                 
        end
    end
    return slopes, pvals_one, pvals_two
end


function sign_number_of_hurst_numeric(ice_core_list::Array{ice_core}, n_surrs::Int64, n_surrs_per_gs::Int64,  p_thresholds::Array{Float64}, sranges::Vector{Tuple{Int64, Int64}};
            mother::String = ["PAUL","MORLET"][1], param::Int64 = -1, 
            sur_method::String = "TFTS", 
            dat_type::String = ["orig", "normed", "filt", "normed_filt"][1], cold_type::String = ["gs", "gs_short"][1],
            filter_indicator::Bool = false, smooth_w::Bool = false,  exclude_coi::Bool = true, normalise::Bool = false, only_full::Bool = true,
            folder::String = "NGRIP5", gs_lengths = gs_lengths, gs_lengths_n = gs_lengths_n, batch_no::Int64 = 1)

    #1. make n_surrs surrogates of δ
    #2. place 17 GS randomly
    #3. calc ews (and significance at different pvals) there
    #4. get number of positive trends per time series
    #5. get distribution of that
    sur_method == uppercase(sur_method)
    if sur_method == "TFTS"
        method = TFTS(0.05)
    elseif sur_method == "TFTD"
        method = TFTDRandomFourier(true, 0.02)
    else
        method = RandomFourier(true)
    end

    filt_type = lowercase(dat_type)
    cold_type = lowercase(cold_type)

    for (i,ice) in enumerate(ice_core_list)
    if filt_type == "orig"
            data = ice.δ
        elseif filt_type == "normed"
            data = ice.δ_normed
        elseif filt_type == "filt"
            data = ice.δ_filt_100
        elseif filt_type == "normed_filt"
            data = ice.δ_normed_filt_100
        else
        error("ERROR: dat_type must be \"orig\", \"normed\", \"filt\" or \"normed_filt\".")
        end
        
        if cold_type == "gs"
            len_gs = gs_lengths[i]
            max_len_gs = maximum(len_gs)
        elseif cold_type == "gs_short"
            len_gs = gs_lengths_n[i]
            max_len_gs = maximum(len_gs)
        else
            error("cold_type must be \"gs\" or \"gs_short\".")
        end

        dj = 0.1
        dt = ice.resolution

        sgen = surrogenerator(data, method)
        siter = (sgen() for _ in 1:n_surrs)

        for (i,(s1,s2)) in enumerate(sranges)
        
            if s1 < s2
                qs_one = zeros(Int32,length(p_thresholds),n_surrs)
                qs_two = zeros(Int32,length(p_thresholds),n_surrs)
                pstring = string.(p_thresholds)
                pstring = join(pstring,"_")
                
                if batch_no>0
                    saveto = "new_surrogate_files/$(folder)/n_sig_dists/hurst/batch_$(batch_no)_s1_$(s1)_s2_$(s2)_$(ice.name)_$(mother)_$(filt_type)_$(cold_type)_FiltInd_$(filter_indicator)_normalise_$(normalise)_nocoi_$(exclude_coi)_onlyfull_$(only_full)_smoothw_$(smooth_w)_n_$(n_surrs)_nGS_$(n_surrs_per_gs)_p_$(pstring)_$(sur_method).jld2"
                    saveto_temp = "new_surrogate_files/$(folder)/n_sig_dists/hurst/TEMP_batch_$(batch_no)_s1_$(s1)_s2_$(s2)_$(ice.name)_$(mother)_$(filt_type)_$(cold_type)_FiltInd_$(filter_indicator)_normalise_$(normalise)_nocoi_$(exclude_coi)_onlyfull_$(only_full)_smoothw_$(smooth_w)_n_$(n_surrs)_nGS_$(n_surrs_per_gs)_p_$(pstring)_$(sur_method).jld2"
                else
                    saveto = "new_surrogate_files/$(folder)/n_sig_dists/hurst/s1_$(s1)_s2_$(s2)_$(ice.name)_$(mother)_$(filt_type)_$(cold_type)_FiltInd_$(filter_indicator)_normalise_$(normalise)_nocoi_$(exclude_coi)_onlyfull_$(only_full)_smoothw_$(smooth_w)_n_$(n_surrs)_nGS_$(n_surrs_per_gs)_p_$(pstring)_$(sur_method).jld2"
                    saveto_temp = "new_surrogate_files/$(folder)/n_sig_dists/hurst/TEMP_s1_$(s1)_s2_$(s2)_$(ice.name)_$(mother)_$(filt_type)_$(cold_type)_FiltInd_$(filter_indicator)_normalise_$(normalise)_nocoi_$(exclude_coi)_onlyfull_$(only_full)_smoothw_$(smooth_w)_n_$(n_surrs)_nGS_$(n_surrs_per_gs)_p_$(pstring)_$(sur_method).jld2"
                end
                
                if !isfile(saveto)
                    println("Calculating $saveto")
                    if !isfile(saveto_temp)
                        obj= distribution_significant_increases(ice.name, qs_one, qs_two, 
                                    0, n_surrs_per_gs, p_thresholds, "hurst", s1, s2, ice.resolution)
                    else 
                        obj = load(saveto_temp)["distribution"]
                    end
                    for (isu,xsur) in enumerate(siter)
                        if isu > obj.n
                            slopes, pvals_one, pvals_two = new_get_hurst_slope_pvals(xsur,s1,s2,
                                    len_gs, max_len_gs, n_surrs_per_gs,
                                    dj, dt,
                                    mother, param,
                                    exclude_coi,
                                    method,
                                    filter_indicator, smooth_w, normalise, only_full)
                            for (ip,p) in enumerate(p_thresholds)
                                obj.num_inc_one[ip,isu] =  count(v -> typeof(slopes[v])!= Nothing && slopes[v]>0 && pvals_one[v] <=p, axes(slopes, 1))
                                obj.num_inc_two[ip,isu] = count(v -> typeof(slopes[v])!= Nothing && slopes[v]>0 && pvals_two[v] <=p, axes(slopes, 1))
                            end
                            #save every 5th iteration
                            if isu %5 == 0
                                println("calculated $(isu) surrogates so far.")
                                obj.n = isu
                                save(saveto_temp, "distribution", obj)
                            end
                        end
                    end

                    obj.n = n_surrs
                    save(saveto, "distribution", obj)
                    if isfile(saveto_temp)
                        rm(saveto_temp)
                    end
                else
                    println("$saveto already exists")
                end
            end
        
        end
    end
    println("Done :)")
end


