using Pkg
Pkg.activate("/SEC/ZebraFish/Tools_Analysis/ZebrafishHMM_EntropyPaper/Julia_code/")
using HiddenMarkovModels: viterbi
using HDF5
using ZebrafishHMM2023: load_hmm, HMM_ARTR_Log
# using Statistics: mean, std, var
# using ZebrafishHMM2023: load_artr_wolf_2023, HMM_ARTR_Log, normalize_transition_matrix, ATol, viterbi_artr,
#     artr_wolf_2023_temperatures, artr_wolf_2023_fishes, split_into_repeated_subsequences, load_hmm, artr_wolf_2023


# loading the data
datapath = "Data/neuro.h5"
@assert isfile(datapath)
data = h5open(datapath, "r")
function load_data(fish::Int, temp::Int)
  grp = data["T$(temp)-fish$(fish)"]
  L = read(grp["L"])
  R = read(grp["R"])
  return L, R
end
all_ARTRs = read(data["fish_temp_combinations"])

# load the models
modelpath = "Models/hmms_ARTR_20240620/"
@assert isdir(modelpath)
function load_model(fish::Int, temp::Int)
  path = joinpath(modelpath, "artr_temperature=$(temp)-fish=$(fish).hdf5")
  model = load_hmm(path, HMM_ARTR_Log)
  return model
end

# output
outpath = "Models/hmms_ARTR_outputs.h5"
out = h5open(outpath, "w")


for i in 1:size(all_ARTRs, 2)
  fish, temp = all_ARTRs[:, i]
  println("Fish: $fish, Temp: $temp")
  L, R = load_data(fish, temp)
  spikes = collect(eachcol(vcat(L, R)))
  hmm = load_model(fish, temp)
  state_seq = viterbi(hmm, spikes)
  gen_states, gen_spikes = rand(hmm, length(spikes))
  gen_spikes = permutedims(mapreduce(permutedims, vcat, gen_spikes))
  # println("size of gen_states: $(size(gen_states)), size of gen_spikes: $(size(gen_spikes)) ; $(size(L)),$(size(R))")
  # println(typeof(gen_spikes))
  gen_L = gen_spikes[1:size(L, 1), :]
  gen_R = gen_spikes[size(L, 1)+1:end, :]
  # println("size of gen_L: $(size(gen_L)), size of gen_R: $(size(gen_R))")

  #saving the output
  grp = create_group(out, "T$(temp)-fish$(fish)")
  write(grp, "viterbi", state_seq .- 1)
  write(grp, "gen_states", gen_states .- 1)
  write(grp, "gen_L", UInt8.(gen_L))
  write(grp, "gen_R", UInt8.(gen_R))
end















close(out)
close(data)

