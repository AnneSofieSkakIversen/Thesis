## Stability test testing
### Import Data ############
begin
using LaTeXStrings
using JuMP, JLD2, Random
include("data_import_function_new.jl")

# Set datapath


S1 = [s for s in 1:1:364]# choose which Ancillary Market scenarios to use in scenario generation
S2 = [s for s in 1:1:364] # choose which Prosumer scenarios to use in scenario generation
l = length(S1)*length(S2)

datapath_new = "G:\\Mit drev\\DTU Sustainable Energy\\F23\\Thesis\\Modelling\\Input Data\\Raw Data\\" 
data_base = data_import_function_new(datapath_new,true) # Put to true to set all market prices to zero (except spot)
data = construct_scenarios(data_base,S1,S2)
#data = load("const_data.jld2")
end

#####################################################
### Stability test settings

Random.seed!(2349);
stb_size = 40 # number of samples in model
no_runs = 100  # 100 # number of optimization runs
risktype = "no_markets" # no_markets (remember 0 prices) "no_risk", "slack", "JCC", "HJCC", "CVaR"
stb_folder = "STB_results\\"
saving = true #false
include("model_functions.jl")
###################################################
#### Stability test Models ############################


stb_obj = zeros(no_runs)
for i in 1:no_runs
        stb_S = rand(1:l,stb_size)  # in sample scenarios
        stb_model = Model(Gurobi.Optimizer)
        Cooperative_Risk_model1(stb_model,data,stb_S,risktype)
        optimize!(stb_model) #solvem(stb_model,stb_S,risktype) 
        stb_obj[i] = objective_value(stb_model)
end


#################################################################
######## Scenario objectives #####################

### STB scenario objectives
begin
        Plots.histogram(stb_obj, label = "$(no_runs) runs")
        Plots.vline!([mean(stb_obj)],linewidth = 4,color = "orange",label ="Mean obj. = $(round(mean(stb_obj),digits = 2))\nSt. dev. = $(round(std(stb_obj; corrected=false),digits = 2)) ")
        title!("Histogram, STB "*risktype)
        xlabel!("obj [DKK]")
        ylabel!("Frequency")
        
    if (saving == true) 
        png(stb_folder*risktype*"_$(no_runs)_histogram_Obj.png") 
        
        stb_results = Dict("stb_obj" => stb_obj)
        save(stb_folder*risktype*"_$(no_runs)_obj.jld2", stb_results) 
     end;
end

#######################################################
### Load data
 #   stb_results = load(stb_folder*risktype*"_$(no_runs)_obj.jld2")







#=
#JuMP.solution_summary(out_of_S_model)
#JuMP.lp_sensitivity_report(out_of_S_model)

# fortesting
#save("const_data.jld2", data)
#data = load("const_data.jld2")
#stb_S = rand(1:10,10)  # in sample scenarios
#out_of_S = rand(1:10,20)
#stb_name = "model_stb_260_"*risktype
#write_to_file(stb_model, stb_name*".nl")
#save(stb_name*".jld2", stb_results)
#stb_model = read_from_file("model_stb_260_"*risktype*".nl")
#stb_results = load(stb_name*".jld2")

#OS_name = "model_OS_260_"*risktype
#write_to_file(out_of_S_model, OS_name*".nl")
#save(OS_name*".jld2", out_of_S_results)
=#