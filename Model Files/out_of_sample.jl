## Out of sample testing
### Import Data ############


using LaTeXStrings
using JuMP, JLD2, Random
include("data_import_function_new.jl")
include("plot_output_SC.jl")
# Set datapath

S1 = [s for s in 1:1:364]# choose which Ancillary Market scenarios to use in scenario generation
S2 = [s for s in 1:1:364] # choose which Prosumer scenarios to use in scenario generation
l = length(S1)*length(S2)

datapath_new = "G:\\Mit drev\\DTU Sustainable Energy\\F23\\Thesis\\Modelling\\Input Data\\Raw Data\\" 
data_base = data_import_function_new(datapath_new,false) # Put to true to set all market prices to zero (except spot)
data = construct_scenarios(data_base,S1,S2)
save("const_data.jld2", data)
#data = load("const_data.jld2")


#####################################################
### Out of sample settings
##

Random.seed!(2349);
inS_size = 260
in_S = rand(1:l,inS_size)  # in sample scenarios
out_of_S = rand(setdiff(data["S"],in_S),5*inS_size) # out of sample scenarios - none in Sample
#risktype = "no_risk" # no_markets (remember 0 prices) "no_risk", "slack", "JCC", "HJCC", "CVaR" "CC"
risktypes = ["CVaR"] #["no_risk", "CC","JCC", "HJCC", "CVaR"] # "slack" "no_markets"["no_markets"] #

OS_folder = "OS_results\\"
saving = true

for risktype in risktypes
###################################################
#### In Sample Model ############################
println("\nINITIATING IS")

include("model_functions.jl")

in_S_model = Model(Gurobi.Optimizer)
Cooperative_Risk_model1(in_S_model,data,in_S,risktype)
in_S_results = solvem(in_S_model,in_S,risktype)


##
#####################################
### out of sample ###################


println("\nINITIATING OS")
out_of_S_model = Model(Gurobi.Optimizer)
Cooperative_Risk_model1(out_of_S_model,data,out_of_S,risktype)


T = data["T"]
C = data["C"]
F = data["F"]
E = data["E"]
S = [i for i in eachindex(out_of_S)] # keep!

@constraint(out_of_S_model, out_S_p_c_FCRN[c in C,t in T], out_of_S_model[:p_c_FCRN][c,t] == value(in_S_model[:p_c_FCRN][c,t]))
@constraint(out_of_S_model, out_S_p_c_mFRR[c in C,t in T], out_of_S_model[:p_c_mFRR][c,t] == value(in_S_model[:p_c_mFRR][c,t]))
@constraint(out_of_S_model, out_S_p_c_FFR[c in C,t in T],  out_of_S_model[:p_c_FFR][c,t] == value(in_S_model[:p_c_FFR][c,t]))
@constraint(out_of_S_model, out_S_p_c_FCRD[c in C,t in T], out_of_S_model[:p_c_FCRD][c,t] == value(in_S_model[:p_c_FCRD][c,t]))

### Replace bid apacity constraints with slack version
# Add slack variables
#@variable(out_of_S_model,slack_TCL_up[F,T,S] >= 0)
#@variable(out_of_S_model,slack_TCL_down[F,T,S] >= 0)
#@variable(out_of_S_model,slack_EV_up[E,T,S] >= 0)
#@variable(out_of_S_model,slack_EV_down[E,T,S] >= 0)
# Replace constraints
delete.(out_of_S_model, out_of_S_model[:TCL_up_reserve_limit][F,T,S])
delete.(out_of_S_model, out_of_S_model[:TCL_down_reserve_limit][F,T,S])
delete.(out_of_S_model, out_of_S_model[:EV_up_reserve_limit][E,T,S])
delete.(out_of_S_model, out_of_S_model[:EV_down_reserve_limit][E,T,S])
@constraint(out_of_S_model, out_of_S_TCL_up_reserve_limit[f in F,t in T,s in S], data["FCRN_accept"][f,t,out_of_S[s]]*out_of_S_model[:p_c_FCRN][f,t] + out_of_S_model[:p_c_mFRR][f,t] + data["FCRD_accept"][f,t,out_of_S[s]]*out_of_S_model[:p_c_FCRD][f,t] + out_of_S_model[:p_c_FFR][f,t] <= data["D_base_TCL"][f,t,out_of_S[s]] + out_of_S_model[:slack_TCL_up][f,t,s]) # Limit on ancillary upward reserve bid
@constraint(out_of_S_model, out_of_S_TCL_down_reserve_limit[f in F,t in T,s in S], data["FCRN_accept"][f,t,out_of_S[s]]*out_of_S_model[:p_c_FCRN][f,t] <= data["Pmax_TCL"] - data["D_base_TCL"][f,t,out_of_S[s]] + out_of_S_model[:slack_TCL_down][f,t,s]) # Limit on ancillary upward reserve bid
@constraint(out_of_S_model, out_of_S_EV_up_reserve_limit[e in E,t in T,s in S], data["FCRN_accept"][e,t,out_of_S[s]]*out_of_S_model[:p_c_FCRN][e,t] + out_of_S_model[:p_c_mFRR][e,t] + data["FCRD_accept"][e,t,out_of_S[s]]*out_of_S_model[:p_c_FCRD][e,t] + out_of_S_model[:p_c_FFR][e,t] <= data["EV_con"][e,t,out_of_S[s]]*data["Pmax_ev"][e] + out_of_S_model[:slack_EV_up][e,t,s]) # Limit on ancillary upward reserve bid
@constraint(out_of_S_model, out_of_S_EV_down_reserve_limit[e in E,t in T,s in S], data["FCRN_accept"][e,t,out_of_S[s]]*out_of_S_model[:p_c_FCRN][e,t] <= data["EV_con"][e,t,out_of_S[s]]*data["Pmax_ev"][e] + out_of_S_model[:slack_EV_down][e,t,s]) # Limit on ancillary upward reserve bid

is_valid.(out_of_S_model, out_of_S_model[:TCL_up_reserve_limit][F,T,S])
is_valid.(out_of_S_model, out_of_S_model[:out_of_S_TCL_up_reserve_limit][F,T,S])

### Objective function: Add slack penalty to objective function


OS_pen = 50 # penalty for violating constraints - no_markets no_risk diff = 84 dkk IS.
                # *1/3 because you get kicked out after 3 times you deliver less than 90%
                # + some, because you would earn more if you allow risk. = 50.

capacity_slack_penalty = @expression(out_of_S_model, OS_pen*sum(1/length(out_of_S)*(sum(out_of_S_model[:slack_TCL_up][f,t,s] + out_of_S_model[:slack_TCL_down][f,t,s] for f in F) + sum(out_of_S_model[:slack_EV_up][e,t,s] + out_of_S_model[:slack_EV_down][e,t,s] for e in E)) for s in S,t in T)); #keep S
old_objective = objective_function(out_of_S_model);
new_objective = @expression(out_of_S_model, old_objective + capacity_slack_penalty);
@objective(out_of_S_model, Min, new_objective);


### Solve out of sample problem
out_of_S_results = solvem(out_of_S_model,out_of_S,risktype);


#################################################################
######## IS / OS Scenario objectives #####################

### IS scenario objectives

        T = data["T"]
        C = data["C"]
        F = data["F"]
        E = data["E"]

    IS_obj_s = zeros(length(in_S))
    for s in eachindex(in_S)
    IS_obj_s[s] = value.(sum((data["grid_in_price"][t,in_S[s]]*in_S_model[:p_grid_in][t,s] - data["grid_out_price"][t,in_S[s]]*in_S_model[:p_grid_out][t,s]
    - sum(data["mFRR_price"][t,in_S[s]]*(in_S_model[:p_c_mFRR][c,t] - in_S_model[:q_c_mFRR][c,t,s])
            + data["FFR_price"][t,in_S[s]]*(in_S_model[:p_c_FFR][c,t]  - in_S_model[:q_c_FFR][c,t,s])
            + data["FCRD_price"][c]*data["FCRD_accept"][c,t,in_S[s]]*(in_S_model[:p_c_FCRD][c,t]  - in_S_model[:q_c_FCRD][c,t,s])
            + data["FCRN_price"][c]*data["FCRN_accept"][c,t,in_S[s]]*(in_S_model[:p_c_FCRN][c,t] - in_S_model[:q_c_FCRN_up][c,t,s] - in_S_model[:q_c_FCRN_down][c,t,s]) for c in C)
    - data["balance_price"][t,in_S[s]]*in_S_model[:p_FCRN_up][t,s] 
    - data["balance_price"][t,in_S[s]]*in_S_model[:p_FCRN_down][t,s] 
    - data["balance_price"][t,in_S[s]]*in_S_model[:p_mFRR_up][t,s]
        ) for t in T)
    + data["beta"]*sum((data["T_a"][t,in_S[s]] < data["T_opt"] ? (in_S_model[:T_TCL][f,t,s] - data["T_opt"]).^2 : 0 ) for f in F, t in T))
    end


    in_S_results["IS_obj_s"] = IS_obj_s

    IS_obj = objective_value(in_S_model)
    mean(IS_obj_s)

       
############################################
######## IS / OS Result analysis & Histograms #################

println("$(risktype)")
println("\nINITIATING HISTOGRAMS\n")
        OS_obj_s = zeros(length(out_of_S))
        for s in eachindex(out_of_S)
                OS_obj_s[s] = value.(sum((data["grid_in_price"][t,out_of_S[s]]*out_of_S_model[:p_grid_in][t,s] - data["grid_out_price"][t,out_of_S[s]]*out_of_S_model[:p_grid_out][t,s]
                        - sum(data["mFRR_price"][t,out_of_S[s]]*(out_of_S_model[:p_c_mFRR][c,t] - out_of_S_model[:q_c_mFRR][c,t,s])
                                + data["FFR_price"][t,out_of_S[s]]*(out_of_S_model[:p_c_FFR][c,t]  - out_of_S_model[:q_c_FFR][c,t,s])
                                + data["FCRD_price"][c]*data["FCRD_accept"][c,t,out_of_S[s]]*(out_of_S_model[:p_c_FCRD][c,t]  - out_of_S_model[:q_c_FCRD][c,t,s])
                                + data["FCRN_price"][c]*data["FCRN_accept"][c,t,out_of_S[s]]*(out_of_S_model[:p_c_FCRN][c,t] - out_of_S_model[:q_c_FCRN_up][c,t,s] - out_of_S_model[:q_c_FCRN_down][c,t,s]) for c in C)
                        - data["balance_price"][t,out_of_S[s]]*out_of_S_model[:p_FCRN_up][t,s] 
                        - data["balance_price"][t,out_of_S[s]]*out_of_S_model[:p_FCRN_down][t,s] 
                        - data["balance_price"][t,out_of_S[s]]*out_of_S_model[:p_mFRR_up][t,s]
                        ) for t in T)
                        + data["beta"]*sum((data["T_a"][t,out_of_S[s]] < data["T_opt"] ? (out_of_S_model[:T_TCL][f,t,s] - data["T_opt"]).^2 : 0) for f in F, t in T)
                        + OS_pen*sum(sum(out_of_S_model[:slack_TCL_up][f,t,s] + out_of_S_model[:slack_TCL_down][f,t,s] for f in F) 
                                + sum(out_of_S_model[:slack_EV_up][e,t,s] + out_of_S_model[:slack_EV_down][e,t,s] for e in E)#
                        for t in T));
        end
  

#######################################################
### OS Constraint violation Analysis

slack_TCL_up = zeros(length(C),length(T),length(S))
slack_TCL_up[F,T,S] = value.(Array(out_of_S_model[:slack_TCL_up]))#[F,T,out_of_S])
v_TCL_up = findall(!iszero, slack_TCL_up)

slack_TCL_down = zeros(length(C),length(T),length(S))
slack_TCL_down[F,T,S] = value.(Array(out_of_S_model[:slack_TCL_down]))#[F,T,out_of_S])
v_TCL_down = findall(!iszero, slack_TCL_down)

slack_EV_up = zeros(length(C),length(T),length(S))
slack_EV_up[E,T,S] = value.(Array(out_of_S_model[:slack_EV_up]))#[F,T,out_of_S])
v_EV_up = findall(!iszero, slack_EV_up)

slack_EV_down = zeros(length(C),length(T),length(S))
slack_EV_down[E,T,S] = value.(Array(out_of_S_model[:slack_EV_down]))#[F,T,out_of_S])
v_EV_down = findall(!iszero, slack_EV_down)

# Max violation value
max1 = maximum(Array(value.(out_of_S_model[:slack_TCL_up])))
max2 = maximum(value.(out_of_S_model[:slack_TCL_down]))
max3 = maximum(value.(out_of_S_model[:slack_EV_up]))
max4 = maximum(value.(out_of_S_model[:slack_EV_down]))

# Total violation
sum1 = sum(value.(out_of_S_model[:slack_TCL_up]))
sum2 = sum(value.(out_of_S_model[:slack_TCL_down]))
sum3 = sum(value.(out_of_S_model[:slack_EV_up]))
sum4 = sum(value.(out_of_S_model[:slack_EV_down]))

S_out_obj =[v_TCL_up ;v_TCL_down ; v_EV_up ;v_EV_down]
#obj_out = mean(out_of_S_results[:obj_s])

OS_violations = Dict(   "v_TCL_up" => v_TCL_up,
                        "v_TCL_down" => v_TCL_down,
                        "v_EV_up" => v_EV_up,
                        "v_EV_down" => v_EV_down,
                        "count" => Float64[size(v_TCL_up)[1] size(v_TCL_down)[1] size(v_EV_up)[1] size(v_EV_down)[1]],
                        "sum" => [sum1 sum2 sum3 sum4],
                        "max" =>[max1 max2 max3 max4])

# Save results
if (saving == true)
# IS
IS_name = "model_IS_$(inS_size)_"*risktype
save(OS_folder*IS_name*"_results.jld2", in_S_results)

# OS
OS_name = "model_OS_$(inS_size)_"*risktype
out_of_S_results["OS_obj_s"] = OS_obj_s # save obj_s in results 
save(OS_folder*OS_name*"_results.jld2", out_of_S_results)

# Violations
save(OS_folder*OS_name*"_violations.jld2", OS_violations)
end
println("Final: $(risktype)")

end

