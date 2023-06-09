using LaTeXStrings
using JuMP, JLD2, Random, Gurobi,Statistics
import DataFrames

include("model_functions.jl")
#include("plot_output_SC.jl")
#include("plot_output_SC_avg.jl")

risktypes = ["JCC"] # "CVaR" "no_markets" "no_risk" "CC" "JCC" "HJCC"

#risktype = "CC" # no_markets "no_risk", "slack", "JCC", "HJCC", "CVaR"
data = load("const_data.jld2")
OS_folder = "OS_results\\"
Fig_folder = "Figures\\"
inS_size = 260

for risktype in risktypes 
IS_name = "model_IS_$(inS_size)_"*risktype


### IS results
IS_results = load(OS_folder*IS_name*"_results.jld2")
IS_obj_s = IS_results["obj_s_base"] #load(OS_folder*IS_name*"_obj_s.jld2")
IS_obj = IS_results["obj"]
mean(IS_obj_s)

T = data["T"]
C = data["C"]
F = data["F"]
E = data["E"]
S = IS_results["S"]
SC = IS_results["SC"]


println("\nINITIATING fix bin model")
fix_bin_model = Model(Gurobi.Optimizer)
Cooperative_Risk_model1(fix_bin_model,data,SC,"slack")

delete.(fix_bin_model, fix_bin_model[:y][S])
unregister.(fix_bin_model, :y)

if (risktype == "HJCC")
    println("model HJCC = $(risktype)")
    
    Pmax_EC = data["Pmax_EC"]
    @constraint(fix_bin_model, HJCC_mFRR[t in T, s in S], fix_bin_model[:q_mFRR][t,s] <= IS_results["y"][t,s]*Pmax_EC)
    @constraint(fix_bin_model, HJCC_FFR[t in T, s in S], fix_bin_model[:q_FFR][t,s] <= IS_results["y"][t,s]*Pmax_EC)
    @constraint(fix_bin_model, HJCC_FCRD[t in T,s in S], fix_bin_model[:q_FCRD][t,s] <= IS_results["y"][t,s]*Pmax_EC)
    @constraint(fix_bin_model, HJCC_FCRN_up[t in T, s in S], fix_bin_model[:q_FCRN_up][t,s] <= IS_results["y"][t,s]*Pmax_EC)
    @constraint(fix_bin_model, HJCC_FCRN_down[t in T, s in S], fix_bin_model[:q_FCRN_down][t,s] <= IS_results["y"][t,s]*Pmax_EC)
    
elseif (risktype == "JCC")
    println("model JCC = $(risktype)")
    
    Pmax_EC = data["Pmax_EC"]
    @constraint(fix_bin_model, JCC_mFRR[s in S], sum(fix_bin_model[:q_mFRR][t,s] for t in T) <= IS_results["y"][s]*24*Pmax_EC)
    @constraint(fix_bin_model, JCC_FFR[s in S], sum(fix_bin_model[:q_FFR][t,s] for t in T) <= IS_results["y"][s]*24*Pmax_EC)
    @constraint(fix_bin_model, JCC_FCRD[s in S], sum(fix_bin_model[:q_FCRD][t,s] for t in T) <= IS_results["y"][s]*24*Pmax_EC)
    @constraint(fix_bin_model, JCC_FCRN_up[s in S], sum(fix_bin_model[:q_FCRN_up][t,s] for t in T) <= IS_results["y"][s]*24*Pmax_EC)
    @constraint(fix_bin_model, JCC_FCRN_down[s in S], sum(fix_bin_model[:q_FCRN_down][t,s] for t in T) <= IS_results["y"][s]*24*Pmax_EC)

elseif (risktype == "CC")
    println("model CC = $(risktype)")
    
    Pmax_EC = data["Pmax_EC"]
    @constraint(fix_bin_model, CC_FFR[s in S], sum(fix_bin_model[:q_FFR][t,s] for t in T) <= IS_results["y_FFR"][s]*24*Pmax_EC)
    @constraint(fix_bin_model, CC_FCRN_up[s in S], sum(fix_bin_model[:q_FCRN_up][t,s] for t in T) <= IS_results["y_FCRN_up"][s]*24*Pmax_EC)
    @constraint(fix_bin_model, CC_FCRN_down[s in S], sum(fix_bin_model[:q_FCRN_down][t,s] for t in T) <= IS_results["y_FCRN_down"][s]*24*Pmax_EC)
    @constraint(fix_bin_model, CC_FCRD[s in S], sum(fix_bin_model[:q_FCRD][t,s] for t in T) <= IS_results["y_FCRD"][s]*24*Pmax_EC)
    @constraint(fix_bin_model, CC_mFRR[s in S], sum(fix_bin_model[:q_mFRR][t,s] for t in T) <= IS_results["y_mFRR"][s]*24*Pmax_EC)
else 
    println("no_risk or CVaR model = $(risktype)")
    fix_bin_model = Model(Gurobi.Optimizer)
    Cooperative_Risk_model1(fix_bin_model,data,SC,risktype)

    delete.(fix_bin_model, fix_bin_model[:y][S])
    unregister.(fix_bin_model, :y)
end


# Solve model with fixed binaries.
optimize!(fix_bin_model)
#fix_bin_results = solvem(fix_bin_model,SC,risktype)

if (dual_status(fix_bin_model) == 0)
    println("Dual solution not found for: $(risktype)")
    break;
else

con_name0 = all_constraints(fix_bin_model, AffExpr, MOI.EqualTo{Float64})
con_name1 = all_constraints(fix_bin_model, AffExpr, MOI.GreaterThan{Float64})
con_name2 = all_constraints(fix_bin_model, AffExpr, MOI.LessThan{Float64})
con_name = [con_name1 ; con_name2]

bind_con1 = findall(dual.(con_name1).> 1e-5)
bind_con2 = findall(dual.(con_name2).< -1e-5)

con_name1[bind_con1]
con_name2[bind_con2]

p1 = sortperm(dual.(con_name1),rev=true)
p2 = sortperm(dual.(con_name2))
dual_list1 = Dict("con_nr" => p1[1:length(bind_con1)+1], "con_eq" => SubString.(string.(con_name1[p1[1:length(bind_con1)+1]]),1,35), "dual" => dual.(con_name1[p1[1:length(bind_con1)+1]]))
dual_list2 = Dict("con_nr" => p2[1:length(bind_con2)+1], "con_eq" => SubString.(string.(con_name2[p2[1:length(bind_con2)+1]]),1,35), "dual" => dual.(con_name2[p2[1:length(bind_con2)+1]]))

duals_bind = Dict("obj" => objective_value(fix_bin_model),"dual_list1" => dual_list1, "dual_list2" => dual_list2)
save("Sensitivity_test\\"*IS_name*"_duals_bind.jld2", duals_bind)
end
end








#####Load results ########

risktype = "CVaR"
inS_size = 260

IS_name = "model_IS_$(inS_size)_"*risktype
df = load("Sensitivity_test\\"*IS_name*"_duals_bind.jld2")

length(df["dual_list1"]["con_eq"])
length(df["dual_list2"]["con_eq"])

println("Binding constraints: "*risktype*", positive duals")
for i in 1:min(800,length(df["dual_list1"]["con_eq"]))
println(df["dual_list1"]["con_eq"][i]*"     :    $(df["dual_list1"]["dual"][i])")
end

println("Binding constraints: "*risktype*", negative duals")
for i in 1001:min(length(df["dual_list2"]["con_eq"]),2000)
println(df["dual_list2"]["con_eq"][i]*"     :    $(df["dual_list2"]["dual"][i])")
end

#solution_summary(fix_bin_model; verbose = true)
#lp_sensitivity_report(fix_bin_model)

#list_of_constraint_types(fix_bin_model)
#(AffExpr, MOI.EqualTo{Float64})
#(AffExpr, MOI.GreaterThan{Float64})
#(AffExpr, MOI.LessThan{Float64})
#(VariableRef, MOI.EqualTo{Float64})
#(VariableRef, MOI.GreaterThan{Float64})

