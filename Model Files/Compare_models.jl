using LaTeXStrings
using JuMP, JLD2, Random

include("plot_output_SC.jl")
include("plot_output_SC_avg.jl")

risktypes = ["no_markets" "no_risk" "CC" "JCC" "HJCC" "CVaR"]
data = load("const_data.jld2")
OS_folder = "OS_results\\"
Fig_folder = "Figures\\"
inS_size = 260

no_markets =Dict([])
no_risk = Dict([])
CC = Dict([])
JCC = Dict([])
HJCC = Dict([])
CVaR = Dict([])


for risktype in risktypes
    dict = Dict([])
    T = [i for i in 1:24]
    
    
    
### IS results
IS_name = "model_IS_$(inS_size)_"*risktype
IS_results = load(OS_folder*IS_name*"_results.jld2")
IS_obj_s = IS_results["obj_s_base"]
IS_obj = IS_results["obj"]
pi_s = IS_results["pi_s"]
S = IS_results["S"]
SC = IS_results["SC"]
C = IS_results["C"]
F = data["F"]
E = data["E"]


### IS results
FCRN_bid = sum(IS_results["p_c_FCRN"][c,t] for c in C,t in T)
FCRD_bid = sum(IS_results["p_c_FCRD"][c,t] for c in C,t in T)
FFR_acc = sum(IS_results["p_c_FFR"][c,t] for c in C,t in T)
FCRN_acc = sum(pi_s[s]*IS_results["p_c_FCRN"][c,t]*data["FCRN_accept"][c,t,SC[s]] for c in C,t in T,s in S)
FCRD_acc = sum(pi_s[s]*IS_results["p_c_FCRD"][c,t]*data["FCRD_accept"][c,t,SC[s]] for c in C,t in T,s in S)
mFRR_acc = sum(IS_results["p_c_mFRR"][c,t] for c in C,t in T)
FFR_act = sum(pi_s[s]*IS_results["p_c_FFR"][c,t]*data["FFR_act"][t,SC[s]] for c in C,t in T,s in S)
FCRN_up_act = sum(pi_s[s]*IS_results["p_c_FCRN"][c,t]*data["FCRN_accept"][c,t,SC[s]]*data["FCRN_up_act"][t,SC[s]] for c in C,t in T,s in S)
FCRN_down_act = sum(pi_s[s]*IS_results["p_c_FCRN"][c,t]*data["FCRN_accept"][c,t,SC[s]]*data["FCRN_down_act"][t,SC[s]] for c in C,t in T,s in S)
FCRD_act = sum(pi_s[s]*IS_results["p_c_FCRD"][c,t]*data["FCRD_accept"][c,t,SC[s]]*data["FCRD_up_act"][t,SC[s]] for c in C,t in T,s in S)
mFRR_act = sum(pi_s[s]*IS_results["p_c_mFRR"][c,t]*data["mFRR_act"][t,SC[s]] for c in C,t in T,s in S)

FFR_nd = sum(pi_s[s]*IS_results["q_c_FFR"][c,t,s] for c in C,t in T,s in S)
FCRN_up_nd = sum(pi_s[s]*IS_results["q_c_FCRN_up"][c,t,s] for c in C,t in T,s in S)
FCRN_down_nd = sum(pi_s[s]*IS_results["q_c_FCRN_down"][c,t,s] for c in C,t in T,s in S)
FCRD_nd = sum(pi_s[s]*IS_results["q_c_FCRD"][c,t,s] for c in C,t in T,s in S)
mFRR_nd = sum(pi_s[s]*IS_results["q_c_mFRR"][c,t,s] for c in C,t in T,s in S)

dict["IS_obj"] = IS_obj
dict["IS_bids"] = [FFR_acc FCRN_bid FCRD_bid mFRR_acc]
dict["IS_acc"] = [FFR_acc FCRN_acc FCRD_acc mFRR_acc]
dict["IS_act"] = [FFR_act FCRN_up_act FCRN_down_act FCRD_act mFRR_act]
dict["IS_nd"] = [FFR_nd FCRN_up_nd FCRN_down_nd FCRD_nd mFRR_nd]
dict["IS_grid_dependancy"] = sum(pi_s[s]*IS_results["p_grid_in"][t,s] for t in T,s in S)
dict["IS_total_bid"] = FFR_acc + FCRN_bid + FCRD_bid + mFRR_acc
dict["IS_total_acc"] = FFR_acc + FCRN_acc + FCRD_acc + mFRR_acc

IS_total_bid = FFR_acc + FCRN_bid + FCRD_bid + mFRR_acc
IS_total_acc = FFR_acc + FCRN_acc + FCRD_acc + mFRR_acc
IS_time = IS_results["solve_time"]

### OS results
OS_name = "model_OS_$(inS_size)_"*risktype
OS_results = load(OS_folder*OS_name*"_results.jld2")
OS_obj_s = OS_results["obj_s_base"]
OS_obj = OS_results["obj"]
pi_s = OS_results["pi_s"]
S = OS_results["S"]
SC = OS_results["SC"]
C = OS_results["C"]

FCRN_bid = sum(OS_results["p_c_FCRN"][c,t] for c in C,t in T)
FCRD_bid = sum(OS_results["p_c_FCRD"][c,t] for c in C,t in T)
FFR_acc = sum(OS_results["p_c_FFR"][c,t] for c in C,t in T)
FCRN_acc = sum(pi_s[s]*OS_results["p_c_FCRN"][c,t]*data["FCRN_accept"][c,t,SC[s]] for c in C,t in T,s in S)
FCRD_acc = sum(pi_s[s]*OS_results["p_c_FCRD"][c,t]*data["FCRD_accept"][c,t,SC[s]] for c in C,t in T,s in S)
mFRR_acc = sum(OS_results["p_c_mFRR"][c,t] for c in C,t in T)
FFR_act = sum(pi_s[s]*OS_results["p_c_FFR"][c,t]*data["FFR_act"][t,SC[s]] for c in C,t in T,s in S)
FCRN_up_act = sum(pi_s[s]*OS_results["p_c_FCRN"][c,t]*data["FCRN_accept"][c,t,SC[s]]*data["FCRN_up_act"][t,SC[s]] for c in C,t in T,s in S)
FCRN_down_act = sum(pi_s[s]*OS_results["p_c_FCRN"][c,t]*data["FCRN_accept"][c,t,SC[s]]*data["FCRN_down_act"][t,SC[s]] for c in C,t in T,s in S)
FCRD_act = sum(pi_s[s]*OS_results["p_c_FCRD"][c,t]*data["FCRD_accept"][c,t,SC[s]]*data["FCRD_up_act"][t,SC[s]] for c in C,t in T,s in S)
mFRR_act = sum(pi_s[s]*OS_results["p_c_mFRR"][c,t]*data["mFRR_act"][t,SC[s]] for c in C,t in T,s in S)

FFR_nd = sum(pi_s[s]*OS_results["q_c_FFR"][c,t,s] for c in C,t in T,s in S)
FCRN_up_nd = sum(pi_s[s]*OS_results["q_c_FCRN_up"][c,t,s] for c in C,t in T,s in S)
FCRN_down_nd = sum(pi_s[s]*OS_results["q_c_FCRN_down"][c,t,s] for c in C,t in T,s in S)
FCRD_nd = sum(pi_s[s]*OS_results["q_c_FCRD"][c,t,s] for c in C,t in T,s in S)
mFRR_nd = sum(pi_s[s]*OS_results["q_c_mFRR"][c,t,s] for c in C,t in T,s in S)

dict["OS_obj"] = OS_obj
dict["OS_bids"] = [FFR_acc FCRN_bid FCRD_bid mFRR_acc]
dict["OS_acc"] = [FFR_acc FCRN_acc FCRD_acc mFRR_acc]
dict["OS_act"] = [FFR_act FCRN_up_act FCRN_down_act FCRD_act mFRR_act]
dict["OS_nd"] = [FFR_nd FCRN_up_nd FCRN_down_nd FCRD_nd mFRR_nd]
dict["OS_grid_dependancy"] = sum(pi_s[s]*OS_results["p_grid_in"][t,s] for t in T,s in S)
dict["OS_total_bid"] = FFR_acc + FCRN_bid + FCRD_bid + mFRR_acc
dict["OS_total_acc"] = FFR_acc + FCRN_acc + FCRD_acc + mFRR_acc


dict["OS_cap_viol_TCL_up"] = (data["FCRN_accept"][F,T,SC[S]].*OS_results["p_c_FCRN"][F,T] .+ OS_results["p_c_mFRR"][F,T] + data["FCRD_accept"][F,T,SC[S]].*OS_results["p_c_FCRD"][F,T] .+ OS_results["p_c_FFR"][F,T]) .> data["D_base_TCL"][F,T,SC[S]] # TCL up capacity limit
dict["OS_cap_viol_TCL_down"] = OS_results["p_c_FCRN"][F,T].> data["Pmax_TCL"] .- data["D_base_TCL"][F,T,SC[S]] # TCL down limit
dict["OS_cap_viol_EV_up"] = (data["FCRN_accept"][E,T,SC[S]].*OS_results["p_c_FCRN"][E,T] .+ OS_results["p_c_mFRR"][E,T] + data["FCRD_accept"][E,T,SC[S]].*OS_results["p_c_FCRD"][E,T] .+ OS_results["p_c_FFR"][E,T]) .> data["EV_con"][E,T,SC[S]].*data["Pmax_ev"][E] # Limit on ancillary upward reserve bid
dict["OS_cap_viol_EV_down"] = data["FCRN_accept"][E,T,SC[S]].*OS_results["p_c_FCRN"][E,T] .> data["EV_con"][E,T,SC[S]].*data["Pmax_ev"][E] # Limit on ancillary upward reserve bid

dict["OS_del_viol_TCL_up"] = [(data["FCRN_accept"][f,t,SC[s]]*OS_results["p_c_FCRN"][f,t]*data["FCRN_up_act"][t,SC[s]] + OS_results["p_c_mFRR"][f,t]*data["mFRR_act"][t,SC[s]] + data["FCRD_accept"][f,t,SC[s]]*OS_results["p_c_FCRD"][f,t]*data["FCRD_up_act"][t,SC[s]] + OS_results["p_c_FFR"][f,t]*data["FFR_act"][t,SC[s]]) > data["D_base_TCL"][f,t,SC[s]] for f in F, t in T,s in S] # TCL up capacity limit
dict["OS_del_viol_TCL_down"] = [OS_results["p_c_FCRN"][f,t]*data["FCRN_down_act"][t,SC[s]] > (data["Pmax_TCL"] - data["D_base_TCL"][f,t,SC[s]]) for f in F, t in T, s in S] # TCL down limit
dict["OS_del_viol_EV_up"] = [(data["FCRN_accept"][e,t,SC[s]]*OS_results["p_c_FCRN"][e,t]*data["FCRN_up_act"][t,SC[s]] + OS_results["p_c_mFRR"][e,t]*data["mFRR_act"][t,SC[s]] + data["FCRD_accept"][e,t,SC[s]]*OS_results["p_c_FCRD"][e,t]*data["FCRD_up_act"][t,SC[s]] + OS_results["p_c_FFR"][e,t]*data["FFR_act"][t,SC[s]]) > data["EV_con"][e,t,SC[s]]*data["Pmax_ev"][e] for e in E, t in T,s in S] # Limit on ancillary upward reserve bid
dict["OS_del_viol_EV_down"] = [data["FCRN_accept"][e,t,SC[s]]*OS_results["p_c_FCRN"][e,t]*data["FCRN_down_act"][t,SC[s]] > data["EV_con"][e,t,SC[s]]*data["Pmax_ev"][e] for e in E, t in T,s in S] # Limit on ancillary upward reserve bid

dict["OS_slack_bin_TCL_up"] = OS_results["slack_TCL_up"] .> 0
dict["OS_slack_bin_TCL_down"] = OS_results["slack_TCL_down"] .> 0
dict["OS_slack_bin_EV_up"] = OS_results["slack_EV_up"] .> 0
dict["OS_slack_bin_EV_down"] = OS_results["slack_EV_down"] .> 0

dict["OS_all_slacks"] = ([OS_results["slack_TCL_up"] + OS_results["slack_TCL_down"] ; OS_results["slack_EV_up"] + OS_results["slack_EV_down"]]) .> 0
sum(sum(dict["OS_all_slacks"][c,:,:] for c in C) .>0)/1300


OS_total_bid = FFR_acc + FCRN_bid + FCRD_bid + mFRR_acc
OS_total_acc = FFR_acc + FCRN_acc + FCRD_acc + mFRR_acc
OS_time = OS_results["solve_time"]

(risktype == "no_risk" ? no_risk = copy(dict) : #
(risktype == "CC" ? CC = copy(dict) : #
(risktype == "JCC" ? JCC = copy(dict) : #
(risktype == "HJCC" ? HJCC = copy(dict) : #
(risktype == "CVaR" ? CVaR = copy(dict) : no_markets = copy(dict))))))


println("IS "*risktype)
println("IS Riskmodel & obj & grid import & bid & accepted bids & activation delivery & CT ")
println(risktype*" & $(IS_obj) & $(dict["IS_grid_dependancy"]) & $(IS_total_bid) & $(IS_total_acc) & $(sum(dict["IS_act"]-dict["IS_nd"])) & $(IS_time) ")
println("")
println("OS "*risktype)
println("OS Riskmodel & obj & grid import & bid & accepted bids & activation delivery & CT ")
println(risktype*" & $(OS_obj) & $(dict["OS_grid_dependancy"]) & $(OS_total_bid) & $(OS_total_acc) & $(sum(dict["OS_act"]-dict["OS_nd"])) & $(OS_time) ")
println("")

end

function compare_performance(title_text,ylabel,no_risk,CC,JCC,HJCC,CVaR)
    p = groupedbar([transpose(no_risk) transpose(CC) transpose(JCC) transpose(HJCC) transpose(CVaR)],
    title = title_text,
    group=repeat(["FFR"; "FCRN"; "FCRD"; "mFRR"],outer=5),
    bar_position = :stack,
    xlabel = "Risk model",
    ylabel = ylabel,
    legend = :bottom,
    ylims = (0, 410),
    ytickfontsize=10,
    xtickfontsize=10,
    guidefontsize = 12,
    c = repeat([3;1;4;2],outer= 5),
    xticks = ([1,2,3,4,5], ["No risk" "CC" "JCC" "HJCC" "CVaR"])
    )
end
function compare_performance_act(title_text,no_risk,CC,JCC,HJCC,CVaR)
    p = groupedbar([transpose(no_risk) transpose(CC) transpose(JCC) transpose(HJCC) transpose(CVaR)],
    title = title_text,
    group=repeat(["FFR"; "FCRN up"; "FCRN down"; "FCRD"; "mFRR"],outer=5),
    bar_position = :stack,
    xlabel = "Risk model",
    ylabel = "Amount Activated bids (kWh)",
    legend = :top,
    ylims = (0, 100),
    ytickfontsize=10,
    xtickfontsize=10,
    guidefontsize = 12,
    c = repeat([3;1;"rgb(0,0,205)";4;2],outer= 5),
    xticks = ([1,2,3,4,5], ["No risk" "CC" "JCC" "HJCC" "CVaR"])
    )
end

function grid_dependancy(title_text,ylabel,no_markets,no_risk,CC,JCC,HJCC,CVaR)
    p = groupedbar(transpose([no_markets no_risk CC JCC HJCC CVaR]),
    title = title_text,
    #group=repeat(["FFR"; "FCRN"; "FCRD"; "mFRR"],outer=5),
    bar_position = :stack,
    xlabel = "Risk model",
    ylabel = ylabel,
    legend = :bottom,
    ylims = (200, 260),
    ytickfontsize=10,
    xtickfontsize=10,
    guidefontsize = 12,
    c = 5,
    label = "",
    xticks = ([1,2,3,4,5,6], ["No markets" "No risk" "CC" "JCC" "HJCC" "CVaR"])
    )
end

function obj_values(title_text,ylabel,no_markets,no_risk,CC,JCC,HJCC,CVaR)
    p = groupedbar(transpose([no_markets no_risk CC JCC HJCC CVaR]),
    title = title_text,
    #group=repeat(["FFR"; "FCRN"; "FCRD"; "mFRR"],outer=5),
    bar_position = :stack,
    xlabel = "Risk model",
    ylabel = ylabel,
    legend = :bottom,
    ylims = (400, 550),
    ytickfontsize=10,
    xtickfontsize=10,
    guidefontsize = 12,
    c = 6,
    label = "",
    xticks = ([1,2,3,4,5,6], ["No markets" "No risk" "CC" "JCC" "HJCC" "CVaR"])
    )
end

#Notes: savings are not happening because of bid-strategy, but in the second stage management of the Energy Community.
compare_performance("Bids","Bid size (kW)",no_risk["IS_bids"],CC["IS_bids"],JCC["IS_bids"],HJCC["IS_bids"],CVaR["IS_bids"])
png(Fig_folder*"Compare_reserve_bids.png")
compare_performance("Accepted Bids","Bid size (kW)",no_risk["IS_acc"],CC["IS_acc"],JCC["IS_acc"],HJCC["IS_acc"],CVaR["IS_acc"])
png(Fig_folder*"Compare_accepted_bids.png")
compare_performance_act("Activated Bids",no_risk["IS_act"],CC["IS_act"],JCC["IS_act"],HJCC["IS_act"],CVaR["IS_act"])
png(Fig_folder*"Compare_activated_bids.png")
compare_performance_act("Activation Delivery",no_risk["IS_act"]-no_risk["IS_nd"],CC["IS_act"]-CC["IS_nd"],JCC["IS_act"]-JCC["IS_nd"],HJCC["IS_act"]-HJCC["IS_nd"],CVaR["IS_act"]-CVaR["IS_nd"])
png(Fig_folder*"Compare_delivered_bids.png")
grid_dependancy("Grid Dependance","Expected grid import [kWh]",no_markets["IS_grid_dependancy"],no_risk["IS_grid_dependancy"],CC["IS_grid_dependancy"],JCC["IS_grid_dependancy"],HJCC["IS_grid_dependancy"],CVaR["IS_grid_dependancy"])
png(Fig_folder*"Compare_grid_dependance.png")
obj_values("Objective value","Expected cost [DKK]",no_markets["IS_obj"],no_risk["IS_obj"],CC["IS_obj"],JCC["IS_obj"],HJCC["IS_obj"],CVaR["IS_obj"])
png(Fig_folder*"Compare_objective.png")


### Does OS comply with the 90% rule as interpreted?
F = data["F"]
no_risk_results = load(OS_folder*"model_OS_260_no_risk_results.jld2")
no_risk_OS_cap_risk = sum(0 .< [sum.(no_risk_results["slack_TCL_up"][F,:,:] .+ no_risk_results["slack_TCL_down"][F,:,:]) ; sum.(no_risk_results["slack_EV_up"] + no_risk_results["slack_EV_down"])])/1300
no_risk_OS_del_risk = sum(sum(no_risk["OS_del_viol_TCL_up"].*no_risk["OS_slack_bin_TCL_up"]) + sum(no_risk["OS_del_viol_TCL_down"].*no_risk["OS_slack_bin_TCL_down"]) + sum(no_risk["OS_del_viol_EV_up"].*no_risk["OS_slack_bin_EV_up"]) + sum(no_risk["OS_del_viol_EV_down"].*no_risk["OS_slack_bin_EV_down"]))/1300


sum(no_risk["OS_slack_bin_TCL_up"])
sum(no_risk["OS_del_viol_TCL_up"].*no_risk["OS_slack_bin_TCL_up"])/1300
sum(no_risk["OS_cap_viol_TCL_up"].*no_risk["OS_slack_bin_TCL_up"])/1300
findall(!iszero,no_risk["OS_del_viol_TCL_up"].*no_risk["OS_slack_bin_TCL_up"])


# CC Rule interpretation
CC_results = load(OS_folder*"model_OS_260_CC_results.jld2")
CC_rule_FFR = sum(CC_results["y_FFR"])/1300
CC_rule_FCRN_up = sum(CC_results["y_FCRN_up"])/1300
CC_rule_FCRN_down = sum(CC_results["y_FCRN_down"])/1300
CC_rule_FCRD = sum(CC_results["y_FCRD"])/1300
CC_rule_mFRR = sum(CC_results["y_mFRR"])/1300


JCC_results = load(OS_folder*"model_OS_260_JCC_results.jld2")
SC = JCC_results["SC"]
JCC_rule_FFR = sum(1e-5 .< (data["FFR_act"][:,SC].*sum(eachrow(CC_results["p_c_FFR"]))))/1300
JCC_rule_FCRN_up = sum(1e-5 .< (data["FCRN_up_act"][:,SC].*sum(eachrow(CC_results["p_c_FCRN"]))))/1300
JCC_rule_FCRN_down = sum(1e-5 .< (data["FCRN_down_act"][:,SC].*sum(eachrow(CC_results["p_c_FCRN"]))))/1300
JCC_rule_FCRD = sum(1e-5 .< (data["FCRD_up_act"][:,SC].*sum(eachrow(CC_results["p_c_FCRD"]))))/1300
JCC_rule_mFRR = sum(1e-5 .< (data["mFRR_act"][:,SC].*sum(eachrow(CC_results["p_c_mFRR"]))))/1300

p_CC = groupedbar(transpose(1 .- [CC_rule_FFR CC_rule_FCRN_up CC_rule_FCRN_down CC_rule_FCRD CC_rule_mFRR]),
title = "CC model: Out of sample performance",
#group = ["FFR"; "FCRN up"; "FCRN down"; "FCRD"; "mFRR"],
bar_position = :stack,
xlabel = "Ancillary service",
ylabel = "Probability of full delivery",
#legend = :lower,
label = "",
ylims = (0, 1),
ytickfontsize=10,
xtickfontsize=10,
guidefontsize = 12,
c = repeat([3;1;"rgb(0,0,205)";4;2],outer= 5),
xticks = ([1,2,3,4,5], ["FFR"; "FCRN up"; "FCRN down"; "FCRD"; "mFRR"])
)
hline!([0.9],color = "black",linewidth = 2,label = "90% Rule: CC interpretation")
png(Fig_folder*"CC_OS_performance.png")

# JCC Rule Interpretation
JCC_rule_JCC = sum(load(OS_folder*"model_OS_260_JCC_results.jld2")["y"])/1300
JCC_rule_CC = sum(((CC_results["y_FFR"] .+ CC_results["y_FCRN_up"] .+ CC_results["y_FCRN_down"] .+ CC_results["y_FCRD"] .+ CC_results["y_mFRR"]) .>0))/1300
JCC_rule_HJCC = sum(sum(eachrow(load(OS_folder*"model_OS_260_HJCC_results.jld2")["y"])) .>0)/1300

CVaR_results = load("OS_results\\model_OS_260_CVaR_results.jld2")
C = data["C"]
T = [i for i in 1:24]
S = CVaR_results["S"]
pi_s = CVaR_results["pi_s"]
JCC_rule_CVaR = 1 - sum(1e-3 .< sum(CVaR_results["q_c_FFR"][c,t,s] .+ CVaR_results["q_c_FCRN_up"][c,t,s] .+ CVaR_results["q_c_FCRN_down"][c,t,s] .+ CVaR_results["q_c_FCRD"][c,t,s] .+ CVaR_results["q_c_mFRR"][c,t,s] for c in C,t in T) for s in S)/1300

p_JCC = groupedbar(transpose(1 .- [JCC_rule_CC JCC_rule_JCC JCC_rule_HJCC JCC_rule_CVaR]),
title = "Out of sample performance on JCC Rule",
bar_position = :stack,
xlabel = "Risk model",
ylabel = "Probability of full delivery in all markets",
ylims = (0, 1),
ytickfontsize=10,
xtickfontsize=10,
guidefontsize = 12,
labels = "",
c = [7;8;9;10],
xticks = ([1,2,3,4], ["JCC"; "CC"; "HJCC"; "CVaR"])
)
hline!([0.9],color = "black",linewidth = 2,label = "90% Rule: JCC Interpretation")
png(Fig_folder*"JCC_OS_performance.png")


# HJCC Rule interpretation
HJCC_rule_HJCC = 1 .- sum(eachcol(load(OS_folder*"model_OS_260_HJCC_results.jld2")["y"]))./1300

CVaR_results = load("OS_results\\model_OS_260_CVaR_results.jld2")
CC_results = load("OS_results\\model_OS_260_CC_results.jld2")
JCC_results = load("OS_results\\model_OS_260_JCC_results.jld2")

C = data["C"]
T = [i for i in 1:24]
S = CVaR_results["S"]
pi_s = CVaR_results["pi_s"]
HJCC_rule_CVaR = 1 .-[sum(1e-3 .< sum(CVaR_results["q_c_FFR"][c,t,s] .+ CVaR_results["q_c_FCRN_up"][c,t,s] .+ CVaR_results["q_c_FCRN_down"][c,t,s] .+ CVaR_results["q_c_FCRD"][c,t,s] .+ CVaR_results["q_c_mFRR"][c,t,s] for c in C) for s in S)/1300 for t in T]
HJCC_rule_CC = 1 .-[sum(1e-3 .< sum(CC_results["q_c_FFR"][c,t,s] .+ CC_results["q_c_FCRN_up"][c,t,s] .+ CC_results["q_c_FCRN_down"][c,t,s] .+ CC_results["q_c_FCRD"][c,t,s] .+ CC_results["q_c_mFRR"][c,t,s] for c in C) for s in S)/1300 for t in T]
HJCC_rule_JCC = 1 .-[sum(1e-3 .< sum(JCC_results["q_c_FFR"][c,t,s] .+ JCC_results["q_c_FCRN_up"][c,t,s] .+ JCC_results["q_c_FCRN_down"][c,t,s] .+ JCC_results["q_c_FCRD"][c,t,s] .+ JCC_results["q_c_mFRR"][c,t,s] for c in C) for s in S)/1300 for t in T]

# HJCC performance
p_HJCC = Plots.bar([i for i in 1:24], HJCC_rule_HJCC,
title = "HJCC Out of sample performance",
bar_position = :stack,
xlabel = "Hour",
ylabel = "Probability of full delivery",
ylims = (0, 1),
ytickfontsize=10,
xtickfontsize=10,
guidefontsize = 12,
labels = "Probability of full delivery on ancillary bid activation",
c = 6,
xticks = [6,12,18,24]
)
hline!([0.9],color = "black",linewidth = 2,label = "90% Rule: HJCC interpretation")

# All models: HJCC performance
p2_HJCC = Plots.plot([HJCC_rule_CC HJCC_rule_JCC HJCC_rule_HJCC HJCC_rule_CVaR],
title = "Out of sample performance: HJCC Rule",
#bar_position = :stack,
xlabel = "Hour",
ylabel = "Probability of full delivery",
labels =  ["CC" "JCC" "HJCC" "CVaR"],
ylims = (0.85, 1.02),
ytickfontsize=10,
xtickfontsize=10,
guidefontsize = 12,
linewidth = 2,
legend = :bottom,
#labels = "Probability of full delivery on ancillary bid activation",
c = [7 8 9 10],
xticks = [6,12,18,24]
)
hline!([0.9],color = "black",linewidth = 2,label = "90% Rule: HJCC interpretation")
png(Fig_folder*"HJCC_OS_performance.png")

CVaR_rule_CVaR = (1 .- CVaR["OS_nd"]./CVaR["OS_act"]).*(CVaR["OS_act"].> 1e-5)
CVaR_rule_CC = (1 .- CC["OS_nd"]./CC["OS_act"]).*(CC["OS_act"].> 1e-5)
CVaR_rule_JCC = (1 .- JCC["OS_nd"]./JCC["OS_act"]).*(JCC["OS_act"].> 1e-5)
CVaR_rule_HJCC = (1 .- HJCC["OS_nd"]./HJCC["OS_act"]).*(HJCC["OS_act"].> 1e-5)

# All models: CVaR performance
p_CVaR = groupedbar(transpose([CVaR_rule_CC ; CVaR_rule_JCC ; CVaR_rule_HJCC ; CVaR_rule_CVaR]),
title = "Out of sample performance: CVaR Rule",
bar_position = :dodge,
bar_width=0.7,
#xlabel = "Hour",
xticks = ([1,2,3,4,5], ["FFR"; "FCRN up"; "FCRN down"; "FCRD"; "mFRR"]),
ylabel = "% of activated ancillary amount delivered",
labels =  ["CC" "JCC" "HJCC" "CVaR"],
ylims = (0, 1.02),
ytickfontsize=10,
xtickfontsize=10,
guidefontsize = 12,
linewidth = 2,
legend = :bottom,
#labels = "Probability of full delivery on ancillary bid activation",
c = [7 8 9 10],
#xticks = [6,12,18,24]
)
hline!([0.9],color = "black",linewidth = 2,label = "90% Rule: CVaR interpretation")
png(Fig_folder*"CVaR_OS_performance.png")


# All models: CC performance
CC_rule_CVaR = (1 .- CVaR["OS_nd"]./CVaR["OS_act"]).*(CVaR["OS_act"].> 1e-5)
CC_rule_CC = (1 .- CC["OS_nd"]./CC["OS_act"]).*(CC["OS_act"].> 1e-5)
CC_rule_JCC = (1 .- JCC["OS_nd"]./JCC["OS_act"]).*(JCC["OS_act"].> 1e-5)
CC_rule_HJCC = (1 .- HJCC["OS_nd"]./HJCC["OS_act"]).*(HJCC["OS_act"].> 1e-5)

# All models: CVaR performance
p_CVaR = groupedbar(transpose([CVaR_rule_CC ; CVaR_rule_JCC ; CVaR_rule_HJCC ; CVaR_rule_CVaR]),
title = "Out of sample performance: CVaR Rule",
bar_position = :dodge,
bar_width=0.7,
#xlabel = "Hour",
xticks = ([1,2,3,4,5], ["FFR"; "FCRN up"; "FCRN down"; "FCRD"; "mFRR"]),
ylabel = "% of activated ancillary amount delivered",
labels =  ["CC" "JCC" "HJCC" "CVaR"],
ylims = (0, 1.02),
ytickfontsize=10,
xtickfontsize=10,
guidefontsize = 12,
linewidth = 2,
legend = :bottom,
#labels = "Probability of full delivery on ancillary bid activation",
c = [7 8 9 10],
#xticks = [6,12,18,24]
)
hline!([0.9],color = "black",linewidth = 2,label = "90% Rule: CVaR interpretation")
png(Fig_folder*"CVaR_OS_performance.png")




# Testing OS results on eachothers 90% criteria
JCC_rule_CC = sum(((CC_results["y_FFR"] + CC_results["y_FCRN_up"] + CC_results["y_FCRN_down"] + CC_results["y_FCRD"] + CC_results["y_mFRR"]) .>0))/1300








### OS Violations
sum(no_risk["OS_slack_bin_TCL_up"])
sum(no_risk["OS_del_viol_TCL_up"].*no_risk["OS_slack_bin_TCL_up"])
sum(no_risk["OS_cap_viol_TCL_up"].*no_risk["OS_slack_bin_TCL_up"])
findall(!iszero,no_risk["OS_del_viol_TCL_up"].*no_risk["OS_slack_bin_TCL_up"])

violations = load(OS_folder*OS_name*"_violations.jld2")
v1 = violations["v_TCL_up"]
v2 = violations["v_TCL_down"]
v3 = violations["v_EV_up"]
v4 = violations["v_EV_down"]

##### Binding Constraints and shadow prices ########

#no_markets =Dict([])
s_no_risk = Dict([])
s_CC = Dict([])
s_JCC = Dict([])
s_HJCC = Dict([])
s_CVaR = Dict([])

risktypes2 = ["no_risk" "CC" "JCC" "HJCC" "CVaR"]
for risktype in risktypes2
IS_name = "model_IS_260_"*risktype
    if (risktype == "no_risk")
        s_no_risk = load("Sensitivity_test\\"*IS_name*"_duals_bind.jld2")
        elseif (risktype == "CC")
        s_CC = load("Sensitivity_test\\"*IS_name*"_duals_bind.jld2")
        elseif (risktype == "JCC")
        s_JCC = load("Sensitivity_test\\"*IS_name*"_duals_bind.jld2")
        elseif (risktype == "HJCC")
        s_HJCC = load("Sensitivity_test\\"*IS_name*"_duals_bind.jld2")
        elseif (risktype == "CVaR")
        s_CVaR = load("Sensitivity_test\\"*IS_name*"_duals_bind.jld2")
        else println("risktype not identified")
    end
end

no_risk_bind = length(s_no_risk["dual_list1"]["dual"]) + length(s_no_risk["dual_list2"]["dual"])
CC_bind = length(s_CC["dual_list1"]["dual"]) + length(s_CC["dual_list2"]["dual"])
JCC_bind = length(s_JCC["dual_list1"]["dual"]) + length(s_JCC["dual_list2"]["dual"])
HJCC_bind = length(s_HJCC["dual_list1"]["dual"]) + length(s_HJCC["dual_list2"]["dual"])
CVaR_bind = length(s_CVaR["dual_list1"]["dual"]) + length(s_CVaR["dual_list2"]["dual"])

no_risk_bind1 = [s_no_risk["dual_list1"]["con_eq"] s_no_risk["dual_list1"]["dual"]]
no_risk_bind2 = [s_no_risk["dual_list2"]["con_eq"] s_no_risk["dual_list2"]["dual"]]

CC_bind1 = s_CC["dual_list1"]["con_eq"]
CC_bind2 = s_CC["dual_list2"]["con_eq"]

JCC_bind1 = s_JCC["dual_list1"]["con_eq"]
JCC_bind2 = s_JCC["dual_list2"]["con_eq"]

HJCC_bind1 = s_HJCC["dual_list1"]["con_eq"]
HJCC_bind2 = s_HJCC["dual_list2"]["con_eq"]

CVaR_bind1 = s_CVaR["dual_list1"]["con_eq"]
CVaR_bind2 = s_CVaR["dual_list2"]["con_eq"]
CVaR_bind2 = s_CVaR["dual_list2"]["dual"]
