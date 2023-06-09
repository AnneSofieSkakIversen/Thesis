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
JCC =Dict([])
HJCC = Dict([])
CVaR =Dict([])


for risktype in risktypes
    dict = Dict([])
    IS_name = "model_IS_$(inS_size)_"*risktype
    OS_name = "model_OS_$(inS_size)_"*risktype
    
### IS results
IS_results = load(OS_folder*IS_name*"_results.jld2")
OS_results = load(OS_folder*OS_name*"_results.jld2")
IS_obj_s = IS_results["obj_s_base"] #load(OS_folder*IS_name*"_obj_s.jld2")
IS_obj = IS_results["obj"]
OS_obj = OS_results["obj"]
pi_s = IS_results["pi_s"]
S = IS_results["S"]
SC = IS_results["SC"]
T = [i for i in 1:24]
C = data["F"]


    
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
dict["OS_obj"] = OS_obj
dict["bids"] = [FFR_acc FCRN_bid FCRD_bid mFRR_acc]
dict["acc"] = [FFR_acc FCRN_acc FCRD_acc mFRR_acc]
dict["act"] = [FFR_act FCRN_up_act FCRN_down_act FCRD_act mFRR_act]
dict["nd"] = [FFR_nd FCRN_up_nd FCRN_down_nd FCRD_nd mFRR_nd]
dict["grid_dependancy"] = sum(pi_s[s]*IS_results["p_grid_in"][t,s] for t in T,s in S)
dict["total_bid"] = FFR_acc + FCRN_bid + FCRD_bid + mFRR_acc
dict["total_acc"] = FFR_acc + FCRN_acc + FCRD_acc + mFRR_acc

(risktype == "no_risk" ? no_risk = copy(dict) : #
(risktype == "CC" ? CC = copy(dict) : #
(risktype == "JCC" ? JCC = copy(dict) : #
(risktype == "HJCC" ? HJCC = copy(dict) : #
(risktype == "CVaR" ? CVaR = copy(dict) : no_markets = copy(dict))))))

total_bid = FFR_acc + FCRN_bid + FCRD_bid + mFRR_acc
total_acc = FFR_acc + FCRN_acc + FCRD_acc + mFRR_acc
IS_time = IS_results["solve_time"]
OS_time = OS_results["solve_time"]
println("Risktype: "*risktype)
println("IS Objective: $(IS_obj)")
println("OS Objective: $(OS_obj)")
println("grid import: $(dict["grid_dependancy"])")
println("Total bid: $(total_bid)")
println("Total accepted: $(total_acc)")
println("IS Computation time: $(IS_time)")
println("OS Computation time: $(OS_time)")

end

function compare_performance(title_text,ylabel,no_risk,CC,JCC,HJCC,CVaR)
    p = groupedbar([transpose(no_risk) transpose(CC) transpose(JCC) transpose(HJCC) transpose(CVaR)],
    title = title_text,
    group=repeat(["FFR"; "FCRN"; "FCRD"; "mFRR"],outer=5),
    bar_position = :stack,
    xlabel = "Risk model",
    ylabel = ylabel,
    legend = :bottom,
    ylims = (0, 4),
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
    ylims = (0, 1),
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
    ylims = (200, 250),
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
    ylims = (300, 450),
    ytickfontsize=10,
    xtickfontsize=10,
    guidefontsize = 12,
    c = 6,
    label = "",
    xticks = ([1,2,3,4,5,6], ["No markets" "No risk" "CC" "JCC" "HJCC" "CVaR"])
    )
end

# Notes: savings are not happening because of bid-strategy, but in the second stage management of the Energy Community.
compare_performance("Bids","Bid size (kW)",no_risk["bids"],CC["bids"],JCC["bids"],HJCC["bids"],CVaR["bids"])
#png(Fig_folder*"Compare_TCL_reserve_bids.png")
compare_performance("Accepted Bids","Bid size (kW)",no_risk["acc"],CC["acc"],JCC["acc"],HJCC["acc"],CVaR["acc"])
#png(Fig_folder*"Compare_TCL_accepted_bids.png")
compare_performance_act("Activated Bids",no_risk["act"],CC["act"],JCC["act"],HJCC["act"],CVaR["act"])
#png(Fig_folder*"Compare_TCL_activated_bids.png")
compare_performance_act("Activation Delivery",no_risk["act"]-no_risk["nd"],CC["act"]-CC["nd"],JCC["act"]-JCC["nd"],HJCC["act"]-HJCC["nd"],CVaR["act"]-CVaR["nd"])
#png(Fig_folder*"Compare_TCL_delivered_bids.png")
#grid_dependancy("Grid Dependance","Expected grid import [kWh]",no_markets["grid_dependancy"],no_risk["grid_dependancy"],CC["grid_dependancy"],JCC["grid_dependancy"],HJCC["grid_dependancy"],CVaR["grid_dependancy"])
#png(Fig_folder*"Compare_grid_dependance.png")
#obj_values("Objective value","Expected cost [DKK]",no_markets["IS_obj"],no_risk["IS_obj"],CC["IS_obj"],JCC["IS_obj"],HJCC["IS_obj"],CVaR["IS_obj"])
#png(Fig_folder*"Compare_objective.png")


### Does OS comply with the 90% rule as interpreted?