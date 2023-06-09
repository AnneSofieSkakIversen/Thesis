
    using LaTeXStrings
    using JuMP, JLD2, Random
    using Plots, Statistics
    data = load("const_data.jld2")

    include("plot_output_SC.jl")
    include("plot_output_SC_avg.jl")
begin
risktype = "JCC" # no_markets (remember 0 prices) "no_risk", "slack", "JCC", "HJCC", "CVaR"
inS_size = 260
OS_folder = "OS_results\\"
OS_name = "model_OS_$(inS_size)_"*risktype
violations = load(OS_folder*OS_name*"_violations.jld2")
results = load(OS_folder*OS_name*"_results.jld2")
saving = false #true

V_folder = "OS_violations\\"

v1 = violations["v_TCL_up"]
v2 = violations["v_TCL_down"]

v3 = violations["v_EV_up"]
v4 = violations["v_EV_down"]
v = unique!([v1; v2 ;v3; v4])

sum = violations["sum"]



SC = results["SC"]
S = results["S"]
F = data["F"]
C = data["C"]
E = data["E"]
T = [i for i in 1:24]
Pmax_ev = [0 0 0 0 10 10 20 7 7 20] #charging power [kW] of the battery. Same models as above. (C), AC charger, one DC charger high each
    
end
data["FCRN_accept"][c,t,s]*data["FCRN_up_act"][t,s]*results["p_c_FCRN"][c,t] + data["FCRD_accept"][c,t,s]*data["FCRD_up_act"][t,s]*results["p_c_FCRD"][c,t] + results["p_c_mFRR"][c,t].*data["mFRR_act"][t,s] + results["p_c_FFR"][c,t].*data["FFR_act"][t,s]


a = findall([data["D_base_TCL"][c,t,SC[s]] .< data["FCRN_accept"][c,t,SC[s]]*data["FCRN_up_act"][t,SC[s]]*results["p_c_FCRN"][c,t] + data["FCRD_accept"][c,t,SC[s]]*data["FCRD_up_act"][t,SC[s]]*results["p_c_FCRD"][c,t] + results["p_c_mFRR"][c,t].*data["mFRR_act"][t,SC[s]] + results["p_c_FFR"][c,t].*data["FFR_act"][t,SC[s]] for c in F,t in T,s in S])
b = findall([data["Pmax_TCL"] .- data["D_base_TCL"][f,t,SC[s]] .< data["FCRN_accept"][f,t,SC[s]]*results["p_c_FCRN"][f,t]*data["FCRN_down_act"][t,SC[s]] for f in F,t in T,s in S])

c = findall([data["EV_con"][c,t,SC[s]]*Pmax_ev[c] .< data["FCRN_accept"][c,t,SC[s]]*data["FCRN_up_act"][t,SC[s]]*results["p_c_FCRN"][c,t] + data["FCRD_accept"][c,t,SC[s]]*data["FCRD_up_act"][t,SC[s]]*results["p_c_FCRD"][c,t] + results["p_c_mFRR"][c,t].*data["mFRR_act"][t,SC[s]] + results["p_c_FFR"][c,t].*data["FFR_act"][t,SC[s]] for c in E,t in T,s in S])
d = findall([data["EV_con"][c,t,SC[s]]*Pmax_ev[c] .< data["FCRN_accept"][c,t,SC[s]]*results["p_c_FCRN"][c,t]*data["FCRN_down_act"][t,SC[s]] for c in E,t in T,s in S])

length(a) + length(b) + length(c) + length(d)
#                                            p_c_FCRN_down[c,t,s] + q_c_FCRN_down[c,t,s] == FCRN_accept[c,t,s]*FCRN_down_act[t,s]*p_c_FCRN[c,t]
# p_c_FCRD_up[c,t,s] + q_c_FCRD[c,t,s] == FCRD_accept[c,t,s]*FCRD_up_act[t,s]*p_c_FCRD[c,t]
#p_c_mFRR_up[c,t,s] + q_c_mFRR[c,t,s] == p_c_mFRR[c,t].*mFRR_act[t,s]
# p_c_FFR_up[c,t,s] + q_c_FFR[c,t,s] == p_c_FFR[c,t].*FFR_act[t,s]


v_TCL = unique!([v1; v2])
v_EV = unique!([v3; v4])

v_t = [v[i][2] for i in eachindex(v)]
v_sc = [v[i][3] for i in eachindex(v)] # scenarios with

v_TCL_c = [v[i][1] for i in eachindex(v_TCL)]
v_TCL_t = [v[i][2] for i in eachindex(v_TCL)]
v_TCL_sc = [v[i][3] for i in eachindex(v_TCL)]
v_TCL_sc_u = unique!(v_TCL_sc)

v_EV_c = [v[i][1] for i in eachindex(v_EV)]
v_EV_t = [v[i][2] for i in eachindex(v_EV)]
v_EV_sc = [v[i][3] for i in eachindex(v_EV)]
v_EV_sc_u = unique!(v_EV_sc)


# Histogram without violations
#v1_sc = unique!([v1[i][3] for i in eachindex(v1)] )
#v2_sc = unique!([v2[i][3] for i in eachindex(v2)] )
#v3_sc = unique!([v3[i][3] for i in eachindex(v3)] )
#v4_sc = unique!([v4[i][3] for i in eachindex(v4)] )
#v_s = unique!([v1_sc ; v2_sc ; v3_sc ;v4_sc]) 
v_s = unique!([v[i][3] for i in eachindex(v)]) # unique scenarios with violation



S = [i for i in eachindex(results["SC"])]
obj_s = results["OS_obj_s"]
v_obj_s = results["OS_obj_s"][v_s]
nv_obj_s = results["OS_obj_s"][setdiff(S,v_s)]



Plots.histogram(obj_s, ylims = (0,4),label = "All OS Scenario objectives")
#Plots.histogram!(nv_obj_s, label = "Not Violating OS Samples")
Plots.histogram!(v_obj_s,   bins = 50,label = "$(length(v_s)) violating OS scenarios")
Plots.vline!([mean(obj_s)],linewidth = 4,color = "blue",linestyle = :dash,label ="Mean OS obj = $(round(mean(obj_s),digits = 2))")#\nStd. dev. = $(round(std(obj_s; corrected=false),digits = 2))")
#Plots.vline!([mean(nv_obj_s)],linewidth = 4,color = "blue",linestyle = :dash,label ="E(nv_obj) = $(round(mean(nv_obj_s),digits = 2))")#\nStd. dev. = $(round(std(nv_obj_s; corrected=false),digits = 2))")
Plots.vline!([mean(v_obj_s)],linewidth = 4,color = 2, linestyle = :dash,label ="Mean viol. obj = $(round(mean(v_obj_s),digits = 2))")#\nStd. dev. = $(round(std(v_obj_s; corrected=false),digits = 2))")
title!("Histogram: violating OS scenarios, "*risktype)
xlabel!("OS obj [DKK]")
ylabel!("Frequency")
if (saving == true) png(V_folder*risktype*"_histogram_Viol_OS_Obj.png") end

####################################
#### Analysis of violating scenarios
v_results = load(OS_folder*OS_name*"_results.jld2")
v_results["pi_s"] .*= 0
v_results["pi_s"][v_s] .+= 1/length(v_s)
v_results["S"] = v_s # including only violating scenarios
SC = v_results["SC"]
SC[v_s]
C = data["C"]
T = data["T"]

function plot_v_histogram(data_name,data,v_data,bins1,bins2,ymax)
 
Plots.histogram(data, ylims = (0,ymax), bins = bins1, label = "Non-violating OS scenarios")
Plots.histogram!(v_data, bins = bins2, label = "violating OS scenarios")
Plots.vline!([mean(data)],linewidth = 3,color = "blue", linestyle = :dash, label ="Mean non-violating scenarios = $(round(mean(data),digits = 2))")
Plots.vline!([mean(v_data)],linewidth = 3,color = 2, linestyle = :dash, label ="Mean violating scenarios = $(round(mean(v_data),digits = 2))")

title!("OS Violation Histogram: "*data_name)
xlabel!(data_name)
ylabel!("Frequency")
end

#plot_v_histogram("OS_obj_s [DKK]",obj_s,v_obj_s,100,50,10)
plot_v_histogram("PV [kW]",mean.(eachcol(data["PV"][:,SC[S]])),mean.(eachcol(data["PV"][:, SC[v_TCL_sc]])),50,5,100) # bin bin yscale
plot_v_histogram("PV [kW]",mean.(eachcol(data["PV"][:,SC[S]])),mean.(eachcol(data["PV"][:, SC[v_EV_sc]])),1,1,10)


plot_v_histogram("D_household [kWh]",mean.(data["D_household"][:,:,s] for s in SC[S]),sum.(data["D_household"][v_TCL] for s in SC[v_s]),50,40,10)
plot_v_histogram("EV connection [1/0]",sum.(data["EV_con"][:,:,s] for s in SC[S]),sum.(data["EV_con"][:,:,s] for s in SC[v_EV_sc]),50,20,10)
plot_v_histogram("Ambient temperature [C]",mean.(data["T_a"][:,s] for s in SC[S]),mean.(data["T_a"][:,s] for s in SC[v_s]),50,40,100)


plot_v_histogram("FCRN Activation [%]",mean.(data["FCRN_up_act"][:,s] for s in SC[S]),mean.(data["T_a"][:,s] for s in SC[v_s]),50,40,100)





include("plot_output_SC.jl")

plot_TCL(v_results)

plot_EV_dynamics(v_results)

#stacked_reserve_bid(v_results,data["C"])
#stacked_reserve_bid(v_results,data["F"])
#stacked_reserve_bid(v_results,data["E"])


stacked_activation(v_results)
#stacked_bar_q(v_results)
stacked_bar_q_mean(v_results)


FFR,FCRD,FCRN_up,FCRN_down,mFRR,all_AM = q_hour_mean(v_results)
Plots.plot(FCRD,FCRN_up,FCRN_down,all_AM, layout = 4)
Plots.plot(FFR,mFRR,all_AM, layout = 3)








################


data_name = "Household demand"
unit ="[kWh]"
data1 = data["D_household"][:,:,SC[S]]
v_c = v1_c
v_s = v1_s

test = unique(v1_c)

function plot_v_hour_avg(data_name,unit,data1,v_c,v_s,v_tu) 
    v_cu = unique(v_c)
    v_su = unique(v_s)
    if (v_c == false)
        p_data = sum(1/length(S)*data1[:,s] for s in S)
        v_data = sum(1/length(v_s)*data1[:,s] for s in v_s)
    elseif (v_c !== false)
        #C = [i for i in 1:size(data1)[1]]
        p_data = sum(1/length(S)*data1[c,:,s] for c in v_cu, s in S)
        v_data = sum(1/length(v_s)*data1[c,:,v_s[i]] for c in v_cu,i in eachindex(v_s))
    end
        Plots.plot(p_data,
        label="Avg. "*data_name,
        title = "Violating OS scenarios: "*data_name,
        xlabel = "Hours",
        ylabel = unit)
        Plots.plot!(v_data,
        label="Avg."*data_name*", $(length(v_s)) violating sc.")
        Plots.vline!([v_tu],linewidth = 1,color = "orangered", label ="$(length(v_tu)) hours with violation")
end



v1_c = [v1[i][1]  for i in eachindex(v1)]
v1_s = [v1[i][3] for i in eachindex(v1)] 
v1_tu = unique!([v1[i][2] for i in eachindex(v1)] )

v2_c = [v2[i][1]  for i in eachindex(v2)]
v2_s = [v2[i][3] for i in eachindex(v2)] 
v2_tu = unique!([v2[i][2] for i in eachindex(v2)])


v3_c = [v3[i][1]  for i in eachindex(v3)]
v3_s = [v3[i][3] for i in eachindex(v3)] 
v3_tu = unique!([v3[i][2] for i in eachindex(v3)])


v4_c = [v4[i][1]  for i in eachindex(v4)]
v4_s = [v4[i][3] for i in eachindex(v4)] 
v4_tu = unique!([v4[i][2] for i in eachindex(v4)])

#########################################################
########## TCL up Violations (v1) #######################


#plot_v_hour_avg("PV","kW",data["PV"][:,SC[S]],false,v1_s,v1_tu)
plot_v_hour_avg("Ambient temperature","[C]",data["T_a"][:,SC[S]],false,v1_s,v1_tu)
#plot_v_hour_avg("Household demand","[kWh]",data["D_household"][:,:,SC[S]],v1_c,v1_s,v1_tu)
plot_v_hour_avg("TCL Base Demand","[kWh]",data["D_base_TCL"][:,:,SC[S]],v1_c,v1_s,v1_tu)

#plot_v_hour_avg("FFR activation","[1/0]",data["FFR_act"][:,SC[S]],false,v1_s,v1_tu)
#plot_v_hour_avg("mFRR activation","[1/0]",data["mFRR_act"][:,SC[S]],false,v1_s,v1_tu)

FCRN_up = zeros(length(C),24,length(S))
FCRN_down = zeros(length(C),24,length(S))
FCRD_up = zeros(length(C),24,length(S))
for c in C
FCRN_up[c,:,:] = reshape(v_results["p_c_FCRN"][c,:].*data["FCRN_accept"][c,:,SC[S]].*data["FCRN_up_act"][:,SC[S]],(24,length(S)))
FCRN_down[c,:,:] = reshape(v_results["p_c_FCRN"][c,:].*data["FCRN_accept"][c,:,SC[S]].*data["FCRN_down_act"][:,SC[S]],(24,length(S)))
FCRD_up[c,:,:] = reshape(v_results["p_c_FCRD"][c,:].*data["FCRD_accept"][c,:,SC[S]].*data["FCRD_up_act"][:,SC[S]],(24,length(S)))
end
plot_v_hour_avg("FCRN up activation","[1/0]",FCRN_up,v1_c,v1_s,v1_tu) # TCL
plot_v_hour_avg("FCRN down activation","[1/0]",FCRN_down,v1_c,v1_s,v1_tu) # TCL
#plot_v_hour_avg("FCRD up activation","[1/0]",FCRD_up,v1_c,v1_s,v1_tu) # TCL




#########################################################
########## TCL up Violations (v2) #######################
# (no violations)
plot_v_hour_avg("PV","kW",data["PV"][:,SC[S]],false,v2_s,v2_tu)
plot_v_hour_avg("Ambient temperature","[C]",data["T_a"][:,SC[S]],false,v2_s,v2_tu)
plot_v_hour_avg("Household demand","[kWh]",data["D_household"][:,:,SC[S]],v2_c,v2_s,v2_tu)
plot_v_hour_avg("TCL Base Demand","[kWh]",data["D_base_TCL"][:,:,SC[S]],v2_c,v2_s,v2_tu)

plot_v_hour_avg("FFR activation","[1/0]",data["FFR_act"][:,SC[S]],false,v2_s,v2_tu)
plot_v_hour_avg("mFRR activation","[1/0]",data["mFRR_act"][:,SC[S]],false,v2_s,v2_tu)

FCRN_up = zeros(length(C),24,length(S))
FCRN_down = zeros(length(C),24,length(S))
FCRD_up = zeros(length(C),24,length(S))
for c in C
FCRN_up[c,:,:] = reshape(v_results["p_c_FCRN"][c,:].*data["FCRN_accept"][c,:,SC[S]].*data["FCRN_up_act"][:,SC[S]],(24,length(S)))
FCRN_down[c,:,:] = reshape(v_results["p_c_FCRN"][c,:].*data["FCRN_accept"][c,:,SC[S]].*data["FCRN_down_act"][:,SC[S]],(24,length(S)))
FCRD_up[c,:,:] = reshape(v_results["p_c_FCRD"][c,:].*data["FCRD_accept"][c,:,SC[S]].*data["FCRD_up_act"][:,SC[S]],(24,length(S)))
end
plot_v_hour_avg("FCRN up activation","[1/0]",FCRN_up,v2_c,v2_s,v2_tu)
plot_v_hour_avg("FCRN down activation","[1/0]",FCRN_down,v2_c,v2_s,v2_tu)
plot_v_hour_avg("FCRD up activation","[1/0]",FCRD_up,v2_c,v2_s,v2_tu)




#########################################################
########## EV up Violations (v3) #######################
# only hour 23, prosumer 6 and 7, where cars are home late.

#plot_v_hour_avg("PV","kW",data["PV"][:,SC[S]],false,v3_s,v3_tu)
#plot_v_hour_avg("Household demand","[kWh]",data["D_household"][:,:,SC[S]],v3_c,v3_s,v3_tu)
EV_con = zeros(length(data["C"]),24,length(S))
for e in data["E"]
EV_con[e,:,:] = reshape(data["EV_con"][e,:,SC[S]],(24,length(S)))
end
plot_v_hour_avg("EV connection","[1/0]",EV_con,v3_c,v3_s,v3_tu)


#plot_v_hour_avg("FFR activation","[1/0]",data["FFR_act"][:,SC[S]],false,v3_s,v3_tu)
#plot_v_hour_avg("mFRR activation","[1/0]",data["mFRR_act"][:,SC[S]],false,v3_s,v3_tu)

FCRN_up = zeros(length(C),24,length(S))
FCRN_down = zeros(length(C),24,length(S))
FCRD_up = zeros(length(C),24,length(S))
for c in C
FCRN_up[c,:,:] = reshape(v_results["p_c_FCRN"][c,:].*data["FCRN_accept"][c,:,SC[S]].*data["FCRN_up_act"][:,SC[S]],(24,length(S)))
FCRN_down[c,:,:] = reshape(v_results["p_c_FCRN"][c,:].*data["FCRN_accept"][c,:,SC[S]].*data["FCRN_down_act"][:,SC[S]],(24,length(S)))
FCRD_up[c,:,:] = reshape(v_results["p_c_FCRD"][c,:].*data["FCRD_accept"][c,:,SC[S]].*data["FCRD_up_act"][:,SC[S]],(24,length(S)))
end
plot_v_hour_avg("FCRN up activation","[1/0]",FCRN_up,v3_c,v3_s,v3_tu) 
plot_v_hour_avg("FCRN down activation","[1/0]",FCRN_down,v3_c,v3_s,v3_tu)
#plot_v_hour_avg("FCRD up activation","[1/0]",FCRD_up,v3_c,v3_s,v3_tu)


#########################################################
########## EV down Violations (v4) #######################
# only hour 23, prosumer 6 and 7, where cars are home late.

#plot_v_hour_avg("PV","kW",data["PV"][:,SC[S]],false,v4_s,v4_tu)
#plot_v_hour_avg("Household demand","[kWh]",data["D_household"][:,:,SC[S]],v4_c,v4_s,v4_tu)
EV_con = zeros(length(data["C"]),24,length(S))
for e in data["E"]
EV_con[e,:,:] = reshape(data["EV_con"][e,:,SC[S]],(24,length(S)))
end
plot_v_hour_avg("EV connection","[1/0]",EV_con,v4_c,v4_s,v4_tu)


#plot_v_hour_avg("FFR activation","[1/0]",data["FFR_act"][:,SC[S]],false,v4_s,v4_tu)
plot_v_hour_avg("mFRR activation","[1/0]",data["mFRR_act"][:,SC[S]],false,v4_s,v4_tu)

FCRN_up = zeros(length(C),24,length(S))
FCRN_down = zeros(length(C),24,length(S))
FCRD_up = zeros(length(C),24,length(S))
for c in C
FCRN_up[c,:,:] = reshape(v_results["p_c_FCRN"][c,:].*data["FCRN_accept"][c,:,SC[S]].*data["FCRN_up_act"][:,SC[S]],(24,length(S)))
FCRN_down[c,:,:] = reshape(v_results["p_c_FCRN"][c,:].*data["FCRN_accept"][c,:,SC[S]].*data["FCRN_down_act"][:,SC[S]],(24,length(S)))
FCRD_up[c,:,:] = reshape(v_results["p_c_FCRD"][c,:].*data["FCRD_accept"][c,:,SC[S]].*data["FCRD_up_act"][:,SC[S]],(24,length(S)))
end
plot_v_hour_avg("FCRN up activation","[1/0]",FCRN_up,v4_c,v4_s,v4_tu) 
plot_v_hour_avg("FCRN down activation","[1/0]",FCRN_down,v4_c,v4_s,v4_tu)
#plot_v_hour_avg("FCRD up activation","[1/0]",FCRD_up,v4_c,v4_s,v4_tu)



#=
include("plot_output_SC_avg.jl")

plot_TCL(v_results)

plot_EV_dynamics(v_results)

#stacked_reserve_bid(v_results,data["C"])
#stacked_reserve_bid(v_results,data["F"])
#stacked_reserve_bid(v_results,data["E"])


stacked_activation(v_results)
#stacked_bar_q(v_results)
stacked_bar_q_mean(v_results)


FFR,FCRD,FCRN_up,FCRN_down,mFRR,all_AM = q_hour_mean(v_results)
Plots.plot(FCRD,FCRN_up,FCRN_down,all_AM, layout = 4)
Plots.plot(FFR,mFRR,all_AM, layout = 3)


plot_prices(v_results)
plot_reserve_prices(v_results)

=#




