using JuMP, Gurobi
#using StochasticPrograms
using JLD2


#############################################
#### BASIC MODEL 2 - WITH MARKETS #############

function Cooperative_Risk_model1(Coop_model,data,SC,AM_risktype)
    S = [i for i in 1:length(SC)]
    ### Load Parameters ###
    #include("parameter_load.jl")
    
    ### Sets ###
    T = data["T"] # Number of time instances
    #N = data["N"] # The inflexible prosumers
    F = data["F"] # The fully flexible prosumers
    E = data["E"] # The prosumers with EV
    C = union(F,E) #union(N,F,E) # All prosumers
    
    #data["prob_s"] .*= 0 # adjusting sc probability to sample size
    #data["prob_s"][S] .+= 1/length(S)
    pi_s = zeros(length(SC)) .+= 1/length(SC)
    
    ### Parameters ###
    grid_in_price = data["grid_in_price"][:,SC]
    grid_out_price = data["grid_out_price"][:,SC] 

    balance_price = data["balance_price"][:,SC]
    
    FFR_price = data["FFR_price"][:,SC] 
    FFR_act = data["FFR_act"][:,SC] 
    FCRN_down_act = data["FCRN_down_act"][:,SC] 
    FCRN_up_act = data["FCRN_up_act"][:,SC]
    FCRN_accept = data["FCRN_accept"][:,:,SC]
    FCRN_price = data["FCRN_price"]###[:,:,SC] 
    FCRD_up_act = data["FCRD_up_act"][:,SC]
    FCRD_accept = data["FCRD_accept"][:,:,SC]
    FCRD_price = data["FCRD_price"]###[:,:,SC]
    mFRR_act = data["mFRR_act"][:,SC]
    mFRR_price = data["mFRR_price"][:,SC]

    PV = data["PV"][:,SC] # Data for PV
    PV_peak = data["PV_peak"] #(C)
    D_household = data["D_household"][:,:,SC] # Data for household demand for each prosumer (C, T, S)
    n_demand = data["n_demand"] #(C)
    EV_con = data["EV_con"][:,:,SC] # Binary indicating if EV is connected (T)
    
    η_c = data["η_c"] # Charging efficiency for EV battery
    η_d = data["η_d"] # Discharging efficiency for EV battery
    η_TCL = data["η_TCL"]   # COP for TCL
    D_base_TCL = data["D_base_TCL"][:,:,SC]

    Pmax_EC = data["Pmax_EC"] # Max power p_in/p_out of energy community
    Pmax_c = data["Pmax_c"] # Max power flow, p_grid of prosumer
    Pmax_TCL = data["Pmax_TCL"] # Max heatpump power d_TCL
    SOC_max = data["SOC_max"] # Max SOC / EV capacity
    Pmax_ev = data["Pmax_ev"] # EV charging/discharging capacity (E)
    
    #TCL Parameters
    R_th = data["R_th"]
    C_th = data["C_th"]
    
    T_a = data["T_a"][:,SC] # Ambient temperature [C]
    T_opt = data["T_opt"] # optimal indoor temperature [C]
    beta = data["beta"] # Weight of Temperature penalty in objective
    
    # Chance Constraint / CVaR Parameters
    epsilon = data["epsilon"] # confidence level (being 90% sure...)
    #gamma = data["gamma"] # (to deliver 95% of activation)
 

### SA Variables ###
    @variable(Coop_model,p_EC[T,S]) # Power flow to community (>0: import, <0: export)
    @variable(Coop_model,p_grid_out[T,S] >= 0) # Export power from whole community to the grid
    @variable(Coop_model,p_grid_in[T,S] >= 0) # Import power from whole community to the grid

    #@variable(Coop_model,p_mFRR[T] >= 0) 
    #@variable(Coop_model,p_FFR[T] >= 0) 
    #@variable(Coop_model,p_FCRD[T,S] >= 0) 
    #@variable(Coop_model,p_FCRN[T,S] >= 0) 

    @variable(Coop_model,p_mFRR_up[T,S] >= 0) 
    @variable(Coop_model,p_FFR_up[T,S] >= 0) 
    @variable(Coop_model,p_FCRD_up[T,S] >= 0) 
    @variable(Coop_model,p_FCRN_up[T,S] >= 0) 
    @variable(Coop_model,p_FCRN_down[T,S]>= 0) 

    @variable(Coop_model,q_mFRR[T,S] >= 0) 
    @variable(Coop_model,q_FFR[T,S] >= 0) 
    @variable(Coop_model,q_FCRD[T,S] >= 0) 
    @variable(Coop_model,q_FCRN_up[T,S] >= 0) 
    @variable(Coop_model,q_FCRN_down[T,S]>= 0)

# OS Slack variables
@variable(Coop_model,slack_TCL_up[F,T,S] >= 0)
@variable(Coop_model,slack_TCL_down[F,T,S] >= 0)
@variable(Coop_model,slack_EV_up[E,T,S] >= 0)
@variable(Coop_model,slack_EV_down[E,T,S] >= 0)


if (AM_risktype == "HJCC")
### Hourly JCC 
    @variable(Coop_model,y[T,S], Bin)
    println("model HJCC = $(AM_risktype)")
else
### JCC binary
    @variable(Coop_model,y[S], Bin)
    println("model otherthanHJCC = $(AM_risktype)")
end

if (AM_risktype == "CVaR")
### CVaR auxillary
    @variable(Coop_model,q_TCL_up[F,T,S]>=0)
    @variable(Coop_model,q_TCL_down[F,T,S]>=0)
    @variable(Coop_model,q_EV_up[E,T,S]>=0)
    @variable(Coop_model,q_EV_down[E,T,S]>=0)

    @variable(Coop_model,η_TCL_up)
    @variable(Coop_model,η_TCL_down)
    @variable(Coop_model,η_EV_up)
    @variable(Coop_model,η_EV_down)

    @variable(Coop_model,delta_TCL_up[S]>=0)
    @variable(Coop_model,delta_TCL_down[S]>=0)
    @variable(Coop_model,delta_EV_up[S]>=0)
    @variable(Coop_model,delta_EV_down[S]>=0)

    @variable(Coop_model,η_mFRR)
    @variable(Coop_model,η_FFR)
    @variable(Coop_model,η_FCRD)
    @variable(Coop_model,η_FCRN_up)
    @variable(Coop_model,η_FCRN_down)

    @variable(Coop_model,delta_mFRR[S]>=0)
    @variable(Coop_model,delta_FFR[S]>=0)
    @variable(Coop_model,delta_FCRD[S]>=0)
    @variable(Coop_model,delta_FCRN_up[S]>=0)
    @variable(Coop_model,delta_FCRN_down[S]>=0)
    println("model CVaR = $(AM_risktype)")
end

### All Prosumers Variables ###
    @variable(Coop_model, p[C,T,S]) 
    @variable(Coop_model,p_c_grid[C,T,S])    # Power for each prosumer, in each scenario

    # All Reserve
    @variable(Coop_model,p_c_mFRR[C,T] >= 0) 
    @variable(Coop_model,p_c_FCRD[C,T] >= 0) 
    @variable(Coop_model,p_c_FFR[C,T] >= 0) 
    @variable(Coop_model,p_c_FCRN[C,T] >= 0)

    # All Activation
    @variable(Coop_model,p_c_mFRR_up[C,T,S] >= 0) 
    @variable(Coop_model,p_c_FCRD_up[C,T,S] >= 0) 
    @variable(Coop_model,p_c_FFR_up[C,T,S] >= 0) 
    @variable(Coop_model,p_c_FCRN_up[C,T,S] >= 0) 
    @variable(Coop_model,p_c_FCRN_down[C,T,S] >= 0) 

    # All Ancillary slacks
    @variable(Coop_model,q_c_mFRR[C,T,S] >= 0) 
    @variable(Coop_model,q_c_FFR[C,T,S] >= 0) 
    @variable(Coop_model,q_c_FCRD[C,T,S] >= 0) 
    @variable(Coop_model,q_c_FCRN_up[C,T,S] >= 0)
    @variable(Coop_model,q_c_FCRN_down[C,T,S] >= 0)

    ### TCL Variables ###
    @variable(Coop_model,d_TCL[F,T,S] >= 0)   # Demand For Coop_modely shiftable demands at each scenario
    @variable(Coop_model,T_TCL[F,T,S])   # Indoor temperature for TCL owner

    ### EV Variables ###
    @variable(Coop_model,ev_ch[E,T,S] >= 0) # Charging power for EV prosumer
    @variable(Coop_model,ev_dch[E,T,S] >= 0) # Discharging power for EV prosumer
    @variable(Coop_model,SOC[E,T,S] >= 0)  # State of Charge variable for each EV battery


### SW Objective ###

    @objective(Coop_model, Min, sum(pi_s[s]*(grid_in_price[t,s]*p_grid_in[t,s] - grid_out_price[t,s]*p_grid_out[t,s]
                                - sum(mFRR_price[t,s]*(p_c_mFRR[c,t] - q_c_mFRR[c,t,s])
                                        + FFR_price[t,s]*(p_c_FFR[c,t]  - q_c_FFR[c,t,s])
                                        + FCRD_price[c]*FCRD_accept[c,t,s]*(p_c_FCRD[c,t]  - q_c_FCRD[c,t,s])
                                        + FCRN_price[c]*FCRN_accept[c,t,s]*(p_c_FCRN[c,t] - q_c_FCRN_up[c,t,s] - q_c_FCRN_down[c,t,s]) for c in C)
                                - balance_price[t,s]*p_FCRN_up[t,s] 
                                - balance_price[t,s]*p_FCRN_down[t,s] 
                                - balance_price[t,s]*p_mFRR_up[t,s]
                        ) for t in T, s in S)
                        + beta*sum(pi_s[s]*((T_a[t,s] < T_opt) ? (T_TCL[f,t,s] - T_opt)^2 : 0) for f in F, t in T, s in S))

### SA Constraints ###    
    # p_EC    Max total import/export
    @constraint(Coop_model, max_power_in[t in T, s in S], p_EC[t,s] <= Pmax_EC) #community import maximum
    @constraint(Coop_model, max_power_out[t in T, s in S], p_EC[t,s] >= - Pmax_EC) #community export maximum
    
    # Balance constraints
    @constraint(Coop_model, SA_Total_Flow[t in T, s in S], p_EC[t,s] == p_grid_in[t,s] - p_grid_out[t,s] - p_mFRR_up[t,s] - p_FFR_up[t,s] - p_FCRD_up[t,s] - p_FCRN_up[t,s] + p_FCRN_down[t,s])
    @constraint(Coop_model, SA_grid_flow[t in T, s in S], p_grid_in[t,s] - p_grid_out[t,s] == sum(p_c_grid[c,t,s] for c in C))
    @constraint(Coop_model, SA_mFRR_flow[t in T, s in S], p_mFRR_up[t,s] == sum(p_c_mFRR_up[c,t,s] for c in C))
    @constraint(Coop_model, SA_FFR_flow[t in T, s in S], p_FFR_up[t,s] == sum(p_c_FFR_up[c,t,s] for c in C))
    @constraint(Coop_model, SA_FCRD_flow[t in T, s in S], p_FCRD_up[t,s] == sum(p_c_FCRD_up[c,t,s] for c in C))
    @constraint(Coop_model, SA_FCRN_up_flow[t in T, s in S], p_FCRN_up[t,s] == sum(p_c_FCRN_up[c,t,s] for c in C))
    @constraint(Coop_model, SA_FCRN_down_flow[t in T, s in S], p_FCRN_down[t,s] == sum(p_c_FCRN_down[c,t,s] for c in C))

    # Ancillary Bid balance
    #@constraint(Coop_model, SA_mFRR_bid[t in T], p_mFRR[t] == sum(p_c_mFRR[c,t] for c in C))
    #@constraint(Coop_model, SA_FFR_bid[t in T], p_FFR[t] == sum(p_c_FFR[c,t] for c in C))
    #@constraint(Coop_model, SA_FCRD_bid[t in T], p_FCRD[t,s] == sum(p_c_FCRD[c,t]*FCRD_accept[c,t,s] for c in C))
    #@constraint(Coop_model, SA_FCRN_up_bid[t in T], p_FCRN[t,s] == sum(p_c_FCRN[c,t]*FCRN_accept[c,t,s] for c in C))

    # SA Reserve slacks
    @constraint(Coop_model, q_reserve_mFRR[t in T, s in S], q_mFRR[t,s] == sum(q_c_mFRR[c,t,s] for c in C)) 
    @constraint(Coop_model, q_reserve_FFR[t in T, s in S], q_FFR[t,s] == sum(q_c_FFR[c,t,s] for c in C)) 
    @constraint(Coop_model, q_reserve_FCRD[t in T, s in S], q_FCRD[t,s] == sum(q_c_FCRD[c,t,s] for c in C)) 
    @constraint(Coop_model, q_reserve_FCRN_up[t in T, s in S], q_FCRN_up[t,s] == sum(q_c_FCRN_up[c,t,s] for c in C)) 
    @constraint(Coop_model, q_reserve_FCRN_down[t in T, s in S], q_FCRN_down[t,s] == sum(q_c_FCRN_down[c,t,s] for c in C))

if (AM_risktype == "no_risk")
    println("Risk_type: no_risk")
    for c in C
        for t in T
            for s in S 
                fix(q_c_FCRN_up[c,t,s], 0; force = true)
                fix(q_c_FCRN_down[c,t,s], 0; force = true)
                fix(q_c_FCRD[c,t,s], 0; force = true)
                fix(q_c_mFRR[c,t,s], 0; force = true)
                fix( q_c_FFR[c,t,s], 0; force = true)
            end
        end
    end

    @constraint(Coop_model, TCL_up_reserve_limit[f in F,t in T,s in S], FCRN_accept[f,t,s]*p_c_FCRN[f,t] + p_c_mFRR[f,t] + FCRD_accept[f,t,s]*p_c_FCRD[f,t] + p_c_FFR[f,t] <= D_base_TCL[f,t,s]) # Limit on ancillary upreg reserve bid
    @constraint(Coop_model, TCL_down_reserve_limit[f in F,t in T,s in S], FCRN_accept[f,t,s]*p_c_FCRN[f,t] <= Pmax_TCL - D_base_TCL[f,t,s]) # Limit on ancillary downreg reserve bid
    @constraint(Coop_model, EV_up_reserve_limit[e in E,t in T,s in S], FCRN_accept[e,t,s]*p_c_FCRN[e,t] + p_c_mFRR[e,t] + FCRD_accept[e,t,s]*p_c_FCRD[e,t] + p_c_FFR[e,t] <= EV_con[e,t,s]*Pmax_ev[e]) # Limit on ancillary upward reserve bid
    @constraint(Coop_model, EV_down_reserve_limit[e in E,t in T,s in S], FCRN_accept[e,t,s]*p_c_FCRN[e,t] <= EV_con[e,t,s]*Pmax_ev[e]) # Limit on ancillary upward reserve bid


elseif (AM_risktype == "CVaR")
    println("Risk_type: CVaR")
### CVaR chance constraints for Ancillary markets
    # Reserve > 95%
    @constraint(Coop_model, TCL_up_reserve_limit[f in F,t in T,s in S], FCRN_accept[f,t,s]*p_c_FCRN[f,t] + p_c_mFRR[f,t] + FCRD_accept[f,t,s]*p_c_FCRD[f,t] + p_c_FFR[f,t] <= D_base_TCL[f,t,s] + q_TCL_up[f,t,s]) # Limit on ancillary upreg reserve bid
    @constraint(Coop_model, TCL_down_reserve_limit[f in F,t in T,s in S], FCRN_accept[f,t,s]*p_c_FCRN[f,t] <= Pmax_TCL - D_base_TCL[f,t,s] + q_TCL_down[f,t,s]) # Limit on ancillary downreg reserve bid
    @constraint(Coop_model, EV_up_reserve_limit[e in E,t in T,s in S], FCRN_accept[e,t,s]*p_c_FCRN[e,t] + p_c_mFRR[e,t] + FCRD_accept[e,t,s]*p_c_FCRD[e,t] + p_c_FFR[e,t] <= EV_con[e,t,s]*Pmax_ev[e] + q_EV_up[e,t,s]) # Limit on ancillary upward reserve bid
    @constraint(Coop_model, EV_down_reserve_limit[e in E,t in T,s in S], FCRN_accept[e,t,s]*p_c_FCRN[e,t] <= EV_con[e,t,s]*Pmax_ev[e] + q_EV_down[e,t,s]) # Limit on ancillary upward reserve bid

    @constraint(Coop_model, CVaR_TCL_up[s in S], delta_TCL_up[s] >= sum(q_TCL_up[f,t,s] for f in F,t in T) - η_TCL_up) # CVaR constraint for each scenario
    @constraint(Coop_model, CVaR_TCL_down[s in S], delta_TCL_down[s] >= sum(q_TCL_down[f,t,s] for f in F,t in T) - η_TCL_down)
    @constraint(Coop_model, CVaR_EV_up[s in S], delta_EV_up[s] >= sum(q_EV_up[e,t,s] for e in E,t in T) - η_EV_up)
    @constraint(Coop_model, CVaR_EV_down[s in S], delta_EV_down[s] >= sum(q_EV_down[e,t,s] for e in E,t in T) - η_EV_down)
    
    @constraint(Coop_model, CVaR_TCL_up2[s in S], η_TCL_up + 1/(1 - epsilon)*sum(pi_s[s]*delta_TCL_up[s] for s in S) <= epsilon*sum(pi_s[s]*(p_c_mFRR[f,t] + FCRN_accept[f,t,s]*p_c_FCRN[f,t] + FCRD_accept[f,t,s]*p_c_FCRD[f,t] + p_c_FFR[f,t]) for f in F,t in T,s in S))
    @constraint(Coop_model, CVaR_TCL_down2[s in S], η_TCL_down + 1/(1 - epsilon)*sum(pi_s[s]*delta_TCL_down[s] for s in S) <= epsilon*sum(pi_s[s]*(FCRN_accept[f,t,s]*p_c_FCRN[f,t]) for f in F,t in T,s in S))
    @constraint(Coop_model, CVaR_EV_up2, η_EV_up + 1/(1 - epsilon)*sum(pi_s[s]*delta_EV_up[s] for s in S) <= epsilon*sum(pi_s[s]*(p_c_mFRR[e,t] + FCRN_accept[e,t,s]*p_c_FCRN[e,t] + FCRD_accept[e,t,s]*p_c_FCRD[e,t] + p_c_FFR[e,t]) for e in E, t in T,s in S))
    @constraint(Coop_model, CVaR_EV_down2, η_FCRN_up + 1/(1 - epsilon)*sum(pi_s[s]*delta_FCRN_up[s] for s in S) <= epsilon*sum(pi_s[s]*(FCRN_accept[e,t,s]*p_c_FCRN[e,t]) for e in E, t in T,s in S))
    
    #@constraint(Coop_model, CVaR_TCL_up2[s in S], η_TCL_up + 1/(1 - epsilon)*sum(pi_s[s]*delta_TCL_up[s] for s in S) <= epsilon*sum(pi_s[s]*D_base_TCL[f,t,s] for f in F,t in T,s in S))
    #@constraint(Coop_model, CVaR_TCL_down2[s in S], η_TCL_down + 1/(1 - epsilon)*sum(pi_s[s]*delta_TCL_down[s] for s in S) <= epsilon*sum(pi_s[s]*(Pmax_TCL - D_base_TCL[f,t,s]) for f in F,t in T,s in S))
    #@constraint(Coop_model, CVaR_EV_up2, η_EV_up + 1/(1 - epsilon)*sum(pi_s[s]*delta_EV_up[s] for s in S) <= epsilon*sum(pi_s[s]*(EV_con[e,t,s]*Pmax_ev[e]) for e in E, t in T,s in S))
    #@constraint(Coop_model, CVaR_EV_down2, η_FCRN_up + 1/(1 - epsilon)*sum(pi_s[s]*delta_FCRN_up[s] for s in S) <= epsilon*sum(pi_s[s]*(EV_con[e,t,s]*Pmax_ev[e]) for e in E, t in T,s in S))

    # Delivery > 95%
    @constraint(Coop_model, CVaR_mFRR[s in S], delta_mFRR[s] >= sum(q_mFRR[t,s] for t in T) - η_mFRR) # CVaR constraint for each scenario
    @constraint(Coop_model, CVaR_FFR[s in S], delta_FFR[s] >= sum(q_FFR[t,s] for t in T) - η_FFR)
    @constraint(Coop_model, CVaR_FCRD[s in S], delta_FCRD[s] >= sum(q_FCRD[t,s] for t in T) - η_FCRD)
    @constraint(Coop_model, CVaR_FCRN_up[s in S], delta_FCRN_up[s] >= sum(q_FCRN_up[t,s] for t in T) - η_FCRN_up)
    @constraint(Coop_model, CVaR_FCRN_down[s in S], delta_FCRN_down[s] >= sum(q_FCRN_down[t,s] for t in T) - η_FCRN_down)

    @constraint(Coop_model, CVaR_mFRR2[s in S], η_mFRR + 1/(1 - epsilon)*sum(pi_s[s]*delta_mFRR[s] for s in S) <= epsilon*sum(pi_s[s]*mFRR_act[t,s]*p_c_mFRR[c,t] for c in C,t in T,s in S))
    @constraint(Coop_model, CVaR_FFR2[s in S], η_FFR + 1/(1 - epsilon)*sum(pi_s[s]*delta_FFR[s] for s in S) <= epsilon*sum(pi_s[s]*FFR_act[t,s]*p_c_FFR[c,t] for c in C,t in T,s in S))
    @constraint(Coop_model, CVaR_FCRD2, η_FCRD + 1/(1 - epsilon)*sum(pi_s[s]*delta_FCRD[s] for s in S) <= epsilon*sum(pi_s[s]*FCRD_accept[c,t,s]*FCRD_up_act[t,s]*p_c_FCRD[c,t] for c in C,t in T,s in S))
    @constraint(Coop_model, CVaR_FCRN_op2, η_FCRN_up + 1/(1 - epsilon)*sum(pi_s[s]*delta_FCRN_up[s] for s in S) <= epsilon*sum(pi_s[s]*FCRN_accept[c,t,s]*FCRN_up_act[t,s]*p_c_FCRN[c,t] for c in C,t in T,s in S))
    @constraint(Coop_model, CVaR_FCRN_down2, η_FCRN_down + 1/(1 - epsilon)*sum(pi_s[s]*delta_FCRN_down[s] for s in S) <= epsilon*sum(pi_s[s]*FCRN_accept[c,t,s]*FCRN_down_act[t,s]*p_c_FCRN[c,t] for c in C,t in T,s in S))
    println("model CVaR = $(AM_risktype)")
elseif (AM_risktype == "HJCC" )
    println("Risk_type: HJCC")
    ### Joint Chance Constraints for individual Hours for Ancillary markets
        @constraint(Coop_model, HJCC_mFRR[t in T, s in S], q_mFRR[t,s] <= y[t,s]*Pmax_EC)
        @constraint(Coop_model, HJCC_FFR[t in T, s in S], q_FFR[t,s] <= y[t,s]*Pmax_EC)
        @constraint(Coop_model, HJCC_FCRD[t in T,s in S], q_FCRD[t,s] <= y[t,s]*Pmax_EC)
        @constraint(Coop_model, HJCC_FCRN_up[t in T, s in S], q_FCRN_up[t,s] <= y[t,s]*Pmax_EC)
        @constraint(Coop_model, HJCC_FCRN_down[t in T, s in S], q_FCRN_down[t,s] <= y[t,s]*Pmax_EC)
    
        @constraint(Coop_model, HJCC[t in T], sum(y[t,s] for s in S) <= length(S)*(epsilon))  
    println("model HJCC = $(AM_risktype)")
elseif (AM_risktype == "JCC" )
    println("Risk_type: JCC")
### Joint Chance Constraints for Ancillary markets
    @constraint(Coop_model, TCL_up_reserve_limit[f in F,t in T,s in S], FCRN_accept[f,t,s]*p_c_FCRN[f,t] + p_c_mFRR[f,t] + FCRD_accept[f,t,s]*p_c_FCRD[f,t] + p_c_FFR[f,t] <= (1 - y[s])*D_base_TCL[f,t,s] + y[s]*Pmax_TCL) # Limit on ancillary upreg reserve bid
    @constraint(Coop_model, TCL_down_reserve_limit[f in F,t in T,s in S], FCRN_accept[f,t,s]*p_c_FCRN[f,t] <= (1 - y[s])*(Pmax_TCL - D_base_TCL[f,t,s]) + y[s]*Pmax_TCL) # Limit on ancillary downreg reserve bid
    @constraint(Coop_model, EV_up_reserve_limit[e in E,t in T,s in S], FCRN_accept[e,t,s]*p_c_FCRN[e,t] + p_c_mFRR[e,t] + FCRD_accept[e,t,s]*p_c_FCRD[e,t] + p_c_FFR[e,t] <= (1 - y[s])*EV_con[e,t,s]*Pmax_ev[e] + y[s]*Pmax_ev[e]) # Limit on ancillary upward reserve bid
    @constraint(Coop_model, EV_down_reserve_limit[e in E,t in T,s in S], FCRN_accept[e,t,s]*p_c_FCRN[e,t] <= (1 - y[s])*EV_con[e,t,s]*Pmax_ev[e]  + y[s]*Pmax_ev[e]) # Limit on ancillary upward reserve bid

    @constraint(Coop_model, JCC_mFRR[s in S], sum(q_mFRR[t,s] for t in T) <= y[s]*24*Pmax_EC)
    @constraint(Coop_model, JCC_FFR[s in S], sum(q_FFR[t,s] for t in T) <= y[s]*24*Pmax_EC)
    @constraint(Coop_model, JCC_FCRD[s in S], sum(q_FCRD[t,s] for t in T) <= y[s]*24*Pmax_EC)
    @constraint(Coop_model, JCC_FCRN_up[s in S], sum(q_FCRN_up[t,s] for t in T) <= y[s]*24*Pmax_EC)
    @constraint(Coop_model, JCC_FCRN_down[s in S], sum(q_FCRN_down[t,s] for t in T) <= y[s]*24*Pmax_EC)

    @constraint(Coop_model, JCC, sum(y[s] for s in S) <= length(S)*(epsilon)) 
    println("model JCC = $(AM_risktype)")
elseif (AM_risktype == "CC" )
    println("Risk_type: CC")
    println("model CC = $(AM_risktype)")
### CC ancillary Variables
    @variable(Coop_model,y_mFRR[S], Bin)
    @variable(Coop_model,y_FFR[S], Bin)
    @variable(Coop_model,y_FCRD[S], Bin)
    @variable(Coop_model,y_FCRN_up[S], Bin)
    @variable(Coop_model,y_FCRN_down[S], Bin)

### Joint Chance Constraints for Ancillary markets
    @constraint(Coop_model, CC_mFRR[s in S], sum(q_mFRR[t,s] for t in T) <= y_mFRR[s]*24*Pmax_EC)
    @constraint(Coop_model, CC_FFR[s in S], sum(q_FFR[t,s] for t in T) <= y_FFR[s]*24*Pmax_EC)
    @constraint(Coop_model, CC_FCRD[s in S], sum(q_FCRD[t,s] for t in T) <= y_FCRD[s]*24*Pmax_EC)
    @constraint(Coop_model, CC_FCRN_up[s in S], sum(q_FCRN_up[t,s] for t in T) <= y_FCRN_up[s]*24*Pmax_EC)
    @constraint(Coop_model, CC_FCRN_down[s in S], sum(q_FCRN_down[t,s] for t in T) <= y_FCRN_down[s]*24*Pmax_EC)

    @constraint(Coop_model, CC1, sum(y_mFRR[s] for s in S) <= length(S)*(epsilon))
    @constraint(Coop_model, CC2, sum(y_FFR[s] for s in S) <= length(S)*(epsilon)) 
    @constraint(Coop_model, CC3, sum(y_FCRD[s] for s in S) <= length(S)*(epsilon)) 
    @constraint(Coop_model, CC4, sum(y_FCRN_up[s] for s in S) <= length(S)*(epsilon)) 
    @constraint(Coop_model, CC5, sum(y_FCRN_down[s] for s in S) <= length(S)*(epsilon)) 

end

### All prosumers constraints ####
# Balance constraints
@constraint(Coop_model, Prosumer_max_power_out[c in C,t in T, s in S], p[c,t,s] >= -Pmax_c) #Maximum prosumer eksport
@constraint(Coop_model, Prosumer_max_power_in[c in C,t in T, s in S], p[c,t,s] <= Pmax_c) #Maximum prosumer import

@constraint(Coop_model, Prosumer_max_bid[c in C,t in T, s in S], FCRN_accept[c,t,s]*p_c_FCRN[c,t] + p_c_mFRR[c,t] + FCRD_accept[c,t,s]*p_c_FCRD[c,t] + p_c_FFR[c,t] <= Pmax_c) #Maximum prosumer bid

# Ancillary Constraints
@constraint(Coop_model, q_max_dump_FCRN_up[c in C,t in T, s in S], p_c_FCRN_up[c,t,s] + q_c_FCRN_up[c,t,s] == FCRN_accept[c,t,s]*FCRN_up_act[t,s]*p_c_FCRN[c,t])
@constraint(Coop_model, q_max_dump_FCRN_down[c in C,t in T, s in S], p_c_FCRN_down[c,t,s] + q_c_FCRN_down[c,t,s] == FCRN_accept[c,t,s]*FCRN_down_act[t,s]*p_c_FCRN[c,t])
@constraint(Coop_model, q_max_dump_FCRD_up[c in C,t in T, s in S], p_c_FCRD_up[c,t,s] + q_c_FCRD[c,t,s] == FCRD_accept[c,t,s]*FCRD_up_act[t,s]*p_c_FCRD[c,t])
@constraint(Coop_model, q_max_dump_mFRR[c in C,t in T, s in S], p_c_mFRR_up[c,t,s] + q_c_mFRR[c,t,s] == p_c_mFRR[c,t].*mFRR_act[t,s])
@constraint(Coop_model, q_max_dump_FFR[c in C,t in T, s in S], p_c_FFR_up[c,t,s] + q_c_FFR[c,t,s] == p_c_FFR[c,t].*FFR_act[t,s])

### TCL ###
@constraint(Coop_model, TCL_balance[f in F, t in T, s in S], p[f,t,s] + PV[t,s]*PV_peak[f]*n_demand[f] - D_household[f,t,s] - d_TCL[f,t,s] == 0) #flexible prosumer energy balance
@constraint(Coop_model, TCL_ancillary[f in F, t in T, s in S], p[f,t,s] == p_c_grid[f,t,s] - p_c_FCRN_up[f,t,s] + p_c_FCRN_down[f,t,s] - p_c_mFRR_up[f,t,s] - p_c_FCRD_up[f,t,s] - p_c_FFR_up[f,t,s])
@constraint(Coop_model, TCL_demand[f in F,t in T,s in S], T_TCL[f,t,s] == (t>1 ? T_TCL[f,t-1,s] : T_TCL[f,24,s]) - 1/(R_th[f]*C_th[f])*((t>1 ? T_TCL[f,t-1,s] : T_TCL[f,24,s])-T_a[t,s]) + η_TCL/C_th[f]*d_TCL[f,t,s]) #TCL adjusted demand 
@constraint(Coop_model, d_TCL_capacity[f in F,t in T,s in S], d_TCL[f,t,s] <= Pmax_TCL)
 

### EV ###
# Balance
@constraint(Coop_model, EV_balance[e in E,t in T, s in S], p[e,t,s] + PV[t,s]*PV_peak[e]*n_demand[e] - D_household[e,t,s] - ev_ch[e,t,s] + ev_dch[e,t,s] == 0) # Balance constraint for EV prosumer
@constraint(Coop_model, EV_ancillary[e in E, t in T, s in S], p[e,t,s] == p_c_grid[e,t,s] - p_c_FCRN_up[e,t,s] + p_c_FCRN_down[e,t,s] - p_c_mFRR_up[e,t,s] - p_c_FCRD_up[e,t,s] - p_c_FFR_up[e,t,s]) # OBS test


# SOC
@constraint(Coop_model, EV_SOC[e in E, t in T, s in S], SOC[e,t,s] == (t>1 ? SOC[e,t-1,s] : SOC[e,24,s]) + η_c*ev_ch[e,t,s] - η_d*ev_dch[e,t,s]
                                                           - 0.6*EV_con[e,t,s]*(EV_con[e,t,s] - (t>1 ? EV_con[e,t-1,s] : EV_con[e,24,s]))*SOC_max[e]) # State of Charge constraint, assuming next day is similar
@constraint(Coop_model, SOC_drive_demand[e in E,t in T, s in S], SOC[e,t,s] >= 0.8*SOC_max[e]*(EV_con[e,t,s] - (t<24 ? EV_con[e,t+1,s] : EV_con[e,1,s]))) # min SOC, when EV disconnects

# Limits
@constraint(Coop_model, EV_SOC_max[e in E,t in T, s in S], SOC[e,t,s] <= SOC_max[e]) # Maximum SOC Capacity 
@constraint(Coop_model, EV_charge_max[e in E,t in T, s in S],  ev_ch[e,t,s] <= Pmax_ev[e]*EV_con[e,t,s]) # Maximum Charging Power
@constraint(Coop_model, EV_discharge_max[e in E,t in T, s in S], ev_dch[e,t,s] <= Pmax_ev[e]*EV_con[e,t,s]) # Maximum Discharging power

# Additions

return Coop_model,SC,pi_s
end

function Cooperative_Risk_model2(Coop_model,data,SC,AM_risktype)
    trade_model,SC,pi_s = Cooperative_Risk_model1(Coop_model,data,SC,AM_risktype)

    T = data["T"]
    C = data["C"]
    F = data["F"]
    E = data["E"]
    S = [i for i in eachindex(SC)] # keep!
    


    # Allowing q_c to take negative values, so prosumers can deliver more than their reserve bid and trade reserves.
    delete.(trade_model, trade_model[:q_c_mFRR][C,T,S])
    delete.(trade_model, trade_model[:q_c_FFR][C,T,S])
    delete.(trade_model, trade_model[:q_c_FCRD][C,T,S])
    delete.(trade_model, trade_model[:q_c_FCRN_up][C,T,S])
    delete.(trade_model, trade_model[:q_c_FCRN_down][C,T,S])


    unregister.(trade_model, :q_c_mFRR) #?
    unregister.(trade_model, trade_model[:q_c_FFR][C,T,S])
    unregister.(trade_model, trade_model[:q_c_FCRD][C,T,S])
    unregister.(trade_model, trade_model[:q_c_FCRN_up][C,T,S])
    unregister.(trade_model, trade_model[:q_c_FCRN_down][C,T,S])


    @variable(trade_model,q_c_mFRR[C,T,S]) 
    @variable(trade_model,q_c_FFR[C,T,S]) 
    @variable(trade_model,q_c_FCRD[C,T,S]) 
    @variable(trade_model,q_c_FCRN_up[C,T,S])
    @variable(trade_model,q_c_FCRN_down[C,T,S])

end

function solvem(Coop_model,SC,AM_risktype)
    S = [i for i in 1:length(SC)]
### Solving problem ###
    optimize!(Coop_model)
    println("Termination status: $(termination_status(Coop_model))")

# Solution
if termination_status(Coop_model) == MOI.OPTIMAL
    println("SA Optimal objective value: $(objective_value(Coop_model))")
    println("Using Risk type: $(AM_risktype)")
###
    T = data["T"]
    F = data["F"]
    E = data["E"]
    C = union(F,E)
    T_opt = data["T_opt"]
    T_a = data["T_a"][:,SC]
    beta = data["beta"]
    pi_s = ones(length(S))./length(S)

    
    T_TCL = [value.(Coop_model[:T_TCL][f,t,s]) for f in F,t in T, s in S]

    p_mFRR = [value.(sum(Coop_model[:p_c_mFRR][C,t])) for t in T]
    p_FFR = [value.(sum(Coop_model[:p_c_FFR][C,t])) for t in T]
    p_FCRD = [value.(sum(Coop_model[:p_c_FCRD][C,t])) for t in T]
    p_FCRN = [value.(sum(Coop_model[:p_c_FCRN][C,t])) for t in T]

    p_grid_out = [value.(Coop_model[:p_grid_out][t,s]) for t in T, s in S]
    p_grid_in = [value.(Coop_model[:p_grid_in][t,s]) for t in T, s in S]
    p_FCRN_up = [value.(Coop_model[:p_FCRN_up][t,s]) for t in T, s in S]
    p_FCRN_down = [value.(Coop_model[:p_FCRN_down][t,s]) for t in T, s in S]
    p_mFRR_up = [value.(Coop_model[:p_mFRR_up][t,s]) for t in T, s in S]
    p_FCRD_up = [value.(Coop_model[:p_FCRD_up][t,s]) for t in T, s in S]
    p_FFR_up = [value.(Coop_model[:p_FFR_up][t,s]) for t in T, s in S]

    p_c_FCRN = [value.(Coop_model[:p_c_FCRN][c,t]) for c in C, t in T]
    p_c_FCRN = [value.(Coop_model[:p_c_FCRN][c,t]) for c in C, t in T]
    p_c_mFRR = [value.(Coop_model[:p_c_mFRR][c,t]) for c in C, t in T]
    p_c_FCRD = [value.(Coop_model[:p_c_FCRD][c,t]) for c in C, t in T]
    p_c_FFR = [value.(Coop_model[:p_c_FFR][c,t]) for c in C, t in T]

    q_c_FCRN_up = [value.(Coop_model[:q_c_FCRN_up][c,t,s]) for c in C, t in T, s in S]
    q_c_FCRN_down = [value.(Coop_model[:q_c_FCRN_down][c,t,s]) for c in C, t in T, s in S]
    q_c_mFRR = [value.(Coop_model[:q_c_mFRR][c,t,s]) for c in C, t in T, s in S]
    q_c_FFR = [value.(Coop_model[:q_c_FFR][c,t,s]) for c in C, t in T, s in S]
    q_c_FCRD = [value.(Coop_model[:q_c_FCRD][c,t,s]) for c in C, t in T, s in S]

    grid_in_price = data["grid_in_price"][:,SC] 
    grid_out_price = data["grid_out_price"][:,SC] 
    balance_price = data["balance_price"][:,SC]
    FFR_price = data["FFR_price"][:,SC] 
    #FFR_act = data["FFR_act"][:,SC] 
    #FCRN_down_act = data["FCRN_down_act"]#[:,SC] 
    #FCRN_up_act = data["FCRN_up_act"]#[:,SC]
    FCRN_accept = data["FCRN_accept"][:,:,SC]
    FCRN_price = data["FCRN_price"]###[:,:,SC] 
    #FCRD_up_act = data["FCRD_up_act"][:,SC]
    FCRD_accept = data["FCRD_accept"][:,:,SC]
    FCRD_price = data["FCRD_price"]###[:,:,SC]
    #mFRR_act = data["mFRR_act"][:,SC]
    mFRR_price = data["mFRR_price"][:,SC]

    
    obj_s_base = zeros(length(S))
    ### Single scenario objective
    for s in S
    obj_s_base[s] = sum((grid_in_price[t,s]*p_grid_in[t,s] - grid_out_price[t,s]*p_grid_out[t,s]
    - sum(mFRR_price[t,s]*(p_c_mFRR[c,t] - q_c_mFRR[c,t,s])
    + FFR_price[t,s]*(p_c_FFR[c,t]  - q_c_FFR[c,t,s])
    + FCRD_price[c]*FCRD_accept[c,t,s]*(p_c_FCRD[c,t]  - q_c_FCRD[c,t,s])
    + FCRN_price[c]*FCRN_accept[c,t,s]*(p_c_FCRN[c,t] - q_c_FCRN_up[c,t,s] - q_c_FCRN_down[c,t,s]) for c in C)
    - balance_price[t,s]*p_FCRN_up[t,s] 
    - balance_price[t,s]*p_FCRN_down[t,s] 
    - balance_price[t,s]*p_mFRR_up[t,s]
    ) for t in T)
    + beta*sum(((T_a[t,s] < T_opt) ? (T_TCL[f,t,s] - T_opt)^2 : 0) for f in F, t in T)
    end

obj_tjeck = sum(pi_s[s]*(grid_in_price[t,s]*p_grid_in[t,s] - grid_out_price[t,s]*p_grid_out[t,s]
                - sum(mFRR_price[t,s]*(p_c_mFRR[c,t] - q_c_mFRR[c,t,s])
                + FFR_price[t,s]*(p_c_FFR[c,t]  - q_c_FFR[c,t,s])
                + FCRD_price[c]*FCRD_accept[c,t,s]*(p_c_FCRD[c,t]  - q_c_FCRD[c,t,s])
                + FCRN_price[c]*FCRN_accept[c,t,s]*(p_c_FCRN[c,t] - q_c_FCRN_up[c,t,s] - q_c_FCRN_down[c,t,s]) for c in C)
                - balance_price[t,s]*p_FCRN_up[t,s] 
                - balance_price[t,s]*p_FCRN_down[t,s] 
                - balance_price[t,s]*p_mFRR_up[t,s]
                ) for t in T, s in S)
                + beta*sum(pi_s[s]*((T_a[t,s] < T_opt) ? (T_TCL[f,t,s] - T_opt)^2 : 0) for f in F, t in T, s in S)


    sol_Coop_model = Dict("obj" => objective_value(Coop_model),
            #
            "obj_s_base" => obj_s_base,
            "obj_tjeck" => obj_tjeck,
            "solve_time" => solve_time(Coop_model),
            "pi_s" => pi_s,
            "C" => C,
            "SC" => SC,
            "S" => [i for i in eachindex(SC)],
            #
            "p_EC" => [value.(Coop_model[:p_EC][t,s]) for t in T, s in S],
            "p_FCRN" => p_FCRN, #[value.(Coop_model[:p_FCRN][T]),
            "p_mFRR" => p_mFRR,#[value.(Coop_model[:p_mFRR][T]),
            "p_FCRD" => p_FCRD,#[value.(Coop_model[:p_FCRD][T]),
            "p_FFR" => p_FFR, #[value.(Coop_model[:p_FFR][T]),
            #
            # SA Activation
            "p_grid_out" => [value.(Coop_model[:p_grid_out][t,s]) for t in T, s in S],
            "p_grid_in" => [value.(Coop_model[:p_grid_in][t,s]) for t in T, s in S],
            "p_FCRN_up" => [value.(Coop_model[:p_FCRN_up][t,s]) for t in T, s in S],
            "p_FCRN_down" => [value.(Coop_model[:p_FCRN_down][t,s]) for t in T, s in S],
            "p_mFRR_up" => [value.(Coop_model[:p_mFRR_up][t,s]) for t in T, s in S],
            "p_FCRD_up" => [value.(Coop_model[:p_FCRD_up][t,s]) for t in T, s in S],
            "p_FFR_up" => [value.(Coop_model[:p_FFR_up][t,s]) for t in T, s in S],
            #
            # All prosumers
            "p_c_FCRN" => [value.(Coop_model[:p_c_FCRN][c,t]) for c in C, t in T],
            "p_c_mFRR" => [value.(Coop_model[:p_c_mFRR][c,t]) for c in C, t in T],
            "p_c_FCRD" => [value.(Coop_model[:p_c_FCRD][c,t]) for c in C, t in T],
            "p_c_FFR" => [value.(Coop_model[:p_c_FFR][c,t]) for c in C, t in T],
            "p_c" => [value.(Coop_model[:p][c,t,s]) for c in C, t in T, s in S],
            "p_c_grid" => [value.(Coop_model[:p_c_grid][c,t,s]) for c in C, t in T, s in S],
            "p_c_FCRN_up" => [value.(Coop_model[:p_c_FCRN_up][c,t,s]) for c in C, t in T, s in S],
            "p_c_FCRN_down" => [value.(Coop_model[:p_c_FCRN_down][c,t,s]) for c in C, t in T, s in S],
            "p_c_mFRR_up" => [value.(Coop_model[:p_c_mFRR_up][c,t,s]) for c in C, t in T, s in S],
            "p_c_FCRD_up" => [value.(Coop_model[:p_c_FCRD_up][c,t,s]) for c in C, t in T, s in S],
            "p_c_FFR_up" => [value.(Coop_model[:p_c_FFR_up][c,t,s]) for c in C, t in T, s in S],
            #
            #TCL
            "d_TCL" => [value.(Coop_model[:d_TCL][f,t,s]) for f in F, t in T, s in S],
            "T_TCL" => [value.(Coop_model[:T_TCL][f,t,s]) for f in F, t in T, s in S],
            "c_TCL" => [(value.(Coop_model[:T_TCL][f,t,s]) .- T_opt).^2 for f in F, t in T, s in S],
            #
            # EV
            "SOC" => [value.(Coop_model[:SOC][e,t,s]) for e in E, t in T, s in S],
            "ev_ch" => [value.(Coop_model[:ev_ch][e,t,s]) for e in E, t in T, s in S],
            "ev_dch" => [value.(Coop_model[:ev_dch][e,t,s]) for e in E, t in T, s in S],
            # Reserve Slacks
            "q_c_FCRN_up" => [value.(Coop_model[:q_c_FCRN_up][c,t,s]) for c in C, t in T, s in S],
            "q_c_FCRN_down" => [value.(Coop_model[:q_c_FCRN_down][c,t,s]) for c in C, t in T, s in S],
            "q_c_mFRR" => [value.(Coop_model[:q_c_mFRR][c,t,s]) for c in C, t in T, s in S],
            "q_c_FFR" => [value.(Coop_model[:q_c_FFR][c,t,s]) for c in C, t in T, s in S],
            "q_c_FCRD" => [value.(Coop_model[:q_c_FCRD][c,t,s]) for c in C, t in T, s in S],
            #
            "q_mFRR" => [value.(Coop_model[:q_mFRR][t,s]) for t in T, s in S],
            "q_FFR" => [value.(Coop_model[:q_FFR][t,s]) for t in T, s in S],
            "q_FCRD" => [value.(Coop_model[:q_FCRD][t,s]) for t in T, s in S],
            "q_FCRN_up" => [value.(Coop_model[:q_FCRN_up][t,s]) for t in T, s in S],
            "q_FCRN_down" => [value.(Coop_model[:q_FCRN_down][t,s]) for t in T, s in S],
            # JCC & HJCC only
            "y" => (AM_risktype =="JCC" || AM_risktype =="HJCC" ? (AM_risktype =="HJCC" ? [value.(Coop_model[:y][t,s]) for t in T, s in S] : [value.(Coop_model[:y][s]) for s in S]) : false),
            "y_FFR" => (AM_risktype == "CC" ?  [value.(Coop_model[:y_FFR][s]) for s in S] : false),
            "y_FCRN_up" => (AM_risktype == "CC" ?  [value.(Coop_model[:y_FCRN_up][s]) for s in S] : false),
            "y_FCRN_down" => (AM_risktype == "CC" ?  [value.(Coop_model[:y_FCRN_down][s]) for s in S] : false),
            "y_FCRD" => (AM_risktype == "CC" ?  [value.(Coop_model[:y_FCRD][s]) for s in S] : false),
            "y_mFRR" => (AM_risktype == "CC" ?  [value.(Coop_model[:y_mFRR][s]) for s in S] : false),
            # CVaR
            #
            # OS Slacks
            # Add slack variables
            "slack_TCL_up" => [value.(Coop_model[:slack_TCL_up][f,t,s]) for f in F, t in T,s in S],
            "slack_TCL_down" => [value.(Coop_model[:slack_TCL_down][f,t,s]) for f in F, t in T,s in S],
            "slack_EV_up" => [value.(Coop_model[:slack_EV_up][e,t,s]) for e in E, t in T,s in S],
            "slack_EV_down" => [value.(Coop_model[:slack_EV_down][e,t,s]) for e in E, t in T,s in S],
            )
    return sol_Coop_model  

else
    println("Problem infeasible or unbounded")
    return "INFEASIBLE_OR_UNBOUNDED"
end  
end
