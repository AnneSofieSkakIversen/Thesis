function load_parameters(data,SC)
    S = [i for i in 1:length(SC)]
    ### Load Parameters ###
    
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
    

    #    M = 10  #OBS# arbitrary big number # OBS! fix correct individual size
    
    #max_grid_in_price = maximum(maximum(eachrow(grid_in_price)))
    #Lambda_grid_max = 1.2*max_grid_in_price # Theoretical limit on community grid price [dkk/kW]
    #Lambda_max = 10*max_grid_in_price

    # Chance Constraint / CVaR Parameters
    epsilon = data["epsilon"] # confidence level (being 90% sure...)
    #gamma = data["gamma"] # (to deliver 95% of activation)
    
    return T,F,E,C,pi_s,grid_in_price,grid_out_price,balance_price,FFR_price,FFR_act,FCRN_down_act,FCRN_up_act,FCRN_accept,FCRN_price,#
    FCRD_up_act, FCRD_accept, FCRD_price, mFRR_act, mFRR_price, PV, PV_peak, D_household, n_demand, EV_con, η_c, η_d, η_TCL,#
    D_base_TCL, Pmax_EC, Pmax_c, Pmax_TCL, SOC_max, Pmax_ev, R_th, C_th, T_a, T_opt, beta, epsilon#
    T,
    F,
    E,
    C ,
    pi_s,
    grid_in_price,
    grid_out_price,
    balance_price,
    FFR_price,
    FFR_act,
    FCRN_down_act,
    FCRN_up_act ,
    FCRN_accept,
    FCRN_price,
    FCRD_up_act,
    FCRD_accept,
    FCRD_price,
    mFRR_act,
    mFRR_price,
    PV,
    PV_peak,
    D_household,
    n_demand,
    EV_con,
    η_c,
    η_d,
    η_TCL,
    D_base_TCL,
    Pmax_EC,
    Pmax_c,
    Pmax_TCL,
    SOC_max,
    Pmax_ev,
    #TCL Parameters
    R_th,
    C_th,
    T_a,
    T_opt,
    beta,
    # Chance Constraint / CVaR Parameters
    epsilon
end


#T,F,E,C,pi_s,grid_in_price,grid_out_price,balance_price,FFR_price,FFR_act,FCRN_down_act,FCRN_up_act,FCRN_accept,FCRN_price,
#FCRD_up_act, FCRD_accept, FCRD_price, mFRR_act, mFRR_price, PV, PV_peak, D_household, n_demand, EV_con, η_c, η_d, η_TCL,
#D_base_TCL, Pmax_EC, Pmax_c, Pmax_TCL, SOC_max, Pmax_ev, R_th, C_th, T_a, T_opt, beta, epsilon = load_parameters(data,SC)














#=
    ### Sets ###
    T = data["T"] # Number of time instances
    #N = data["N"] # The inflexible prosumers
    F = data["F"] # The fully flexible prosumers
    E = data["E"] # The prosumers with EV
    C = union(F,E) #union(N,F,E) # All prosumers
    #S = data["S"] # Scenarios

    #data["prob_s"] .*= 0 # adjusting sc probability to sample size
    #data["prob_s"][S] .+= 1/length(S)
    pi_s = zeros(length(data["prob_s"]))
    pi_s[S] .+= 1/length(S)
    
    ### Parameters ###
    grid_in_price = data["grid_in_price"]#[:,S] 
    grid_out_price = data["grid_out_price"]#[:,S] 

    #balance_up_price = data["balance_up_price"]#[:,S] 
    #balance_down_price = data["balance_down_price"]#[:,S] 
    balance_price = data["balance_price"]#[:,S]
    
    FFR_price = data["FFR_price"]#[:,S] 
    FFR_act = data["FFR_act"]#[:,S] 
    FCRN_down_act = data["FCRN_down_act"]#[:,S] 
    FCRN_up_act = data["FCRN_up_act"]#[:,S]
    FCRN_accept = data["FCRN_accept"]#[:,:,S]
    FCRN_price = data["FCRN_price"]#[:,:,S] 
    #FCRD_down_act = data["FCRD_down_act"][:,S]
    FCRD_up_act = data["FCRD_up_act"]#[:,S]
    FCRD_accept = data["FCRD_accept"]#[:,:,S]
    FCRD_price = data["FCRD_price"]#[:,:,S]
    mFRR_act = data["mFRR_act"]#[:,S]
    mFRR_price = data["mFRR_price"]#[:,S]

    PV = data["PV"]#[:,S] # Data for PV
    PV_peak = data["PV_peak"] #(C)
    D_household = data["D_household"]#[:,:,S] # Data for household demand for each prosumer (C, T, S)
    n_demand = data["n_demand"] #(C)
    EV_con = data["EV_con"]#[:,:,S] # Binary indicating if EV is connected (T)
    
    η_c = data["η_c"] # Charging efficiency for EV battery
    η_d = data["η_d"] # Discharging efficiency for EV battery
    η_TCL = data["η_TCL"]   # COP for TCL
    D_base_TCL = data["D_base_TCL"]#[:,:,S]

    Pmax_EC = data["Pmax_EC"] # Max power p_in/p_out of energy community
    Pmax_c = data["Pmax_c"] # Max power flow, p_grid of prosumer
    Pmax_TCL = data["Pmax_TCL"] # Max heatpump power d_TCL
    SOC_max = data["SOC_max"] # Max SOC / EV capacity
    Pmax_ev = data["Pmax_ev"] # EV charging/discharging capacity (E)
    
    #TCL Parameters
    R_th = data["R_th"]
    C_th = data["C_th"]
    
    T_a = data["T_a"]#[:,S] # Ambient temperature [C]
    T_opt = data["T_opt"] # optimal indoor temperature [C]
    beta = data["beta"] # Weight of Temperature penalty in objective
    

#    M = 10  #OBS# arbitrary big number # OBS! fix correct individual size
    
    max_grid_in_price = maximum(maximum(eachrow(grid_in_price)))
    Lambda_grid_max = 1.2*max_grid_in_price # Theoretical limit on community grid price [dkk/kW]
    Lambda_max = 10*max_grid_in_price

# Chance Constraint / CVaR Parameters
    epsilon = data["epsilon"] # confidence level (being 90% sure...)
    #gamma = data["gamma"] # (to deliver 95% of activation)

=#