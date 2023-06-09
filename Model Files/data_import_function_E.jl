## Flexibility Service Integration for Citizen Energy Communities
# Functions to import data from "processed data"

using CSV, DataFrames,LinearAlgebra,PowerModels, Statistics, Dates, Random, Distributions,Plots
#using JLD2  #JLD,

#datapath = "G:\\Mit drev\\DTU Sustainable Energy\\F23\\Thesis\\Modelling\\Input Data\\Raw Data\\" 
#data_import_new(datapath,false)

function data_import_function_E(datapath,E,test=false) #Import processed data
    
### Community/General Data ###
    # Model settings
    T = [i for i in 1:24] #Hours in a day
    S = [i for i in 1:364] #Number of scenarios
    #N = [] #1; 2; 3] #Set of inflexble prosumers
    F = [1; 2; 3; 4] #Set of fully flexible prosumers
    #E = [5; 6; 7; 8; 9; 10] #Set of Batteries
    C = union(F,E) # All prosumers

    startIndex = 1 # "2021-01-01 00:00"
    endIndex = startIndex+24*length(S)-1 # S days later
 
    prob_s = ones(length(S))*1/length(S)

    Pmax_EC = 70 #70 [kW] Max Community import/export
    Pmax_c = 20 #20 # [kW]  Max prosumer import/export
    n_demand = ones(length(C))*1 # Number of community members for each consumer type


### Solar profile (normalized) ###
    solar = CSV.read(string(datapath, "PV_peak.csv"), DataFrame; delim=",") # OBS PV data from 2019
    PV = reshape(solar[startIndex:endIndex,3],(24,length(S)))

    PV_peak = [2; 0; 0; 2; 0; 0; 6.4; 6.4; 0; 6.4] # Capacity for each consumer type
    
### Household Demand profile kWh ###
    #demand = load(string(datapath,"DemandScenarios.jld"))["demand"]
    demand_data = CSV.read(string(datapath, "ConsumptionDE35_121.csv"), DataFrame; delim=";",select=[2,5]) 
    demand = reshape(demand_data[startIndex:endIndex,2],(24,length(S))).*2/10e5
    Random.seed!(395);
    D_household = zeros(length(C),24,length(S))
    for c in 1:length(C)
        D_household[c,:,:] = demand[:,rand(1:length(S),length(S))]
    end
### TCL Data ###
    T_a_data = CSV.read(string(datapath,"temperature_1.csv"), dateformat = "yyyy-mm-dd hh", DataFrame;delim=";", decimal='.',select=[1,2]) 
    for m in 2:12
                
        T_a_add = CSV.read(string(datapath,"temperature_$(m).csv"), dateformat = "yyyy-mm-dd hh", DataFrame;delim=";", decimal='.',select=[1,2]) 
        #append!(T_a_data,T_a_add)
        T_a_data = vcat(T_a_data, T_a_add)#, on = :DateTime)
    end
    # Add missing data
    new_row1 = ["01-07-2021 11:00", Float64(18)] #index 4356 = "01-07-2021 12:00"
    insert!.(eachcol(T_a_data), 4356, new_row1)
    
    new_row2 = ["13-07-2021 15:00", Float64(25.3)] 
    insert!.(eachcol(T_a_data), 4648, new_row2)
    
    new_row3 = ["14-09-2021 08:00", Float64(14.5)] 
    insert!.(eachcol(T_a_data), 6153, new_row3)
    
    new_row4 = ["14-09-2021 09:00", Float64(15)] 
    insert!.(eachcol(T_a_data), 6154, new_row4)
    
    new_row5 = ["14-09-2021 10:00", Float64(16)] 
    insert!.(eachcol(T_a_data), 6155, new_row5)
    
    new_row6 = ["24-09-2021 14:00", Float64(14.3)]
    insert!.(eachcol(T_a_data), 6399, new_row6)

    new_row7 = ["30-12-2021 01:00", Float64(2.5)] 
    insert!.(eachcol(T_a_data), 8714, new_row7)

    new_row8 = ["30-12-2021 02:00", Float64(3)]
    insert!.(eachcol(T_a_data), 8715, new_row8)
    
    #T_a_data[startIndex,1]
    #T_a_data[endIndex,1]
    T_a = T_a_data[startIndex:endIndex,2] # Ambient temperature

        Pmax_TCL = 9 #[kW] heatpump power limit (rated electrical power)
    T_opt = 21 # [C]  Optimal indoor temperature
    R_th = [9.655; 5.199; 9.655; 5.199; 0; 0; 0; 0; 0; 0] # Thermal resistance of TCL houses [C/kW]
    C_th = [11.2; 11.2; 16.8; 16.8; 0; 0; 0; 0; 0; 0] # Thermal capacitance of TCL houses [kWh/C]
    η_TCL = 4   # COP for TCL
    beta = 1e-4 # Weight of Temperature penalty in objective
    D_base_TCL = reshape((1 ./(η_TCL.*R_th)).*((η_TCL.*R_th).>0)*transpose(max.(0,T_opt .- T_a)),(length(C),24,length(S))) # CxTxS TCL Baseline demand [kW]
    T_a = reshape(T_a,(24,length(S)))
    
### EV Data ###
    SOC_max = zeros(length(C)) #[0; 0; 0; 0; 50; 50; 50; 50; 50; 50] #maximum capacity of batteries in the problem [kWh]. Must be in same order as listed. The numbers you see are 3 tesla model 3 and 3 vw id.4 pure
    SOC_max[E] .+= 50
    Pmax_ev = 10*ones(length(C)) #[0 0 0 0 10 10 250 7 7 118] #charging power [kW] of the battery. Same models as above. (C), AC charger, one DC charger high each
    Pmax_ev[1:10] = [0 0 0 0 10 10 250 7 7 118]

    EV_con = zeros(length(C),length(T),length(S))
    EV_con[E,T,S] .= 1

    dist_disc = Poisson(3)
    dist_rec = Poisson(4)
 
    Random.seed!(1234);
    for e in E
        hour_disc = 5 .+ rand(dist_disc,length(S))
        hour_rec = 13 .+ rand(dist_rec,length(S))
    for s in S
            if hour_disc[s] < hour_rec[s]
                EV_con[e,hour_disc[s]:min(23,hour_rec[s]),s] *= 0
            end
    end
    end
#plot([hour_disc, hour_rec])
   
    η_c = 0.975 # EV charging efficiency
    η_d = 1.026 # EV discharging efficiency


### Market Data ###
    df_spot_price =  CSV.read(string(datapath,"Elspotprices.csv"), DataFrame;delim=';', decimal=',',select=[2,4]) # Balanceprice for FCRN activation
    spot_price = reshape(Array(df_spot_price[:,2]),(24,365))./1000
    
    ### Tariffs and grid prices in DKK/kWh
    elafgift = 0.9
    tso = 0.049 + 0.061 + 0.00229
    dso_radius_winter = [0.2363,0.2363,0.2363,0.2363,0.2363,0.2363,0.2363,0.2363,0.2363,0.2363,0.2363,0.2363,0.2363,0.2363,0.2363,0.2363,0.2363,0.6307,0.6307,0.6307,0.2363,0.2363,0.2363,0.2363] # https://radiuselnet.dk/wp-content/uploads/Radius-priser-siden-2012_november_2022.pdf?_gl=1*cp6kso*_up*MQ..*_ga*MjExMTI0NzY3OC4xNjgyNjc1OTg1*_ga_TRQXB95J5R*MTY4MjY3NTk4NC4xLjEuMTY4MjY3NTk5MC4wLjAuMA..
    dso_radius_summer = 0.2363
    vindstoed = 0.00375 + 0.000875 + 0.01
    
    #grid out price
    grid_out_price = reshape(spot_price[startIndex:endIndex],(24,length(S))) .- vindstoed

    
    #grid in price
    grid_in_price2 = spot_price
    winter1 = [i for i in 1:(31+28+30)] # jan feb mar
    winter2 = 273 .+[i for i in 1:(31+30+31)] #oct nov dec
    winter = [winter1; winter2]
    summer = 31+28+30 .+ [i for i in 1:3*30+3*31] #apr to sep
    grid_in_price2[:,winter] .+= elafgift .+ tso .+ dso_radius_winter #).*1.25 # tariff added https://energinet.dk/El/Elmarkedet/Tariffer/Aktuelle-tariffer
    grid_in_price2[:,summer] .+=  elafgift + tso + dso_radius_summer #).*1.25 # tariff added https://energinet.dk/El/Elmarkedet/Tariffer/Aktuelle-tariffer
    
    grid_in_price = 1.25.*reshape(grid_in_price2[startIndex:endIndex],(24,length(S)))

    #grid_in_price = grid_in_price2[startIndex:endIndex].*1.25
    #grid_in_price = (spot_price .+ elafgift .+ tso .+ dso_radius_summer).*1.25 # tariff added https://energinet.dk/El/Elmarkedet/Tariffer/Aktuelle-tariffer
    
    
    df_balance_prices =  CSV.read(string(datapath,"RegulatingBalancePowerdata.csv"), DataFrame;delim=';', decimal=',',select=[2,4,9,11,13]) # Balanceprice for FCRN activation
    # Add missing data
    #new_row =  Array(df_balance_prices[2067,:]) 
    #new_row[1,1] = "28-03-2021 02:00" 
    #insert!.(eachcol(df_balance_prices), 2067, new_row)
    
    mFRR_act = (reshape(Array(df_balance_prices[:,2])[startIndex:endIndex],(24,length(S))).>0) #[0~1]
    balance_price = reshape(Array(df_balance_prices[:,3])[startIndex:endIndex],(24,length(S)))./1000 
    balance_up_price = reshape(Array(df_balance_prices[:,4])[startIndex:endIndex],(24,length(S)))./1000 # not used
    balance_down_price = reshape(Array(df_balance_prices[:,5])[startIndex:endIndex],(24,length(S)))./1000  # not used

    df_mFRR =  CSV.read(string(datapath,"MfrrReservesDK2.csv"), DataFrame;delim=';', decimal=',',select=[2,9])
    mFRR_price = reshape(Array(df_mFRR[:,2])[startIndex:endIndex],(24,length(S)))./1000 #2=9 #  [DKK/kWh]
    
    df_FFR_price =  CSV.read(string(datapath,"FFR_price2.csv"), DataFrame;delim=";", decimal='.',select=[3,5]) # FFR reserve price
    # Adding missing data
    for i in 23:-1:0
        new_row = (i>9 ? ["2021-10-12 $(i):00:00", 0] : ["2021-10-12 0$(i):00:00", 0])
        insert!.(eachcol(df_FFR_price), 8617, new_row)
    end
    deleteat!(df_FFR_price,1)
    UTC2_row = ["2022-01-01 00:00:00", 0]
    insert!.(eachcol(df_FFR_price), 8760, UTC2_row)
    FFR_price = reshape(Array(df_FFR_price[:,2])[startIndex:endIndex],(24,length(S)))./1000 # Marginal reserve price for FFR

### Individual prosumer FCR reserve pay-as-bid prices
    FCRN_EV_marginal_price = ones(length(E))*0.30 #0.0561 # [dkk/kW] marginal cost for EV in FCRN
    FCRD_EV_marginal_price = ones(length(E))*0.10 #0.0536 # [dkk/kW] marginal cost for EV in FCRD

    FCRN_TCL_marginal_price = ones(length(F))*0.15 # mean=0.275, 25q = 0.17 75q = 0.34 [dkk/kW] marginal cost for TCL in FCRN # OBS! just a guess...
    FCRD_TCL_marginal_price = ones(length(F))*0.30 # mean = 0.2944, 25 q = 0.20[dkk/kW] marginal cost for TCL in FCRD # OBS! just a guess...

    FCR_prices = CSV.read(string(datapath,"FcrReservesDK2.csv"), DataFrame;delim=";", decimal=',',select=[2,3,5])
    FCRN_price_data = FCR_prices[startIndex:endIndex,2]./1000
    FCRD_price_data = FCR_prices[startIndex:endIndex,3]./1000

    FCRN_price = zeros(length(C))
    FCRN_price[E] = FCRN_EV_marginal_price
    FCRN_price[F] = FCRN_TCL_marginal_price
    FCRD_price = zeros(length(C))
    FCRD_price[E] = FCRD_EV_marginal_price
    FCRD_price[F] = FCRD_TCL_marginal_price

    ### Individual prosumer FCR bid acceptance, based on individual bid prices
    FCRN_EV_accept = FCRN_price_data .> transpose(FCRN_EV_marginal_price)
    FCRD_EV_accept = FCRN_price_data .> transpose(FCRD_EV_marginal_price)

    FCRN_TCL_accept = FCRN_price_data .> transpose(FCRN_TCL_marginal_price)
    FCRD_TCL_accept = FCRD_price_data .> transpose(FCRD_TCL_marginal_price)

    FCRN_accept = zeros(length(C),length(T),length(S))
    FCRD_accept = zeros(length(C),length(T),length(S))

    i = 1
    for e in E
    FCRN_accept[e,:,:] = reshape(Array(FCRN_EV_accept[:,i]),(24,length(S)))
    FCRD_accept[e,:,:] = reshape(Array(FCRD_EV_accept[:,i]),(24,length(S)))
    i+= 1
    end
    i = 1
    for f in F
        FCRN_accept[f,:,:] = reshape(Array(FCRN_TCL_accept[:,i]),(24,length(S)))
        FCRD_accept[f,:,:] = reshape(Array(FCRD_TCL_accept[:,i]),(24,length(S)))
        i+= 1
    end

   
### FCR, FFR activation data based on frequency messurements
    df_act = CSV.read(string(datapath,"freq_data_1.csv"), dateformat = "yyyy-mm-dd hh", DataFrame;delim=",", decimal='.') 
    for m in 2:12
        df_add = CSV.read(string(datapath,"freq_data_$(m).csv"), dateformat = "yyyy-mm-dd hh", DataFrame;delim=",", decimal='.') 
        append!(df_act,df_add)
    end
    new_row =  Array(df_act[2068,:]) 
    new_row[1,1] = "28-03-2021 03:00" 
    insert!.(eachcol(df_act), 2068, new_row)
    
    # UTC+2 -> UTC+1 data
    UTC2_row2 = Array(df_act[1,:])
    UTC2_row2[1,1] = "2022-01-01T00:00:00.0" 
    deleteat!(df_act,1)
    insert!.(eachcol(df_act), 8760, UTC2_row2)
    

    #df_act[startIndex,:Hour]
    #df_act[endIndex,:Hour]
    FCRN_up_act = reshape(df_act[startIndex:endIndex,:FCRN_up_mean],(24,length(S)))
    FCRN_down_act = reshape(df_act[startIndex:endIndex,:FCRN_down_mean],(24,length(S)))
    FCRD_up_act = reshape(df_act[startIndex:endIndex,:FCRD_mean],(24,length(S)))
    FFR_act = reshape(df_act[startIndex:endIndex,:FFR_mean],(24,length(S)))
    
### Choice to remove ancillary services from model
    if test == true
        FFR_price, FCRN_price, mFRR_price, FCRD_price, balance_up_price, balance_down_price = dataload_prices_to_zero(FFR_price, FCRN_price, mFRR_price, FCRD_price, balance_up_price, balance_down_price)
        balance_price = balance_price.*0
        FCRN_accept = FCRN_accept.*0
        FCRN_down_act = FCRN_down_act.*0
        FCRD_accept = FCRD_accept.*0
        FFR_act = FFR_act.*0
    end

### Chance constraint parameters
    epsilon = 0.05 #0.95 # confidence level: 1 - epsilon
    
    # For CVaR
    #data["gamma"] = 0.95 # (to deliver 95% of activation)

    data = Dict(
        "T" => T, #"T_flex_periode" => T_flex_periode, "T_inflex_periode" => T_inflex_periode, 
        "C" => C,"F" => F, "E" => E, "S" => S, "prob_s" => prob_s,
        "Pmax_EC" => Pmax_EC, "SOC_max" => SOC_max, "Pmax_ev" => Pmax_ev, "EV_con" => EV_con,
        "η_c" => η_c, "η_d" => η_d,
        "PV_peak" => PV_peak, "n_demand" => n_demand,
        "PV" => PV, "D_household" => D_household,
        "Pmax_EC" => Pmax_EC, "Pmax_c" => Pmax_c, 
        # TCL
        "Pmax_TCL" => Pmax_TCL, 
        "T_a" => T_a, "T_opt" => T_opt, "R_th" => R_th, "C_th" => C_th, 
        "η_TCL" => η_TCL, "beta" => beta, "D_base_TCL" => D_base_TCL,
        # Prices 
        "grid_in_price" => grid_in_price, "grid_out_price" => grid_out_price,
        "balance_up_price" => balance_up_price, "balance_down_price" => balance_down_price, "balance_price" => balance_price,
        "FFR_price" => FFR_price, "FFR_act" => FFR_act,
        "FCRN_down_act" => FCRN_down_act, "FCRN_up_act" => FCRN_up_act, "FCRN_accept" => FCRN_accept, "FCRN_price" => FCRN_price,
        "FCRN_market_price" => reshape(FCRN_price_data,24,length(S)), "FCRD_market_price" => reshape(FCRD_price_data,24,length(S)),
        #"FCRD_down_act" => FCRD_down_act, 
        "FCRD_up_act" => FCRD_up_act, "FCRD_accept" => FCRD_accept, "FCRD_price" => FCRD_price,
        "mFRR_act" => mFRR_act, "mFRR_price" => mFRR_price,
        "epsilon" => epsilon)
    return(data)

end