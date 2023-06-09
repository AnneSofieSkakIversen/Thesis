## Flexibility Service Integration for Citizen Energy Communities
# Functions to import data from "processed data"

### Add to main: ###
using CSV, JLD, DataFrames, Dates, Statistics
datapath = "G:\\Mit drev\\DTU Sustainable Energy\\F23\\Thesis\\Input Data\\Raw Data\\"
#matrix_construct(datapath)
#data_import(datapath)

function data_import(datapath) #Import processed data
    # construct and load scenario matrix
    data_matrix = string(datapath,"data_matrix.csv")
    data_all = CSV.read(data_matrix, DataFrame; delim=",")

    #PV [kWh] and prosumer demand [kWh]
    solar = reshape(data_all[:,"PV"],(24,183*20))
    demand_reshaped = zeros(10,24,183*20) # (C,T,S)
    for c in 1:10
    demand_reshaped[c,:,:] = reshape(data_all[:,string("demand_C",c)],(24,183*20))
    end

    # Spot and balance prices [DKK/kWh]
    spot_price = reshape(data_all[:,"spot_price"],(24,183*20))
    balance_up_price = reshape(data_all[:,"balance_up_price"],(24,183*20)) # tjeck MWh / kWh, use new  tariffs?
    balance_down_price = reshape(data_all[:,"balance_down_price"],(24,183*20))

    # Market data
    FFR_price = reshape(data_all[:,"FFR_price"],(24,183*20))
    FFR_act = reshape(data_all[:,"FFR_act"],(24,183*20))
    FCRN_down_act = reshape(data_all[:,"FCRN_down_act"],(24,183*20))
    FCRN_up_act = reshape(data_all[:,"FCRN_up_act"],(24,183*20))
    FCRN_accept = reshape(data_all[:,"FCRN_accept"],(24,183*20))
    FCRN_price = reshape(data_all[:,"FCRN_price"],(24,183*20))
    FCRD_down_act = reshape(data_all[:,"FCRD_down_act"],(24,183*20))
    FCRD_up_act = reshape(data_all[:,"FCRD_up_act"],(24,183*20))
    FCRD_accept = reshape(data_all[:,"FCRD_accept"],(24,183*20))
    FCRD_price = reshape(data_all[:,"FCRD_price"],(24,183*20))
    mFRR_act = reshape(data_all[:,"mFRR_act"],(24,183*20))
    mFRR_price = reshape(data_all[:,"mFRR_price"],(24,183*20))

    #Model settings
    T = [i for i in 1:24] #Hours in a day
    S = [i for i in 1:183*20] #Number of scenarios (20 prosumer scenarios and 183 market scenarios)
    N = [1; 2; 3] #Set of inflexble prosumers
    F = [4; 5; 6] #Set of fully flexible prosumers
    E = [7; 8; 9; 10] #Set of Batteries
    M = union(F,E) # Prosumers that can participate in mFRR
    IN = union(N,E) #Set of prosumers with inflexible demand
    C = union(N,F,E) # All prosumers
    T_flex_periode = [15 ; 16 ; 17 ; 18 ; 19 ; 20 ; 21]
    T_inflex_periode = [1 ; 2 ; 3 ; 4 ; 5 ; 6 ; 7 ; 8 ; 9 ; 10 ; 11 ; 12 ; 13 ; 14 ; 22 ; 23 ; 24]
    prob_s = ones(length(S))*1/length(S)

    #Grid connection data
    Pmax = 70 # [kW]

    #This data need to be in the the same order as the battery prosumers, E
    E_max = [0; 0; 0; 0; 0; 0; 5; 5; 5; 5] #maximum capacity of batteries in the problem. Must be in same order as listed. The numbers you see are tesla model 3 and vw id.4 pure
    bat_peak = [0; 0; 0; 0; 0; 0; 3; 3; 3; 3] #charging power of the battery. Same models as above.
    η_c = 0.975
    η_d = 1.026

    PV_peak = [2; 0; 0; 2; 0; 0; 6.4; 6.4; 0; 0] # Capacity for each consumer type
    n_demand = ones(10)*1 # Number of community members for each consumer type

    #Tariffs in DKK/kWh
    elafgift = 0.9
    tso = 0.049 + 0.061 + 0.00229
    dso_radius_winter = [0.2363,0.2363,0.2363,0.2363,0.2363,0.2363,0.2363,0.2363,0.2363,0.2363,0.2363,0.2363,0.2363,0.2363,0.2363,0.2363,0.2363,0.6307,0.6307,0.6307,0.2363,0.2363,0.2363,0.2363]
    dso_radius_summer = 0.2363
    vindstoed = 0.00375 + 0.000875 + 0.01

    grid_in_price = (spot_price .+ elafgift .+ tso .+ dso_radius_summer).*1.25 # tariff added https://energinet.dk/El/Elmarkedet/Tariffer/Aktuelle-tariffer
    grid_out_price = spot_price .- vindstoed
    
    data = Dict(
        "S" => S, "T" => T, "T_flex_periode" => T_flex_periode, "T_inflex_periode" => T_inflex_periode,
        "N" => N, "F" => F, "E" => E, "S" => S, "M" => M, "IN" => IN, "prob_s" => prob_s,
        "Pmax" => Pmax, "E_max" => E_max, "bat_peak" => bat_peak, "η_c" => η_c, "η_d" => η_d,
        "PV_peak" => PV_peak, "n_demand" => n_demand,
        "PV" => solar, "D_base" => demand_reshaped,
        "grid_in_price" => grid_in_price, "grid_out_price" => grid_out_price,
        "balance_up_price" => balance_up_price, "balance_down_price" => balance_down_price,
        "FFR_price" => FFR_price, "FFR_act" => FFR_act,
        "FCRN_down_act" => FCRN_down_act, "FCRN_up_act" => FCRN_up_act, "FCRN_accept" => FCRN_accept, "FCRN_price" => FCRN_price,
        "FCRD_down_act" => FCRD_down_act, "FCRD_up_act" => FCRD_up_act, "FCRD_accept" => FCRD_accept, "FCRD_price" => FCRD_price,
        "mFRR_act" => mFRR_act, "mFRR_price" => mFRR_price)
        
    return(data)
end

#function matrix_construct(datapath)
    startIndex = 1 #5088 # "2021-08-01 00:00"
    endIndex = 48 #startIndex+24*31-1 # "2021-08-31 23:00"


    # Solar profile (normalized) and demand [kWh]
    solar = CSV.read(string(datapath, "PV_peak.csv"), DataFrame; delim=",")[:,[1,3]]
    PV = solar[startIndex:endIndex,2]
    # Ambient temperature avg, max and min [C]
    temp = CSV.read(string(datapath, "temperature_1.csv"), DataFrame; delim=";") #temp 2021
    for m in 2:12
        temp_m = CSV.read(string(datapath, "temperature_",m,".csv"), DataFrame; delim=";")
        temp = vcat(temp,temp_m, cols=:union)
    end
     Temp =  temp[startIndex:endIndex,2]
    # Prosumer demand excl. heatpump OBS! adjust
    demand = load(string(datapath,"DemandScenarios.jld"))["demand"]
    demand_reshaped = Array(zeros(10*24*20))
    demand_reshaped = reshape(demand_reshaped,(10,24,20))
    for i in 1:20
        demand_reshaped[:,:,i] = transpose(demand[:,:,i]) #[C,T,S]
    end
    demand_reshaped
    
    #Read demand data from demand.csv file
    df_Spot = CSV.read(string(datapath,"Elspotprices.csv"), DataFrame; delim=';', decimal=',',select=[2,4]) # Spotprice DK2 [DKK/MWh]
    df_Reserve = CSV.read(string(datapath,"MfrrReservesDK2.csv"), DataFrame; delim=';', decimal=',',select=[2,8,9]) # MFRR DK2 Reserve amount [MW] and Price [DKK/MW] b(marginal)
    df_Balance = CSV.read(string(datapath,"RegulatingBalancePowerdata.csv"), DataFrame; delim=';', decimal=',',select=[2,4,5,11,13]) # Activated mFRR Up DK2 [MWh]

    #### Market Data ####
    startIndex = 1 #5088 # "2021-08-01 00:00"
    endIndex = 48 #startIndex+24*31-1 # "2021-08-31 23:00"

    # Spotprice [dkk/kWh]
    spot_price = Array(df_Spot[:,2])[startIndex:endIndex]./1000

    # Balanceprices [dkk/kWh]
    balUp_price = Array(df_Balance[:,4])[startIndex:endIndex]./1000
    balDwn_price = Array(df_Balance[:,5])[startIndex:endIndex]./1000

    # mFRR Activation [kWh/h] and Reserveprice [dkk/kW/h]
    #data_amount_RegUp = Array(df_Balance[:,2])[startIndex:endIndex]*1000
    #data_amount_RegDwn = Array(df_Balance[:,3])[startIndex:endIndex]*1000
    mFRR_price = Array(df_Reserve[:,3])[startIndex:endIndex]./1000           # [DKK/kWh]    
    mFRR_act = Array(df_Balance[:,2])[startIndex:endIndex]*1000 .> 0 # mFRR Up activation [0,1]
    # data_RegDwn = data_amount_RegDwn .< 0 # mFRR Dwn activation [0,1]

    #FCRN, FCRD FFR reserve prices
    df_FCR = CSV.read(string(datapath,"FcrReservesDK2.csv"), DataFrame; delim=';', decimal=',',select=[2,3,5])
    FCRN_price = Array(df_FCR[:,2])[startIndex:endIndex]./1000       # [DKK/kWh]
    FCRD_price = Array(df_FCR[:,3])[startIndex:endIndex]./1000       # FCRD Up reserve[DKK/kWh]

    #OBS! FFR Reserve prices
    #df_FFR_price = CSV.read(string(datapath,"FFR_price.csv"), DataFrame; delim='"', decimal='.')#,select=[2,3,5])
    df_FFR_price = DataFrame(FFR_price=[0]) #OBS! missing data
    FFR_price = zeros(endIndex-startIndex+1)
    
    
    # FCRN, FCRD and FFR Activation from frequency
    freq = CSV.read(string(datapath,"frequency_Q1.csv"), dateformat="yyyy-mm-dd HH:MM:SS", DataFrame; delim=',', decimal='.',select=[3,5])
    for q in 2:4
        freq_q = CSV.read(string(datapath,"frequency_Q",q,".csv"), dateformat="yyyy-mm-dd HH:MM:SS", DataFrame; delim=',', decimal='.',select=[3,5])
        freq = vcat(freq,freq_q, cols=:union)
    end
    rename!(freq,[:hour,:freq])
    df_FCRN_RegUp = DataFrame(freq)
    df_FCRN_RegDwn = DataFrame(freq)
    df_FCRD_RegUp = DataFrame(freq)
    df_FCRD_RegDwn = DataFrame(freq)
    df_FFR_RegUp = DataFrame(freq)

    rename!(df_FCRN_RegUp,[:hour,:FCRN_RegUp])
    rename!(df_FCRN_RegDwn,[:hour,:FCRN_RegDwn])
    rename!(df_FCRD_RegUp,[:hour,:FCRD_RegUp])
    rename!(df_FCRD_RegDwn,[:hour,:FCRD_RegDwn])
    rename!(df_FFR_RegUp,[:hour,:FFR_RegUp])

    act_data = innerjoin(df_FCRN_RegUp,df_FCRN_RegDwn,df_FCRD_RegUp,df_FCRD_RegDwn, df_FFR_RegUp, on=:hour)

    act_data[!,:hour] = DateTime.(floor.(act_data[!,:hour], Dates.Hour))
    act_data[!,:FCRN_RegUp] = act_data[!,:FCRN_RegUp].<50
    act_data[!,:FCRN_RegDwn] = act_data[!,:FCRN_RegDwn].>50
    act_data[!,:FCRD_RegUp] = act_data[!,:FCRD_RegUp].<49.9
    act_data[!,:FCRD_RegDwn] = act_data[!,:FCRD_RegDwn].>50.1
    act_data[!,:FFR_RegUp] = act_data[!,:FFR_RegUp].<49.7
    act_data2 = combine(groupby(act_data,1),[:FCRN_RegUp,:FCRN_RegDwn,:FCRD_RegUp,:FCRD_RegDwn,:FFR_RegUp].=> mean)
    mean.(eachcol(act_data2[!,2:6]))

    FCRN_RegUp = act_data2[!,:FCRN_RegUp_mean][startIndex:endIndex]
    FCRN_RegDwn = act_data2[!,:FCRN_RegDwn_mean][startIndex:endIndex]
    FCRD_RegUp = act_data2[!,:FCRD_RegUp_mean][startIndex:endIndex]
    FCRD_RegDwn = act_data2[!,:FCRD_RegDwn_mean][startIndex:endIndex]
    FFR_RegUp = act_data2[!,:FFR_RegUp_mean][startIndex:endIndex]
    
    #### Matrix construction ####
    T = [i for i in 1:24] #Hours in a day
    S1 = [i for i in 1:20] #Number of prosumer scenarios
    S2 = [i for i in 1:(endIndex-startIndex+1)/24] #Number of energy market scenarios
    S = [i for i in 1:length(S1)*length(S2)] #Total number of scenarios
    
    
    
    
    S2_hour = repeat(T,length(S2)) # S2 hour index
    S2_sc = reshape(transpose(S2).*ones(24),(length(T)*length(S2),1))

    S1_hour = repeat(T,length(S1))
    S1_sc = reshape(transpose(S1).*ones(24),(length(T)*length(S1),1))

    #S2
    #T
    FCRN_accept= ones(length(T)*length(S2))
    FCRD_accept= ones(length(T)*length(S2))
    
    #data_S2 = DataFrame(markets_scenario=[], hour=[],PV=[],temp=[], spot_price=[], balance_up_price=[], balance_down_price=[], FFR_price=[], FFR_act=[], FCRN_down_act=[], FCRN_up_act=[],
    #                    FCRN_accept=[], FCRN_price=[], FCRD_down_act=[], FCRD_up_act=[], FCRD_accept=[], FCRD_price=[], mFRR_act=[], mFRR_price=[])
    S2_col_names = ["markets_scenario", "hour", "PV","Temp", "spot_price", "balance_up_price", "balance_down_price", "FFR_price", "FFR_act", "FCRN_down_act", "FCRN_up_act","FCRN_accept", "FCRN_price", "FCRD_down_act", "FCRD_up_act", "FCRD_accept", "FCRD_price", "mFRR_act", "mFRR_price"]
    data_for_S2 = [S2_sc S2_hour PV Temp spot_price balUp_price balDwn_price FFR_price FCRN_accept FFR_RegUp FCRN_RegDwn FCRN_RegUp FCRN_price FCRD_RegDwn FCRD_RegUp FCRD_price FCRD_price mFRR_act mFRR_price]
    data_S2 = DataFrame(data_for_S2,S2_col_names)

    data_S1 = DataFrame(prosumer_scenario=[],hour=[],demand_C1=[],demand_C2=[],demand_C3=[],demand_C4=[],demand_C5=[],demand_C6=[],demand_C7=[],demand_C8=[],demand_C9=[],demand_C10=[])
    S1_col_names = ["prosumer_scenario","hour","demand_C1","demand_C2","demand_C3","demand_C4","demand_C5","demand_C6","demand_C7","demand_C8","demand_C9","demand_C10"]
    data_for_S1 = [S1_sc S1_hour ]
    data_S1 = DataFrame(data_for_S1,S2_col_names)


    #=
    data_all = DataFrame(markets_scenario=[], prosumer_scenario=[],hour=[],PV=[],temp=[],demand_C1=[],demand_C2=[],demand_C3=[],demand_C4=[],demand_C5=[],demand_C6=[],demand_C7=[],demand_C8=[],demand_C9=[],demand_C10=[],
                                spot_price=[], balance_up_price=[], balance_down_price=[], FFR_price=[], FFR_act=[], FCRN_down_act=[], FCRN_up_act=[],
                                FCRN_accept=[], FCRN_price=[], FCRD_down_act=[], FCRD_up_act=[], FCRD_accept=[], FCRD_price=[], mFRR_act=[], mFRR_price=[])
    
    for s2 in S2
        for s1 in S1
            for t in T
                push!(data_all,
                [s2 s1 t solar[t] temp[t] transpose(demand_reshaped[:,t,s1]) spot_price[t,s2] balance_up_price[t,s2] balance_down_price[t,s2] FFR_price[1,1] FFR_act[t,s2] FCRN_down_act[t,s2] FCRN_up_act[t,s2] FCRN_accept[t,s2] FCRN_price[t,s2] FCRD_down_act[t,s2] FCRD_up_act[t,s2] FCRD_accept[t,s2] FCRD_price[t,s2] mFRR_act[t,s2] mFRR_price[t,s2] ])
            end
        end
    end
    =#
    # every 24 rows contain all information on one scenario. Every datatype has a dedicated column.
    # Data is structured in [row,column] like this: data_all[(s2-1)*20*24 + (s1-1)*24 + t,datatype] 
    # t: hour of the day (1-24)
    # s1: prosumer scenario (1-20)
    # s2: market scenario (1-183)
      
    file = string(datapath,"data_matrix.csv")
    CSV.write(file, data_all)
#end
