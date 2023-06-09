### Import Data ############
begin
#begin
    using LaTeXStrings
    using JuMP, JLD2, Random
    include("data_import_function_new.jl")
    
    # Set datapath
    S1 = [s for s in 1:1:364]# choose which Ancillary Market scenarios to use in scenario generation
    S2 = [s for s in 1:1:364] # choose which Prosumer scenarios to use in scenario generation
    l = length(S1)*length(S2)
    
    datapath_new = "G:\\Mit drev\\DTU Sustainable Energy\\F23\\Thesis\\Modelling\\Input Data\\Raw Data\\" 
    data_base = data_import_function_new(datapath_new,false) # Put to true to set all market prices to zero (except spot)
#end

##########################################
### Functions ##########################
function change_FCRN(data_copy,FCRN_price_data,FCRN_price)
    C = data_copy["C"]
    T = data_copy["T"]
    S = data_copy["S"]

    FCRN_accept = zeros(length(C),length(T),length(S))
    for c in C
    FCRN_accept[c,:,:] = FCRN_price_data .> FCRN_price[c]
    end
    data_copy["FCRN_accept"] = FCRN_accept 
    data_copy["FCRN_price"] = FCRN_price
end

#FCRD reserve market price and pay-as-bid prices
function change_FCRD(data_copy,FCRD_price_data,FCRD_price)
    C = data_copy["C"]
    T = data_copy["T"]
    S = data_copy["S"]
    FCRD_accept = zeros(length(C),length(T),length(S))
        for c in C
            FCRD_accept[c,:,:] = FCRD_price_data .> FCRD_price[c]
        end
    data_copy["FCRD_accept"] = FCRD_accept 
    data_copy["FCRD_price"] = FCRD_price
end

function change_balance_price(data_copy,balance_price)
    data_copy["balance_price"] = balance_price
end

function change_temperature(data_copy,T_a)
        η_TCL = data_base["η_TCL"]
        R_th = data_base["R_th"]
        T_opt = data_base["T_opt"]
        C = data_base["C"]
        S = data_base["S"]

        dummy1 = (1 ./(η_TCL.*R_th)).*((η_TCL.*R_th).>0)
        dummy2 = max.(0,T_opt .- T_a)
        D_base_TCL = zeros(length(C),24,length(S))
        for c in C
        D_base_TCL[c,:,:] = dummy1[c].*dummy2
        end
        data_copy["T_a"] = T_a
        data_copy["D_base_TCL"] = D_base_TCL
end

function change_R_th(data_copy,R_th)
        η_TCL = data_base["η_TCL"]
        T_opt = data_base["T_opt"]
        C = data_base["C"]
        S = data_base["S"]
        T_a = data_base["T_a"]
        
        dummy1 = (1 ./(η_TCL.*R_th)).*((η_TCL.*R_th).>0)
        dummy2 = max.(0,T_opt .- T_a)
        D_base_TCL = zeros(length(C),24,length(S))
        for c in C
        D_base_TCL[c,:,:] = dummy1[c].*dummy2
        end
        data_copy["R_th"] = R_th
        data_copy["D_base_TCL"] = D_base_TCL
end

function change_η_TCL(data_copy,η_TCL)
    T_opt = data_base["T_opt"]
    C = data_base["C"]
    S = data_base["S"]
    T_a = data_base["T_a"]
    R_th = data_base["R_th"]
    
    dummy1 = (1 ./(η_TCL.*R_th)).*((η_TCL.*R_th).>0)
    dummy2 = max.(0,T_opt .- T_a)
    D_base_TCL = zeros(length(C),24,length(S))
    for c in C
    D_base_TCL[c,:,:] = dummy1[c].*dummy2
    end
    data_copy["η_TCL"] = η_TCL
    data_copy["D_base_TCL"] = D_base_TCL
end


global function plot_obj_sensitivity(results,name,unit)
    ymin = -1 + minimum(results[name])
    ymax = +1 + maximum(results[name])

    Plots.plot(results["level"],results[name],
        title = risktype*" Sensitivity: "*results["variation"],
        xlabel = "Scalar variation",
        ylabel = unit,
        label = name, 
        ylims = (ymin,ymax),
        linewidth = 2,
        ytickfontsize=10,
        xtickfontsize=10,
        guidefontsize = 12
        )
end

global function plot_bid_sensitivity(results,unit)
    #"Grid import"
    ymax = 1 + maximum([results["FCRN"] results["FFR"] results["mFRR"] results["FCRD"] results["Grid import"]])
    Plots.plot(results["level"],[results["FCRN"] results["FFR"] results["mFRR"] results["FCRD"] results["Grid import"]],
        title = risktype*" Sensitivity: "*results["variation"],
        xlabel = "Scalar variation",
        ylabel = unit,
        legend =:right,
        linewidth = 2,
        labels = ["FCRN" "FFR" "mFRR" "FCRD" "Grid import"],
        color = [1 3 2 4 5],
        ylims = (0,ymax),
        ytickfontsize=10,
        xtickfontsize=10,
        guidefontsize = 12
        )
end

##############################################
#### Data Changes ############################
Random.seed!(2349);
    M = l #length(data["S"]) # number of scenarios to sample from
    N = 260             # number of samples in model
    my_S = rand(1:M,N) # sample of N scenarios
folder = "Sensitivity_test\\"
include("model_functions.jl")

risktype = "CVaR"    # "no markets" "slack" "no_risk" "CC" "JCC" "HJCC" "CVaR" "HCVaR"
saving = false
variation = "Ambient temperature" # "FCRN bid price", "FCRD market price","FCRD bid price","Balance price"

if (variation == "FCRN market price")
        level = [i for i in 0.4:0.1:1.2] # OK
elseif (variation == "FCRN bid price")
        level = [i for i in 0:0.2:1.4] # OK
elseif (variation == "FCRD market price")
        level = [i for i in 0:1:3] # OK - no change
elseif (variation == "FCRD bid price")
        level = [i for i in 1:1:6] # OK
elseif (variation == "FFR price")
        level = [i for i in 6:2:12] # OK
elseif (variation == "mFRR price")
        level = [i for i in 0:2:8] # OK
elseif (variation == "Balance price")
        level = [i for i in 0.5:0.5:1.5] # correlated with spot price!
elseif (variation == "Ambient temperature")
        level = [i for i in 0.6:0.2:1.4]
elseif (variation == "R_th")
        level = [i for i in 0.6:0.1:1.4]
elseif (variation == "C_th")
        level = [i for i in 0.6:0.1:1.4]
elseif (variation == "C_th")
        level = [i for i in 0.6:0.1:1.4]
elseif (variation == "EV cars")
        level = [i for i in 2:2:18] # number of cars
end

#volatility = [i for i in 0.6:0.1:1.4]

#begin
no_runs = length(level)
C = data_base["C"]
objectives = zeros(no_runs)
FCRN_bids = zeros(no_runs)
FCRD_bids = zeros(no_runs)
FFR_bids = zeros(no_runs)
mFRR_bids = zeros(no_runs)
grid_import = zeros(no_runs)

for i in eachindex(level)
    data_copy = copy(data_base)
    if variation == "FCRN market price"
        FCRN_market_price = level[i].*data_base["FCRN_market_price"]
        change_FCRN(data_copy,FCRN_market_price,data_base["FCRN_price"])
    elseif variation == "FCRN bid price"
        FCRN_price = level[i].*data_base["FCRN_price"]
        change_FCRN(data_copy,data_base["FCRN_market_price"],FCRN_price)
    #
    elseif variation == "FCRD market price"
            FCRD_market_price = level[i].*data_base["FCRD_market_price"]
            change_FCRD(data_copy,FCRD_market_price,data_base["FCRD_price"])
    elseif variation == "FCRD bid price"
            FCRD_price = level[i].*data_base["FCRD_price"]
            change_FCRD(data_copy,data_base["FCRD_market_price"],FCRD_price)
    #
    elseif variation == "FFR price"
    data_copy["FFR_price"] = level[i].*data_base["FFR_price"]
    #
    elseif variation == "mFRR price"
    data_copy["mFRR_price"] = level[i].*data_base["mFRR_price"]
    #
    elseif variation == "Balance price"
        balance_price = level[i].*data_base["balance_price"]
        change_balance_price(data_copy,balance_price)
    #
    elseif variation == "Ambient temperature"
        T_a = level[i].*data_base["T_a"]
        change_temperature(data_copy,T_a)
    elseif variation == "R_th"
        R_th = level[i].*data_base["R_th"]
        change_R_th(data_copy,R_th)
    elseif variation == "C_th"
        data_copy["C_th"] = level[i].*data_base["C_th"]
    elseif variation == "η_TCL"
        η_TCL = level[i].*data_base["η_TCL"]
        change_η_TCL(data_copy,η_TCL)
    elseif variation == "EV cars"
        include("data_import_function_E.jl")
        E = [e for e in 5:(4+i)]
        data_copy = data_import_function_E(datapath_new,E,false)
    else
        break
    end 
    global data = construct_scenarios(data_copy,S1,S2)
    
    my_model = Model(Gurobi.Optimizer)
    dummy,my_sc = Cooperative_Risk_model1(my_model,data,my_S,risktype) 
    my_results = solvem(my_model,my_S,risktype)

    objectives[i] = objective_value(my_model)
    FCRN_bids[i] = sum(my_results["p_c_FCRN"][c,t] for c in C, t in 1:24)
    FCRD_bids[i] = sum(my_results["p_c_FCRD"][c,t] for c in C, t in 1:24)
    FFR_bids[i] = sum(my_results["p_c_FFR"][c,t] for c in C, t in 1:24)
    mFRR_bids[i] = sum(my_results["p_c_mFRR"][c,t] for c in C, t in 1:24)
    grid_import[i] = sum(1/N*my_results["p_grid_in"][t,s] for t in 1:24,s in eachindex(my_S))
end
global sensitivity_results = Dict("variation" => variation,
                            "level" => level,
                            "Objective" => objectives,
                            "FCRN" => FCRN_bids,
                            "FCRD" => FCRD_bids,
                            "FFR" => FFR_bids,
                            "mFRR" => mFRR_bids,
                            "Grid import" => grid_import
                            )
end

if saving == true
save(folder*risktype*"_"*variation*"_sensitivity.jld2",sensitivity_results)

plot_obj_sensitivity(sensitivity_results,"Objective","DKK") 
png(folder*risktype*"_"*variation*"_obj.png")

plot_bid_sensitivity(sensitivity_results,"kW")
png(folder*risktype*"_"*variation*"_bids.png")
end
