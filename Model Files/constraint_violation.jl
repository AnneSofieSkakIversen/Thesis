## Basic model 1: Equilibrium model including: TCL, EV, PS and CM, spot(grid) market
using LaTeXStrings
using JuMP, JLD2, Random

### Import Data ############
    begin
    #include("data_import_function_new.jl")
    
    # Set datapath
    #datapath_new = "G:\\Mit drev\\DTU Sustainable Energy\\F23\\Thesis\\Modelling\\Input Data\\Raw Data\\" 
    #data = data_import_function_new(datapath_new,false) # Put to true to set all market prices to zero (except spot)
    data = load("const_data.jld2")
    end
    
    ###################################################
    #### Model ############################
    #begin
    begin
    include("model_functions.jl")
    my_S = [i for i in 1:8:364] # 75 scenarios
    my_model = Model(Gurobi.Optimizer)
    risktype = "JCC"                      # "slack" "no_risk" "CC" "JCC" "HJCC" "CVaR"
    dummy,my_sc,pi_s = Cooperative_Risk_model1(my_model,data,my_S,risktype) 
    mr = solvem(my_model,my_S,risktype)
    end

objective_value(mr)
    T = data["T"]
    C = data["C"]
    



# number of [t,s] where constraint is violated
# 1 means violated constraint

# Ancillary Bid balance p_SA !== sum p_c
sum([mr[:p_mFRR][t] .!== sum(mr[:p_c_mFRR][C,t]) for t in T])
sum([mr[:p_FFR][t] .!== sum(mr[:p_c_FFR][C,t]) for t in T])
sum([mr[:p_FCRD][t] .!== sum(mr[:p_c_FCRD][C,t]) for t in T])
sum([mr[:p_FCRN][t] .!== sum(mr[:p_c_FCRN][C,t]) for t in T])


### SA Reserve slacks q_SA !>= sum q_c
maximum([maximum([norm.(mr[:q_mFRR][t,s] .- sum(mr[:q_c_mFRR][C,t,s])) for t in T]) for s in my_S])
maximum([maximum([norm.(mr[:q_FFR][t,s] .- sum(mr[:q_c_FFR][C,t,s])) for t in T]) for s in my_S])
maximum([maximum([norm.(mr[:q_FCRD][t,s] .- sum(mr[:q_c_FCRD][C,t,s])) for t in T]) for s in my_S])
maximum([maximum([norm.(mr[:q_FCRN_up][t,s] .- sum(mr[:q_c_FCRN_up][C,t,s])) for t in T]) for s in my_S])
maximum([maximum([norm.(mr[:q_FCRN_down][t,s] .- sum(mr[:q_c_FCRN_down][C,t,s])) for t in T]) for s in my_S])


v1 = sum([mr[:q_FFR][t,s] .!== sum(mr[:q_c_FFR][C,t,s]) for t in T, s in my_S])
v2 = sum([mr[:q_FCRD][t,s] .!== sum(mr[:q_c_FCRD][C,t,s]) for t in T, s in my_S])
v3 = sum([mr[:q_FCRN_up][t,s] .!== sum(mr[:q_c_FCRN_up][C,t,s]) for t in T, s in my_S])
v4 = sum([mr[:q_FCRN_down][t,s] .!== sum(mr[:q_c_FCRN_down][C,t,s]) for t in T, s in my_S])
v5 = sum([mr[:q_mFRR][t,s] .!== sum(mr[:q_c_mFRR][C,t,s]) for t in T, s in my_S])

v_index_FFR = findall(!iszero, [mr[:q_FFR][t,s] .< sum(mr[:q_c_FFR][C,t,s]) for t in T, s in my_S]) 
v_index_mFRR = findall(!iszero, [mr[:q_mFRR][t,s] .< sum(mr[:q_c_mFRR][C,t,s]) for t in T, s in my_S]) 
v_index_FCRD = findall(!iszero, [mr[:q_FCRD][t,s] .< sum(mr[:q_c_FCRD][C,t,s]) for t in T, s in my_S]) 
v_index_FCRN_up = findall(!iszero,[mr[:q_FCRN_up][t,s] .< sum(mr[:q_c_FCRN_up][C,t,s]) for t in T, s in my_S]) 
v_index_FCRN_down = findall(!iszero,[mr[:q_FCRN_down][t,s] .< sum(mr[:q_c_FCRN_down][C,t,s]) for t in T, s in my_S])

#a = mr[:q_FCRN_up]
#b = [sum(mr[:q_c_FCRN_up][C,t,s]) for t in T, s in my_S]
#sum((a .- b).^2)

#a = mr[:q_FCRD]
#b = [sum(mr[:q_c_FCRD][C,t,s]) for t in T, s in my_S]
#sum((a .- b).^2)



# JCC problem scenarios: my_S(2,7,9,10)
# CC problem scenarios: no problems? my_S(2,7,9,10,13,34,38,40)

include("plot_output.jl")


stacked_bar_bid(mr)
#png("JCC_activation_s273.png")
stacked_bar_activation(mr,s)
stacked_bar_q(mr,s)
stacked_bar_q_mean(mr,my_S)

plot_TCL_demand(mr,s) 
plot_TCL_Temp(mr,s)
plot_EV_dynamics(mr,s)

#risktype = "JCC"  
begin
    y = mr[:y]
    y_epsilon = sum(y)/size(y)[1]
    no_y_S = my_S[findall(iszero, [mr[:y][s] for s in my_S])]    

stacked_bar_q_mean(mr,my_S)
#png("avg_delivery_amount_JCC_with_Y.png")

FFR,FCRD,FCRN_up,FCRN_down,mFRR,all_AM = q_hour_mean(mr,my_S,false)
#png(Plots.plot(FCRD,FCRN_up,FCRN_down,all_AM, layout = 4),"hour_delivery_JCC_with_Y.png")

FFR2,FCRD2,FCRN_up2,FCRN_down2,mFRR2,all_AM2 = q_hour_mean(mr,no_y_S,false)
#png(Plots.plot(FCRD2,FCRN_up2,FCRN_down2,all_AM2, layout = 4),"hour_delivery_JCC_no_Y.png")
end

#risktype = "HJCC"
begin
    y = mr[:y][T,my_S]
    y_epsilon = sum(y)/(size(y)[1]*size(y)[2])
    
    
    no_y_S = my_S[findall(iszero, [sum(Array(y[T,s])) for s in my_S])]
    no_y_T = findall(iszero, [sum(Array(y[t,my_S])) for t in T])


    stacked_bar_q_mean(mr,my_S)

FFR,FCRD,FCRN_up,FCRN_down,mFRR,all_AM = q_hour_mean(mr,my_S,false)
#png(
    Plots.plot(FCRD,FCRN_up,FCRN_down,all_AM, layout = 4)
#,"hour_delivery_HJCC_with_Y.png")

FFR2,FCRD2,FCRN_up2,FCRN_down2,mFRR2,all_AM2 = q_hour_mean(mr,no_y_S,no_y_T)
Plots.plot(FCRD2,FCRN_up2,FCRN_down2,all_AM2, layout = 4)
end  



Gurobi.GRBgeterrormsg(my_model)

#solution_summary(my_model)
value.(my_model[:q_reserve_FCRD])
dual.(my_model[:q_reserve_FCRD]) # no dual values available
report = lp_sensitivity_report(my_model; atol = 1e-7); #Unable to compute LP sensitivity because model is not a linear program (or it contains interval constraints).

### JCC for Ancillary markets
Pmax_EC = data["Pmax_EC"]
sum([sum(mr[:q_mFRR][T,s]) .> y[s].*3*Pmax_EC for s in my_S])#*3 not 24
sum([sum(mr[:q_FFR][T,s]) .> y[s]*3*Pmax_EC for s in my_S])
sum([sum(mr[:q_FCRD][T,s]) .> y[s]*3*Pmax_EC for s in my_S])
sum([sum(mr[:q_FCRN_up][T,s]) .> y[s]*3*Pmax_EC for s in my_S])
sum([sum(mr[:q_FCRN_down][T,s]) .> y[s]*3*Pmax_EC for s in my_S])

sum(mr[:y][s] for s in my_S) > length(my_S)*(1 - data["epsilon"])
