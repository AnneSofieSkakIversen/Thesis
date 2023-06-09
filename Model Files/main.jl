## Basic model 1: Equilibrium model including: TCL, EV, PS and CM, spot(grid) market


### Import Data ############
begin
    using LaTeXStrings
    using JuMP, JLD2, Random
    include("data_import_function_new.jl")
    
    # Set datapath
    
    
    S1 = [s for s in 1:1:364]# choose which Ancillary Market scenarios to use in scenario generation
    S2 = [s for s in 1:1:364] # choose which Prosumer scenarios to use in scenario generation
    l = length(S1)*length(S2)
    
    #datapath_new = "G:\\Mit drev\\DTU Sustainable Energy\\F23\\Thesis\\Modelling\\Input Data\\Raw Data\\" 
    #data_base = data_import_function_new(datapath_new,false) # Put to true to set all market prices to zero (except spot)
    #data = construct_scenarios(data_base,S1,S2)
    #save("const_data.jld2", data)

    data = load("const_data.jld2")
end

##########################results[:#########################
#### Model ############################
begin
    Random.seed!(2349);
    M = length(data["S"]) # number of scenarios to sample from
    N = 260  #260           # number of samples in model
    my_S = rand(1:M,N) # sample of N scenarios
    #data["beta"] = 1

    include("model_functions_bid.jl")
    my_model = Model(Gurobi.Optimizer)
    risktype = "JCC"    # "slack" "no_risk" "CC" "JCC" "HJCC" "CVaR" "HCC"
    dummy,my_sc = Cooperative_Risk_model1(my_model,data,my_S,risktype) 
    my_results = solvem(my_model,my_S,risktype)

end
my_results["obj"]
my_results["obj_tjeck"]
mean(my_results["obj_s_base"])
objective_value(my_model)

###################################################
###--------------- Print outputs ---------------###
###################################################

begin
T = [i for i in 1:24]
obj = my_results["obj"]
grid_in = sum(my_results["pi_s"][s]*my_results["p_grid_in"][t,s] for t in T,s in eachindex(my_S))
grid_out = sum(my_results["pi_s"][s]*my_results["p_grid_out"][t,s] for t in T,s in eachindex(my_S))

println("obj: $(obj), grid_in: $(grid_in), grid_out: $(grid_out)")
end
include("plot_output_SC_avg.jl")


folder = "Risk Bid Results\\"
save(folder*"Risk_bid_"*risktype*"_results.jld2",my_results)

plot_TCL(my_results) #; png(folder*risktype*"_TCL_dynamics_$(data["beta"]).png") ;
plot_EV_dynamics(my_results) #png(folder*risktype*"_EV_dynamics.png") ;

stacked_reserve_bid(my_results,data["C"])
png(folder*risktype*"_ALL_reserve_bid.png") ;
stacked_reserve_bid(my_results,data["F"])
png(folder*risktype*"_TCL_reserve_bid.png") ;
stacked_reserve_bid(my_results,data["E"])
png(folder*risktype*"_EV_reserve_bid.png") ;


stacked_activation(my_results); png(folder*risktype*"_bid_activation.png") ;
#stacked_bar_q(my_results)
stacked_bar_q_mean(my_results,risktype); png(folder*risktype*"_missing_delivery.png") ;


FFR,FCRD,FCRN_up,FCRN_down,mFRR,all_AM = q_hour_mean(my_results)
Plots.plot(FCRD,FCRN_up,FCRN_down,all_AM, layout = 4); png(folder*risktype*"_missing_delivery_hour1.png") ;
Plots.plot(FFR,mFRR,all_AM, layout = 3); png(folder*risktype*"_missing_delivery_hour2.png") ;


plot_prices(my_results); #png(folder*risktype*"_balance_prices.png") ; # grid and balance
plot_reserve_prices(my_results); #png(folder*risktype*"_reserve_prices.png") ; # ancillary reserve prices


#Prosumer overview:
# N: 
#1  PV
#2
#3

# TCL:
#4  PV  TCL
#5      TCL
#6      TCL

# EV:
#7  PV*  EV
#8  PV*  EV
#9      EV
#10     EV




### Saving figures
#mod.write("myModel.lp")

#foldername = "CVaR"
#savefig(plot_ref, filename_string)
#png(plot_ref, filename_string)
# plot(,fmt = :png)

