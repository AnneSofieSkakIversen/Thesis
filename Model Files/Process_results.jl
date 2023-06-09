using LaTeXStrings
using JuMP, JLD2, Random

include("plot_output_SC.jl")
include("plot_output_SC_avg.jl")

risktype = "no_risk" # []"no_markets" "CC" "JCC" "HJCC" "CVaR"]

#risktype = "no_markets" # no_markets "no_risk", "slack", "JCC", "HJCC", "CVaR"
data = load("const_data.jld2")
OS_folder = "OS_results\\"
Fig_folder = "Figures\\"
inS_size = 260

for risktype in risktypes
IS_name = "model_IS_$(inS_size)_"*risktype
OS_name = "model_OS_$(inS_size)_"*risktype




### IS results
IS_results = load(OS_folder*IS_name*"_results.jld2")
IS_obj_s = IS_results["obj_s_base"] #load(OS_folder*IS_name*"_obj_s.jld2")
IS_obj = IS_results["obj"]
mean(IS_obj_s)
   
 
    Plots.histogram(IS_obj_s, 
    label = "$(inS_size) Scenarios",
    title = "Histogram: "*risktype*" scenario objectives",
    xlabel = "Scenario objective [DKK]",
    xlims = (-300,1300),
    bins = 20,
    ylabel = "Occurances")
    Plots.vline!([IS_obj],linewidth = 4,color = "orange",label ="Mean = $(round(IS_obj,digits = 2))\nSt. dev. = $(round(std(IS_obj_s; corrected=false),digits = 2))")
    png(Fig_folder*risktype*"_$(inS_size)_IS_histogram_Obj_s.png")
    #stacked_bar_bid(IS_results)
    # png(Fig_folder*risktype*"_$(inS_size)_reserve_bid.png")
    stacked_activation(IS_results)
    png(Fig_folder*risktype*"_$(inS_size)_activation.png")
    stacked_bar_q_mean(IS_results,risktype)
    png(Fig_folder*risktype*"_$(inS_size)_q_mean.png")
    
    FFR,FCRD,FCRN_up,FCRN_down,mFRR,all_AM = q_hour_mean(IS_results)
    Plots.plot(FCRD,FCRN_up,FCRN_down,all_AM, layout = 4)
    png(Fig_folder*risktype*"_IS_missing_delivery_hour1.png") ;
    Plots.plot(FFR,mFRR,all_AM, layout = 3)
    png(Fig_folder*risktype*"_IS_missing_delivery_hour2.png") ;

    #plot_prices(IS_results);  png(Fig_folder*risktype*"_balance_prices.png") ; # grid and balance
    #plot_reserve_prices(IS_results);  png(Fig_folder*risktype*"_reserve_prices.png") ; # ancillary reserve prices
    # png(Fig_folder*risktype*"$(inS_size)_reserve_bid.png")

###


### OS results
        OS_results = load(OS_folder*OS_name*"_results.jld2")
        OS_obj_s = OS_results["OS_obj_s"]
        OS_obj = OS_results["obj"]
        mean(OS_obj_s)

               
        Plots.histogram(OS_obj_s, label = "1300 OS Samples",
        title = "OS Histogram: "*risktype*" scenario objectives",
        xlabel = "Scenario objective [DKK]",
        xlims = (-300,2800),
        bins = 100,
        ylabel = "Occurances")
        Plots.vline!([OS_obj],linewidth = 4,color = "orange",label ="OS Mean = $(round(OS_obj,digits = 2)), St. dev. = $(round(std(OS_obj_s;corrected=false),digits = 2))")
        Plots.vline!([IS_obj],linewidth = 3,linestyle = :dash ,color = "orange",label ="IS Mean = $(round(IS_obj,digits = 2)), std. = $(round(std(IS_obj_s; corrected=false),digits = 2))")
        png(Fig_folder*risktype*"_OS_histogram_Obj_s.png")

end       
        

#JuMP.solution_summary(out_of_S_model)
#JuMP.lp_sensitivity_report(out_of_S_model)

#=
### Violations
OS_violations = load(OS_folder*OS_name*"_violations.jld2")
        
total_reserve = sum(IS_results["p_FCRN"] + IS_results["p_FCRD"] + IS_results["p_mFRR"] + IS_results["p_FFR"])
grid_dependancy = sum(IS_results["p_grid_in"])/inS_size

IS_non_delivery_frequency = sum(eachcol(0 .<(IS_results["q_FCRN_up"] .+ IS_results["q_FCRN_down"] .+ IS_results["q_FCRD"] .+ IS_results["q_mFRR"] .+ IS_results["q_FFR"])))
violation_quantity
=#



#### Print information
println("IS "*risktype*":           obj = $(objective_value(in_S_model))")
println("Out of S:                  obj = $(objective_value(out_of_S_model))")
println("Violation info:            TCL up     TCL down    EV up   EV down")
println("no. violations:            $(OS_violations["count"])")
println("total capacity violated:   $(OS_violations["sum"])")
println("Max violation in one hour: $(OS_violations["max"])")




### 

stacked_reserve_bid(IS_results,data["C"])
png(Fig_folder*risktype*"_ALL_reserve_bid.png") ;
stacked_reserve_bid(IS_results,data["F"])
png(Fig_folder*risktype*"_TCL_reserve_bid.png") ;
stacked_reserve_bid(IS_results,data["E"])
#png(folder*risktype*"_EV_reserve_bid.png") ;









#=


OS_name = "model_OS_$(inS_size)_"*risktype
#write_to_file(out_of_S_model, OS_folder*OS_name*".nl")
OS_results = load(OS_folder*OS_name*"_results.jld2")
violations = load(OS_folder*OS_name*"_violations.jld2")

=#