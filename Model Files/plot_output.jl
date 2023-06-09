## Functions to plot outputs

using Plots.PlotMeasures
using LaTeXStrings
using PlotlyJS
using StatsPlots, Plots
using Statistics

## To save image, use savefig("destination_path_and_file_name.png")
#include("data_import_function_old.jl")

function plot_EV_dynamics(model,s)
    # model = in_S_results
    # s = 40
    F = data["F"]
    T = data["T"]
    E = data["E"]
    C = data["C"]

    ev_ch_sum = value.([sum(model[:ev_ch][E,t,s]) for t in T])
    #ev_ch_sum = sum(results[:ev_ch][e,:,1] for e in 1:length(E))
    ev_dch_sum = value.([sum(model[:ev_dch][E,t,s]) for t in T])#sum(results[:ev_dch][e,:,1] for e in 1:length(E))
    
    EV_con = sum(data["EV_con"][e,:,s] for e in E)/6*100
    totalB = sum(data["SOC_max"][e] for e in E) 

    SOC_sum2 = value.([sum(model[:SOC][E,t,s]) for t in T])./totalB * 100
    SOC_available = value.([sum(data["EV_con"][E,t,s].*model[:SOC][E,t,s]) for t in T])./totalB * 100
    
    
    Plots.plot([i for i in 1:24], ev_ch_sum,
        title = "Battery charge/discharge (scenario $(s))",
        xlabel = "Hours",
        linewidth =2,
        ylabel = "Power [kW]",
        labels = "Battery Charge",
        xticks = [1, 6, 12, 18, 24],
        legend = :topleft,
        right_margin = 15mm,
        linetype = :steppost) 
    
        Plots.plot!(twinx(),[i for i in 1:24], [SOC_sum2,SOC_available],
        labels = ["SOC" "Available SOC"],
        linewidth = [1 2],
        ylabel = "State of Charge [%]",
        legend = :topright,
        color = ["gray" "black"],
        ylims = (0,100),
        xticks = [1, 6, 12, 18, 24],
        linetype = :steppost)

        Plots.plot!([i for i in 1:24], -ev_dch_sum,
        linewidth =2,
        labels = "Battery Discharge",
        linetype = :steppost)
#=
        Plots.plot!(twinx(),[i for i in 1:24], EV_con,
        labels = "EV con",
        linewidth =2,
        #ylabel = "State of Charge [%]",
        #legend = :topright,
        color = "green",
        ylims = (0,100),
        linetype = :steppost
        )
=#
end

function plot_TCL_demand(model,s)
    #model = out_of_S_results #model_results
    #s = 40
    F = data["F"]
    T = data["T"]
    E = data["E"]
    C = data["C"]
    
    beta = data["beta"]
    T_opt = ones(24).*data["T_opt"]
    T_TCL = 1/length(F)*value.(sum(model[:T_TCL][F,t,s]) for t in T)
    T_a = data["T_a"][:,s]
    #c_TCL = sum(model[:c_TCL][f,:,s] for f in F)
    d_TCL = value.(sum(model[:d_TCL][F,t,s]) for t in T)
    D_base_TCL = sum(data["D_base_TCL"][f,:,s]*data["n_demand"][f] for f in F)
    PV = sum(data["PV_peak"][f]*data["PV"][:,s]*data["n_demand"][f] for f in F)
    
    Plots.plot([i for i in 1:24], (d_TCL),
        legs = (:topleft),
        xlabel = "Hours",
        ylabel = "Power [kW]",
        color = "black",
        label = "TCL Consumption",
        linewidth = 2,
        xticks = [1, 6, 12, 18, 24],
        legend = :topleft)
    Plots.plot!([i for i in 1:24], PV,
        color = "orange",
        label = "PV production",
        linewidth = 1,
        xticks = [1, 6, 12, 18, 24])
    Plots.plot!([i for i in 1:24], D_base_TCL,
        color = "black",
        label = "TCL base demand",
        linewidth = 1,
        linestyle = :dot,
        xticks = [1, 6, 12, 18, 24])
    Plots.plot(twinx(),[i for i in 1:24], T_a,
        xlabel = "",
        ylabel = "Temperature [C]",   
        ylims =  (-5,22),
        color = "lightblue",
        label = "Ambient Temperature",
        linewidth = 1,
        xticks = [1, 6, 12, 18, 24],
        legend = :topright,
        )
end

function plot_TCL_Temp(model,s)
    #s = 40
    F = data["F"]
    T = data["T"]
    E = data["E"]
    C = data["C"]

    beta = data["beta"]
    T_a = data["T_a"][:,s]
    T_opt = ones(24).*data["T_opt"]
    T_TCL = 1/length(F)*value.(sum(model[:T_TCL][F,t,s]) for t in T)
    c_TCL = value.(sum(model[:c_TCL][F,t,s]) for t in T)
    PV = sum(data["PV_peak"][f]*data["PV"][:,s]*data["n_demand"][f] for f in F)

    Plots.plot([i for i in 1:24], T_a,
        xlabel = "Hours",
        ylabel = "Temperature [C]",   
        ylims =  (-5,22),
        color = "lightblue",
        label = "Ambient Temperature",
        linewidth = 1,
        xticks = [1, 6, 12, 18, 24],
        legend = :left,
        right_margin = 18mm)
    Plots.plot!([i for i in 1:24], T_opt,
        color = "black",
        label = "T opt",
        linewidth = 1,
        xticks = [1, 6, 12, 18, 24])
    Plots.plot!([i for i in 1:24], T_TCL,
        color = "black",
        label = "Indoor Temperature",
        linewidth = 2,
        xticks = [1, 6, 12, 18, 24])
    Plots.plot!(twinx(),[i for i in 1:24], c_TCL,
        color = "red",
        label = "Penalty",
        linewidth = 1,
        xticks = [1, 6, 12, 18, 24],
        legend = :right,
        ylabel = "Temperature penalty [DKK]")
    end

function plot_CM_flow(model,s)
    # s = 2
    F = data["F"]
    T = data["T"]
    E = data["E"]
    C = data["C"] 
    
    p_EC = [value.(model[:p_EC][t,s]) for t in T]         #T,S
     p_in = p_EC.*(p_EC.>= 0)
     p_out = p_EC.*(p_EC.<= 0)       #T,S
     p_grid_out = -[value.(model[:p_grid_out][t,s]) for t in T]   #T,S
     p_grid_in = [value.(model[:p_grid_in][t,s]) for t in T]    #T,S

 # Activated Bids 
    p_FFR_up    = -[value.(model[:p_FFR_up][t,s]) for t in T]
    p_FCRD_up   = -[value.(model[:p_FCRD_up][t,s]) for t in T]
    p_FCRN_up   = -[value.(model[:p_FCRN_up][t,s]) for t in T]
    p_FCRN_down = [value.(model[:p_FCRN_down][t,s]) for t in T]
    p_mFRR_up   = -[value.(model[:p_mFRR_up][t,s]) for t in T]

    ymin = -1 + minimum(p_grid_out .+ p_FFR_up .+ p_FCRD_up .+ p_FCRN_up .+ p_mFRR_up)
    ymax = 1 + maximum(p_grid_in .+ p_FCRN_down)

    groupedbar([transpose(p_grid_in) transpose(p_FCRN_down) transpose(p_grid_out) transpose(p_FFR_up) transpose(p_FCRD_up) transpose(p_FCRN_up) transpose(p_mFRR_up)],
        group=repeat(["Import from grid" , "FCR-N down regulation" , "Export to grid" , "FFR up regulation" , "FCR-D up regulation" , "FCR-N up regulation" , "mFRR up regulation"], inner = 24),
        bar_position = :stack,
        xlabel = "Hour (h)",
        ylabel = "Power (kW)",
        legend = :bottomleft,
        ylims = (ymin, ymax),
        ytickfontsize=10,
        xtickfontsize=10,
        guidefontsize = 12,
        c = repeat([5, "rgb(0,0,205)", "rgb(173,255,47)", 3, 4, 1, 2],inner=24),
        xticks = xticks = [1, 6, 12, 18, 24]
    )
    Plots.plot!(p_EC,lw=3,color="black",label="EC Power flow")

end

function plot_prosumer_flow(model,s,c)
    #s = 38     # C,T,S
    #c = 1
    F = data["F"]
    T = data["T"]
    E = data["E"]
    C = data["C"]

    t = [i for i in 1:24]

    p_grid = [value.(model[:p_c_grid][c,t,s]) for t in T]       #C,T,S
    p_c = [model[:p_c][c,t,s] for t in T]         #T,S
      
    # Activated Bids 
   p_FFR_up    = -[value.(model[:p_c_FFR_up][c,t,s] ) for t in T]
   p_FCRD_up   = -[value.(model[:p_c_FCRD_up][c,t,s]) for t in T]
   p_FCRN_up   = -[value.(model[:p_c_FCRN_up][c,t,s]) for t in T]
   p_FCRN_down = [value.(model[:p_c_FCRN_down][c,t,s] ) for t in T]
   p_mFRR_up   = -[value.(model[:p_c_mFRR_up][c,t,s]) for t in T]

   ymin = -1 + minimum(p_grid .+ p_FFR_up .+ p_FCRD_up .+ p_FCRN_up .+ p_mFRR_up)
   ymax = 1 + maximum([p_grid p_FCRN_down])

   groupedbar([transpose(p_grid) transpose(p_FCRN_down) transpose(p_FFR_up) transpose(p_FCRD_up) transpose(p_FCRN_up) transpose(p_mFRR_up)],
       group=repeat(["Import/export from grid" , "FCR-N down regulation" , "FFR up regulation" , "FCR-D up regulation" , "FCR-N up regulation" , "mFRR up regulation"], inner = 24),
       bar_position = :stack,
       xlabel = "Hour (h)",
       ylabel = "Power (kW)",
       legend = :bottomleft,
       ylims = (ymin, ymax),
       ytickfontsize=10,
       xtickfontsize=10,
       guidefontsize = 12,
       #c = repeat([6, "rgb(0,0,205)", "rgb(173,255,47)", 3, 4, 1, 2,5],inner=24),
       xticks = xticks = [1, 6, 12, 18, 24]
   )
   plot!(p_c,lw=3,color="black",label="p_c Power flow")
end



function plot_prices(model,s,S)
# model = out_of_S_results
F = data["F"]
T = data["T"]
E = data["E"]
C = data["C"]
    grid_in_price = data["grid_in_price"][:,s]
    grid_out_price = data["grid_out_price"][:,s]
    balance_price = data["balance_price"][:,s]
    
    #balance_up_price = data["balance_up_price"][:,s]
    #balance_down_price = data["balance_down_price"][:,s]
    
    grid_in_price_avg = sum(model[:pi_s][s]*data["grid_in_price"][:,s] for s in S)
    grid_out_price_avg = sum(model[:pi_s][s]*data["grid_out_price"][:,s] for s in S)
    balance_price_avg = sum(model[:pi_s][s]*data["balance_price"][:,s] for s in S)
    #balance_up_price_avg = sum(model[:pi_s][s]*data["balance_up_price"][:,s] for s in S)
    #balance_down_price_avg = sum(model[:pi_s][s]*data["balance_down_price"][:,s] for s in S)
    
    Plots.plot([i for i in 1:24], grid_in_price,
        xlabel = "Hours",
        ylabel = "Price [DKK/kWh]",
        labels = "Import price",
        linewidth =2,
        color = 5,
        xticks = xticks = [1, 6, 12, 18, 24])
    Plots.plot!([i for i in 1:24], grid_out_price,
        labels = "Export price", 
        linewidth =2, color = "rgb(173,255,47)")
    Plots.plot!([i for i in 1:24], balance_price,
        labels = "Balancing price",
        linewidth =2,
        color = 2)
    Plots.plot!([i for i in 1:24], grid_in_price_avg,
        labels = "Avg. Export price", 
        linewidth =1, color = 5,
        linestyle = :dash)
    Plots.plot!([i for i in 1:24], grid_out_price_avg,
        labels = "Avg. Import price",
        linewidth =1,
        color = "rgb(173,255,47)",
        linestyle = :dash)
    Plots.plot!([i for i in 1:24], balance_price_avg,
        labels = "Avg. Import price",
        linewidth =1,
        color = 2,
        linestyle = :dash)
end

function plot_prices_market(model,s)
    e = 9 # 
    f = 1

    F = data["F"]
    T = data["T"]
    E = data["E"]
    C = data["C"]

    FFR_price = data["FFR_price"][:,s]
    FCRD_price = data["FCRD_price"]#[:,s]
    FCRN_price = data["FCRN_price"]#[:,s]
    mFRR_price = data["mFRR_price"][:,s]
    
    Plots.plot([i for i in 1:24], FFR_price,
        xlabel = "Hours",
        ylabel = "Price [DKK/kWh]",
        labels = "FFR price",
        linewidth =2,
        color = 3,
        xticks = xticks = [1, 6, 12, 18, 24])
    Plots.plot!([i for i in 1:24], data["FCRD_accept"][f,:,s]*FCRD_price[f],
        labels = "FCR-D price TCL", 
        linewidth =2)
    Plots.plot!([i for i in 1:24], data["FCRD_accept"][e,:,s]*FCRD_price[e],
        labels = "FCR-D price EV", 
        linewidth =2)
    Plots.plot!([i for i in 1:24], data["FCRN_accept"][f,:,s]*FCRN_price[f],
        labels = "FCR-N price TCL",
        linewidth =2,
        )
    Plots.plot!([i for i in 1:24], data["FCRN_accept"][e,:,s]*FCRN_price[e],
        labels = "FCR-N price EV",
        linewidth =2,
        )
    Plots.plot!([i for i in 1:24], mFRR_price,
        labels = "mFRR price",
        linewidth =2,
        legend = :topleft)
end

#=
function plot_prices_market_avg(model,s)

    FFR_price_avg = sum(data["FFR_price"][:,s]*model[:pi_s][s] for s in S)
    FCRD_price_avg = data["FCRD_price"] #sum(data["FCRD_price"][:,s]*model[:pi_s][s] for s in S)
    FCRN_price_avg = data["FCRN_price"] #sum(data["FCRN_price"][:,s]*model[:pi_s][s] for s in S)
    mFRR_price_avg = sum(data["mFRR_price"][:,s]*model[:pi_s][s] for s in S)
    
    Plots.plot([i for i in 1:24], FFR_price_avg,
        xlabel = "Hours",
        ylabel = "Price [DKK/kWh]",
        labels = "FFR price",
        linewidth =2,
        color = 3,
        xticks = xticks = [1, 6, 12, 18, 24])
    Plots.plot!([i for i in 1:24], FCRD_price_avg,
        labels = "FCR-D price", 
        linewidth =2, color = 4)
    Plots.plot!([i for i in 1:24], FCRN_price_avg,
        labels = "FCR-N price",
        linewidth =2,
        color = 1)
    Plots.plot!([i for i in 1:24], mFRR_price_avg,
        labels = "mFRR price",
        linewidth =2,
        color = 2,
        legend = :topleft)
end
=#

# Plot the stacked area of the full model - Does not work so far 
function stacked_bar_activation(model, s)
    #model = mr  
    #s = s     # C,T,S
        
    F = data["F"]
    T = data["T"]
    E = data["E"]
    C = data["C"]
    # Time
    #t = [i for i in 1:24]

    # Power flow
    #p_sum = -sum(model[:p_c][c,:,s] for c in C)
    p_sum = -[model[:p_EC][t,s] for t in T]

    # Import / Export from grid
    p_grid_in = -model[:p_grid_in][:,s] 
    p_grid_out = model[:p_grid_out][:,s] 

    # Activated Bids 
    p_FFR_up    = [sum(model[:p_c_FFR_up][C,t,s]) for t in T]
    p_FCRD_up   = [sum(model[:p_c_FCRD_up][C,t,s]) for t in T]
    p_FCRN_up   = [sum(model[:p_c_FCRN_up][C,t,s]) for t in T]
    p_FCRN_down = -[sum(model[:p_c_FCRN_down][C,t,s]) for t in T]
    p_mFRR_up   = [sum(model[:p_c_mFRR_up][C,t,s]) for t in T]

    ymax = 1 + maximum(p_grid_out .+ p_FFR_up .+ p_FCRD_up .+ p_FCRN_up .+ p_mFRR_up)
    ymin = - 1 + minimum(p_grid_in .+ p_FCRN_down)


    groupedbar([transpose(p_grid_in) transpose(p_FCRN_down) transpose(p_grid_out) transpose(p_FFR_up) transpose(p_FCRD_up) transpose(p_FCRN_up) transpose(p_mFRR_up)],
        group=repeat(["Import from grid" , "FCR-N down regulation" ,"Export to grid" , "FFR up regulation" , "FCR-D up regulation" , "FCR-N up regulation" , "mFRR up regulation"], inner = 24),
        bar_position = :stack,
        xlabel = "Hour (h)",
        ylabel = "Power (kW)",
        legend = :bottom,
        ylims = (ymin, ymax),
        ytickfontsize=10,
        xtickfontsize=10,
        guidefontsize = 12,
        c = repeat([5, "rgb(0,0,205)", "rgb(173,255,47)", 3, 4, 1, 2],inner=24),
        xticks = xticks = [1, 6, 12, 18, 24]
        )
    plot!(p_sum,lw=3,color="black",label="Power flow")

end

function stacked_bar_activation_mean(model,S)
    #model = model_results  
        
    F = data["F"]
    T = data["T"]
    E = data["E"]
    C = data["C"]

    # Time
    t = [i for i in 1:24]

    # Power flow
    pi = model[:pi_s]
    p_sum = -[sum(pi[s]*model[:p_EC][t,s] for s in S) for t in T]

    # Import / Export from grid
    p_grid_in = -[sum(pi[s]*model[:p_grid_in][t,s]  for s in S) for t in T]
    p_grid_out = [sum(pi[s]*model[:p_grid_out][t,s]  for s in S) for t in T]

    # Activated Bids 
    p_FFR_up    = [sum(pi[s]*model[:p_c_FFR_up][c,t,s] for c in 1:length(C), s in S) for t in T]
    p_FCRD_up   = [sum(pi[s]*model[:p_c_FCRD_up][c,t,s] for c in 1:length(C), s in S) for t in T]
    p_FCRN_up   = [sum(pi[s]*model[:p_c_FCRN_up][c,t,s] for c in 1:length(C), s in S) for t in T]
    p_FCRN_down = -[sum(pi[s]*model[:p_c_FCRN_down][c,t,s] for c in 1:length(C), s in S) for t in T]
    p_mFRR_up   = [sum(pi[s]*model[:p_c_mFRR_up][c,t,s] for c in 1:length(C), s in S) for t in T]

    ymax = 1 + maximum(p_grid_out .+ p_FFR_up .+ p_FCRD_up .+ p_FCRN_up .+ p_mFRR_up)
    ymin = - 1 + minimum(p_grid_in .+ p_FCRN_down)


    p1 = groupedbar([transpose(p_grid_in) transpose(p_FCRN_down) transpose(p_grid_out) transpose(p_FFR_up) transpose(p_FCRD_up) transpose(p_FCRN_up) transpose(p_mFRR_up)],
        group=repeat(["Import from grid" , "FCR-N down regulation" ,"Export to grid" , "FFR up regulation" , "FCR-D up regulation" , "FCR-N up regulation" , "mFRR up regulation"], inner = 24),
        bar_position = :stack,
        title = "Avg. activation",
        xlabel = "Hour (h)",
        ylabel = "Power (kW)",
        legend = :bottom,
        ylims = (ymin, ymax),
        ytickfontsize=10,
        xtickfontsize=10,
        guidefontsize = 12,
        legendfontsize = 6,
        c = repeat([5, "rgb(0,0,205)", "rgb(173,255,47)", 3, 4, 1, 2],inner=24),
        xticks = xticks = [1, 6, 12, 18, 24]
        )
    plot = plot!(p1,p_sum,lw=3,color="black",label="Power flow")
    return plot
end

# Plot the stacked area of the full model
function stacked_bar_bid(model)
    #model = model_results  
    #s = s     # C,T,S
    # Time
    t = [i for i in 1:24]
    F = data["F"]
    T = data["T"]
    E = data["E"]
    C = data["C"]
    # Activated Bids 
    p_FFR   = model[:p_FFR][t]
    p_FCRD  = model[:p_FCRD][t]
    p_FCRN  = model[:p_FCRN][t]
    p_mFRR   = model[:p_mFRR][t]

    ymin = 0
    ymax = maximum(p_FFR .+ p_FCRD .+ p_FCRN .+ p_mFRR)

    #model[:p_FFR]
    groupedbar([transpose(p_FFR) transpose(p_FCRN) transpose(p_FCRD) transpose(p_mFRR)],

        #bar(cats,vals,colour=colours,group=labels, legend=:topleft)
        
        group=repeat([ "FFR Bid" , "FCR-N Bid" , "FCR-D Bid" , "mFRR Bid"], inner = 24),
        bar_position = :stack,
        xlabel = "Hour (h)",
        ylabel = "Power (kW)",
        legend = :topleft,
        ylims = (ymin,ymax),
        ytickfontsize=10,
        xtickfontsize=10,
        guidefontsize = 12,
        #c = repeat([5, "rgb(0,0,205)", "rgb(173,255,47)", 3, 4, 1, 2],inner=24),
        c = repeat([3,1,4,2],inner=24),
        #c = repeat([3, 1, 4,  "rgb(220,20,60)",2],inner=24),
        xticks = xticks = [1, 6, 12, 18, 24]
        )
end

function stacked_bar_q(model, s)
    #model = model_results  
    #s = 1   # C,T,S
        
    # Time
    T = [i for i in 1:24]
    F = data["F"]
    T = data["T"]
    E = data["E"]
    C = data["C"]
 
    # Unmet bid activation
    q_FFR_up    = sum(model[:q_FFR][t,s] for  t in T)
    q_FCRD_up   = sum(model[:q_FCRD][t,s] for  t in T)
    q_FCRN_up   = sum(model[:q_FCRN_up][t,s] for  t in T)
    #q2_FCRN_up   = sum(model[:q_FCRN_up][t,s] for  t in T)
    q_FCRN_down = sum(model[:q_FCRN_down][t,s] for  t in T)
    #q2_FCRN_down   = sum(model[:q_FCRN_down][t,s] for  t in T)
    q_mFRR_up   = sum(model[:q_mFRR][t,s] for  t in T)

    # Met bid activation
    p_FFR_up    = sum(model[:p_FFR_up][t,s] for  t in T)
    p_FCRD_up   = sum(model[:p_FCRD_up][t,s] for  t in T)
    p_FCRN_up   = sum(model[:p_FCRN_up][t,s] for  t in T)
    p_FCRN_down = sum(model[:p_FCRN_down][t,s] for  t in T)
    p_mFRR_up   = sum(model[:p_mFRR_up][t,s] for  t in T)

    # Activated and Accepted Bids 
    act_FFR    = sum(sum(model[:p_c_FFR][C,t])*data["FFR_act"][t,s] for t in T)
    act_mFRR   = sum(model[:p_mFRR][t]*data["mFRR_act"][t,s] for  t in T)
    act_FCRD   = sum(model[:p_c_FCRD][c,t]*data["FCRD_up_act"][t,s]*data["FCRD_accept"][c,t,s] for c in C,t in T)
    act_FCRN_up   = sum(model[:p_c_FCRN][c,t]*data["FCRN_up_act"][t,s]*data["FCRN_accept"][c,t,s] for c in C,t in T)
    act_FCRN_down = sum(model[:p_c_FCRN][c,t]*data["FCRN_down_act"][t,s]*data["FCRN_accept"][c,t,s] for  c in C,t in T)
    
    #act_FCRD   = sum(model[:p_FCRD][t]*data["FCRD_up_act"][t,s]*data["FCRD_accept"][c,t,s] for c in C,t in T)
    #act_FCRN_up   = sum(model[:p_c_FCRN][c,t]*data["FCRN_up_act"][t,s]*data["FCRN_accept"][C[c],t,s] for c in 1:length(C),t in T)
    #act_FCRN_down = sum(model[:p_c_FCRN][c,t]*data["FCRN_down_act"][t,s]*data["FCRN_accept"][C[c],t,s] for  c in 1:length(C),t in T)
    

    act_bids = [act_FFR,0,
                act_FCRD,0,
                act_FCRN_up,0,
                act_FCRN_down,0,
                act_mFRR,0]
    delivery = [0,p_FFR_up,
                0,p_FCRD_up,
                0,p_FCRN_up,
                0,p_FCRN_down,
                0,p_mFRR_up]
    act_unmet = [0,q_FFR_up,
                 0,q_FCRD_up,
                 0,q_FCRN_up,
                 0,q_FCRN_down,
                 0,q_mFRR_up]
    

    pdata = [delivery act_unmet act_bids] #
    
    ctg = repeat(["Delivered","Missing delivery","Activated bid"], inner = 10)
    nam = repeat(["01","01 FFR_up",
                  "02","02 FCRD_up",
                  "03","03 FCRN_up",
                  "04","04 FCRN_down",
                  "05","05 mFRR_up"], outer = 3)
    groupedbar(nam,pdata, group = ctg, bar_position = :stack, xlabel = "Ancillary Markets", ylabel = "Energy [kWh]",
        title = "CM delivery on activated ancillary bids in sc $(s)", bar_width = 0.67,
        lw = 0, framestyle = :box,xtickfontsize=5,legend=:right)
end

function stacked_bar_q_mean(model,S)
    #model = my_model_results  
       # C,T,S
       F = data["F"]
       T = data["T"]
       E = data["E"]
       C = data["C"]
    # Time
    T = [i for i in 1:24]

    pi_s = model[:pi_s]
    # Unmet bid activation
    q_FFR_up    = sum(pi_s[s]*model[:q_FFR][t,s] for  t in T, s in S)
    q_FCRD_up   = sum(pi_s[s]*model[:q_FCRD][t,s] for  t in T, s in S)
    q_FCRN_up   = sum(pi_s[s]*model[:q_FCRN_up][t,s] for  t in T, s in S)
    q_FCRN_down = sum(pi_s[s]*model[:q_FCRN_down][t,s] for  t in T, s in S)
    q_mFRR_up   = sum(pi_s[s]*model[:q_mFRR][t,s] for  t in T, s in S)
    
    q_sum = q_FFR_up + q_FCRD_up + q_FCRN_up + q_FCRN_down + q_mFRR_up
    
    # Met bid activation
    p_FFR_up    = sum(pi_s[s]*model[:p_FFR_up][t,s] for  t in T, s in S)
    p_FCRD_up   = sum(pi_s[s]*model[:p_FCRD_up][t,s] for  t in T, s in S)
    p_FCRN_up   = sum(pi_s[s]*model[:p_FCRN_up][t,s] for  t in T, s in S)
    p_FCRN_down = sum(pi_s[s]*model[:p_FCRN_down][t,s] for  t in T, s in S)
    p_mFRR_up   = sum(pi_s[s]*model[:p_mFRR_up][t,s] for  t in T, s in S)
    p_sum = p_FFR_up + p_FCRD_up + p_FCRN_up + p_FCRN_down + p_mFRR_up
    

    # Activated and Accepted Bids 
    act_FFR    = sum(pi_s[s]*model[:p_c_FFR][c,t]*data["FFR_act"][t,s] for c in 1:length(C),t in T, s in S)
    act_FCRD   = sum(pi_s[s]*model[:p_c_FCRD][c,t]*data["FCRD_up_act"][t,s]*data["FCRD_accept"][c,t,s] for  c in C,t in T, s in S)
    act_FCRN_up   = sum(pi_s[s]*model[:p_c_FCRN][c,t]*data["FCRN_up_act"][t,s]*data["FCRN_accept"][c,t,s] for  c in C,t in T, s in S)
    act_FCRN_down = sum(pi_s[s]*model[:p_c_FCRN][c,t]*data["FCRN_down_act"][t,s]*data["FCRN_accept"][c,t,s] for  c in C,t in T, s in S)
    act_mFRR   = sum(pi_s[s]*model[:p_c_mFRR][c,t]*data["mFRR_act"][t,s] for  c in C,t in T, s in S)
    act_sum = act_FFR + act_FCRD + act_FCRN_up + act_FCRN_down + act_mFRR
    

    act_bids = [act_FFR, act_FCRD, act_FCRN_up, act_FCRN_down, act_mFRR, act_sum]
    act_bids = act_bids .+ iszero.(act_bids)

    delivery = [p_FFR_up, p_FCRD_up, p_FCRN_up, p_FCRN_down, p_mFRR_up, p_sum]
    delivery = delivery .* (delivery .>10e-5)

    act_unmet = [q_FFR_up, q_FCRD_up, q_FCRN_up, q_FCRN_down, q_mFRR_up, q_sum]
    act_unmet = act_unmet .* (act_unmet .>10e-5)

    # percent delivery and unmet bid of activated bid
    pdata = [act_unmet./act_bids delivery./act_bids]
   # pdata = isnan.(pdata) pdata


    ctg = repeat(["1 unmet act","2 delivery"], inner = 6)
    nam = repeat(["01 FFR",
                  "02 FCRD",
                  "03 FCRN_up",
                  "04 FCRN_down",
                  "05 mFRR_",
                  "06 all_markets"], outer = 2)
    groupedbar(nam,pdata, group = ctg, bar_position = :stack, xlabel = "Ancillary market", ylabel = "Percent of acc. and act. bid",
        title = "CM delivery on activated bids", bar_width = 0.5,
        lw = 0, framestyle = :box,xtickfontsize=9,legend=:right)
    plot!(data["epsilon"].*ones(length(act_unmet)),lw=3,color="black",label="90% delivery")

end

function q_hour_mean(model,S,TT)
    #model = model_results  
       # C,T,S
        
    # Time
    
    F = data["F"]
    E = data["E"]
    C = data["C"]
    T = data["T"]
    no_y = ones(24,length(data["S"]))

    if (TT !== false) # TT == true for HJCC
        no_y = zeros(24,length(data["S"]))
        no_y[TT,S] .+= 1
    end

    pi_s = model[:pi_s]
    # Not delivered
    q_FFR_up    = [sum(pi_s[s]*model[:q_FFR][t,s]*no_y[t,s] for  s in S) for t in T]
    q_FCRD_up   = [sum(pi_s[s]*model[:q_FCRD][t,s]*no_y[t,s] for  s in S) for t in T]
    q_FCRN_up   = [sum(pi_s[s]*model[:q_FCRN_up][t,s]*no_y[t,s] for  s in S) for t in T]
    q_FCRN_down = [sum(pi_s[s]*model[:q_FCRN_down][t,s]*no_y[t,s] for  s in S) for t in T]
    q_mFRR_up   = [sum(pi_s[s]*model[:q_mFRR][t,s]*no_y[t,s] for  s in S) for t in T]
    
    q_sum = q_FFR_up + q_FCRD_up + q_FCRN_up + q_FCRN_down + q_mFRR_up
    
    # Delivered
    p_FFR_up   = [sum([pi_s[s]*sum(model[:p_c_FFR_up][C,t,s])*no_y[t,s] for s in S]) for  t in T]
    p_FCRD_up   = [sum([pi_s[s]*sum(model[:p_c_FCRD_up][C,t,s])*no_y[t,s] for s in S]) for  t in T]
    p_FCRN_up   = [sum([pi_s[s]*sum(model[:p_c_FCRN_up][C,t,s])*no_y[t,s] for s in S]) for  t in T]
    p_FCRN_down   = [sum([pi_s[s]*sum(model[:p_c_FCRN_down][C,t,s])*no_y[t,s] for s in S]) for  t in T]
    p_mFRR_up   = [sum([pi_s[s]*sum(model[:p_c_mFRR_up][C,t,s])*no_y[t,s] for s in S]) for  t in T]
    

    p_sum = p_FFR_up + p_FCRD_up + p_FCRN_up + p_FCRN_down + p_mFRR_up
    

    # Activated and Accepted Bids 
    act_FFR    = [sum(pi_s[s]*model[:p_FFR][t].*data["FFR_act"][t,s]*no_y[t,s] for c in C,s in S) for t in T]
    act_FCRD   = [sum(pi_s[s]*model[:p_c_FCRD][c,t].*data["FCRD_up_act"][t,s].*data["FCRD_accept"][c,t,s]*no_y[t,s] for  c in C,s in S) for t in T]
    act_FCRN_up   = [sum(pi_s[s]*model[:p_c_FCRN][c,t].*data["FCRN_up_act"][t,s].*data["FCRN_accept"][c,t,s]*no_y[t,s] for  c in C,s in S) for t in T]
    act_FCRN_down = [sum(pi_s[s]*model[:p_c_FCRN][c,t].*data["FCRN_down_act"][t,s].*data["FCRN_accept"][c,t,s]*no_y[t,s] for  c in C,s in S) for t in T]
    act_mFRR   = [sum(pi_s[s]*model[:p_mFRR][t].*data["mFRR_act"][t,s]*no_y[t,s] for  c in C,s in S) for t in T]
    act_sum = act_FFR + act_FCRD + act_FCRN_up + act_FCRN_down + act_mFRR
    

    act_bids = [act_FFR, act_FCRD, act_FCRN_up, act_FCRN_down, act_mFRR, act_sum]
    
    delivery = [p_FFR_up, p_FCRD_up, p_FCRN_up, p_FCRN_down, p_mFRR_up, p_sum]
    
    act_unmet = [q_FFR_up, q_FCRD_up, q_FCRN_up, q_FCRN_down, q_mFRR_up, q_sum]
    
    # percent delivery and unmet bid of activated bid
    FFR = Plots.plot([i for i in 1:24], [act_FFR,p_FFR_up,q_FFR_up],
        title = "FFR",    
        xlabel = "Hours",
        ylabel = "kW",
        labels = ["Activation" "Delivered" "Not delivered"],
        legendfontsize = 7,
        titlefontsize = 9,
        xtickfontsize=7,
        ytickfontsize=7,
        guidefontsize=7,
        legend = :left,
        xticks = [1, 6, 12, 18, 24],
        linetype = :steppost)
    FCRD = Plots.plot([i for i in 1:24], [act_FCRD,p_FCRD_up,q_FCRD_up],
        title = "FCRD",    
        xlabel = "Hours",
        ylabel = "kW",
        labels = ["Activation" "Delivered" "Not delivered"],
        legendfontsize = 5,
        titlefontsize = 9,
        xtickfontsize=7,
        ytickfontsize=7,
        guidefontsize=7,
        
        legend = :topright,
        xticks = [1, 6, 12, 18, 24],
        linetype = :steppost)
    FCRN_up = Plots.plot([i for i in 1:24], [act_FCRN_up,p_FCRN_up,q_FCRN_up],
        title = "FCRN up",    
        xlabel = "Hours",
        ylabel = "kW",
        labels = ["Activation" "Delivered" "Not delivered"],
        legendfontsize = 5,
        titlefontsize = 9,
        xtickfontsize=7,
        ytickfontsize=7,
        guidefontsize=7,
        
        legend = :topright,
        xticks = [1, 6, 12, 18, 24],
        linetype = :steppost)
    FCRN_down = Plots.plot([i for i in 1:24], [act_FCRN_down,p_FCRN_down,q_FCRN_down],
        title = "FCRN down",    
        xlabel = "Hours",
        ylabel = "kW",
        labels = ["Activation" "Delivered" "Not delivered"],
        legendfontsize = 5,
        titlefontsize = 9,
        xtickfontsize=7,
        ytickfontsize=7,
        guidefontsize=7,
        
        legend = :topright,
        xticks = [1, 6, 12, 18, 24],
        linetype = :steppost)
    mFRR = Plots.plot([i for i in 1:24], [act_mFRR,p_mFRR_up,q_mFRR_up],
        title = "mFRR",    
        xlabel = "Hours",
        ylabel = "kW",
        labels = ["Activation" "Delivered" "Not delivered"],
        legendfontsize = 5,
        titlefontsize = 9,
        xtickfontsize=7,
        ytickfontsize=7,
        guidefontsize=7,
        
        legend = :topright,
        xticks = [1, 6, 12, 18, 24],
        linetype = :steppost)
    
    all_AM = Plots.plot([i for i in 1:24], [act_sum,p_sum,q_sum],
        title = "all",
        titlefontsize = 9,
        xtickfontsize=7,
        ytickfontsize=7,
        guidefontsize=7,
        
        xlabel = "Hours",
        ylabel = "kW",
        labels = ["Activation" "Delivered" "Not delivered"],
        legendfontsize = 5,
        legend = :topright,
        xticks = [1, 6, 12, 18, 24],
        linetype = :steppost)
    
    Plots.plot(FFR,FCRD,FCRN_up,FCRN_down,mFRR,all_AM, layout = 6)
    
    return FFR,FCRD,FCRN_up,FCRN_down,mFRR,all_AM
end
