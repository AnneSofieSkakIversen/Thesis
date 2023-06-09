## Functions to plot outputs

using Plots.PlotMeasures
using LaTeXStrings
using PlotlyJS
using StatsPlots, Plots
using Statistics

## To save image, use savefig("destination_path_and_file_name.png")
#include("data_import_function_old.jl")

function plot_EV_dynamics_s(model_results,s)
    SC = model_results["SC"]
    S = [i for i in eachindex(SC)]
    sc = SC[s]
    pi_s = model_results["pi_s"]
    E = data["E"]
    
    #ev_ch_sum = [value.(sum(model_results["ev_ch][E,t,SC[s]])) for t in T]#eachindex(E))
    ev_ch_sum = sum(model_results["ev_ch"][e,:,s] for e in eachindex(E))
    ev_dch_sum = sum(model_results["ev_dch"][e,:,s] for e in eachindex(E))
    
    #EV_con = sum(data["EV_con"][e,:,SC[s]] for e in E)/6*100
    totalB = sum(data["SOC_max"][e] for e in E) 

    SOC_sum2 = sum(model_results["SOC"][e,:,s] for e in eachindex(E))./totalB * 100
    SOC_available = sum(data["EV_con"][E[e],:,SC[s]].*model_results["SOC"][e,:,s] for e in eachindex(E))./totalB * 100
    
    
    Plots.plot([i for i in 1:24], ev_ch_sum,
        title = "Battery charge/discharge (scenario $(sc))",
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
end

function plot_TCL_demand_s(model_results,s)
    SC = model_results["SC"]
    sc = SC[s]
    
    #model_results = model_results
    #s = 25
    
    
    beta = data["beta"]
    T_opt = ones(24).*data["T_opt"]
    T_TCL = 1/length(F)*sum(model_results["T_TCL"][f,:,s] for f in F)
    T_a = data["T_a"][:,sc]
    #c_TCL = sum(model_results["c_TCL][f,:,s] for f in F)
    d_TCL = sum(model_results["d_TCL"][f,:,s] for f in F)
    D_base_TCL = sum(data["D_base_TCL"][f,:,sc]*data["n_demand"][f] for f in F)
    PV = sum(data["PV_peak"][f]*data["PV"][:,sc]*data["n_demand"][f] for f in F)
    
    Plots.plot([i for i in 1:24], (d_TCL),
        legs = (:topleft),
        xlabel = "Hours",
        ylabel = "Power [kW]",
        color = "black",
        label = "TCL Consumption",
        linewidth = 2,
        xticks = [1, 6, 12, 18, 24],
        legend = :topleft,
        right_margin = 18mm)
    Plots.plot!(twinx(),[i for i in 1:24], T_a,
        xlabel = "",
        ylabel = "Temperature [C]",   
        ylims =  (-5,22),
        color = "lightblue",
        label = "Ambient Temperature",
        linewidth = 1,
        xticks = [1, 6, 12, 18, 24],
        legend = :topright,
        )
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
end

function plot_TCL_Temp_s(model_results,s)
    SC = model_results["SC"]
    sc = SC[s]
    
    #model_results = JCC_results
    #s = 25
    
    beta = data["beta"]
    T_a = data["T_a"][:,sc]
    T_opt = ones(24).*data["T_opt"]
    T_TCL = 1/length(F)*sum(model_results["T_TCL"][f,:,s] for f in F)
    c_TCL = sum(model_results["c_TCL"][f,:,s] for f in F)
    PV = sum(data["PV_peak"][f]*data["PV"][:,sc]*data["n_demand"][f] for f in F)

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

function plot_CM_flow_s(model_results,s)
    SC = model_results["SC"]
    sc = SC[s]
    
     #model_results["p_grid]       # C,T,S
     p_EC = - model_results["p_EC"][:,s]         #T,S
     p_in = p_EC.*(p_EC.>= 0)
     p_out = p_EC.*(p_EC.<= 0)       #T,S
     p_grid_out = -model_results["p_grid_out"][:,s]   #T,S
     p_grid_in = model_results["p_grid_in"][:,s]    #T,S

     #p_sum = sum(model_results["p_c][:,s] for c in C)
    # Activated Bids 
    p_FFR_up    = -model_results["p_FFR_up"][:,s]
    p_FCRD_up   = -model_results["p_FCRD_up"][:,s]
    p_FCRN_up   = -model_results["p_FCRN_up"][:,s]
    p_FCRN_down = model_results["p_FCRN_down"][:,s]
    p_mFRR_up   = -model_results["p_mFRR_up"][:,s]

    ymin = 1 + minimum(p_grid_out .+ p_FFR_up .+ p_FCRD_up .+ p_FCRN_up .+ p_mFRR_up)
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
    plot!(p_EC,lw=3,color="black",label="EC Power flow")

end

function plot_prosumer_flow_s(model_results,sc,c)
    SC = model_results["SC"]
    sc = SC[s]
    
    #model_results = model_results  
    #sc = sc     # C,T,S
    t = [i for i in 1:24]

    p_grid = model_results["p_c_grid"][c,:,s]       #C,T,S
    p_c = model_results["p_c"][c,:,s]        #T,S
      
    # Activated Bids 
   p_FFR_up    = -model_results["p_c_FFR_up"][c,:,s] 
   p_FCRD_up   = -model_results["p_c_FCRD_up"][c,:,s]
   p_FCRN_up   = -model_results["p_c_FCRN_up"][c,:,s]
   p_FCRN_down = model_results["p_c_FCRN_down"][c,:,s] 
   p_mFRR_up   = -model_results["p_c_mFRR_up"][c,:,s]

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

function plot_prices_s(model_results,s)
    SC = model_results["SC"]
    sc = SC[s]
    pi_s = model_results["pi_s"]
    
    #
    grid_in_price = data["grid_in_price"][:,sc]
    grid_out_price = data["grid_out_price"][:,sc]
    balance_price = data["balance_price"][:,sc]
    
    #balance_up_price = data["balance_up_price"][:,sc]
    #balance_down_price = data["balance_down_price"][:,sc]
    
    grid_in_price_avg = sum(pi_s[s]*data["grid_in_price"][:,SC[s]] for s in S)
    grid_out_price_avg = sum(pi_s[s]*data["grid_out_price"][:,SC[s]] for s in S)
    #balance_up_price_avg = sum(pi_s[s]*data["balance_up_price"][:,SC[s]] for s in S)
    #balance_down_price_avg = sum(pi_s[s]*data["balance_down_price"][:,SC[s]] for s in S)
    
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
end

function plot_prices_market_s(model_results,s)
    SC = model_results["SC"]
    sc = SC[s]
    
    e = 9 # to get pay-as-bid prices
    f = 1 # to get pay-as-bid prices
    FFR_price = data["FFR_price"][:,sc]
    FCRD_price = data["FCRD_price"]
    FCRN_price = data["FCRN_price"]
    mFRR_price = data["mFRR_price"][:,sc]

    FCRD_accept = data["FCRD_accept"]
    FCRN_accept = data["FCRN_accept"]
    
    Plots.plot([i for i in 1:24], FFR_price,
        xlabel = "Hours",
        ylabel = "Price [DKK/kWh]",
        labels = "FFR price",
        linewidth =2,
        color = 3,
        xticks = xticks = [1, 6, 12, 18, 24])
    Plots.plot!([i for i in 1:24], FCRD_accept[f,:,sc]*FCRD_price[f],
        labels = "FCR-D price TCL", 
        linewidth =2)
    Plots.plot!([i for i in 1:24], FCRD_accept[e,:,sc]*FCRD_price[e],
        labels = "FCR-D price EV", 
        linewidth =2)
    Plots.plot!([i for i in 1:24], FCRN_accept[f,:,sc]*FCRN_price[f],
        labels = "FCR-N price TCL",
        linewidth =2,
        )
    Plots.plot!([i for i in 1:24], FCRN_accept[e,:,sc]*FCRN_price[e],
        labels = "FCR-N price EV",
        linewidth =2,
        )
    Plots.plot!([i for i in 1:24], mFRR_price,
        labels = "mFRR price",
        linewidth =2,
        legend = :topleft)
end

function plot_prices_market_avg_s(model_results,s)
    SC = model_results["SC"]
    S = [i for i in eachindex(SC)]
    sc = SC[s]
    pi_s = model_results["pi_s"]
    
    FFR_price_avg = sum(data["FFR_price"][:,sc]*pi_s[s] for s in SC)
    FCRD_price_avg = data["FCRD_price"] #sum(data["FCRD_price"][:,sc]*pi_s[s] for s in SC)
    FCRN_price_avg = data["FCRN_price"] #sum(data["FCRN_price"][:,sc]*pi_s[s] for s in SC)
    mFRR_price_avg = sum(data["mFRR_price"][:,sc]*pi_s[s] for s in SC)
    
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

function stacked_bar_activation_s(model_results,s)
    SC = model_results["SC"]
    sc = SC[s]
    
    
    #model_results = model_results  
    #sc = sc     # C,T,S
        
    # Time
    t = [i for i in 1:24]

    # Power flow
    #p_sum = -sum(model_results["p_c][c,:,s] for c in C)
    p_sum = -model_results["p_EC"][:,s]

    # Import / Export from grid
    p_grid_in = -model_results["p_grid_in"][:,s] 
    p_grid_out = model_results["p_grid_out"][:,s] 

    # Activated Bids 
    p_FFR_up    = sum(model_results["p_c_FFR_up"][c,:,s] for c in eachindex(C))
    p_FCRD_up   = sum(model_results["p_c_FCRD_up"][c,:,s] for c in eachindex(C))
    p_FCRN_up   = sum(model_results["p_c_FCRN_up"][c,:,s] for c in eachindex(C))
    p_FCRN_down = -sum(model_results["p_c_FCRN_down"][c,:,s] for c in eachindex(C))
    p_mFRR_up   = sum(model_results["p_c_mFRR_up"][c,:,s] for c in eachindex(C))

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

function stacked_bar_q(model_results,s)
    SC = model_results["SC"]
    sc = SC[s]
    S = model_results["S"]
    #model_results = model_results  
    #sc = 1   # C,T,S
        
    # Time
    T = [i for i in 1:24]
    C = model_results["C"]
 
    # Unmet bid activation
    q_FFR_up    = sum(model_results["q_FFR"][t,s] for  t in T)
    q_FCRD_up   = sum(model_results["q_FCRD"][t,s] for  t in T)
    q_FCRN_up   = sum(model_results["q_FCRN_up"][t,s] for  t in T)
    #q2_FCRN_up   = sum(model_results["q_FCRN_up"][t,s] for  t in T)
    q_FCRN_down = sum(model_results["q_FCRN_down"][t,s] for  t in T)
    #q2_FCRN_down   = sum(model_results["q_FCRN_down"][t,s] for  t in T)
    
    q_mFRR_up   = sum(model_results["q_mFRR"][t,s] for  t in T)

    # Met bid activation
    p_FFR_up    = sum(model_results["p_FFR_up"][t,s] for  t in T)
    p_FCRD_up   = sum(model_results["p_FCRD_up"][t,s] for  t in T)
    p_FCRN_up   = sum(model_results["p_FCRN_up"][t,s] for  t in T)
    p_FCRN_down = sum(model_results["p_FCRN_down"][t,s] for  t in T)
    p_mFRR_up   = sum(model_results["p_mFRR_up"][t,s] for  t in T)

    # Activated and Accepted Bids 
    act_FFR    = sum(model_results["p_FFR"][t]*data["FFR_act"][t,sc] for t in T)
    act_FCRD   = sum(model_results["p_FCRD"][t]*data["FCRD_up_act"][t,sc]*data["FCRD_accept"][c,t,sc] for c in C,t in T)
    act_FCRN_up   = sum(model_results["p_c_FCRN"][c,t]*data["FCRN_up_act"][t,sc]*data["FCRN_accept"][C[c],t,sc] for c in eachindex(C),t in T)
    act_FCRN_down = sum(model_results["p_c_FCRN"][c,t]*data["FCRN_down_act"][t,sc]*data["FCRN_accept"][C[c],t,sc] for  c in eachindex(C),t in T)
    act_mFRR   = sum(model_results["p_mFRR"][t]*data["mFRR_act"][t,sc] for  t in T)

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
        title = "CM delivery on activated ancillary bids in sc $(sc)", bar_width = 0.67,
        lw = 0, framestyle = :box,xtickfontsize=5,legend=:right)
end