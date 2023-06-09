## Functions to plot outputs

using Plots.PlotMeasures
using LaTeXStrings
using PlotlyJS
using StatsPlots, Plots
using Statistics

## To save image, use savefig("destination_path_and_file_name.png")
#include("data_import_function_old.jl")

function plot_EV_dynamics(model_results)
    #model_results = my_results
    SC = model_results["SC"]
    S = model_results["S"]
    pi_s = model_results["pi_s"]
    E = data["E"]
    
    #ev_ch_sum = [value.(sum(model_results["ev_ch][E,t,SC[s]])) for t in T]#eachindex(E))
    ev_ch_sum = sum(pi_s[s]*model_results["ev_ch"][e,:,s] for e in eachindex(E),s in S)
    ev_dch_sum = sum(pi_s[s]*model_results["ev_dch"][e,:,s] for e in eachindex(E),s in S)
    
    totalB = sum(data["SOC_max"][e] for e in E) 
    EV_cap = sum(pi_s[s]*data["EV_con"][e,:,SC[s]]*data["SOC_max"][e] for e in E,s in S)./totalB*100
    
    SOC_sum2 = sum(pi_s[s]*model_results["SOC"][e,:,s] for e in eachindex(E), s in S)./totalB * 100
    SOC_available = sum(pi_s[s]*data["EV_con"][E[e],:,SC[s]].*model_results["SOC"][e,:,s] for e in eachindex(E), s in S)./totalB * 100
    
    
    Plots.plot([i for i in 1:24], ev_ch_sum,
        title = "EV Battery Dynamics (avg.)",
        xlabel = "Hours",
        linewidth =2,
        ylabel = "Power [kW]",
        labels = "EV Charge",
        xticks = [1, 6, 12, 18, 24],
        legend = :topleft,
        right_margin = 15mm,
        linetype = :steppost) 
    
        Plots.plot!(twinx(),[i for i in 1:24], [EV_cap,SOC_available],
        labels = ["Available capacity" "SOC"],
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

function plot_TCL(model_results)
     #model_results = v_results
   
    SC = model_results["SC"]
    S = model_results["S"]
    pi_s = model_results["pi_s"]
    F = data["F"]

    T_max = [maximum(b) for b in eachrow(data["T_a"][:,SC])]
    T_min = [minimum(b) for b in eachrow(data["T_a"][:,SC])]
    y2max = maximum(T_max)
    y2min = minimum(T_min)

    
    beta = data["beta"]
    T_avg = sum(pi_s[s]*data["T_a"][:,SC[s]] for s in S)
    T_opt = ones(24).*data["T_opt"]
    T_TCL = 1/length(F)*sum(pi_s[s]*model_results["T_TCL"][f,:,s] for f in F,s in S)
    c_TCL = sum(pi_s[s]*model_results["c_TCL"][f,:,s] for f in F,s in S)
    PV = sum(pi_s[s]*data["PV_peak"][f]*data["PV"][:,SC[s]]*data["n_demand"][f] for f in F,s in S)

    d_TCL = sum(pi_s[s]*model_results["d_TCL"][f,:,s] for f in F,s in S)
    D_base_TCL = sum(pi_s[s]*data["D_base_TCL"][f,:,SC[s]]*data["n_demand"][f] for f in F,s in S)
 
    Plots.plot([i for i in 1:24], [D_base_TCL,d_TCL],
        title = "TCL avg. Consumption and Temperature",
        legs = (:topleft),
        xlabel = "Hours",
        ylabel = "Power [kW]",
        color = "black",
        label = ["TCL base demand" "TCL Consumption"],
        linewidth = [1 2],
        linestyle = [:dash :solid],
        xticks = [1, 6, 12, 18, 24],
        legend = :topleft,
        right_margin = 2mm)
    Plots.plot!(twinx(),[i for i in 1:24], [T_opt,T_TCL],
        color = "blue",
        label = ["T opt" "Indoor Temperature"],
        ylabel = "Temperature [C]",
        linewidth = [1 2],
        ylims = (18,23),
        linestyle = [:dash :solid],
        legend = :topright,
        xticks = [1, 6, 12, 18, 24])
end

function plot_prices(model_results)
    SC = model_results["SC"]
    pi_s = model_results["pi_s"]
    S = model_results["S"]
    #
    #grid_in_price = data["grid_in_price"][:,SC[s]]
    #grid_out_price = data["grid_out_price"][:,SC[s]]
    #balance_price = data["balance_price"][:,SC[s]]
        
    grid_in_price_avg = sum(pi_s[s]*data["grid_in_price"][:,SC[s]] for s in S)
    grid_out_price_avg = sum(pi_s[s]*data["grid_out_price"][:,SC[s]] for s in S)
    balance_price_avg = sum(pi_s[s]*data["balance_price"][:,SC[s]] for s in S)
    

    Plots.plot([i for i in 1:24], [grid_in_price_avg],
        title = "Avg. spot and balance prices",
        xlabel = "Hours",
        ylabel = "Price [DKK/kWh]",
        labels = "Avg. Import price",
        linewidth =2,
        color = 5,
        #linestyle = :solid,
        xticks = xticks = [1, 6, 12, 18, 24])
    Plots.plot!([i for i in 1:24], grid_out_price_avg,
        labels = "Export price", 
        linewidth =2, 
        #linestyle = :dash,
        color = "rgb(173,255,47)")
    Plots.plot!([i for i in 1:24], balance_price_avg,
        labels = "Balancing price",
        linewidth =2,
        #linestyle = :dash,
        color = "black")
end

function plot_reserve_prices(model_results)
    #model_results = my_results
    SC = model_results["SC"]
    S = model_results["S"]
    pi_s = model_results["pi_s"]
    #sc = SC[s]
    
    e = 9 # to get pay-as-bid prices
    f = 1 # to get pay-as-bid prices
    avg_FFR_price = sum(pi_s[s]*data["FFR_price"][:,SC[s]] for s in S)
    avg_mFRR_price = sum(pi_s[s]*data["mFRR_price"][:,SC[s]] for s in S)
    avg_FCRD_price_TCL = sum(pi_s[s]*data["FCRD_accept"][f,:,SC[s]]*data["FCRD_price"][f] for s in S)
    avg_FCRN_price_TCL = sum(pi_s[s]*data["FCRN_accept"][f,:,SC[s]]*data["FCRN_price"][f] for s in S)
    avg_FCRD_price_EV = sum(pi_s[s]*data["FCRD_accept"][e,:,SC[s]]*data["FCRD_price"][e] for s in S)
    avg_FCRN_price_EV = sum(pi_s[s]*data["FCRN_accept"][e,:,SC[s]]*data["FCRN_price"][e] for s in S)
    
    avg_FCRD_market_price = sum(pi_s[s]*data["FCRD_market_price"][:,SC[s]] for s in S)
    avg_FCRN_market_price = sum(pi_s[s]*data["FCRN_market_price"][:,SC[s]] for s in S)
    
    Plots.plot([i for i in 1:24], avg_FFR_price,
        xlabel = "Hours",
        ylabel = "Price [DKK/kW]",
        labels = "FFR price",
        linewidth =2,
        color = 3,
        linetype = :steppost,
        xticks = xticks = [1, 6, 12, 18, 24])
    Plots.plot!([i for i in 1:24], [avg_FCRD_market_price,avg_FCRD_price_EV,avg_FCRD_price_TCL],
        labels = ["FCR-D Market price" "FCR-D price EV" "FCR-D price TCL"], 
        linetype = :steppost,
        linestyle = [:dash :solid :solid],
        color = [4 4 "pink"],
        linewidth =2)
    Plots.plot!([i for i in 1:24], [avg_FCRN_market_price,avg_FCRN_price_EV,avg_FCRN_price_TCL],
        labels = ["FCR-N Market price" "FCR-N price EV" "FCR-N price TCL"],
        linestyle = [:dash :solid :solid],
        linewidth =2,
        color = [1 1 "lightblue"],
        linetype = :steppost)
    Plots.plot!([i for i in 1:24], avg_mFRR_price,
        labels = "mFRR price",
        linewidth =2,
        color = 2,
        linetype = :steppost,
        legend = :topleft)
end

# Plot the stacked area of the full model - Does not work so far 
function stacked_activation(model_results)
    #model_results = my_results
    SC = model_results["SC"]
    S = model_results["S"]
    
    pi_s = model_results["pi_s"]
    
    #model_results = model_results  
    #sc = sc     # C,T,S
        
    # Time
    t = [i for i in 1:24]
    C = model_results["C"]
    # Power flow
    p_sum = -sum(pi_s[s]*model_results["p_EC"][:,s] for s in S)

    # Import / Export from grid
    p_grid_in = -sum(pi_s[s]*model_results["p_grid_in"][:,s] for s in S)
    p_grid_out = sum(pi_s[s]*model_results["p_grid_out"][:,s] for s in S)

    # Activated Bids 
    p_FFR_up    = sum(pi_s[s]*model_results["p_c_FFR_up"][c,:,s] for c in eachindex(C),s in S)
    p_FCRD_up   = sum(pi_s[s]*model_results["p_c_FCRD_up"][c,:,s] for c in eachindex(C),s in S)
    p_FCRN_up   = sum(pi_s[s]*model_results["p_c_FCRN_up"][c,:,s] for c in eachindex(C),s in S)
    p_FCRN_down = -sum(pi_s[s]*model_results["p_c_FCRN_down"][c,:,s] for c in eachindex(C),s in S)
    p_mFRR_up   = sum(pi_s[s]*model_results["p_c_mFRR_up"][c,:,s] for c in eachindex(C),s in S)

    ymax = 1 + maximum(p_grid_out .+ p_FFR_up .+ p_FCRD_up .+ p_FCRN_up .+ p_mFRR_up)
    ymin = - 1 + minimum(p_grid_in .+ p_FCRN_down)


    groupedbar([transpose(p_grid_in) transpose(p_FCRN_down) transpose(p_grid_out) transpose(p_FFR_up) transpose(p_FCRD_up) transpose(p_FCRN_up) transpose(p_mFRR_up)],
        title = "Avg. Activation",
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

# Plot the stacked area of the full model
function stacked_reserve_bid(model_results,C) # choose C=F or C=E for TCL or EV only.
    #SC = model_results["SC"]
    #S = model_results["S"]
    #model_results = model_results  
    #sc = sc     # C,T,S
    # Time
    t = [i for i in 1:24]
    # Activated Bids 
    #p_FFR   = model_results["p_FFR]
    #p_FCRD  = model_results["p_FCRD]
    #p_FCRN  = model_results["p_FCRN]
    #p_mFRR   = model_results["p_mFRR]

    p_FFR   = sum(model_results["p_c_FFR"][c,:] for c in C)
    p_FCRD  = sum(model_results["p_c_FCRD"][c,:] for c in C)
    p_FCRN  = sum(model_results["p_c_FCRN"][c,:] for c in C)
    p_mFRR   = sum(model_results["p_c_mFRR"][c,:] for c in C)

    ymin = 0
    ymax = 1 + maximum(p_FFR .+ p_FCRD .+ p_FCRN .+ p_mFRR)

    #model_results["p_FFR]
    groupedbar([transpose(p_FFR) transpose(p_FCRN) transpose(p_FCRD) transpose(p_mFRR)],
        title = (length(C)<10 ? (length(C)<5 ? "TCL Reserve bids" : "EV Reserve bids") : "Community Reserve bids"),
        #bar(cats,vals,colour=colours,group=labels, legend=:topleft),
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

function stacked_bar_q_mean(model_results,risktype)
    # model_results = IS_results  
    SC = model_results["SC"]
    S = model_results["S"]
    pi_s = model_results["pi_s"]
   
        
    # Time
    T = [i for i in 1:24]
    C = model_results["C"]
 
    # Unmet bid activation
    q_FFR_up    = sum(pi_s[s]*model_results["q_FFR"][t,s] for  t in T,s in S)
    q_FCRD_up   = sum(pi_s[s]*model_results["q_FCRD"][t,s] for  t in T,s in S)
    q_FCRN_up   = sum(pi_s[s]*model_results["q_FCRN_up"][t,s] for  t in T,s in S)
    q_FCRN_down = sum(pi_s[s]*model_results["q_FCRN_down"][t,s] for  t in T,s in S)
    q_mFRR_up   = sum(pi_s[s]*model_results["q_mFRR"][t,s] for  t in T,s in S)
    q_sum = q_FFR_up + q_FCRD_up + q_FCRN_up + q_FCRN_down + q_mFRR_up

    # Met bid activation
    p_FFR_up    = sum(pi_s[s]*model_results["p_FFR_up"][t,s] for  t in T,s in S)
    p_FCRD_up   = sum(pi_s[s]*model_results["p_FCRD_up"][t,s] for  t in T,s in S)
    p_FCRN_up   = sum(pi_s[s]*model_results["p_FCRN_up"][t,s] for  t in T,s in S)
    p_FCRN_down = sum(pi_s[s]*model_results["p_FCRN_down"][t,s] for  t in T,s in S)
    p_mFRR_up   = sum(pi_s[s]*model_results["p_mFRR_up"][t,s] for  t in T,s in S)
    p_sum = p_FFR_up + p_FCRD_up + p_FCRN_up + p_FCRN_down + p_mFRR_up
    

    # Activated and Accepted Bids 
    act_FFR    = sum(pi_s[s]*model_results["p_FFR"][t]*data["FFR_act"][t,SC[s]] for t in T,s in S)
    act_FCRD   = sum(pi_s[s]*model_results["p_FCRD"][t]*data["FCRD_up_act"][t,SC[s]]*data["FCRD_accept"][c,t,SC[s]] for c in C,t in T,s in S)
    act_FCRN_up   = sum(pi_s[s]*model_results["p_c_FCRN"][c,t]*data["FCRN_up_act"][t,SC[s]]*data["FCRN_accept"][C[c],t,SC[s]] for c in eachindex(C),t in T,s in S)
    act_FCRN_down = sum(pi_s[s]*model_results["p_c_FCRN"][c,t]*data["FCRN_down_act"][t,SC[s]]*data["FCRN_accept"][C[c],t,SC[s]] for  c in eachindex(C),t in T,s in S)
    act_mFRR   = sum(pi_s[s]*model_results["p_mFRR"][t]*data["mFRR_act"][t,SC[s]] for  t in T,s in S)
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


    ctg = repeat(["1 Not delivered","2 Delivered"], inner = 6)
    nam = repeat(["FFR",
                  "FCRD",
                  "FCRN up",
                  "FCRN down",
                  "mFRR",
                  "All"], outer = 2)
    groupedbar(nam,pdata, group = ctg, bar_position = :stack, 
        xlabel = "Ancillary service", 
        ylabel = "Percent of activated bid",
        title = risktype*": Expected delivery", bar_width = 0.5,
        lw = 0, framestyle = :box,xtickfontsize=9,legend=:right,
        yticks = ([0, 0.25, 0.5, 0.75, 1],["0%", "25%", "50%", "75%", "100%"]),
        c = repeat(["rgb(177,223,243)" 6],outer = 6))
    hline!([1-data["epsilon"]],lw=2,color="black",label="95% delivery")

    #=

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
        title = "CM delivery on activated ancillary bids", bar_width = 0.67,
        lw = 0, framestyle = :box,xtickfontsize=5,legend=:right)
=#
    end

function q_hour_mean(model_results)
    #model_results = in_S_results  
       # C,T,S
        
    # Time
    
    F = data["F"]
    E = data["E"]
    C = data["C"]
    T = data["T"]
    SC = model_results["SC"]
    S = model_results["S"]
    no_y = ones(24,length(S))
    pi_s = model_results["pi_s"]

    # Not delivered
    q_FFR_up    = [sum(pi_s[s]*model_results["q_FFR"][t,s]*no_y[t,s] for  s in S) for t in T]
    q_FCRD_up   = [sum(pi_s[s]*model_results["q_FCRD"][t,s]*no_y[t,s] for  s in S) for t in T]
    q_FCRN_up   = [sum(pi_s[s]*model_results["q_FCRN_up"][t,s]*no_y[t,s] for  s in S) for t in T]
    q_FCRN_down = [sum(pi_s[s]*model_results["q_FCRN_down"][t,s]*no_y[t,s] for  s in S) for t in T]
    q_mFRR_up   = [sum(pi_s[s]*model_results["q_mFRR"][t,s]*no_y[t,s] for  s in S) for t in T]
    
    q_sum = q_FFR_up + q_FCRD_up + q_FCRN_up + q_FCRN_down + q_mFRR_up
    
    # Delivered
    p_FFR_up   = [sum([pi_s[s]*sum(model_results["p_c_FFR_up"][C,t,s])*no_y[t,s] for s in S]) for  t in T]
    p_FCRD_up   = [sum([pi_s[s]*sum(model_results["p_c_FCRD_up"][C,t,s])*no_y[t,s] for s in S]) for  t in T]
    p_FCRN_up   = [sum([pi_s[s]*sum(model_results["p_c_FCRN_up"][C,t,s])*no_y[t,s] for s in S]) for  t in T]
    p_FCRN_down   = [sum([pi_s[s]*sum(model_results["p_c_FCRN_down"][C,t,s])*no_y[t,s] for s in S]) for  t in T]
    p_mFRR_up   = [sum([pi_s[s]*sum(model_results["p_c_mFRR_up"][C,t,s])*no_y[t,s] for s in S]) for  t in T]
    

    p_sum = p_FFR_up + p_FCRD_up + p_FCRN_up + p_FCRN_down + p_mFRR_up
    

    # Activated and Accepted Bids 
    act_FFR    = [sum(pi_s[s]*model_results["p_FFR"][t].*data["FFR_act"][t,SC[s]]*no_y[t,s] for c in C,s in S) for t in T]
    act_FCRD   = [sum(pi_s[s]*model_results["p_c_FCRD"][c,t].*data["FCRD_up_act"][t,SC[s]].*data["FCRD_accept"][c,t,SC[s]]*no_y[t,s] for  c in C,s in S) for t in T]
    act_FCRN_up   = [sum(pi_s[s]*model_results["p_c_FCRN"][c,t].*data["FCRN_up_act"][t,SC[s]].*data["FCRN_accept"][c,t,SC[s]]*no_y[t,s] for  c in C,s in S) for t in T]
    act_FCRN_down = [sum(pi_s[s]*model_results["p_c_FCRN"][c,t].*data["FCRN_down_act"][t,SC[s]].*data["FCRN_accept"][c,t,SC[s]]*no_y[t,s] for  c in C,s in S) for t in T]
    act_mFRR   = [sum(pi_s[s]*model_results["p_mFRR"][t].*data["mFRR_act"][t,SC[s]]*no_y[t,s] for  c in C,s in S) for t in T]
    act_sum = act_FFR + act_FCRD + act_FCRN_up + act_FCRN_down + act_mFRR
    

    #act_bids = [act_FFR, act_FCRD, act_FCRN_up, act_FCRN_down, act_mFRR, act_sum]
    
    #delivery = [p_FFR_up, p_FCRD_up, p_FCRN_up, p_FCRN_down, p_mFRR_up, p_sum]
    
    #act_unmet = [q_FFR_up, q_FCRD_up, q_FCRN_up, q_FCRN_down, q_mFRR_up, q_sum]
    
    groupedbar([transpose(p_FFR_up) transpose(q_FFR_up)],
    group=repeat(["2 Delivered" , "1 Not Delivered"], inner = 24),
    bar_position = :stack,
    title = "Mean FFR",    
    xlabel = "Hours",
    ylabel = "kW",
    #labels = ["Delivered" "Not Delivered"],
    legendfontsize = 7,
    titlefontsize = 9,
    xtickfontsize=7,
    ytickfontsize=7,
    guidefontsize=7,
    legend = :left,
    xticks = [1, 6, 12, 18, 24],
    linetype = :steppost,
    c = repeat([3, "rgb(198,227,171)"],inner=24), #, 5, "rgb(0,0,205)" "rgb(173,255,47)", 3, 4, 1, 2
    )
    FFR = Plots.plot!([i for i in 1:24], act_FFR,
        label = "Activation", 
        color = "black",
        linewidth = 2,
        linetype = :stepmid)

    groupedbar([transpose(p_FCRD_up) transpose(q_FCRD_up)],
    group=repeat(["2 Delivered" , "1 Not Delivered"], inner = 24),
    bar_position = :stack,
    title = "Mean FCRD",    
    xlabel = "Hours",
    ylabel = "kW",
    #labels = ["Delivered" "Not Delivered"],
    legendfontsize = 7,
    titlefontsize = 9,
    xtickfontsize=7,
    ytickfontsize=7,
    guidefontsize=7,
    legend = :topright,
    xticks = [1, 6, 12, 18, 24],
    linetype = :steppost,
    c = repeat([4, "rgb(232,163,178)"],inner=24),
    )
    FCRD = Plots.plot!([i for i in 1:24], act_FCRD,
    label = "Activation", 
    color = "black",
    linewidth = 2,
    linetype = :stepmid)
    
    groupedbar([transpose(p_FCRN_up) transpose(q_FCRN_up)],
    group=repeat(["2 Delivered","1 Not Delivered"], inner = 24),
    bar_position = :stack,
    title = "Mean FCRN up",    
    xlabel = "Hours",
    ylabel = "kW",
    #labels = ["Delivered" "Not Delivered"],
    legendfontsize = 7,
    titlefontsize = 9,
    xtickfontsize=7,
    ytickfontsize=7,
    guidefontsize=7,
    legend = :topright,
    xticks = [1, 6, 12, 18, 24],
    linetype = :steppost,
    c = repeat([1, "rgb(177,223,243)"],inner=24),
    )
    FCRN_up = Plots.plot!([i for i in 1:24], act_FCRN_up,
    label = "Activation", 
    color = "black",
    linewidth = 2,
    linetype = :stepmid)
    #
    groupedbar([transpose(p_FCRN_down) transpose(q_FCRN_down)],
    group=repeat(["2 Delivered" , "1 Not Delivered"], inner = 24),
    bar_position = :stack,
    title = "Mean FCRN down",    
    xlabel = "Hours",
    ylabel = "kW",
    #labels = ["Delivered" "Not Delivered"],
    legendfontsize = 7,
    titlefontsize = 9,
    xtickfontsize=7,
    ytickfontsize=7,
    guidefontsize=7,
    legend = :topright,
    xticks = [1, 6, 12, 18, 24],
    linetype = :steppost,
    c = repeat(["rgb(0,0,205)", "rgb(177,223,243)"],inner=24),
    )
    FCRN_down = Plots.plot!([i for i in 1:24], act_FCRN_down,
    label = "Activation", 
    color = "black",
    linewidth = 2,
    linetype = :stepmid)
    #
    groupedbar([transpose(p_mFRR_up) transpose(q_mFRR_up)],
    group=repeat(["2 Delivered" , "1 Not Delivered"], inner = 24),
    bar_position = :stack,
    title = "Mean mFRR up",    
    xlabel = "Hours",
    ylabel = "kW",
    #labels = ["Delivered" "Not Delivered"],
    legendfontsize = 7,
    titlefontsize = 9,
    xtickfontsize=7,
    ytickfontsize=7,
    guidefontsize=7,
    legend = :topright,
    xticks = [1, 6, 12, 18, 24],
    linetype = :steppost,
    c = repeat([2, "rgb(255,199,174)"],inner=24),
    )
    mFRR = Plots.plot!([i for i in 1:24], act_mFRR,
    label = "Activation", 
    color = "black",
    linewidth = 2,
    linetype = :stepmid)


    #groupedbar([transpose(p_sum) transpose(q_sum)],
    groupedbar([transpose(p_FFR_up) transpose(p_FCRD_up) transpose(p_FCRN_up) transpose(p_FCRN_down) transpose(p_mFRR_up) transpose(q_FFR_up) transpose(q_FCRD_up) transpose(q_FCRN_up) transpose(q_FCRN_down) transpose(q_mFRR_up)],
    group=repeat(["2 FFR Delivered" , "2 FCRD Delivered", "2 FCRN up Delivered", "2 FCRN down Delivered", "2 mFRR Delivered", "1 FFR Not Delivered", "1 FCRD Not Delivered", "1 FCRN up Not Delivered", "1 FCRN down Not Delivered", "1 mFRR Not Delivered"], inner = 24),
    bar_position = :stack,
    title = "All markets",    
    xlabel = "Hours",
    ylabel = "kW",
    legendfontsize = 4,
    titlefontsize = 9,
    xtickfontsize=7,
    ytickfontsize=7,
    guidefontsize=7,
    legend = :top,
    xticks = [1, 6, 12, 18, 24],
    linetype = :steppost,
    c = repeat([3, 4, 1, "rgb(0,0,205)", 2,#
        "rgb(198,227,171)", "rgb(232,163,178)", "rgb(177,223,243)", "rgb(177,223,243)", "rgb(255,199,174)"],inner=24)
        )
    all_AM = Plots.plot!([i for i in 1:24], act_sum,
    label = "Activation", 
    color = "black",
    linewidth = 2,
    linetype = :stepmid)
                    
    Plots.plot(FFR,FCRD,FCRN_up,FCRN_down,mFRR,all_AM, layout = 6)
    
    return FFR,FCRD,FCRN_up,FCRN_down,mFRR,all_AM
end