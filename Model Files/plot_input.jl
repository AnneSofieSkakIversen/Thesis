
using StatsPlots, Plots
using Statistics
using JuMP, JLD2, Random
# using DataFrames
# using LaTeXStrings
# using Plots.PlotMeasures

## To save image, use savefig("destination_path_and_file_name.png")
include("data_import_function_new.jl")

#       Choose scenario
#S = [1; 2 ; 5; 12; 19 ; 20]
# Set datapath
datapath = "G:\\Mit drev\\DTU Sustainable Energy\\F23\\Thesis\\Modelling\\Input Data\\Raw Data\\" 
data = data_import_function_new(datapath,false) # Put to true to set all market prices to zero (except spot)

F = data["F"] # The fully flexible prosumers
E = data["E"] # The prosumers with EV
S = data["S"] # Scenarios
C = union(F,E) # All prosumers

# Chance of activation
mean(data["FFR_act"])
mean(data["FCRN_up_act"])
mean(data["FCRN_down_act"])
mean(data["FCRD_up_act"])
mean(data["mFRR_act"])

mean(data["FFR_price"])
mean((data["FCRN_price"][1] .+ data["balance_price"].*(data["FCRN_up_act"] + data["FCRN_down_act"])).*data["FCRN_accept"][1,:,:]) # TCL
mean((data["FCRN_price"][10] .+ data["balance_price"].*(data["FCRN_up_act"] + data["FCRN_down_act"])).*data["FCRN_accept"][10,:,:]) # EV

mean((data["FCRD_price"][1] .+ data["balance_price"].*data["FCRD_up_act"]).*data["FCRD_accept"][1,:,:]) # TCL
mean((data["FCRD_price"][10] .+ data["balance_price"].*data["FCRD_up_act"]).*data["FCRD_accept"][10,:,:]) #EV

mean((data["mFRR_price"].+ data["balance_price"].*data["mFRR_act"]))




# Solar production profile
Plots.plot(Matrix(data["PV"][:,S])*data["PV_peak"][1]*data["n_demand"][1],
    label="",
    title = "Solar production profile - $(length(S)) scenarios",
    xlabel = "Hours",
    ylabel = "kW")

# Total Solar production profile
Plots.plot(sum(Matrix(data["PV"][:,S])*data["PV_peak"][c]*data["n_demand"][c] for c in C),
    label="",
    title = "Total Solar production profile - Scenario $(S)",
    xlabel = "Hours",
    ylabel = "kW")


# Community member 10 demand profile
Plots.plot(Matrix(data["D_household"][10,:,:])*data["n_demand"][10],
    label="",
    title = "Demand profile - prosumer 10",
    xlabel = "Hours",
    ylabel = "kW")
    
# Spot price
Plots.plot(Matrix(data["grid_in_price"][:,S]),
    label="",
    title = "Grid_in_price - Scenarios",
    xlabel = "Hours",
    ylabel = "DKK/kWh")

Plots.plot(Matrix(data["grid_out_price"][:,S]),
    label="",
    title = "Grid_out_price - Scenarios",
    xlabel = "Hours",
    ylabel = "DKK/kWh")

# Up regulation price
Plots.plot(Matrix(data["balance_up_price"]),
    label="",
    title = "Balance_up_price",
    xlabel = "Hours",
    ylabel = "DKK/kWh")

# Down regulation price
Plots.plot(Matrix(data["balance_down_price"]),
    label="",
    title = "Balance_down_price",
    xlabel = "Hours",
    ylabel = "DKK/kWh")

# FFR price
Plots.plot(Matrix(data["FFR_price"]),
    label="",
    title = "FFR_price",
    xlabel = "Hours",
    ylabel = "DKK/kWh")

# FFR activation profile
Plots.plot(Matrix(data["FFR_act"]),
    label="",
    title = "FFR_act",
    xlabel = "Hours",
    ylabel = "Activation")

    findmax(data["FFR_act"])
    data["FFR_act"][8,38]

# FCRN down activation
Plots.plot(Matrix(data["FCRN_down_act"]),
    label="",
    title = "FCRN_down_act",
    xlabel = "Hours",
    ylabel = "Activation")

# FCRN up activation
plot(Matrix(data["FCRN_up_act"]),
    label="",
    title = "FCRN_up_act",
    xlabel = "Hours",
    ylabel = "Activation")

# FCRN acceptance
Plots.plot(Matrix(data["FCRN_accept"][5,:,:]),
    label="",
    title = "FCRN_acceptance",
    xlabel = "Hours",
    ylabel = "Acceptance")

# FCRN price
plot(Matrix(data["FCRN_price"]),
    label="",
    title = "FCRN reserve price",
    xlabel = "Hours",
    ylabel = "DKK/kWh")

# FCRD acceptance
scatter(Matrix(data["FCRD_accept"]),
    label="",
    title = "FCRD_acceptance",
    xlabel = "Hours",
    ylabel = "Acceptance")

# FCRD price
plot(Matrix(data["FCRD_price"]),
    label="",
    title = "FCRD reserve price",
    xlabel = "Hours",
    ylabel = "DKK/kWh")

# mFRR up activation
scatter(Matrix(data["mFRR_act"]),
    label="",
    title = "mFRR_activation",
    xlabel = "Hours",
    ylabel = "Activation")

# mFRR price
plot(Matrix(data["mFRR_price"]),
    label="",
    title = "mFRR reserve price",
    xlabel = "Hours",
    ylabel = "DKK/kWh")

# Plot mean of all the market prices
grid_in_price_mean = mean(eachcol(data["grid_in_price"]))
grid_out_price_mean = mean(eachcol(data["grid_out_price"]))
balance_up_price_mean = mean(eachcol(data["balance_up_price"]))
balance_down_price_mean = mean(eachcol(data["balance_down_price"]))
FFR_price_mean = mean(eachcol(data["FFR_price"]))
FCRN_price_mean = mean(eachcol(data["FCRN_price"]))
FCRD_price_mean = mean(eachcol(data["FCRD_price"]))
mFRR_price_mean = mean(eachcol(data["mFRR_price"]))

plot([i for i in 1:24],
    [grid_in_price_mean grid_out_price_mean balance_up_price_mean balance_down_price_mean],
    label = ["grid_in_price" "grid_out_price" "balance_up_price" "balance_down_price"],
    xlabel = "Hour (h)",
    ylabel = "Price [DKK/kWh]",
    legend = :outerright,
    lw=2,
    legendfontsize = 10,
    ytickfontsize=10,
    xtickfontsize=10,
    guidefontsize = 12)

plot([i for i in 1:24],
    [FFR_price_mean FCRN_price_mean FCRD_price_mean mFRR_price_mean],
    label = ["FFR_price" "FCRN_price" "FCRD_price" "mFRR_price"],
    xlabel = "Hour (h)",
    ylabel = "Price [DKK/kWh]",
    legend = :outerright,
    lw=2,
    legendfontsize = 10,
    ytickfontsize=10,
    xtickfontsize=10,
    guidefontsize = 12)


# Frequency plot
    freq_folder = "C:\\Users\\Anne Sofie\\Desktop\\RawFrequency2021\\2021_2\\"
    freq_file = "2021-02-07.csv"

    using Dates
    start = 320980-599 #287759 #320980-33221 #320980
    slut = start + 599*2

    freq = CSV.read(string(freq_folder,freq_file), dateformat = "yyyy-mm-dd hh.HHH",DataFrame;delim=",", decimal='.')
    findmin(freq[:,2])
   
    #freq = freq[start:slut,:]

time = freq[:,1]
    
p = Plots.plot(freq[start:slut,2],
    label = "Frequency",
    xlabel = "Seconds",
    ylabel = "Frequency [Hz]",
    legend = :topright,
    lw=1,
    color = "black",
    ylims = (49.45,50.55),
    xticks = (0:100:2*600,string.(0:10:120)),
    legendfontsize = 8,
    ytickfontsize=10,
    xtickfontsize=10,
    guidefontsize = 12)
  Plots.hline!([50 49.9 50.1 49.5 50.5], 
               labels = "",#50 Hz" "49.9 Hz" "50.1 Hz" "49.5 Hz" "50.5 Hz"]
               color = ["black" 1 1 4 4] )
    Plots.hspan!(p,[49.9, 50], alpha = 0.2, 
                label = "FCRN up",
                color = 1)
    Plots.hspan!(p,[50,50.1], alpha = 0.4, 
                label = "FCRN down",
                color = 1)
    Plots.hspan!(p,[49.5, 49.9], alpha = 0.2, 
                label = "FCRD up",
                color = 4)
    Plots.hspan!(p,[50.1,50.5], alpha = 0.4, 
                label = "FCRN down",
                color = 4)
    Plots.hspan!(p,[49.0, 49.7], alpha = 0.2, 
                label = "FFR & mFRR up",
                color = 3)
    
    #png("Behavior\\freq_example.png")
hist_freq = CSV.read("G:\\Mit drev\\DTU Sustainable Energy\\F23\\Thesis\\Modelling\\Input Data\\hist_freq.csv", DataFrame; header=1, delim=",")
histogram(hist_freq[!,1],bins = 200,fillalpha = 0,
            title = "Histogram: Frequencies and Activation zones", 
            label = "Frequency data",
            ylabel = "Occurances",
            xlabel = "Frequency [Hz]",
            legendfontsize = 10,
            ytickfontsize=10,
            xtickfontsize=10,
            guidefontsize = 12,
            color = "black",
            legend = :topleft,
            linewidth = 2,
            xticks = ([49.5, 49.7, 49.9, 50, 50.1], ["49.5" "49.7" "49.9" "50.0" "50.1"]),
            xlims = (minimum(hist_freq[!,1]),maximum(hist_freq[!,1])))
Plots.vline!([50 49.9 50.1 49.5 50.5], 
           labels = "",#50 Hz" "49.9 Hz" "50.1 Hz" "49.5 Hz" "50.5 Hz"]
           color = ["black" 1 1 4 4] )
Plots.vspan!([49.9, 50], alpha = 0.2, 
            label = "FCRN up",
            color = 1)
Plots.vspan!([50,50.1], alpha = 0.4, 
            label = "FCRN down",
            color = 1)
Plots.vspan!([49.5, 49.9], alpha = 0.2, 
            label = "FCRD up",
            color = 4)
Plots.vspan!([50.1,50.5], alpha = 0.4, 
            label = "FCRN down",
            color = 4)
Plots.vspan!([49.0, 49.7], alpha = 0.2, 
            label = "FFR up",
            color = 3)
#png("Behavior\\freq_histogram.png")

f = transpose([49.4 49.5 49.9 50.1 50.2])
FFR = transpose([1 0 0 0 0])
FCRN = transpose([1 1 1 -1 -1])
FCRD = transpose([1 1 0 0 0])

Plots.plot(f,[FFR FCRN FCRD],
        title = "Activation levels from Frequency", 
        label = ["FFR" "FCRN" "FCRD"],
        ylabel = "Activation level [%]",
        xlabel = "Frequency [Hz]",
        legendfontsize = 10,
        ytickfontsize=10,
        xtickfontsize=10,
        color = [3 1 4],
        legend = :bottomleft,
        linewidth = 2,
        ylims = (-1.1,1.1),
        xticks = ([49.4, 49.5, 49.9, 50, 50.1], ["49.4" "49.5" "49.9" "50" "50.1"]),
        yticks = ([0, 0.5, 1], ["0%" "50%" "100%"]),
        seriestype = [:steppost :path :path ],
        guidefontsize = 12)
Plots.vline!([50], 
        label = "50 Hz",
        color = "black",
        linestyle = :dash,
        linewidth = 3)
annotate!(49.4, 0.1, text("Up regulation", :black, :left, 10, rotation = 90))
annotate!(50.2, -0.1, text("Down regulation", :black, :left, 10, rotation = 270))
Plots.plot!(twinx(),
yticks = ([-1, -0.5, 0], 
        ["100%" "50%" "0%"]),
        ylims = (-1.1,1.1),
        ytickfontsize=10)
#Plots.hspan!([-2, 0], alpha = 0.2, 
#            label = "",
#            color = "yellow")
#png("Behavior\\Freq_activation.png")


# EV mean connection
T = [i for i in 1:24]
C = [i for i in 1:10]
    Plots.plot([mean(sum(data["EV_con"][c,t,:] for c in C)) for t in T],
    label="",
    title = "Mean number of EV's connected",
    xlabel = "Hours",
    xticks = [0,6,12,18, 24],
    lw = 2,
    legendfontsize = 10,
    ytickfontsize=10,
    xtickfontsize=10,
    guidefontsize = 12,
    ylabel = "Number of EV's")
    png("Behavior\\EV_con_example.png")


S = 260
E = [i for i in 5:10]
    EV_con = zeros(length(C),length(T),S)
    EV_con[E,T,S] .= 1
    dist_disc = Poisson(3)
    dist_rec = Poisson(4)
 
bar([5 .+ Poisson(3),13 .+ Poisson(4)],
    xlims = (0,24.5),
    ylabel = "Probability",
    title = "Distributions of EV departures and arrivals",
    xlabel = "Hours",
    linestyle = :solid,
    linewidth = 1,
    alpha = 0.7,
    legendfontsize = 10,
    ytickfontsize=10,
    xtickfontsize=10,
    guidefontsize = 12,
    legend = :right,
    yticks = ([0, 0.05, 0.1,0.15,0.20], ["0%" "5%" "10%" "15%" "20%"]),
    label= ["Departure time" "Arrival time"],
    xticks = [0,6,12,18, 24],)
    png("Behavior\\EV_Poissons.png")



### Baseline for ambient temperature
TCL_results = load("OS_results\\model_IS_260_no_risk_results.jld2")
data = load("const_data.jld2")
begin
s = 34
h = [i for i in 1:24]
    SC = TCL_results["SC"]
    T_a = data["T_a"][h,SC[s]]
    T_opt = ones(length(h)).*data["T_opt"]
    T_TCL = TCL_results["T_TCL"][2,h,s]
    D_base_TCL = data["D_base_TCL"][2,h,SC[s]]
    d_TCL = TCL_results["d_TCL"][2,h,s]
    Plots.plot(D_base_TCL, c = 1,fillrange = 9,alpha = 0.3,label = "Down-regulation Capacity")
    Plots.plot!(D_base_TCL, c = 2, fillrange = 0,alpha = 0.3,label = "Up-regulation Capacity")
    Plots.plot!(max.(D_base_TCL,d_TCL), c = 1,fillrange = D_base_TCL,alpha = 0.8,label = "Delivered Down-regulation")
    Plots.plot!(min.(D_base_TCL,d_TCL), c = 2,fillrange = D_base_TCL,alpha = 0.8,label = "Delivered Up-regulation")
    Plots.plot!([i for i in h], [D_base_TCL,d_TCL],
    title = "TCL avg. Consumption and Temperature",
    legs = (:topleft),
    xlabel = "Hours",
    ylabel = "Power [kW]",
    color = "black",
    label = ["TCL baseline consumption" "TCL Consumption"],
    linewidth = [1 2],
    linestyle = [:dash :solid],
    xticks = [1,6, 12, 18, 24],
    xlims = (1,24),
    legend = :topright,
    ylims = (-0.01,2),
    legendfontsize = 10,
    ytickfontsize=10,
    xtickfontsize=10,
    guidefontsize = 12,
    #right_margin = 2mm
    )
end
png("Behavior\\TCL_capacity.png")

maximum(data["FFR_price"])
   