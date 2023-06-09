using StatsPlots, Plots
using Statistics
# using DataFrames
# using LaTeXStrings
# using Plots.PlotMeasures

## To save image, use savefig("destination_path_and_file_name.png")
include("data_import_function.jl")
using CSV, JLD, DataFrames
datapath = "G:\\Mit drev\\DTU Sustainable Energy\\E22\\XXXXX Special Course\\Paper - Scenario Matrix\\Processed data\\"
matrix_construct(datapath)
data = data_import(datapath)

#       Choose scenario
S = data["S"][100:183]

N = data["N"] # The inflexible prosumers
F = data["F"] # The fully flexible prosumers
E = data["E"] # The prosumers with EV
S = data["S"] # Scenarios
C = union(N,F,E) # All prosumers
IN = union(N,E) # All prosumers that are inflexible


# Solar production profile
plot(Matrix(data["PV"][:,S])*data["PV_peak"][1]*data["n_demand"][1],
    label="",
    title = "Solar production profile - 20 scenarios",
    xlabel = "Hours",
    ylabel = "kW")

# Total Solar production profile
#plot(sum(Matrix(data["PV"][:,S])*data["PV_peak"][c]*data["n_demand"][c] for c in C),
#    label="",
#    title = "Total Solar production profile - Scenario $(S)",
#    xlabel = "Hours",
#    ylabel = "kW")


# Community member 10 demand profile
plot(Matrix(data["D_base"][10,:,:])*data["n_demand"][10],
    label="",
    title = "Demand profile - prosumer 10",
    xlabel = "Hours",
    ylabel = "kW")
    
# Spot price
plot(Matrix(data["grid_in_price"][:,S]),
    label="",
    title = "Grid_in_price - Scenario $(S)",
    xlabel = "Hours",
    ylabel = "DKK/kWh")

plot(Matrix(data["grid_out_price"][:,S]),
    label="",
    title = "Grid_out_price - Scenario $(S)",
    xlabel = "Hours",
    ylabel = "DKK/kWh")

# Up regulation price
plot(Matrix(data["balance_up_price"]),
    label="",
    title = "Balance_up_price",
    xlabel = "Hours",
    ylabel = "DKK/kWh")

# Down regulation price
plot(Matrix(data["balance_down_price"]),
    label="",
    title = "Balance_down_price",
    xlabel = "Hours",
    ylabel = "DKK/kWh")

# FFR price
plot(Matrix(data["FFR_price"]),
    label="",
    title = "FFR_price",
    xlabel = "Hours",
    ylabel = "DKK/kWh")

# FFR activation profile
plot(Matrix(data["FFR_act"]),
    label="",
    title = "FFR_act",
    xlabel = "Hours",
    ylabel = "Activation")

# FCRN down activation
plot(Matrix(data["FCRN_down_act"]),
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
scatter(Matrix(data["FCRN_accept"]),
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
