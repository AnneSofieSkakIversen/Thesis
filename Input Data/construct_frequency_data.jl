using CSV, JLD2, DataFrames,LinearAlgebra,PowerModels, Statistics, Dates,Plots

    datapath_freq = "C:\\Users\\Anne Sofie\\Desktop\\RawFrequency2021\\" 
        #G:\\Mit drev\\DTU Sustainable Energy\\F23\\Thesis\\Modelling\\Input Data\\FrequencyData2021\\"

    hist_freq = DataFrame(Freq = Float64[])
for m in 1:12  #"1" # Month of folder to process

    files_and_dirs = readdir(datapath_freq*"2021_$(m)") # readdir(pwd())  # reading files and directory
    hour_act_data = DataFrame(Hour = DateTime[], FFR_mean = Float64[], FCRD_mean = Float64[], FCRN_up_mean = Float64[],  FCRN_down_mean = Float64[]) 
   
    for i in 1:length(files_and_dirs)
            #if isfile(i)
                df = CSV.read(string(datapath_freq*"2021_$(m)\\",files_and_dirs[i]), dateformat = "yyyy-mm-dd hh:mm:ss", DataFrame;delim=",", decimal='.')
            #end
        rename!(df,[:string_time,:Frequency])
    
        myFormat = Dates.DateFormat("yyyy-mm-dd HH:MM:SS.s")
        df.string_time = Dates.DateTime.(df[!,:string_time], myFormat)
    
        Frequency = Array(df[!,:Frequency]) #[startIndex:endIndex]
        #=
        df.FFR = Frequency.*(Frequency .< 49.7)
        df.FCRD = (49.9 .- max.(49.5,Frequency)).*(Frequency .< 49.9) ./(49.9 - 49.5) # Linear up-activation % of FCRD bid from 0 to 1 between frequencies 49.9 and 49.5
        df.FCRN_up = (50 .- max.(49.9,Frequency)).*(Frequency .< 50) ./(50 - 49.9) # Linear up-activation % of FCRN bid from 0 to 1 between frequencies 50.0 and 49.9
        df.FCRN_down = (min.(50.1,Frequency) .-50).*(Frequency .> 50) ./(50.1 - 50) # Linear down-activation % of FCRN bid from 0 to 1 between frequencies 50.0 and 50.1

        df.Hour = floor.(df.string_time, Dates.Hour)  
        gd = groupby(transform(df, :string_time => x-> floor.(x,Dates.Hour)),:Hour)
        cd = combine(gd, [:FFR, :FCRD,:FCRN_up,:FCRN_down] .=> mean)
    
        append!(hour_act_data, cd)
        =#
        append!(hist_freq, DataFrame(Freq = Frequency[rand(1:length(Frequency),10000)]))
        
    end
    #CSV.write(datapath_freq*"freq_data_$(m).csv", hour_act_data)
end
CSV.write(datapath_freq*"hist_freq.csv", hist_freq)
histogram(hist_freq[!,1],bins = 60,fillalpha = 0,title = "Histogram: Frequencies 2021", label = "Raw data",ylabel = "Occurances",xlabel = "Frequency [Hz]")



    #N_samples = Int(floor(size(cd)[1]/24))
    #local FCRN_up_act1 = reshape(cd[1:24*N_samples,:FCRN_up_mean],(24,N_samples))
    #local FCRN_down_act1 = reshape(cd[1:(24*N_samples),:FCRN_down_mean],(24,N_samples))
    #local FCRD_act1 = reshape(cd[1:24*N_samples,:FCRD_mean],(24,N_samples))
    #local FFR_act1 = reshape(cd[1:24*N_samples,:FFR_mean],(24,N_samples))
    
#    return cd #FCRN_up_act1,FCRN_down_act1,FCRD_act1,FFR_act1
#end
