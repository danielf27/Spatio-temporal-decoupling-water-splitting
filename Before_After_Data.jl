using DelimitedFiles, Clustering
using Distances
#Before and after thing
filename_nyiso = "caiso_data.csv"

# Eprice_year = readdlm(filename, ',',header=true)[1]
temp_Eprice = readdlm(filename_nyiso, ',',header=true)[1]
# For nyiso_data: 1=Bal, 2=Low, 3=Solar, 4=Wind
Eprice_year = temp_Eprice[:,4]

BA_time = 24 #hours
tot_time = 24+BA_time*2
Scen = "Solar"
Tot_time = BA_time*2+24
daysoff = Int(floor((BA_time-1)/24)+1)
Nperiod = Int(365 - (daysoff*2))

NewPrices = zeros(Tot_time, Nperiod)
k = 1
for j in (daysoff-1)*24+2:length(Eprice_year)-24*daysoff
    if j % 24 == 1
        global NewPrices[:,k] = Eprice_year[j-BA_time:j+23+BA_time]
        global k = k+1
    else
        continue
    end
end

clust_size = 10

# cluster X into 20 clusters using K-means
R = kmeans(NewPrices, clust_size; maxiter=200, display=:iter)

@assert nclusters(R) == clust_size # verify the number of clusters

a = assignments(R) # get the assignments of points to clusters
c = counts(R) # get the cluster sizes
M = R.centers # get the cluster centers
Prob = c/Nperiod

writedlm("Caiso_Prices_$tot_time" * "hr_$Scen" * ".csv", M, ',')
writedlm("Caiso_Probs_$tot_time" * "hr_$Scen" * ".csv", Prob, ',')
