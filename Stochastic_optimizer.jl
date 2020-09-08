using JuMP, Ipopt
using KNITRO
using DelimitedFiles

method = "Reg" # CS, Binary, Reg
RegGap = 0.05
hours =  72
# hours =  1
Scen = "Solar"
temp_Eprice = readdlm("Prices_$hours" * "hr_$Scen" * ".csv", ',',header=false)
Prob = readdlm("Probs_$hours" * "hr_$Scen" * ".csv", ',', header=false)

temp_Eprice = temp_Eprice

Nscen = size(temp_Eprice, 2)

Price_data = temp_Eprice
m = Model(with_optimizer(KNITRO.Optimizer, maxit = 100_000, ms_enable=1, ms_maxsolves=100, par_blasnumthreads=8))

timeset = 1:hours
scenset = 1:Nscen
@variable(m, 0 <= j_ch[t in timeset, s in scenset] <= 1500) # current density of charging
@variable(m, 0 <= j_dis[t in timeset, s in scenset] <= 200) #
# @variable(m, 0 <= j_ch[t in timeset, s in scenset] <= 3000) # 2x current density
# @variable(m, 0 <= j_dis[t in timeset, s in scenset] <= 400) # 2x current density

@variable(m, 0.1 <= a_ch <= 20) # area of charge cell
@variable(m, 0.1 <= a_dis <= 20) # area of discharge cell
@variable(m, 10 <= v <= 1e4) # volume of tank [liters]
# When controlling area and volume ratios, use below
# @constraint(m, ach_adis, a_ch == 1.51*a_dis )
# @constraint(m, v_adis, v == a_dis*0.27*1000 )


@variable(m, c[t in timeset, s in scenset] >= 0)
@variable(m, d[t in timeset, s in scenset] >= 0)
@variable(m, soc[t in timeset, s in scenset] >= 0)


@variable(m, eta_FE[t in timeset, s in scenset] >= 0)

@variable(m, CostElec[t in timeset, s in scenset] >= maximum(Price_data/1e6)*((-0.0005*(178)^2 + 0.1796*(178)) * 20)*(-1))
@variable(m, CostCap >= 0)
@variable(m, CostStor >= 0)
@variable(m, CostH >= 0)

Lambda_E = Price_data / 1e6
Lambda_P = 2000
Lambda_V = 0.56
F = 96485 #Faraday's constant
K = 3600/(F*1000)
r = 0.05
h = 10 # lifetime
# Gamma = (r*(1+r)^h)/(((1+r)^h - 1)*365) # h=10years / 5%
Gamma = 3.5481e-4 * length(timeset) / 24

Cbar = 1500 * (4e-7(1500)^2 + 0.0006*1500 +1.58) * 20
# Cbar = 3000 * (4e-7(3000)^2 + 0.0006*3000 +1.58) * 20
Dbar = (-0.0005 * (178)^2 + 0.1796*178) * 20
# 0.0001 *j_dis is for double current density
# Dbar = (-0.0001 * (356)^2 + 0.0898*356) * 20

@constraint(m, iniSOC[s in scenset], soc[1,s] == F * 0.6 * v /2 )
@constraint(m, finalSOC[s in scenset], soc[end,s] == F * 0.6 * v /2 )
# added constraint below (lowSOC) in order to keep SOC below 0.95
@constraint(m, chargeUB[t in timeset, s in scenset], c[t,s] <= Cbar)
@constraint(m, dischargeUB[t in timeset, s in scenset], d[t,s] <= Dbar)

if method == "CS"
    @variable(m, CSL[z in 1:length(timeset)*length(scenset)] >= 0)
    @variable(m, CSR[z in 1:length(timeset)*length(scenset)] >= 0)
    @constraint(m, CSL_def[t in timeset, s in scenset], CSL[Nscen*(s-1)+t] == j_ch[t,s]) #113.650379...
    @constraint(m, CSR_def[t in timeset, s in scenset], CSR[Nscen*(s-1)+t] == j_dis[t,s])
    @constraint(m, CS, [CSL; CSR] in MOI.Complements(length(timeset)*length(scenset)))
elseif method == "Binary"
    @variable(m, u[t in timeset, s in scenset], Bin)
    @constraint(m, CS1[t in timeset, s in scenset], j_ch[t,s] <= Cbar*u[t,s])
    @constraint(m, CS2[t in timeset, s in scenset], j_dis[t,s] <= Dbar*(1-u[t,s]))
elseif method == "Reg"
    @NLconstraint(m, CS1[t in timeset, s in scenset], (j_ch[t,s] + RegGap) * (j_dis[t,s] + RegGap) - RegGap^2 >= 0)
    @NLconstraint(m, CS2[t in timeset, s in scenset], (j_ch[t,s]) * (j_dis[t,s]) - RegGap^2 <= 0)
end

@constraint(m, tank[t in timeset, s in scenset], soc[t,s] <= F * 0.6 * v * 0.95)

# This is the normal one
@NLconstraint(m, chDef[t in timeset, s in scenset], c[t,s] == j_ch[t,s] * (4e-7*(j_ch[t,s])^2 + 6e-4*j_ch[t,s] + 1.58) * a_ch)
# This is for half overpotential (divide current density in voltage def by 2 for double current density)
# @NLconstraint(m, chDef[t in timeset, s in scenset], c[t,s] == j_ch[t,s] * ((4e-7*(j_ch[t,s])^2 + 6e-4*(j_ch[t,s]) + 1.58)-((4e-7*(j_ch[t,s])^2 + 6e-4*(j_ch[t,s]) + 1.58) - 1.6)/2) * a_ch)
# This is the normal one
@NLconstraint(m, disDef[t in timeset, s in scenset], d[t,s] == (-0.0005*(j_dis[t,s])^2 + 0.1796*j_dis[t,s]) * a_dis)
# This is for half overptential and double current density
# @NLconstraint(m, disDef[t in timeset, s in scenset], d[t, s] == j_dis[t,s]*(0.37 - (0.37 - (-0.0001*(j_dis[t,s])^2 + 0.0898*j_dis[t,s])/(j_dis[t,s]))/2) * a_dis)
# This is for half overpotential
# @NLconstraint(m, disDef[t in timeset, s in scenset], d[t, s] == j_dis[t,s]*(0.37 - (0.37 - (-0.0005*(j_dis[t,s])^2 + 0.1796*j_dis[t, s])/(j_dis[t,s]))/2) * a_dis)
@NLconstraint(m, SOC[t in setdiff(timeset,[1]), s in scenset], soc[t,s] == soc[t-1,s] + 3600*j_ch[t,s]*a_ch*eta_FE[t,s] - 3600*j_dis[t,s]*a_dis)
@NLconstraint(m, etaDef[t in setdiff(timeset,[1]), s in scenset], eta_FE[t,s] == -0.6857*((soc[t-1,s]/(F*0.6*v))^2) - 0.3143 * (soc[t-1,s]/(F*0.6*v)) + 0.9343)
@NLconstraint(m, iniFE[s in scenset], eta_FE[1,s] == -0.6857*((soc[1,s]/(F*0.6*v))^2) - 0.3143 * (soc[1,s]/(F*0.6*v)) + 0.9343)
# This is for improved FE
# @NLconstraint(m, etaDef[t in setdiff(timeset,[1]), s in scenset], eta_FE[t, s] == 1 - (1 - (-0.6857*((soc[t-1,s]/(F*0.6*v))^2) - 0.3143 * (soc[t-1,s]/(F*0.6*v)) + 0.9343))/2)
# @NLconstraint(m, iniFE[s in scenset], eta_FE[1, s] == 1 - (1 - (-0.6857*((soc[1,s]/(F*0.6*v))^2) - 0.3143 * (soc[1,s]/(F*0.6*v)) + 0.9343))/2)

@NLconstraint(m, CostElec_Def[t in timeset, s in scenset], CostElec[t,s] == (Lambda_E[t,s] * (c[t,s] - d[t,s])))
@NLconstraint(m, CostCap == (Gamma * Lambda_P *(a_ch + a_dis)))
@NLconstraint(m, CostStor == (Gamma * Lambda_V * v))

# @NLconstraint(m, objDef, CostH == (sum(Prob[s] * CostElec[t,s] / sum(Prob[ss] for ss in scenset) for t in setdiff(timeset, [1]), s in scenset) + CostCap + CostStor) /
#     (K * sum(Prob[s] * j_ch[t,s] * a_ch  / sum(Prob[ss] for ss in scenset) for t in setdiff(timeset,[1]), s in scenset)))
@NLconstraint(m, objDef, CostH == (sum(Prob[s] * CostElec[t,s] for t in setdiff(timeset, [1]), s in scenset) + CostCap + CostStor) /
    (K * sum(Prob[s] * j_ch[t,s] * a_ch for t in setdiff(timeset,[1]), s in scenset)))
denom = K * sum(Prob[s] * j_ch[t,s] * a_ch for t in setdiff(timeset,[1]), s in scenset)
@NLobjective(m, Min, CostH)
optimize!(m)


@show objective_value(m)
@show value.(CostElec)
@show value.(CostCap)
@show value.(CostStor)
@show value.(a_ch)
@show value.(a_dis)
@show value.(v)
@show value.(soc)/(F*0.6*value.(v))


SOC_result = value.(soc)/(F*0.6*value.(v))
# writedlm("SOC.csv", SOC_result, ',')

Preamb = "ach_adis_up"
writedlm("$Preamb" * "Second_Cost_$hours" * "hr_$Scen" * ".csv", objective_value(m), ',')
writedlm("$Preamb" * "Second_a_ch_$hours" * "hr_$Scen" * ".csv", value.(a_ch), ',')
writedlm("$Preamb" * "Second_a_dis_$hours" * "hr_$Scen" * ".csv", value.(a_dis), ',')
writedlm("$Preamb" * "Second_v_$hours" * "hr_$Scen" * ".csv", value.(v), ',')
writedlm("$Preamb" * "Second_j_ch_$hours" * "hr_$Scen" * ".csv", value.(j_ch), ',')
writedlm("$Preamb" * "Second_j_dis_$hours" * "hr_$Scen" * ".csv", value.(j_dis), ',')
writedlm("$Preamb" * "Second_FE_$hours" * "hr_$Scen" * ".csv", value.(eta_FE), ',')
writedlm("$Preamb" * "Second_s_$hours" * "hr_$Scen" * ".csv", SOC_result, ',')
writedlm("$Preamb" * "Second_CostElec_$hours" * "hr_$Scen" * ".csv", value.(CostElec), ',')
writedlm("$Preamb" * "Second_CostCap_$hours" * "hr_$Scen" * ".csv", value.(CostCap), ',')
writedlm("$Preamb" * "Second_CostStor_$hours" * "hr_$Scen" * ".csv", value.(CostStor), ',')
writedlm("$Preamb" * "Second_denom_$hours" * "hr_$Scen" * ".csv", value.(denom), ',')


print("status: ",termination_status(m))
