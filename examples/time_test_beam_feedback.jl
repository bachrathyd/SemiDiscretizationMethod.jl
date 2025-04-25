######Product-Free full mapping
######The mapping is repsented by a left and a right matrix
## Time-complexity measurements for a publication

5 + 5
###import Pkg
########
###Pkg.activate("")
######### Pkg.activate()# ez az eredeti, mindent tartalamzo
#########] dev SemiDiscretizationMethod

using Revise
using SemiDiscretizationMethod
using StaticArrays
using Plots

using CSV
using DataFrames
using BenchmarkTools

using LaTeXStrings
Threads.nthreads()


using SemiDiscretizationMethod
using StaticArrays

#include("beam_delay_feedback.jl")
## --------------- time test -------------------------
## -----------------------------------------------------
#Journal
BenchmarkTools.DEFAULT_PARAMETERS.samples = 1000.0
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 2.0


#Nv = ceil.(10 .^ (1.2:0.05:5.61)) #5.61  #100:100:3000
Nv = ceil.(10 .^ (1.2:0.05:6))
Nv = ceil.(10 .^ (1.2:0.1:5.5))
Nv = ceil.(10 .^ (1.2:0.05:4.5))
#Nv = ceil.(10 .^ (1.2:0.2:6))
Twaitfor_SH = 100.0;

#Fast(er) test
BenchmarkTools.DEFAULT_PARAMETERS.samples = 20.0
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1.0
#


#Nv = ceil.(10 .^ (1.0:0.025:3.00))
#Twaitfor_SH = 1.0;

kpow = 1e-4
NEV = 1
#kpow=1e-20  #has stron effect in T<tau
#NEV=20

tmake_SH_PRi = zeros(Float64, length(Nv));
tmake_SH_Ci = zeros(Float64, length(Nv));
tmake_SH_PhiALL = zeros(Float64, length(Nv));
teig_SH = zeros(Float64, length(Nv));
tfixP_SH = zeros(Float64, length(Nv));
μSH = zeros(Float64, length(Nv));
fixSH = zeros(Float64, length(Nv));

#deviation
tmake_SH_PRi_S = zeros(Float64, length(Nv));
tmake_SH_Ci_S = zeros(Float64, length(Nv));
tmake_SH_PhiALL_S = zeros(Float64, length(Nv));
teig_SH_S = zeros(Float64, length(Nv));
tfixP_SH_S = zeros(Float64, length(Nv));
fixSH_S = zeros(Float64, length(Nv));


domoreSH = true;
domoreSH_PR = true;
#@show domoreSH=false;

tmake_LR = zeros(Float64, length(Nv));
teig_LR = zeros(Float64, length(Nv));
tfixP_LR = zeros(Float64, length(Nv));
μLR = zeros(Float64, length(Nv));
fixLR = zeros(Float64, length(Nv));
#deviation
tmake_LR_S = zeros(Float64, length(Nv));
teig_LR_S = zeros(Float64, length(Nv));
tfixP_LR_S = zeros(Float64, length(Nv));
fixLR_S = zeros(Float64, length(Nv));


#f1 = (x) -> log(x) ./ log(10)
#f2 = (x) -> log(x) ./ log(10)

@time begin
    @warn "Beam modell is used!!! see the beam_delay_deedbacl.jl file"
    Tper = TperpTspeed * Tspeed
    T=Tper
    prob_lddep = delay_beam_feedback(E, A, ρ, L, η, N, P, τpTspeed, Tspeed, TperpTspeed, Val(N - 1))
    τmax = τpTspeed * Tspeed
    method = SemiDiscretization(2, Tper / Ndisc) # 3rd order semi discretization with Δt=0.1
    Nsteps = Int((Tper + 100eps(Tper)) ÷ method.Δt)
end

#kNdisc=10
#@profview  
for kNdisc in vcat([1, 1], 1:length(Nv)) #the first is repated to get read of the first compliation time
    println("=====================================================")
    if kNdisc == 1
        tmake_SH_PRi[kNdisc] = 0.0
        tmake_SH_Ci[kNdisc] = 0.0
        tmake_SH_PhiALL[kNdisc] = 0.0
        teig_SH[kNdisc] = 0.0
        tfixP_SH[kNdisc] = 0.0
        μSH[kNdisc] = 0.0
        fixSH[kNdisc] = 0.0

        domoreSH = true
        domoreSH_PR = true
        tmake_LR[kNdisc] = 0.0
        teig_LR[kNdisc] = 0.0
        tfixP_LR[kNdisc] = 0.0
        μLR[kNdisc] = 0.0
        fixLR[kNdisc] = 0.0
    end
    Ndisc = Nv[kNdisc]
    println([kNdisc / length(Nv), Ndisc])

    #@show Naver = maximum([ceil(4 - (log(Ndisc) / log(10))) * 2, 1])
    Naver = 1
    method = SemiDiscretization(1, T / Ndisc) # 3rd order semi discretization with Δt=0.1
    Nsteps = Int((T + 100eps(T)) ÷ method.Δt)





    for _ = 1:Naver


        t = @benchmark SemiDiscretizationMethod.calculateResults($prob_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true)
        print("PR calc: ")
        @show t
        tmake_SH_PRi[kNdisc] += BenchmarkTools.median(t).time / 1e9
        tmake_SH_PRi_S[kNdisc] += BenchmarkTools.std(t).time / 1e9
        resultPR = SemiDiscretizationMethod.calculateResults(prob_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true)

        if (kNdisc < 5 || tmake_SH_Ci[kNdisc-1] < Twaitfor_SH * Naver) #domoreSH_PR
            #tmake_SH_Ci[kNdisc] += @elapsed mmpp2 = SemiDiscretizationMethod.DiscreteMappingSteps(resultPR)
            
            print("Creating Ci:  ")
            t = @benchmark SemiDiscretizationMethod.DiscreteMappingSteps($resultPR)
            @show t
            tmake_SH_Ci[kNdisc] += BenchmarkTools.median(t).time / 1e9
            tmake_SH_Ci_S[kNdisc] += BenchmarkTools.std(t).time / 1e9
            mmpp2 = SemiDiscretizationMethod.DiscreteMappingSteps(resultPR)


            #domoreSH_PR = tmake_SH_Ci[kNdisc] < Twaitfor_SH * Naver
        else
            tmake_SH_Ci[kNdisc] = NaN
        end
        if (kNdisc < 5 || tmake_SH_PhiALL[kNdisc-1] < Twaitfor_SH * Naver)#domoreSH
            #println("domoreSH")
            #tmake_SH_PhiALL[kNdisc] += @elapsed mappingFull = DiscreteMapping_1step(mathieu_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true) #The discrete mapping of the system
            t = @benchmark DiscreteMapping_1step($prob_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true) #The discrete mapping of the system

            print("DiscreteMapping_1step - products:  ")
            @show t
            tmake_SH_PhiALL[kNdisc] += BenchmarkTools.median(t).time / 1e9
            tmake_SH_PhiALL_S[kNdisc] += BenchmarkTools.std(t).time / 1e9
            mappingFull = DiscreteMapping_1step(prob_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true) #The discrete mapping of the system

            print("spectralRadiusOfMapping full Phi:  ")
            #teig_SH[kNdisc] += @elapsed μSH[kNdisc] = spectralRadiusOfMapping(mappingFull,nev=NEV,tol=10.0^kpow)
            t = @benchmark spectralRadiusOfMapping($mappingFull, nev=NEV, tol=10.0^kpow)
            @show t
            teig_SH[kNdisc] += BenchmarkTools.median(t).time / 1e9
            teig_SH_S[kNdisc] += BenchmarkTools.std(t).time / 1e9
            μSH[kNdisc] = spectralRadiusOfMapping(mappingFull, nev=NEV, tol=10.0^kpow)

            # fixSH[kNdisc]  += @elapsed μSH[kNdisc] = fixPointOfMapping(mappingFull)
            t = @benchmark fixPointOfMapping($mappingFull)
            
            print("fixPointOfMapping:  ")
            @show t
            fixSH[kNdisc] += BenchmarkTools.median(t).time / 1e9
            fixSH_S[kNdisc] += BenchmarkTools.std(t).time / 1e9


            #domoreSH = tmake_SH_PhiALL[kNdisc] < Twaitfor_SH * Naver
        else
            tmake_SH_PhiALL[kNdisc] = NaN
            teig_SH[kNdisc] = NaN
            fixSH[kNdisc] = NaN
        end

        t = @benchmark DiscreteMapping_LR($prob_lddep, method, τmax, n_steps=Int((T + 100eps(T)) ÷ method.Δt), calculate_additive=true)#The discrete mapping of the system
        println("---- LR ---")
        print("DiscreteMapping_LR:  ")
        @show t
        tmake_LR[kNdisc] += BenchmarkTools.median(t).time / 1e9
        tmake_LR_S[kNdisc] += BenchmarkTools.std(t).time / 1e9
        mappingLR = DiscreteMapping_LR(prob_lddep, method, τmax, n_steps=Int((T + 100eps(T)) ÷ method.Δt), calculate_additive=true)#The discrete mapping of the system



        #teig_LR[kNdisc] += @elapsed μLR[kNdisc] = spectralRadiusOfMapping(mappingLR,nev=NEV,tol=10.0^kpow)
        t = @benchmark spectralRadiusOfMapping($mappingLR, nev=NEV, tol=10.0^kpow)
        
        print("DiscreteMapping_LR - spectralRadiusOfMapping  ")
        @show t
        teig_LR[kNdisc] += BenchmarkTools.median(t).time / 1e9
        teig_LR_S[kNdisc] += BenchmarkTools.std(t).time / 1e9
        μLR[kNdisc] = spectralRadiusOfMapping(mappingLR, nev=NEV, tol=10.0^kpow)

        # tfixP_LR[kNdisc] += @elapsed xPLR = fixPointOfMapping(mappingLR)
        t = @benchmark fixPointOfMapping($mappingLR)
        
        print("DiscreteMapping_LR - fixPointOfMapping  ")
        @show t
        fixLR[kNdisc] += BenchmarkTools.median(t).time / 1e9
        fixLR_S[kNdisc] += BenchmarkTools.std(t).time / 1e9

    end

    tmake_SH_PRi[kNdisc] /= Naver
    tmake_SH_Ci[kNdisc] /= Naver
    tmake_SH_PhiALL[kNdisc] /= Naver
    teig_SH[kNdisc] /= Naver
    fixSH[kNdisc] /= Naver

    tmake_LR[kNdisc] /= Naver
    teig_LR[kNdisc] /= Naver
    fixLR[kNdisc] /= Naver


    tmake_SH_PRi_S[kNdisc] /= Naver
    tmake_SH_Ci_S[kNdisc] /= Naver
    tmake_SH_PhiALL_S[kNdisc] /= Naver
    teig_SH_S[kNdisc] /= Naver
    fixSH_S[kNdisc] /= Naver

    tmake_LR_S[kNdisc] /= Naver
    teig_LR_S[kNdisc] /= Naver
    fixLR_S[kNdisc] /= Naver


    #if T=τmax
    #@show norm(xP-xPLR)
    #plot(xP[1:2:end])
    #plot!(xPLR[1:2:end])

    # @show norm(μSH[kNdisc]-μLR[kNdisc])

    #if (kNdisc>5000  ||  mod(kNdisc,10)==0)

    f1 = (x) -> log(x) ./ log(10)
    f2 = (x) -> log(x) ./ log(10)
    df = DataFrame(
        Nv_data=Nv,
        tmake_SH_data=tmake_SH_PhiALL,
        teig_SH_data=teig_SH,
        tfixP_SH_data=tfixP_SH,
        tmake_LR_data=tmake_LR,
        teig_LR_data=teig_LR,
        tfixP_LR_data=tfixP_LR,
        Nv_data_log=f2.(Nv),
        tmake_SH_data_log=f1.(tmake_SH_PhiALL),
        teig_SH_data_log=f1.(teig_SH),
        tfixP_SH_data_log=f1.(tfixP_SH),
        tmake_LR_data_log=f1.(tmake_LR),
        teig_LR_data_log=f1.(teig_LR),
        tfixP_LR_data_log=f1.(tfixP_LR))

    CSV.write("CPU_Timev__T_beam_delay_feedback" * string(T) * "_tau_" * string(τmax) * ".csv", df)
    println("data saved")
    #end




    tmake_SH_PRi_S[isnan.(tmake_SH_PRi_S)] .= 0.0
    tmake_SH_Ci_S[isnan.(tmake_SH_Ci_S)] .= 0.0
    tmake_SH_PhiALL_S[isnan.(tmake_SH_PhiALL_S)] .= 0.0
    teig_SH_S[isnan.(teig_SH_S)] .= 0.0
    fixSH_S[isnan.(fixSH_S)] .= 0.0

    tmake_LR_S[isnan.(tmake_LR_S)] .= 0.0
    teig_LR_S[isnan.(teig_LR_S)] .= 0.0
    fixLR_S[isnan.(fixLR_S)] .= 0.0


    #f1 = (x) -> log(x) ./ log(10)
    #f2 = (x) -> log(x) ./ log(10)

    f1 = (x) -> x == 0.0 ? NaN : x
    f2 = (x) -> x == 0.0 ? NaN : x

    ScaleSTD = 0.5

    #fplot=scatter
    #fplot! = scatter!
    fplot = plot
    fplot! = plot!

    #fplot(f2.(Nv), f1.(tmake_SH_PRi .+tmake_SH_PRi_S), labels="tmake_SH_PRi")
    #fplot!(f2.(Nv), f1.(tmake_SH_PRi ), labels="tmake_SH_PRi")
    #fplot!(f2.(Nv), f1.(tmake_SH_PRi .-tmake_SH_PRi_S), labels="tmake_SH_PRi")
    fplot(f2.(Nv), f1.(tmake_SH_PRi), ribbon=tmake_SH_PRi_S .* ScaleSTD, labels="tmake_SH_PRi")
    fplot!(f2.(Nv), f1.(tmake_SH_Ci), ribbon=tmake_SH_Ci_S .* ScaleSTD, labels="tmake_SH_Ci")
    #fplot!(f2.(Nv), f1.(tmake_SH_Ci),ribbon=tmake_SH_Ci_S .* ScaleSTD, labels="tmake_SH_Ci")
    fplot!(f2.(Nv), f1.(tmake_SH_PhiALL), ribbon=tmake_SH_PhiALL_S .* ScaleSTD, labels="tmake_SH_Phi_1step...prod(...)")
    fplot!(f2.(Nv), f1.(teig_SH), ribbon=teig_SH_S .* ScaleSTD, labels="teig_SH: eigs(PHI) only")
    #fplot!(f2.(Nv), f1.(tfixP_SH),ribbon=tfixP_SH_S .* ScaleSTD, labels="teig_SH: fix Point")
    fplot!(f2.(Nv), f1.(tmake_LR), ribbon=tmake_LR_S .* ScaleSTD, labels="tmake_LR", linewidth=4)
    fplot!(f2.(Nv), f1.(teig_LR), ribbon=teig_LR_S .* ScaleSTD, labels="teig_LR: eigs(ΦR,ΦL) only", linewidth=4)

    fplot!(f2.(Nv), f1.(fixSH), ribbon=fixSH_S .* ScaleSTD, labels="fixSH", linewidth=1)
    fplot!(f2.(Nv), f1.(fixLR), ribbon=fixLR_S .* ScaleSTD, labels="fixLR", linewidth=4)

    #fplot!(xticks =collect( 10 .^(1:0.25:5))) 
    #fplot!(yticks =collect( 10.0.^ (-5:1:5))) 
    fplot!(yaxis=(:log10, [0.00001, :auto]))
    fplot!(xaxis=:log10)
    fplot!(gridlinewidth=2)
    display(
        fplot!(labels="teig_LR: fix Point",
            xlabel=L"number of steps, (N)", ylabel=L"CPU-time", legend=:bottomright,
            xticks=10.0 .^ (1:6), yticks=10.0 .^ (-5:3))
    )

    savefig("myplot_beam_feedback_N100.svg")
end


tmake_SH_PRi_S .= [-maximum([-a * 0.99, -b]) for (a, b) in zip(tmake_SH_PRi, tmake_SH_PRi_S)]
tmake_SH_Ci_S .= [-maximum([-a * 0.99, -b]) for (a, b) in zip(tmake_SH_Ci, tmake_SH_Ci_S)]
tmake_SH_PhiALL_S .= [-maximum([-a * 0.99, -b]) for (a, b) in zip(tmake_SH_PhiALL, tmake_SH_PhiALL_S)]
teig_SH_S .= [-maximum([-a * 0.99, -b]) for (a, b) in zip(teig_SH, teig_SH_S)]
fixSH_S .= [-maximum([-a * 0.99, -b]) for (a, b) in zip(fixSH, fixSH_S)]

tmake_LR_S .= [-maximum([-a * 0.99, -b]) for (a, b) in zip(tmake_LR, tmake_LR_S)]
teig_LR_S .= [-maximum([-a * 0.99, -b]) for (a, b) in zip(teig_LR, teig_LR_S)]
fixLR_S .= [-maximum([-a * 0.99, -b]) for (a, b) in zip(fixLR, fixLR_S)]



#tmake_SH_PRi_S[isnan.(tmake_SH_PRi_S)] .= 0.0
#tmake_SH_Ci_S[isnan.(tmake_SH_Ci_S)] .= 0.0
#tmake_SH_PhiALL_S[isnan.(tmake_SH_PhiALL_S)] .= 0.0
#teig_SH_S[isnan.(teig_SH_S)] .= 0.0
#tfixP_SH_S[isnan.(tfixP_SH_S)] .= 0.0
#tmake_LR_S[isnan.(tmake_LR_S)] .= 0.0
#teig_LR_S[isnan.(teig_LR_S)] .= 0.0
#tfixP_LR_S[isnan.(tfixP_LR_S)] .= 0.0


f1 = (x) -> x == 0.0 ? NaN : x
f2 = (x) -> x == 0.0 ? NaN : x

ScaleSTD = 0.5

#fplot=scatter
#fplot! = scatter!
fplot = plot
fplot! = plot!

#fplot(f2.(Nv), f1.(tmake_SH_PRi .+tmake_SH_PRi_S), labels="tmake_SH_PRi")
#fplot!(f2.(Nv), f1.(tmake_SH_PRi ), labels="tmake_SH_PRi")
#fplot!(f2.(Nv), f1.(tmake_SH_PRi .-tmake_SH_PRi_S), labels="tmake_SH_PRi")
fplot(f2.(Nv), f1.(tmake_SH_PRi), ribbon=tmake_SH_PRi_S .* ScaleSTD, labels="tmake_SH_PRi")
fplot!(f2.(Nv), f1.(tmake_SH_Ci), ribbon=tmake_SH_Ci_S .* ScaleSTD, labels="tmake_SH_Ci")
#fplot!(f2.(Nv), f1.(tmake_SH_Ci),ribbon=tmake_SH_Ci_S .* ScaleSTD, labels="tmake_SH_Ci")
fplot!(f2.(Nv), f1.(tmake_SH_PhiALL), ribbon=tmake_SH_PhiALL_S .* ScaleSTD, labels="tmake_SH_Phi_1step...prod(...)")
fplot!(f2.(Nv), f1.(teig_SH), ribbon=teig_SH_S .* ScaleSTD, labels="teig_SH: eigs(PHI) only")
#fplot!(f2.(Nv), f1.(tfixP_SH),ribbon=tfixP_SH_S .* ScaleSTD, labels="teig_SH: fix Point")
fplot!(f2.(Nv), f1.(tmake_LR), ribbon=tmake_LR_S .* ScaleSTD, labels="tmake_LR", linewidth=2)
fplot!(f2.(Nv), f1.(teig_LR), ribbon=teig_LR_S .* ScaleSTD, labels="teig_LR: eigs(ΦR,ΦL) only", linewidth=2)


fplot!(f2.(Nv), f1.(fixSH), ribbon=fixSH_S .* ScaleSTD, labels="fixSH", linewidth=1)
fplot!(f2.(Nv), f1.(fixLR), ribbon=fixLR_S .* ScaleSTD, labels="fixLR", linewidth=2)



#fplot!(xticks =collect( 10 .^(1:0.25:5))) 
#fplot!(yticks =collect( 10.0.^ (-5:1:5))) 
fplot!(yaxis=(:log10, [0.00001, :auto]))
fplot!(xaxis=:log10)
fplot!(gridlinewidth=2)



function power_fit(x::Vector{Float64}, y::Vector{Float64}, N)
    ynew = deepcopy(y[x.>N])
    xnew = deepcopy(x[x.>N])
    xnew = (xnew[.!isnan.(ynew)])
    ynew = (ynew[.!isnan.(ynew)])
    logx = log.(xnew)
    logy = log.(ynew)

    # Linear fit to the transformed problem
    X = hcat(ones(length(logx)), logx)
    β = X \ logy  # Solves for log(c) and p

    c = exp(β[1])  # Convert back from log(c) to c
    p = β[2]

    # Calculate fitted y values using the original x values
    fitted_y = c .* xnew .^ p
    println([c, p])
    return xnew, fitted_y #c, p, 
end

fplot!(power_fit(Nv, tmake_SH_PRi, 300)..., labels="", linewidth=1, linecolor=:black, linestyle=:dash)
fplot!(power_fit(Nv, tmake_SH_Ci, 300)..., labels="", linewidth=1, linecolor=:black, linestyle=:dash)
fplot!(power_fit(Nv, tmake_SH_PhiALL, 300)..., labels="", linewidth=1, linecolor=:black, linestyle=:dash)
fplot!(power_fit(Nv, teig_SH, 300)..., labels="", linewidth=1, linecolor=:black, linestyle=:dash)

fplot!(power_fit(Nv, tmake_LR, 300)..., labels="", linewidth=1, linecolor=:black, linestyle=:dash)
fplot!(power_fit(Nv, teig_LR, 300)..., labels="", linewidth=1, linecolor=:black, linestyle=:dash)

fplot!(power_fit(Nv, fixSH, 300)..., labels="", linewidth=1, linecolor=:black, linestyle=:dash)
fplot!(power_fit(Nv, fixLR, 300)..., labels="", linewidth=1, linecolor=:black, linestyle=:dash)


display(
    fplot!(labels="teig_LR: fix Point",
        xlabel=L"number of steps, (N)", ylabel=L"CPU-time", legend=:bottomright,
        xticks=10.0 .^ (1:6), yticks=10.0 .^ (-5:3))
)



savefig("myplot_beam_feedback_eigs_test_N100.svg")