using LinearAlgebra
using StaticArrays
using SemiDiscretizationMethod

function beam_matrices(E, A, ρ, L, η, N)
    dx = L / N
    m = ρ * A * dx
    k = E * A / dx

    # Mass matrix
    M = Diagonal(fill(m, N))

    # Stiffness matrix
    K = zeros(N, N)
    for i in 1:N
        if i > 1
            K[i, i] += k
            K[i, i-1] -= k
        end
        if i < N
            K[i, i] += k
            K[i, i+1] -= k
        end
    end

    # Apply fixed boundary condition at one end (first node)
    K = K[2:end, 2:end]
    M = M[2:end, 2:end]
    C = K * η
    return M, C, K
end

function first_order_system(M, C, K)
    n = size(M, 1)
    Z = zeros(n, n)
    In = Matrix(I, n, n)

    A = [Z In;
        -M\K -M\C]

    return A
end

function delay_beam_feedback(E, A, ρ, L, η, N, P, τpTspeed, Tspeed, TperpTspeed, ::Val{Nm1}) where {Nm1}
    T = TperpTspeed * Tspeed
    M, C, K = beam_matrices(E, A, ρ, L, η, N)
    A = first_order_system(M, C, K)

    B = zeros(2 * (Nm1), 2 * (Nm1))
    B[2*(N-1), 1] = P * E / (L / N)  # feedback from displacement of the first node to velocity of the last node

    F = zeros(Nm1)
    node = Nm1
    F[node] = 1.0
    Ffirst = vcat(zeros(Nm1), F)

    #AMx = ProportionalMX(t -> SMatrix{2 * (Nm1),2 * (Nm1)}(A)::SMatrix{2 * (Nm1),2 * (Nm1),Float64})
    #τ1 = t -> τpTspeed * Tspeed::Float64
    #BMx1 = DelayMX(τ1, t -> SMatrix{2 * (Nm1),2 * (Nm1)}(B.* (t<T*0.8))::SMatrix{2 * (Nm1),2 * (Nm1),Float64})
    #cVec = Additive(t -> SVector{2 * (Nm1)}((1.0 + 0.0 * cos(2π / T * t * 4)) .* Ffirst)::SVector{2 * (Nm1),Float64})


 #   AMx = ProportionalMX(t -> MMatrix{2 * (Nm1),2 * (Nm1)}(A)::MMatrix{2 * (Nm1),2 * (Nm1),Float64})
 #   τ1 = t -> τpTspeed * Tspeed::Float64
 #   BMx1 = DelayMX(τ1, t -> MMatrix{2 * (Nm1),2 * (Nm1)}(B .* (t<T*0.8))::MMatrix{2 * (Nm1),2 * (Nm1),Float64})
 #   cVec = Additive(t -> MVector{2 * (Nm1)}((0.0 + 1.0 * cos(2π / T * t )) .* Ffirst)::MVector{2 * (Nm1),Float64})



    AMx = ProportionalMX(t ->A)
    τ1 = t -> τpTspeed * Tspeed::Float64
    BMx1 = DelayMX(τ1, t ->B.* (t<T*0.8))
    cVec = Additive(t -> (1.0 + 0.0 * cos(2π / T * t * 4)) .* Ffirst)

    #LDDEProblem(AMx, [BMx1, BMx2], cVec)
    LDDEProblem(AMx, BMx1, cVec)
end

function foo(τpTspeed, P)
    # Example parameters
    E = Float64(210e9)           # Young's modulus (Pa)
    A = 1e-4            # Cross-sectional area (m²)
    ρ = 7800.0            # Density (kg/m³)
    L = 100.0             # Length of beam (m)
    η = 0.01
    N = 15              # Number of discrete elements

    c = sqrt(E / ρ)

    Tspeed = L / c
    τpTspeed = 0.2#.955#τ = τpT * Tspeed
    P = 0.2

    TperpTspeed = 0.4
    Tper = TperpTspeed * Tspeed

    prob_lddep = delay_beam_feedback(E, A, ρ, L, η, N, P, τpTspeed, Tspeed, TperpTspeed, Val(N - 1))

    Ndisc = 10
    method = SemiDiscretization(0, Tper / Ndisc) # 3rd order semi discretization with Δt=0.1

   @show  Nsteps = Int((Tper + 100eps(Tper)) ÷ method.Δt)
   @show r_approx= τpTspeed * Tspeed / method.Δt
    #mappingLR = DiscreteMapping_LR(mathieu_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true)#The discrete mapping of the system
    #μLR = spectralRadiusOfMapping(mappingLR)

    mapping = DiscreteMapping_LR(prob_lddep, method, τpTspeed * Tspeed, n_steps=Nsteps, calculate_additive=true)#The discrete mapping of the system
    return μ = spectralRadiusOfMapping(mapping)

end
println("-----------------------")
@time foo(0.999995, -0.7)
@time foo(0.999995, -0.7)
@time foo(1.0, 0.02)


 # Example parameters
 E = Float64(210e9)            # Young's modulus (Pa)
 A = 1e-4            # Cross-sectional area (m²)
 ρ = 7800.0            # Density (kg/m³)
 L = 100.0             # Length of beam (m)
 η = 0.01
 N = 15              # Number of discrete elements

 c = sqrt(E / ρ)

 Tspeed = L / c
   τpTspeed = 0.2#0.955#τ = τpT * Tspeed
    P = 0.2

 TperpTspeed = 0.4
 Tper = TperpTspeed * Tspeed
 
 
 println("modell creation:")
 @time prob_lddep = delay_beam_feedback(E, A, ρ, L, η, N, P, τpTspeed, Tspeed, TperpTspeed, Val(N - 1))
 
 @code_warntype prob_lddep.A(1.2)
 @code_warntype prob_lddep.Bs[1](1.2)
 @code_warntype prob_lddep.Bs[1].τ(1.2)
 
 Ndisc = 40
 method = SemiDiscretization(0, Tper / Ndisc) # 3rd order semi discretization with Δt=0.1
 
 Nsteps = Int((Tper + 100eps(Tper)) ÷ method.Δt)
 τmax= τpTspeed * Tspeed
 
 @show  Nsteps = Int((Tper + 100eps(Tper)) ÷ method.Δt)
 @show r_approx= τmax / method.Δt
 #mappingLR = DiscreteMapping_LR(mathieu_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true)#The discrete mapping of the system
 #μLR = spectralRadiusOfMapping(mappingLR)

 @time mapping = DiscreteMapping_LR(prob_lddep, method, τmax, n_steps=Nsteps, calculate_additive=true);#The discrete mapping of the system
 @show @time μ = spectralRadiusOfMapping(mapping)
 
 