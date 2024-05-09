begin
	import Pkg
	Pkg.add("Distributions")
	Pkg.add("DataFrames")
	Pkg.add("DifferentialEquations")
	using DataFrames
	using Distributions
	using DifferentialEquations
end

begin
	# lower bound of mas bin (g)
	m_lower = [] 
	# upper bound of mass bin (g)
	m_upper = []
	
	mlow = 10^(-18)/2 # (g)
	mup = 10^(-18) # (g)
	# number of bins
	bin = 5
		
	for i in 1:bin
	    mlow = mlow*2
	    push!(m_lower,mlow)
	end
	
	for i in 1:bin
	    mup = mup*2
	    push!(m_upper,mup)
	end
	
	# avg mass of one particle in size bin k (g)
	m_k = []
	mass = 0
	for i in 1:length(m_lower)
	    mass = (m_lower[i] + m_upper[i])/2
	    push!(m_k,mass)
	end
	# print(m_lower)
	# print(m_upper)
	# print(m_k)
end

begin
	# particle density (g cm-3)
	ρ = 1.8 
	# lower limit of particle diamter (cm)
	Dp_lower = []
	# upper limit of particle diameter (cm)
	Dp_upper = [] 
	
	Dplow = 0
	Dpup = 0
	
	for x in m_lower
	    Dplow = ((3*x/(4*ρ*pi))^(1/3))*2 
	    push!(Dp_lower,Dplow)
	end
	
	for x in m_upper
	    Dpup = ((3*x/(4*ρ*pi))^(1/3))*2 
	    push!(Dp_upper,Dpup)
	end

end

begin
	# diamter bin width (cm)
	dDp = [] 
	sizerange = 0
	for i in 1:length(Dp_lower)
	    sizerange = Dp_upper[i] - Dp_lower[i]
	    push!(dDp,sizerange)
	end

	# mass bin width (g)
	dm = [] 
	msizerange = 0
	for i in 1:length(m_lower)
	    msizerange = m_upper[i] - m_lower[i]
	    push!(dm,msizerange)
	end

	# average diameter of size bin k (cm)
	Dp = Float64[] 
	avgdia = 0
	for i in 1:length(Dp_lower)
	    avgdia = (Dp_upper[i] + Dp_lower[i])/2
	    push!(Dp,avgdia)
	end
	print(Dp)
end

begin
    # viscosity of air at T = 298 K
    μ = 1.83*10^(-4) # (g cm-1 s-1)
    # mean free path of air at T = 298 K
    λ = 0.0686*10^(-4) # (cm)
    # Temperature
    T = 298 # (K)
    # Boltzmann constant
    k = 1.38*10^(-16) #(cm2 g s-2 K-1)
    
    # parameters needed to calculate coagulation coefficient
        
    # slip correction
    C_c = Float64[]
    cc = 0
    for i in 1:bin
        if Dp[i] < 10^(-5) # (cm) # free molecular regime (0.1 um)
            cc = 1 + 2*λ / Dp[i]*(1.257+0.4*exp(-1.1 * Dp[i]/(2*λ)))
        else # continuum regime
            cc = 1
        end
        push!(C_c,cc)
    end
    print(C_c)
        
    D = Float64[] # diffusivity cm^2 s^-1
    Di = 0
    for i in 1:bin
        Di = k*T*C_c[i]/(3*pi*μ*Dp[i])
        push!(D,Di)
    end
    c = Float64[]
    ci = 0
    for i in 1:bin
        ci = sqrt(8*k*T/(pi*m_k[i]))
        push!(c,ci)
    end
    
    l = Float64[]
    li = 0
    for i in 1:bin
        li = 8*D[i]/(pi*c[i])
        push!(l,li)
    end
    
    g = Float64[]
    gi = 0
    for i in 1:bin
        gi = sqrt(2)/(3*Dp[i]*l[i])*((Dp[i]+l[i])^3-(Dp[i]^2+l[i]^2)^(3/2))-Dp[i]
        push!(g,gi)
    end
end    

# coagulation coef K_jk
function coagcoef(i,j) # (cm^3 s^-1) 
    if Dp[i] < 10^(-5) && Dp[j] < 10^(-5) # free molecular regime
        # Funchs Form of the Brownian Coagulation Coef
        return 2*pi*(Dp[i]+Dp[j])*(D[i]+D[j])*((Dp[i]+Dp[j])/(Dp[i]+Dp[j]+2*sqrt((g[i])^2+(g[j])^2))+8*(D[i]+D[j])/(sqrt((c[i])^2+(c[j])^2)*(Dp[i]+Dp[j])))^(-1)
    else # continuum regime
        return  2*k*T/(3*μ)*(Dp[i]+Dp[j])^2/(Dp[i]*Dp[j])
    end
end

using Plots

begin
	ξ = 1.0625 # depends on the bin spacing (for mass-doubling spacing)
	I = bin # tot number of size bin
end

Kvec_eq = reduce(vcat, coagcoef(i,j) for i in 1:bin, j in 1:bin)

begin
	Kvec = zeros(I*I)
	for i in 1:I*I
		Kvec[i] = Kvec_eq[i]
		i += 1
	end
	Kvec
end

Kmat = zeros(I,I)

# Fill the matrix from the vector
for j in 1:I
    Kmat[:, j] = Kvec[(j-1)*I+1:j*I]
end

begin
	using ModelingToolkit
	@variables t 
	@variables N(t)[1:I]
	@variables M(t)[1:I]
	@variables f(t)[1:I]
	@variables ψ(t)[1:I]
	# @parameters Kmat
	# @constants ξ = 1.0625 
	N = collect(N)
	f = collect(f)
	ψ = collect(ψ)
	M = collect(M)
	ddt = Differential(t)
end

eqs = [
    [f[i] ~ 2 * N[i] / m_lower[i] * (2 - m_k[i] / m_lower[i]) for i in 1:I]
    [ψ[i] ~ 2 * N[i] / m_lower[i] * (m_k[i] / m_lower[i] - 1) for i in 1:I]
    [ddt(N[k]) ~ (k > 1 ? 1/2 * Kmat[k-1,k-1] * N[k-1]^2 : 0) -
                Kmat[k,k] * N[k]^2 - N[k] * sum(Kmat[k, k+1:I] .* N[k+1:I]) +
                (k > 2 ? ψ[k-1] * sum(Kmat[k-1, 1:k-2] .* M[1:k-2]) : 0) +
                (k > 1 ? -ψ[k] * sum(Kmat[k, 1:k-1] .* M[1:k-1]) : 0) +
                (k > 1 ? (ψ[k] - f[k]) / (2 * m_lower[k]) * ξ * sum(Kmat[k, 1:k-1] .* M[1:k-1] .* m_k[1:k-1]) : 0)+(k > 2 ? -(ψ[k-1] - f[k-1]) / (2 * m_lower[k-1]) * ξ * sum(Kmat[k-1, 1:k-2] .* M[1:k-2] .* m_k[1:k-2]) : 0) for k in 1:I]
    [ddt(M[k]) ~ (k > 1 ? Kmat[k-1,k-1] * N[k-1] * M[k-1] : 0) -
                Kmat[k,k] * N[k] * M[k] + 
                (k > 1 ? N[k] * sum(Kmat[k, 1:k-1] .* M[1:k-1]) : 0) + 
                (k <= I-1 ? -M[k] * sum(Kmat[k, k+1:I] .* N[k+1:I]) : 0) +
                (k > 2 ? ψ[k-1] * m_lower[k] * sum(Kmat[k-1, 1:k-2] .* M[1:k-2]) : 0) + (k > 1 && k <= I-1 ? -ψ[k] * m_lower[k+1] * sum(Kmat[k, 1:k-1] .* M[1:k-1]) : 0) +
                (k > 2 ? f[k-1] / 2 * ξ * sum(Kmat[k-1, 1:k-2] .* M[1:k-2] .* m_k[1:k-2]) : 0) +
                (k > 1 ? -f[k] / 2 * ξ * sum(Kmat[k, 1:k-1] .* M[1:k-1] .* m_k[1:k-1]) : 0) +
                (k > 1 ? (ψ[k] - f[k]) / (2 * m_lower[k]) * ξ^3 * sum(Kmat[k, 1:k-1] .* M[1:k-1] .* (m_k[1:k-1] .^ 2)) : 0) +
                (k > 2 ? -(ψ[k-1] - f[k-1]) / (2 * m_lower[k-1]) * ξ^3 * sum(Kmat[k-1, 1:k-2] .* M[1:k-2] .* (m_k[1:k-2] .^ 2)) : 0) for k in 1:I
]
]
@named tomas = ODESystem(eqs)
states(tomas)
tomas_simplified = structural_simplify(tomas)
states(tomas_simplified)
begin
    # initial conditions
    initial_valsN = [N[i] => pdf(Normal(15,5),i)*10^10 for i in 1:I] 
    initial_valsM = [M[i] => m_k[i]*(pdf(Normal(15,5),i)*10^10) for i in 1:I]
    ps = [T => 298]
    # time span
    tspan = (0,86400) # run for a day
    
    prob = ODEProblem(tomas_simplified,[initial_valsN;initial_valsM],tspan,ps)
end

sol = solve(prob,Tsit5())

time = sol[t]
timesteps = length(time) # number of time steps
N_sol = sol[N]
M_sol = sol[M]
begin
	N_1= []
	M_1= []
	for i in 1:length(time)
		push!(N_1, N_sol[i][1])
		push!(M_1, M_sol[i][1])
	end
	N_2= []
	M_2= []
	for i in 1:length(time)
		push!(N_2, N_sol[i][2])
		push!(M_2, M_sol[i][2])
	end
	N_3= []
	M_3= []
	for i in 1:length(time)
		push!(N_3, N_sol[i][3])
		push!(M_3, M_sol[i][3])
	end
	N_4= []
	M_4= []
	for i in 1:length(time)
		push!(N_4, N_sol[i][4])
		push!(M_4, M_sol[i][4])
	end
	N_5= []
	M_5= []
	for i in 1:length(time)
		push!(N_5, N_sol[i][5])
		push!(M_5, M_sol[i][5])
	end
    N_6= []
	M_6= []
	for i in 1:length(time)
		push!(N_6, N_sol[i][6])
		push!(M_6, M_sol[i][6])
	end
	N_7= []
	M_7= []
	for i in 1:length(time)
		push!(N_7, N_sol[i][7])
		push!(M_7, M_sol[i][7])
	end
	N_8= []
	M_8= []
	for i in 1:length(time)
		push!(N_8, N_sol[i][8])
		push!(M_8, M_sol[i][8])
	end
	N_9= []
	M_9= []
	for i in 1:length(time)
		push!(N_9, N_sol[i][9])
		push!(M_9, M_sol[i][9])
	end
	N_10= []
	M_10= []
	for i in 1:length(time)
		push!(N_10, N_sol[i][10])
		push!(M_10, M_sol[i][10])
	end
	N_11= []
	M_11= []
	for i in 1:length(time)
		push!(N_11, N_sol[i][11])
		push!(M_11, M_sol[i][11])
	end
    N_12= []
	M_12= []
	for i in 1:length(time)
		push!(N_12, N_sol[i][12])
		push!(M_12, M_sol[i][12])
	end
	N_13= []
	M_13= []
	for i in 1:length(time)
		push!(N_13, N_sol[i][13])
		push!(M_13, M_sol[i][13])
	end
	N_14= []
	M_14= []
	for i in 1:length(time)
		push!(N_14, N_sol[i][14])
		push!(M_14, M_sol[i][14])
	end
	N_15= []
	M_15= []
	for i in 1:length(time)
		push!(N_15, N_sol[i][15])
		push!(M_15, M_sol[i][15])
	end
	
end

begin
    plot(time,N_1, xlabel="Time(s)", ylabel="Number of aerosols",yaxis=:log10, label="N1(t)")
    plot!(time,N_2, label="N2(t)")
    plot!(time,N_3, label="N3(t)")
    plot!(time,N_4, label="N4(t)")
    plot!(time,N_5, label="N5(t)")
    plot!(time,N_6, label="N6(t)")
    plot!(time,N_7, label="N7(t)")
    plot!(time,N_8, label="N8(t)")
    plot!(time,N_9, label="N9(t)")
    plot!(time,N_10, label="N10(t)")
    plot!(time,N_11, label="N11(t)")
    plot!(time,N_12, label="N12(t)")
    plot!(time,N_13, label="N13(t)")
    plot!(time,N_14, label="N14(t)")
    plot!(time,N_15, label="N15(t)")
end

begin
    plot(time,M_1, xlabel="Time (s)", ylabel="Mass of aerosols (kg)",yaxis=:log10, label="M1(t)")
    plot!(time,M_2, label="M2(t)")
    plot!(time,M_3, label="M3(t)")
    plot!(time,M_4, label="M4(t)")
    plot!(time,M_5, label="M5(t)")
    plot!(time,M_6, label="M6(t)")
    plot!(time,M_7, label="M7(t)")
    plot!(time,M_8, label="M8(t)")
    plot!(time,M_9, label="M9(t)")
    plot!(time,M_10, label="M10(t)")
    plot!(time,M_11, label="M11(t)")
    plot!(time,M_12, label="M12(t)")
    plot!(time,M_13, label="M13(t)")
    plot!(time,M_14, label="M14(t)")
    plot!(time,M_15, label="M15(t)")
end
