m_lower = [] # lower bound of mass bin (kg)
m_upper = [] # upper bounf of mass bin (kg)

mlow = 10^(-21)/2
mup = 10^(-21)

for i in 1:30
    mlow = mlow*2
    push!(m_lower,mlow)
end

for i in 1:30
    mup = mup*2
    push!(m_upper,mup)
end

m_k = [] # avg mass of one particle in size bin k (kg)
mass = 0
for i in 1:length(m_lower)
    mass = (m_lower[i] + m_upper[i])/2
    push!(m_k,mass)
end

ρ = 1.8*10^(-3) # # particle density (kg/cm^3)

Dp_lower = [] # lower limit of particle diamter (cm)
Dp_upper = [] # upper limit of particle diameter (cm)

Dplow = 0
Dpup = 0

for x in m_lower
    Dplow = ((3*x/(4*ρ*pi))^(1/3))*2 #*10^4 
    push!(Dp_lower,Dplow)
end

for x in m_upper
    Dpup = ((3*x/(4*ρ*pi))^(1/3))*2 #*10^4 
    push!(Dp_upper,Dpup)
end
dDp = [] # size range (cm)
sizerange = 0
for i in 1:length(Dp_lower)
    sizerange = Dp_upper[i] - Dp_lower[i]
    push!(dDp,sizerange)
end

Dp = [] # avg diameter of size bin k (cm)
avgdia = 0
for i in 1:length(Dp_lower)
    avgdia = (Dp_upper[i] + Dp_lower[i])/2
    push!(Dp,avgdia)
end

μ = 1.83*10^(-7) # viscosity of air at T = 298 K (kg cm-1 s-1)
λ = 0.0686*10^(-4) # mean free path of air at T = 298 K (cm)
T = 298 # Temperature (K)
k = 1.38*10^(-19) # Boltzmann constant(cm^2 kg s-2 K-1)

# parameters needed to calculate coagulation coefficient
C_c = []
cc = 0
for i in 1:30
    if Dp[i] < 0.1 # free molecular regime
        cc = 1+1.657*2*λ/Dp[i]
    else # continuum regime
        cc = 1
    end
    push!(C_c,cc)
end

D = [] # diffusivity cm^2 s^-1
Di = 0
for i in 1:30
    Di = k*T*C_c[i]/(3*pi*μ*Dp[i])
    push!(D,Di)
end

c = []
ci = 0
for i in 1:30
    ci = sqrt(8*k*T/(pi*m_k[i]))
    push!(c,ci)
end

l = []
li = 0
for i in 1:30
    li = 8*D[i]/(pi*c[i])
    push!(l,li)
end

g = []
gi = 0
for i in 1:30
    gi = sqrt(2)/(3*Dp[i]*l[i])*((Dp[i]+l[i])^3-(Dp[i]^2+l[i]^2)^(3/2))-Dp[i]
    push!(g,gi)
end

# coagulation coef K_jk
function coagcoef(i,j) # (cm^3 s^-1) 
    if Dp[i] < 0.1 && Dp[j] < 0.1 # free molecular regime
        # Funchs Form of the Brownian Coagulation Coef
        return 2*pi*(Dp[i]+Dp[j])*(D[i]+D[j])*((Dp[i]+Dp[j])/(Dp[i]+Dp[j]+2*sqrt((g[i])^2+(g[j])^2))+8*(D[i]+D[j])/(sqrt((c[i])^2+(c[j])^2)*(Dp[i]+Dp[j])))^(-1)
    else # continuum regime
        return  2*pi*(Dp[i]+Dp[j])*(D[i]+D[j])
    end
end
  
import Pkg
Pkg.add("Distributions")

ξ = 1.0625
I = 30 # tot number of size bin

using ModelingToolkit
using Distributions
using DifferentialEquations
using Plots

@variables t 
@variables N(t)[1:I]
@variables M(t)[1:I]
@variables K(t)[1:I,1:I]
@variables f(t)[1:I]
@variables ψ(t)[1:I]
@parameters T = 298
@constants ξ = 1.0625 
N = collect(N)
K = collect(K)
f = collect(f)
ψ = collect(ψ)
M = collect(M)
ddt = Differential(t)

eqs = [
        reduce(vcat,[K[i,j] ~ coagcoef(i,j) for i in 1:I,j in 1:I])
       [f[i] ~ 2*N[i]/m_lower[i]*(2-m_k[i]/m_lower[i]) for i in 1:I]
       [ψ[i] ~ 2*N[i]/m_lower[i]*(m_k[i]/m_lower[i]-1) for i in 1:I]
       [ddt(N[1]) ~ - K[1,1]*(N[1])^2 - N[1]*sum(K[1,2:I].*N[2:I])]
       [ddt(N[2]) ~ 1/2*K[1,1]*N[1]^2 - K[2,2]*N[2]^2 - N[2]*sum(K[2,2:30].*N[2:30]) - ψ[2]*K[2,1]*M[1] + (ψ[2]-f[2])/(2*m_lower[2])*ξ*K[2,1]*M[1]*m_k[1]]
       [ddt(N[i]) ~ 1/2*K[i-1,i-1]*(N[i-1])^2 - K[i,i]*(N[i])^2 - N[i]*sum(K[i,i+1:I].*N[i+1:I]) + ψ[i-1]*sum(K[i-1,1:i-2].*M[1:i-2]) - ψ[i]*sum(K[i,1:i-1].*M[1:i-1]) + (ψ[i]-f[i])/(2*m_lower[i])*ξ*sum(K[i,1:i-1].*M[1:i-1].*m_k[1:i-1]) - (ψ[i-1]-f[i-1])/(2*m_lower[i-1])*ξ*sum(K[i-1,1:i-2].*M[1:i-2].*m_k[1:i-2]) for i in 3:29]
       [ddt(N[30]) ~ 1/2*K[29,29]*N[29]^2 - K[30,30]*N[30]^2 + ψ[29]*sum(K[29,1:28].*M[1:28]) - ψ[30]*sum(K[30,1:29].*M[1:29]) + (ψ[30]-f[30])/(2*m_lower[30])*ξ*sum(K[30,1:29].*M[1:29].*m_k[1:29]) - (ψ[29]-f[29])/(2*m_lower[29])*ξ*sum(K[29,1:28].*M[1:28].*m_k[1:28])]
       [ddt(M[1]) ~ - K[1,1]*N[1]*M[1] - M[1]*sum(K[1,2:I].*N[2:I])]
       [ddt(M[2]) ~ K[1,1]*N[1]*M[1] - K[2,2]*N[2]*M[2] + N[2]*K[2,1]*M[1] - M[2]*sum(K[2,3:30].*N[3:30]) - ψ[2]*m_lower[3]*K[2,1]*M[1] - f[2]/2*ξ*K[2,1]*M[1]*m_k[1] + (ψ[2]-f[2])/(2*m_lower[2])*ξ^3*K[2,1]*M[1]*m_k[1]^2]
       [ddt(M[i]) ~ K[i-1,i-1]*N[i-1]*M[i-1] - K[i,i]*N[i]*M[i] + N[i]*sum(K[i,1:i-1].*M[1:i-1]) - M[i]*sum(K[i,i+1:I].*N[i+1:I]) + ψ[i-1]*m_lower[i]*sum(K[i-1,1:i-2].*M[1:i-2]) - ψ[i]*m_lower[i+1]*sum(K[i,1:i-1].*M[1:i-1]) + f[i-1]/2*ξ*sum(K[i-1,1:i-2].*M[1:i-2].*m_k[1:i-2]) - f[i]/2*ξ*sum(K[i,1:i-1].*M[1:i-1].*m_k[1:i-1]) + (ψ[i]-f[i])/(2*m_lower[i])*ξ^3*sum(K[i,1:i-1].*M[1:i-1].*m_k[1:i-1].*m_k[1:i-1]) - (ψ[i-1]-f[i-1])/(2*m_lower[i-1])*ξ^3*sum(K[i-1,1:i-2].*M[1:i-2].*m_k[1:i-2].*m_k[1:i-2]) for i in 3:29]
       [ddt(M[30]) ~ K[29,29]*N[29]*M[29] - K[30,30]*N[30]*M[30] + N[30]*sum(K[30,1:29].*M[1:29]) + ψ[29]*m_lower[30]*sum(K[29,1:28].*M[1:28]) + f[29]/2*ξ*sum(K[29,1:28].*M[1:28].*m_k[1:28]) - f[30]/2*ξ*sum(K[30,1:29].*M[1:29].*m_k[1:29]) + (ψ[30]-f[30])/(2*m_lower[30])*ξ^3*sum(K[30,1:29].*M[1:29].*m_k[1:29].*m_k[1:29]) - (ψ[29]-f[29])/(2*m_lower[29])*ξ^3*sum(K[29,1:28].*M[1:28].*m_k[1:28].*m_k[1:28])]
    ] 

@named tomas = ODESystem(eqs)

states(tomas)

tomas_simplified = structural_simplify(tomas)

states(tomas_simplified)

# initial conditions
initial_valsN = [N[i] => pdf(LogNormal(15,5),i)*10^8 for i in 1:I] # unit!
initial_valsM = [M[i] => m_k[i]*(pdf(LogNormal(15,5),i)*10^8) for i in 1:I]
ps = [T => 298]
# time span
tspan = (0,86400) # run for a day

prob = ODEProblem(tomas_simplified,[initial_valsN;initial_valsM],tspan,ps)

sol = solve(prob,Tsit5())

plot(sol,xlabel="time", ylabel="Number/Mass distribution")
