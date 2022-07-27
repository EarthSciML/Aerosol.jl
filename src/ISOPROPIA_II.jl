using Markdown, OrdinaryDiffEq, DiffEqBase, Plots, StochasticDiffEq, DifferentialEquations, ModelingToolkit, NonlinearSolve
export ISOPROPIA

@variables NH4_aq Na_aq H_aq Cl_aq NO3_aq SO4_aq HNO3_aq NH3_aq HCl_aq HSO4_aq OH_aq H2O_aq Ca_aq K_aq Mg_aq 
# unit: mol/kg(water)

@variables P_NH3 P_HNO3 P_HCl P_H2O
# unit: atm

@variables NH42SO4_s NH4HSO4_s NH43HSO42_s NH4NO3_s NH4Cl_s NaCl_s NaNO3_s NaHSO4_s Na2SO4_s CaSO4_s CaNO32_s CaCl2_s K2SO4_s KHSO4_s KNO3_s KCl_s MgSO4_s MgNO32_s MgCl2_s 
# unit: mol/m3(air)

@variables I W_
@parameters T a_w Na_i H2SO4_i NH3_i HNO3_i HCl_i Ca_i K_i Mg_i # m_all

begin
K1 = 6.067*10^5
K2 = 7.974*10^11
K3 = 4.319*10^-5
K4 = 1.569*10^-2
K5 = 24.016
K6 = 0.872
K7 = 8.680
K8 = 1.079*10^5
K9 = 2.507*10^15
K10 = 9.557*10^21
K11 = 1.015*10^-2
K12 = 5.764*10^1
K13 = 1.805*10^-5
K14 = 2.511*10^6
K15 = 2.1*10^5
K16 = 1.971*10^6
K17 = 2.5*10^3
K18 = 1.01.*10^-14
K19 = 4.799*10^-1
K20 = 1.817
K21 = 1.086*10^-16
K22 = 1.197*10
K23 = 3.766*10
K24 = 2.413*10^4
K25 = 4.199*10^-17
K26 = 1.383
K27 = 2.972*10

H4 = -9.589
H5 = -8.423
H6 = 14.075
H7 = -6.167
H8 = 36.798
H9 = -8.754
H10 = -1.347
H11 = 8.85
H12 = 13.79
H13 = -1.5
H14 = 29.17
H15 = 29.17
H16 = 30.20
H17 = 30.20
H18 = -22.52
H19 = 0.98
H20 = -2.65
H21 = -71
H22 = -8.22
H23 = -1.56
H24 = 0.79
H25 = -74.735
H26 = -2.87
H27 = -5.19

C4 = 45.807
C5 = 17.964
C6 = 19.388
C7 = 19.953
C11 = 25.14
C12 = -5.39
C13 = 26.92
C14 = 16.83
C15 = 16.83
C16 = 19.91
C17 = 19.91
C18 = 26.92
C19 = 39.75
C20 = 38.57
C21 = 2.40
C22 = 16.01
C23 = 16.90
C24 = 14.75
C25 = 6.025
C26 = 15.83
C27 = 54.40
end

function Equilibriumconst(K₀,T,H_group,c_group) # to calculate equilibrium constant at given temperature
	# input variables：
	# K₀ --> equilibrium constant at 298.15K(Tₒ), T --> given temperature,K, H_group --> ΔH₀/(R*T₀), c_group --> Δcₒ/R 

	# Equation: 
	# K = K₀*exp(-ΔH₀/(R*T₀)*(T₀/T-1)-Δcₒ/R*(1+log(T₀/T)-T₀/T))

	# K₀: equilibrium constant at 298.15K(Tₒ)
	# T: given temperature, K
	# R: universal gas constant 
	# ΔH₀:  enthalpy change of the reaction at 298.15K(Tₒ)
	# Δc₀: change of molar heat capacity of products minus reactants

	K = K₀*exp(-H_group*(298.15/T-1)-c_group*(1+log(298.15/T)-298.15/T))
end

function Equilibriumconst(K₀, T, H_group, c_group) # to calculate equilibrium constant at given temperature
    # input variables：
    # K₀ --> equilibrium constant at 298.15K(Tₒ), T --> given temperature,K, H_group --> ΔH₀/(R*T₀), c_group --> Δcₒ/R 

    # Equation: 
    # K = K₀*exp(-ΔH₀/(R*T₀)*(T₀/T-1)-Δcₒ/R*(1+log(T₀/T)-T₀/T))

    # K₀: equilibrium constant at 298.15K(Tₒ)
    # T: given temperature, K
    # R: universal gas constant 
    # ΔH₀:  enthalpy change of the reaction at 298.15K(Tₒ)
    # Δc₀: change of molar heat capacity of products minus reactants

    K = K₀ * exp(-H_group * (298.15 / T - 1) - c_group * (1 + log(298.15 / T) - 298.15 / T))
end

function cal_I(z_all, m_all) # to compute the ionic strength of the multicomponent solution
    # input variables: 
    # mi --> molalities of ions, zi --> valence of ions
    I = 0
    for i in 1:length(z_all)
        I += 1 / 2 * (m_all[i] * z_all[i]^2)
    end
    return I
end

function Watercontent(M, m) # to compute the water uptake of aerosols using ZSR correlation
    # input variables:
    # M --> vector of molar concentration of species i (mol* m^-3 air), m --> vector of molality of an aqueous binary solution of the i-th electrolyte with the same relative humidity
    W = sum(M ./ m) # mass concentration of aerosol water (kg* m^-3 air)
end

function m_aw(aw, k0, k1, k2, k3, k4, k5) # to compute the molality of binary solutions as a function of water activity 
    # input variables:
    # aw --> water activity, ki --> coefficients of m(aw) from the polynomial fit
    m = k0
    n = [k1, k2, k3, k4, k5]
    for i in 1:5
        m += n[i] * aw^i
    end
    return m
end

function DRH(DRH₀, g1, T) # to cpmpute the deliquescence relative humidity at given temperature
    # input variables:
    # DRH₀ --> deliquescence relative humidity at 298.15K, g1 --> -Mₛ*mₛ*Lₛ/1000R, T --> given temperature
    DRH = exp(DRH₀ * g1 * (1 / T - 1 / 298.15))
end

function log_gamma_0(z1, z2, q, I)
    C = 1 + 0.055 * q * exp(-0.023 * I^3)
    Γ_ = exp(-0.5107 * I^0.5 / (1 + C * I^0.5))
    B = 0.75 - 0.065 * q
    Γ0 = (1 + B * (1 + 0.1 * I)^q - B) * Γ_
    log_gamma_12_0 = z1 * z2 * log(Γ0)
    return log_gamma_12_0
end

# function F1(z1,z2,z2_,m2,m2_,q,q2_,I)
# 	Aγ = 0.511 # kg^0.5 mol^-0.5
# 	n = length(z2_)
# 	Y_21 = ((z1+z2)/2)^2*m2/I
# 	log_γ_12 = log_gamma_0(z1,z2,q,I)
# 	A_group = Aγ*I^0.5/(1+I^0.5)
# 	F1 = Y_21*log_γ_12+A_group*z1*z2*Y_21
# 	for i in 1:n
# 		Y_i1 = ((z1+z2_[i])/2)^2*m2_[i]/I
# 		log_γ_1i = log_gamma_0(z1,z2_[i],q2_[i],I)
# 		F1 += Y_i1*log_γ_1i+A_group*z1*z2_[i]*Y_i1
# 	end
# 	return F1

# end
# function F2(z1,z2,z1_,m1,m1_,q,q1_,I)
# 	Aγ = 0.511 # kg^0.5 mol^-0.5
# 	m = length(z1_)
# 	X_12 = ((z1+z2)/2)^2*m1/I
# 	log_γ_12 = log_gamma_0(z1,z2,q,I)
# 	A_group = Aγ*I^0.5/(1+I^0.5)
# 	F2 = X_12*log_γ_12+A_group*z1*z2*X_12
# 	for i in 1:m
# 		X_i2 = ((z1_[i]+z2)/2)^2*m1_[i]/I
# 		log_γ_i2 = log_gamma_0(z1_[i],z2,q1_[i],I)
# 		F2 += X_i2*log_γ_i2+A_group*z1_[i]*z2*X_i2
# 	end
# 	return F2
# end


function F(cation, anions, I)
    z1_ = F_input(cation, anions)[1]
    m1_ = F_input(cation, anions)[2]
    q1_ = F_input(cation, anions)[3]
    z2_ = F_input(cation, anions)[4]
    m2_ = F_input(cation, anions)[5]
    q2_ = F_input(cation, anions)[6]
    z1 = F_input(cation, anions)[7]
    z2 = F_input(cation, anions)[8]
    m1 = F_input(cation, anions)[9]
    m2 = F_input(cation, anions)[10]
    q = F_input(cation, anions)[11]

    Aγ = 0.511 # kg^0.5 mol^-0.5
    n = length(z2_)
    Y_21 = ((z1 + z2) / 2)^2 * m2 / I
    log_γ_12 = log_gamma_0(z1, z2, q, I)
    A_group = Aγ * I^0.5 / (1 + I^0.5)
    F1 = Y_21 * log_γ_12 + A_group * z1 * z2 * Y_21
    for i in 1:n
        Y_i1 = ((z1 + z2_[i]) / 2)^2 * m2_[i] / I
        log_γ_1i = log_gamma_0(z1, z2_[i], q2_[i], I)
        F1 += Y_i1 * log_γ_1i + A_group * z1 * z2_[i] * Y_i1
    end
    m = length(z1_)
    X_12 = ((z1 + z2) / 2)^2 * m1 / I
    F2 = X_12 * log_γ_12 + A_group * z1 * z2 * X_12
    for i in 1:m
        X_i2 = ((z1_[i] + z2) / 2)^2 * m1_[i] / I
        log_γ_i2 = log_gamma_0(z1_[i], z2, q1_[i], I)
        F2 += X_i2 * log_γ_i2 + A_group * z1_[i] * z2 * X_i2
    end
    return F1, F2
end

function F_input(cation, anions)
    # z1,z2,z1_,z2_,m1,m2,m1_,m2_,q,q1_,q2_,I
    z2_NO3 = [1, 1, 1, 2, 2, 1]
    m2_NO3 = [H_aq, Na_aq, K_aq, Ca_aq, Mg_aq, NH4_aq]
    q2_NO3 = [2.6, -0.39, -2.33, 0.93, 2.32, -1.15]
    z2_Cl = [1, 1, 1, 2, 2, 1]
    m2_Cl = [H_aq, Na_aq, K_aq, Ca_aq, Mg_aq, NH4_aq]
    q2_Cl = [6.0, 2.23, 0.92, 2.4, 2.90, 0.82]
    z2_SO4 = [1, 1, 1, 1, 2]
    m2_SO4 = [H_aq, NH4_aq, Na_aq, K_aq, Mg_aq]
    q2_SO4 = [-0.1, -0.25, -0.19, -0.25, 0.15]
    if cation == "Ca"
        z1 = 2
        m1 = Ca_aq
        z1_Ca = [1, 1]
        m1_Ca = [NO3_aq, Cl_aq]
        q1_Ca = [0.93, 2.4]
        if anions == "NO3"
            z2 = 1
            m2 = NO3_aq
            q = 0.93
            z1_ = [1]
            m1_ = [Cl_aq]
            q1_ = [2.4]
            z2_ = [1, 1, 1, 2]
            m2_ = [H_aq, Na_aq, K_aq, Mg_aq]
            q2_ = [2.6, -0.39, -2.33, 2.32]
        elseif anions == "Cl"
            z2 = 1
            m2 = Cl_aq
            q = 2.4
            z1_ = [1]
            m1_ = [NO3_aq]
            q1_ = [0.93]
            z2_ = [1, 1, 1, 2]
            m2_ = [H_aq, Na_aq, K_aq, Mg_aq]
            q2_ = [6.0, 2.23, 0.92, 2.90]
        else
            print("Wrong F anion input")
        end

    elseif cation == "Na"
        z1 = 1
        m1 = Na_aq
        z1_Na = [1, 1, 2]
        m1_Na = [NO3_aq, Cl_aq, SO4_aq]
        q1_Na = [-0.39, 2.23, -0.19]
        if anions == "NO3"
            z2 = 1
            m2 = NO3_aq
            q = q1_Na[1]
            z1_ = deleteat!(z1_Na, 1)
            m1_ = deleteat!(m1_Na, 1)
            q1_ = deleteat!(q1_Na, 1)
            z2_ = deleteat!(z2_NO3, 2)
            m2_ = deleteat!(m2_NO3, 2)
            q2_ = deleteat!(q2_NO3, 2)
        elseif anions == "Cl"
            z2 = 2
            m2 = Cl_aq
            q = q1_Na[2]
            z1_ = deleteat!(z1_Na, 2)
            m1_ = deleteat!(m1_Na, 2)
            q1_ = deleteat!(q1_Na, 2)
            z2_ = deleteat!(z2_Cl, 2)
            m2_ = deleteat!(m2_Cl, 2)
            q2_ = deleteat!(q2_Cl, 2)
        elseif anions == "SO4"
            z2 = 2
            m2 = SO4_aq
            q = q1_Na[3]
            z1_ = deleteat!(z1_Na, 3)
            m1_ = deleteat!(m1_Na, 3)
            q1_ = deleteat!(q1_Na, 3)
            z2_ = deleteat!(z2_SO4, 3)
            m2_ = deleteat!(m2_SO4, 3)
            q2_ = deleteat!(q2_SO4, 3)
        else
            print("Wrong F anion input")
        end

    elseif cation == "K"
        z1 = 1
        m1 = K_aq
        z1_K = [1, 1, 2]
        m1_K = [NO3_aq, Cl_aq, SO4_aq]
        q1_K = [-2.33, 0.92, -0.25]
        if anions == "NO3"
            z2 = 1
            m2 = NO3_aq
            q = q1_K[1]
            z1_ = deleteat!(z1_K, 1)
            m1_ = deleteat!(m1_K, 1)
            q1_ = deleteat!(q1_K, 1)
            z2_ = deleteat!(z2_NO3, 3)
            m2_ = deleteat!(m2_NO3, 3)
            q2_ = deleteat!(q2_NO3, 3)
        elseif anions == "Cl"
            z2 = 1
            m2 = Cl_aq
            q = q1_K[2]
            z1_ = deleteat!(z1_K, 2)
            m1_ = deleteat!(m1_K, 2)
            q1_ = deleteat!(q1_K, 2)
            z2_ = deleteat!(z2_Cl, 3)
            m2_ = deleteat!(m2_Cl, 3)
            q2_ = deleteat!(q2_Cl, 3)
        elseif anions == "SO4"
            z2 = 2
            m2 = SO4_aq
            q = q1_K[3]
            z1_ = deleteat!(z1_K, 3)
            m1_ = deleteat!(m1_K, 3)
            q1_ = deleteat!(q1_K, 3)
            z2_ = deleteat!(z2_SO4, 4)
            m2_ = deleteat!(m2_SO4, 4)
            q2_ = deleteat!(q2_SO4, 4)
        else
            print("Wrong F anion input")
        end

    elseif cation == "Mg"
        z1 = 2
        m1 = Mg_aq
        z1_Mg = [1, 1, 2]
        m1_Mg = [NO3_aq, Cl_aq, SO4_aq]
        q1_Mg = [2.32, 2.90, 0.15]
        if anions == "NO3"
            z2 = 1
            m2 = NO3_aq
            q = q1_Mg[1]
            z1_ = deleteat!(z1_Mg, 1)
            m1_ = deleteat!(m1_Mg, 1)
            q1_ = deleteat!(q1_Mg, 1)
            z2_ = deleteat!(z2_NO3, 5)
            m2_ = deleteat!(m2_NO3, 5)
            q2_ = deleteat!(q2_NO3, 5)
        elseif anions == "Cl"
            z2 = 1
            m2 = Cl_aq
            q = q1_Mg[2]
            z1_ = deleteat!(z1_Mg, 2)
            m1_ = deleteat!(m1_Mg, 2)
            q1_ = deleteat!(q1_Mg, 2)
            z2_ = deleteat!(z2_Cl, 5)
            m2_ = deleteat!(m2_Cl, 5)
            q2_ = deleteat!(q2_Cl, 5)
        elseif anions == "SO4"
            z2 = 2
            m2 = SO4_aq
            q = q1_Mg[3]
            z1_ = deleteat!(z1_Mg, 3)
            m1_ = deleteat!(m1_Mg, 3)
            q1_ = deleteat!(q1_Mg, 3)
            z2_ = deleteat!(z2_SO4, 5)
            m2_ = deleteat!(m2_SO4, 5)
            q2_ = deleteat!(q2_SO4, 5)
        else
            print("Wrong F anion input")
        end

    elseif cation == "H"
        z1 = 1
        m1 = H_aq
        z1_H = [1, 1, 1, 2]
        m1_H = [NO3_aq, Cl_aq, HSO4_aq, SO4_aq]
        q1_H = [2.6, 6.0, 8.0, -0.1]
        if anions == "NO3"
            z2 = 1
            m2 = NO3_aq
            q = q1_H[1]
            z1_ = deleteat!(z1_H, 1)
            m1_ = deleteat!(m1_H, 1)
            q1_ = deleteat!(q1_H, 1)
            z2_ = deleteat!(z2_NO3, 1)
            m2_ = deleteat!(m2_NO3, 1)
            q2_ = deleteat!(q2_NO3, 1)
        elseif anions == "Cl"
            z2 = 1
            m2 = Cl_aq
            q = q1_H[2]
            z1_ = deleteat!(z1_H, 2)
            m1_ = deleteat!(m1_H, 2)
            q1_ = deleteat!(q1_H, 2)
            z2_ = deleteat!(z2_Cl, 1)
            m2_ = deleteat!(m2_Cl, 1)
            q2_ = deleteat!(q2_Cl, 1)
        elseif anions == "HSO4"
            z2 = 1
            m2 = HSO4_aq
            q = q1_H[3]
            z1_ = deleteat!(z1_H, 3)
            m1_ = deleteat!(m1_H, 3)
            q1_ = deleteat!(q1_H, 3)
            z2_ = 0
            m2_ = 0
            q2_ = 0
        elseif anions == "SO4"
            z2 = 2
            m2 = SO4_aq
            q = q1_H[4]
            z1_ = deleteat!(z1_H, 4)
            m1_ = deleteat!(m1_H, 4)
            q1_ = deleteat!(q1_H, 4)
            z2_ = deleteat!(z2_SO4, 1)
            m2_ = deleteat!(m2_SO4, 1)
            q2_ = deleteat!(q2_SO4, 1)
        else
            print("Wrong F anion input")
        end

    elseif cation == "NH4"
        z1 = 1
        m1 = NH4_aq
        z1_NH4 = [1, 1, 2]
        m1_NH4 = [NO3_aq, Cl_aq, SO4_aq]
        q1_NH4 = [-1.15, 0.82, -0.25]
        if anions == "SO4"
            z2 = 2
            m2 = SO4_aq
            q = q1_NH4[3]
            z1_ = deleteat!(z1_NH4, 3)
            m1_ = deleteat!(m1_NH4, 3)
            q1_ = deleteat!(q1_NH4, 3)
            z2_ = deleteat!(z2_SO4, 2)
            m2_ = deleteat!(m2_SO4, 2)
            q2_ = deleteat!(q2_SO4, 2)
        elseif anions == "Cl"
            z2 = 1
            m2 = Cl_aq
            q = q1_NH4[2]
            z1_ = deleteat!(z1_NH4, 2)
            m1_ = deleteat!(m1_NH4, 2)
            q1_ = deleteat!(q1_NH4, 2)
            z2_ = deleteat!(z2_Cl, 6)
            m2_ = deleteat!(m2_Cl, 6)
            q2_ = deleteat!(q2_Cl, 6)
        elseif anions == "NO3"
            z2 = 1
            m2 = NO3_aq
            q = q1_NH4[1]
            z1_ = deleteat!(z1_NH4, 1)
            m1_ = deleteat!(m1_NH4, 1)
            q1_ = deleteat!(q1_NH4, 1)
            z2_ = deleteat!(z2_NO3, 6)
            m2_ = deleteat!(m2_NO3, 6)
            q2_ = deleteat!(q2_NO3, 6)
        else
            print("Wrong F anion input")
        end
    else
        print("Wrong F cation input")
    end
    return z1_, m1_, q1_, z2_, m2_, q2_, z1, z2, m1, m2, q
end

function ActivityCoefficient(T, cation, anions, I)
    Aγ = 0.511 # kg^0.5 mol^-0.5
    z1 = F_input(cation, anions)[7]
    z2 = F_input(cation, anions)[8]
    F1 = F(cation, anions, I)[1]
    F2 = F(cation, anions, I)[2]

    log_γ₁₂_T0 = exp(-Aγ * (z1 * z2 * I^0.5 / (1 + I^0.5)) + z1 * z2 / (z1 + z2) * (F1 / z1 + F2 / z2))
    A = -0.41 * I^0.5 / (1 + I^0.5) + 0.0391^0.92
    γ₁₂ = exp(log_γ₁₂_T0 * (1.125 - 0.005 * (T - 273.15)) - (0.125 - 0.005 * (T - 273.15)) * A)
    return γ₁₂
end

function Wateruptake_denominator(w)
    r = zeros(18)
    m = ones(18, 7)
    m[1, :] = [36.356, -165.66, 447.46, -673.55, 510.91, -155.56, 0]
    m[2, :] = [20.847, -97.599, 273.220, -422.120, 331.160, -105.450, 0]
    m[3, :] = [1.061, -0.101, 1.579 * 10^-2, -1.950 * 10^-3, 9.515 * 10^-5, -1.547 * 10^-6, 0]
    m[4, :] = [1061.51, -4748.97, 8096.16, -6166.16, 1757.47, 0, 0]
    m[5, :] = [1.2141 * 10^4, -5.1173 * 10^4, 8.1252 * 10^4, -5.7527 * 10^4, 1.5305 * 10^4, 0, 0]
    m[6, :] = [179.721, -721.266, 1161.03, -841.479, 221 / 943, 0, 0]
    m[7, :] = [-0.778, 177.74, -719.79, 1174.6, -863.44, 232.31, 0]
    m[8, :] = [12.166, -16.154, 0, 10.886, 0, -6.815, 0]
    m[9, :] = [11.505, -26.518, 34.937, -19.829, 0, 0, 0]
    m[10, :] = [0.9988, -2.6947 * 10^-2, 1.9610 * 10^-4, 2.8154 * 10^-5, -6.1359 * 10^-7, 0, 0]
    m[11, :] = [1.0614, -0.1014, 1.5796 * 10^-2, -1.9501 * 10^-3, 9.5147 * 10^-5, -1.5473 * 10^-6, 0]
    m[12, :] = [1.0084, -4.9390 * 10^-2, 8.888 * 10^-3, -2.1570 * 10^-3, 1.6170 * 10^-4, -1.99 * 10^-6, -1.142 * 10^-7]
    m[13, :] = [1.0052, -6.4840 * 10^-2, 3.519 * 10^-2, -1.319 * 10^-2, 1.9250 * 10^-3, -1.224 * 10^-4, 2.87 * 10^-6]
    m[14, :] = [0.9968, -2.969 * 10^-2, 1.735 * 10^-5, -3.253 * 10^-4, 3.571 * 10^-5, -9.787 * 10^-7, 0]
    m[15, :] = [0.9968, -2.611 * 10^-2, -1.599 * 10^-3, 1.355 * 10^-4, -2.3170 * 10^-6, -1.113 * 10^-8, 0]
    m[16, :] = [1.0053, -2.4991 * 10^-2, 4.4688 * 10^-4, 1.6453 * 10^-5, -3.8940 * 10^-7, -4.7668 * 10^-8, 1.3753 * 10^-9]
    m[17, :] = [1.0261, -4.9766 * 10^-2, 3.2757 * 10^-3, -2.4477 * 10^-4, 1.0766 * 10^-5, -1.8329 * 10^-7, 0]
    m[18, :] = [1.0088, -5.3730 * 10^-2, 1.4201 * 10^-3, -9.2484 * 10^-4, 2.2796 * 10^-4, -1.5445 * 10^-5, 0]
    for i in 1:18
        for j in 1:7
            r[i] += m[i, j] * w^(j - 1)
        end
    end
    return r
end
@register_symbolic Wateruptake_denominator(w)

function Wateruptake(numerator::Vector{Float64},a_w::Float64) 
	denominator = Wateruptake_denominator(a_w) #unit: mol/kg
	return sum(numerator ./ denominator)
end
	
@register_symbolic Wateruptake(numerator, a_w)

z_all = [1, 1, 1, 1, 1, 2, 1, 2, 1, 2, 1]

eqs = [
	# Ionic strength
	I ~ cal_I(z_all,[NH4_aq, Na_aq, H_aq, Cl_aq, NO3_aq, SO4_aq, OH_aq, Ca_aq, K_aq, Mg_aq, HSO4_aq])
	
	# Equilibrium Equations
	K1 ~ Ca_aq*(NO3_aq^2)*ActivityCoefficient(T,"Ca","NO3",I)^3
	K2 ~ Ca_aq*(Cl_aq^2)*ActivityCoefficient(T,"Ca","Cl",I)^3
	K3 ~ Ca_aq*SO4_aq
	Equilibriumconst(K4,T,H4,C4) ~ (K_aq^2)*SO4_aq*ActivityCoefficient(T,"K","SO4",I)^3
	Equilibriumconst(K5,T,H5,C5) ~ K_aq*HSO4_aq*((ActivityCoefficient(T,"H","HSO4",I)*ActivityCoefficient(T,"K","Cl",I)/ActivityCoefficient(T,"H","Cl",I))^0.5)^2
	Equilibriumconst(K6,T,H6,C6) ~ K_aq*NO3_aq*ActivityCoefficient(T,"K","NO3",I)^2
	Equilibriumconst(K7,T,H7,C7) ~ K_aq*Cl_aq*ActivityCoefficient(T,"K","Cl",I)^2
	K8 ~ Mg_aq*SO4_aq*ActivityCoefficient(T,"Mg","SO4",I)^2
	K9 ~ Mg_aq*(NO3_aq^2)*ActivityCoefficient(T,"Mg","NO3",I)^3
	K10 ~ Mg_aq*(Cl_aq^2)*ActivityCoefficient(T,"Mg","Cl",I)^3
	Equilibriumconst(K11,T,H11,C11) ~ (H_aq*SO4_aq/HSO4_aq)*(ActivityCoefficient(T,"H","SO4",I)^3)/(ActivityCoefficient(T,"H","HSO4",I)^2)
	Equilibriumconst(K12,T,H12,C12) ~ NH3_aq/P_NH3 # No activity coefficient for NH3_aq
	Equilibriumconst(K13,T,H13,C13) ~ NH4_aq*OH_aq/NH3_aq # No activity coefficient for NH4+ and NH3_aq
	
	Equilibriumconst(K14,T,H14,C14) ~ H_aq*NO3_aq/P_HNO3*(ActivityCoefficient(T,"H","NO3",I)^2)
	# Equilibriumconst(K15,T,H15,C15) ~ HNO3_aq/P_HNO3 # No activity coefficient for HNO3_aq
	
	# Equilibriumconst(K16,T,H16,C16) ~ H_aq*Cl_aq/P_HCl*ActivityCoefficient(T,"H","Cl",I)^2
	Equilibriumconst(K17,T,H17,C17) ~ HCl_aq/P_HCl # No activity coefficient for HCl_aq
	Equilibriumconst(K18,T,H18,C18) ~ H_aq*OH_aq/a_w
	# Equilibriumconst(K19,T,H19,C19) ~ (Na_aq^2)*SO4_aq*ActivityCoefficient(T,"Na","SO4",I)^3
	Equilibriumconst(K20,T,H20,C20) ~ (NH4_aq^2)*SO4_aq*ActivityCoefficient(T,"NH4","SO4",I)^3
	Equilibriumconst(K21,T,H21,C21) ~ P_NH3*P_HCl

	Equilibriumconst(K22,T,H22,C22) ~ Na_aq*NO3_aq*ActivityCoefficient(T,"Na","NO3",I)^2
	# Equilibriumconst(K23,T,H23,C23) ~ Na_aq*Cl_aq*ActivityCoefficient(T,"Na","Cl",I)^2
	# Equilibriumconst(K24,T,H24,C24) ~ Na_aq*HSO4_aq*ActivityCoefficient(T,"H","HSO4",I)*ActivityCoefficient(T,"Na","Cl",I)/ActivityCoefficient(T,"H","Cl",I)
	# Equilibriumconst(K25,T,H25,C25) ~ P_NH3*P_HNO3
	# Equilibriumconst(K26,T,H26,C26) ~ NH4_aq*HSO4_aq*ActivityCoefficient(T,"H","HSO4",I)*ActivityCoefficient(T,"NH4","Cl",I)/ActivityCoefficient(T,"H","Cl",I)
	# Equilibriumconst(K27,T,H27,C27) ~ (NH4_aq^3)*SO4_aq*HSO4_aq*(ActivityCoefficient(T,"NH4","SO4",I)^3)*(ActivityCoefficient(T,"H","HSO4",I)*ActivityCoefficient(T,"NH4","Cl",I)/ActivityCoefficient(T,"H","Cl",I))^0.5 
	
	# Charge Conservation
	0 ~ Na_aq + H_aq + NH4_aq + K_aq + 2*Ca_aq + 2*Mg_aq - OH_aq - NO3_aq - Cl_aq - HSO4_aq - 2*SO4_aq 
	
	# Water uptake kg/m3
	# W_ ~ Wateruptake([CaNO32_s, CaCl2_s, KHSO4_s, K2SO4_s, KNO3_s, KCl_s, MgSO4_s, MgNO32_s, MgCl2_s, NaNO3_s, NaHSO4_s, NaCl_s, Na2SO4_s, NH42SO4_s, NH4Cl_s, NH4NO3_s, NH4HSO4_s, NH43HSO42_s],a_w)
	0.5 ~ W_
	P_H2O ~ W_/18*T*0.082*10^-5
	
	# Mass conservation
	Na_i ~ Na_aq*W_ + Na2SO4_s + NaHSO4_s + NaNO3_s + NaCl_s # unit: mol/m3
	Ca_i ~ Ca_aq*W_ + CaSO4_s + CaNO32_s + CaCl2_s
	K_i ~ K_aq*W_ + K2SO4_s + KHSO4_s + KNO3_s + KCl_s
	Mg_i ~ Mg_aq*W_ + MgSO4_s + MgNO32_s + MgCl2_s
	H2SO4_i ~ SO4_aq*W_ + HSO4_aq*W_ + NH42SO4_s + Na2SO4_s + CaSO4_s + K2SO4_s + MgSO4_s
	HCl_i ~ P_HCl/(T*0.082*10^-5) + HCl_aq*W_ + Cl_aq*W_ + NH4Cl_s + NaCl_s + CaCl2_s + KCl_s + MgCl2_s
	HNO3_i ~ P_HNO3/(T*0.082*10^-5) + NO3_aq*W_ + HNO3_aq*W_ + NH4NO3_s + NaNO3_s + CaNO32_s + KNO3_s + MgNO32_s
	NH3_i ~ P_NH3/(T*0.082*10^-5) + NH4_aq+W_ + NH3_aq*W_ + NH42SO4_s + NH4HSO4_s + NH4NO3_s + NH4Cl_s
]

@named ns = NonlinearSystem(eqs,[NH4_aq, Na_aq, H_aq, Cl_aq, NO3_aq, SO4_aq, HNO3_aq, NH3_aq, HCl_aq, HSO4_aq, OH_aq, H2O_aq, Ca_aq, K_aq, Mg_aq, P_NH3, P_HNO3, P_HCl, P_H2O, NH42SO4_s, NH4HSO4_s, NH43HSO42_s, NH4NO3_s, NH4Cl_s, NaCl_s, NaNO3_s, NaHSO4_s, Na2SO4_s, CaSO4_s, CaNO32_s, CaCl2_s, K2SO4_s, KHSO4_s, KNO3_s, KCl_s, MgSO4_s, MgNO32_s, MgCl2_s, I, W_],[T, a_w, Na_i, H2SO4_i, NH3_i, HNO3_i, HCl_i, Ca_i, K_i, Mg_i])

nlsys_func = generate_function(ns)[2]
f = eval(nlsys_func)

# [CaNO32_s, CaCl2_s, KHSO4_s, K2SO4_s, KNO3_s, KCl_s, MgSO4_s, MgNO32_s, MgCl2_s, NaNO3_s, NaHSO4_s, NaCl_s, Na2SO4_s, NH42SO4_s, NH4Cl_s, NH4NO3_s, NH4HSO4_s, NH43HSO42_s]
# M_mass(g/mol) = [164.09, 110.98, 136.169, 174.259, 101.1032, 74.5513, 120.366, 148.3, 95.211, 84.9947, 120.06, 58.44, 142.04, 132.14, 53.491, 80.043, 115.11, 247.25] 

M_mass_input = (1, 1, 22.99*10^6, 98.079*10^6, 17.031*10^6, 63.01*10^6, 36.458*10^6, 40.08*10^6, 39.1*10^6, 24.31*10^6)
params = (298.15, 0.55, 0.0, 10.0, 3.4, 2.0, 0.0, 0.4, 0.33, 0.0)./M_mass_input
# T a_w Na_i H2SO4_i NH3_i HNO3_i HCl_i Ca_i K_i Mg_i

j_func = generate_jacobian(ns)[2] # second is in-place
#	takes 20 minutes to run!!!
j! = eval(j_func)

using NLsolve
guess = ones(40)

nlsolve((out, x) -> f(out, x, params), (out, x) -> j!(out, x, params), ones(40))