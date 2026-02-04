export WaterProperties, KelvinEffect, KohlerTheory, DropletGrowth, CloudDynamics,
       IcePhysics, RainFormation, AerosolScavenging, CloudPhysics

"""
    WaterProperties(; name=:WaterProperties)

Thermodynamic properties of water and water solutions, including specific heat,
latent heats, surface tension, and saturation vapor pressure.

Based on Chapter 17 of Seinfeld & Pandis (2006) "Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change", 2nd Edition, John Wiley & Sons.

Equations 17.1-17.6, Table 17.2.
"""
@component function WaterProperties(; name = :WaterProperties)
    @constants begin
        T_ref = 273.15, [description = "Reference temperature (0°C)", unit = u"K"]
        T_unit = 1.0,
        [description = "Unit temperature for nondimensionalization", unit = u"K"]
        # Saturation vapor pressure polynomial coefficients (Table 17.2)
        # Polynomial: p°(mbar) = a0 + a1*T + a2*T^2 + ... (T in °C)
        # Converted to SI: p°(Pa) = 100*(a0 + a1*T_C + a2*T_C^2 + ...)
        # where T_C is temperature difference in K (numerically same as °C difference).
        # Coefficients multiplied by 100 (mbar->Pa conversion factor).
        # After nondimensionalization by T_unit, the effective units are Pa.
        a0_w = 610.7799961, [description = "Water sat. vap. press. coeff. a0 (Table 17.2)", unit = u"Pa"]
        a1_w = 44.36518521, [description = "Water sat. vap. press. coeff. a1 (Table 17.2)", unit = u"Pa"]
        a2_w = 1.428945805, [description = "Water sat. vap. press. coeff. a2 (Table 17.2)", unit = u"Pa"]
        a3_w = 2.650648471e-2,
        [description = "Water sat. vap. press. coeff. a3 (Table 17.2)", unit = u"Pa"]
        a4_w = 3.031240396e-4,
        [description = "Water sat. vap. press. coeff. a4 (Table 17.2)", unit = u"Pa"]
        a5_w = 2.034080948e-6,
        [description = "Water sat. vap. press. coeff. a5 (Table 17.2)", unit = u"Pa"]
        a6_w = 6.136820929e-9,
        [description = "Water sat. vap. press. coeff. a6 (Table 17.2)", unit = u"Pa"]
        # Ice coefficients (also with 100x for mbar->Pa)
        a0_i = 610.9177956, [description = "Ice sat. vap. press. coeff. a0 (Table 17.2)", unit = u"Pa"]
        a1_i = 50.34698970, [description = "Ice sat. vap. press. coeff. a1 (Table 17.2)", unit = u"Pa"]
        a2_i = 1.886013408, [description = "Ice sat. vap. press. coeff. a2 (Table 17.2)", unit = u"Pa"]
        a3_i = 4.176223716e-2,
        [description = "Ice sat. vap. press. coeff. a3 (Table 17.2)", unit = u"Pa"]
        a4_i = 5.824720280e-4,
        [description = "Ice sat. vap. press. coeff. a4 (Table 17.2)", unit = u"Pa"]
        a5_i = 4.838803174e-6,
        [description = "Ice sat. vap. press. coeff. a5 (Table 17.2)", unit = u"Pa"]
        a6_i = 1.838826904e-8,
        [description = "Ice sat. vap. press. coeff. a6 (Table 17.2)", unit = u"Pa"]
        # Eq. 17.2 constants - Specific heat of liquid water
        # c_pw = 4.175 + 1.3e-5*(T-308)^2 + 1.6e-8*(T-308)^4 in J/(g*K)
        # Converted to J/(kg*K) by multiplying by 1000
        cpw_base = 4175.0, [description = "Base specific heat of water (Eq. 17.2)", unit = u"J/(kg*K)"]
        cpw_c2 = 0.013,
        [description = "Quadratic coefficient for c_pw (Eq. 17.2: 1.3e-5 J/(g*K) * 1000)", unit = u"J/(kg*K)"]
        cpw_c4 = 1.6e-5,
        [description = "Quartic coefficient for c_pw (Eq. 17.2: 1.6e-8 J/(g*K) * 1000)", unit = u"J/(kg*K)"]
        T_cpw_ref = 308.0, [description = "Reference T for c_pw correlation (Eq. 17.2)", unit = u"K"]
        # Eq. 17.3 constants - Latent heat of vaporization
        # ΔH_v(kJ/g) = 2.5*(273.15/T)^(0.167 + 3.67e-4*T)
        # Converted to J/kg: multiply by 1e6
        Hv_ref = 2.5e6,
        [description = "Latent heat reference value (Eq. 17.3: 2.5 kJ/g = 2.5e6 J/kg)", unit = u"J/kg"]
        Hv_exp_const = 0.167,
        [description = "Exponent constant for latent heat (Eq. 17.3, dimensionless)", unit = u"1"]
        Hv_exp_coeff = 3.67e-4,
        [description = "Exponent coefficient for latent heat (Eq. 17.3)", unit = u"1/K"]
        # Eq. 17.4 constants - Latent heat of melting
        # ΔH_m(J/g) = 333.5 + 2.03*T_C - 0.0105*T_C^2, T_C in Celsius
        # Converted to J/kg by multiplying by 1000
        Hm_base = 333500.0,
        [description = "Base latent heat of melting (Eq. 17.4: 333.5 J/g = 333500 J/kg)", unit = u"J/kg"]
        Hm_c1 = 2030.0,
        [description = "Linear coeff. for latent heat of melting (Eq. 17.4: 2.03 J/g = 2030 J/kg)",
            unit = u"J/kg"]
        Hm_c2 = -10.5,
        [description = "Quadratic coeff. for latent heat of melting (Eq. 17.4: -0.0105 J/g = -10.5 J/kg)",
            unit = u"J/kg"]
        # Eq. 17.5 constants - Surface tension of pure water
        # σ_w0 = 0.0761 - 1.55e-4*(T-273) in J/m^2 = N/m
        sigma_base = 0.0761, [description = "Base surface tension (Eq. 17.5)", unit = u"N/m"]
        sigma_coeff = 1.55e-4,
        [description = "Surface tension temperature coefficient (Eq. 17.5)", unit = u"N/(m*K)"]
    end

    @parameters begin
        T = 293.15, [description = "Temperature", unit = u"K"]
    end

    @variables begin
        T_C(t), [description = "Temperature difference from 273.15 K", unit = u"K"]
        c_pw(t), [description = "Specific heat of liquid water", unit = u"J/(kg*K)"]
        ΔH_v(t), [description = "Latent heat of vaporization", unit = u"J/kg"]
        ΔH_m(t), [description = "Latent heat of melting", unit = u"J/kg"]
        σ_w0(t), [description = "Surface tension of pure water", unit = u"N/m"]
        p_sat_water(t), [description = "Saturation vapor pressure over water", unit = u"Pa"]
        p_sat_ice(t), [description = "Saturation vapor pressure over ice", unit = u"Pa"]
    end

    # T_C has units of K (it is a temperature difference in kelvin).
    # The polynomial uses dimensionless powers of T_C/T_unit.
    eqs = [
        # Temperature offset from 273.15 K (in kelvin)
        T_C ~ T - T_ref,

        # Eq. 17.2 - Specific heat of liquid water (273 < T <= 308 K)
        c_pw ~
        cpw_base + cpw_c2 * ((T - T_cpw_ref) / T_unit)^2 +
        cpw_c4 * ((T - T_cpw_ref) / T_unit)^4,

        # Eq. 17.3 - Latent heat of vaporization (J/kg)
        # Hv_exp_coeff has units 1/K, so Hv_exp_coeff * T is dimensionless
        ΔH_v ~ Hv_ref * (T_ref / T)^(Hv_exp_const + Hv_exp_coeff * T),

        # Eq. 17.4 - Latent heat of melting (J/g)
        # T_C/T_unit makes the polynomial argument dimensionless
        ΔH_m ~ Hm_base + Hm_c1 * (T_C / T_unit) + Hm_c2 * (T_C / T_unit)^2,

        # Eq. 17.5 - Surface tension of pure water (N/m)
        σ_w0 ~ sigma_base - sigma_coeff * T_C,

        # Table 17.2 - Saturation vapor pressure over water (Pa)
        # Polynomial in dimensionless T_C/T_unit
        p_sat_water ~
        a0_w + a1_w * (T_C / T_unit) + a2_w * (T_C / T_unit)^2 +
        a3_w * (T_C / T_unit)^3 + a4_w * (T_C / T_unit)^4 +
        a5_w * (T_C / T_unit)^5 + a6_w * (T_C / T_unit)^6,

        # Table 17.2 - Saturation vapor pressure over ice (Pa)
        p_sat_ice ~
        a0_i + a1_i * (T_C / T_unit) + a2_i * (T_C / T_unit)^2 +
        a3_i * (T_C / T_unit)^3 + a4_i * (T_C / T_unit)^4 +
        a5_i * (T_C / T_unit)^5 + a6_i * (T_C / T_unit)^6
    ]

    return System(eqs, t; name)
end

"""
    KelvinEffect(; name=:KelvinEffect)

Kelvin (curvature) effect on vapor pressure over curved droplet surfaces.
The vapor pressure over a curved surface is enhanced relative to a flat surface
due to surface energy effects.

Based on Equations 17.9 and 17.28 from Seinfeld & Pandis (2006).

# Equations
- Eq. 17.9: p_w(D_p)/p° = exp(4M_w σ_w / (R T ρ_w D_p))
- Eq. 17.28: A = 4M_w σ_w / (R T ρ_w) ≈ 0.66/T µm (for T in K)
"""
@component function KelvinEffect(; name = :KelvinEffect)
    @constants begin
        M_w = 0.018015, [description = "Molecular weight of water", unit = u"kg/mol"]
        R = 8.314, [description = "Universal gas constant", unit = u"J/(mol*K)"]
        ρ_w = 997.0, [description = "Density of liquid water at 25°C", unit = u"kg/m^3"]
    end

    @parameters begin
        T = 293.15, [description = "Temperature", unit = u"K"]
        σ_w = 0.0728, [description = "Surface tension of water (Table 17.1)", unit = u"N/m"]
        D_p = 1.0e-7, [description = "Droplet diameter", unit = u"m"]
    end

    @variables begin
        kelvin_ratio(t),
        [description = "Kelvin ratio p_w(D_p)/p° (Eq. 17.9, dimensionless)", unit = u"1"]
        A_kelvin(t), [description = "Kelvin parameter A (Eq. 17.28)", unit = u"m"]
    end

    eqs = [
        # Eq. 17.28 - Kelvin parameter A = 4*M_w*σ_w / (R*T*ρ_w)
        A_kelvin ~ 4 * M_w * σ_w / (R * T * ρ_w),

        # Eq. 17.9 - Kelvin equation: p_w(D_p)/p° = exp(A/D_p)
        kelvin_ratio ~ exp(A_kelvin / D_p)
    ]

    return System(eqs, t; name)
end

"""
    KohlerTheory(; name=:KohlerTheory)

Köhler theory for equilibrium vapor pressure over solution droplets,
combining the Kelvin (curvature) effect and Raoult (solute) effect.

Includes standard Köhler equation, critical diameter and supersaturation,
and extension for partially insoluble particles.

Based on Equations 17.26-17.40 from Seinfeld & Pandis (2006).

# Key Equations
- Eq. 17.27: ln(p_w(D_p)/p°) = A/D_p - B/D_p³
- Eq. 17.28: A = 4M_w σ_w / (R T ρ_w)
- Eq. 17.29: B = 6n_s M_w / (π ρ_w)
- Eq. 17.30: D_pc = √(3B/A) (critical diameter)
- Eq. 17.32: ln(S_c) = √(4A³/27B) (critical saturation)
- Eq. 17.33: n_s = νπd_s³ρ_s / (6M_s) (moles of solute)
- Eq. 17.38-17.40: Extensions for insoluble material
"""
@component function KohlerTheory(; name = :KohlerTheory)
    @constants begin
        M_w = 0.018015, [description = "Molecular weight of water", unit = u"kg/mol"]
        R = 8.314, [description = "Universal gas constant", unit = u"J/(mol*K)"]
        ρ_w = 997.0, [description = "Density of liquid water", unit = u"kg/m^3"]
        π_val = 3.141592653589793, [description = "Pi (dimensionless)", unit = u"1"]
        vol_eps = 1e-60,
        [description = "Small volume for numerical stability", unit = u"m^3"]
    end

    @parameters begin
        T = 293.15, [description = "Temperature", unit = u"K"]
        σ_w = 0.0728, [description = "Surface tension of water", unit = u"N/m"]
        D_p = 1.0e-7, [description = "Wet droplet diameter", unit = u"m"]
        d_s = 5.0e-8, [description = "Dry particle diameter", unit = u"m"]
        ν_s = 2.0,
        [description = "Van't Hoff factor (number of ions) (dimensionless)", unit = u"1"]
        M_s = 0.05844, [description = "Molecular weight of solute", unit = u"kg/mol"]
        ρ_s = 2165.0, [description = "Density of solute", unit = u"kg/m^3"]
        ε_m = 1.0,
        [description = "Mass fraction of soluble material (dimensionless)", unit = u"1"]
        ρ_u = 2000.0, [description = "Density of insoluble material", unit = u"kg/m^3"]
    end

    @variables begin
        A(t), [description = "Kelvin parameter", unit = u"m"]
        B(t), [description = "Solute parameter", unit = u"m^3"]
        n_s(t), [description = "Moles of solute", unit = u"mol"]
        d_u(t), [description = "Diameter of insoluble core", unit = u"m"]
        ln_S(t), [description = "Log of saturation ratio (dimensionless)", unit = u"1"]
        S(t), [description = "Saturation ratio (dimensionless)", unit = u"1"]
        D_pc(t), [description = "Critical droplet diameter", unit = u"m"]
        S_c(t), [description = "Critical saturation ratio (dimensionless)", unit = u"1"]
        ln_S_c(t), [description = "Log of critical saturation (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Eq. 17.28 - Kelvin parameter A (m)
        A ~ 4 * M_w * σ_w / (R * T * ρ_w),

        # Eq. 17.40 - Moles of solute for mixed particle
        n_s ~ (ε_m / M_s) * ν_s * π_val * d_s^3 / (6 * (ε_m / ρ_s + (1 - ε_m) / ρ_u)),

        # Eq. 17.29 - Solute parameter B (m^3)
        B ~ 6 * n_s * M_w / (π_val * ρ_w),

        # Eq. 17.39 - Insoluble core diameter (m)
        d_u ~ d_s * (1 / ((ρ_u / ρ_s) * (ε_m / (1 - ε_m + 1e-30)) + 1))^(1 / 3),

        # Eq. 17.38 - Köhler equation with insoluble core (log saturation)
        ln_S ~ A / D_p - B / (D_p^3 - d_u^3 + vol_eps),

        # Saturation ratio S
        S ~ exp(ln_S),

        # Eq. 17.30 - Critical droplet diameter
        D_pc ~ sqrt(3 * B / A),

        # Eq. 17.32 - Critical saturation
        ln_S_c ~ sqrt(4 * A^3 / (27 * B)),

        # Critical saturation ratio
        S_c ~ exp(ln_S_c)
    ]

    return System(eqs, t; name)
end

"""
    DropletGrowth(; name=:DropletGrowth)

Diffusional growth rate of cloud droplets by condensation of water vapor.
Accounts for the thermal gradient between droplet surface and environment,
gas-phase diffusion, and modified diffusivity for kinetic effects.

Based on Equations 17.60-17.72 from Seinfeld & Pandis (2006).

# Key Equations
- Eq. 17.61: D_v = 0.211(p_ref/p)(T/273)^1.94 cm²/s (water vapor diffusivity)
- Eq. 17.62: D'_v = D_v / [1 + 2D_v/(α_c D_p)√(2πM_w/RT)] (kinetic correction)
- Eq. 17.70: D_p dD_p/dt = [S_v,∞ - S_eq] / [ρ_w RT/(4p° D'_v M_w) + ΔH_v ρ_w/(4k'_a T)(ΔH_v M_w/(TR) - 1)]
- Eq. 17.71: k_a = 10⁻³(4.39 + 0.071T) J/(m·s·K) (thermal conductivity)
- Eq. 17.72: k'_a correction for kinetic effects
"""
@component function DropletGrowth(; name = :DropletGrowth)
    @constants begin
        M_w = 0.018015, [description = "Molecular weight of water", unit = u"kg/mol"]
        M_a = 0.02897, [description = "Molecular weight of air", unit = u"kg/mol"]
        R = 8.314, [description = "Universal gas constant", unit = u"J/(mol*K)"]
        ρ_w = 997.0, [description = "Density of liquid water", unit = u"kg/m^3"]
        π_val = 3.141592653589793, [description = "Pi (dimensionless)", unit = u"1"]
        # Eq. 17.61 constants (D_v in m^2/s)
        Dv_ref = 0.211e-4,
        [description = "Reference diffusivity at 273 K, 1 atm", unit = u"m^2/s"]
        p_ref = 101325.0, [description = "Reference pressure (1 atm)", unit = u"Pa"]
        T_Dv_ref = 273.0, [description = "Reference temperature for D_v", unit = u"K"]
        # Eq. 17.71 constants (k_a in W/(m*K))
        ka_c0 = 4.39e-3,
        [description = "Thermal conductivity constant", unit = u"J/(m*s*K)"]
        ka_c1 = 0.071e-3,
        [description = "Thermal conductivity slope", unit = u"J/(m*s*K^2)"]
        # Eq. 17.72 constant (air density at STP)
        ρ_air_stp = 1.293, [description = "Air density at STP", unit = u"kg/m^3"]
        # Saturation vapor pressure (Magnus formula)
        p_sat_ref = 610.94, [description = "Reference saturation pressure", unit = u"Pa"]
        magnus_a = 17.625,
        [description = "Magnus coefficient a (dimensionless)", unit = u"1"]
        T_magnus_offset = 273.15, [description = "Magnus temperature offset", unit = u"K"]
        magnus_b = 243.04, [description = "Magnus coefficient b", unit = u"K"]
        T_unit = 1.0,
        [description = "Unit temperature for nondimensionalization", unit = u"K"]
        # Latent heat constants (Eq. 17.3, converted to J/kg)
        Hv_ref_kg = 2.5e6,
        [description = "Reference latent heat of vaporization", unit = u"J/kg"]
        T_Hv_ref = 273.15,
        [description = "Reference temperature for latent heat", unit = u"K"]
        Hv_exp_const = 0.167,
        [description = "Exponent constant (dimensionless)", unit = u"1"]
        Hv_exp_coeff = 3.67e-4, [description = "Exponent coefficient", unit = u"1/K"]
        one = 1.0, [description = "Unit value (dimensionless)", unit = u"1"]
        vol_eps = 1e-60,
        [description = "Small volume for numerical stability", unit = u"m^3"]
    end

    @parameters begin
        T = 293.15, [description = "Temperature", unit = u"K"]
        p_a = 101325.0, [description = "Atmospheric pressure", unit = u"Pa"]
        S_v = 1.01,
        [description = "Environmental saturation ratio (dimensionless)", unit = u"1"]
        σ_w = 0.0728, [description = "Surface tension of water", unit = u"N/m"]
        D_p = 1.0e-6, [description = "Droplet diameter", unit = u"m"]
        n_s = 0.0, [description = "Moles of solute in droplet", unit = u"mol"]
        d_u = 0.0, [description = "Insoluble core diameter", unit = u"m"]
        α_c = 1.0,
        [description = "Mass accommodation coefficient (dimensionless)", unit = u"1"]
        α_T = 1.0,
        [description = "Thermal accommodation coefficient (dimensionless)", unit = u"1"]
        c_p_air = 1005.0, [description = "Specific heat of air", unit = u"J/(kg*K)"]
    end

    @variables begin
        D_v(t), [description = "Diffusivity of water vapor in air", unit = u"m^2/s"]
        D_v_prime(t), [description = "Modified diffusivity of water vapor", unit = u"m^2/s"]
        k_a(t), [description = "Thermal conductivity of air", unit = u"J/(m*s*K)"]
        k_a_prime(t),
        [description = "Modified thermal conductivity of air", unit = u"J/(m*s*K)"]
        p_sat(t), [description = "Saturation vapor pressure", unit = u"Pa"]
        ΔH_v(t), [description = "Latent heat of vaporization", unit = u"J/kg"]
        S_eq(t),
        [
            description = "Equilibrium saturation ratio over droplet (dimensionless)", unit = u"1"]
        dDp_dt(t), [description = "Droplet diameter growth rate", unit = u"m/s"]
    end

    eqs = [
        # Eq. 17.61 - Diffusivity of water vapor in air (m^2/s)
        D_v ~ Dv_ref * (p_ref / p_a) * (T / T_Dv_ref)^1.94,

        # Eq. 17.62 - Modified diffusivity (kinetic correction)
        D_v_prime ~ D_v / (one + (2 * D_v / (α_c * D_p)) * sqrt(2 * π_val * M_w / (R * T))),

        # Eq. 17.71 - Thermal conductivity of air (J/(m*s*K))
        k_a ~ ka_c0 + ka_c1 * T,

        # Eq. 17.72 - Modified thermal conductivity (kinetic correction)
        k_a_prime ~
        k_a / (one +
         (2 * k_a / (α_T * D_p * ρ_air_stp * c_p_air)) * sqrt(2 * π_val * M_a / (R * T))),

        # Saturation vapor pressure (Magnus formula)
        p_sat ~
        p_sat_ref *
        exp(magnus_a * (T - T_magnus_offset) / ((T - T_magnus_offset) + magnus_b)),

        # Eq. 17.3 - Latent heat of vaporization (J/kg)
        ΔH_v ~ Hv_ref_kg * (T_Hv_ref / T)^(Hv_exp_const + Hv_exp_coeff * T),

        # Equilibrium saturation ratio over the droplet (Köhler, Eq. 17.27)
        S_eq ~ exp(4 * M_w * σ_w / (R * T * ρ_w * D_p) -
            6 * n_s * M_w / (π_val * ρ_w * (D_p^3 - d_u^3 + vol_eps))),

        # Eq. 17.70 - Simplified droplet growth equation
        dDp_dt ~
        (one / D_p) * (S_v - S_eq) /
        (ρ_w * R * T / (4 * p_sat * D_v_prime * M_w) +
         ΔH_v * ρ_w / (4 * k_a_prime * T) * (ΔH_v * M_w / (T * R) - one))
    ]

    return System(eqs, t; name)
end

"""
    CloudDynamics(; name=:CloudDynamics)

Adiabatic cooling, cloud formation, and supersaturation dynamics.
Includes dry adiabatic lapse rate, lifting condensation level (LCL),
dew point temperature, and saturation mixing ratio.

Based on Equations 17.43-17.49 from Seinfeld & Pandis (2006).

# Key Equations
- Eq. 17.43: T_d ≈ T_0 + (RT_0²/ΔH_v M_w)ln(RH) (dew point temperature)
- Eq. 17.44: Γ_d = g/c_p ≈ 9.76 K/km (dry adiabatic lapse rate)
- Eq. 17.46: w_v = M_w p_w / (M_a p_a) (water vapor mixing ratio)
- Eq. 17.49: h_LCL = (T_0 - T_L)/Γ (lifting condensation level)
"""
@component function CloudDynamics(; name = :CloudDynamics)
    @constants begin
        M_w = 0.018015, [description = "Molecular weight of water", unit = u"kg/mol"]
        M_a = 0.02897, [description = "Molecular weight of air", unit = u"kg/mol"]
        R = 8.314, [description = "Universal gas constant", unit = u"J/(mol*K)"]
        g = 9.81, [description = "Gravitational acceleration", unit = u"m/s^2"]
        c_p_air = 1005.0,
        [description = "Specific heat of air at constant pressure", unit = u"J/(kg*K)"]
        # Magnus formula constants
        p_sat_ref = 610.94, [description = "Reference saturation pressure", unit = u"Pa"]
        magnus_a = 17.625,
        [description = "Magnus coefficient a (dimensionless)", unit = u"1"]
        T_magnus_offset = 273.15, [description = "Magnus temperature offset", unit = u"K"]
        magnus_b = 243.04, [description = "Magnus coefficient b", unit = u"K"]
    end

    @parameters begin
        T_0 = 293.15, [description = "Initial temperature", unit = u"K"]
        RH = 0.8, [description = "Relative humidity (dimensionless)", unit = u"1"]
        p_a = 101325.0, [description = "Atmospheric pressure", unit = u"Pa"]
        ΔH_v = 2.5e6, [description = "Latent heat of vaporization", unit = u"J/kg"]
    end

    @variables begin
        Γ_d(t), [description = "Dry adiabatic lapse rate", unit = u"K/m"]
        T_d(t), [description = "Dew point temperature", unit = u"K"]
        h_LCL(t), [description = "Lifting condensation level height", unit = u"m"]
        p_sat(t), [description = "Saturation vapor pressure at T_0", unit = u"Pa"]
        w_vs(t), [description = "Saturation mixing ratio (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Eq. 17.44 - Dry adiabatic lapse rate
        Γ_d ~ g / c_p_air,

        # Saturation vapor pressure at T_0 (Magnus formula)
        p_sat ~
        p_sat_ref *
        exp(magnus_a * (T_0 - T_magnus_offset) / ((T_0 - T_magnus_offset) + magnus_b)),

        # Saturation mixing ratio
        w_vs ~ M_w * p_sat / (M_a * p_a),

        # Eq. 17.43 - Dew point temperature
        T_d ~ T_0 + (R * T_0^2 / (ΔH_v * M_w)) * log(RH),

        # Eq. 17.49 - Lifting condensation level height
        h_LCL ~ (T_0 - T_d) / Γ_d
    ]

    return System(eqs, t; name)
end

"""
    IcePhysics(; name=:IcePhysics)

Ice nucleation and ice-water equilibrium thermodynamics.
Includes freezing point depression for solutions, Kelvin equation for ice particles,
and empirical ice nuclei concentration.

Based on Equations 17.100-17.103 from Seinfeld & Pandis (2006).

# Key Equations
- Eq. 17.100: T_0 - T_e = (RT_0 T_e/ΔH_m)(-ln a_w) (freezing point from water activity)
- Eq. 17.101: ΔT_f = (RT_0² M_w)/(1000 ΔH_m) m ≈ 1.86 K·kg/mol × m (freezing point depression)
- Eq. 17.102: p_i = p_sat,i exp(4M_w σ_ia / (R T ρ_i D_p)) (Kelvin effect for ice)
- Eq. 17.103: IN(L⁻¹) = exp[0.6(253-T)] (ice nuclei concentration)
"""
@component function IcePhysics(; name = :IcePhysics)
    @constants begin
        M_w = 0.018015, [description = "Molecular weight of water", unit = u"kg/mol"]
        R = 8.314, [description = "Universal gas constant", unit = u"J/(mol*K)"]
        ρ_i = 917.0, [description = "Density of ice", unit = u"kg/m^3"]
        T_freeze = 273.15, [description = "Freezing point of pure water", unit = u"K"]
        ΔH_m_mass = 333500.0,
        [description = "Latent heat of fusion of water", unit = u"J/kg"]
        ΔH_s = 2.83e6, [description = "Latent heat of sublimation", unit = u"J/kg"]
        # Eq. 17.101: ΔT_f = R*T_0^2*M_w / (1000*ΔH_m) * m
        # 1000 g/kg factor: convert ΔH_m from J/kg to J/g effectively
        # The formula is ΔT_f = (R*T_0^2*M_w)/(ΔH_m_molar) * m with ΔH_m_molar = ΔH_m*M_w
        # Actually: ΔT_f = K_f * m, where K_f = R*T_0^2 / (1000*ΔH_m_g) = R*T_0^2 / (ΔH_m_mass*M_w*(1000/M_w))
        # For simplicity, use the cryoscopic constant for water = 1.86 K·kg/mol
        K_f = 1.86, [description = "Cryoscopic constant of water", unit = u"K*kg/mol"]
        # Eq. 17.103 constants
        IN_ref = 1000.0,
        [
            description = "Ice nuclei concentration reference (per L -> per m^3)", unit = u"1/m^3"]
        IN_coeff = 0.6, [description = "Ice nuclei exponential coefficient", unit = u"1/K"]
        T_IN_ref = 253.0,
        [description = "Reference temperature for IN parameterization", unit = u"K"]
    end

    @parameters begin
        T = 263.15, [description = "Temperature", unit = u"K"]
        σ_ia = 0.106, [description = "Ice-air surface tension", unit = u"N/m"]
        D_p = 1.0e-6, [description = "Ice particle diameter", unit = u"m"]
        a_w = 1.0, [description = "Water activity (dimensionless)", unit = u"1"]
        molality = 0.0, [description = "Solution molality", unit = u"mol/kg"]
    end

    @variables begin
        ΔT_f(t), [description = "Freezing point depression", unit = u"K"]
        T_e(t), [description = "Equilibrium freezing temperature", unit = u"K"]
        kelvin_ice(t), [description = "Kelvin factor for ice (dimensionless)", unit = u"1"]
        IN(t), [description = "Ice nuclei concentration", unit = u"1/m^3"]
    end

    eqs = [
        # Eq. 17.101 - Freezing point depression
        ΔT_f ~ K_f * molality,

        # Eq. 17.100 - Equilibrium freezing temperature from water activity
        # Uses approximation T_0*T_e ≈ T_0^2, where ΔH_m_molar = ΔH_m_mass * M_w
        T_e ~ T_freeze + R * T_freeze^2 / (ΔH_m_mass * M_w) * log(a_w + 1e-30),

        # Eq. 17.102 - Kelvin equation for ice particles
        kelvin_ice ~ exp(4 * M_w * σ_ia / (R * T * ρ_i * D_p)),

        # Eq. 17.103 - Ice nuclei concentration (per liter, converted to per m^3)
        IN ~ IN_ref * exp(IN_coeff * (T_IN_ref - T))
    ]

    return System(eqs, t; name)
end

"""
    RainFormation(; name=:RainFormation)

Rain formation by collision and coalescence processes, plus raindrop
size distributions.

Based on Equations 17.104-17.108 from Seinfeld & Pandis (2006).

# Key Equations
- Eq. 17.104: E = y²/(D_p + d_p)² (collision efficiency)
- Eq. 17.105: E_c = (D_p/(D_p + d_p))² (coalescence efficiency, for D_p > 400 µm)
- Eq. 17.106: dm/dt = (π/4)E_t(D_p + d_p)²w_L(v_D - v_d) (continuous accretion)
- Eq. 17.107: F(D_p) = 1 - exp[-(D_p/(1.3p_0^0.232 mm))^2.25] (Best CDF)
- Eq. 17.108: n(D_p) = n_0 exp(-ψD_p) (Marshall-Palmer distribution)
"""
@component function RainFormation(; name = :RainFormation)
    @constants begin
        π_val = 3.141592653589793, [description = "Pi (dimensionless)", unit = u"1"]
        # Eq. 17.107 constants - Best (1950) raindrop distribution
        # F(D_p) = 1 - exp[-(D_p / (1.3*p_0^0.232 mm))^2.25]
        # where p_0 is rainfall intensity in mm/h
        F_diam_ref = 1.3e-3,
        [description = "Best distribution diameter scale (Eq. 17.107: 1.3 mm)", unit = u"m"]
        F_exp = 2.25,
        [description = "Best distribution exponent (Eq. 17.107, dimensionless)", unit = u"1"]
        one = 1.0, [description = "Unity (dimensionless)", unit = u"1"]
        # For p_0 exponent nondimensionalization
        p0_unit = 1.0,
        [description = "Unit rainfall intensity for nondimensionalization", unit = u"mm/hr"]
    end

    @parameters begin
        D_p = 1.0e-3, [description = "Collector drop diameter", unit = u"m"]
        d_p = 1.0e-5, [description = "Collected droplet diameter", unit = u"m"]
        E_coll = 0.5, [description = "Collision efficiency (Eq. 17.104, dimensionless)", unit = u"1"]
        E_coal = 1.0, [description = "Coalescence efficiency (Eq. 17.105, dimensionless)", unit = u"1"]
        w_L = 1.0e-3,
        [description = "Liquid water content (mass mixing ratio, dimensionless)", unit = u"1"]
        v_D = 5.0, [description = "Fall velocity of collector drop", unit = u"m/s"]
        v_d = 0.01, [description = "Fall velocity of collected droplet", unit = u"m/s"]
        ρ_a = 1.225, [description = "Air density", unit = u"kg/m^3"]
        p_0 = 1.0, [description = "Rainfall intensity (Eq. 17.107)", unit = u"mm/hr"]
        n_0 = 8.0e6,
        [description = "Marshall-Palmer intercept parameter (Eq. 17.108: 8000 m^-3 mm^-1)", unit = u"1/m^4"]
        ψ = 4100.0, [description = "Marshall-Palmer slope parameter (Eq. 17.108)", unit = u"1/m"]
    end

    @variables begin
        E_t(t), [description = "Collection efficiency (Eq. 17.106, dimensionless)", unit = u"1"]
        dm_dt(t), [description = "Mass accretion rate (Eq. 17.106)", unit = u"kg/s"]
        F_dist(t), [description = "Best raindrop CDF (Eq. 17.107, dimensionless)", unit = u"1"]
        n_MP(t), [description = "Marshall-Palmer number density (Eq. 17.108)", unit = u"1/m^4"]
    end

    eqs = [
        # Collection efficiency (product of collision and coalescence)
        E_t ~ E_coll * E_coal,

        # Eq. 17.106 - Continuous accretion equation
        # dm/dt = (π/4) * E_t * (D_p + d_p)^2 * w_L * ρ_a * |v_D - v_d|
        dm_dt ~ (π_val / 4) * E_t * (D_p + d_p)^2 * w_L * ρ_a * abs(v_D - v_d),

        # Eq. 17.107 - Best raindrop cumulative distribution
        # F(D_p) = 1 - exp[-(D_p / (F_diam_ref * p_0^0.232))^F_exp]
        # p_0/p0_unit makes the base of the power dimensionless
        F_dist ~ one - exp(-(D_p / (F_diam_ref * (p_0 / p0_unit)^0.232))^F_exp),

        # Eq. 17.108 - Marshall-Palmer raindrop size distribution
        # n(D_p) = n_0 * exp(-ψ * D_p)
        n_MP ~ n_0 * exp(-ψ * D_p)
    ]

    return System(eqs, t; name)
end

"""
    AerosolScavenging(; name=:AerosolScavenging)

Aerosol scavenging by cloud droplets, including mass and number
scavenging ratios and scavenging coefficients.

Based on Equations 17.84-17.92 from Seinfeld & Pandis (2006).

# Key Equations
- Eq. 17.84: F_i = (C_{i,0} - C_{i,int})/C_{i,0} (mass scavenging ratio)
- Eq. 17.85: F_N = (N_0 - N_int)/N_0 (number scavenging ratio)
- Eq. 17.91: Λ(D_p,t) = ∫K(D_p,x)n_d(x,t)dx (scavenging coefficient)
- Eq. 17.92: Λ(D_p) = N_d K(D_p, 10µm) (monodisperse approximation)
"""
@component function AerosolScavenging(; name = :AerosolScavenging)
    @parameters begin
        C_0 = 1.0e9, [description = "Initial aerosol mass concentration", unit = u"1/m^3"]
        C_int = 5.0e8,
        [description = "Interstitial aerosol mass concentration", unit = u"1/m^3"]
        N_0 = 1.0e9, [description = "Initial aerosol number concentration", unit = u"1/m^3"]
        N_int = 3.0e8,
        [description = "Interstitial aerosol number concentration", unit = u"1/m^3"]
        K_coll = 1.0e-12, [description = "Collision coefficient", unit = u"m^3/s"]
        N_d = 1.0e8, [description = "Cloud droplet number concentration", unit = u"1/m^3"]
    end

    @variables begin
        F_mass(t), [description = "Mass scavenging ratio (dimensionless)", unit = u"1"]
        F_number(t), [description = "Number scavenging ratio (dimensionless)", unit = u"1"]
        Λ(t), [description = "Scavenging coefficient", unit = u"1/s"]
    end

    eqs = [
        # Eq. 17.84 - Mass scavenging ratio
        F_mass ~ (C_0 - C_int) / C_0,

        # Eq. 17.85 - Number scavenging ratio
        F_number ~ (N_0 - N_int) / N_0,

        # Eq. 17.92 - Scavenging coefficient for monodisperse cloud droplets
        Λ ~ N_d * K_coll
    ]

    return System(eqs, t; name)
end

"""
    CloudPhysics(; name=:CloudPhysics)

Complete cloud physics system composing all subsystems: water properties,
Kelvin effect, Köhler theory, droplet growth, cloud dynamics, ice physics,
rain formation, and aerosol scavenging.

Based on Chapter 17 of Seinfeld & Pandis (2006) "Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change", 2nd Edition, John Wiley & Sons.
"""
@component function CloudPhysics(; name = :CloudPhysics)
    water = WaterProperties(; name = :water)
    kelvin = KelvinEffect(; name = :kelvin)
    kohler = KohlerTheory(; name = :kohler)
    growth = DropletGrowth(; name = :growth)
    cloud = CloudDynamics(; name = :cloud)
    ice = IcePhysics(; name = :ice)
    rain = RainFormation(; name = :rain)
    scav = AerosolScavenging(; name = :scav)

    return System(
        Equation[], t; systems = [
            water, kelvin, kohler, growth, cloud, ice, rain, scav], name)
end
