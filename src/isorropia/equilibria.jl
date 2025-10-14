@mtkmodel EqConst begin
    @description """An equilibrium constant based on Equation 5 in Fountoukis and Nenes (2007),
parameterized by values from Table 2 in the same paper."""
    @constants begin
        logK⁰, [description = "Log of the equilibrium constant at 298.15 K"]
        H⁺, [description = "ΔH⁰ / (R * T₀) (unitless)"]
        C⁺, [description = "ΔC⁰ₚ / R (unitless)"]
        T₀ = 293.15, [unit = u"K", description = "Standard temperature"]
    end
    @parameters begin
        T, [unit = u"K", description = "Temperature"]
    end
    @variables begin
        logK_eq(t), [description = "Log of the equilibrium constant", guess=logK⁰]
    end
    @equations begin
        logK_eq ~ logK⁰ + (-H⁺ * (T₀ / T - 1) - C⁺ * (1 + log(T₀ / T) - T₀ / T))
    end
end

@mtkmodel EquilibriumConstants begin
    @parameters begin
        T, [unit = u"K", description = "Temperature"]
    end
    @components begin
        # Equilibrium constants from Table 2 of Fountoukis and Nenes (2007).
        # NOTE: Assuming that H⁺ and C⁺ are zero when they are left out of Table 2.
        r1 = EqConst(logK⁰ = log(6.067e5), H⁺ = -11.299, C⁺ = 0.0)
        r2 = EqConst(logK⁰ = log(7.974e11), H⁺ = -14.087, C⁺ = 0.0)
        r3 = EqConst(logK⁰ = log(4.319e-5), H⁺ = 0.0, C⁺ = 0.0)
        r4 = EqConst(logK⁰ = log(1.569e-2), H⁺ = -9.589, C⁺ = 45.807)
        r5 = EqConst(logK⁰ = log(24.016), H⁺ = -8.423, C⁺ = 17.964)
        r6 = EqConst(logK⁰ = log(0.872), H⁺ = 14.075, C⁺ = 19.388)
        r7 = EqConst(logK⁰ = log(8.680), H⁺ = -6.167, C⁺ = 19.953)
        r8 = EqConst(logK⁰ = log(1.079e5), H⁺ = 36.798, C⁺ = 0.0)
        r9 = EqConst(logK⁰ = log(2.507e15), H⁺ = -8.754, C⁺ = 0.0)
        r10 = EqConst(logK⁰ = log(9.557e21), H⁺ = -1.347, C⁺ = 0.0)
        r11 = EqConst(logK⁰ = log(1.015e-2), H⁺ = 8.85, C⁺ = 25.14)
        r12 = EqConst(logK⁰ = log(5.764e1), H⁺ = 13.79, C⁺ = -5.39)
        r13 = EqConst(logK⁰ = log(1.805e-5), H⁺ = -1.50, C⁺ = 26.92)
        r14 = EqConst(logK⁰ = log(2.511e6), H⁺ = 29.17, C⁺ = 16.83)
        r15 = EqConst(logK⁰ = log(2.1e5), H⁺ = 29.17, C⁺ = 16.83)
        r16 = EqConst(logK⁰ = log(1.971e6), H⁺ = 30.20, C⁺ = 19.91)
        r17 = EqConst(logK⁰ = log(2.5e3), H⁺ = 30.20, C⁺ = 19.91)
        r18 = EqConst(logK⁰ = log(1.010e-14), H⁺ = -22.52, C⁺ = 26.92)
        r19 = EqConst(logK⁰ = log(4.799e-1), H⁺ = 0.98, C⁺ = 39.75)
        r20 = EqConst(logK⁰ = log(1.87e0), H⁺ = -2.65, C⁺ = 38.57)
        r21 = EqConst(logK⁰ = log(1.086e-16), H⁺ = -71.00, C⁺ = 2.40)
        r22 = EqConst(logK⁰ = log(1.197e1), H⁺ = -8.22, C⁺ = 16.01)
        r23 = EqConst(logK⁰ = log(3.766e1), H⁺ = -1.56, C⁺ = 16.90)
        r24 = EqConst(logK⁰ = log(2.413e4), H⁺ = 0.79, C⁺ = 14.75)
        r25 = EqConst(logK⁰ = log(4.199e-17), H⁺ = -74.375, C⁺ = 6.025)
        r26 = EqConst(logK⁰ = log(1.383e0), H⁺ = -2.87, C⁺ = 15.83)
        r27 = EqConst(logK⁰ = log(2.972e1), H⁺ = -5.19, C⁺ = 54.40)
    end
    @constants begin
        k1_unit = 1, [unit = u"mol^3/kg^3", description = "Unit conversion factor"]
        k2_unit = 1, [unit = u"mol^3/kg^3", description = "Unit conversion factor"]
        k3_unit = 1, [unit = u"mol^2/kg^2", description = "Unit conversion factor"]
        k4_unit = 1, [unit = u"mol^3/kg^3", description = "Unit conversion factor"]
        k5_unit = 1, [unit = u"mol^2/kg^2", description = "Unit conversion factor"]
        k6_unit = 1, [unit = u"mol^2/kg^2", description = "Unit conversion factor"]
        k7_unit = 1, [unit = u"mol^2/kg^2", description = "Unit conversion factor"]
        k8_unit = 1, [unit = u"mol^2/kg^2", description = "Unit conversion factor"]
        k9_unit = 1, [unit = u"mol^3/kg^3", description = "Unit conversion factor"]
        k10_unit = 1, [unit = u"mol^3/kg^3", description = "Unit conversion factor"]
        k11_unit = 1, [unit = u"mol/kg", description = "Unit conversion factor"]
        k12_unit = 1, [unit = u"mol/kg/Constants.atm", description = "Unit conversion"]
        k13_unit = 1, [unit = u"mol/kg", description = "Unit conversion factor"]
        k14_unit = 1, [unit = u"mol^2/kg^2/Constants.atm", description = "Unit conversion"]
        k1b_unit = 1, [unit = u"mol/kg", description = "Unit conversion"]
        k15_unit = 1, [unit = u"mol^2/kg^2/Constants.atm", description = "Unit conversion"]
        k16_unit = 1, [unit = u"mol^2/kg^2/Constants.atm", description = "Unit conversion"]
        k17_unit = 1, [unit = u"mol/kg/Constants.atm", description = "Unit conversion"]
        k2b_unit = 1, [unit = u"mol/kg", description = "Unit conversion factor"]
        k18_unit = 1, [unit = u"mol^2/kg^2", description = "Unit conversion factor"]
        k19_unit = 1, [unit = u"mol^3/kg^3", description = "Unit conversion factor"]
        k20_unit = 1, [unit = u"mol^3/kg^3", description = "Unit conversion factor"]
        k21_unit = 1, [unit = u"Constants.atm^2", description = "Unit conversion factor"]
        k22_unit = 1, [unit = u"mol^2/kg^2", description = "Unit conversion factor"]
        k23_unit = 1, [unit = u"mol^2/kg^2", description = "Unit conversion factor"]
        k24_unit = 1, [unit = u"mol^2/kg^2", description = "Unit conversion factor"]
        k25_unit = 1, [unit = u"Constants.atm^2", description = "Unit conversion factor"]
        k26_unit = 1, [unit = u"mol^2/kg^2", description = "Unit conversion factor"]
        k27_unit = 1, [unit = u"mol^5/kg^5", description = "Unit conversion factor"]
    end
    @equations begin
        r1.T ~ T
        r2.T ~ T
        r3.T ~ T
        r4.T ~ T
        r5.T ~ T
        r6.T ~ T
        r7.T ~ T
        r8.T ~ T
        r9.T ~ T
        r10.T ~ T
        r11.T ~ T
        r12.T ~ T
        r13.T ~ T
        r14.T ~ T
        r15.T ~ T
        r16.T ~ T
        r17.T ~ T
        r18.T ~ T
        r19.T ~ T
        r20.T ~ T
        r21.T ~ T
        r22.T ~ T
        r23.T ~ T
        r24.T ~ T
        r25.T ~ T
        r26.T ~ T
        r27.T ~ T
    end
end
