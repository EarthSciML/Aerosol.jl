@mtkmodel Solid begin
    @description "A solid with the given concentration."
    @variables begin
     #   m(t), [description = "Concentration of the solid in water", unit = u"mol/kg"]
        M(t), [description = "Molarity of the solid in air", unit = u"mol/m^3", guess=1]
        W(t), [description = "Aerosol water mass in air", unit = u"kg/m^3"]
    end
    @equations begin
      #  M ~ m * W
    end
end

@mtkmodel Solids begin
    @description "Solids in Isorropia II."
    @components begin
        CaNO32 = Solid()
        CaCl2 = Solid()
        CaSO4 = Solid()
        KHSO4 = Solid()
        K2SO4 = Solid()
        KNO3 = Solid()
        KCl = Solid()
        MgSO4 = Solid()
        MgNO32 = Solid()
        MgCl2 = Solid()
        NaCl = Solid()
        NaNO3 = Solid()
        Na2SO4 = Solid()
        NaHSO4 = Solid()
        NH4Cl = Solid()
        NH4NO3 = Solid()
        NH42SO4 = Solid()
        NH4HSO4 = Solid()
        NH43HSO42 = Solid()
    end
    @variables begin
        W(t), [description = "Aerosol water mass in air", unit = u"kg/m^3"]
    end
    @equations begin
        CaNO32.W ~ W
        CaCl2.W ~ W
        CaSO4.W ~ W
        KHSO4.W ~ W
        K2SO4.W ~ W
        KNO3.W ~ W
        KCl.W ~ W
        MgSO4.W ~ W
        MgNO32.W ~ W
        MgCl2.W ~ W
        NaCl.W ~ W
        NaNO3.W ~ W
        Na2SO4.W ~ W
        NaHSO4.W ~ W
        NH4Cl.W ~ W
        NH4NO3.W ~ W
        NH42SO4.W ~ W
        NH4HSO4.W ~ W
        NH43HSO42.W ~ W
    end
end


sys = mtkcompile(Solids(name=:x), fully_determined=false)

equations(sys)
unknowns(sys)
