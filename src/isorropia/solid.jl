@mtkmodel Solid begin
    @description "A solid with the given concentration."
    @variables begin
        M(t), [description = "Molarity of the solid in air", unit = u"mol/m^3", guess=1e-8]
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
end

@named xxx = Solids()
sys = mtkcompile(Solids(name=:x), fully_determined=false)

equations(sys)
unknowns(sys)
