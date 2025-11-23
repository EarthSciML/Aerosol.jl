@mtkmodel Solid begin
    @description "A solid with the given concentration."
    @variables begin
        M(t), [description = "Molarity of the solid in air", unit = u"mol/m^3", guess=1e-20]
    end
end

@mtkmodel Solids begin
    @description "Solids in Isorropia II."
    # @components begin
    #     CaNO32 = Solid()
    #     CaCl2 = Solid()
    #     CaSO4 = Solid()
    #     KHSO4 = Solid()
    #     K2SO4 = Solid()
    #     KNO3 = Solid()
    #     KCl = Solid()
    #     MgSO4 = Solid()
    #     MgNO32 = Solid()
    #     MgCl2 = Solid()
    #     NaCl = Solid()
    #     NaNO3 = Solid()
    #     Na2SO4 = Solid()
    #     NaHSO4 = Solid()
    #     NH4Cl = Solid()
    #     NH4NO3 = Solid()
    #     NH42SO4 = Solid()
    #     NH4HSO4 = Solid()
    #     NH43HSO42 = Solid()
    # end
    @variables begin
        #! format: off
        NH4(t), [description = "Molarity of NH4+ in the solid", unit = u"mol/m^3", guess=1e-14]
        Na(t), [description = "Molarity of Na+ in the solid", unit = u"mol/m^3", guess=1e-14]
        Ca(t), [description = "Molarity of Ca2+ in the solid", unit = u"mol/m^3", guess=1e-14]
        K(t), [description = "Molarity of K+ in the solid", unit = u"mol/m^3", guess=1e-14]
        Mg(t), [description = "Molarity of Mg2+ in the solid", unit = u"mol/m^3", guess=1e-14]
        Cl(t), [description = "Molarity of Cl- in the solid", unit = u"mol/m^3", guess=1e-14]
        NO3(t), [description = "Molarity of NO3- in the solid", unit = u"mol/m^3", guess=1e-14]
        #SO4(t), [description = "Molarity of SO4-- in the solid", unit = u"mol/m^3", guess=1e-14]
        #HSO4(t), [description = "Molarity of HSO4- in the solid", unit = u"mol/m^3", guess=1e-14]
        #! format: on
    end
    @equations begin
        #0 ~ min(NH4, Na, Ca, K, Mg, Cl, NO3, SO4)#, HSO4) # Non-negativity constraints
#        SO4 ~ HSO4
    end
end
