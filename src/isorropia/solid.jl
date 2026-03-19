"""
    Solid(; name=:Solid)

A solid with the given concentration.
"""
@component function Solid(; name = :Solid)
    @variables begin
        M(t),
            [description = "Molarity of the solid in air", unit = u"mol/m^3", guess = 1.0e-20]
    end

    eqs = []

    return System(eqs, t; name)
end

"""
    Solids(; name=:Solids)

Solids in Isorropia II.
"""
@component function Solids(; name = :Solids)
    # CaNO32 = Solid(; name=:CaNO32)
    # CaCl2 = Solid(; name=:CaCl2)
    # CaSO4 = Solid(; name=:CaSO4)
    # KHSO4 = Solid(; name=:KHSO4)
    # K2SO4 = Solid(; name=:K2SO4)
    # KNO3 = Solid(; name=:KNO3)
    # KCl = Solid(; name=:KCl)
    # MgSO4 = Solid(; name=:MgSO4)
    # MgNO32 = Solid(; name=:MgNO32)
    # MgCl2 = Solid(; name=:MgCl2)
    # NaCl = Solid(; name=:NaCl)
    # NaNO3 = Solid(; name=:NaNO3)
    # Na2SO4 = Solid(; name=:Na2SO4)
    # NaHSO4 = Solid(; name=:NaHSO4)
    # NH4Cl = Solid(; name=:NH4Cl)
    # NH4NO3 = Solid(; name=:NH4NO3)
    # NH42SO4 = Solid(; name=:NH42SO4)
    # NH4HSO4 = Solid(; name=:NH4HSO4)
    # NH43HSO42 = Solid(; name=:NH43HSO42)

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

    eqs = [
        #0 ~ min(NH4, Na, Ca, K, Mg, Cl, NO3, SO4)#, HSO4) # Non-negativity constraints
        #        SO4 ~ HSO4
    ]

    return System(eqs, t; name)
end
