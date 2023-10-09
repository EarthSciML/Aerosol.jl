module ISORROPIA
using EarthSciMLBase
using ModelingToolkit, Catalyst, Unitful
using IfElse

export Isorropia

include("units.jl")

"""
    Isorropia(t)

An implementation of ISORROPIA II, a model for the thermodynamic equilibrium of gas-aerosol interactions, as described in:

> Fountoukis, C. and Nenes, A., 2007. ISORROPIA II: a computationally efficient thermodynamic equilibrium model for K+–Ca 2+–Mg 2+–NH 4+–Na+–SO 4 2−–NO 3−–Cl−–H 2 O aerosols. Atmospheric Chemistry and Physics, 7(17), pp.4639-4659.

"""
struct Isorropia <: EarthSciMLODESystem
    sys::ODESystem
    rxn_sys::ReactionSystem
    Isorropia(sys::ModelingToolkit.ODESystem, rxn_sys::ReactionSystem) = new(sys, rxn_sys)
    function Isorropia(tt; name=:Isorropia)
        global t = tt # Make t global so that it can be used when defining MTK variables.

        # Load the other files into the global namespace. Ideally we would use functions 
        # instead, but there is so much systemic complexity it's difficult to encapsulate
        # everything into functions.
        include(joinpath(@__DIR__, "species.jl"))
        include(joinpath(@__DIR__, "aqueous.jl"))
        include(joinpath(@__DIR__, "deliquescence.jl"))
        include(joinpath(@__DIR__, "solid.jl"))
        include(joinpath(@__DIR__, "gas.jl"))
        include(joinpath(@__DIR__, "water.jl"))
        include(joinpath(@__DIR__, "reactions.jl"))

        statevars = [all_solids; all_ions; all_gases; I; W]
        ps = [T, RH, H2O_aq, metastable]

        @parameters so4rate=100 [unit = u"s^-1", description = "Rate of SO4 to conversion aerosol phase (pseudo-instantaneous)"] # should be a @constants
        eqs = vcat(
            # All SO4 immediately goes to aerosol phase as per Section 3.3 (item 1) of Fountoukis and Nenes (2007).
            Reaction(ifelse(SO4_g > tiny_conc, so4rate, tinyrate/unit_conc), [SO4_g], [SO4_aq]),    
            all_rxn_eqs...)

        rxn_sys = ReactionSystem(eqs, t, statevars, [so4rate; ps]; 
            systems=[IW, DRH, [rxn.sys for rxn ∈ all_rxns]...], 
            checks=true, combinatoric_ratelaws=false, name=name)
        sys = convert(ODESystem, rxn_sys)
        new(sys, rxn_sys)
    end
end
end