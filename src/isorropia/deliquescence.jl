"""
Create a system of equations for the deliquescence relative humidity of each salt.
"""
function deliquescence(t, RH, active_salts, salts, all_ions)
    @parameters metastable = 0 [description = "Whether the solution is 'metastable' (i.e. RH is decreasing over time). Value should be 1 for 'true' or 0 for 'false'"]
    metastable = ParentScope(metastable)
    @variables MDRH(t) = 0.3 [description = "Minimum deliquescence relative humidity (unitless)"]
    
    drh_eqs = []
    drh_vars = []
    for s ∈ active_salts
        saltname = nameof(s)
        drhname = Symbol("DRH_", saltname)
        eval(quote
            v = only(@variables $drhname($t) [description = "Deliquescence relative humidity of $(nameof($s)) (unitless)"])
            eq = v ~ $drh($s)
            push!($drh_eqs, eq)
            push!($drh_vars, v)
        end)
    end
    
    f_deliquescence = Dict()
    for s ∈ active_salts
        saltname = nameof(s)
        fname = Symbol("f_", saltname)
        drhname = Symbol("DRH_", saltname)
        eval(quote
            v = only(@variables $fname($t) [description = "Deliquescence factor for $(nameof($s)) (unitless)"])
            v = ParentScope(v)
            # Fountoukis and Nenes (2007) Eq. 22.
            # If the solution is "metastable" (i.e. RH is decreasing over time), 
            # then the deliquescence factor (one minus the fraction of solids that become liquid) is always 0.
            # Otherwise the solution is "stable" and 
            #   - Solids deliquesce (become liquid) completely when RH > DRH.
            #   - Solids do not deliquesce at all when when RH < MDRH.
            #   - Solids deliquesce partially when MDRH < RH < DRH.
            eq = v ~ ifelse($(metastable) > 0.5,
                0.0,
                min(1.0, max(0.0, ($RH - $drhname) / (min($MDRH, $drhname) - $drhname))),
                #ifelse(RH > $drhname, 0.0, 1.0)
            )
            push!($drh_eqs, eq)
            push!($drh_vars, v)
            $f_deliquescence[$s] = v
        end)
    end
    
    # TODO(CT): We don't account for the temperature dependency of MDRH. The ISORROPIA II paper 
    # doesn't say anything about this, but presumbly MDRH is temperature-dependent,
    # and the ISORROPIA I paper gives a temperature-dependent MDRH for each solution included in that 
    # model.
    @named deliquescence = ODESystem(Equation[
            drh_eqs
            MDRH ~ solution_mdrh_recurrent(1, active_salts, salts, all_ions)
        ], t, [drh_vars; MDRH], [metastable])
    return deliquescence, f_deliquescence
end

@constants T₀₃ = 298.15 [unit = u"K", description = "Standard temperature 3"]
@constants unit_T = 1 [unit = u"K", description = "Unit temperature"]    

function drh(s::SaltLike)
    if ismissing(s.l_term)
        return s.drh # If l_term doesn't exist, then DRH doesn't vary with temperature.
    end
    s.drh * exp(-s.l_term * (1.0 / T - 1.0 / T₀₃) * unit_T)
end


mdrhs = [
    ([:CaNO32, :CaCl2, :K2SO4, :KNO3, :KCl, :MgSO4, :MgNO32, :MgCl2, :NaNO3, :NaCl, :NH4NO3, :NH4Cl], 0.200),
    ([:NH42SO4, :NH4NO3, :NH4Cl, :Na2SO4, :K2SO4, :MgSO4], 0.460),
    ([:CaNO32, :K2SO4, :KNO3, :KCl, :MgSO4, :MgNO32, :MgCl2, :NaNO3, :NaCl, :NH4NO3, :NH4Cl], 0.240),
    ([:NH42SO4, :NH4Cl, :Na2SO4, :K2SO4, :MgSO4], 0.691),
    ([:CaNO32, :K2SO4, :KNO3, :KCl, :MgSO4, :MgNO32, :NaNO3, :NaCl, :NH4NO3, :NH4Cl], 0.240),
    ([:NH42SO4, :Na2SO4, :K2SO4, :MgSO4], 0.697),
    ([:K2SO4, :MgSO4, :KHSO4, :NH4HSO4, :NaHSO4, :NH42SO4, :Na2SO4, :NH43HSO42], 0.240),
    ([:NH42SO4, :NH4NO3, :Na2SO4, :K2SO4, :MgSO4], 0.494),
    ([:K2SO4, :KNO3, :KCl, :MgSO4, :MgNO32, :NaNO3, :NaCl, :NH4NO3, :NH4Cl], 0.240),
    ([:K2SO4, :MgSO4, :KHSO4, :NaHSO4, :NH42SO4, :Na2SO4, :NH43HSO42], 0.363),
    ([:K2SO4, :KNO3, :KCl, :MgSO4, :NaNO3, :NaCl, :NH4NO3, :NH4Cl], 0.596),
    ([:K2SO4, :MgSO4, :KHSO4, :NH42SO4, :Na2SO4, :NH43HSO42], 0.610),
    ([:CaNO32, :K2SO4, :KNO3, :KCl, :MgSO4, :MgNO32, :NaNO3, :NaCl, :NH4NO3, :NH4Cl], 0.240),
    ([:K2SO4, :KNO3, :KCl, :MgSO4, :MgNO32, :NaNO3, :NaCl, :NH4NO3, :NH4Cl], 0.240),
    # TODO(CT): Fountoukis and Nenes (2007) imply that this table is missing all of the mixtures that were included in ISORROPIA I.
]

@constants conc_threshold = 5.0e-10 [unit = u"mol/m_air^3", description = "Concentration threshold for whether a salt is considered present in a solution"]

"""
Return whether the given index of `mdrhs` represents the composition of the current solution
based on whether all of the ions in the given list of salts are present at
concentrations greater than `conc_threshold`, and all the ions not in the given list of salts
are present at concentrations less than `conc_threshold`.
"""
function is_solution(i, active_salts, salts, all_ions)
    active_salt_dict = filter(((k,v),) -> !isnothing(findfirst(isequal(v), active_salts)), salts)
    if length(setdiff(mdrhs[i][1], keys(active_salt_dict))) > 0
        return false # The equation system doesn't contain all of the salts for this MDRH.
    end
    ions_present = unique(vcat([[salts[s].cation.m, salts[s].anion.m] for s ∈ mdrhs[i][1]]...))
    ions_not_present = setdiff([i.m for i ∈ values(all_ions)], ions_present)
    sum(vcat([ion > conc_threshold for ion ∈ ions_present],
        [ion < conc_threshold for ion ∈ ions_not_present])) > length(all_ions) - 0.5
end

function solution_mdrh_recurrent(i, active_salts, salts, all_ions)
    if i == length(mdrhs)
        return ifelse(is_solution(i, active_salts, salts, all_ions), mdrhs[i][2],
            sum([mdrh[2] for mdrh ∈ mdrhs]) / length(mdrhs) # Average of all MDRHs
        )
    end
    ifelse(is_solution(i, active_salts, salts, all_ions), mdrhs[i][2], solution_mdrh_recurrent(i + 1, active_salts, salts, all_ions))
end