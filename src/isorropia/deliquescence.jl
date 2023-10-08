@parameters metastable = 0 [description = "Whether the solution is 'metastable' (i.e. RH is decreasing over time). Value should be 1 for 'true' or 0 for 'false'"]
metastable = ParentScope(metastable)
@variables MDRH(t) = 0.3 [description = "Minimum deliquescence relative humidity (unitless)"]

@constants T₀₃ = 298.15 [unit = u"K", description = "Standard temperature 3"]
@constants unit_T = 1 [unit = u"K", description = "Unit temperature"]
function drh(s::SaltLike)
    if ismissing(s.l_term)
        return s.drh # If l_term doesn't exist, then DRH doesn't vary with temperature.
    end
    s.drh * exp(-s.l_term * (1.0 / T - 1.0 / T₀₃) * unit_T)
end

drh_eqs = []
drh_vars = []
for s ∈ all_salts
    saltname = nameof(s)
    drhname = Symbol("DRH_", saltname)
    eval(quote
        v = (@variables $drhname($t) [description = "Deliquescence relative humidity of $(nameof($s)) (unitless)"])[1]
        eq = v ~ $drh($s)
        push!(drh_eqs, eq)
        push!(drh_vars, v)
    end)
end

f_deliquescence = Dict()
for s ∈ all_salts
    saltname = nameof(s)
    fname = Symbol("f_", saltname)
    drhname = Symbol("DRH_", saltname)
    eval(quote
        v = (@variables $fname($t) [description = "Deliquescence factor for $(nameof($s)) (unitless)"])[1]
        v = ParentScope(v)
        # Fountoukis and Nenes (2007) Eq. 22.
        # If the solution is "metastable" (i.e. RH is decreasing over time), 
        # then the deliquescence factor (one minus the fraction of solids that become liquid) is always 0.
        # Otherwise the solution is "stable" and 
        #   - Solids deliquesce (become liquid) completely when RH > DRH.
        #   - Solids do not deliquesce at all when when RH < MDRH.
        #   - Solids deliquesce partially when MDRH < RH < DRH.
        eq = v ~ ifelse(metastable > 0.5,
            0.0,
            min(1.0, max(0.0, (RH - $drhname) / (min(MDRH, $drhname) - $drhname))),
            #ifelse(RH > $drhname, 0.0, 1.0)
        )
        push!(drh_eqs, eq)
        push!(drh_vars, v)
        f_deliquescence[$s] = v
    end)
end

mdrhs = [
    ([CaNO32_aqs, CaCl2_aqs, K2SO4_aqs, KNO3_aqs, KCl_aqs, MgSO4_aqs, MgNO32_aqs, MgCl2_aqs, NaNO3_aqs, NaCl_aqs, NH4NO3_aqs, NH4Cl_aqs], 0.200),
    ([NH42SO4_aqs, NH4NO3_aqs, NH4Cl_aqs, Na2SO4_aqs, K2SO4_aqs, MgSO4_aqs], 0.460),
    ([CaNO32_aqs, K2SO4_aqs, KNO3_aqs, KCl_aqs, MgSO4_aqs, MgNO32_aqs, MgCl2_aqs, NaNO3_aqs, NaCl_aqs, NH4NO3_aqs, NH4Cl_aqs], 0.240),
    ([NH42SO4_aqs, NH4Cl_aqs, Na2SO4_aqs, K2SO4_aqs, MgSO4_aqs], 0.691),
    ([CaNO32_aqs, K2SO4_aqs, KNO3_aqs, KCl_aqs, MgSO4_aqs, MgNO32_aqs, NaNO3_aqs, NaCl_aqs, NH4NO3_aqs, NH4Cl_aqs], 0.240),
    ([NH42SO4_aqs, Na2SO4_aqs, K2SO4_aqs, MgSO4_aqs], 0.697),
    ([K2SO4_aqs, MgSO4_aqs, KHSO4_aqs, NH4HSO4_aqs, NaHSO4_aqs, NH42SO4_aqs, Na2SO4_aqs, NH43HSO42_aqs], 0.240),
    ([NH42SO4_aqs, NH4NO3_aqs, Na2SO4_aqs, K2SO4_aqs, MgSO4_aqs], 0.494),
    ([K2SO4_aqs, KNO3_aqs, KCl_aqs, MgSO4_aqs, MgNO32_aqs, NaNO3_aqs, NaCl_aqs, NH4NO3_aqs, NH4Cl_aqs], 0.240),
    ([K2SO4_aqs, MgSO4_aqs, KHSO4_aqs, NaHSO4_aqs, NH42SO4_aqs, Na2SO4_aqs, NH43HSO42_aqs], 0.363),
    ([K2SO4_aqs, KNO3_aqs, KCl_aqs, MgSO4_aqs, NaNO3_aqs, NaCl_aqs, NH4NO3_aqs, NH4Cl_aqs], 0.596),
    ([K2SO4_aqs, MgSO4_aqs, KHSO4_aqs, NH42SO4_aqs, Na2SO4_aqs, NH43HSO42_aqs], 0.610),
    ([CaNO32_aqs, K2SO4_aqs, KNO3_aqs, KCl_aqs, MgSO4_aqs, MgNO32_aqs, NaNO3_aqs, NaCl_aqs, NH4NO3_aqs, NH4Cl_aqs], 0.240),
    ([K2SO4_aqs, KNO3_aqs, KCl_aqs, MgSO4_aqs, MgNO32_aqs, NaNO3_aqs, NaCl_aqs, NH4NO3_aqs, NH4Cl_aqs], 0.240),
    # TODO(CT): Fountoukis and Nenes (2007) imply that this table is missing all of the mixtures that were included in ISORROPIA I.
]

@constants conc_threshold=1.0e-15 [unit = u"mol/m_air^3", description = "Concentration threshold for whether a salt is considered present in a solution"]

"""
Return whether the given index of `mdrhs` represents the composition of the current solution
based on whether all of the ions in the given list of salts are present at
concentrations greater than `conc_threshold`, and all the ions not in the given list of salts
are present at concentrations less than `conc_threshold`.
"""
function is_solution(i)
    ions_present = unique(vcat([[s.cation.m, s.anion.m] for s ∈ mdrhs[i][1]]...))
    ions_not_present = setdiff(all_ions, ions_present)
    sum(vcat([ion > conc_threshold for ion ∈ ions_present],
        [ion < conc_threshold for ion ∈ ions_not_present])) > length(all_ions) - 0.5
end

function solution_mdrh_recurrent(i)
    if i == length(mdrhs)
        return ifelse(is_solution(i), mdrhs[i][2], 
            sum([mdrh[2] for mdrh ∈ mdrhs]) / length(mdrhs) # Average of all MDRHs
        )
    end
    ifelse(is_solution(i), mdrhs[i][2], solution_mdrh_recurrent(i + 1))
end

# TODO(CT): We don't account for the temperature dependency of MDRH. The ISORROPIA II paper 
# doesn't say anything about this, but presumbly MDRH is temperature-dependent,
# and the ISORROPIA I paper gives a temperature-dependent MDRH for each solution included in that 
# model.
@named DRH = NonlinearSystem([
    drh_eqs;
    MDRH ~ solution_mdrh_recurrent(1); 
], [drh_vars; MDRH], [metastable])
