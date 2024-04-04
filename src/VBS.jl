using ModelingToolkit, Unitful
export VBS

# The function to calculate Ci_star at different temperatures
function calc_Ci_star(T, Ci_star_standard)
    T1 = 300 # K 
    ΔH = 100000 # J/mol
    R = 8.314 # J/K/mol
    return Ci_star_standard / exp(ΔH/R*(1/(T)-1/T1))
end
@register_symbolic calc_Ci_star(T, Ci_star_standard)

function VBS(Ci)
    @parameters T = 300 [unit = u"K", description = "temperature"]
    @parameters T_unit = 1 [unit = u"K"]
    @variables Ci_star[1:8] [description = "saturation concentrations"]

    Ci_star_standard = [0.01,0.1,1,10,100,1000,10000,100000] # standard saturation concentrations
    @parameters C_OA_unit = 1 [unit = u"m^3/μg", description = "make C_OA unitless"]
    @variables C_OA [unit = u"μg/(m^3)", description = "organic compound in aerosol phase"]
    @variables ξ[1:8] [description = "partitioning coefficient at corresponding Ci_star"]
    
    eqs = [
        [ξ[i] ~ 1/(1+calc_Ci_star(T/T_unit,Ci_star_standard[i])/(C_OA*C_OA_unit)) for i ∈ 1:8]
        C_OA*C_OA_unit ~ sum(ξ.*Ci)
    ]
    
    NonlinearSystem(eqs, [C_OA, collect(ξ)...], [T,C_OA_unit,T_unit]; name=:VBS)
end

function VBS(Ci,Ci_star_standard)
    @parameters T = 300 [unit = u"K", description = "temperature"]
    @parameters T_unit = 1 [unit = u"K"]
    n = length(Ci_star_standard)
    @variables Ci_star[1:n] [description = "saturation concentrations"]

    @parameters C_OA_unit = 1 [unit = u"m^3/μg", description = "make C_OA unitless"]
    @variables C_OA [unit = u"μg/(m^3)", description = "organic compound in aerosol phase"]
    @variables ξ[1:n] [description = "partitioning coefficient at corresponding Ci_star"]
    
    eqs = [
        [ξ[i] ~ 1/(1+calc_Ci_star(T/T_unit,Ci_star_standard[i])/(C_OA*C_OA_unit)) for i ∈ 1:n]
        C_OA*C_OA_unit ~ sum(ξ.*Ci)
    ]
    
    NonlinearSystem(eqs, [C_OA, collect(ξ)...], [T,C_OA_unit,T_unit]; name=:VBS)
end
