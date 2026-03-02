# =============================================================================
# ISORROPIA II Thermodynamic Equilibrium Model (Metastable State)
#
# Implements the inorganic aerosol thermodynamic equilibrium model for the
# K⁺–Ca²⁺–Mg²⁺–NH₄⁺–Na⁺–SO₄²⁻–NO₃⁻–Cl⁻–H₂O system.
#
# Reference: Fountoukis, C. and Nenes, A.: ISORROPIA II: a computationally
#   efficient thermodynamic equilibrium model for
#   K⁺–Ca²⁺–Mg²⁺–NH₄⁺–Na⁺–SO₄²⁻–NO₃⁻–Cl⁻–H₂O aerosols,
#   Atmos. Chem. Phys., 7, 4639–4659, 2007.
#
# Numerical values from: USEPA CMAQ isocom.f (CCTM/src/aero/aero6/)
# =============================================================================

# =============================================================================
# Part 1: Thermodynamic Data
# =============================================================================

# Reference temperature for helper functions (needed at module level for @register_symbolic)
# Main thermodynamic constants moved to @constants inside component
const _ISO2_T0 = 298.15  # K

# --------------------------------------------------------------------------
# Equilibrium constants at 298.15K with temperature dependence parameters
# Van't Hoff equation (Eq. 5):
#   K(T) = K₀ × exp(A×(T₀/T - 1) + B×(1 + ln(T₀/T) - T₀/T))
# where A = ΔH°(T₀)/(R·T₀), B = Δcₚ°/R
# Format: (K0, A, B)
# --------------------------------------------------------------------------

# van't Hoff equilibrium constant parameters: (K0, A, B) where A = ΔH°/(R·T₀), B = Δcₚ°/R
# Values are moved to @constants inside component following ModelingToolkit standards

# Solid salt dissolution constants (not used in metastable implementation)
# These would be needed for stable solution branches - included in constants below

# --------------------------------------------------------------------------
# Kusik-Meissner Q parameters for binary activity coefficients (from CMAQ isocom.f)
# Format: (q_parameter, charge_product_z)
# Note: Module-level const needed for @register_symbolic functions
# --------------------------------------------------------------------------
const _ISO2_KM_PARAMS = Dict(
    :NaCl => (2.23, 1),
    :Na2SO4 => (-0.19, 2),
    :NaNO3 => (-0.39, 1),
    :NH42SO4 => (-0.25, 2),
    :NH4NO3 => (-1.15, 1),
    :NH4Cl => (0.82, 1),
    :H2SO4 => (-0.1, 2),
    :HHSO4 => (8.0, 1),
    :HNO3 => (2.6, 1),
    :HCl => (6.0, 1),
    :CaNO32 => (0.93, 2),
    :CaCl2 => (2.4, 2),
    :K2SO4 => (-0.25, 2),
    :KNO3 => (-2.33, 1),
    :KCl => (0.92, 1),
    :MgSO4 => (0.15, 4),
    :MgNO32 => (2.32, 2),
    :MgCl2 => (2.9, 2),
)

# --------------------------------------------------------------------------
# ZSR Binary Molality Lookup Tables
# 100 entries from aw = 0.01 to 1.00 (index i → aw ≈ i/100)
# Values are molality (mol/kg water) at the given water activity
# Source: CMAQ isocom.f BLOCK DATA BLKISO (AIM Model III database)
# Note: These are module-level constants for performance (used in registered functions)
# --------------------------------------------------------------------------

# (NH4)2SO4 - AWAS
const _ISO2_ZSR_NH42SO4 = [
    187.72, 187.72, 187.72, 187.72, 187.72, 187.72, 187.72, 187.72, 187.72, 187.72,
    158.13, 134.41, 115.37, 100.1, 87.86, 78.0, 70.0, 63.45, 58.02, 53.46,
    49.59, 46.26, 43.37, 40.84, 38.59, 36.59, 34.79, 33.16, 31.67, 30.31,
    29.07, 27.91, 26.84, 25.84, 24.91, 24.03, 23.21, 22.44, 21.7, 21.01,
    20.34, 19.71, 19.11, 18.54, 17.99, 17.46, 16.95, 16.46, 15.99, 15.54,
    15.1, 14.67, 14.26, 13.86, 13.47, 13.09, 12.72, 12.36, 12.01, 11.67,
    11.33, 11.0, 10.68, 10.37, 10.06, 9.75, 9.45, 9.15, 8.86, 8.57,
    8.29, 8.01, 7.73, 7.45, 7.18, 6.91, 6.64, 6.37, 6.1, 5.83,
    5.56, 5.29, 5.02, 4.74, 4.47, 4.19, 3.91, 3.63, 3.34, 3.05,
    2.75, 2.45, 2.14, 1.83, 1.51, 1.19, 0.87, 0.56, 0.26, 0.1,
]

# Na2SO4 - AWSS
const _ISO2_ZSR_Na2SO4 = [
    24.1, 24.1, 24.1, 24.1, 24.1, 24.1, 24.1, 24.1, 24.1, 24.1,
    23.17, 22.34, 21.58, 20.9, 20.27, 19.69, 19.15, 18.64, 18.17, 17.72,
    17.3, 16.9, 16.52, 16.16, 15.81, 15.48, 15.16, 14.85, 14.55, 14.27,
    13.99, 13.73, 13.47, 13.21, 12.97, 12.73, 12.5, 12.27, 12.05, 11.84,
    11.62, 11.42, 11.21, 11.01, 10.82, 10.63, 10.44, 10.25, 10.07, 9.89,
    9.71, 9.53, 9.36, 9.19, 9.02, 8.85, 8.68, 8.51, 8.35, 8.19,
    8.02, 7.86, 7.7, 7.54, 7.38, 7.22, 7.06, 6.9, 6.74, 6.58,
    6.42, 6.26, 6.1, 5.94, 5.78, 5.61, 5.45, 5.28, 5.11, 4.93,
    4.76, 4.58, 4.39, 4.2, 4.01, 3.81, 3.6, 3.39, 3.16, 2.93,
    2.68, 2.41, 2.13, 1.83, 1.52, 1.19, 0.86, 0.54, 0.25, 0.1,
]

# NaNO3 - AWSN
const _ISO2_ZSR_NaNO3 = [
    394.54, 394.54, 394.54, 394.54, 394.54, 394.54, 394.54, 394.54, 394.54, 394.54,
    338.91, 293.01, 254.73, 222.61, 195.56, 172.76, 153.53, 137.32, 123.65, 112.08,
    102.26, 93.88, 86.68, 80.45, 75.02, 70.24, 66.02, 62.26, 58.89, 55.85,
    53.09, 50.57, 48.26, 46.14, 44.17, 42.35, 40.65, 39.06, 37.57, 36.17,
    34.85, 33.6, 32.42, 31.29, 30.22, 29.2, 28.22, 27.28, 26.39, 25.52,
    24.69, 23.89, 23.12, 22.37, 21.65, 20.94, 20.26, 19.6, 18.96, 18.33,
    17.72, 17.12, 16.53, 15.96, 15.4, 14.85, 14.31, 13.78, 13.26, 12.75,
    12.25, 11.75, 11.26, 10.77, 10.29, 9.82, 9.35, 8.88, 8.42, 7.97,
    7.52, 7.07, 6.62, 6.18, 5.75, 5.32, 4.89, 4.47, 4.05, 3.64,
    3.24, 2.84, 2.45, 2.07, 1.7, 1.34, 0.99, 0.65, 0.31, 0.1,
]

# NaCl - AWSC
const _ISO2_ZSR_NaCl = [
    28.16, 28.16, 28.16, 28.16, 28.16, 28.16, 28.16, 28.16, 28.16, 28.16,
    27.17, 26.27, 25.45, 24.69, 23.98, 23.33, 22.72, 22.14, 21.59, 21.08,
    20.58, 20.12, 19.67, 19.24, 18.82, 18.43, 18.04, 17.67, 17.32, 16.97,
    16.63, 16.31, 15.99, 15.68, 15.38, 15.08, 14.79, 14.51, 14.24, 13.97,
    13.7, 13.44, 13.18, 12.93, 12.68, 12.44, 12.2, 11.96, 11.73, 11.5,
    11.27, 11.05, 10.82, 10.6, 10.38, 10.16, 9.95, 9.74, 9.52, 9.31,
    9.1, 8.89, 8.69, 8.48, 8.27, 8.07, 7.86, 7.65, 7.45, 7.24,
    7.04, 6.83, 6.62, 6.42, 6.21, 6.0, 5.79, 5.58, 5.36, 5.15,
    4.93, 4.71, 4.48, 4.26, 4.03, 3.8, 3.56, 3.32, 3.07, 2.82,
    2.57, 2.3, 2.04, 1.76, 1.48, 1.2, 0.91, 0.61, 0.3, 0.1,
]

# NH4NO3 - AWAN
const _ISO2_ZSR_NH4NO3 = [
    960.19, 960.19, 960.19, 960.19, 960.19, 960.19, 960.19, 960.19, 960.19, 960.19,
    853.15, 763.85, 688.2, 623.27, 566.92, 517.54, 473.91, 435.06, 400.26, 368.89,
    340.48, 314.63, 291.01, 269.36, 249.46, 231.11, 214.17, 198.5, 184.0, 170.58,
    158.15, 146.66, 136.04, 126.25, 117.24, 108.97, 101.39, 94.45, 88.11, 82.33,
    77.06, 72.25, 67.85, 63.84, 60.16, 56.78, 53.68, 50.81, 48.17, 45.71,
    43.43, 41.31, 39.32, 37.46, 35.71, 34.06, 32.5, 31.03, 29.63, 28.3,
    27.03, 25.82, 24.67, 23.56, 22.49, 21.47, 20.48, 19.53, 18.61, 17.72,
    16.86, 16.02, 15.2, 14.41, 13.64, 12.89, 12.15, 11.43, 10.73, 10.05,
    9.38, 8.73, 8.09, 7.47, 6.86, 6.27, 5.7, 5.15, 4.61, 4.09,
    3.6, 3.12, 2.66, 2.23, 1.81, 1.41, 1.03, 0.67, 0.32, 0.1,
]

# NH4Cl - AWAC
const _ISO2_ZSR_NH4Cl = [
    1209.0, 1209.0, 1209.0, 1209.0, 1209.0, 1209.0, 1209.0, 1209.0, 1209.0, 1209.0,
    1067.6, 949.27, 848.62, 761.82, 686.04, 619.16, 559.55, 505.92, 457.25, 412.69,
    371.55, 333.21, 297.13, 262.81, 229.78, 197.59, 165.98, 135.49, 108.57, 88.29,
    74.4, 64.75, 57.69, 52.25, 47.9, 44.3, 41.27, 38.65, 36.36, 34.34,
    32.52, 30.88, 29.39, 28.02, 26.76, 25.6, 24.51, 23.5, 22.55, 21.65,
    20.8, 20.0, 19.24, 18.52, 17.83, 17.17, 16.54, 15.93, 15.35, 14.79,
    14.25, 13.73, 13.22, 12.73, 12.26, 11.8, 11.35, 10.92, 10.49, 10.08,
    9.67, 9.28, 8.89, 8.51, 8.14, 7.77, 7.42, 7.06, 6.72, 6.37,
    6.03, 5.7, 5.37, 5.05, 4.72, 4.4, 4.08, 3.77, 3.45, 3.14,
    2.82, 2.51, 2.2, 1.89, 1.57, 1.26, 0.94, 0.62, 0.31, 0.1,
]

# NH4HSO4 - AWAB
const _ISO2_ZSR_NH4HSO4 = [
    312.84, 312.84, 312.84, 312.84, 312.84, 312.84, 312.84, 312.84, 312.84, 312.84,
    271.43, 237.19, 208.52, 184.28, 163.64, 145.97, 130.79, 117.72, 106.42, 96.64,
    88.16, 80.77, 74.33, 68.67, 63.7, 59.3, 55.39, 51.89, 48.76, 45.93,
    43.38, 41.05, 38.92, 36.97, 35.18, 33.52, 31.98, 30.55, 29.22, 27.98,
    26.81, 25.71, 24.67, 23.7, 22.77, 21.9, 21.06, 20.27, 19.52, 18.8,
    18.11, 17.45, 16.82, 16.21, 15.63, 15.07, 14.53, 14.01, 13.51, 13.02,
    12.56, 12.1, 11.66, 11.24, 10.82, 10.42, 10.04, 9.66, 9.29, 8.93,
    8.58, 8.24, 7.91, 7.58, 7.26, 6.95, 6.65, 6.35, 6.05, 5.76,
    5.48, 5.2, 4.92, 4.64, 4.37, 4.09, 3.82, 3.54, 3.27, 2.99,
    2.7, 2.42, 2.12, 1.83, 1.52, 1.22, 0.9, 0.59, 0.28, 0.1,
]

# H2SO4 - AWSA (also used for H-HSO4)
const _ISO2_ZSR_H2SO4 = [
    34.0, 33.56, 29.22, 26.55, 24.61, 23.11, 21.89, 20.87, 19.99, 18.45,
    17.83, 17.26, 16.73, 16.25, 15.8, 15.38, 14.98, 14.61, 14.26, 13.93,
    13.61, 13.3, 13.01, 12.73, 12.47, 12.21, 11.96, 11.72, 11.49, 11.26,
    11.04, 10.83, 10.62, 10.42, 10.23, 10.03, 9.85, 9.67, 9.49, 9.31,
    9.14, 8.97, 8.81, 8.65, 8.49, 8.33, 8.18, 8.02, 7.87, 7.73,
    7.58, 7.44, 7.29, 7.15, 7.01, 6.88, 6.74, 6.61, 6.47, 6.34,
    6.21, 6.07, 5.94, 5.81, 5.68, 5.55, 5.43, 5.3, 5.17, 5.04,
    4.91, 4.78, 4.65, 4.52, 4.39, 4.26, 4.13, 4.0, 3.86, 3.73,
    3.59, 3.45, 3.31, 3.17, 3.02, 2.87, 2.71, 2.56, 2.39, 2.22,
    2.05, 1.87, 1.68, 1.48, 1.27, 1.04, 0.8, 0.55, 0.28, 0.1,
]

# NaHSO4 - AWSB (also used for KHSO4)
const _ISO2_ZSR_NaHSO4 = [
    55.99, 55.99, 55.99, 55.99, 55.99, 55.99, 55.99, 55.99, 55.99, 55.99,
    53.79, 51.81, 49.99, 48.31, 46.75, 45.28, 43.91, 42.62, 41.39, 40.22,
    39.1, 38.02, 36.99, 36.0, 35.04, 34.11, 33.21, 32.34, 31.49, 30.65,
    29.84, 29.04, 28.27, 27.5, 26.75, 26.01, 25.29, 24.57, 23.87, 23.17,
    22.49, 21.81, 21.15, 20.49, 19.84, 19.21, 18.58, 17.97, 17.37, 16.77,
    16.19, 15.63, 15.08, 14.54, 14.01, 13.51, 13.01, 12.53, 12.07, 11.62,
    11.19, 10.77, 10.36, 9.97, 9.59, 9.23, 8.87, 8.53, 8.2, 7.88,
    7.57, 7.27, 6.97, 6.69, 6.41, 6.14, 5.88, 5.62, 5.36, 5.11,
    4.87, 4.63, 4.39, 4.15, 3.92, 3.68, 3.45, 3.21, 2.98, 2.74,
    2.49, 2.24, 1.98, 1.72, 1.44, 1.16, 0.87, 0.57, 0.28, 0.1,
]

# Letovicite (NH4)3H(SO4)2 - AWLC
const _ISO2_ZSR_LC = [
    125.37, 125.37, 125.37, 125.37, 125.37, 125.37, 125.37, 125.37, 125.37, 125.37,
    110.1, 97.5, 86.98, 78.08, 70.49, 63.97, 58.33, 53.43, 49.14, 45.36,
    42.03, 39.07, 36.44, 34.08, 31.97, 30.06, 28.33, 26.76, 25.32, 24.01,
    22.81, 21.7, 20.67, 19.71, 18.83, 18.0, 17.23, 16.5, 15.82, 15.18,
    14.58, 14.01, 13.46, 12.95, 12.46, 11.99, 11.55, 11.13, 10.72, 10.33,
    9.96, 9.6, 9.26, 8.93, 8.61, 8.3, 8.0, 7.72, 7.44, 7.17,
    6.91, 6.66, 6.42, 6.19, 5.96, 5.74, 5.52, 5.31, 5.11, 4.91,
    4.71, 4.53, 4.34, 4.16, 3.99, 3.81, 3.64, 3.48, 3.31, 3.15,
    2.99, 2.84, 2.68, 2.53, 2.37, 2.22, 2.06, 1.91, 1.75, 1.6,
    1.44, 1.28, 1.12, 0.95, 0.79, 0.62, 0.45, 0.29, 0.14, 0.1,
]

# Ca(NO3)2 - AWCN
const _ISO2_ZSR_CaNO32 = [
    32.89, 31.46, 30.12, 28.84, 27.64, 26.51, 25.44, 24.44, 23.49, 22.59,
    21.75, 20.96, 20.22, 19.51, 18.85, 18.23, 17.64, 17.09, 16.56, 16.07,
    15.61, 15.17, 14.75, 14.36, 13.99, 13.63, 13.3, 12.98, 12.68, 12.39,
    12.11, 11.84, 11.59, 11.35, 11.11, 10.88, 10.66, 10.45, 10.24, 10.04,
    9.84, 9.65, 9.46, 9.28, 9.1, 8.92, 8.74, 8.57, 8.4, 8.23,
    8.06, 7.9, 7.73, 7.57, 7.41, 7.25, 7.1, 6.94, 6.79, 6.63,
    6.48, 6.33, 6.18, 6.03, 5.89, 5.74, 5.6, 5.46, 5.32, 5.17,
    5.04, 4.9, 4.76, 4.62, 4.49, 4.35, 4.22, 4.08, 3.94, 3.8,
    3.66, 3.52, 3.38, 3.23, 3.08, 2.93, 2.77, 2.6, 2.43, 2.25,
    2.07, 1.87, 1.67, 1.45, 1.22, 0.97, 0.72, 0.44, 0.14, 0.1,
]

# CaCl2 - AWCC
const _ISO2_ZSR_CaCl2 = [
    19.9, 19.0, 18.15, 17.35, 16.6, 15.89, 15.22, 14.58, 13.99, 13.43,
    12.9, 12.41, 11.94, 11.5, 11.09, 10.7, 10.34, 9.99, 9.67, 9.37,
    9.09, 8.83, 8.57, 8.34, 8.12, 7.91, 7.71, 7.53, 7.35, 7.19,
    7.03, 6.88, 6.74, 6.6, 6.47, 6.35, 6.23, 6.12, 6.01, 5.9,
    5.8, 5.7, 5.61, 5.51, 5.42, 5.33, 5.24, 5.16, 5.07, 4.99,
    4.91, 4.82, 4.74, 4.66, 4.58, 4.5, 4.42, 4.34, 4.26, 4.19,
    4.11, 4.03, 3.95, 3.87, 3.79, 3.72, 3.64, 3.56, 3.48, 3.41,
    3.33, 3.25, 3.17, 3.09, 3.01, 2.93, 2.85, 2.76, 2.68, 2.59,
    2.5, 2.41, 2.32, 2.23, 2.13, 2.03, 1.93, 1.82, 1.71, 1.59,
    1.47, 1.35, 1.22, 1.07, 0.93, 0.77, 0.61, 0.44, 0.25, 0.1,
]

# K2SO4 - AWPS
const _ISO2_ZSR_K2SO4 = [
    1014.82, 969.72, 926.16, 884.11, 843.54, 804.41, 766.68, 730.32, 695.3, 661.58,
    629.14, 597.93, 567.92, 539.09, 511.41, 484.83, 459.34, 434.89, 411.47, 389.04,
    367.58, 347.05, 327.43, 308.69, 290.8, 273.73, 257.47, 241.98, 227.24, 213.22,
    199.9, 187.26, 175.27, 163.91, 153.15, 142.97, 133.36, 124.28, 115.73, 107.66,
    100.08, 92.95, 86.26, 79.99, 74.12, 68.63, 63.5, 58.73, 54.27, 50.14,
    46.3, 42.74, 39.44, 36.4, 33.59, 31.0, 28.63, 26.45, 24.45, 22.62,
    20.95, 19.43, 18.05, 16.79, 15.64, 14.61, 13.66, 12.81, 12.03, 11.33,
    10.68, 10.09, 9.55, 9.06, 8.6, 8.17, 7.76, 7.38, 7.02, 6.66,
    6.32, 5.98, 5.65, 5.31, 4.98, 4.64, 4.31, 3.96, 3.62, 3.27,
    2.92, 2.57, 2.22, 1.87, 1.53, 1.2, 0.87, 0.57, 0.28, 0.1,
]

# KNO3 - AWPN (first 44 entries capped at 1000)
const _ISO2_ZSR_KNO3 = [
    1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0,
    1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0,
    1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0,
    1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0,
    1000.0, 1000.0, 1000.0, 953.05, 881.09, 813.39, 749.78, 690.09, 634.14, 581.77,
    532.83, 487.16, 444.61, 405.02, 368.26, 334.18, 302.64, 273.51, 246.67, 221.97,
    199.31, 178.56, 159.6, 142.33, 126.63, 112.4, 99.54, 87.96, 77.55, 68.24,
    59.92, 52.53, 45.98, 40.2, 35.11, 30.65, 26.75, 23.35, 20.4, 17.85,
    15.63, 13.72, 12.06, 10.61, 9.35, 8.24, 7.25, 6.37, 5.56, 4.82,
    4.12, 3.47, 2.86, 2.28, 1.74, 1.24, 0.79, 0.4, 0.2, 0.1,
]

# KCl - AWPC
const _ISO2_ZSR_KCl = [
    172.62, 165.75, 159.1, 152.67, 146.46, 140.45, 134.64, 129.03, 123.61, 118.38,
    113.34, 108.48, 103.79, 99.27, 94.93, 90.74, 86.71, 82.84, 79.11, 75.53,
    72.09, 68.79, 65.63, 62.59, 59.68, 56.9, 54.23, 51.68, 49.24, 46.91,
    44.68, 42.56, 40.53, 38.6, 36.76, 35.0, 33.33, 31.75, 30.24, 28.81,
    27.45, 26.16, 24.94, 23.78, 22.68, 21.64, 20.66, 19.74, 18.86, 18.03,
    17.25, 16.51, 15.82, 15.16, 14.54, 13.96, 13.41, 12.89, 12.4, 11.94,
    11.5, 11.08, 10.69, 10.32, 9.96, 9.62, 9.3, 8.99, 8.69, 8.4,
    8.12, 7.85, 7.59, 7.33, 7.08, 6.83, 6.58, 6.33, 6.08, 5.84,
    5.59, 5.34, 5.09, 4.83, 4.57, 4.31, 4.04, 3.76, 3.48, 3.19,
    2.9, 2.6, 2.29, 1.98, 1.66, 1.33, 0.99, 0.65, 0.3, 0.1,
]

# MgSO4 - AWMS
const _ISO2_ZSR_MgSO4 = [
    0.93, 2.5, 3.94, 5.25, 6.45, 7.54, 8.52, 9.4, 10.19, 10.89,
    11.5, 12.04, 12.51, 12.9, 13.23, 13.5, 13.72, 13.88, 13.99, 14.07,
    14.1, 14.09, 14.05, 13.98, 13.88, 13.75, 13.6, 13.43, 13.25, 13.05,
    12.83, 12.61, 12.37, 12.13, 11.88, 11.63, 11.37, 11.12, 10.86, 10.6,
    10.35, 10.09, 9.85, 9.6, 9.36, 9.13, 8.9, 8.68, 8.47, 8.26,
    8.07, 7.87, 7.69, 7.52, 7.35, 7.19, 7.03, 6.89, 6.75, 6.62,
    6.49, 6.37, 6.26, 6.15, 6.04, 5.94, 5.84, 5.75, 5.65, 5.56,
    5.47, 5.38, 5.29, 5.2, 5.11, 5.01, 4.92, 4.82, 4.71, 4.6,
    4.49, 4.36, 4.24, 4.1, 3.96, 3.81, 3.65, 3.48, 3.3, 3.11,
    2.92, 2.71, 2.49, 2.26, 2.02, 1.76, 1.5, 1.22, 0.94, 0.64,
]

# Mg(NO3)2 - AWMN
const _ISO2_ZSR_MgNO32 = [
    12.0, 11.84, 11.68, 11.52, 11.36, 11.2, 11.04, 10.88, 10.72, 10.56,
    10.4, 10.25, 10.09, 9.93, 9.78, 9.63, 9.47, 9.32, 9.17, 9.02,
    8.87, 8.72, 8.58, 8.43, 8.29, 8.15, 8.01, 7.87, 7.73, 7.59,
    7.46, 7.33, 7.2, 7.07, 6.94, 6.82, 6.69, 6.57, 6.45, 6.33,
    6.21, 6.01, 5.98, 5.87, 5.76, 5.65, 5.55, 5.44, 5.34, 5.24,
    5.14, 5.04, 4.94, 4.84, 4.75, 4.66, 4.56, 4.47, 4.38, 4.29,
    4.21, 4.12, 4.03, 3.95, 3.86, 3.78, 3.69, 3.61, 3.53, 3.45,
    3.36, 3.28, 3.19, 3.11, 3.03, 2.94, 2.85, 2.76, 2.67, 2.58,
    2.49, 2.39, 2.3, 2.2, 2.1, 1.99, 1.88, 1.77, 1.66, 1.54,
    1.42, 1.29, 1.16, 1.02, 0.88, 0.73, 0.58, 0.42, 0.25, 0.1,
]

# MgCl2 - AWMC
const _ISO2_ZSR_MgCl2 = [
    11.24, 10.99, 10.74, 10.5, 10.26, 10.03, 9.81, 9.59, 9.38, 9.18,
    8.98, 8.79, 8.6, 8.42, 8.25, 8.07, 7.91, 7.75, 7.59, 7.44,
    7.29, 7.15, 7.01, 6.88, 6.75, 6.62, 6.5, 6.38, 6.27, 6.16,
    6.05, 5.94, 5.85, 5.75, 5.65, 5.56, 5.47, 5.38, 5.3, 5.22,
    5.14, 5.06, 4.98, 4.91, 4.84, 4.77, 4.7, 4.63, 4.57, 4.5,
    4.44, 4.37, 4.31, 4.25, 4.19, 4.13, 4.07, 4.01, 3.95, 3.89,
    3.83, 3.77, 3.71, 3.65, 3.58, 3.52, 3.46, 3.39, 3.33, 3.26,
    3.19, 3.12, 3.05, 2.98, 2.9, 2.82, 2.75, 2.67, 2.58, 2.49,
    2.41, 2.32, 2.22, 2.13, 2.03, 1.92, 1.82, 1.71, 1.6, 1.48,
    1.36, 1.24, 1.11, 0.98, 0.84, 0.7, 0.56, 0.41, 0.25, 0.1,
]

# Collect all ZSR tables indexed by pair ID matching CMAQ convention
const _ISO2_ZSR_TABLES = Dict(
    1 => _ISO2_ZSR_NaCl,
    2 => _ISO2_ZSR_Na2SO4,
    3 => _ISO2_ZSR_NaNO3,
    4 => _ISO2_ZSR_NH42SO4,
    5 => _ISO2_ZSR_NH4NO3,
    6 => _ISO2_ZSR_NH4Cl,
    7 => _ISO2_ZSR_H2SO4,
    8 => _ISO2_ZSR_H2SO4,      # H-HSO4 uses same table as H2SO4
    9 => _ISO2_ZSR_NH4HSO4,
    12 => _ISO2_ZSR_NaHSO4,
    13 => _ISO2_ZSR_LC,
    15 => _ISO2_ZSR_CaNO32,
    16 => _ISO2_ZSR_CaCl2,
    17 => _ISO2_ZSR_K2SO4,
    18 => _ISO2_ZSR_NaHSO4,    # KHSO4 uses NaHSO4 table
    19 => _ISO2_ZSR_KNO3,
    20 => _ISO2_ZSR_KCl,
    21 => _ISO2_ZSR_MgSO4,
    22 => _ISO2_ZSR_MgNO32,
    23 => _ISO2_ZSR_MgCl2,
)

# --------------------------------------------------------------------------
# DRH values at 298.15K and temperature correction coefficients (Table 4, Eq. 17)
# DRH(T) = DRH(T₀) × exp(coeff × (1/T - 1/T₀))
# Format: (DRH_at_298K, temperature_correction_coefficient_in_K)
# --------------------------------------------------------------------------
const _ISO2_DRH = Dict(
    :NH42SO4 => (0.7997, 80.0),
    :NH4HSO4 => (0.4, 384.0),
    :LC => (0.69, 186.0),
    :NH4NO3 => (0.6183, 852.0),
    :NH4Cl => (0.771, 239.0),
    :NaCl => (0.7528, 25.0),
    :Na2SO4 => (0.93, 80.0),
    :NaNO3 => (0.7379, 304.0),
    :NaHSO4 => (0.52, -45.0),
    :CaSO4 => (0.97, 0.0),
    :CaNO32 => (0.4906, 509.4),
    :CaCl2 => (0.283, 551.1),
    :K2SO4 => (0.9751, 35.6),
    :KHSO4 => (0.86, 0.0),
    :KNO3 => (0.9248, 0.0),
    :KCl => (0.8426, 158.9),
    :MgSO4 => (0.8613, -714.5),
    :MgNO32 => (0.54, 230.2),
    :MgCl2 => (0.3284, 42.23),
)

# --------------------------------------------------------------------------
# Mutual Deliquescence Relative Humidity (MDRH) for salt mixtures (Table 5)
# Indexed by composition regime identifier
# --------------------------------------------------------------------------
const _ISO2_MDRH = Dict(
    # Sulfate-poor, crustal+Na rich (R₁ > 2, crustal+Na present)
    :A => 0.2,  # Full system with crustal species
    :B => 0.363,  # Sulfate-poor, Na present, no crustal
    :C => 0.54,  # Sulfate-poor, NH4 only (no Na, no crustal)
    # Sulfate-rich (R₁ < 2)
    :D => 0.4,  # NH4HSO4 + LC + (NH4)2SO4
    :E => 0.363,  # Na + NH4 sulfate rich
    :F => 0.2,  # Full system sulfate rich with crustal
)

# =============================================================================
# Part 2: Computational Functions
# =============================================================================

"""
    _iso2_eq_const(K0, A, B, T)

Compute equilibrium constant at temperature T using the van't Hoff equation (Eq. 5).
K(T) = K₀ × exp(A×(T₀/T - 1) + B×(1 + ln(T₀/T) - T₀/T))
"""
function _iso2_eq_const(K0, A, B, T)
    T0T = _ISO2_T0 / T
    return K0 * exp(A * (T0T - 1.0) + B * (1.0 + log(T0T) - T0T))
end
@register_symbolic _iso2_eq_const(K0, A, B, T)

"""
    _iso2_km_gamma(q, z, I)

Compute binary mean activity coefficient using the Kusik-Meissner relationship (Eqs. 9-13).
Returns γ± (not log γ±).
- q: salt-specific interaction parameter
- z: charge product z₁×z₂
- I: ionic strength (mol/kg)
"""
function _iso2_km_gamma(q, z, I)
    I_safe = max(I, 1.0e-20)
    B_param = 0.75 - 0.065 * q                           # Eq. 11
    sI = sqrt(I_safe)
    C_param = I_safe < 6.0 ? 1.0 + 0.055 * q * exp(-0.023 * I_safe^3) : 1.0  # Eq. 13
    log_Gamma_star = -0.5107 * sI / (1.0 + C_param * sI) # Eq. 12
    Gamma0 = max((1.0 + B_param * (1.0 + 0.1 * I_safe)^q - B_param) * 10.0^log_Gamma_star, 1.0e-30)  # Eq. 10
    return Gamma0^z                                       # γ± = (Γ°)^z (Eq. 9)
end
@register_symbolic _iso2_km_gamma(q, z, I)

"""
    _iso2_gamma_T(gamma298, I, T)

Apply temperature correction to activity coefficient (Eq. 14).
Returns γ(T) given γ(298K).
"""
function _iso2_gamma_T(gamma298, I, T)
    TI = T - 273.15  # Eq. 14: T in Celsius
    I_safe = max(I, 1.0e-20)
    log_g298 = log10(max(gamma298, 1.0e-30))
    # A parameter in Eq. 14
    A_corr = -0.41 * sqrt(I_safe) / (1.0 + sqrt(I_safe)) + 0.039 * I_safe^0.92
    CF1 = 1.125 - 0.005 * TI
    CF2 = (0.125 - 0.005 * TI) * A_corr
    log_g_T = CF1 * log_g298 - CF2  # Eq. 14: no z factor
    return 10.0^log_g_T
end
@register_symbolic _iso2_gamma_T(gamma298, I, T)

"""
    _iso2_zsr_m0(table, RH)

Interpolate binary molality from ZSR lookup table at given RH.
Returns molality in mol/kg water.
"""
function _iso2_zsr_m0(table::Vector{Float64}, RH)
    RH_safe = clamp(RH, 0.01, 0.99)
    # Tables are indexed from 1 (aw=0.01) to 100 (aw=1.00)
    idx_real = RH_safe * 100.0
    idx_lo = max(1, min(99, floor(Int, idx_real)))
    idx_hi = idx_lo + 1
    frac = idx_real - idx_lo
    return table[idx_lo] * (1.0 - frac) + table[idx_hi] * frac
end

"""
    _iso2_zsr_water(RH, c_Na, c_NH4, c_SO4, c_HSO4, c_NO3, c_Cl, c_Ca, c_K, c_Mg)

Compute aerosol water content using ZSR correlation (Eq. 16).
Ion concentrations in mol/m³ air. Returns water content in kg/m³ air.
Uses hierarchical ion pairing following ISORROPIA conventions.
"""
function _iso2_zsr_water(RH, c_Na, c_NH4, c_SO4, c_HSO4, c_NO3, c_Cl, c_Ca, c_K, c_Mg)
    RH_safe = clamp(RH, 0.01, 0.99)

    # Working copies of available ions (mol/m³)
    Na_r = max(c_Na, 0.0)
    NH4_r = max(c_NH4, 0.0)
    SO4_r = max(c_SO4, 0.0)
    HSO4_r = max(c_HSO4, 0.0)
    NO3_r = max(c_NO3, 0.0)
    Cl_r = max(c_Cl, 0.0)
    Ca_r = max(c_Ca, 0.0)
    K_r = max(c_K, 0.0)
    Mg_r = max(c_Mg, 0.0)

    W = 0.0

    # --- Crustal cation pairs (prioritize with sulfate, then nitrate, then chloride) ---
    # Ca(NO3)2
    pair = min(Ca_r, NO3_r / 2.0)
    if pair > 1.0e-30
        W += pair / _iso2_zsr_m0(_ISO2_ZSR_CaNO32, RH_safe)
    end
    Ca_r -= pair
    NO3_r -= 2.0 * pair

    # CaCl2
    pair = min(Ca_r, Cl_r / 2.0)
    if pair > 1.0e-30
        W += pair / _iso2_zsr_m0(_ISO2_ZSR_CaCl2, RH_safe)
    end
    Ca_r -= pair
    Cl_r -= 2.0 * pair

    # K2SO4
    pair = min(K_r / 2.0, SO4_r)
    if pair > 1.0e-30
        W += pair / _iso2_zsr_m0(_ISO2_ZSR_K2SO4, RH_safe)
    end
    K_r -= 2.0 * pair
    SO4_r -= pair

    # KNO3
    pair = min(K_r, NO3_r)
    if pair > 1.0e-30
        W += pair / _iso2_zsr_m0(_ISO2_ZSR_KNO3, RH_safe)
    end
    K_r -= pair
    NO3_r -= pair

    # KCl
    pair = min(K_r, Cl_r)
    if pair > 1.0e-30
        W += pair / _iso2_zsr_m0(_ISO2_ZSR_KCl, RH_safe)
    end
    K_r -= pair
    Cl_r -= pair

    # MgSO4
    pair = min(Mg_r, SO4_r)
    if pair > 1.0e-30
        W += pair / _iso2_zsr_m0(_ISO2_ZSR_MgSO4, RH_safe)
    end
    Mg_r -= pair
    SO4_r -= pair

    # Mg(NO3)2
    pair = min(Mg_r, NO3_r / 2.0)
    if pair > 1.0e-30
        W += pair / _iso2_zsr_m0(_ISO2_ZSR_MgNO32, RH_safe)
    end
    Mg_r -= pair
    NO3_r -= 2.0 * pair

    # MgCl2
    pair = min(Mg_r, Cl_r / 2.0)
    if pair > 1.0e-30
        W += pair / _iso2_zsr_m0(_ISO2_ZSR_MgCl2, RH_safe)
    end
    Mg_r -= pair
    Cl_r -= 2.0 * pair

    # --- Sodium pairs ---
    # Na2SO4
    pair = min(Na_r / 2.0, SO4_r)
    if pair > 1.0e-30
        W += pair / _iso2_zsr_m0(_ISO2_ZSR_Na2SO4, RH_safe)
    end
    Na_r -= 2.0 * pair
    SO4_r -= pair

    # NaHSO4
    pair = min(Na_r, HSO4_r)
    if pair > 1.0e-30
        W += pair / _iso2_zsr_m0(_ISO2_ZSR_NaHSO4, RH_safe)
    end
    Na_r -= pair
    HSO4_r -= pair

    # NaNO3
    pair = min(Na_r, NO3_r)
    if pair > 1.0e-30
        W += pair / _iso2_zsr_m0(_ISO2_ZSR_NaNO3, RH_safe)
    end
    Na_r -= pair
    NO3_r -= pair

    # NaCl
    pair = min(Na_r, Cl_r)
    if pair > 1.0e-30
        W += pair / _iso2_zsr_m0(_ISO2_ZSR_NaCl, RH_safe)
    end
    Na_r -= pair
    Cl_r -= pair

    # --- Ammonium pairs ---
    # (NH4)2SO4
    pair = min(NH4_r / 2.0, SO4_r)
    if pair > 1.0e-30
        W += pair / _iso2_zsr_m0(_ISO2_ZSR_NH42SO4, RH_safe)
    end
    NH4_r -= 2.0 * pair
    SO4_r -= pair

    # NH4HSO4
    pair = min(NH4_r, HSO4_r)
    if pair > 1.0e-30
        W += pair / _iso2_zsr_m0(_ISO2_ZSR_NH4HSO4, RH_safe)
    end
    NH4_r -= pair
    HSO4_r -= pair

    # NH4NO3
    pair = min(NH4_r, NO3_r)
    if pair > 1.0e-30
        W += pair / _iso2_zsr_m0(_ISO2_ZSR_NH4NO3, RH_safe)
    end
    NH4_r -= pair
    NO3_r -= pair

    # NH4Cl
    pair = min(NH4_r, Cl_r)
    if pair > 1.0e-30
        W += pair / _iso2_zsr_m0(_ISO2_ZSR_NH4Cl, RH_safe)
    end
    NH4_r -= pair
    Cl_r -= pair

    # --- Remaining acid ---
    # H2SO4 (remaining SO4)
    if SO4_r > 1.0e-30
        W += SO4_r / _iso2_zsr_m0(_ISO2_ZSR_H2SO4, RH_safe)
    end

    # H-HSO4 (remaining HSO4)
    if HSO4_r > 1.0e-30
        W += HSO4_r / _iso2_zsr_m0(_ISO2_ZSR_H2SO4, RH_safe)
    end

    return max(W, 1.0e-20)
end
@register_symbolic _iso2_zsr_water(RH, c_Na, c_NH4, c_SO4, c_HSO4, c_NO3, c_Cl, c_Ca, c_K, c_Mg)

"""
    _iso2_fb_ncp(a, b)

Fischer-Burmeister NCP (nonlinear complementarity problem) function.
Returns √(a² + b² + ε) - a - b, which equals zero when either a=0 or b=0 (or both).
Used for solid precipitation equilibria: either solid amount is zero, or solution is saturated.
"""
function _iso2_fb_ncp(a, b)
    return sqrt(a^2 + b^2 + 1.0e-20) - a - b
end
@register_symbolic _iso2_fb_ncp(a, b)

"""
    _iso2_smooth_min(a, b)

Smooth differentiable approximation to min(a, b).
Uses the identity: min(a,b) = 0.5*(a + b - √((a-b)² + ε)).
"""
function _iso2_smooth_min(a, b)
    return 0.5 * (a + b - sqrt((a - b)^2 + 1.0e-20))
end
@register_symbolic _iso2_smooth_min(a, b)

"""
    _iso2_drh_T(drh0, coeff, T)

Compute temperature-dependent deliquescence relative humidity (Eq. 17, Table 4).
DRH(T) = DRH(T₀) × exp(coeff × (1/T - 1/T₀))
Result is clamped to [0, 1].
"""
function _iso2_drh_T(drh0, coeff, T)
    return clamp(drh0 * exp(coeff * (1.0 / T - 1.0 / _ISO2_T0)), 0.0, 1.0)
end
@register_symbolic _iso2_drh_T(drh0, coeff, T)

"""
    _iso2_smooth_step(x)

Smooth logistic step function: σ(x) = 1/(1 + exp(-clamp(x, -20, 20))).
Returns ≈0 for x ≪ 0 and ≈1 for x ≫ 0, with smooth transition near x=0.
"""
function _iso2_smooth_step(x)
    x_safe = clamp(x, -20.0, 20.0)
    return 1.0 / (1.0 + exp(-x_safe))
end
@register_symbolic _iso2_smooth_step(x)

"""
    _iso2_compute_mdrh(Na, SO4, NH3, NO3, Cl, Ca, K, Mg)

Compute the Mutual Deliquescence Relative Humidity (MDRH) from total composition
(all in mol/m³). Uses sulfate ratio R₁ and presence of crustal/Na species to
select the appropriate MDRH from Table 5.
"""
function _iso2_compute_mdrh(Na, SO4, NH3, NO3, Cl, Ca, K, Mg)
    # Sulfate ratio R₁ = (NH3 + Na + 2Ca + K + 2Mg) / SO4
    SO4_safe = max(SO4, 1.0e-30)
    R1 = (NH3 + Na + 2.0 * Ca + K + 2.0 * Mg) / SO4_safe

    has_crustal = (Ca + K + Mg) > 1.0e-20
    has_Na = Na > 1.0e-20

    return if R1 >= 2.0
        # Sulfate-poor regime
        if has_crustal
            return 0.2  # Table 5: full system with crustal
        elseif has_Na
            return 0.363  # Table 5: Na present, no crustal
        else
            return 0.54  # Table 5: NH4-only system
        end
    else
        # Sulfate-rich regime
        if has_crustal
            return 0.2  # Table 5: sulfate-rich with crustal
        elseif has_Na
            return 0.363  # Table 5: Na + NH4 sulfate rich
        else
            return 0.4  # Table 5: NH4HSO4 + LC + (NH4)2SO4
        end
    end
end
@register_symbolic _iso2_compute_mdrh(Na, SO4, NH3, NO3, Cl, Ca, K, Mg)

# =============================================================================
# Part 3: ModelingToolkit Component
# =============================================================================

"""
    IsorropiaEquilibrium(; name=:IsorropiaEquilibrium)

ISORROPIA II inorganic aerosol thermodynamic equilibrium model for the
K⁺–Ca²⁺–Mg²⁺–NH₄⁺–Na⁺–SO₄²⁻–NO₃⁻–Cl⁻–H₂O system.

Supports both **metastable** (no solids) and **stable** (with solid precipitation)
solution modes, controlled by the `stable` parameter (0=metastable, 1=stable).

The model solves 51 coupled algebraic equations:
- 8 mass balance equations (with solid stoichiometric contributions in stable mode)
- 1 charge balance (electroneutrality)
- 5 aqueous equilibrium expressions (HSO₄⁻, NH₃, HNO₃, HCl, H₂O)
- 1 ZSR water content equation with MDRH cutoff (Eqs. 16-18, 22)
- 1 ionic strength definition
- 16 activity coefficient equations (Kusik-Meissner, Eqs. 9-13)
- 19 solid precipitation equations (Fischer-Burmeister NCP complementarity)

In stable mode, solid salts precipitate when the solution becomes supersaturated.
Each solid equilibrium uses a smooth complementarity formulation:
either the solid amount is zero, or the solution is exactly at saturation.

**Deliquescence**: Temperature-dependent DRH values are computed for each salt
using Eq. 17 and Table 4. The mutual deliquescence relative humidity (MDRH)
from Table 5 is enforced in stable mode: below the MDRH, the aerosol water
content is driven to zero, producing a completely dry (crystalline) aerosol.
In metastable mode, the MDRH cutoff is not applied.

Activity coefficients are computed using the Kusik-Meissner binary relationship
(Eqs. 9-13) with temperature correction (Eq. 14).

**Reference**: Fountoukis, C. and Nenes, A.: ISORROPIA II: a computationally
efficient thermodynamic equilibrium model for
K⁺–Ca²⁺–Mg²⁺–NH₄⁺–Na⁺–SO₄²⁻–NO₃⁻–Cl⁻–H₂O aerosols,
Atmos. Chem. Phys., 7, 4639–4659, 2007.
"""
@component function IsorropiaEquilibrium(; name = :IsorropiaEquilibrium)

    # ------------------------------------------------------------------
    # Constants
    # ------------------------------------------------------------------
    @constants begin
        # Physical constants
        p_ref = 101325.0, [description = "Reference pressure (1 atm)", unit = u"Pa"]
        R_gas_si = 8.314462, [description = "Gas constant (SI)", unit = u"Pa*m^3/mol/K"]

        # Thermodynamic reference temperature (Table 2)
        T_0 = 298.15, [description = "Reference temperature for thermodynamic data", unit = u"K"]

        # Reference molality for unit conversions
        m_ref = 1.0, [description = "Reference molality for dimensionless ratios", unit = u"mol/kg"]

        # Equilibrium constant parameters (K0, A, B) from Table 2, Eq. 5
        # Van't Hoff: K(T) = K₀ × exp[A×(T₀/T - 1) + B×(1 + ln(T₀/T) - T₀/T)]
        # where A = ΔH°/(R·T₀) and B = Δcₚ°/R are dimensionless
        # HSO4⁻ ↔ H⁺ + SO4²⁻ (Table 2: mol/kg)
        K1_0 = 1.015e-2, [description = "HSO4 dissociation equilibrium constant at T0", unit = u"mol/kg"]
        K1_A = 8.85, [description = "HSO4 dissociation enthalpy parameter (dimensionless)"]
        K1_B = 25.14, [description = "HSO4 dissociation heat capacity parameter (dimensionless)"]

        # NH3(g) ↔ NH3(aq) (Table 2: mol/(kg·atm) converted to SI)
        # Standard atmosphere = 101325 Pa, so K21_SI = K21_atm / 101325
        K21_0_atm = 5.764e1, [description = "NH3 dissolution equilibrium constant at T0 in atm units", unit = u"mol/kg/Constants.atm"]
        K21_0 = K21_0_atm * 9.86923e-6, [description = "NH3 dissolution equilibrium constant at T0 in SI units", unit = u"mol/(kg*Pa)"]
        K21_A = 13.79, [description = "NH3 dissolution enthalpy parameter (dimensionless)"]
        K21_B = -5.393, [description = "NH3 dissolution heat capacity parameter (dimensionless)"]

        # NH3(aq) ↔ NH4⁺ + OH⁻ (Table 2: mol/kg)
        K22_0 = 1.805e-5, [description = "NH3 dissociation equilibrium constant at T0", unit = u"mol/kg"]
        K22_A = -1.5, [description = "NH3 dissociation enthalpy parameter (dimensionless)"]
        K22_B = 26.92, [description = "NH3 dissociation heat capacity parameter (dimensionless)"]

        # HCl(g) ↔ H⁺ + Cl⁻ (Table 2: mol²/(kg²·atm) converted to SI)
        # K3_SI = K3_atm / (1 atm) = K3_atm / 101325 Pa
        K3_0_atm = 1.971e6, [description = "HCl gas-liquid equilibrium constant at T0 in atm units", unit = u"mol^2/kg^2/Constants.atm"]
        K3_0 = K3_0_atm * 9.86923e-6, [description = "HCl gas-liquid equilibrium constant at T0 in SI units", unit = u"mol^2/(kg^2*Pa)"]
        K3_A = 30.2, [description = "HCl gas-liquid enthalpy parameter (dimensionless)"]
        K3_B = 19.91, [description = "HCl gas-liquid heat capacity parameter (dimensionless)"]

        # HNO3(g) ↔ H⁺ + NO3⁻ (Table 2: mol²/(kg²·atm) converted to SI)
        # K4_SI = K4_atm / (1 atm) = K4_atm / 101325 Pa
        K4_0_atm = 2.511e6, [description = "HNO3 gas-liquid equilibrium constant at T0 in atm units", unit = u"mol^2/kg^2/Constants.atm"]
        K4_0 = K4_0_atm * 9.86923e-6, [description = "HNO3 gas-liquid equilibrium constant at T0 in SI units", unit = u"mol^2/(kg^2*Pa)"]
        K4_A = 29.17, [description = "HNO3 gas-liquid enthalpy parameter (dimensionless)"]
        K4_B = 16.83, [description = "HNO3 gas-liquid heat capacity parameter (dimensionless)"]

        # H2O ↔ H⁺ + OH⁻ (Table 2: mol²/kg²)
        Kw_0 = 1.01e-14, [description = "Water dissociation equilibrium constant at T0", unit = u"mol^2/kg^2"]
        Kw_A = -22.52, [description = "Water dissociation enthalpy parameter (dimensionless)"]
        Kw_B = 26.92, [description = "Water dissociation heat capacity parameter (dimensionless)"]

        # --- Solid dissolution equilibrium constants (Table 2) ---
        # Format: K0 [units at 298.15K], A = ΔH°/(R·T₀), B = Δcₚ°/R
        # All K0 values converted to SI units where needed

        # (NH4)2SO4(s) ↔ 2NH4+ + SO4²⁻ (K5, Table 2: mol³/kg³)
        K5_0 = 1.817, [description = "(NH4)2SO4 dissolution K0", unit = u"mol^3/kg^3"]
        K5_A = -2.65, [description = "(NH4)2SO4 dissolution A (dimensionless)"]
        K5_B = 38.57, [description = "(NH4)2SO4 dissolution B (dimensionless)"]

        # NH4HSO4(s) ↔ NH4+ + H+ + SO4²⁻ (K6, Table 2: mol³/kg³)
        K6_0 = 1.383e2, [description = "NH4HSO4 dissolution K0", unit = u"mol^3/kg^3"]
        K6_A = -2.87, [description = "NH4HSO4 dissolution A (dimensionless)"]
        K6_B = 15.83, [description = "NH4HSO4 dissolution B (dimensionless)"]

        # (NH4)3H(SO4)2(s) ↔ 3NH4+ + HSO4⁻ + SO4²⁻ (K7 letovicite, Table 2: mol⁵/kg⁵)
        K7_0 = 2.93e1, [description = "Letovicite dissolution K0", unit = u"mol^5/kg^5"]
        K7_A = -5.19, [description = "Letovicite dissolution A (dimensionless)"]
        K7_B = 54.4, [description = "Letovicite dissolution B (dimensionless)"]

        # NH4NO3(s) ↔ NH3(g) + HNO3(g) (K8, Table 2: atm² converted to Pa²)
        # K8_SI = K8_atm² × (101325 Pa/atm)² = K8_atm² × 1.0269e10
        K8_0_atm2 = 4.3559e-17, [description = "NH4NO3 dissolution K0 in atm²", unit = u"Constants.atm^2"]
        K8_0 = K8_0_atm2 * 1.0269e10, [description = "NH4NO3 dissolution K0 in SI", unit = u"Pa^2"]
        K8_A = -74.38, [description = "NH4NO3 dissolution A (dimensionless)"]
        K8_B = 6.12, [description = "NH4NO3 dissolution B (dimensionless)"]

        # NH4Cl(s) ↔ NH3(g) + HCl(g) (K9, Table 2: atm² converted to Pa²)
        # K9_SI = K9_atm² × (101325 Pa/atm)² = K9_atm² × 1.0269e10
        K9_0_atm2 = 1.086e-16, [description = "NH4Cl dissolution K0 in atm²", unit = u"Constants.atm^2"]
        K9_0 = K9_0_atm2 * 1.0269e10, [description = "NH4Cl dissolution K0 in SI", unit = u"Pa^2"]
        K9_A = -71.0, [description = "NH4Cl dissolution A (dimensionless)"]
        K9_B = 2.4, [description = "NH4Cl dissolution B (dimensionless)"]

        # Na2SO4(s) ↔ 2Na+ + SO4²⁻ (K10, Table 2: mol³/kg³)
        K10_0 = 4.799e-1, [description = "Na2SO4 dissolution K0", unit = u"mol^3/kg^3"]
        K10_A = 0.98, [description = "Na2SO4 dissolution A (dimensionless)"]
        K10_B = 39.75, [description = "Na2SO4 dissolution B (dimensionless)"]

        # NaHSO4(s) ↔ Na+ + H+ + SO4²⁻ (K11, Table 2: mol³/kg³)
        K11_0 = 2.413e4, [description = "NaHSO4 dissolution K0", unit = u"mol^3/kg^3"]
        K11_A = 0.79, [description = "NaHSO4 dissolution A (dimensionless)"]
        K11_B = 14.75, [description = "NaHSO4 dissolution B (dimensionless)"]

        # NaNO3(s) ↔ Na+ + NO3⁻ (K12, Table 2: mol²/kg²)
        K12_0 = 1.197e1, [description = "NaNO3 dissolution K0", unit = u"mol^2/kg^2"]
        K12_A = -8.22, [description = "NaNO3 dissolution A (dimensionless)"]
        K12_B = 16.01, [description = "NaNO3 dissolution B (dimensionless)"]

        # NaCl(s) ↔ Na+ + Cl⁻ (K13, Table 2: mol²/kg²)
        K13_0 = 3.766e1, [description = "NaCl dissolution K0", unit = u"mol^2/kg^2"]
        K13_A = -1.56, [description = "NaCl dissolution A (dimensionless)"]
        K13_B = 16.9, [description = "NaCl dissolution B (dimensionless)"]

        # CaSO4(s) ↔ Ca²+ + SO4²⁻ (K14, Table 2: mol²/kg²)
        K14_0 = 4.319e-5, [description = "CaSO4 dissolution K0", unit = u"mol^2/kg^2"]
        K14_A = 0.0, [description = "CaSO4 dissolution A (dimensionless)"]
        K14_B = 0.0, [description = "CaSO4 dissolution B (dimensionless)"]

        # Ca(NO3)2(s) ↔ Ca²+ + 2NO3⁻ (K15, Table 2: mol³/kg³)
        K15_0 = 6.067e5, [description = "Ca(NO3)2 dissolution K0", unit = u"mol^3/kg^3"]
        K15_A = -11.299, [description = "Ca(NO3)2 dissolution A (dimensionless)"]
        K15_B = 0.0, [description = "Ca(NO3)2 dissolution B (dimensionless)"]

        # CaCl2(s) ↔ Ca²+ + 2Cl⁻ (K16, Table 2: mol³/kg³)
        K16_0 = 7.974e11, [description = "CaCl2 dissolution K0", unit = u"mol^3/kg^3"]
        K16_A = -8.698, [description = "CaCl2 dissolution A (dimensionless)"]
        K16_B = 0.0, [description = "CaCl2 dissolution B (dimensionless)"]

        # K2SO4(s) ↔ 2K+ + SO4²⁻ (K17, Table 2: mol³/kg³)
        K17_0 = 1.569e-2, [description = "K2SO4 dissolution K0", unit = u"mol^3/kg^3"]
        K17_A = -9.585, [description = "K2SO4 dissolution A (dimensionless)"]
        K17_B = 45.81, [description = "K2SO4 dissolution B (dimensionless)"]

        # KHSO4(s) ↔ K+ + H+ + SO4²⁻ (K18, Table 2: mol³/kg³)
        K18_0 = 24.016, [description = "KHSO4 dissolution K0", unit = u"mol^3/kg^3"]
        K18_A = -8.423, [description = "KHSO4 dissolution A (dimensionless)"]
        K18_B = 17.96, [description = "KHSO4 dissolution B (dimensionless)"]

        # KNO3(s) ↔ K+ + NO3⁻ (K19, Table 2: mol²/kg²)
        K19_0 = 0.872, [description = "KNO3 dissolution K0", unit = u"mol^2/kg^2"]
        K19_A = -14.08, [description = "KNO3 dissolution A (dimensionless)"]
        K19_B = 19.39, [description = "KNO3 dissolution B (dimensionless)"]

        # KCl(s) ↔ K+ + Cl⁻ (K20, Table 2: mol²/kg²)
        K20_0 = 8.68, [description = "KCl dissolution K0", unit = u"mol^2/kg^2"]
        K20_A = -6.902, [description = "KCl dissolution A (dimensionless)"]
        K20_B = 19.95, [description = "KCl dissolution B (dimensionless)"]

        # MgSO4(s) ↔ Mg²+ + SO4²⁻ (K21s, Table 2: mol²/kg²)
        # Note: using K21s to avoid conflict with K21 (NH3 dissolution)
        K21s_0 = 1.079e5, [description = "MgSO4 dissolution K0", unit = u"mol^2/kg^2"]
        K21s_A = 36.798, [description = "MgSO4 dissolution A (dimensionless)"]
        K21s_B = 0.0, [description = "MgSO4 dissolution B (dimensionless)"]

        # Mg(NO3)2(s) ↔ Mg²+ + 2NO3⁻ (K22s, Table 2: mol³/kg³)
        K22s_0 = 2.507e15, [description = "Mg(NO3)2 dissolution K0", unit = u"mol^3/kg^3"]
        K22s_A = -8.754, [description = "Mg(NO3)2 dissolution A (dimensionless)"]
        K22s_B = 0.0, [description = "Mg(NO3)2 dissolution B (dimensionless)"]

        # MgCl2(s) ↔ Mg²+ + 2Cl⁻ (K23, Table 2: mol³/kg³)
        K23_0 = 9.557e21, [description = "MgCl2 dissolution K0", unit = u"mol^3/kg^3"]
        K23_A = -1.347, [description = "MgCl2 dissolution A (dimensionless)"]
        K23_B = 0.0, [description = "MgCl2 dissolution B (dimensionless)"]
    end

    # ------------------------------------------------------------------
    # Parameters (inputs)
    # ------------------------------------------------------------------
    @parameters begin
        T_env = 298.15, [description = "Ambient temperature", unit = u"K"]
        RH = 0.8, [description = "Ambient relative humidity (dimensionless)"]
        stable = 0, [description = "Solution mode: 0=metastable (no solids), 1=stable (with solids) (dimensionless)"]
        W_Na_total = 0.0, [description = "Total sodium concentration", unit = u"mol/m^3"]
        W_SO4_total = 1.0e-7, [description = "Total sulfate concentration", unit = u"mol/m^3"]
        W_NH3_total = 2.0e-7, [description = "Total ammonia + ammonium concentration", unit = u"mol/m^3"]
        W_NO3_total = 1.0e-7, [description = "Total nitrate + nitric acid concentration", unit = u"mol/m^3"]
        W_Cl_total = 0.0, [description = "Total chloride + HCl concentration", unit = u"mol/m^3"]
        W_Ca_total = 0.0, [description = "Total calcium concentration", unit = u"mol/m^3"]
        W_K_total = 0.0, [description = "Total potassium concentration", unit = u"mol/m^3"]
        W_Mg_total = 0.0, [description = "Total magnesium concentration", unit = u"mol/m^3"]
    end

    # ------------------------------------------------------------------
    # Variables (unknowns)
    # ------------------------------------------------------------------
    @variables begin
        # Aqueous ion concentrations (mol/m³ air)
        c_H(t), [description = "Aqueous H⁺ concentration", unit = u"mol/m^3", guess = 1.0e-10]
        c_Na(t), [description = "Aqueous Na⁺ concentration", unit = u"mol/m^3", guess = 1.0e-10]
        c_NH4(t), [description = "Aqueous NH₄⁺ concentration", unit = u"mol/m^3", guess = 1.0e-7]
        c_Cl(t), [description = "Aqueous Cl⁻ concentration", unit = u"mol/m^3", guess = 1.0e-10]
        c_SO4(t), [description = "Aqueous SO₄²⁻ concentration", unit = u"mol/m^3", guess = 1.0e-7]
        c_HSO4(t), [description = "Aqueous HSO₄⁻ concentration", unit = u"mol/m^3", guess = 1.0e-10]
        c_NO3(t), [description = "Aqueous NO₃⁻ concentration", unit = u"mol/m^3", guess = 1.0e-10]
        c_Ca(t), [description = "Aqueous Ca²⁺ concentration", unit = u"mol/m^3", guess = 1.0e-10]
        c_K(t), [description = "Aqueous K⁺ concentration", unit = u"mol/m^3", guess = 1.0e-10]
        c_Mg(t), [description = "Aqueous Mg²⁺ concentration", unit = u"mol/m^3", guess = 1.0e-10]
        c_OH(t), [description = "Aqueous OH⁻ concentration", unit = u"mol/m^3", guess = 1.0e-10]

        # Gas-phase concentrations (mol/m³ air)
        g_NH3(t), [description = "Gas-phase NH₃ concentration", unit = u"mol/m^3", guess = 1.0e-8]
        g_HNO3(t), [description = "Gas-phase HNO₃ concentration", unit = u"mol/m^3", guess = 1.0e-8]
        g_HCl(t), [description = "Gas-phase HCl concentration", unit = u"mol/m^3", guess = 1.0e-10]

        # Aerosol water content (kg/m³ air)
        W_w(t), [description = "Aerosol liquid water content", unit = u"kg/m^3", guess = 1.0e-6]

        # Ionic strength (mol/kg) - intermediate variable
        I_s(t), [description = "Ionic strength of the solution", unit = u"mol/kg", guess = 1.0]

        # Activity coefficient products (dimensionless) - intermediate variables
        γ_H2SO4(t), [description = "Mean activity coeff. H₂SO₄ (dimensionless)", guess = 1.0]
        γ_HHSO4(t), [description = "Mean activity coeff. H-HSO₄ (dimensionless)", guess = 1.0]
        γ_HNO3(t), [description = "Mean activity coeff. HNO₃ (dimensionless)", guess = 1.0]
        γ_HCl(t), [description = "Mean activity coeff. HCl (dimensionless)", guess = 1.0]
        γ_NH4Cl(t), [description = "Mean activity coeff. NH₄Cl (dimensionless)", guess = 1.0]

        # Additional activity coefficients for solid dissolution equilibria
        γ_NaCl(t), [description = "Mean activity coeff. NaCl (dimensionless)", guess = 1.0]
        γ_NaNO3(t), [description = "Mean activity coeff. NaNO₃ (dimensionless)", guess = 1.0]
        γ_Na2SO4(t), [description = "Mean activity coeff. Na₂SO₄ (dimensionless)", guess = 1.0]
        γ_NH42SO4(t), [description = "Mean activity coeff. (NH₄)₂SO₄ (dimensionless)", guess = 1.0]
        γ_CaNO32(t), [description = "Mean activity coeff. Ca(NO₃)₂ (dimensionless)", guess = 1.0]
        γ_CaCl2(t), [description = "Mean activity coeff. CaCl₂ (dimensionless)", guess = 1.0]
        γ_K2SO4(t), [description = "Mean activity coeff. K₂SO₄ (dimensionless)", guess = 1.0]
        γ_KNO3(t), [description = "Mean activity coeff. KNO₃ (dimensionless)", guess = 1.0]
        γ_KCl(t), [description = "Mean activity coeff. KCl (dimensionless)", guess = 1.0]
        γ_MgNO32(t), [description = "Mean activity coeff. Mg(NO₃)₂ (dimensionless)", guess = 1.0]
        γ_MgCl2(t), [description = "Mean activity coeff. MgCl₂ (dimensionless)", guess = 1.0]

        # Solid species concentrations (mol/m³ air) — zero in metastable mode
        s_CaSO4(t), [description = "Solid CaSO₄ concentration", unit = u"mol/m^3", guess = 0.0]
        s_CaNO32(t), [description = "Solid Ca(NO₃)₂ concentration", unit = u"mol/m^3", guess = 0.0]
        s_CaCl2(t), [description = "Solid CaCl₂ concentration", unit = u"mol/m^3", guess = 0.0]
        s_K2SO4(t), [description = "Solid K₂SO₄ concentration", unit = u"mol/m^3", guess = 0.0]
        s_KHSO4(t), [description = "Solid KHSO₄ concentration", unit = u"mol/m^3", guess = 0.0]
        s_KNO3(t), [description = "Solid KNO₃ concentration", unit = u"mol/m^3", guess = 0.0]
        s_KCl(t), [description = "Solid KCl concentration", unit = u"mol/m^3", guess = 0.0]
        s_MgSO4(t), [description = "Solid MgSO₄ concentration", unit = u"mol/m^3", guess = 0.0]
        s_MgNO32(t), [description = "Solid Mg(NO₃)₂ concentration", unit = u"mol/m^3", guess = 0.0]
        s_MgCl2(t), [description = "Solid MgCl₂ concentration", unit = u"mol/m^3", guess = 0.0]
        s_Na2SO4(t), [description = "Solid Na₂SO₄ concentration", unit = u"mol/m^3", guess = 0.0]
        s_NaHSO4(t), [description = "Solid NaHSO₄ concentration", unit = u"mol/m^3", guess = 0.0]
        s_NaNO3(t), [description = "Solid NaNO₃ concentration", unit = u"mol/m^3", guess = 0.0]
        s_NaCl(t), [description = "Solid NaCl concentration", unit = u"mol/m^3", guess = 0.0]
        s_NH42SO4(t), [description = "Solid (NH₄)₂SO₄ concentration", unit = u"mol/m^3", guess = 0.0]
        s_NH4HSO4(t), [description = "Solid NH₄HSO₄ concentration", unit = u"mol/m^3", guess = 0.0]
        s_LC(t), [description = "Solid letovicite (NH₄)₃H(SO₄)₂ concentration", unit = u"mol/m^3", guess = 0.0]
        s_NH4NO3(t), [description = "Solid NH₄NO₃ concentration", unit = u"mol/m^3", guess = 0.0]
        s_NH4Cl(t), [description = "Solid NH₄Cl concentration", unit = u"mol/m^3", guess = 0.0]
    end

    # ------------------------------------------------------------------
    # Intermediate computations
    # ------------------------------------------------------------------

    # Aqueous equilibrium constants at ambient temperature (Eq. 5, Table 2)
    # Use proper units throughout per project standards
    K1 = _iso2_eq_const(K1_0, K1_A, K1_B, T_env)
    K21 = _iso2_eq_const(K21_0, K21_A, K21_B, T_env)
    K22 = _iso2_eq_const(K22_0, K22_A, K22_B, T_env)
    K3 = _iso2_eq_const(K3_0, K3_A, K3_B, T_env)
    K4 = _iso2_eq_const(K4_0, K4_A, K4_B, T_env)
    Kw = _iso2_eq_const(Kw_0, Kw_A, Kw_B, T_env)
    K2 = K21 * K22

    # Solid dissolution equilibrium constants at ambient temperature
    Ksp5 = _iso2_eq_const(K5_0, K5_A, K5_B, T_env)       # (NH4)2SO4 [mol³/kg³]
    Ksp6 = _iso2_eq_const(K6_0, K6_A, K6_B, T_env)       # NH4HSO4 [mol³/kg³]
    Ksp7 = _iso2_eq_const(K7_0, K7_A, K7_B, T_env)       # Letovicite [mol⁵/kg⁵]
    Ksp8 = _iso2_eq_const(K8_0, K8_A, K8_B, T_env)       # NH4NO3 [Pa²]
    Ksp9 = _iso2_eq_const(K9_0, K9_A, K9_B, T_env)       # NH4Cl [Pa²]
    Ksp10 = _iso2_eq_const(K10_0, K10_A, K10_B, T_env) / m_ref^3   # Na2SO4 [mol³/kg³]
    Ksp11 = _iso2_eq_const(K11_0, K11_A, K11_B, T_env) / m_ref^3   # NaHSO4 [mol³/kg³]
    Ksp12 = _iso2_eq_const(K12_0, K12_A, K12_B, T_env) / m_ref^2   # NaNO3 [mol²/kg²]
    Ksp13 = _iso2_eq_const(K13_0, K13_A, K13_B, T_env) / m_ref^2   # NaCl [mol²/kg²]
    Ksp14 = _iso2_eq_const(K14_0, K14_A, K14_B, T_env) / m_ref^2   # CaSO4 [mol²/kg²]
    Ksp15 = _iso2_eq_const(K15_0, K15_A, K15_B, T_env) / m_ref^3   # Ca(NO3)2 [mol³/kg³]
    Ksp16 = _iso2_eq_const(K16_0, K16_A, K16_B, T_env) / m_ref^3   # CaCl2 [mol³/kg³]
    Ksp17 = _iso2_eq_const(K17_0, K17_A, K17_B, T_env) / m_ref^3   # K2SO4 [mol³/kg³]
    Ksp18 = _iso2_eq_const(K18_0, K18_A, K18_B, T_env) / m_ref^3   # KHSO4 [mol³/kg³]
    Ksp19 = _iso2_eq_const(K19_0, K19_A, K19_B, T_env) / m_ref^2   # KNO3 [mol²/kg²]
    Ksp20 = _iso2_eq_const(K20_0, K20_A, K20_B, T_env) / m_ref^2   # KCl [mol²/kg²]
    Ksp21s = _iso2_eq_const(K21s_0, K21s_A, K21s_B, T_env) / m_ref^2 # MgSO4 [mol²/kg²]
    Ksp22s = _iso2_eq_const(K22s_0, K22s_A, K22s_B, T_env) / m_ref^3 # Mg(NO3)2 [mol³/kg³]
    Ksp23 = _iso2_eq_const(K23_0, K23_A, K23_B, T_env) / m_ref^3   # MgCl2 [mol³/kg³]

    # Molalities (mol/kg) - proper units maintained per project standards
    m_H = c_H / W_w
    m_Na = c_Na / W_w
    m_NH4 = c_NH4 / W_w
    m_SO4 = c_SO4 / W_w
    m_HSO4 = c_HSO4 / W_w
    m_NO3 = c_NO3 / W_w
    m_Cl = c_Cl / W_w
    m_Ca = c_Ca / W_w
    m_K = c_K / W_w
    m_Mg = c_Mg / W_w
    m_OH = c_OH / W_w

    # Partial pressures (Pa) - from ideal gas law: p = c * R * T
    p_NH3 = g_NH3 * R_gas_si * T_env
    p_HNO3 = g_HNO3 * R_gas_si * T_env
    p_HCl = g_HCl * R_gas_si * T_env

    # Compute MDRH from composition (Table 5) — uses total inputs
    MDRH = _iso2_compute_mdrh(
        W_Na_total, W_SO4_total, W_NH3_total,
        W_NO3_total, W_Cl_total,
        W_Ca_total, W_K_total, W_Mg_total
    )

    # MDRH factor: smooth transition from dry (0) to wet (1) around MDRH
    # Scale factor 100 gives transition width ~0.02 in RH
    mdrh_factor = _iso2_smooth_step(100.0 * (RH - MDRH))

    # Derived activity coefficients (from cross-relations)
    # γ_KHSO4 = √(γ_HHSO4 × γ_KCl / γ_HCl)
    γ_KHSO4_derived = sqrt(γ_HHSO4 * γ_KCl / γ_HCl)
    # γ_NaHSO4 = √(γ_HHSO4 × γ_NaCl / γ_HCl)
    γ_NaHSO4_derived = sqrt(γ_HHSO4 * γ_NaCl / γ_HCl)
    # γ_NH4HSO4 = √(γ_HHSO4 × γ_NH4Cl / γ_HCl)
    γ_NH4HSO4_derived = sqrt(γ_HHSO4 * γ_NH4Cl / γ_HCl)
    # γ_LC = (γ_NH42SO4³ × γ_NH4HSO4)^0.2
    γ_LC_derived = (γ_NH42SO4^3 * γ_NH4HSO4_derived)^0.2
    # γ_MgSO4 = γ_Na2SO4 (approximate, from ISORROPIA convention)
    γ_MgSO4_derived = γ_Na2SO4

    # Saturation ratios for each solid (IAP / Ksp)
    # Aqueous-phase dissolution: IAP = product of molalities × activity coefficients
    SR_NH42SO4 = m_NH4^2 * m_SO4 * γ_NH42SO4^3 / Ksp5
    SR_NH4HSO4 = m_NH4 * m_H * m_SO4 * γ_NH4HSO4_derived^3 / Ksp6
    SR_LC = m_NH4^3 * m_HSO4 * m_SO4 * γ_LC_derived^5 / Ksp7
    SR_Na2SO4 = m_Na^2 * m_SO4 * γ_Na2SO4^3 / Ksp10
    SR_NaHSO4 = m_Na * m_H * m_SO4 * γ_NaHSO4_derived^3 / Ksp11
    SR_NaNO3 = m_Na * m_NO3 * γ_NaNO3^2 / Ksp12
    SR_NaCl = m_Na * m_Cl * γ_NaCl^2 / Ksp13
    SR_CaSO4 = m_Ca * m_SO4 * γ_CaNO32^2 / Ksp14  # Use γ_CaNO32 as proxy for Ca activity
    SR_CaNO32 = m_Ca * m_NO3^2 * γ_CaNO32^3 / Ksp15
    SR_CaCl2 = m_Ca * m_Cl^2 * γ_CaCl2^3 / Ksp16
    SR_K2SO4 = m_K^2 * m_SO4 * γ_K2SO4^3 / Ksp17
    SR_KHSO4 = m_K * m_H * m_SO4 * γ_KHSO4_derived^3 / Ksp18
    SR_KNO3 = m_K * m_NO3 * γ_KNO3^2 / Ksp19
    SR_KCl = m_K * m_Cl * γ_KCl^2 / Ksp20
    SR_MgSO4 = m_Mg * m_SO4 * γ_MgSO4_derived^2 / Ksp21s
    SR_MgNO32 = m_Mg * m_NO3^2 * γ_MgNO32^3 / Ksp22s
    SR_MgCl2 = m_Mg * m_Cl^2 * γ_MgCl2^3 / Ksp23
    # Gas-phase equilibria
    SR_NH4NO3 = p_NH3 * p_HNO3 / Ksp8
    SR_NH4Cl = p_NH3 * p_HCl / Ksp9

    # ------------------------------------------------------------------
    # Equations (50 total)
    # ------------------------------------------------------------------
    eqs = [
        # === Mass balances with solid stoichiometric contributions (Eqs. 1-8) ===
        # 1. Sodium: c_Na + (solid terms only if stable mode) = W_Na_total
        0 ~ (c_Na + stable * (2 * s_Na2SO4 + s_NaHSO4 + s_NaNO3 + s_NaCl)) -
            W_Na_total,

        # 2. Sulfate: c_SO4 + c_HSO4 + (solid terms only if stable mode) = W_SO4_total
        0 ~ (
            c_SO4 + c_HSO4 + stable * (s_CaSO4 + s_K2SO4 + s_KHSO4 + s_Na2SO4 + s_NaHSO4 +
                s_NH42SO4 + s_NH4HSO4 + 2 * s_LC + s_MgSO4)
        ) -
            W_SO4_total,

        # 3. Ammonia/ammonium: c_NH4 + g_NH3 + (solid terms only if stable mode) = W_NH3_total
        0 ~ (
            c_NH4 + g_NH3 + stable * (2 * s_NH42SO4 + s_NH4HSO4 + 3 * s_LC +
                s_NH4NO3 + s_NH4Cl)
        ) -
            W_NH3_total,

        # 4. Nitrate: c_NO3 + g_HNO3 + (solid terms only if stable mode) = W_NO3_total
        0 ~ (
            c_NO3 + g_HNO3 + stable * (2 * s_CaNO32 + s_KNO3 + s_NaNO3 +
                2 * s_MgNO32 + s_NH4NO3)
        ) -
            W_NO3_total,

        # 5. Chloride: c_Cl + g_HCl + (solid terms only if stable mode) = W_Cl_total
        0 ~ (
            c_Cl + g_HCl + stable * (2 * s_CaCl2 + s_KCl + s_NaCl +
                2 * s_MgCl2 + s_NH4Cl)
        ) -
            W_Cl_total,

        # 6. Calcium: c_Ca + (solid terms only if stable mode) = W_Ca_total
        0 ~ (c_Ca + stable * (s_CaSO4 + s_CaNO32 + s_CaCl2)) -
            W_Ca_total,

        # 7. Potassium: c_K + (solid terms only if stable mode) = W_K_total
        0 ~ (c_K + stable * (2 * s_K2SO4 + s_KHSO4 + s_KNO3 + s_KCl)) -
            W_K_total,

        # 8. Magnesium: c_Mg + (solid terms only if stable mode) = W_Mg_total
        0 ~ (c_Mg + stable * (s_MgSO4 + s_MgNO32 + s_MgCl2)) -
            W_Mg_total,

        # === Charge balance (Eq. 9) ===
        0 ~ (c_H + c_Na + c_NH4 + 2 * c_Ca + c_K + 2 * c_Mg) -
            (c_Cl + 2 * c_SO4 + c_HSO4 + c_NO3 + c_OH),

        # === Aqueous equilibrium expressions (Eqs. 10-14) ===
        # 10. HSO4⁻ ↔ H⁺ + SO4²⁻ (K1)
        0 ~ m_H * m_SO4 * γ_H2SO4^3 - K1 * m_HSO4 * γ_HHSO4^2,

        # 11. NH₃(g) + H₂O ↔ NH₄⁺ + OH⁻ (K2 = K21×K22)
        # Note: Should use water activity (aw) instead of RH for concentrated solutions
        # Activity coefficient ratio approximates the mean ionic activity coefficient
        0 ~ m_NH4 * m_OH * (γ_NH4Cl / γ_HCl)^2 - K2 * RH * p_NH3,

        # 12. HNO₃(g) ↔ H⁺ + NO₃⁻ (K4)
        0 ~ m_H * m_NO3 * γ_HNO3^2 - K4 * p_HNO3,

        # 13. HCl(g) ↔ H⁺ + Cl⁻ (K3)
        0 ~ m_H * m_Cl * γ_HCl^2 - K3 * p_HCl,

        # 14. H₂O ↔ H⁺ + OH⁻ (Kw)
        # Note: Should use water activity (aw) instead of RH for concentrated solutions
        # For dilute solutions: aw ≈ RH, but for concentrated aerosols this approximation has limitations
        0 ~ m_H * m_OH - Kw * RH,

        # === Water content with MDRH cutoff (ZSR, Eqs. 15-18, 22) ===
        # Factor (1 - stable*(1 - mdrh_factor)):
        #   Metastable (stable=0): factor = 1.0 → ZSR water unchanged
        #   Stable, RH > MDRH: factor ≈ 1.0 → normal water content
        #   Stable, RH < MDRH: factor ≈ 0.0 → dry particle (W_w → tiny)
        0 ~ W_w -
            (1.0 - stable * (1.0 - mdrh_factor)) *
            _iso2_zsr_water(
            RH, c_Na, c_NH4, c_SO4,
            c_HSO4, c_NO3, c_Cl,
            c_Ca, c_K, c_Mg
        ),

        # === Ionic strength (Eq. 16) ===
        0 ~ I_s - 0.5 * (
            m_H + m_Na + m_NH4 + m_Cl + 4 * m_SO4 + m_HSO4 +
                m_NO3 + 4 * m_Ca + m_K + 4 * m_Mg + m_OH
        ),

        # === Activity coefficients — original 5 (Eqs. 17-21) ===
        γ_H2SO4 ~ _iso2_gamma_T(_iso2_km_gamma(-0.1, 2, I_s), I_s, T_env),
        γ_HHSO4 ~ _iso2_gamma_T(_iso2_km_gamma(8.0, 1, I_s), I_s, T_env),
        γ_HNO3 ~ _iso2_gamma_T(_iso2_km_gamma(2.6, 1, I_s), I_s, T_env),
        γ_HCl ~ _iso2_gamma_T(_iso2_km_gamma(6.0, 1, I_s), I_s, T_env),
        γ_NH4Cl ~ _iso2_gamma_T(_iso2_km_gamma(0.82, 1, I_s), I_s, T_env),

        # === Activity coefficients — additional 11 (Eqs. 22-32) ===
        γ_NaCl ~ _iso2_gamma_T(_iso2_km_gamma(2.23, 1, I_s), I_s, T_env),
        γ_NaNO3 ~ _iso2_gamma_T(_iso2_km_gamma(-0.39, 1, I_s), I_s, T_env),
        γ_Na2SO4 ~ _iso2_gamma_T(_iso2_km_gamma(-0.19, 2, I_s), I_s, T_env),
        γ_NH42SO4 ~ _iso2_gamma_T(_iso2_km_gamma(-0.25, 2, I_s), I_s, T_env),
        γ_CaNO32 ~ _iso2_gamma_T(_iso2_km_gamma(0.93, 2, I_s), I_s, T_env),
        γ_CaCl2 ~ _iso2_gamma_T(_iso2_km_gamma(2.4, 2, I_s), I_s, T_env),
        γ_K2SO4 ~ _iso2_gamma_T(_iso2_km_gamma(-0.25, 2, I_s), I_s, T_env),
        γ_KNO3 ~ _iso2_gamma_T(_iso2_km_gamma(-2.33, 1, I_s), I_s, T_env),
        γ_KCl ~ _iso2_gamma_T(_iso2_km_gamma(0.92, 1, I_s), I_s, T_env),
        γ_MgNO32 ~ _iso2_gamma_T(_iso2_km_gamma(2.32, 2, I_s), I_s, T_env),
        γ_MgCl2 ~ _iso2_gamma_T(_iso2_km_gamma(2.9, 2, I_s), I_s, T_env),

        # === Solid precipitation NCP equations (Eqs. 33-50) ===
        # Each: stable * FB(s, 1 - SR) + (1-stable) * s = 0
        # When stable=0 (metastable): s = 0 for all solids
        # When stable=1: either s=0 or SR=1 (at saturation)

        # 33. CaSO4 — always insoluble (Table 4 footnote a), handled as smooth min constraint
        0 ~ stable * _iso2_fb_ncp(s_CaSO4, 1 - SR_CaSO4) +
            (1 - stable) * s_CaSO4,

        # 34. Ca(NO3)2
        0 ~ stable * _iso2_fb_ncp(s_CaNO32, 1 - SR_CaNO32) +
            (1 - stable) * s_CaNO32,

        # 35. CaCl2
        0 ~ stable * _iso2_fb_ncp(s_CaCl2, 1 - SR_CaCl2) +
            (1 - stable) * s_CaCl2,

        # 36. K2SO4
        0 ~ stable * _iso2_fb_ncp(s_K2SO4, 1 - SR_K2SO4) +
            (1 - stable) * s_K2SO4,

        # 37. KHSO4
        0 ~ stable * _iso2_fb_ncp(s_KHSO4, 1 - SR_KHSO4) +
            (1 - stable) * s_KHSO4,

        # 38. KNO3
        0 ~ stable * _iso2_fb_ncp(s_KNO3, 1 - SR_KNO3) +
            (1 - stable) * s_KNO3,

        # 39. KCl
        0 ~ stable * _iso2_fb_ncp(s_KCl, 1 - SR_KCl) +
            (1 - stable) * s_KCl,

        # 40. MgSO4
        0 ~ stable * _iso2_fb_ncp(s_MgSO4, 1 - SR_MgSO4) +
            (1 - stable) * s_MgSO4,

        # 41. Mg(NO3)2
        0 ~ stable * _iso2_fb_ncp(s_MgNO32, 1 - SR_MgNO32) +
            (1 - stable) * s_MgNO32,

        # 42. MgCl2
        0 ~ stable * _iso2_fb_ncp(s_MgCl2, 1 - SR_MgCl2) +
            (1 - stable) * s_MgCl2,

        # 43. Na2SO4
        0 ~ stable * _iso2_fb_ncp(s_Na2SO4, 1 - SR_Na2SO4) +
            (1 - stable) * s_Na2SO4,

        # 44. NaHSO4
        0 ~ stable * _iso2_fb_ncp(s_NaHSO4, 1 - SR_NaHSO4) +
            (1 - stable) * s_NaHSO4,

        # 45. NaNO3
        0 ~ stable * _iso2_fb_ncp(s_NaNO3, 1 - SR_NaNO3) +
            (1 - stable) * s_NaNO3,

        # 46. NaCl
        0 ~ stable * _iso2_fb_ncp(s_NaCl, 1 - SR_NaCl) +
            (1 - stable) * s_NaCl,

        # 47. (NH4)2SO4
        0 ~ stable * _iso2_fb_ncp(s_NH42SO4, 1 - SR_NH42SO4) +
            (1 - stable) * s_NH42SO4,

        # 48. NH4HSO4
        0 ~ stable * _iso2_fb_ncp(s_NH4HSO4, 1 - SR_NH4HSO4) +
            (1 - stable) * s_NH4HSO4,

        # 49. Letovicite (NH4)3H(SO4)2
        0 ~ stable * _iso2_fb_ncp(s_LC, 1 - SR_LC) +
            (1 - stable) * s_LC,

        # 50. NH4NO3 (gas-solid equilibrium)
        0 ~ stable * _iso2_fb_ncp(s_NH4NO3, 1 - SR_NH4NO3) +
            (1 - stable) * s_NH4NO3,

        # 51 → renumbered as Eq. 50 in the 50-equation system
        # NH4Cl (gas-solid equilibrium)
        0 ~ stable * _iso2_fb_ncp(s_NH4Cl, 1 - SR_NH4Cl) +
            (1 - stable) * s_NH4Cl,
    ]

    return System(eqs, t; name)
end

export IsorropiaEquilibrium
