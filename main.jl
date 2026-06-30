using CubeTools
using LinearAlgebra
using StaticArrays
using Statistics
using Printf

TelurioData = open_cube("Telurio_bulk_model.vhart.cube")
Display(TelurioData)

N = length(TelurioData.atoms)


EField = ElectricField(TelurioData)
E = []

# Obtém o módulo do campo elétrico para cada átomo
for i in 1:N 
    Ei = norm(EField[i])
    push!(E, Ei)
end

# Alguma Estatística
E_mean = mean(E)
E_std = std(E)
E_uncertainty = E_std/√N

ExpBase = floor(log10(abs(E_mean)))
E_mean_scaled = E_mean / 10^ExpBase
E_uncertainty_scaled = E_uncertainty / 10^ExpBase
println("[ElectricField] Electric Field Statistics (Bohr units)    ")
@printf("E    = %6.4E Ha/Bohr\n", E_mean)
@printf("σ    = %6.4E Ha/Bohr\n", E_std)
@printf("σ/√N = %6.4E Ha/Bohr\n", E_uncertainty)
println()
println("[ElectricField] Electric Field with Uncertainty    ")

@printf("E = (%6.4f ± %1.4f) × 10^(%+3d) Ha/Bohr\n", E_mean_scaled, E_uncertainty_scaled, ExpBase)
