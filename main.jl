using CubeTools
using Printf
using LinearAlgebra

TelurioData = open_cube("Telurio.vhart.cube")

T = ElectricField(TelurioData, :angstrom)

# Imprime o resumo do arquivo .cube
Display(TelurioData)

# Imprime o campo elétrico em cada átomo em V/A
for (i, atom) in enumerate(TelurioData.atoms)
    @printf("Atom %d: E = (% 9.4f, % 9.4f, % 9.4f, % 9.4f) V/Å\n", i, T[i][1], T[i][2], T[i][3], norm(T[i]))
end 
