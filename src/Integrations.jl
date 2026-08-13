# Integrations.jl
# Integração numérica de campos volumétricos (cube.data) sobre regiões
# localizadas em torno de átomos.
#
# --- Passo 1 (validação): integral sobre uma esfera simples ---------------
# Estratégia: contagem de voxels (soma de Riemann). Para cada átomo, soma-se
# cube.data em todo ponto de grade cuja posição cartesiana cai dentro da
# esfera de raio R, multiplicado pelo volume da célula de grade
#   dV = |det(cube.dl)|
# (o paralelepípedo formado pelos três vetores de passo, válido também para
# malhas não ortogonais).
#
# NOTA (limitação conhecida desta primeira versão): não há wrap de imagens
# periódicas. Átomos próximos da borda da célula, com R maior que a distância
# até a borda, terão a esfera truncada artificialmente. Fica para um próximo
# passo, junto com testes de convergência em função da resolução da malha.

"""
    SphereIntegral(cube::CubeFile, R::Real; units::Symbol = :bohr) -> Vector{Float64}

Integra `cube.data` numa esfera de raio `R` centrada em cada átomo de
`cube.atoms`, por contagem de voxels.

`R` é interpretado nas unidades indicadas por `units` (`:bohr` ou
`:angstrom`); internamente é convertido para Bohr, a unidade nativa do
arquivo .cube. Retorna um vetor com um valor por átomo, na mesma ordem de
`cube.atoms`.

# Exemplo
```julia
cube = open_cube("Telurio.vhart.cube")
Q = SphereIntegral(cube, 1.0; units = :angstrom)  # integral em torno de cada átomo, R = 1 Å
```
"""
function SphereIntegral(cube::CubeFile, R::Real; units::Symbol = :bohr)::Vector{Float64}
    units in (:bohr, :angstrom) || error("units deve ser :bohr ou :angstrom, recebido: $(units)")
    R_bohr = units == :angstrom ? R * AngstromToBohr : R
    R_bohr > 0 || error("O raio da esfera deve ser positivo.")

    nx, ny, nz = cube.npoints
    dV = abs(det(cube.dl))   # volume da célula de grade (paralelepípedo dl1 x dl2 x dl3)

    integrals = zeros(Float64, length(cube.atoms))

    for (a, atom) in enumerate(cube.atoms)
        acc = 0.0
        for ix in 0:nx-1, iy in 0:ny-1, iz in 0:nz-1
            r = cube.origin .+ ix .* cube.dl[1, :] .+ iy .* cube.dl[2, :] .+ iz .* cube.dl[3, :]
            if norm(r - atom.position) <= R_bohr
                acc += cube.data[ix+1, iy+1, iz+1]
            end
        end
        integrals[a] = acc * dV
    end

    return integrals
end
