# Interpola o potencial eletrostático armazenado no CubeFile na posição de
# cada átomo, por síntese de Fourier (interpolação trigonométrica/band-limited),
# a partir da FFT dos valores amostrados na malha.
#
# (mesma justificativa de malha não ortogonal / coordenadas fracionárias u
#  do cabeçalho original — nada muda aqui, pois é só sobre AMOSTRAGEM)
#
# g(u) = (1/N) Σ_m V̂(m) · exp(i·q·u),   q_a = 2π·m_a
#
# Diferente da versão de campo elétrico, aqui não há gradiente nem regra da
# cadeia: T/T_inv só são usados para converter a posição cartesiana do átomo
# em coordenada fracionária u, não para transformar nenhuma derivada.
function ValueAtAtomicPositions(cube::CubeFile) :: Vector{Float64}
    V_r = cube.data
    nx = cube.npoints[1]; ny = cube.npoints[2]; nz = cube.npoints[3]
    N  = nx * ny * nz

    T     = SMatrix{3,3,Float64}(cube.dl .* cube.npoints)
    T_inv = inv(T)

    qx = fftfreq(nx, nx) * 2π
    qy = fftfreq(ny, ny) * 2π
    qz = fftfreq(nz, nz) * 2π

    V_k = fft(V_r)

    Value = Float64[]

    for atom in cube.atoms
        u = T_inv' * (atom.position - cube.origin)

        φx = reshape(exp.(im .* qx .* u[1]), nx,  1,  1)
        φy = reshape(exp.(im .* qy .* u[2]),  1, ny,  1)
        φz = reshape(exp.(im .* qz .* u[3]),  1,  1, nz)

        v = real(sum(V_k .* φx .* φy .* φz)) / N

        push!(Value, v)
    end

    return Value
end