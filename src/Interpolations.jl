# Calcula o campo elétrico na posição de cada átomo a partir da densidade eletrônica.
#
# O campo elétrico é o gradiente negativo do potencial eletrostático: E = -∇V.
#
# O gradiente espectral é calculado no espaço de Fourier como  i·k·ρ̂(k)
# e avaliado nas posições atômicas via síntese de Fourier direta:
#   ∇ρ(r) = (1/N) Σ_k  i·k·ρ̂(k) · exp(i·k·(r − origin))
#
# Usar (r − origin) em vez de r é essencial porque a FFT assume que os dados
# começam na posição zero; a origem da célula é o deslocamento físico real.
function ElectricField(cube::CubeFile) :: Vector{SVector{3, Float64}}
    V_r = cube.data
    nx = cube.npoints[1]; ny = cube.npoints[2]; nz = cube.npoints[3]
    N  = nx * ny * nz

    fsx = 1/cube.dl[1, 1]; fsy = 1/cube.dl[2, 2]; fsk = 1 / cube.dl[3, 3] 
    EField = SVector{3, Float64}[]

    # Frequências de Fourier em rad/Bohr para cada eixo.
    # Julia: fftfreq(n, fs) espera fs = taxa de amostragem = 1/espaçamento (Bohr⁻¹).
    # Passar o espaçamento d diretamente daria k d² vezes menor → fases erradas.
    kx = fftfreq(nx, fsx) * 2π
    ky = fftfreq(ny, fsy) * 2π
    kz = fftfreq(nz, fsk) * 2π

    V_k = fft(V_r)

    # Gradiente no espaço de Fourier: multiplicação por i·k
    Êx = im .* reshape(kx, nx,  1,  1) .* V_k
    Êy = im .* reshape(ky,  1, ny,  1) .* V_k
    Êz = im .* reshape(kz,  1,  1, nz) .* V_k

    ox, oy, oz = cube.origin

    for (i, atom) in enumerate(cube.atoms)
        # Posição do átomo relativa à origem da célula
        rx = atom.position[1] - ox
        ry = atom.position[2] - oy
        rz = atom.position[3] - oz

        # Fatores de fase exp(i·k·r) para síntese de Fourier em r
        φx = reshape(exp.(im .* kx .* rx), nx,  1,  1)
        φy = reshape(exp.(im .* ky .* ry),  1, ny,  1)
        φz = reshape(exp.(im .* kz .* rz),  1,  1, nz)

        # O campo elétrico é o gradiente negativo do potencial eletrostático.
        # A FFT normaliza por N, então precisamos dividir pelo número total de pontos.
        # A função real() é usada para descartar a parte imaginária residual da FFT.
        Ex = -real(sum(Êx .* φx .* φy .* φz)) / N
        Ey = -real(sum(Êy .* φx .* φy .* φz)) / N
        Ez = -real(sum(Êz .* φx .* φy .* φz)) / N

        push!(EField, SVector{3, Float64}(Ex, Ey, Ez))
    end

    return EField
end
