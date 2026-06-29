# Calcula o campo elétrico na posição de cada átomo a partir da densidade eletrônica.
#
# O campo elétrico é aproximado pelo gradiente negativo do potencial, que por sua vez
# é proporcional ao gradiente da densidade: E ≈ -∇ρ.
#
# O gradiente espectral é calculado no espaço de Fourier como  i·k·ρ̂(k)
# e avaliado nas posições atômicas via síntese de Fourier direta:
#   ∇ρ(r) = (1/N) Σ_k  i·k·ρ̂(k) · exp(i·k·(r − origin))
#
# Usar (r − origin) em vez de r é essencial porque a FFT assume que os dados
# começam na posição zero; a origem da célula é o deslocamento físico real.
function ElectricField(cube::CubeFile)
    ρ = cube.data
    nx = cube.npoints[1]; ny = cube.npoints[2]; nz = cube.npoints[3]
    N  = nx * ny * nz

    # Frequências de Fourier em rad/Bohr para cada eixo
    # fftfreq(n, 1/d) retorna frequências em 1/Bohr; × 2π converte para rad/Bohr
    kx = fftfreq(nx, 1 / cube.dl[1, 1]) * 2π
    ky = fftfreq(ny, 1 / cube.dl[2, 2]) * 2π
    kz = fftfreq(nz, 1 / cube.dl[3, 3]) * 2π

    ρ̂ = fft(ρ)

    # Gradiente no espaço de Fourier: multiplicação por i·k
    Êx = im .* reshape(kx, nx,  1,  1) .* ρ̂
    Êy = im .* reshape(ky,  1, ny,  1) .* ρ̂
    Êz = im .* reshape(kz,  1,  1, nz) .* ρ̂

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

        Ex = real(sum(Êx .* φx .* φy .* φz)) / N
        Ey = real(sum(Êy .* φx .* φy .* φz)) / N
        Ez = real(sum(Êz .* φx .* φy .* φz)) / N

        @printf("Campo elétrico no átomo %0d (Z=%d): (%14.6f, %14.6f, %14.6f)\n",
                i, atom.Z, Ex, Ey, Ez)
    end
end
