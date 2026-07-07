using LinearAlgebra

# Calcula o campo elétrico na posição de cada átomo a partir do potencial
# eletrostático armazenado no CubeFile, por diferenciação espectral (FFT).
#
# Válida tanto para malhas ORTOGONAIS quanto NÃO ORTOGONAIS: a malha do .cube
# é definida por três vetores de passo `cube.dl[i, :]` (não necessariamente
# perpendiculares entre si), e o ponto de índice (ix,iy,iz) fica em
#   r(ix,iy,iz) = origin + ix·dl[1,:] + iy·dl[2,:] + iz·dl[3,:].
#
# --- Ideia central ---------------------------------------------------------
# A FFT assume uma função amostrada em pontos igualmente espaçados ao longo
# de eixos independentes. Isso NÃO é verdade em r (Cartesiano) quando a malha
# é não ortogonal, mas é sempre verdade nas coordenadas FRACIONÁRIAS u ∈ [0,1)³
# associadas à malha, já que por construção
#   r(ix,iy,iz) − origin = T' · (ix/nx, iy/ny, iz/nz) = T' · u,
# onde T é a matriz cujas linhas são os vetores de rede completos da célula:
#   T[i, :] = npoints[i] · dl[i, :]     (dl escalado do passo para a aresta inteira)
#
# Logo, tratamos g(u) := V(r(u)) como função periódica de período 1 em cada
# u_a e aplicamos a FFT padrão nela — o array cube.data já está amostrado
# exatamente nos pontos u = (ix/nx, iy/ny, iz/nz), então basta usar `fft`
# diretamente sobre cube.data, sem qualquer referência a distâncias físicas.
#
# A série de Fourier de g em u usa frequências adimensionais (modos inteiros)
#   q_a = 2π·m_a,  m_a ∈ {0, 1, ..., ⌈n_a/2⌉−1, −⌊n_a/2⌋, ..., −1}
# obtidas com `fftfreq(n, n)` (passar fs = n devolve exatamente os inteiros m).
# O gradiente espectral em u é então ∇_u g = i·q·ĝ(m), avaliado num ponto u
# arbitrário (não necessariamente um nó da malha) por síntese de Fourier direta:
#   ∂g/∂u_a (u) = (1/N) Σ_m  i·q_a·V̂(m) · exp(i·q·u)
#
# --- Voltando para coordenadas cartesianas ---------------------------------
# Pela regra da cadeia, com r − origin = T'·u:
#   ∂r_j/∂u_a = T[a,j]   ⇒   ∇_u g = T · ∇_r V   ⇒   ∇_r V = T⁻¹ · ∇_u g
# (produto matriz-vetor com a matriz inversa de T, e NÃO sua transposta).
# Essa relação se reduz exatamente ao caso ortogonal simples (k físico = q/L)
# quando dl é diagonal, então esta é a ÚNICA função de campo elétrico
# necessária — ela substitui as antigas `ElectricField`/`GeneralElectricField`.
#
# O campo elétrico é o gradiente negativo do potencial: E = -∇_r V.
function ElectricField(cube::CubeFile, units::Symbol = :bohr) :: Vector{SVector{3, Float64}}
    V_r = cube.data
    nx = cube.npoints[1]; ny = cube.npoints[2]; nz = cube.npoints[3]
    N  = nx * ny * nz

    # Matriz de rede da célula completa (linhas = vetores de rede a_i) e sua inversa.
    T     = SMatrix{3,3,Float64}(cube.dl .* cube.npoints)
    T_inv = inv(T)

    # Frequências de Fourier adimensionais (modos inteiros × 2π), conjugadas às
    # coordenadas fracionárias u — não confundir com números de onda físicos.
    qx = fftfreq(nx, nx) * 2π
    qy = fftfreq(ny, ny) * 2π
    qz = fftfreq(nz, nz) * 2π

    V_k = fft(V_r)

    # Gradiente espectral em relação a u: multiplicação por i·q
    Ûx = im .* reshape(qx, nx,  1,  1) .* V_k
    Ûy = im .* reshape(qy,  1, ny,  1) .* V_k
    Ûz = im .* reshape(qz,  1,  1, nz) .* V_k

    EField = SVector{3, Float64}[]

    for atom in cube.atoms
        # Coordenada fracionária do átomo: u = (T')⁻¹ · (r − origin) = T_inv' · (r − origin)
        u = T_inv' * (atom.position - cube.origin)

        # Fatores de fase exp(i·q·u) para síntese de Fourier na posição u
        φx = reshape(exp.(im .* qx .* u[1]), nx,  1,  1)
        φy = reshape(exp.(im .* qy .* u[2]),  1, ny,  1)
        φz = reshape(exp.(im .* qz .* u[3]),  1,  1, nz)

        # A FFT normaliza por N, então dividimos pelo número total de pontos.
        # real() descarta a parte imaginária residual (erro numérico de arredondamento).
        gu = SVector{3, Float64}(
            real(sum(Ûx .* φx .* φy .* φz)) / N,
            real(sum(Ûy .* φx .* φy .* φz)) / N,
            real(sum(Ûz .* φx .* φy .* φz)) / N,
        )

        # Regra da cadeia: ∇_r V = T⁻¹ · ∇_u g. Campo elétrico = -∇_r V.
        push!(EField, -(T_inv * gu))
    end

    if units == :angstrom
        # V está em Hartree e r em Bohr, logo E está em Hartree/(e·Bohr).
        # Conversão para eV/Å: Hartree→eV multiplica, Bohr→Å divide.
        EField .= EField .* (HartreeToEV / BohrToAngstrom)
    end

    return EField
end

