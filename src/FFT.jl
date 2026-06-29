# Calcula o gradiente de um campo escalar armazenado em um CubeFile
# usando diferenciação espectral (via FFT).
#
# O gradiente é calculado no espaço de Fourier como:
#   ∇f(r) = FFT⁻¹( i·k · F̂(k) )
# onde F̂(k) é a transformada de Fourier de f(r) e k é o vetor de onda.
#
# Retorna um vetor com três CubeData: [∇x, ∇y, ∇z]
function Gradient(cube::CubeFile) :: Vector{CubeData}
    nx = cube.npoints[1]; ny = cube.npoints[2]; nz = cube.npoints[3]

    data_fft = fft(cube.data)

    # Comprimento físico de cada passo de malha em Bohr (norma do vetor dl[i,:])
    d = Float64[norm(row) for row in eachrow(cube.dl)]

    # Frequências de Fourier em rad/Bohr para cada eixo.
    # Julia: fftfreq(n, fs) espera fs = taxa de amostragem = 1/espaçamento (Bohr⁻¹).
    kx = 2π .* fftfreq(nx, 1/d[1])
    ky = 2π .* fftfreq(ny, 1/d[2])
    kz = 2π .* fftfreq(nz, 1/d[3])

    # Reshape para broadcast 3D: kx ao longo do eixo 1, ky do 2, kz do 3
    KX = reshape(kx, nx, 1, 1)
    KY = reshape(ky, 1, ny, 1)
    KZ = reshape(kz, 1, 1, nz)

    # Diferenciação espectral: multiplicação por i·k no espaço de Fourier
    ∇x = real(ifft(im .* KX .* data_fft))
    ∇y = real(ifft(im .* KY .* data_fft))
    ∇z = real(ifft(im .* KZ .* data_fft))

    return [∇x, ∇y, ∇z]
end
