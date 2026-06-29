function Gradient(cube::CubeFile) :: Vector{CubeData}
    nx = cube.npoints[1], ny = cube.npoints[2], nz = cube.npoints[3]

    data_fft = fft(cube.data)

    d = Float64[norm(col) for col in eachcol(cube.dl)]

    kx = 2π .* fftfreq(nx, d[1])
    ky = 2π .* fftfreq(ny, d[2])
    kz = 2π .* fftfreq(nz, d[3])

    # Meshgrid de frequências
    KX = reshape(kx, nx, 1, 1)
    KY = reshape(ky, 1, ny, 1)
    KZ = reshape(kz, 1, 1, nz)

    ∇x = real(ifft(im .* KX .* data_fft))
    ∇y = real(ifft(im .* KY .* data_fft))
    ∇z = real(ifft(im .* KZ .* data_fft))

    return [∇x, ∇y, ∇z]
end

