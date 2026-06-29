# Lê um arquivo .cube e retorna um objeto CubeFile.
# O formato .cube armazena os dados com x variando mais lentamente e z mais rapidamente.
# O array resultante tem shape (nx, ny, nz) com indexação [ix, iy, iz].
function open_cube(FileName::String, periodicity=(true, true, true))::CubeFile
    open(FileName) do file
        title1 = readline(file)
        title2 = readline(file)

        # Linha 3: número de átomos e origem da malha
        origin_line = split(readline(file))
        n_atoms = parse(Int, origin_line[1])
        origin  = SVector{3, Float64}(parse.(Float64, origin_line[2:4]))

        # Linhas 4–6: número de pontos e vetor de passo para cada eixo
        npts   = zeros(Int, 3)
        dl_buf = zeros(Float64, 3, 3)
        for i in 1:3
            axis_line    = split(readline(file))
            npts[i]      = parse(Int, axis_line[1])
            dl_buf[i, :] = parse.(Float64, axis_line[2:4])
        end
        npoints = SVector{3, Int}(npts)
        dl      = SMatrix{3, 3, Float64, 9}(dl_buf)

        # Linhas seguintes: átomos (Z, carga, x, y, z)
        atoms = Atom[]
        for _ in 1:n_atoms
            atom_line = split(readline(file))
            Z         = parse(Int,     atom_line[1])
            charge    = parse(Float64, atom_line[2])
            pos       = SVector{3, Float64}(parse.(Float64, atom_line[3:5]))
            push!(atoms, Atom(Z, charge, pos))
        end

        # Resto do arquivo: dados volumétricos em ordem x-outer, z-inner
        data = Float64[]
        for line in eachline(file)
            append!(data, parse.(Float64, split(line)))
        end

        if length(data) != prod(npoints)
            error("O número de pontos lidos ($(length(data))) não corresponde ao esperado ($(prod(npoints))).")
        end

        # O arquivo .cube armazena com x variando mais devagar e z mais rápido.
        # Julia é column-major: reshape com nz primeiro captura z como índice externo.
        # permutedims(_, (3,2,1)) reordena para a semântica [ix, iy, iz].
        data = permutedims(reshape(data, npoints[3], npoints[2], npoints[1]), (3, 2, 1))

        CubeFile(title1, title2, origin, dl, npoints, periodicity, atoms, data)
    end
end


# Escreve um objeto CubeFile em disco no formato .cube padrão.
# Os dados são escritos com x variando mais lentamente e z mais rapidamente,
# conforme a especificação do formato.
function save_cube(cube::CubeFile, FileName::String)
    open(FileName, "w") do file
        println(file, cube.LineOne)
        println(file, cube.LineTwo)

        @printf(file, "%5d %12.6f %12.6f %12.6f\n",
            length(cube.atoms), cube.origin[1], cube.origin[2], cube.origin[3])

        for i in 1:3
            @printf(file, "%5d %12.6f %12.6f %12.6f\n",
                cube.npoints[i], cube.dl[i, 1], cube.dl[i, 2], cube.dl[i, 3])
        end

        for atom in cube.atoms
            @printf(file, "%5d %12.6f %12.6f %12.6f %12.6f\n",
                atom.Z, atom.charge, atom.position[1], atom.position[2], atom.position[3])
        end

        # permutedims (3,2,1) transforma (nx,ny,nz) → (nz,ny,nx):
        # iteração column-major resulta em z variando mais rápido, conforme o formato
        data = permutedims(cube.data, (3, 2, 1))
        for (i, value) in enumerate(data)
            @printf(file, " %12.5E", value)
            if i % 6 == 0
                println(file)
            end
        end
        println(file)
    end
end
