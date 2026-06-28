function open_cube(FileName::String, periodicity=(true, true, true))::CubeFile
    open(FileName) do file
        title1 = readline(file)
        title2 = readline(file)

        origin_line = split(readline(file))
        n_atoms = parse(Int, origin_line[1])
        origin  = SVector{3, Float64}(parse.(Float64, origin_line[2:4]))

        npts = zeros(Int, 3)
        axes_buf = zeros(Float64, 3, 3)
        for i in 1:3
            axis_line      = split(readline(file))
            npts[i]        = parse(Int, axis_line[1])
            axes_buf[i, :] = parse.(Float64, axis_line[2:4])
        end
        npoints = SVector{3, Int}(npts)
        axes = SMatrix{3, 3, Float64, 9}(axes_buf)

        atoms = Atom[]
        for _ in 1:n_atoms
            atom_line = split(readline(file))
            Z         = parse(Int,     atom_line[1])
            charge    = parse(Float64, atom_line[2])
            pos       = SVector{3, Float64}(parse.(Float64, atom_line[3:5]))
            push!(atoms, Atom(Z, charge, pos))
        end

        data = Float64[]
        for line in eachline(file)
            append!(data, parse.(Float64, split(line)))
        end

        # Verifica se o número de pontos lidos corresponde ao esperado
        expected_points = prod(npoints)
        if length(data) != expected_points
            error("O número de pontos lidos ($length(data)) não corresponde ao esperado ($expected_points).")
        end 

        # Cria o objeto CubeFile com os dados lidos
        CubeFile(
            (title1, title2),
            origin,
            axes,
            npoints,
            SVector{3, Bool}(periodicity),
            atoms,
            data,
        )
    end
end

function save_cube(cube::CubeFile, FileName::String)
    open(FileName, "w") do file
        println(file, cube.title[1])
        println(file, cube.title[2])

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

        for (i, value) in enumerate(cube.data)
            @printf(file, " %12.5E", value)
            if i % 6 == 0
                println(file)
            end
        end
        println(file)
    end
end
