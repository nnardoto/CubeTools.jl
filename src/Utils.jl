function Display(cube :: CubeFile)
    # Calcula os vetores de dimensão do cubo
    a = cube.npoints[1] * cube.dl[1,:]
    b = cube.npoints[2] * cube.dl[2,:]
    c = cube.npoints[3] * cube.dl[3,:]


    # Desenha uma Janela usando Printf para exibir informações do arquivo .cube
    if cube.units == :bohr
        println("[Summary] Units: Bohr (length) and Hartree (energy)")
    elseif cube.units == :angstrom
        println("[Summary] Units: Angstrom (length) and eV (energy)")
    else
        println("[Summary] Units: Unknown")
    end
    println()  
    println("+" * "-"^50 * "+")
    println(cube.LineOne)
    println(cube.LineTwo)
    println("+" * "-"^50 * "+")
    println("[Summary] Number of atoms: ", length(cube.atoms))
    @printf("[Summary] Data size:       %d x %d x %d\n", size(cube.data)...)
    println()
    if cube.units == :bohr
        println("[Summary] Cell Vectors (Bohr):")
        @printf("  a: % 9.4f % 9.4f % 9.4f\n", a[1], a[2], a[3])
        @printf("  b: % 9.4f % 9.4f % 9.4f\n", b[1], b[2], b[3])
        @printf("  c: % 9.4f % 9.4f % 9.4f\n", c[1], c[2], c[3])
    elseif cube.units == :angstrom
        println("[Summary] Cell Vectors (Angstrom):")
        @printf("  a: % 9.4f % 9.4f % 9.4f\n", a[1], a[2], a[3])
        @printf("  b: % 9.4f % 9.4f % 9.4f\n", b[1], b[2], b[3])
        @printf("  c: % 9.4f % 9.4f % 9.4f\n", c[1], c[2], c[3])
    end
    println()
    if cube.units == :bohr
        println("[Summary] Atom Positions: (Bohr):")
    elseif cube.units == :angstrom
        println("[Summary] Atom Positions: (Angstrom):")
    end
    for (i, atom) in enumerate(cube.atoms)
        @printf("  Atom %d: % 9.4f % 9.4f % 9.4f\n", i, atom.position[1], atom.position[2], atom.position[3])
    end
    println("+" * "-"^50 * "+")

    println()
end


function Bohr2Angstrom(cube::CubeFile) :: CubeFile
    buffer = copy(cube)
    buffer.origin .= buffer.origin .* BohrToAngstrom
    buffer.dl .= buffer.dl .* BohrToAngstrom
    buffer.atoms = [Atom(atom.Z, atom.charge, atom.position * BohrToAngstrom) for atom in buffer.atoms]
    buffer.data .= buffer.data .* (HartreeToEV / BohrToAngstrom)
    buffer.units = :angstrom
    return buffer
end

function Angstrom2Bohr(cube::CubeFile) :: CubeFile
    buffer = copy(cube)
    buffer.origin .= buffer.origin .* AngstromToBohr
    buffer.dl .= buffer.dl .* AngstromToBohr
    buffer.atoms = [Atom(atom.Z, atom.charge, atom.position * AngstromToBohr) for atom in buffer.atoms]
    buffer.data .= buffer.data .* (EVToHartree * BohrToAngstrom)
    buffer.units = :bohr
    return buffer
end