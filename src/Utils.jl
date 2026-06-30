function Display(cube :: CubeFile)
    # Calcula os vetores de dimensão do cubo
    a = cube.npoints[1] * cube.dl[1,:]
    b = cube.npoints[2] * cube.dl[2,:]
    c = cube.npoints[3] * cube.dl[3,:]
   
    # Desenha uma Janela usando Printf para exibir informações do arquivo .cube
    println("+" * "-"^50 * "+")
    println(cube.LineOne)
    println(cube.LineTwo)
    println("+" * "-"^50 * "+")
    println("[Summary] Number of atoms: ", length(cube.atoms))
    @printf("[Summary] Data size:       %d x %d x %d\n", size(cube.data)...)
    println()
    println("[Summary] Cell Vectors (Bohr):")
    @printf("  a: % 9.4f % 9.4f % 9.4f\n", a[1], a[2], a[3])
    @printf("  b: % 9.4f % 9.4f % 9.4f\n", b[1], b[2], b[3])
    @printf("  c: % 9.4f % 9.4f % 9.4f\n", c[1], c[2], c[3])
    println()
    println("[Summary] Atom Positions: (Bohr)")
    for (i, atom) in enumerate(cube.atoms)
        @printf("  Atom %d: % 9.4f % 9.4f % 9.4f\n", i, atom.position[1], atom.position[2], atom.position[3])
    end
    println("+" * "-"^50 * "+")

    println()
end