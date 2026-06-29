using CubeTools

H2 = open_cube("H2.cube")
println("Arquivo Aberto com sucesso!")

println("Calculando o gradiente do arquivo H2.cube...")
∇_H2 = Gradient(H2)
println("Gradiente calculado com sucesso!")