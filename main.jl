push!(LOAD_PATH, joinpath(@__DIR__, "src"))

using CubeTools
using StaticArrays

H2 = read_cube("H2.cube")


println("Title: ", H2.title[1])
println("Title: ", H2.title[2])
println("Origin: ", H2.origin)
println("Number of points: ", H2.npoints)

write_cube(H2, "H2_out.cube")


#write_cube(H2, "H2.cube")
