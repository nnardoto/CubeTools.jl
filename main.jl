push!(LOAD_PATH, joinpath(@__DIR__, "src"))

using CubeTools
using StaticArrays

H2 = open_cube("H2.cube")

copyH2 = -3*H2 + 3*H2

save_cube(copyH2, "copyH2_out.cube")


#write_cube(H2, "H2.cube")
