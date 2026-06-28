module CubeTools
    using StaticArrays
    using Printf

    include("Types.jl")
    include("IO.jl")

    export CubeFile, Atom, read_cube, write_cube
end
