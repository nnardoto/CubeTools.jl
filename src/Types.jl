const BohrToAngstrom = 0.529177210903
const AngstromToBohr = 1.0 / BohrToAngstrom
const HartreeToEV = 27.211386245988
const EVToHartree = 1.0 / HartreeToEV


struct Atom
    Z        :: Int
    charge   :: Float64
    position :: SVector{3, Float64}
end


mutable struct CubeFile
    title       :: Tuple{String, String}
    origin      :: SVector{3, Float64}
    dl          :: SMatrix{3, 3, Float64, 9}
    npoints     :: SVector{3, Int}
    periodicity :: SVector{3, Bool}
    atoms       :: Vector{Atom}
    data        :: Vector{Float64}
end
