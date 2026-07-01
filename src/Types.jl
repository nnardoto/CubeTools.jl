# Constantes de conversão de unidades atômicas
const BohrToAngstrom = 0.529177210903   # 1 Bohr em Ångström (CODATA 2018)
const AngstromToBohr = 1.0 / BohrToAngstrom
const HartreeToEV    = 27.211386245988  # 1 Hartree em elétron-volt (CODATA 2018)
const EVToHartree    = 1.0 / HartreeToEV


# Alias para o array 3D de dados volumétricos (em unidades atômicas)
const CubeData = Array{Float64, 3}


# Representa um átomo no arquivo .cube
# Z        : número atômico
# charge   : carga nuclear efetiva (campo auxiliar do formato .cube)
# position : posição em unidades atômicas (Bohr)
struct Atom
    Z        :: Int
    charge   :: Float64
    position :: SVector{3, Float64}
end


# Representa um arquivo .cube completo
# LineOne/LineTwo : linhas de título/comentário do cabeçalho
# origin          : origem da malha em Bohr
# dl              : vetores de passo da malha (linhas = eixos x, y, z) em Bohr
# npoints         : número de pontos em cada direção
# periodicity     : flags de periodicidade por eixo
# atoms           : lista de átomos
# data            : densidade eletrônica em cada ponto da malha (nx × ny × nz)
mutable struct CubeFile
    LineOne     :: String
    LineTwo     :: String
    origin      :: SVector{3, Float64}
    dl          :: SMatrix{3, 3, Float64, 9}
    npoints     :: SVector{3, Int64}
    periodicity :: NTuple{3, Bool}
    atoms       :: Vector{Atom}
    data        :: CubeData
    units       :: Symbol  # :bohr ou :angstrom
end
