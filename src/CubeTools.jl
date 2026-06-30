# Módulo principal do CubeTools.jl
# Ferramentas para leitura, escrita e análise de arquivos .cube (formato Gaussian/VASP)
module CubeTools
    using StaticArrays   # Vetores e matrizes de tamanho fixo em tempo de compilação
    using Printf         # Formatação de saída numérica
    using FFTW           # Transformada de Fourier rápida
    using LinearAlgebra  # Norma de vetores e álgebra linear
    using Interpolations # Interpolação trilinear para avaliação em pontos arbitrários

    include("Types.jl")        # Tipos de dados: Atom, CubeFile e constantes físicas
    include("IO.jl")           # Leitura e escrita de arquivos .cube
    include("Arithmetics.jl")  # Operações aritméticas entre CubeFiles
    include("FFT.jl")          # Gradiente via transformada de Fourier
    include("Interpolations.jl") # Campo elétrico por interpolação espectral
    include("Utils.jl")        # Funções utilitárias, como Display

    export CubeFile, Atom, open_cube, save_cube, Gradient, ElectricField, Display
end
