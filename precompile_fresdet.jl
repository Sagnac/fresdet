#=

using PackageCompiler

create_sysimage([:Images, :FFTW, :StatsBase, :GLMakie];
                sysimage_path = "fresdet.ext",
                precompile_execution_file = "precompile_fresdet.jl")

.ext = .so on Linux, .dll on Windows, .dylib on macOS

When the plots appear interact with the program so that all of the functions are
precompiled and then close the windows to continue the precompilation process.

The sysimage can be loaded at start-up using the -J or --sysimage switch.

=#

include("fresdet.jl")
HTG, FFT, M = fresdet("images/image.png");
save("images/Histogram.png", HTG)
save("images/FFT.png", FFT)
savefft("images/barefft.png", M)
wait(fftscreen)
wait(Hscreen)
