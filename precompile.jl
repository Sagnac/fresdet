#=

using PackageCompiler

create_sysimage([:Images, :FFTW, :StatsBase, :GLMakie];
                sysimage_path = "fresdet.ext",
                precompile_execution_file = "precompile.jl")

.ext = .so on Linux, .dll on Windows, .dylib on macOS

When the plots appear interact with the figure so that all of the functions are
precompiled and then close the window to continue the precompilation process.

The sysimage can be loaded at start-up using the -J or --sysimage switch.

=#

include("fresdet.jl")
fig, S = fresdet("images/image.png", true);
save("images/fresdet.png", fig)
savefft("images/barefft.png", S)
