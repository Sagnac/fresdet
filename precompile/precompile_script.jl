include("fresdet.jl")
fig, S = fresdet("images/image.png"; script = true);
save("images/fresdet.png", fig)
savefft("images/barefft.png", S)
