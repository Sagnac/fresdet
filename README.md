# [fresdet.jl](fresdet.jl?raw=true)

Simple analysis tool for estimating the original resolution of standard upscaled images using ffts and basic statistical metrics.

![FFT.png](FFT.png) ![Histogram.png](Histogram.png)

----

## Caveats

Note that the fft will not definitively tell you whether an image has been upscaled; only under the assumption that it has been can you use it to estimate by how much. In general spectral evaluation can be used to assess the relative amount of detail present in an image and there are various reasons why a spectrum can be limited. Namely, the image could simply be inherently low resolution lacking in fine detail, or it could have been manipulated and processed in various ways. For instance, the application of a low pass filter like a Gaussian blur will attenuate the high frequency amplitudes without changing the image size; in that case, if you know that only that filter has been applied you could possibly estimate the radius or variance parameters.

Upsampling with standard interpolation algorithms which convolve the image with a specific kernel can sometimes trace distinctive features in the Fourier transform, but concluding that this operation took place usually requires closer inspection and knowledge of how the various kernels behave. Since upscaling in this manner is essentially low pass filtering after zero insertion it can be tricky discerning the difference between this and other low pass operations; if a variety of operations have been performed on the image then it becomes increasingly difficult to know which ones as well as estimate the magnitude of each.

The specific application of this tool involves situations where you suspect an image has been upscaled, usually by a standard degree, and you want to see if there's evidence supporting this. For example, if you want to quickly check if a 1080p video has been upscaled from a 720p source you could take a screenshot, run it through the program, and if the measurement yields scaling values close to 1.5 there's a good chance it has been. Generally, the program allows you to quickly approximate the equivalent underlying resolution of an image relative to its dimensions (concretely in relation to the Nyquist limit).

----

## Loading an image

Installation of Julia is required as well as the following packages:
```
Images
FFTW
StatsBase
GLMakie
```

You can run the program from the Julia REPL:
```
include("fresdet.jl")
fresdet("image.png");
```

Calling fresdet returns three values: the histogram figure, the fft figure, and the fft matrix allowing you to save three figures:
```
HTG, FFT, M = fresdet("image.png");
save("Histogram.png", HTG)
save("FFT.png", FFT)
savefft("barefft.png", M) # saves the fft without the axis in its original size
```

Alternatively, you can run the program as a script, passing the image as an argument:
```
julia fresdet.jl image.png
```

If you don't want to wait for the packages to load (Julia is a JIT compiled language so loading some packages like GLMakie can be slow) or the first plot latency you can create a custom sysimage with the PackageCompiler.jl package which will precompile all of the necessary functions and speed things up. A precompilation script is provided in [precompile_fresdet.jl](precompile_fresdet.jl). Note that this sysimage can be quite large in size, but at the moment there isn't a way to provide a small relocatable binary of the program which works across all operating systems.

----

## Using the program

If the image has been significantly upsized it will be immediately apparent on the fft as many of the high frequency pixels away from the center will have values close to zero. Clicking and dragging on the fft figure will allow you to select an area which should correspond to the bandwidth of the original image in each direction. The dimensions of this rectangle as well as the scaling factors in both directions will be printed. If well selected, the dimensions of the rectangle should correspond to the original resolution of the resized image and the scaling factors will tell you by how much the image has been resized.

The histogram along with the statistics are provided in order to help fine tune the selection, allowing for a more accurate estimate of the resolution; the figure will update after a new area is selected on the fft. The general idea here is that sampling from any spurious frequency components near the edges of the main part of the spectrum should induce noticeable changes in the distribution. Ideally the values should be approximately normally distributed; the distribution should look unimodal with small differences between the mean, median, and mode as well as exhibit a low skewness and kurtosis. Any deviation from this should be apparent, giving you an approximate bound on the dimensions of the original image.
