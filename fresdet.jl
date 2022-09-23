# fresdet.jl
# Version 0.1.0
# 2022-9-23
# https://github.com/Sagnac/fresdet

# Simple analysis tool for estimating the original resolution
# of standard upscaled images using ffts and basic statistical metrics.

# Julia v1.8.0

using Images # v0.25.2
using FFTW # v1.5.0
using StatsBase # v0.33.21
using GLMakie # v0.6.13

function fftimage(file)
    image1 = file |> load |> rotr90
    grayfloat = image1 .|> Gray .|> Float64
    F = grayfloat |> fft |> fftshift
    F2 = F .|> abs2
    S = F2 .+ 1 .|> log
    S * 255 / maximum(S) .|> round
end

function stats(v)
    n::Int = maximum(v) - minimum(v) + 1
    c = (mean(v), median(v), mode(v), std(v), skewness(v), kurtosis(v))
    m, d, o, s, g, k = round.(c; digits = 3)
    s6 = "mean: $m\nmedian: $d\nmode: $o\nst.dev.: $s\nskew: $g\nkurtosis: $k"
    n, replace(s6, "-" => "\u2212")
end

Hscreen = nothing
fftscreen = nothing

function fresdet(file)

    S = fftimage(file)
    w, h = size(S)
    Hscreen === nothing || GLMakie.destroy!(Hscreen)
    fftscreen === nothing || GLMakie.destroy!(fftscreen)
    t = 25 # text size
    res = (1000, 1000)
    Hfig = Figure(resolution = res)
    local Haxis

    # computes the measures and generates the histogram
    function H(Sv, x, y)
        n, s6 = stats(Sv)
        Haxis = GLMakie.Axis(Hfig[1,1], xlabel = "pixel value [0-255]",
                             ylabel = "pixel count")
        hist!(Haxis, Sv, bins = n)
        # abuse the legend since it's easier than making a custom text! box
        Legend(Hfig[1,1], [MarkerElement(marker = "", color = :transparent)], [s6],
               "$x x $y", tellwidth = false, tellheight = false, titlesize = t,
               labelsize = t, halign = :right, valign = :center)
        DataInspector(Hfig)
        return
    end

    # set the fft image
    fftfig = Figure(resolution = res)
    fftaxis = GLMakie.Axis(fftfig[1,1], aspect = DataAspect(),
                           title = "Resolution: $w x $h",
                           subtitle = "x-scale: 1, y-scale: 1",
                           titlesize = t, subtitlesize = t)
    deregister_interaction!(fftaxis, :rectanglezoom)
    image!(fftaxis, S, colormap = :tokyo, interpolate = false)
    rect = select_rectangle(fftaxis)
    DataInspector(fftfig)

    # render
    H(vec(S), w, h)
    set_window_config!(title = "Histogram", focus_on_show = false)
    global Hscreen = display(Hfig)
    set_window_config!(title = "FFT", focus_on_show = true)
    global fftscreen = display(GLMakie.Screen(), fftfig)
    set_window_config!(title = "Histogram", focus_on_show = false)

    # rectangle selection callback
    R = nothing
    on(rect) do v
        # bounds
        x01, y01 = round.(Int, v.origin)
        Lx0, Ly0 = round.(Int, v.widths)
        x02, y02 = Lx0 + x01, Ly0 + y01
        x01 = clamp(x01, 0, w - 1)
        y01 = clamp(y01, 0, h - 1)
        x02 = clamp(x02, x01 + 1, w)
        y02 = clamp(y02, y01 + 1, h)
        Lx0, Ly0 = x02 - x01, y02 - y01
        xs, ys = round.((w / Lx0, h / Ly0); digits = 2)

        R === nothing || delete!(fftaxis, R)
        # draw rectangle
        R = lines!(fftaxis, Rect(x01, y01, Lx0, Ly0), linestyle = :dot,
                   color = RGBf(0, 1, 0.38), linewidth = 7)
        fftaxis.title = "Effective Resolution: $Lx0 x $Ly0"
        fftaxis.subtitle = "x-scale: $xs, y-scale: $ys"
        println("\nEffective Resolution: $Lx0 x $Ly0")
        println("x-scale: $xs\ny-scale: $ys")

        # wipe the figure to clear the legend and recreate the axis
        empty!(Haxis)
        empty!(Hfig)
        ss = S[x01+1:x02, y01+1:y02]
        H(vec(ss), Lx0, Ly0)

        # erase the lines if the histogram window closes
        on(events(Hfig).window_open) do Hwindow
            if !Hwindow && events(fftfig).window_open[]
                delete!(fftaxis, R)
                fftaxis.title = "Resolution: $w x $h"
                fftaxis.subtitle = "x-scale: 1, y-scale: 1"
                R = nothing
            end
            nothing
        end

        # re-open the histogram/stats window on request
        !events(Hfig).window_open[] && (Hscreen = display(Hfig))

        nothing

    end

    # script support
    if !isinteractive()
        wait(fftscreen)
        wait(Hscreen)
    end

    return Hfig, fftfig, S

end

# saves the fft 1:1 without the axis
function savefft(filename::String, S::Matrix)
    r = size(S)
    Z = image(S, colormap = :tokyo, interpolate = false,
              figure = (resolution = r,),
              axis = (width = r[1], height = r[2]),)
    hidespines!(Z.axis)
    hidedecorations!(Z.axis)
    save(filename, Z)
    return nothing
end

# script support
if !isempty(ARGS)
    HTG, FFT, M = fresdet(ARGS[1]);
    empty!(ARGS)
end

print()
