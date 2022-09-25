# fresdet.jl
# Version 0.2.0
# 2022-9-25
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
    fftfig = Figure(resolution = res)
    Sv = vec(S)
    n, s6 = stats(Sv)
    local fftaxis, rect

    # compute the measures and generate the histogram
    Haxis = GLMakie.Axis(Hfig[1,1], title = "$w x $h pixels", titlesize = t,
                         xlabel = "pixel value [0-255]", ylabel = "pixel count")
    H = hist!(Haxis, Sv, bins = n, color = RGBf(0, 0.5, 0.8))
    # abuse the legend since it's easier than making a custom text! box
    legend = axislegend(Haxis, [MarkerElement(marker = "", color = :transparent)],
                        [s6], titlesize = t, labelsize = t, position = :rc)
    DataInspector(Hfig)

    # sets the fft image
    function I()
        fftaxis = GLMakie.Axis(fftfig[1,1], aspect = DataAspect(),
                               title = "Resolution: $w x $h",
                               titlesize = t, subtitlesize = t)
        deregister_interaction!(fftaxis, :rectanglezoom)
        image!(fftaxis, S, colormap = :tokyo, interpolate = false)
        rect = select_rectangle(fftaxis)
        DataInspector(fftfig)
        return
    end

    # render
    I()
    set_window_config!(title = "Histogram", focus_on_show = false)
    global Hscreen = display(Hfig)
    set_window_config!(title = "FFT", focus_on_show = true)
    global fftscreen = display(GLMakie.Screen(), fftfig)
    set_window_config!(title = "Histogram", focus_on_show = false)

    # rectangle selection callback
    R = nothing
    local sliderx, slidery, Lx, Ly, xs, ys
    listen() = on(rect) do v
        # bounds
        x01, y01 = round.(Int, v.origin)
        Lx0, Ly0 = round.(Int, v.widths)
        x02, y02 = Lx0 + x01, Ly0 + y01
        x01 = clamp(x01, 0, w - 1)
        y01 = clamp(y01, 0, h - 1)
        x02 = clamp(x02, x01 + 1, w)
        y02 = clamp(y02, y01 + 1, h)
        Lx0, Ly0 = x02 - x01, y02 - y01

        # sliders
        if R === nothing
            sliderx = IntervalSlider(fftfig[2,1], range = 0:w,
                                     startvalues = (x01, x02), snap = false)
            slidery = IntervalSlider(fftfig[1,2], range = 0:h,
                                     startvalues = (y01, y02), snap = false,
                                     horizontal = false)
            x1 = @lift(clamp($(sliderx.interval)[1], 0, w - 1))
            x2 = @lift(clamp($(sliderx.interval)[2], $x1 + 1, w))
            y1 = @lift(clamp($(slidery.interval)[1], 0, h - 1))
            y2 = @lift(clamp($(slidery.interval)[2], $y1 + 1, h))
            Lx, Ly = @lift($x2 - $x1), @lift($y2 - $y1)
            A = @lift($Lx * $Ly)

            # draw rectangle
            r4 = @lift(Rect($x1, $y1, $Lx, $Ly))
            R = lines!(fftaxis, r4, linestyle = :dot,
                       color = RGBf(0, 1, 0.38), linewidth = 7)

            # create new histogram on selection
            Sv = @lift(vec(S[$x1+1:$x2, $y1+1:$y2]))
            delete!(Haxis, H)
            H = hist!(Haxis, Sv, color = RGBf(0, 0.5, 0.8))

            # interval change callback
            on(A; update = true) do _
                xs, ys = round.((w / Lx[], h / Ly[]); digits = 2)
                fftaxis.title = "Effective Resolution: $(Lx[]) x $(Ly[])"
                fftaxis.subtitle = "x-scale: $xs, y-scale: $ys"
                H.bins, legend.entrygroups[][1][2][1].label = stats(Sv[])
                Haxis.title = "$(Lx[]) x $(Ly[]) pixels"
                reset_limits!(Haxis)
                nothing
            end

        else

            set_close_to!(sliderx, x01, x02)
            sliderx.startvalues = (x01, x02)
            set_close_to!(slidery, y01, y02)
            slidery.startvalues = (y01, y02)

        end

        println("\nEffective Resolution: $(Lx[]) x $(Ly[])")
        println("x-scale: $xs\ny-scale: $ys")

        # re-open the histogram/stats window on request
        !events(Hfig).window_open[] && (Hscreen = display(Hfig))

        nothing

    end

    listen()

    # erase the lines and remove the sliders if the histogram window closes
    on(events(Hfig).window_open) do Hwindow
        if !Hwindow && events(fftfig).window_open[] && R !== nothing
            empty!(fftaxis)
            empty!(fftfig)
            I()
            listen()
            R = nothing
        end
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
