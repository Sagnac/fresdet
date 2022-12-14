# fresdet.jl
# Version 0.7.1
# 2022-10-25
# https://github.com/Sagnac/fresdet

# Simple analysis tool for estimating the original resolution
# of standard upscaled images using ffts and basic statistical metrics.

# Julia v1.8.0

using Images # v0.25.2
using FFTW # v1.5.0
using StatsBase # v0.33.21
using GLMakie # v0.7.0
import GLMakie: Axis, Makie.async_latest as async_latest
import GLMakie: Makie.position2string as p2s
import GLMakie: Makie.color2text as c2t
import GLMakie: Makie.StaticVector as SV

p2s(p::SV{2}) = "x: $(round(Int, p[1]))\ny: $(round(Int, p[2]))"
c2t(name, x::Integer, y::Integer, z) = "$x, $y = $(Int(z))"

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

function fresdet(file, script = !isinteractive())

    S = fftimage(file)
    w, h = size(S)
    t = 25 # text size
    res = (2000, 1000)
    fig = Figure(resolution = res)

    local of1, of2, ov, AR

    # sliders
    xc = Observable(div(w, 2))
    yc = Observable(div(h, 2))
    sliderx = Slider(fig[2,2], range = 1:w,
                     startvalue = w, snap = false)
    slidery = Slider(fig[1,3], range = 1:h,
                     startvalue = h, snap = false, horizontal = false)
    Lx = sliderx.value
    Ly = slidery.value

    # point setting
    function xpoints(xc, Lx)
        if xc > div(w, 2)
            x2 = clamp(xc + cld(Lx, 2), xc + 1, w)
            x1 = x2 - Lx
        else
            x1 = clamp(xc - div(Lx, 2), 0, xc)
            x2 = Lx + x1
        end
        return x1, x2
    end

    function ypoints(yc, Ly)
        if yc > div(h, 2)
            y2 = clamp(yc + cld(Ly, 2), yc + 1, h)
            y1 = y2 - Ly
        else
            y1 = clamp(yc - div(Ly, 2), 0, yc)
            y2 = Ly + y1
        end
        return y1, y2
    end

    x = lift(xpoints, xc, Lx; ignore_equal_values = true)
    x1, x2 = [lift(x -> x[i], x; ignore_equal_values = true) for i in 1:2]

    y = lift(ypoints, yc, Ly; ignore_equal_values = true)
    y1, y2 = [lift(y -> y[i], y; ignore_equal_values = true) for i in 1:2]

    xs, ys = [@lift(round.((w / $Lx, h / $Ly); digits = 2)[i]) for i in 1:2]

    # add center button
    centre = Button(fig[2,3], label = "Center", textsize = 0.76t)
    on(centre.clicks) do _
        xc[] = div(w, 2)
        yc[] = div(h, 2)
        nothing
    end

    # slider callback functions under lock
    function f1(x)
        off(of2)
        set_close_to!(slidery, round(Int, Lx[] / AR[]))
        of2 = on(f2, Ly)
        return
    end

    function f2(y)
        off(of1)
        set_close_to!(sliderx, round(Int, Ly[] * AR[]))
        of1 = on(f1, Lx)
        return
    end

    # set lock toggle
    arlock = Toggle(fig[3,3])
    Label(fig[3,2], "Lock AR", halign = :right, tellwidth = false,
          textsize = 0.76t)

    on(arlock.active) do active
        if active
            AR = Lx[] / Ly[]
            of1 = on(f1, Lx)
            of2 = on(f2, Ly)
        else
            off(of1)
            off(of2)
        end
        nothing
    end

    # set the fft image
    fftaxis = Axis(fig[1,2], aspect = DataAspect(),
                   title = @lift("Effective Resolution: $($Lx) x $($Ly)"),
                   subtitle = @lift("x-scale: $($xs), y-scale: $($ys)"),
                   titlesize = t, subtitlesize = t)
    deregister_interaction!(fftaxis, :rectanglezoom)
    image!(fftaxis, S, colormap = :tokyo, interpolate = false)
    rect = select_rectangle(fftaxis)
    DataInspector(fig)

    # draw rectangle
    r4 = @lift(Rect($x1, $y1, $Lx, $Ly))
    R = lines!(fftaxis, r4, linestyle = :dot,
               color = RGBf(0, 1, 0.38), linewidth = 7)

    # selection observables
    Sv = @lift(vec(S[$x1+1:$x2, $y1+1:$y2]))
    on(Sv) do vs
        if allequal(vs)
            push!(Sv[], isempty(vs) ? 0 : vs[1] + 1)
        end
        nothing
    end
    n, s6 = stats(Sv[])
    n = Observable(n)
    Svs = Observable(Sv[])
    Svlt = async_latest(Sv)

    # add interactive elements to histogram
    HGrid = fig[1,1] = GridLayout(tellwidth = false, tellheight = false,
                                  valign = :top, halign = :right)
    HGrid[1,1] = Label(fig, "")
    HGrid[2,1] = Label(fig, "Live", textsize = t)
    HGrid[2,2] = live = Toggle(fig)
    HGrid[3,1] = Label(fig, @lift($(live.active) || $Svs == $Sv ? "\u2713" : ""),
                       textsize = 2t)
    HGrid[3,2] = update = Button(fig, label = "Update", textsize = t)
    HGrid[3,3] = Label(fig, "")

    # set the histogram
    Haxis = Axis(fig[1,1], title = "$w x $h pixels", titlesize = t,
                 xlabel = "pixel value [0-255]", ylabel = "pixel count")
    deregister_interaction!(Haxis, :rectanglezoom)
    H = hist!(Haxis, Svs, bins = n, color = RGBf(0, 0.5, 0.8))
    # abuse the legend since it's easier than making a custom text! box
    legend = axislegend(Haxis, [MarkerElement(marker = "", color = :transparent)],
                        [s6], titlesize = t, labelsize = t, position = :rc)

    # interval change callback for related histogram attributes
    function refresh(Svlt = Sv[])
        n[], legend.entrygroups[][1][2][1].label[] = stats(Svlt)
        Haxis.title[] = "$(Lx[]) x $(Ly[]) pixels"
        reset_limits!(Haxis)
        return
    end

    # console printing
    function info()
        println("\nEffective Resolution: $(Lx[]) x $(Ly[])")
        println("x-scale: $(xs[])\ny-scale: $(ys[])")
        return
    end

    # histogram callbacks
    on(update.clicks) do _
        Svs[] = Sv[]
        refresh()
        info()
        nothing
    end

    on(live.active) do active
        delete!(Haxis, H)
        if active
            H = hist!(Haxis, Svlt, bins = n, color = RGBf(0, 0.8, 0.3))
            ov = on(refresh, Svlt)
        else
            off(ov)
            H = hist!(Haxis, Svs, bins = n, color = RGBf(0, 0.5, 0.8))
            Svs[] = Sv[]
        end
        refresh()
    end

    # rectangle selection callback
    on(rect) do v
        # bounds
        x01, y01 = round.(Int, v.origin)
        Lx0, Ly0 = round.(Int, v.widths)
        x01 = clamp(x01, 0, w - 1)
        y01 = clamp(y01, 0, h - 1)
        x02 = clamp(Lx0 + x01, x01 + 1, w)
        y02 = clamp(Ly0 + y01, y01 + 1, h)
        Lx0, Ly0 = x02 - x01, y02 - y01

        # set
        arlock.active[] && off(of1) && off(of2)
        xc[] = x01 + div(Lx0, 2)
        set_close_to!(sliderx, Lx0)
        yc[] = y01 + div(Ly0, 2)
        set_close_to!(slidery, Ly0)
        AR = Lx[] / Ly[]
        Svs[] = Sv[]
        refresh()
        if arlock.active[]
            of1 = on(f1, Lx)
            of2 = on(f2, Ly)
        end
        info()
        nothing
    end

    # render
    set_window_config!(title = "fresdet", focus_on_show = true)
    screen = display(fig)

    script && wait(screen)

    return fig, S

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
    fig, S = fresdet(ARGS[1]);
    empty!(ARGS)
end

print()
