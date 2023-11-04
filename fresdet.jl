# fresdet.jl
# Version 0.10.1
# 2023-11-4
# https://github.com/Sagnac/fresdet

# Simple analysis tool for estimating the original resolution
# of standard upscaled images using ffts and basic statistical metrics.

using Images
using FFTW
using StatsBase
using Printf
using GLMakie
using GLMakie: GLFW.GetPrimaryMonitor, MonitorProperties
using .Makie: Axis, async_latest

function fftimage(file)
    image = file |> load |> rotr90
    image = image .|> Gray .|> Float64
    F2 = image |> fft |> fftshift .|> abs2
    map!(x -> log(x + 1), F2, F2)
    F2 * 255 / maximum(F2) .|> round
end

function stats(v)
    n::Int = maximum(v) - minimum(v) + 1
    c = [f(v) for f in (mean, median, mode, std, skewness, kurtosis)]
    s6 = @sprintf(
             """
             mean: %.3f
             median: %.3f
             mode: %.3f
             st.dev.: %.3f
             skew: %.3f
             kurtosis: %.3f
             """,
             c...
         )
    n, replace(s6, "-" => "\u2212")
end

function fresdet(file; script = !isinteractive())

    S = fftimage(file)
    w, h = size(S)
    monitor_properties = MonitorProperties(GetPrimaryMonitor())
    (; height) = monitor_properties.videomode
    res = (height, height/2)
    t = 0.1 * monitor_properties.dpi[1] # text size
    t2 = 0.76t # smaller text size used for some of the outer widgets
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
    y = lift(ypoints, yc, Ly; ignore_equal_values = true)

    scale = @lift(@sprintf("x-scale: %.2f, y-scale: %.2f", w / $Lx, h / $Ly))

    # add center button
    centre = Button(fig[2,3], label = "Center", fontsize = t2)
    on(centre.clicks) do _
        xc[] = div(w, 2)
        yc[] = div(h, 2)
        Svs[] = Sv[]
        refresh()
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
          fontsize = t2)

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

    # set the strings for the inspector labels
    image_inspector(_, i, p) = @sprintf "%u, %u = %u" i[1] i[2] p[3]
    histogram_inspector(_, _, p) = @sprintf "x: %.0f\ny: %.0f" p[1] p[2]

    # set the fft image
    fftaxis = Axis(fig[1,2], aspect = DataAspect(),
                   title = @lift("Effective Resolution: $($Lx) x $($Ly)"),
                   subtitle = scale,
                   titlesize = t, subtitlesize = t)
    deregister_interaction!(fftaxis, :rectanglezoom)
    image!(fftaxis, S, colormap = :tokyo, interpolate = false,
           inspector_label = image_inspector)
    rect = select_rectangle(fftaxis)
    DataInspector(fig)

    # draw rectangle
    selection = lift((x, y, Lx, Ly) -> (x[1], y[1], Lx, Ly), x, y, Lx, Ly;
                     ignore_equal_values = true)
    r4 = @lift(Rect($selection...))
    R = lines!(fftaxis, r4, linestyle = :dot, inspectable = false,
               color = RGBf(0, 1, 0.38), linewidth = 7)

    # selection observables
    Sv = lift((x, y) -> vec(@view S[x[1]+1:x[2], y[1]+1:y[2]]), x, y;
              ignore_equal_values = true)
    n, s6 = stats(Sv[])
    n = Observable(n)
    Svs = Observable(Sv[])
    Svlt = async_latest(Sv)

    # add interactive elements to histogram
    HGrid = fig[1,1] = GridLayout(tellwidth = false, tellheight = false,
                                  valign = :top, halign = :right)
    HGrid[1,1] = Label(fig, "")
    HGrid[2,1] = Label(fig, "Live", fontsize = t)
    HGrid[2,2] = live = Toggle(fig)
    HGrid[3,1] = Label(fig, @lift($(live.active) || $Svs == $Sv ? "\u2713" : ""),
                       fontsize = 2t)
    HGrid[3,2] = update = Button(fig, label = "Update", fontsize = t)
    HGrid[3,3] = Label(fig, "")

    # set the histogram
    Haxis = Axis(fig[1,1], title = @lift("$($Lx) x $($Ly) pixels"), titlesize = t,
                 xlabel = "pixel value [0-255]", ylabel = "pixel count")
    deregister_interaction!(Haxis, :rectanglezoom)
    H = hist!(Haxis, Svs, bins = n, color = RGBf(0, 0.5, 0.8),
              inspector_label = histogram_inspector)
    # abuse the legend since it's easier than making a custom text! box
    legend = axislegend(Haxis, [MarkerElement(marker = '!', color = :transparent)],
                        [s6], titlesize = t, labelsize = t, position = :rc)

    # interval change callback for related histogram attributes
    function refresh(Svlt = Sv[])
        n[], legend.entrygroups[][1][2][1].label[] = stats(Svlt)
        reset_limits!(Haxis)
        return
    end

    # console printing
    function info()
        println("\nEffective Resolution: $(Lx[]) x $(Ly[])")
        println(scale[])
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
            H = hist!(Haxis, Svlt, bins = n, color = RGBf(0, 0.8, 0.3),
                      inspector_label = histogram_inspector)
            ov = on(refresh, Svlt)
        else
            off(ov)
            H = hist!(Haxis, Svs, bins = n, color = RGBf(0, 0.5, 0.8),
                      inspector_label = histogram_inspector)
            Svs[] = Sv[]
        end
        refresh()
    end

    # rectangle selection callback
    on(rect) do v
        # bounds
        x01, y01 = trunc.(Int, v.origin)
        Lx0, Ly0 = ceil.(Int, v.widths)
        x01 = clamp(x01, 0, w - 1)
        y01 = clamp(y01, 0, h - 1)
        x02 = clamp(Lx0 + x01, x01 + 1, w)
        y02 = clamp(Ly0 + y01, y01 + 1, h)
        Lx0, Ly0 = x02 - x01, y02 - y01

        # set
        arlock.active[] && off(of1) && off(of2)
        xc.val = x01 + div(Lx0, 2)
        yc.val = y01 + div(Ly0, 2)
        set_close_to!(sliderx, Lx0)
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

    # fft save button
    savefft = Button(fig[2,1], label = "Save FFT", halign = :right,
                     tellwidth = false, fontsize = t2)
    on(savefft.clicks) do _
        save("FFT.png", clamp01nan!(rotl90(S / 255)))
        println("FFT saved as FFT.png in the present working directory.\n")
        nothing
    end

    # render
    set_window_config!(title = "fresdet", focus_on_show = true)
    screen = display(fig)

    script && wait(screen)

    return fig

end

# script support
if !isempty(ARGS)
    fig, S = fresdet(ARGS[1]);
    empty!(ARGS)
end

print()
