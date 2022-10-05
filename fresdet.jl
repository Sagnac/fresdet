# fresdet.jl
# Version 0.4.0
# 2022-10-4
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
    deregister_interaction!(Haxis, :rectanglezoom)
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
    local sliderx, slidery, Lx, Ly, xs, ys, Svs, H2, refresh, status, HGrid, live
    local Lxp, xc, Lyp, yc, AR, f1, f2, of1, of2, ov2
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

            # draw rectangle
            r4 = @lift(Rect($x1, $y1, $Lx, $Ly))
            R = lines!(fftaxis, r4, linestyle = :dot,
                       color = RGBf(0, 1, 0.38), linewidth = 7)

            # set attribute observables
            xs, ys = [@lift(round.((w / $Lx, h / $Ly); digits = 2)[i]) for i in 1:2]
            connect!(fftaxis.title, @lift("Effective Resolution: $($Lx) x $($Ly)"))
            connect!(fftaxis.subtitle, @lift("x-scale: $($xs), y-scale: $($ys)"))

            # create new histogram on selection
            Sv = @lift(vec(S[$x1+1:$x2, $y1+1:$y2]))
            on(Sv) do vs
                if length(vs) === 1
                    push!(Sv[], vs[1] + 1)
                end
                nothing
            end
            Svs = Observable(Sv[])
            link = @lift($Svs == $Sv ? "\u2713" : "")
            !@isdefined(status) ? (status = link) : connect!(status, link)
            if !@isdefined(live) || !live.active[]
                delete!(Haxis, H)
                H = hist!(Haxis, Svs, color = RGBf(0, 0.5, 0.8))
            else
                delete!(Haxis, H2)
                H2 = hist!(Haxis, Sv, color = RGBf(0, 0.8, 0.3))
                ov2 = onany(refresh, sliderx.interval, slidery.interval)
            end

            # interval change callback for related histogram attributes
            function refresh(_...)
                Svs[] = Sv[]
                H.bins, legend.entrygroups[][1][2][1].label = stats(Sv[])
                @isdefined(H2) && (H2.bins[] = H.bins[])
                Haxis.title = "$(Lx[]) x $(Ly[]) pixels"
                reset_limits!(Haxis)
                return
            end

            refresh()

            # add center button
            centre = Button(fftfig[2,2], label = "Center", textsize = 0.76t)
            on(centre.clicks) do _
                x1c = div(w, 2) - div(Lx[], 2)
                y1c = div(h, 2) - div(Ly[], 2)
                set_close_to!(sliderx, x1c, x1c + Lx[])
                set_close_to!(slidery, y1c, y1c + Ly[])
                xc = x1[] + div(Lx[], 2)
                yc = y1[] + div(Ly[], 2)
                nothing
            end

            # slider callback functions under lock
            function f1(x)
                off(of1)
                if Lx[] != Lxp
                    set_close_to!(sliderx, 2 * xc - x[2], x[2])
                    if ylock.active[]
                        Δy = Int(div(Lx[] / AR, 2))
                        off(of2)
                        set_close_to!(slidery, yc - Δy, yc + Δy)
                        of2 = on(f2, slidery.interval)
                    else
                        AR = Lx[] / Ly[]
                    end
                end
                Lxp, xc = Lx[], x1[] + div(Lx[], 2)
                Lyp, yc = Ly[], y1[] + div(Ly[], 2)
                of1 = on(f1, sliderx.interval)
                return
            end

            function f2(y)
                off(of2)
                if Ly[] != Lyp
                    set_close_to!(slidery, 2 * yc - y[2], y[2])
                    if xlock.active[]
                        Δx = Int(div(Ly[] * AR, 2))
                        off(of1)
                        set_close_to!(sliderx, xc - Δx, xc + Δx)
                        of1 = on(f1, sliderx.interval)
                    else
                        AR = Lx[] / Ly[]
                    end
                end
                Lxp, xc = Lx[], x1[] + div(Lx[], 2)
                Lyp, yc = Ly[], y1[] + div(Ly[], 2)
                of2 = on(f2, slidery.interval)
                return
            end

            # set lock toggles
            xlock, ylock = [Toggle(fftfig[i+2,2]) for i in 1:2]
            [Label(fftfig[i+2,1], s, halign = :right, tellwidth = false,
             textsize = 0.76t) for (i, s) in pairs(("Lock X", "Lock Y"))]

            on(xlock.active) do active
                if active
                    Lxp, xc = Lx[], x1[] + div(Lx[], 2)
                    Lyp, yc = Ly[], y1[] + div(Ly[], 2)
                    AR = Lx[] / Ly[]
                    of1 = on(f1, sliderx.interval)
                else
                    off(of1)
                end
                nothing
            end

            on(ylock.active) do active
                if active
                    Lxp, xc = Lx[], x1[] + div(Lx[], 2)
                    Lyp, yc = Ly[], y1[] + div(Ly[], 2)
                    AR = Lx[] / Ly[]
                    of2 = on(f2, slidery.interval)
                else
                    off(of2)
                end
                nothing
            end

        else

            @isdefined(of1) && (flag1 = off(of1))
            @isdefined(of2) && (flag2 = off(of2))
            set_close_to!(sliderx, x01, x02)
            sliderx.startvalues = (x01, x02)
            set_close_to!(slidery, y01, y02)
            slidery.startvalues = (y01, y02)
            refresh()
            Lxp, xc = Lx[], x01 + div(Lx[], 2)
            Lyp, yc = Ly[], y01 + div(Ly[], 2)
            AR = Lx[] / Ly[]
            @isdefined(flag1) && flag1 && (of1 = on(f1, sliderx.interval))
            @isdefined(flag2) && flag2 && (of2 = on(f2, slidery.interval))

        end

        # add interactive elements to histogram
        if !@isdefined(HGrid)
            HGrid = Hfig[1,1] = GridLayout(tellwidth = false, tellheight = false,
                                           valign = :top, halign = :right)
            HGrid[1,1] = Label(Hfig, "")
            HGrid[2,1] = Label(Hfig, "Live", textsize = t)
            HGrid[2,2] = live = Toggle(Hfig)
            HGrid[3,1] = Label(Hfig, status, textsize = 2t)
            HGrid[3,2] = update = Button(Hfig, label = "Update", textsize = t)
            HGrid[3,3] = Label(Hfig, "")

            on(update.clicks) do _
                refresh()
                info()
                nothing
            end

            on(live.active) do active
                if active
                    delete!(Haxis, H)
                    H2 = hist!(Haxis, Sv, color = RGBf(0, 0.8, 0.3))
                    ov2 = onany(refresh, sliderx.interval, slidery.interval)
                else
                    delete!(Haxis, H2)
                    H = hist!(Haxis, Svs, color = RGBf(0, 0.5, 0.8))
                    off.(ov2)
                end
                refresh()
            end
        end

        # re-open the histogram/stats window on request
        !events(Hfig).window_open[] && (Hscreen = display(Hfig))

        info()

    end

    function info()
        println("\nEffective Resolution: $(Lx[]) x $(Ly[])")
        println("x-scale: $(xs[])\ny-scale: $(ys[])")
        return
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
