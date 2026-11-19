#=
make_screenshots.jl – regenerate the annotated screenshots used in guide.tex

    julia make_screenshots.jl 0      # one example (0 … 6)
    julia make_screenshots.jl all    # all examples, each in its own Julia process

Every example defines the same global names, so each one is loaded in a separate
Julia process. The program builds the GUI without showing it, puts it into a
useful state, draws numbered red markers on the widgets (positions are taken
from the widgets themselves) and saves the picture to figures/.
=#

const ROOT = normpath(joinpath(@__DIR__, ".."))
const FIGS = joinpath(@__DIR__, "figures")

if length(ARGS) == 1 && ARGS[1] == "all"
    for k in 0:6
        run(`$(Base.julia_cmd()) $(@__FILE__) $k`)
    end
    exit()
end
const EX = parse(Int, ARGS[1])
mkpath(FIGS)

file = Dict(0 => "Example0_fixed_composition.jl", 1 => "Example1_solution_playground.jl",
            2 => "Example2_eutectic_TX.jl",       3 => "Example3_peritectic_TX.jl",
            4 => "Example4_olivine_loop_TX.jl",   5 => "Example5_ternary_DiAnFo.jl",
            6 => "Example6_ternary_CAS.jl")[EX]
include(joinpath(ROOT, file))                    # defines build_app(), loads GLMakie

# ----------------------------------------------------------------------------
# annotation helpers (pixel coordinates of the figure, origin bottom left)
# ----------------------------------------------------------------------------
const MARK = RGBf(0.85, 0.05, 0.05)

blocks(fig, T) = [b for b in fig.content if b isa T]
txt(b) = (s = b isa Button ? b.label[] : b isa Union{Axis,Axis3} ? b.title[] : b.text[];
          s isa AbstractString ? s : string(s))
function findblock(fig, T, s)
    for b in blocks(fig, T)
        occursin(s, txt(b)) && return b
    end
    error("no $T containing \"$s\"")
end
bbox(b) = b.layoutobservables.computedbbox[]
function axbox(ax)                                   # plotting area of an axis
    v = ax.scene.viewport[]
    Rect2f(v.origin..., v.widths...)
end
"Union of several rectangles."
function union_rect(rs...)
    lo = reduce((a, b) -> min.(a, b), [r.origin for r in rs])
    hi = reduce((a, b) -> max.(a, b), [r.origin .+ r.widths for r in rs])
    Rect2f(lo..., (hi .- lo)...)
end
"Pixel position of the data point (x, y) of a 2D axis."
function data2px(ax, x, y)
    v = ax.scene.viewport[]; l = ax.finallimits[]
    Point2f(v.origin .+ (Point2f(x, y) .- l.origin) ./ l.widths .* v.widths)
end

overlay(fig) = (sc = Scene(fig.scene; camera = campixel!, clear = false); sc)

function badge!(sc, p, n)
    s = scatter!(sc, [p]; markersize = 34, color = MARK, strokecolor = :white, strokewidth = 2.5)
    t = text!(sc, [p]; text = string(n), color = :white, font = :bold, fontsize = 20, align = (:center, :center))
    translate!(s, 0, 0, 20); translate!(t, 0, 0, 21)                   # always on top of the frames
end
"Red frame around a rectangle and a numbered badge at its top-left corner."
function mark!(sc, r::Rect2, n; pad = 4, corner = :lt)
    r2 = Rect2f((r.origin .- pad)..., (r.widths .+ 2pad)...)
    o, w = r2.origin, r2.widths
    lines!(sc, Point2f[o, o .+ (w[1], 0), o .+ w, o .+ (0, w[2]), o]; color = MARK, linewidth = 3)
    p = corner === :lt ? Point2f(o[1] - 18, o[2] + w[2]) : corner === :rt ? Point2f(o[1] + w[1] + 18, o[2] + w[2]) :
        corner === :lb ? Point2f(o[1] - 18, o[2]) : Point2f(o[1] + w[1] + 18, o[2])
    fw, fh = sc.viewport[].widths                     # keep the badge inside the figure
    badge!(sc, Point2f(clamp(p[1], 18, fw - 18), clamp(p[2], 18, fh - 18)), n)
end
mark!(sc, b, n; kw...) = mark!(sc, bbox(b), n; kw...)
"Mark a click position inside an axis: ring + badge."
function click!(sc, ax, x, y, n)
    p = data2px(ax, x, y)
    scatter!(sc, [p]; marker = :circle, markersize = 30, color = :transparent, strokecolor = MARK, strokewidth = 3.5)
    badge!(sc, p .+ Point2f(26, 26), n)
end
"Toggle (or checkbox) together with the label to its right."
function with_label(fig, s, T = Toggle)
    lab = findblock(fig, Label, s)
    lb  = bbox(lab)
    c   = lb.origin .+ lb.widths ./ 2
    cands = blocks(fig, T)
    t = cands[argmin([abs(bbox(b).origin[2] + bbox(b).widths[2] / 2 - c[2]) +
                      (bbox(b).origin[1] < lb.origin[1] ? 0.0 : 1e6) +
                      0.2 * abs(lb.origin[1] - bbox(b).origin[1]) for b in cands])]
    union_rect(bbox(t), lb)
end
sliderrow(sg, i) = union_rect(bbox(sg.labels[i]), bbox(sg.sliders[i]), bbox(sg.valuelabels[i]))

function shot(fig, name; annotate = sc -> nothing)
    Makie.update_state_before_display!(fig)
    sc = overlay(fig)
    annotate(sc)
    save(joinpath(FIGS, name), fig; px_per_unit = 1.4)
    empty!(sc)                                        # remove the markers: later screenshots are clean
    println("saved ", name)
end

# ----------------------------------------------------------------------------
# the screenshots of every example
# ----------------------------------------------------------------------------
if EX == 0
    fig, app = build_app()
    app.add_point(500.0, 20.0); app.add_point(900.0, 10.0); app.add_point(700.0, 12.0)
    axPT = findblock(fig, Axis, "P–T diagram")
    sg = only(blocks(fig, SliderGrid))
    shot(fig, "ex0_overview.png"; annotate = sc -> begin
        click!(sc, axPT, 1000.0, 30.0, 1)
        mark!(sc, findblock(fig, Button, "Compute"), 2)
        mark!(sc, union_rect(axbox.(blocks(fig, Axis)[2:end])...), 3; corner = :lt)
        mark!(sc, sliderrow(sg, 1), 4)
        mark!(sc, with_label(fig, "show fields"), 5)
        mark!(sc, findblock(fig, Button, "Initialise"), 6)
        mark!(sc, findblock(fig, Button, "Clear"), 7)
    end)
    app.set_bulk(0.75)
    shot(fig, "ex0_bulk075.png")

elseif EX == 1
    fig, app = build_app()
    ax = findblock(fig, Axis, "g–x")
    sg = only(blocks(fig, SliderGrid))
    shot(fig, "ex1_overview.png"; annotate = sc -> begin
        mark!(sc, sliderrow(sg, 1), 1)
        mark!(sc, sliderrow(sg, 2), 2)
        mark!(sc, sliderrow(sg, 3), 3)
        mark!(sc, sliderrow(sg, 4), 4)
        click!(sc, ax, 0.62, -1.2, 4)
        cbs = [with_label(fig, s, Checkbox) for s in ("α  solid", "A  pure", "A2B", "AB2", "B  pure", "L  ideal", "all phases")]
        mark!(sc, union_rect(cbs...), 5)
        mark!(sc, axbox(findblock(fig, Axis, "stable assemblage")), 6)
        mark!(sc, union_rect([with_label(fig, s, Checkbox) for s in ("G of the system", "tangent at", "show the pseudo")]...), 7)
    end)
    app.set_T(850.0); app.set_X(0.4)
    shot(fig, "ex1_850K.png")
    for i in 2:6; app.checkboxes[i].checked[] = false; end      # α only
    app.set_T(1000.0); app.set_X(0.2)
    shot(fig, "ex1_alpha_only.png")
    for i in 2:6; app.checkboxes[i].checked[] = true; end
    app.set_T(1300.0); app.set_X(0.5)
    shot(fig, "ex1_1300K.png")

elseif EX in (2, 3)
    # add_point(X, T): X is the weight fraction on the axis of the diagram (An in Example 2, SiO2 in Example 3)
    fig, app = build_app()
    if EX == 2
        app.add_point(0.30, 1250.0); app.add_point(0.30, 1300.0); app.add_point(0.60, 1450.0)
    else
        app.add_point(X_EN, 1550.0); app.add_point(X_EN, 1570.0); app.add_point(0.85, 1750.0)
    end
    axTX = findblock(fig, Axis, "T–X diagram")
    sg = only(blocks(fig, SliderGrid))
    gx = [b for b in blocks(fig, Axis) if b !== axTX]
    shot(fig, "ex$(EX)_overview.png"; annotate = sc -> begin
        EX == 2 ? click!(sc, axTX, 0.20, 1500.0, 1) : click!(sc, axTX, 0.62, 1860.0, 1)
        mark!(sc, findblock(fig, Button, "Compute"), 2)
        mark!(sc, union_rect(axbox.(gx)...), 3)
        mark!(sc, union_rect([sliderrow(sg, i) for i in 1:length(sg.sliders)]...), 4)
        mark!(sc, findblock(fig, Button, "Initialise"), 5)
        mark!(sc, union_rect(with_label(fig, "show fields"), with_label(fig, "compute on click")), 6)
    end)

elseif EX == 4
    fig, app = build_app()
    app.add_point(0.5, 1300.0); app.add_point(0.65, 1495.0); app.add_point(0.5, 1700.0)   # X = weight fraction Fa
    app.experiments.checked[] = true
    axTX = findblock(fig, Axis, "T–X diagram")
    sg = only(blocks(fig, SliderGrid))
    shot(fig, "ex4_overview.png"; annotate = sc -> begin
        click!(sc, axTX, 0.25, 1600.0, 1)
        mark!(sc, findblock(fig, Button, "Compute"), 2)
        mark!(sc, union_rect(sliderrow(sg, 2), sliderrow(sg, 3)), 3)
        mark!(sc, with_label(fig, "experiments at 1 bar", Checkbox), 4)
        mark!(sc, with_label(fig, "exact loop"), 5)
    end)
    app.experiments.checked[] = false
    app.set_Ws(80.0); app.add_point(0.5, 1150.0)
    shot(fig, "ex4_W80.png")

elseif EX == 5
    fig, app = build_app()
    axT = findblock(fig, Axis, "isothermal section")
    sg = only(blocks(fig, SliderGrid))
    ax3 = only(blocks(fig, Axis3))
    shot(fig, "ex5_overview.png"; annotate = sc -> begin
        click!(sc, axT, 0.62, 0.30, 1)
        mark!(sc, only(blocks(fig, Textbox)), 2)
        mark!(sc, sliderrow(sg, 1), 3)
        mark!(sc, sliderrow(sg, 2), 4)
        mark!(sc, union_rect(with_label(fig, "isothermal section"), with_label(fig, "reference")), 5)
        mark!(sc, findblock(fig, Button, "Update section"), 6)
        mark!(sc, axbox(ax3), 7)
    end)
    app.set_T(1290.0); app.update_section(); app.toggles.reference.active[] = true
    shot(fig, "ex5_1290_reference.png")
    app.toggles.section.active[] = false
    shot(fig, "ex5_projection.png")

elseif EX == 6
    fig, app = build_app()
    axT = findblock(fig, Axis, "isothermal section")
    sg = only(blocks(fig, SliderGrid))
    ax3 = only(blocks(fig, Axis3))
    shot(fig, "ex6_overview.png"; annotate = sc -> begin
        click!(sc, axT, 0.50, 0.30, 1)                   # tocart(w) = (w_Al2O3 + w_SiO2/2, √3/2 w_SiO2)
        mark!(sc, only(blocks(fig, Textbox)), 2)
        mark!(sc, sliderrow(sg, 1), 3)
        mark!(sc, union_rect(with_label(fig, "isothermal section"), with_label(fig, "reference")), 4)
        mark!(sc, findblock(fig, Button, "Update section"), 5)
        mark!(sc, axbox(ax3), 6)
    end)
    app.toggles.section.active[] = false
    shot(fig, "ex6_projection.png")
    app.toggles.reference.active[] = true
    shot(fig, "ex6_reference.png")
    app.toggles.reference.active[] = false; app.toggles.section.active[] = true
    app.set_T(1300.0); app.set_bulk([0.40, 0.30, 0.30]); app.update_section()
    shot(fig, "ex6_1300.png")
end
