#=
Example 0 – Interactive phase equilibria with three phases of fixed composition
================================================================================

Three phases built from two components (oxides) A and B (cf. LinProg1.ipynb):

    A  (X_B = 0)      AB (X_B = 0.5)      B  (X_B = 1)

with Gibbs energies  g = g0 - s T + v P.  At every (P,T) the stable assemblage
is found by minimising the Gibbs energy of the system with linear programming
(see gibbs_lp.jl).

How to use
----------
* At start-up the program is *initialised*: the LP is solved on a P–T grid and
  the field boundaries are located by bisection between neighbouring grid
  nodes that have different assemblages (black dots).
* Left-click anywhere in the P–T diagram to pick a point (a cross appears),
  then press "Compute g–x". The g–x diagram for that point appears in the
  next free panel on the right (up to 3; the oldest is dropped after that).
* Turn on "compute on click" to skip the button.
* Move the "Bulk X_B" slider to change the bulk composition of the system.
* "Show fields" / "analytic line" can be switched off so that students first
  have to discover the topology by themselves.

Run from a terminal:      julia Example0_fixed_composition.jl
or from the Julia REPL:   include("Example0_fixed_composition.jl")
The window opens smaller if it would not fit on the screen (window_size.jl).
=#

using Pkg
Pkg.activate(@__DIR__)

include(joinpath(@__DIR__, "gibbs_lp.jl"))
using GLMakie
include(joinpath(@__DIR__, "window_size.jl"))       # fit_to_screen!: a window that fits on the screen

# =============================================================================
# USER INPUT
# =============================================================================
const COMPONENTS = ["A", "B"]                        # the two oxides
#                   name   oxide amounts    g0      s       v
const PHASES = [Phase("A",  [1.0, 0.0],   10.0,  0.020,  1.00),
                Phase("B",  [0.0, 1.0],   14.0,  0.016,  0.80),
                Phase("AB", [0.5, 0.5],   16.0,  0.028,  1.15)]
# units: g0 [kJ/mol], s [kJ/mol/K], v [kJ/kbar/mol = J/bar]
# Reaction A + B = AB (per mole of oxide: ½A + ½B = AB):
#   Δg0 = 4, Δs = 0.010, Δv = 0.25  ->  P = 0.04 T - 16   (AB stable at high T, low P)

const TLIM      = (200.0, 1200.0)                    # temperature range [K]
const PLIM      = (0.0, 40.0)                        # pressure range [kbar]
const XB_START  = 0.25                               # initial bulk composition X_B
const NGRID     = 25                                 # initial grid nodes per axis
const REACTIONS = [(1, 2, 3)]                        # analytic lines: phase 3 on tie-line 1–2
const NSLOTS    = 3                                  # number of g–x panels

const SLOT_COLORS  = [:royalblue3, :darkorange2, :forestgreen]
const FIELD_COLORS = ["#a6cee3", "#fdbf6f", "#b2df8a", "#fb9a99",
                      "#cab2d6", "#ffff99", "#8dd3c7", "#d9d9d9"]

# =============================================================================
# APPLICATION
# =============================================================================
function build_app()
    names = [p.name for p in PHASES]
    xm    = xcomp(PHASES)
    Acomp = comp_matrix(PHASES)
    lp    = GibbsLP(PHASES, [1 - XB_START, XB_START])

    # ---------------------------------------------------------------- state
    xb        = Ref(XB_START)
    points    = Tuple{Float64,Float64}[]             # (T,P) of the computed slots
    pending   = Ref{Union{Nothing,Tuple{Float64,Float64}}}(nothing)
    color_of  = Dict{Vector{Int},Int}()              # assemblage -> field colour (kept stable)
    fieldplts = Any[]                                # plots of the current grid (re-created)

    # ---------------------------------------------------------------- layout
    fig = Figure(size = (1500, 900), fontsize = 15)
    Label(fig[0, 1:2], "Example 0 – three phases of fixed composition (A, AB, B): " *
          "Gibbs energy minimisation by linear programming", fontsize = 20, font = :bold, word_wrap = true, tellwidth = false)

    left  = fig[1, 1] = GridLayout()
    right = fig[1, 2] = GridLayout()
    colsize!(fig.layout, 1, Relative(0.52))

    axPT = Axis(left[1, 1]; xlabel = "T [K]", ylabel = "P [kbar]",
                title = "P–T diagram  (left-click to pick a point)")
    limits!(axPT, TLIM..., PLIM...)
    for name in (:rectanglezoom, :scrollzoom, :dragpan)   # keep the diagram fixed
        name in keys(interactions(axPT)) && deregister_interaction!(axPT, name)
    end

    # controls
    ctrl = left[2, 1] = GridLayout(tellwidth = false)
    sg = SliderGrid(ctrl[1, 1:6],
        (label = "Bulk X_B",   range = 0.0:0.01:1.0, startvalue = XB_START, format = "{:.2f}"),
        (label = "Grid nodes", range = 5:1:100,      startvalue = NGRID,    format = n -> "$n × $n"))
    sl_xb, sl_grid = sg.sliders

    btn_init  = Button(ctrl[2, 1], label = "Initialise grid")
    btn_comp  = Button(ctrl[2, 2], label = "Compute g–x", buttoncolor = RGBf(0.80, 0.88, 1.0))
    btn_clear = Button(ctrl[2, 3], label = "Clear points")
    tg_fields = Toggle(ctrl[3, 1], active = true);  Label(ctrl[3, 2], "show fields",       halign = :left)
    tg_line   = Toggle(ctrl[3, 3], active = true);  Label(ctrl[3, 4], "analytic line",     halign = :left)
    tg_auto   = Toggle(ctrl[3, 5], active = false); Label(ctrl[3, 6], "compute on click",  halign = :left)

    status = Label(left[3, 1], "", halign = :left, tellwidth = false, color = :gray25)
    cursor = Label(left[4, 1], " ", halign = :left, tellwidth = false, color = :gray25)
    rowsize!(left, 1, Relative(0.72))

    # g–x panels
    axGX = [Axis(right[k, 1]; ylabel = "g [kJ/mol]",
                 xticks = (xm, names), xlabel = k == NSLOTS ? "X_B" : "")
            for k in 1:NSLOTS]

    # ------------------------------------------------- static P–T decorations
    Tline = collect(range(TLIM..., length = 200))
    for (i, j, k) in REACTIONS
        lines!(axPT, Tline, reaction_line(PHASES, i, j, k, Tline);
               color = :black, linestyle = :dash, linewidth = 1.5, visible = tg_line.active)
    end

    # pending point (cross) and computed points (numbered markers)
    pend_pos = Observable(Point2f(NaN, NaN)); pend_vis = Observable(false)
    scatter!(axPT, pend_pos; marker = :xcross, markersize = 20, color = :black, visible = pend_vis)
    slot_pos = [Observable(Point2f(NaN, NaN)) for _ in 1:NSLOTS]
    slot_vis = [Observable(false) for _ in 1:NSLOTS]
    for k in 1:NSLOTS
        scatter!(axPT, slot_pos[k]; marker = :xcross, markersize = 20,
                 color = SLOT_COLORS[k], visible = slot_vis[k])
        text!(axPT, slot_pos[k]; text = string(k), offset = (8, 6), fontsize = 20,
              font = :bold, color = SLOT_COLORS[k], visible = slot_vis[k])
    end

    # ------------------------------------------------------------ functions
    "Initialise: scan the grid, refine the boundaries, draw the fields."
    function init_grid!()
        n   = round(Int, sl_grid.value[])
        Ts  = range(TLIM..., length = n)
        Ps  = range(PLIM..., length = n)
        t0  = time()
        labels, asm = scan_grid(lp, Ts, Ps)
        Tb, Pb      = refine_boundaries(lp, Ts, Ps, labels)
        dt  = time() - t0

        # stable colour per assemblage, independent of discovery order
        for a in asm
            get!(color_of, a, mod1(length(color_of) + 1, length(FIELD_COLORS)))
        end
        cidx = [color_of[asm[l]] for l in labels]

        foreach(p -> delete!(axPT, p), fieldplts); empty!(fieldplts)
        cmap = cgrad(FIELD_COLORS, length(FIELD_COLORS), categorical = true)
        hm = heatmap!(axPT, Ts, Ps, cidx; colormap = cmap, alpha = 0.55,
                      colorrange = (0.5, length(FIELD_COLORS) + 0.5), visible = tg_fields.active)
        translate!(hm, 0, 0, -10)
        nodes = scatter!(axPT, [Point2f(T, P) for T in Ts for P in Ps];
                         markersize = 2.5, color = (:black, 0.25), visible = tg_fields.active)
        translate!(nodes, 0, 0, -9)
        bnd = scatter!(axPT, Tb, Pb; markersize = 6, color = :black, visible = tg_fields.active)
        translate!(bnd, 0, 0, -8)
        push!(fieldplts, hm, nodes, bnd)

        # label every field at the centre of its grid nodes
        for (id, a) in enumerate(asm)
            idx = findall(==(id), labels)
            pos = Point2f(sum(Ts[I[1]] for I in idx) / length(idx),
                          sum(Ps[I[2]] for I in idx) / length(idx))
            push!(fieldplts, text!(axPT, pos; text = assemblage_name(PHASES, a),
                  align = (:center, :center), fontsize = 20, font = :bold,
                  color = :gray20, visible = tg_fields.active))
        end
        nlp = n^2 + 15 * length(Tb)
        status.text[] = @sprintf("Initialised: %d × %d grid, %d boundary points, ≈ %d LPs in %.2f s  |  bulk X_B = %.2f",
                                 n, n, length(Tb), nlp, dt, xb[])
    end

    "Draw the g–x diagram of slot k for the point (T,P)."
    function draw_slot!(k, T, P)
        ax = axGX[k]
        empty!(ax)
        c   = SLOT_COLORS[k]
        r   = solve!(lp, P, T)
        g   = r.g
        hx, hg = lower_hull(xm, g)
        μ   = r.μ

        vlines!(ax, [xb[]]; color = :gray50, linestyle = :dot)                           # bulk
        lines!(ax, [0.0, 1.0], μ; color = c, linestyle = :dash, linewidth = 1.5)          # tangent (μ)
        lines!(ax, hx, hg; color = c, linewidth = 3)                                    # G of system
        st  = falses(length(g)); st[r.stable] .= true
        scatter!(ax, xm, g; markersize = 16, strokewidth = 2, strokecolor = c,
                 color = [s ? c : :white for s in st])
        text!(ax, xm, g; text = names, offset = (0, 12), align = (:center, :bottom), fontsize = 15)
        scatter!(ax, [xb[]], [r.G]; marker = :star5, markersize = 22, color = :gold, strokewidth = 1.5, strokecolor = :black)
        text!(ax, [0.0, 1.0], μ; text = [rich("μ", subscript(COMPONENTS[1]), @sprintf(" = %.2f", μ[1])),
                                          rich("μ", subscript(COMPONENTS[2]), @sprintf(" = %.2f", μ[2]))],
              align = [(:left, :top), (:right, :top)], offset = [(4, -6), (-4, -6)], fontsize = 13, color = c)

        lo, hi = extrema(vcat(g, μ)); pad = max(hi - lo, 1e-3)
        limits!(ax, -0.06, 1.06, lo - 0.2pad, hi + 0.25pad)
        amounts = join([@sprintf("%s (%.2f)", names[i], r.α[i]) for i in r.stable], " + ")
        note = ""
        if length(r.stable) < length(COMPONENTS)
            # fewer phases than components: the tangent (μ) is not unique, the LP shows one of them
            note = "   |   1 phase: μ not unique (any line through it below the others)"
        else
            # metastable phases touching the (unique) tangent: we sit on a univariant reaction
            on_line = [i for i in eachindex(g) if !st[i] && abs(g[i] - dot(Acomp[:, i], μ)) < 1e-6 * (1 + abs(g[i]))]
            isempty(on_line) ||
                (note = "   |   on reaction: " * join(names[sort(vcat(r.stable, on_line))], ", ") * " coexist")
        end
        ax.title    = @sprintf("%d)  T = %.0f K,  P = %.1f kbar", k, T, P)
        ax.subtitle = "stable: " * amounts * @sprintf("   |   G = %.2f kJ/mol", r.G) * note
        ax.titlecolor = c
    end

    "Re-draw all g–x panels and the numbered markers."
    function redraw_slots!()
        for k in 1:NSLOTS
            if k <= length(points)
                T, P = points[k]
                draw_slot!(k, T, P)
                slot_pos[k][] = Point2f(T, P); slot_vis[k][] = true
            else
                empty!(axGX[k])
                axGX[k].title = "$k)  – empty –"; axGX[k].titlecolor = SLOT_COLORS[k]
                axGX[k].subtitle = "pick a point in the P–T diagram, then press Compute g–x"
                limits!(axGX[k], -0.06, 1.06, 0, 1)
                slot_vis[k][] = false
            end
        end
    end

    "Move the pending point into the next g–x slot."
    function compute!()
        if pending[] === nothing
            status.text[] = "First left-click on the P–T diagram to choose a point."
            return
        end
        push!(points, pending[])
        length(points) > NSLOTS && popfirst!(points)
        pending[] = nothing; pend_vis[] = false
        redraw_slots!()
    end

    "Programmatic access (used for testing / screenshots)."
    function add_point!(T, P)
        pending[] = (T, P); compute!()
    end

    # ------------------------------------------------------------ callbacks
    register_interaction!(axPT, :pick_point) do event::MouseEvent, axis
        if event.type === MouseEventTypes.leftclick
            T = clamp(Float64(event.data[1]), TLIM...)
            P = clamp(Float64(event.data[2]), PLIM...)
            pending[] = (T, P)
            pend_pos[] = Point2f(T, P); pend_vis[] = true
            r = solve!(lp, P, T)
            status.text[] = @sprintf("Picked T = %.0f K, P = %.1f kbar  (%s)  – press Compute g–x",
                                     T, P, assemblage_name(PHASES, r.stable))
            tg_auto.active[] && compute!()
        end
        return Consume(false)
    end

    on(events(fig).mouseposition) do _
        if is_mouseinside(axPT.scene)
            T, P = mouseposition(axPT)
            if TLIM[1] <= T <= TLIM[2] && PLIM[1] <= P <= PLIM[2]
                r = solve!(lp, P, T)
                cursor.text[] = @sprintf("cursor: T = %.0f K, P = %.1f kbar  →  %s",
                                         T, P, assemblage_name(PHASES, r.stable))
            end
        end
        return Consume(false)
    end

    on(_ -> compute!(), btn_comp.clicks)
    on(_ -> init_grid!(), btn_init.clicks)
    on(btn_clear.clicks) do _
        empty!(points); pending[] = nothing; pend_vis[] = false
        redraw_slots!()
    end
    on(sl_xb.value) do v                # new bulk composition: everything is recomputed
        xb[] = v
        set_bulk!(lp, [1 - v, v])
        init_grid!()
        redraw_slots!()
    end

    # ------------------------------------------------------- initialisation
    init_grid!()
    redraw_slots!()
    set_bulk_slider!(v) = set_close_to!(sl_xb, v)
    return fig, (add_point = add_point!, init_grid = init_grid!, set_bulk = set_bulk_slider!)
end

# =============================================================================
# RUN
# =============================================================================
if abspath(PROGRAM_FILE) == (@__FILE__) || isinteractive()
    fig, app = build_app()
    fit_to_screen!(fig)                 # smaller than designed if the screen is small (see window_size.jl)
    screen = display(fig)
    isinteractive() || wait(screen)     # keep the window open when run from a terminal
end
