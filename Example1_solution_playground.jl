#=
Example 1 – Interactive g–x diagram: how the phases and the stable assemblage change with P and T
=================================================================================================

Two components A and B, one composition axis X = X_B, and six phases that can
be switched on and off:

    α      solid solution between αA and αB with a Margules (regular-solution)
           parameter W:  g = (1-X) g_αA + X g_αB + RT[X lnX + (1-X) ln(1-X)] + W X(1-X)
           W = W_H - T W_S + P W_V.  If W is large enough (W > 2RT) the g curve
           gets a hump and the solution unmixes: α + α  (a solvus).
    A      pure phase, X = 0   (a low-T polymorph of the α end-member αA)
    A2B    compound of fixed composition, X = 1/3
    AB2    compound of fixed composition, X = 2/3
    B      pure phase, X = 1   (a low-T polymorph of αB)
    L      ideal solution (a "melt") of L_A and L_B (W = 0)

Every phase has g = g0 - s T + v P (end-members and compounds), so moving the
T and P sliders moves every curve and point. The stable assemblage is found
exactly as in Examples 0 and 2–5: by minimising G with linear programming
(gibbs_lp.jl), the solutions being discretised into pseudocompounds. Geometrically
the result is the lower convex hull of everything that is drawn (thick black
line) and the common tangent at the bulk composition (dashed, its intercepts
are μ_A and μ_B).

How to use
----------
* Sliders: T, P, the Margules parameter W_H of α, the bulk composition X_B and
  the number of pseudocompounds per solution.
* Left-click in the g–x diagram to set the bulk composition.
* Check boxes: switch every phase on or off (also all at once). A switched-off
  phase is removed from the LP and from the plot.
* The strip below the g–x diagram shows the stable assemblage for every bulk
  composition at the current P and T (one LP per column), i.e. one horizontal
  line through a T–X diagram.

Suggested exercises: (1) only α on: raise W or lower T until the solvus
appears; (2) α + compounds: find the T where A2B and AB2 stop being stable;
(3) everything on: heat X_B = 0.5 from 800 K to 1400 K and name every
assemblage; (4) repeat at 15 kbar; (5) X_B = 0.05, cool from 900 K: where does
pure A appear, and why does it appear at a lower T than for X_B = 0?

Default data (P = 0): pure A and B below ≈ 750 / 700 K (A ⇄ αA, B ⇄ αB),
A2B + AB2 (with α, A or B) below ≈ 910 K, α + α solvus up to ≈ 1085 K
(crest at X_B = 0.5), one α up to ≈ 1285 K, then α + L and L.

Run from a terminal:      julia Example1_solution_playground.jl
or from the Julia REPL:   include("Example1_solution_playground.jl")
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
const COMPONENTS = ["A", "B"]

"Binary solution: end-members e1 (X = 0) and e2 (X = 1), W = WH - T WS + P WV."
struct Solution
    name  :: String
    e1    :: Phase
    e2    :: Phase
    WH    :: Float64      # [kJ/mol] (for α: start value of the slider)
    WS    :: Float64      # [kJ/mol/K]
    WV    :: Float64      # [kJ/kbar/mol]
    color :: Symbol
end

#                              name   oxide amounts   g0     s       v
const ALPHA = Solution("α", Phase("αA", [1.0, 0.0],    0.0,  0.0,    0.0),
                            Phase("αB", [0.0, 1.0],    0.0,  0.0,    0.0),
                       18.0, 0.0, 0.05, :royalblue3)
const MELT  = Solution("L", Phase("L_A", [1.0, 0.0],  35.0,  35 / 1600, 0.20),   # T_m(A) = 1600 K
                            Phase("L_B", [0.0, 1.0],  30.0,  30 / 1400, 0.10),   # T_m(B) = 1400 K
                       0.0, 0.0, 0.0, :firebrick)                                # ideal
# phases of fixed composition (pure phases and compounds), in the order of the check boxes
const COMPOUNDS = [(Phase("A",   [1.0, 0.0], -0.9, -0.0012, -0.02), :sienna),
                   (Phase("A2B", [2/3, 1/3], -2.0, -0.0012, -0.03), :darkorange2),
                   (Phase("AB2", [1/3, 2/3], -2.0, -0.0012, -0.03), :forestgreen),
                   (Phase("B",   [0.0, 1.0], -0.7, -0.0010, -0.02), :teal)]
# g of the compounds: -2.0 + 0.0012 T - 0.03 P  (stable at low T and high P)
# pure A: -0.9 + 0.0012 T - 0.02 P,  pure B: -0.7 + 0.0010 T - 0.02 P, relative to the
# α end-members (g = 0): low-T polymorphs, A ⇄ αA at 750 K and B ⇄ αB at 700 K (P = 0)

const T_START  = 1000.0                      # [K]
const P_START  = 0.0                         # [kbar]
const X_START  = 0.40
const NPC      = 201                         # pseudocompounds per solution
const TLIM     = (500.0, 1800.0)
const NSTRIP   = 200                         # columns of the assemblage strip

const FIELD_COLORS = ["#a6cee3", "#fdbf6f", "#b2df8a", "#fb9a99", "#cab2d6",
                      "#ffff99", "#8dd3c7", "#d9d9d9", "#bc80bd", "#ccebc5"]

# =============================================================================
# THERMODYNAMICS (independent of the plotting)
# =============================================================================
const SOLUTIONS = [ALPHA, MELT]
# phase "switches", in the order of the check boxes
const PHASE_NAMES = [ALPHA.name; [p.name for (p, _) in COMPOUNDS]; MELT.name]

Wof(s::Solution, WH, T, P) = WH - T * s.WS + P * s.WV

"Continuous g of a solution at X [kJ/mol]."
function gsol(s::Solution, WH, X, T, P)
    g1, g2 = gibbs(s.e1, P, T), gibbs(s.e2, P, T)
    (1 - X) * g1 + X * g2 + RGAS * T * (xlogx(X) + xlogx(1 - X)) + Wof(s, WH, T, P) * X * (1 - X)
end

"Tangent to the solution curve at X: [μ_A, μ_B]."
function musol(s::Solution, WH, X, T, P)
    W = Wof(s, WH, T, P); RT = RGAS * T
    [gibbs(s.e1, P, T) + RT * log(1 - X) + W * X^2, gibbs(s.e2, P, T) + RT * log(X) + W * (1 - X)^2]
end

"""
Build the phase list of the LP from the active phases. Returns the phases and,
for every phase, (:sol, solution index, X) or (:cmp, compound index, X).
"""
function build_phases(active, WH, n)
    phases = Phase[]; kind = Tuple{Symbol,Int,Float64}[]
    if active[1]
        append!(phases, solution_pseudocompounds(ALPHA.name, ALPHA.e1, ALPHA.e2, WH, n; WS = ALPHA.WS, WV = ALPHA.WV))
        append!(kind, [(:sol, 1, x) for x in range(0, 1, length = n)])
    end
    for (c, (p, _)) in enumerate(COMPOUNDS)
        active[1+c] || continue
        push!(phases, p); push!(kind, (:cmp, c, p.comp[2] / sum(p.comp)))
    end
    if active[end]
        append!(phases, solution_pseudocompounds(MELT.name, MELT.e1, MELT.e2, MELT.WH, n; WS = MELT.WS, WV = MELT.WV))
        append!(kind, [(:sol, 2, x) for x in range(0, 1, length = n)])
    end
    return phases, kind
end

"""
Group the stable phases into physical phases: stable pseudocompounds of the same
solution are one phase unless the continuous g curve between them has a hump
(then they are two coexisting phases, e.g. α + α across a solvus).
Returns (name, amount, X, solution index or 0) sorted by X.
"""
function phase_list(phases, kind, r, WH, T, P)
    out = NamedTuple{(:name, :amount, :x, :sol),Tuple{String,Float64,Float64,Int}}[]
    for (i, k) in enumerate(kind)
        k[1] === :cmp && i in r.stable && push!(out, (name = phases[i].name, amount = r.α[i], x = k[3], sol = 0))
    end
    for s in eachindex(SOLUTIONS)
        idx = sort([i for i in r.stable if kind[i][1] === :sol && kind[i][2] == s], by = i -> kind[i][3])
        isempty(idx) && continue
        sol, W = SOLUTIONS[s], s == 1 ? WH : SOLUTIONS[s].WH
        groups = [[idx[1]]]
        for i in idx[2:end]
            a, b = kind[groups[end][end]][3], kind[i][3]
            ga, gb = gsol(sol, W, a, T, P), gsol(sol, W, b, T, P)
            hump = maximum(gsol(sol, W, (1 - t) * a + t * b, T, P) - ((1 - t) * ga + t * gb) for t in 0.1:0.1:0.9)
            hump <= 1e-9 ? push!(groups[end], i) : push!(groups, [i])
        end
        for g in groups
            am = sum(r.α[g])
            push!(out, (name = sol.name, amount = am, x = sum(r.α[i] * kind[i][3] for i in g) / am, sol = s))
        end
    end
    return sort(out, by = p -> p.x)
end

phase_label(p) = p.sol > 0 ? @sprintf("%s(%.2f)", p.name, p.x) : p.name
"Key of an assemblage in a fixed phase order (α, compounds, L), e.g. \"α + L\" (used for the strip colours)."
assemblage_key(pl) = join(sort([p.name for p in pl], by = n -> findfirst(==(n), PHASE_NAMES)), " + ")

# =============================================================================
# APPLICATION
# =============================================================================
function build_app()
    # ---------------------------------------------------------------- state
    T, P, X, WH, n = Ref(T_START), Ref(P_START), Ref(X_START), Ref(ALPHA.WH), Ref(NPC)
    active = fill(true, length(PHASE_NAMES))
    built  = Ref(build_phases(active, WH[], n[]))
    lp     = Ref{Any}(GibbsLP(built[][1], [1 - X[], X[]]))
    color_of = Dict{String,Int}()

    function rebuild_lp!()
        built[] = build_phases(active, WH[], n[])
        lp[] = isempty(built[][1]) ? nothing : GibbsLP(built[][1], [1 - X[], X[]])
    end
    "LP at bulk x; `nothing` if the active phases cannot make that composition."
    function equil(x)
        lp[] === nothing && return nothing
        set_bulk!(lp[], [1 - x, x])
        try
            solve!(lp[], P[], T[])
        catch
            nothing
        end
    end

    # ---------------------------------------------------------------- layout
    fig = Figure(size = (1500, 900), fontsize = 15)
    Label(fig[0, 1:2], "Example 1 – g–x diagram: phases, solutions, solvus and melt as functions of P and T",
          fontsize = 20, font = :bold, word_wrap = true, tellwidth = false)
    left  = fig[1, 1] = GridLayout()
    right = fig[1, 2] = GridLayout()
    colsize!(fig.layout, 2, Fixed(470))

    ax = Axis(left[1, 1]; ylabel = "g [kJ/mol]", xticks = 0:0.1:1,
              title = "g–x diagram  (left-click to set the bulk composition)")
    axS = Axis(left[2, 1]; xlabel = "X_B", yticksvisible = false, yticklabelsvisible = false,
               xticks = 0:0.1:1, title = "stable assemblage at this P and T, for every bulk composition",
               titlesize = 14, titlefont = :regular)
    linkxaxes!(ax, axS)
    rowsize!(left, 2, Fixed(70))
    for a in (ax, axS), name in (:rectanglezoom, :scrollzoom, :dragpan)
        name in keys(interactions(a)) && deregister_interaction!(a, name)
    end
    info = Label(left[3, 1], " ", halign = :left, tellwidth = false, justification = :left)

    sg = SliderGrid(right[1, 1],
        (label = "T [K]",          range = TLIM[1]:1.0:TLIM[2],  startvalue = T_START,  format = "{:.0f}"),
        (label = "P [kbar]",       range = 0.0:0.5:30.0,         startvalue = P_START,  format = "{:.1f}"),
        (label = "W_H of α",       range = -10.0:0.5:40.0,       startvalue = ALPHA.WH, format = "{:.1f} kJ/mol"),
        (label = "bulk X_B",       range = 0.0:0.005:1.0,        startvalue = X_START,  format = "{:.3f}"),
        (label = "pseudocomp.",    range = 11:2:401,             startvalue = NPC,      format = n -> "$n"),
        width = 450, tellheight = true)
    sl_T, sl_P, sl_W, sl_X, sl_n = sg.sliders

    boxes = right[2, 1] = GridLayout(tellwidth = false, halign = :left)
    Label(boxes[1, 1:2], "phases present", font = :bold, halign = :left)
    phase_colors = [ALPHA.color, last.(COMPOUNDS)..., MELT.color]
    descr = ["α  solid solution (Margules W)",
             [(x = p.comp[2] / sum(p.comp); x in (0, 1) ? @sprintf("%s  pure phase, X_B = %.0f", p.name, x) :
                                                         @sprintf("%s  compound, X_B = %.2f", p.name, x))
              for (p, _) in COMPOUNDS]...,
             "L  ideal melt"]
    cbs = map(eachindex(PHASE_NAMES)) do i
        cb = Checkbox(boxes[i+1, 1], checked = true, checkboxcolor_checked = phase_colors[i])
        Label(boxes[i+1, 2], descr[i], halign = :left, color = phase_colors[i])
        cb
    end
    cb_all = Checkbox(boxes[length(PHASE_NAMES)+2, 1], checked = true)
    Label(boxes[length(PHASE_NAMES)+2, 2], "all phases", halign = :left, font = :bold)
    Label(boxes[length(PHASE_NAMES)+3, 1:2], "display", font = :bold, halign = :left)
    cb_hull = Checkbox(boxes[length(PHASE_NAMES)+4, 1], checked = true)
    Label(boxes[length(PHASE_NAMES)+4, 2], "G of the system (lower convex hull)", halign = :left)
    cb_tan = Checkbox(boxes[length(PHASE_NAMES)+5, 1], checked = true)
    Label(boxes[length(PHASE_NAMES)+5, 2], "tangent at the bulk (μ_A, μ_B)", halign = :left)
    cb_pts = Checkbox(boxes[length(PHASE_NAMES)+6, 1], checked = false)
    Label(boxes[length(PHASE_NAMES)+6, 2], "show the pseudocompounds", halign = :left)
    status = Label(right[3, 1], " ", halign = :left, tellwidth = false, color = :gray30, justification = :left)
    Label(right[4, 1], " ", tellheight = false)                         # spacer

    # ------------------------------------------------------------ drawing
    "Assemblage strip: one LP per column."
    function draw_strip!()
        empty!(axS)
        edges = range(0.0, 1.0, length = NSTRIP + 1)
        xs = (edges[1:end-1] .+ edges[2:end]) ./ 2
        t0 = time()
        keys_ = map(xs) do x
            r = equil(x)
            r === nothing ? "not possible" : assemblage_key(phase_list(built[][1], built[][2], r, WH[], T[], P[]))
        end
        dt = time() - t0
        for k in keys_
            get!(color_of, k, mod1(length(color_of) + 1, length(FIELD_COLORS)))
        end
        cidx = [k == "not possible" ? 0 : color_of[k] for k in keys_]
        cmap = vcat([RGBf(1, 1, 1)], [parse(RGBf, c) for c in FIELD_COLORS])
        heatmap!(axS, collect(edges), [0.0, 1.0], reshape(Float64.(cidx), :, 1);
                 colormap = cgrad(cmap, length(cmap), categorical = true), colorrange = (-0.5, length(FIELD_COLORS) + 0.5))
        # labels on runs that are wide enough
        i = 1
        while i <= NSTRIP
            j = i
            while j < NSTRIP && keys_[j+1] == keys_[i]; j += 1; end
            if (j - i + 1) / NSTRIP >= 0.07
                text!(axS, (xs[i] + xs[j]) / 2, 0.5; text = keys_[i], align = (:center, :center),
                      fontsize = 13, font = :bold, color = :gray15)
            end
            i = j + 1
        end
        vlines!(axS, [X[]]; color = :black, linewidth = 2)
        limits!(axS, -0.04, 1.04, 0, 1)
        return dt
    end

    "The g–x diagram at the current P, T, W and bulk composition."
    function draw_gx!()
        empty!(ax)
        phases, kind = built[]
        lines!(ax, [0.0, 1.0], [0.0, 0.0]; color = (:gray50, 0.4), linewidth = 0.8)
        if isempty(phases)
            info.text[] = "No phase is active."
            return
        end
        g = gibbs(phases, P[], T[])
        allg = Float64[]
        # solution curves
        for (s, sol) in enumerate(SOLUTIONS)
            idx = [i for (i, k) in enumerate(kind) if k[1] === :sol && k[2] == s]
            isempty(idx) && continue
            xs = [kind[i][3] for i in idx]
            lines!(ax, xs, g[idx]; color = sol.color, linewidth = 2.5)
            cb_pts.checked[] && scatter!(ax, xs, g[idx]; markersize = 4, color = sol.color)
            append!(allg, g[idx])
            m = idx[clamp(round(Int, (s == 1 ? 0.15 : 0.85) * length(idx)), 1, length(idx))]   # label position
            text!(ax, kind[m][3], g[m]; text = sol.name, color = sol.color, fontsize = 20, font = :bold,
                  offset = (0, 10), align = (:center, :bottom))
        end
        # compounds
        for (i, k) in enumerate(kind)
            k[1] === :cmp || continue
            c = COMPOUNDS[k[2]][2]
            scatter!(ax, [k[3]], [g[i]]; marker = :rect, markersize = 15, color = :white, strokecolor = c, strokewidth = 2.5)
            text!(ax, k[3], g[i]; text = phases[i].name, color = c, fontsize = 16, font = :bold,
                  offset = (0, 12), align = (k[3] == 0 ? :left : k[3] == 1 ? :right : :center, :bottom))
            push!(allg, g[i])
        end
        # lower convex hull = G of the system for every bulk composition
        xk = [k[3] for k in kind]
        hx, hg = lower_hull(xk, g)
        cb_hull.checked[] && lines!(ax, hx, hg; color = :black, linewidth = 4.5, alpha = 0.55)

        r = equil(X[])
        lo, hi = extrema(allg); pad = 0.12 * max(hi - lo, 1.0)
        if r === nothing
            vlines!(ax, [X[]]; color = :gray50, linestyle = :dot)
            limits!(ax, -0.04, 1.04, lo - pad, hi + pad)
            info.text[] = @sprintf("T = %.0f K, P = %.1f kbar, bulk X_B = %.3f:  the active phases cannot make this bulk composition",
                                   T[], P[], X[])
            return
        end
        pl = phase_list(phases, kind, r, WH[], T[], P[])
        # μ: a single solution phase has a unique tangent; otherwise the LP duals
        μ = if length(pl) == 1 && pl[1].sol > 0 && 0 < pl[1].x < 1
            s = pl[1].sol; musol(SOLUTIONS[s], s == 1 ? WH[] : SOLUTIONS[s].WH, pl[1].x, T[], P[])
        else
            r.μ
        end
        if cb_tan.checked[]
            lines!(ax, [0.0, 1.0], μ; color = :black, linestyle = :dash, linewidth = 1.5)
            text!(ax, [0.0, 1.0], μ; text = [rich("μ", subscript("A"), @sprintf(" = %.2f", μ[1])),
                                              rich("μ", subscript("B"), @sprintf(" = %.2f", μ[2]))],
                  align = [(:left, :top), (:right, :top)], offset = [(4, -6), (-4, -6)], fontsize = 14)
        end
        # stable phases
        for p in pl
            c = p.sol > 0 ? SOLUTIONS[p.sol].color : COMPOUNDS[findfirst(q -> q[1].name == p.name, COMPOUNDS)][2]
            gp = (1 - p.x) * μ[1] + p.x * μ[2]
            scatter!(ax, [p.x], [gp]; markersize = p.sol > 0 ? 16 : 17, marker = p.sol > 0 ? :circle : :rect,
                     color = c, strokecolor = :black, strokewidth = 1.5)
        end
        vlines!(ax, [X[]]; color = :gray50, linestyle = :dot)
        scatter!(ax, [X[]], [r.G]; marker = :star5, markersize = 24, color = :gold, strokewidth = 1.5, strokecolor = :black)
        lo = min(lo, minimum(μ)); hi = max(hi, maximum(μ)); pad = 0.12 * max(hi - lo, 1.0)
        limits!(ax, -0.04, 1.04, lo - pad, hi + pad)
        ax.title = @sprintf("g–x at T = %.0f K, P = %.1f kbar, W(α) = %.1f kJ/mol  (left-click to set the bulk)",
                            T[], P[], Wof(ALPHA, WH[], T[], P[]))
        crit = Wof(ALPHA, WH[], T[], P[]) / (2RGAS)
        info.text[] = "stable: " * join([@sprintf("%s %.0f%%", phase_label(p), 100p.amount) for p in pl], " + ") *
                      @sprintf("     G = %.3f kJ/mol,  μ_A = %.3f,  μ_B = %.3f", r.G, μ...) *
                      (length(pl) == 1 && pl[1].sol == 0 ? "   (one compound: μ not unique)" : "") *
                      (active[1] ? @sprintf("\nα: W = %.1f kJ/mol → W/2R = %.0f K (the solvus crest if α were alone)", Wof(ALPHA, WH[], T[], P[]), crit) : "")
    end

    function update!(; strip = true)
        draw_gx!()
        dt = strip ? draw_strip!() : 0.0
        np = length(built[][1])
        strip && (status.text[] = @sprintf("%d phases in the LP\nstrip: %d LPs in %.2f s", np, NSTRIP, dt))
    end

    # ------------------------------------------------------------ callbacks
    on(v -> (T[] = v; update!()), sl_T.value)
    on(v -> (P[] = v; update!()), sl_P.value)
    on(sl_W.value) do v
        WH[] = v; rebuild_lp!(); update!()
    end
    on(sl_n.value) do v
        n[] = round(Int, v); rebuild_lp!(); update!()
    end
    on(v -> (X[] = v; update!(strip = false); draw_strip!()), sl_X.value)
    syncing = Ref(false)
    for (i, cb) in enumerate(cbs)
        on(cb.checked) do v
            active[i] = v
            syncing[] && return
            rebuild_lp!(); update!()
        end
    end
    on(cb_all.checked) do v
        syncing[] = true
        foreach(cb -> cb.checked[] = v, cbs)
        syncing[] = false
        rebuild_lp!(); update!()
    end
    for cb in (cb_hull, cb_tan, cb_pts)
        on(_ -> update!(strip = false), cb.checked)
    end
    register_interaction!(ax, :pick_x) do event::MouseEvent, axis
        event.type === MouseEventTypes.leftclick && set_close_to!(sl_X, clamp(event.data[1], 0.0, 1.0))
        return Consume(false)
    end

    # ------------------------------------------------------- initialisation
    update!()
    return fig, (set_T = v -> set_close_to!(sl_T, v), set_P = v -> set_close_to!(sl_P, v),
                 set_W = v -> set_close_to!(sl_W, v), set_X = v -> set_close_to!(sl_X, v),
                 checkboxes = cbs, all = cb_all, show_points = cb_pts)
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
