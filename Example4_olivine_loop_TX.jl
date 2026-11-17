#=
Example 4 – Interactive T–X diagram of a binary solid solution that melts (olivine, Fo–Fa)
==========================================================================================

Two components, forsterite Mg2SiO4 (Fo) and fayalite Fe2SiO4 (Fa). Two phases,
both SOLUTIONS:

    Ol   olivine (Mg,Fe)2SiO4, solid solution of Fo(s) and Fa(s)
    L    melt, solution of Fo(l) and Fa(l)

    g(x) = (1-x) g_Fo + x g_Fa + n R T [x ln x + (1-x) ln(1-x)] + W x (1-x)

MOLES IN THE CALCULATION, WEIGHT IN THE PLOTS. x = X_Fa is the MOLE fraction of
fayalite, g is in kJ per mole and the linear programme balances moles. The
plots are in WEIGHT per cent, as the published diagram (Bowen & Schairer 1935)
and the other examples:

    composition   wt fraction  w = x M_Fa / [(1-x) M_Fo + x M_Fa]
    phase amounts mass fractions (lever rule in the wt% diagram)
    g, μ          per gram:  g [J/g] = g [kJ/mol] / M̄,  μ_i [J/g] = μ_i [kJ/mol] / M_i

A straight line in (x, g per mole) is also a straight line in (w, g per gram),
so the common tangent stays a straight line in the wt% plot.

with n = 2 mixing sites per formula unit (Mg and Fe mix on the two octahedral
sites), the same for the melt. With W = 0 for both phases (ideal) the result is
the classic melting LOOP: the solidus and the liquidus meet only at the melting
points of the pure end-members, and at every T between them a Mg-rich olivine
coexists with a more Fe-rich melt (common tangent to the two g curves).

Both curves are discretised into pseudocompounds and the same linear programme
as in Examples 0–3 (gibbs_lp.jl) finds the stable assemblage.

Data (illustrative, 1 atm): melting points Fo 1890 °C, Fa 1205 °C;
enthalpies of fusion Fo 142 kJ/mol (Richet et al. 1993), Fa 89.3 kJ/mol
(Stebbins & Carmichael 1984). The solids are the reference (g = 0); the pure
melts follow from  g(l) - g(s) = ΔH (1 - T/T_m) + ΔV P.

For W = 0 the loop has a closed form in mole fractions (dashed black lines, "exact loop"):
    k_i = exp(-[g_i(l) - g_i(s)] / (n R T)),   X_Ol = (1 - k_Fo) / (k_Fa - k_Fo),   X_L = k_Fa X_Ol

Things to try: W of olivine > 0 (non-ideal solid: the loop develops a minimum,
and for large W a solvus Ol + Ol appears below it, turning the loop into a
eutectic-like diagram); W of the melt; pressure; few pseudocompounds.

How to use (as in Example 2)
----------------------------
* Left-click in the T–X diagram to pick a bulk composition and temperature,
  then press "Compute g–x" (or turn on "compute on click"). The g–x diagrams
  are sorted by T (T₁ lowest, at the bottom); the isotherm and the tie-line
  (olivine – melt) are drawn in the T–X diagram.

Run from a terminal:      julia Example4_olivine_loop_TX.jl
or from the Julia REPL:   include("Example4_olivine_loop_TX.jl")
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
const COMPONENTS = ["Fo", "Fa"]
#                       name     comp         g0     s      v
const SOL_ENDS = [Phase("Fo(s)", [1.0, 0.0],  0.0,   0.0,   0.0),     # solids = reference
                  Phase("Fa(s)", [0.0, 1.0],  0.0,   0.0,   0.0)]
# melting data:  ΔH [kJ/mol]   T_m [°C] at 1 atm   dT_m/dP [K/kbar] (illustrative)
const MELTING  = [(ΔH = 142.0, Tm = 1890.0, dTdP = 4.7),               # forsterite
                  (ΔH =  89.3, Tm = 1205.0, dTdP = 6.0)]               # fayalite
const SITES    = 2.0                          # mixing sites per formula unit (both phases)

const SOL_NAME, LIQ_NAME = "Ol", "L"
const SOL_COLOR, LIQ_COLOR = :royalblue3, :firebrick
const WS_START = 0.0                          # Margules W of olivine [kJ/mol]
const WL_START = 0.0                          # Margules W of the melt [kJ/mol]
const NPC      = 201                          # pseudocompounds per solution
const P_START  = 0.0                          # [kbar]
const TLIM     = (1100.0, 2000.0)             # [°C]
const NGRID    = 25
const NSLOTS   = 3

const SLOT_COLORS  = [:darkorange2, :forestgreen, :purple]
const FIELD_COLORS = ["#b2df8a", "#d9cbe6", "#a6cee3", "#fdbf6f",
                      "#fb9a99", "#ffff99", "#8dd3c7", "#d9d9d9"]

# pure melts from the melting data (ΔS = ΔH/T_m, ΔV = ΔS dT_m/dP)
const LIQ_ENDS = [let s = SOL_ENDS[k], m = MELTING[k], Tm = m.Tm + 273.15
                      Phase(COMPONENTS[k] * "(l)", s.comp, s.g0 + m.ΔH, s.s + m.ΔH / Tm,
                            s.v + m.ΔH / Tm * m.dTdP)
                  end for k in 1:2]
const MOLAR_MASS = [140.69, 203.77]                    # Mg2SiO4, Fe2SiO4 [g/mol]: moles <-> weight
const XTICKS = (0:0.2:1, ["Fo", "20", "40", "60", "80", "Fa"])
const XLABEL = "wt% " * COMPONENTS[2]

# experimental data at 1 bar (see the file for the format)
const EXP_FILE   = joinpath(@__DIR__, "data", "olivine_1bar_experiments.csv")

# =============================================================================
# THERMODYNAMICS (independent of the plotting)
# =============================================================================
K(T_C) = T_C + 273.15

# ----- moles <-> weight. The calculation uses the mole fraction x of Fa,
#       the plots use its weight fraction w.
"weight fraction -> mole fraction of the 2nd component"
w2x(w) = (w / MOLAR_MASS[2]) / ((1 - w) / MOLAR_MASS[1] + w / MOLAR_MASS[2])
"mole fraction -> weight fraction of the 2nd component"
x2w(x) = x * MOLAR_MASS[2] / ((1 - x) * MOLAR_MASS[1] + x * MOLAR_MASS[2])
"mass [g] of one mole of a phase (or bulk) of mole fraction x"
molar_mass(x) = (1 - x) * MOLAR_MASS[1] + x * MOLAR_MASS[2]

"All phases: olivine pseudocompounds followed by melt pseudocompounds (compositions in mole fractions)."
system_phases(Ws, Wl, n) = vcat(solution_pseudocompounds(SOL_NAME, SOL_ENDS..., Ws, n; sites = SITES),
                                solution_pseudocompounds(LIQ_NAME, LIQ_ENDS..., Wl, n; sites = SITES))

"""
Solve at bulk composition w (WEIGHT fraction of Fa, as on the diagram), temperature T [°C]
and pressure P [kbar]. The bulk is converted to moles; the result is in moles and kJ/mol.
"""
function equilibrium!(lp::GibbsLP, w, T_C, P)
    x = w2x(w)
    set_bulk!(lp, [1 - x, x])
    return solve!(lp, P, K(T_C))
end

"""
Physical phases: name, amount (MASS fraction of the system), weight fraction w and mole
fraction x of Fa. Neighbouring pseudocompounds are one phase.
"""
function phase_list(phases, r)
    xm = xcomp(phases)
    pl = map(assemblage_groups(phases, r.stable)) do g
        n = sum(r.α[g])                                        # moles of this phase
        x = sum(r.α[i] * xm[i] for i in g) / n
        (name = phases[g[1]].name, mass = n * molar_mass(x), x = x)
    end
    mtot = sum(p.mass for p in pl)
    return [(name = p.name, amount = p.mass / mtot, w = x2w(p.x), x = p.x) for p in pl]
end
phase_label(p) = @sprintf("%s(%.0f wt%% %s)", p.name, 100p.w, COMPONENTS[2])

"Tangent to the continuous g curve of `ends` at MOLE fraction X: [μ_Fo, μ_Fa] in kJ/mol."
function sol_mu(ends, W, X, T_C, P)
    nRT = SITES * RGAS * K(T_C); g = gibbs(ends, P, K(T_C))
    [g[1] + nRT * log(1 - X) + W * X^2, g[2] + nRT * log(X) + W * (1 - X)^2]
end

"""
    exact_loop(P, Ts)

Closed-form solidus and liquidus of the IDEAL loop (W = 0 in both phases):
k_i = exp(-Δg_i / (nRT)) with Δg_i = g_i(l) - g_i(s);  X_Ol = (1-k_Fo)/(k_Fa-k_Fo),
X_L = k_Fa X_Ol. Returns (X_Ol, X_L, T) for every T between the melting points;
X_Ol and X_L are MOLE fractions (use x2w to plot them).
"""
function exact_loop(P, Ts)
    Xs, Xl, To = Float64[], Float64[], Float64[]
    for T in Ts
        nRT = SITES * RGAS * K(T)
        dg  = gibbs(LIQ_ENDS, P, K(T)) .- gibbs(SOL_ENDS, P, K(T))
        kFo, kFa = exp(-dg[1] / nRT), exp(-dg[2] / nRT)
        xs = (1 - kFo) / (kFa - kFo); xl = kFa * xs
        (0 < xs < 1 && 0 < xl < 1) || continue
        push!(Xs, xs); push!(Xl, xl); push!(To, T)
    end
    return Xs, Xl, To
end

"""
    read_experiments(file)

Read the 1-bar experiments (T_C, type, Fa, unit, source). Compositions given in mol%
are converted to the weight fraction of Fa used on the diagram. Returns a vector of
NamedTuples (T, type, X, source) with X = weight fraction Fa; an empty vector if the
file is missing.
"""
function read_experiments(file)
    out = NamedTuple{(:T, :type, :X, :source),Tuple{Float64,String,Float64,String}}[]
    isfile(file) || return out
    header = true
    for line in eachline(file)
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        if header                                   # first non-comment line = column names
            header = false; continue
        end
        f = strip.(split(s, ',', limit = 5))
        length(f) >= 4 || continue
        T, w = tryparse(Float64, f[1]), tryparse(Float64, f[3])
        (T === nothing || w === nothing) && continue
        w /= 100
        X = lowercase(f[4]) == "wt" ? w : x2w(w)
        push!(out, (T = T, type = lowercase(f[2]), X = X, source = length(f) >= 5 ? f[5] : ""))
    end
    return out
end

# =============================================================================
# APPLICATION
# =============================================================================
function build_app()
    # ---------------------------------------------------------------- state
    Ws, Wl = Ref(WS_START), Ref(WL_START)
    P      = Ref(P_START)
    phases = Ref(system_phases(Ws[], Wl[], NPC))
    lp     = Ref(GibbsLP(phases[], [0.5, 0.5]))
    npc    = Ref(NPC)

    points    = NamedTuple{(:X, :T, :c),Tuple{Float64,Float64,Int}}[]
    pending   = Ref{Union{Nothing,Tuple{Float64,Float64}}}(nothing)
    color_of  = Dict{Any,Int}()
    fieldplts = Any[]
    slotplts  = Any[]

    # X is always the WEIGHT fraction of Fa (the coordinate of the diagram)
    equil(X, T) = equilibrium!(lp[], X, T, P[])
    key(X, T)   = assemblage_key(phases[], equil(X, T).stable)

    # ---------------------------------------------------------------- layout
    fig = Figure(size = (1500, 950), fontsize = 15)
    Label(fig[0, 1:2], "Example 4 – melting of a binary solid solution (olivine Fo–Fa, melt L): " *
          "Gibbs energy minimisation by linear programming", fontsize = 20, font = :bold, word_wrap = true, tellwidth = false)

    left  = fig[1, 1] = GridLayout()
    right = fig[1, 2] = GridLayout()
    colsize!(fig.layout, 2, Relative(0.5))

    axTX = Axis(right[1, 1]; xlabel = XLABEL, ylabel = "T [°C]", xticks = XTICKS,
                title = "T–X diagram  (left-click to pick a point)")
    limits!(axTX, 0, 1, TLIM...)
    for name in (:rectanglezoom, :scrollzoom, :dragpan)
        name in keys(interactions(axTX)) && deregister_interaction!(axTX, name)
    end

    # one sub-grid per row, so that the widths of the rows do not add up (keeps the controls inside the window)
    ctrl = right[2, 1] = GridLayout(tellwidth = false)
    sg = SliderGrid(ctrl[1, 1],
        (label = "P [kbar]",         range = 0.0:0.5:30.0,   startvalue = P_START,  format = "{:.1f}"),
        (label = "W of Ol [kJ/mol]", range = -20.0:0.5:80.0, startvalue = WS_START, format = "{:.1f}"),
        (label = "W of L [kJ/mol]",  range = -20.0:0.5:80.0, startvalue = WL_START, format = "{:.1f}"),
        (label = "pseudocomp.",      range = 3:2:401,        startvalue = NPC,      format = n -> "$n per phase"),
        (label = "Grid nodes",       range = 5:1:100,        startvalue = NGRID,    format = n -> "$n × $n"))
    sl_P, sl_Ws, sl_Wl, sl_npc, sl_grid = sg.sliders

    cb = ctrl[2, 1] = GridLayout(halign = :left)
    btn_init  = Button(cb[1, 1], label = "Initialise grid")
    btn_comp  = Button(cb[1, 2], label = "Compute g–x", buttoncolor = RGBf(0.80, 0.88, 1.0))
    btn_clear = Button(cb[1, 3], label = "Clear points")
    ct = ctrl[3, 1] = GridLayout(halign = :left)
    tg_fields = Toggle(ct[1, 1], active = true);  Label(ct[1, 2], "show fields",       halign = :left)
    tg_line   = Toggle(ct[1, 3], active = true);  Label(ct[1, 4], "exact loop (W = 0)", halign = :left)
    tg_auto   = Toggle(ct[1, 5], active = false); Label(ct[1, 6], "compute on click",  halign = :left)
    ce = ctrl[4, 1] = GridLayout(halign = :left, tellwidth = false)
    cb_exp    = Checkbox(ce[1, 1], checked = false)
    Label(ce[1, 2], "experiments at 1 bar (data/olivine_1bar_experiments.csv; compare with P = 0)",
          halign = :left, justification = :left, word_wrap = true, tellwidth = false)
    colsize!(ce, 2, Auto(1))

    status = Label(right[3, 1], "", halign = :left, tellwidth = false, color = :gray25)
    cursor = Label(right[4, 1], " ", halign = :left, tellwidth = false, color = :gray25)
    # rowsize!(right, 1, Relative(0.63))   (removed: the diagram now takes the height the controls leave)

    axGX = [Axis(left[k, 1]; ylabel = "g [J/g]", xticks = XTICKS, xlabel = k == NSLOTS ? XLABEL : "")
            for k in 1:NSLOTS]
    linkxaxes!(axGX...)

    pend_pos = Observable(Point2f(NaN, NaN)); pend_vis = Observable(false)
    scatter!(axTX, pend_pos; marker = :xcross, markersize = 20, color = :black, visible = pend_vis)

    # experimental points (re-read from the file every time the box is ticked)
    expplts = Any[]
    function draw_experiments!(show)
        foreach(p -> p isa Legend ? delete!(p) : delete!(axTX, p), expplts); empty!(expplts)
        show || return
        data = read_experiments(EXP_FILE)
        styles = Dict("liquidus" => (:circle, LIQ_COLOR, "liquidus (exp.)"),
                      "solidus"  => (:rect, SOL_COLOR, "solidus (exp.)"),
                      "melting"  => (:diamond, :black, "melting point (exp.)"))
        for (type, (marker, col, lab)) in styles
            d = filter(r -> r.type == type, data)
            isempty(d) && continue
            push!(expplts, scatter!(axTX, [r.X for r in d], [r.T for r in d]; marker, markersize = 13,
                                    color = :white, strokecolor = col, strokewidth = 2.5, label = lab))
        end
        foreach(p -> translate!(p, 0, 0, 5), expplts)
        isempty(expplts) || push!(expplts, axislegend(axTX; position = :lt, framevisible = true,
                                                      labelsize = 13, padding = (6, 6, 4, 4)))
        n = length(data)
        status.text[] = @sprintf("%d experimental point%s read from %s", n, n == 1 ? "" : "s",
                                 basename(EXP_FILE)) * (isfile(EXP_FILE) ? "" : "  (file not found)")
    end

    # ------------------------------------------------------------ functions
    "Initialise: scan the grid, refine the boundaries, draw the fields and the exact loop."
    function init_grid!()
        n   = round(Int, sl_grid.value[])
        Xs  = [(i - 0.5) / n for i in 1:n]            # cell centres (wt fraction)
        Ts  = range(TLIM..., length = n)
        t0  = time()
        labels, asm = scan_grid(key, Xs, Ts)
        Xb, Tb      = refine_boundaries(key, Xs, Ts, labels)
        dt  = time() - t0

        for a in asm
            get!(color_of, a, mod1(length(color_of) + 1, length(FIELD_COLORS)))
        end
        cidx = [color_of[asm[l]] for l in labels]

        foreach(p -> delete!(axTX, p), fieldplts); empty!(fieldplts)
        cmap = cgrad(FIELD_COLORS, length(FIELD_COLORS), categorical = true)
        hm = heatmap!(axTX, Xs, Ts, cidx; colormap = cmap, alpha = 0.55,
                      colorrange = (0.5, length(FIELD_COLORS) + 0.5), visible = tg_fields.active)
        translate!(hm, 0, 0, -10)
        nodes = scatter!(axTX, [Point2f(X, T) for X in Xs for T in Ts];
                         markersize = 2.5, color = (:black, 0.25), visible = tg_fields.active)
        translate!(nodes, 0, 0, -9)
        bnd = scatter!(axTX, Xb, Tb; markersize = 6, color = :black, visible = tg_fields.active)
        translate!(bnd, 0, 0, -8)
        push!(fieldplts, hm, nodes, bnd)
        # exact loop, only meaningful for the ideal case
        if Ws[] == 0 && Wl[] == 0
            xs, xl, te = exact_loop(P[], range(TLIM..., length = 600))
            for (x, lab) in ((x2w.(xs), "solidus"), (x2w.(xl), "liquidus"))     # mole -> weight fraction
                l = lines!(axTX, x, te; color = :black, linestyle = :dash, linewidth = 1.5, visible = tg_line.active)
                translate!(l, 0, 0, -7); push!(fieldplts, l)
            end
        end
        for (id, idx) in field_regions(labels)
            length(idx) >= 3 || continue
            pos = Point2f(sum(Xs[I[1]] for I in idx) / length(idx), sum(Ts[I[2]] for I in idx) / length(idx))
            push!(fieldplts, text!(axTX, pos; text = join(asm[id], " + "), align = (:center, :center),
                                   fontsize = 20, font = :bold, color = :gray20, visible = tg_fields.active))
        end
        nlp = n^2 + 15 * length(Xb)
        tm = [MELTING[k].Tm + MELTING[k].dTdP * P[] for k in 1:2]
        status.text[] = @sprintf("Initialised: %d × %d grid, %d boundary points, ≈ %d LPs in %.2f s\n",
                                 n, n, length(Xb), nlp, dt) *
                        @sprintf("P = %.1f kbar (T_m: Fo %.0f, Fa %.0f °C), W(Ol) = %.1f, W(L) = %.1f kJ/mol",
                                 P[], tm..., Ws[], Wl[]) *
                        (Ws[] == 0 && Wl[] == 0 ? "" : "  (exact loop only for W = 0)")
    end

    "Draw the g–x diagram of panel k for bulk composition X (wt fraction Fa) at temperature T [°C]."
    function draw_slot!(k, X, T, c, rank)
        ax  = axGX[k]
        empty!(ax)
        ph  = phases[]
        r   = equil(X, T)                              # solved in moles: r.g, r.μ, r.G in kJ/mol
        xmol = xcomp(ph)                               # mole fractions of the pseudocompounds
        n   = npc[]
        iS, iL = 1:n, n+1:2n
        pl  = phase_list(ph, r)
        # a single solution phase has a unique tangent (μ): the tangent to its continuous curve
        one = length(pl) == 1 && 0 < pl[1].x < 1
        μmol = one ? (pl[1].name == SOL_NAME ? sol_mu(SOL_ENDS, Ws[], pl[1].x, T, P[]) :
                                                sol_mu(LIQ_ENDS, Wl[], pl[1].x, T, P[])) : r.μ      # [kJ/mol]

        # ---- display in weight: wt fraction on the x axis, g and μ per gram [J/g]
        xm  = x2w.(xmol)
        g   = 1000 .* r.g ./ molar_mass.(xmol)
        μ   = 1000 .* μmol ./ MOLAR_MASS               # intercepts of the tangent in this plot
        G   = 1000 * r.G / molar_mass(w2x(X))

        vlines!(ax, [X]; color = :gray50, linestyle = :dot)
        lines!(ax, xm[iS], g[iS]; color = SOL_COLOR, linewidth = 2.5)
        lines!(ax, xm[iL], g[iL]; color = LIQ_COLOR, linewidth = 2.5)
        if n <= 41
            scatter!(ax, xm[iS], g[iS]; markersize = 6, color = SOL_COLOR)
            scatter!(ax, xm[iL], g[iL]; markersize = 6, color = LIQ_COLOR)
        end
        hx, hg = lower_hull(xm, g)
        lines!(ax, hx, hg; color = c, linewidth = 3, linestyle = (:dot, :dense))
        lines!(ax, [0.0, 1.0], μ; color = c, linestyle = :dash, linewidth = 1.5)
        for p in pl                                     # stable phases on the tangent
            scatter!(ax, [p.w], [(1 - p.w) * μ[1] + p.w * μ[2]]; markersize = 14,
                     color = p.name == SOL_NAME ? SOL_COLOR : LIQ_COLOR, strokewidth = 2, strokecolor = c)
        end
        for (idx, name, col, f) in ((iS, SOL_NAME, SOL_COLOR, 0.12), (iL, LIQ_NAME, LIQ_COLOR, 0.88))
            m = idx[clamp(round(Int, f * n), 1, n)]
            text!(ax, [xm[m]], [g[m]]; text = name, offset = (0, 8), align = (:center, :bottom),
                  fontsize = 18, font = :bold, color = col)
        end
        scatter!(ax, [X], [G]; marker = :star5, markersize = 22, color = :gold, strokewidth = 1.5, strokecolor = :black)
        text!(ax, [0.0, 1.0], μ;
              text = [rich("μ", subscript(COMPONENTS[i]), @sprintf(" = %.1f J/g", round(μ[i], digits = 1) + 0.0)) for i in 1:2],
              align = [(:left, :top), (:right, :top)], offset = [(4, -6), (-4, -6)], fontsize = 13, color = c)

        lo, hi = extrema(vcat(g, μ)); pad = max(hi - lo, 1e-3)
        limits!(ax, -0.04, 1.04, lo - 0.2pad, hi + 0.25pad)

        amounts = join([@sprintf("%s %.0f wt%%", phase_label(p), 100p.amount) for p in pl], " + ")
        note = one ? "   |   μ: tangent to the " * pl[1].name * " curve" : ""
        ax.title    = rich(rich("T", subscript(string(rank)), color = c),
                           @sprintf(" = %.0f °C,  bulk %.0f wt%% %s,  P = %.1f kbar", T, 100X, COMPONENTS[2], P[]))
        ax.subtitle = "stable: " * amounts * @sprintf("   |   G = %.1f J/g", G) * note *       # 2nd line: molar values
                      @sprintf("\nμ(%s) = %.2f kJ/mol,   μ(%s) = %.2f kJ/mol", COMPONENTS[1], round(μmol[1], digits = 2) + 0.0,
                               COMPONENTS[2], round(μmol[2], digits = 2) + 0.0)
        ax.titlecolor = c

        xs = [p.w for p in pl]                          # tie-line in the T–X diagram (wt fractions)
        push!(slotplts,
              hlines!(axTX, [T]; color = c, linewidth = 1, linestyle = :dash),
              lines!(axTX, [minimum(xs), maximum(xs)], [T, T]; color = c, linewidth = 3),
              scatter!(axTX, xs, fill(T, length(xs)); markersize = 11, color = c),
              scatter!(axTX, [X], [T]; marker = :xcross, markersize = 18, color = c),
              text!(axTX, Point2f(1.0, T); text = rich("T", subscript(string(rank))), color = c,
                    align = (:right, :bottom), offset = (-4, 2), fontsize = 18, font = :bold))
    end

    function redraw_slots!()
        foreach(p -> delete!(axTX, p), slotplts); empty!(slotplts)
        order = sortperm(points; by = p -> p.T)
        np    = length(points)
        for k in 1:NSLOTS
            rank = NSLOTS - k + 1
            if rank <= np
                p = points[order[rank]]
                draw_slot!(k, p.X, p.T, SLOT_COLORS[p.c], rank)
            else
                empty!(axGX[k])
                axGX[k].title = "– empty –"; axGX[k].titlecolor = :gray40
                axGX[k].subtitle = "pick a point in the T–X diagram, then press Compute g–x\n "
                limits!(axGX[k], -0.04, 1.04, 0, 1)
            end
        end
    end

    function compute!()
        if pending[] === nothing
            status.text[] = "First left-click on the T–X diagram to choose a point."
            return
        end
        length(points) >= NSLOTS && popfirst!(points)
        used = [p.c for p in points]
        X, T = pending[]
        push!(points, (X = X, T = T, c = findfirst(c -> !(c in used), 1:NSLOTS)))
        pending[] = nothing; pend_vis[] = false
        redraw_slots!()
    end

    add_point!(X, T) = (pending[] = (X, T); compute!())      # X = weight fraction Fa

    function rebuild!()
        Ws[] = sl_Ws.value[]; Wl[] = sl_Wl.value[]; npc[] = round(Int, sl_npc.value[])
        phases[] = system_phases(Ws[], Wl[], npc[])
        lp[]     = GibbsLP(phases[], [0.5, 0.5])
        init_grid!()
        redraw_slots!()
    end

    # ------------------------------------------------------------ callbacks
    register_interaction!(axTX, :pick_point) do event::MouseEvent, axis
        if event.type === MouseEventTypes.leftclick
            X = clamp(Float64(event.data[1]), 0.0, 1.0)
            T = clamp(Float64(event.data[2]), TLIM...)
            pending[] = (X, T)
            pend_pos[] = Point2f(X, T); pend_vis[] = true
            pl = phase_list(phases[], equil(X, T))
            status.text[] = @sprintf("Picked %.0f wt%% %s, T = %.0f °C  (%s)  – press Compute g–x",
                                     100X, COMPONENTS[2], T, join(phase_label.(pl), " + "))
            tg_auto.active[] && compute!()
        end
        return Consume(false)
    end
    on(events(fig).mouseposition) do _
        if is_mouseinside(axTX.scene)
            X, T = mouseposition(axTX)
            if 0 <= X <= 1 && TLIM[1] <= T <= TLIM[2]
                pl = phase_list(phases[], equil(X, T))
                cursor.text[] = @sprintf("cursor: %.0f wt%% %s, T = %.0f °C  →  %s", 100X, COMPONENTS[2], T,
                                         join(phase_label.(pl), " + "))
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
    on(sl_P.value) do v
        P[] = v; init_grid!(); redraw_slots!()
    end
    for sl in (sl_Ws, sl_Wl, sl_npc)
        on(_ -> rebuild!(), sl.value)
    end
    on(draw_experiments!, cb_exp.checked)

    # ------------------------------------------------------- initialisation
    init_grid!()
    redraw_slots!()
    return fig, (add_point = add_point!, init_grid = init_grid!,
                 set_P  = v -> set_close_to!(sl_P, v),
                 set_Ws = v -> set_close_to!(sl_Ws, v),
                 set_Wl = v -> set_close_to!(sl_Wl, v), experiments = cb_exp)
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
