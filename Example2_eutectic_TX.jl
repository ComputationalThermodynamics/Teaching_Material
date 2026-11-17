#=
Example 2 – Interactive T–X diagram of a binary eutectic: diopside – anorthite
==============================================================================

Two components, diopside CaMgSi2O6 (Di) and anorthite CaAl2Si2O8 (An). Phases:

    Di(s)  pure solid, x = 0
    An(s)  pure solid, x = 1
    L      liquid (melt) solution of Di(l) and An(l), any x

MOLES IN THE CALCULATION, WEIGHT IN THE PLOTS. The thermodynamics is written in
moles: x is the MOLE fraction of An, g is in kJ per mole, and the linear
programme balances moles. The Di–An diagram of the course (part 4) is drawn in
WEIGHT per cent, so everything that is shown is converted for the display:

    composition   wt fraction  w = x M_An / [(1-x) M_Di + x M_An]
    phase amounts mass fractions (lever rule in the wt% diagram)
    g, μ          per gram:  g [J/g] = g [kJ/mol] / M̄,  μ_i [J/g] = μ_i [kJ/mol] / M_i

A straight line in (x, g per mole) is also a straight line in (w, g per gram),
so the common tangent stays a straight line in the wt% plot; its intercepts
are the chemical potentials per gram (the molar values are printed as well).

The liquid is a solution: its Gibbs energy is a curve in g–x space,

    g_L(x) = (1-x) g_Di(l) + x g_An(l) + RT [x ln x + (1-x) ln(1-x)] + W x (1-x)

To keep the *same* linear programme as in Example 0 (see gibbs_lp.jl), the curve
is discretised into many "pseudocompounds": liquids of fixed composition
x = 0, 1/(n-1), …, 1. Each of them is an ordinary phase with g = g0 - s T + v P
(the configurational entropy is part of s). The LP then picks the stable phases;
two neighbouring pseudocompounds together stand for one liquid whose
composition lies between them.

The pure liquids follow from the melting data of the solids:
    g_Di(l) - g_Di(s) = ΔH_Di (1 - T/T_m,Di) + ΔV_Di P          (and the same for An)

CALIBRATION. The melting points are those of the diagram (Di 1391 °C, An
1553 °C). The two liquidus curves of the diagram were digitised (117 points, in
wt%, data/DiAn_liquidus_digitized.txt), converted to mole fractions, and ΔH of
both end-members and W of the liquid were fitted to them
(calibration/fit_Example2_DiAn.jl):

    ΔH_Di = 128.1 kJ/mol,  ΔH_An = 154.7 kJ/mol,  W = -25.5 kJ/mol
    model eutectic 41.9 wt% An (35.9 mol%), 1274 °C   (diagram: 42 wt% An, 1274 °C)
    liquidus misfit: rms 2.2 K

ΔH and W are fitted to the liquidus, they are not calorimetric values. The
pressure dependence (melting slopes dT_m/dP) is illustrative only.

How to use
----------
* At start-up the program is *initialised*: the LP is solved on a T–X grid and
  the field boundaries are located by bisection between neighbouring grid
  nodes that have different assemblages (black dots).
* Left-click anywhere in the T–X diagram to pick a bulk composition and a
  temperature (a cross appears), then press "Compute g–x". The g–x diagram at
  that temperature appears on the left; the panels are sorted by temperature
  (T₁ lowest, at the bottom). Up to 3 panels; the oldest point is dropped after
  that. The isotherm and the tie-line of the coexisting phases are drawn in the
  T–X diagram.
* Turn on "compute on click" to skip the button.
* Sliders: pressure P, Margules parameter W of the liquid (calibrated value
  -25.5 kJ/mol; W = 0: ideal; W > 0 large enough: liquid immiscibility L + L), the
  number of liquid pseudocompounds (discretisation of the curve), and the grid size.
* "Show fields" / "exact liquidus" can be switched off so that students first
  have to discover the topology by themselves. The exact liquidus is computed
  for the continuous liquid (no pseudocompounds) as a check of the LP.

Run from a terminal:      julia Example2_eutectic_TX.jl
or from the Julia REPL:   include("Example2_eutectic_TX.jl")
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
const COMPONENTS = ["Di", "An"]                      # the two components
const MOLAR_MASS = [216.55, 278.21]                  # CaMgSi2O6, CaAl2Si2O8 [g/mol]: moles <-> weight
#                    name   amounts (mol)  g0    s     v
const SOLIDS = [Phase("Di", [1.0, 0.0],   0.0,  0.0,  0.0),
                Phase("An", [0.0, 1.0],   0.0,  0.0,  0.0)]
# g of the solids is the reference (= 0 at every P,T), so the g–x diagrams are easy
# to read. Any other values work as well (units as in Example 0).

# melting data of the pure end-members (in the order of SOLIDS)
#   ΔH [kJ/mol]   T_m [°C] at P = 0   dT_m/dP [K/kbar]  ->  ΔV = ΔS dT_m/dP,  ΔS = ΔH/T_m
const MELTING = [(ΔH = 128.1, Tm = 1391.0, dTdP = 9.0),            # diopside
                 (ΔH = 154.7, Tm = 1553.0, dTdP = 5.0)]            # anorthite
# ΔH calibrated against the Di–An diagram in mole fractions (see header); dT_m/dP is illustrative

const LIQ_NAME  = "L"
const W_START   = -25.5                              # Margules parameter of L [kJ/mol] (calibrated)
const WS_L      = 0.0                                # its T dependence: W(T) = W - T WS_L
const NPC       = 201                                # number of liquid pseudocompounds
const P_START   = 0.0                                # pressure [kbar]
const TLIM      = (1200.0, 1600.0)                   # temperature range [°C]
const NGRID     = 25                                 # initial grid nodes per axis
const NSLOTS    = 3                                  # number of g–x panels

const SLOT_COLORS  = [:royalblue3, :darkorange2, :forestgreen]
const LIQ_COLOR    = :firebrick
const FIELD_COLORS = ["#a6cee3", "#fdbf6f", "#b2df8a", "#fb9a99",
                      "#cab2d6", "#ffff99", "#8dd3c7", "#d9d9d9"]

const XTICKS = (0:0.2:1, [COMPONENTS[1], "20", "40", "60", "80", COMPONENTS[2]])
const XLABEL = "wt% " * COMPONENTS[2]

# pure liquids Di(l), An(l) from the melting data (T_m converted to K; ΔS = ΔH/T_m, ΔV = ΔS dT_m/dP)
const LIQ_ENDS = [let Tm = m.Tm + 273.15
                      Phase(s.name * "(l)", s.comp, s.g0 + m.ΔH, s.s + m.ΔH / Tm, s.v + m.ΔH / Tm * m.dTdP)
                  end for (s, m) in zip(SOLIDS, MELTING)]

# =============================================================================
# THERMODYNAMICS (independent of the plotting)
# =============================================================================
K(T_C) = T_C + 273.15                                # temperatures are entered in °C

# ----- moles <-> weight. The calculation uses the mole fraction x of the 2nd component,
#       the plots use its weight fraction w.
"weight fraction -> mole fraction of the 2nd component"
w2x(w) = (w / MOLAR_MASS[2]) / ((1 - w) / MOLAR_MASS[1] + w / MOLAR_MASS[2])
"mole fraction -> weight fraction of the 2nd component"
x2w(x) = x * MOLAR_MASS[2] / ((1 - x) * MOLAR_MASS[1] + x * MOLAR_MASS[2])
"mass [g] of one mole of a phase (or bulk) of mole fraction x"
molar_mass(x) = (1 - x) * MOLAR_MASS[1] + x * MOLAR_MASS[2]

"All phases: the two solids followed by the n liquid pseudocompounds (compositions in mole fractions)."
system_phases(W, n) = vcat(SOLIDS, solution_pseudocompounds(LIQ_NAME, LIQ_ENDS..., W, n; WS = WS_L))

"""
Solve at bulk composition w (WEIGHT fraction of the 2nd component, as on the diagram),
T [°C] and P [kbar]. The bulk is converted to moles; the result is in moles and kJ/mol.
"""
function equilibrium!(lp::GibbsLP, w, T, P)
    x = w2x(w)
    set_bulk!(lp, [1 - x, x])
    return solve!(lp, P, K(T))
end

"""
The physical phases of a solution r: name, amount (MASS fraction of the system),
weight fraction w and mole fraction x of the 2nd component.
Neighbouring liquid pseudocompounds are merged into one liquid.
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

phase_label(p) = p.name == LIQ_NAME ? @sprintf("%s(%.0f wt%% %s)", p.name, 100p.w, COMPONENTS[2]) : p.name

"""
Chemical potentials [kJ/mol] of the two components in the (continuous) liquid of
MOLE fraction X: the tangent to the g_L curve at X (T in °C).
    μ_1 = g_1(l) + RT ln(1-X) + W X²,   μ_2 = g_2(l) + RT ln X + W (1-X)²
"""
function liquid_mu(X, T, P, W)
    TK = K(T); RT = RGAS * TK; Wt = W - TK * WS_L
    gl = gibbs(LIQ_ENDS, P, TK)
    return [gl[1] + RT * log(1 - X) + Wt * X^2, gl[2] + RT * log(X) + Wt * (1 - X)^2]
end

"""
    exact_liquidus(W, P, Ts)

Liquidus of the continuous liquid (no pseudocompounds): for every T [°C] find the
liquid compositions X (mole fraction) with  μ_1(L) = g_1(s)  or  μ_2(L) = g_2(s).
A root is kept only if it is a true equilibrium: the other solid is not
supersaturated and no other liquid lies below the tangent (no unmixing).
Returns the points as vectors (X, T); X is a MOLE fraction (use x2w to plot it).
"""
function exact_liquidus(W, P, Ts; ny = 1001)
    ys = range(1e-9, 1 - 1e-9, length = ny)
    Xo, To = Float64[], Float64[]
    for T in Ts
        TK = K(T); RT = RGAS * TK; Wt = W - TK * WS_L
        gs = gibbs(SOLIDS, P, TK)
        gl = gibbs(LIQ_ENDS, P, TK)
        μA(x) = liquid_mu(x, T, P, W)[1]
        μB(x) = liquid_mu(x, T, P, W)[2]
        gL(y) = (1 - y) * gl[1] + y * gl[2] + RT * (xlogx(y) + xlogx(1 - y)) + Wt * y * (1 - y)
        for h in (x -> μA(x) - gs[1], x -> μB(x) - gs[2])
            hv = h.(ys)
            for i in 1:ny-1
                sign(hv[i]) == sign(hv[i+1]) && continue
                a, b = ys[i], ys[i+1]
                for _ in 1:40
                    m = (a + b) / 2
                    sign(h(m)) == sign(hv[i]) ? (a = m) : (b = m)
                end
                x  = (a + b) / 2
                mA, mB = μA(x), μB(x)
                ok = gs[1] - mA >= -1e-6 && gs[2] - mB >= -1e-6 &&
                     all(gL(y) - ((1 - y) * mA + y * mB) >= -1e-6 for y in ys)
                ok && (push!(Xo, x); push!(To, T))
            end
        end
    end
    return Xo, To
end


# =============================================================================
# APPLICATION
# =============================================================================
function build_app()
    # ---------------------------------------------------------------- state
    W      = Ref(W_START)
    P      = Ref(P_START)
    phases = Ref(system_phases(W[], NPC))
    lp     = Ref(GibbsLP(phases[], [0.5, 0.5]))

    points    = NamedTuple{(:X, :T, :c),Tuple{Float64,Float64,Int}}[]   # c = colour index
    pending   = Ref{Union{Nothing,Tuple{Float64,Float64}}}(nothing)
    color_of  = Dict{Any,Int}()                       # assemblage -> field colour (kept stable)
    fieldplts = Any[]                                 # plots of the current grid (re-created)
    slotplts  = Any[]                                 # isotherms / tie-lines   (re-created)

    # X is always the WEIGHT fraction of An (the coordinate of the diagram)
    equil(X, T) = equilibrium!(lp[], X, T, P[])
    key(X, T)   = assemblage_key(phases[], equil(X, T).stable)

    # ---------------------------------------------------------------- layout
    fig = Figure(size = (1500, 950), fontsize = 15)
    Label(fig[0, 1:2], "Example 2 – binary eutectic (diopside Di, anorthite An, liquid L): " *
          "Gibbs energy minimisation by linear programming", fontsize = 20, font = :bold, word_wrap = true, tellwidth = false)

    left  = fig[1, 1] = GridLayout()
    right = fig[1, 2] = GridLayout()
    colsize!(fig.layout, 2, Relative(0.5))

    axTX = Axis(right[1, 1]; xlabel = XLABEL, ylabel = "T [°C]",
                xticks = XTICKS,
                title = "T–X diagram  (left-click to pick a point)")
    limits!(axTX, 0, 1, TLIM...)
    for name in (:rectanglezoom, :scrollzoom, :dragpan)   # keep the diagram fixed
        name in keys(interactions(axTX)) && deregister_interaction!(axTX, name)
    end

    # controls
    # one sub-grid per row, so that the widths of the rows do not add up (keeps the controls inside the window)
    ctrl = right[2, 1] = GridLayout(tellwidth = false)
    sg = SliderGrid(ctrl[1, 1],
        (label = "P [kbar]",        range = 0.0:0.5:30.0,  startvalue = P_START, format = "{:.1f}"),
        (label = "W of L [kJ/mol]", range = -40.0:0.5:40.0, startvalue = W_START, format = "{:.1f}"),
        (label = "L pseudocomp.",   range = 3:2:401,       startvalue = NPC,     format = n -> "$n"),
        (label = "Grid nodes",      range = 5:1:100,       startvalue = NGRID,   format = n -> "$n × $n"))
    sl_P, sl_W, sl_npc, sl_grid = sg.sliders

    cb = ctrl[2, 1] = GridLayout(halign = :left)
    btn_init  = Button(cb[1, 1], label = "Initialise grid")
    btn_comp  = Button(cb[1, 2], label = "Compute g–x", buttoncolor = RGBf(0.80, 0.88, 1.0))
    btn_clear = Button(cb[1, 3], label = "Clear points")
    ct = ctrl[3, 1] = GridLayout(halign = :left)
    tg_fields = Toggle(ct[1, 1], active = true);  Label(ct[1, 2], "show fields",      halign = :left)
    tg_line   = Toggle(ct[1, 3], active = true);  Label(ct[1, 4], "exact liquidus",   halign = :left)
    tg_auto   = Toggle(ct[1, 5], active = false); Label(ct[1, 6], "compute on click", halign = :left)

    status = Label(right[3, 1], "", halign = :left, tellwidth = false, color = :gray25)
    cursor = Label(right[4, 1], " ", halign = :left, tellwidth = false, color = :gray25)
    # rowsize!(right, 1, Relative(0.72))   (removed: the diagram now takes the height the controls leave)

    # g–x panels (top = highest T, as in figure 2)
    axGX = [Axis(left[k, 1]; ylabel = "g [J/g]",
                 xticks = XTICKS,
                 xlabel = k == NSLOTS ? XLABEL : "")
            for k in 1:NSLOTS]
    linkxaxes!(axGX...)

    # pending point (cross)
    pend_pos = Observable(Point2f(NaN, NaN)); pend_vis = Observable(false)
    scatter!(axTX, pend_pos; marker = :xcross, markersize = 20, color = :black, visible = pend_vis)

    # ------------------------------------------------------------ functions
    "Initialise: scan the grid, refine the boundaries, draw the fields and the exact liquidus."
    function init_grid!()
        n   = round(Int, sl_grid.value[])
        Xs  = [(i - 0.5) / n for i in 1:n]            # cell centres (wt fraction): X = 0 and 1 are not nodes
        Ts  = range(TLIM..., length = n)
        t0  = time()
        labels, asm = scan_grid(key, Xs, Ts)
        Xb, Tb      = refine_boundaries(key, Xs, Ts, labels)
        dt  = time() - t0
        Xl, Tl      = exact_liquidus(W[], P[], range(TLIM..., length = 600))
        Xl          = x2w.(Xl)                        # mole fraction -> weight fraction for the plot

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
        liq = scatter!(axTX, Xl, Tl; markersize = 2.5, color = :black, visible = tg_line.active)
        translate!(liq, 0, 0, -7)
        push!(fieldplts, hm, nodes, bnd, liq)

        # label every connected field region at the centre of its grid nodes
        for (id, idx) in field_regions(labels)
            length(idx) >= 3 || continue            # skip slivers of one or two nodes
            pos = Point2f(sum(Xs[I[1]] for I in idx) / length(idx),
                          sum(Ts[I[2]] for I in idx) / length(idx))
            push!(fieldplts, text!(axTX, pos; text = join(asm[id], " + "),
                  align = (:center, :center), fontsize = 20, font = :bold,
                  color = :gray20, visible = tg_fields.active))
        end
        nlp = n^2 + 15 * length(Xb)
        status.text[] = @sprintf("Initialised: %d × %d grid, %d boundary points, ≈ %d LPs in %.2f s\n",
                                 n, n, length(Xb), nlp, dt) *
                        @sprintf("P = %.1f kbar, W = %.1f kJ/mol, %d L pseudocompounds",
                                 P[], W[], length(phases[]) - length(SOLIDS))
    end

    "Draw the g–x diagram of panel k for bulk composition X (wt fraction) at temperature T."
    function draw_slot!(k, X, T, c, rank)
        ax  = axGX[k]
        empty!(ax)
        ph  = phases[]
        r   = equil(X, T)                              # solved in moles: r.g, r.μ, r.G in kJ/mol
        xmol = xcomp(ph)                               # mole fractions of the phases
        iS  = 1:length(SOLIDS)
        iL  = length(SOLIDS)+1:length(ph)
        st  = falses(length(ph)); st[r.stable] .= true
        pl  = phase_list(ph, r)
        # A single liquid has a unique tangent (μ). When the bulk composition hits a
        # pseudocompound exactly, the LP sees only that one point and returns an
        # arbitrary supporting line; use the tangent to the continuous L curve instead.
        liq_only = length(pl) == 1 && pl[1].name == LIQ_NAME && 0 < pl[1].x < 1
        μmol = liq_only ? liquid_mu(pl[1].x, T, P[], W[]) : r.μ          # [kJ/mol]

        # ---- display in weight: wt fraction on the x axis, g and μ per gram [J/g]
        xm  = x2w.(xmol)
        g   = 1000 .* r.g ./ molar_mass.(xmol)
        μ   = 1000 .* μmol ./ MOLAR_MASS               # intercepts of the tangent in this plot
        G   = 1000 * r.G / molar_mass(w2x(X))

        vlines!(ax, [X]; color = :gray50, linestyle = :dot)                              # bulk
        lines!(ax, xm[iL], g[iL]; color = LIQ_COLOR, linewidth = 2)                        # g of L
        length(iL) <= 41 && scatter!(ax, xm[iL], g[iL]; markersize = 6, color = LIQ_COLOR) # pseudocompounds
        hx, hg = lower_hull(xm, g)
        lines!(ax, hx, hg; color = c, linewidth = 3, linestyle = (:dot, :dense))           # G of system (all X)
        lines!(ax, [0.0, 1.0], μ; color = c, linestyle = :dash, linewidth = 1.5)          # tangent (μ)
        scatter!(ax, xm[iS], g[iS]; markersize = 16, strokewidth = 2, strokecolor = c,
                 color = [st[i] ? c : :white for i in iS])
        text!(ax, xm[iS], g[iS]; text = [p.name for p in SOLIDS], offset = (0, 12),
              align = (:center, :bottom), fontsize = 15)
        sL = filter(in(iL), r.stable)                                                    # stable liquid(s)
        isempty(sL) || scatter!(ax, xm[sL], g[sL]; markersize = 13, color = LIQ_COLOR,
                                strokewidth = 2, strokecolor = c)
        im = iL[cld(length(iL), 2)]                                                       # middle of the curve
        text!(ax, [xm[im]], [g[im]]; text = LIQ_NAME, offset = (0, 8), align = (:center, :bottom),
              fontsize = 16, font = :bold, color = LIQ_COLOR)
        scatter!(ax, [X], [G]; marker = :star5, markersize = 22, color = :gold, strokewidth = 1.5, strokecolor = :black)
        text!(ax, [0.0, 1.0], μ;
              text = [rich("μ", subscript(COMPONENTS[i]), @sprintf(" = %.1f J/g", round(μ[i], digits = 1) + 0.0)) for i in 1:2],
              align = [(:left, :top), (:right, :top)], offset = [(4, -6), (-4, -6)], fontsize = 13, color = c)

        lo, hi = extrema(vcat(g, μ)); pad = max(hi - lo, 1e-3)
        limits!(ax, -0.06, 1.06, lo - 0.25pad, hi + 0.25pad)

        amounts = join([@sprintf("%s %.0f wt%%", phase_label(p), 100p.amount) for p in pl], " + ")
        note = ""
        if liq_only
            note = "   |   μ: tangent to L at its composition"
        elseif length(r.stable) < length(COMPONENTS)
            note = "   |   1 phase: μ not unique (any line through it below the others)"
        else
            # a solid that is not stable but touches the tangent: we sit on the eutectic
            on_line = [SOLIDS[i].name for i in iS
                       if !st[i] && abs(g[i] - (μ[1] * (1 - xm[i]) + μ[2] * xm[i])) < 1e-6 * (1 + abs(g[i]))]
            isempty(on_line) || (note = "   |   on eutectic: " * join(vcat(phase_label.(pl), on_line), ", ") * " coexist")
        end
        ax.title    = rich(rich("T", subscript(string(rank)), color = c),
                           @sprintf(" = %.0f °C,  bulk %.0f wt%% %s,  P = %.1f kbar", T, 100X, COMPONENTS[2], P[]))
        ax.subtitle = "stable: " * amounts * @sprintf("   |   G = %.1f J/g", G) * note *       # 2nd line: molar values
                      @sprintf("\nμ(%s) = %.2f kJ/mol,   μ(%s) = %.2f kJ/mol", COMPONENTS[1], round(μmol[1], digits = 2) + 0.0,
                               COMPONENTS[2], round(μmol[2], digits = 2) + 0.0)
        ax.titlecolor = c

        # isotherm, tie-line and bulk point in the T–X diagram (wt fractions)
        xs = [p.w for p in pl]
        push!(slotplts,
              hlines!(axTX, [T]; color = c, linewidth = 1, linestyle = :dash),
              lines!(axTX, [minimum(xs), maximum(xs)], [T, T]; color = c, linewidth = 3),
              scatter!(axTX, xs, fill(T, length(xs)); markersize = 11, color = c),
              scatter!(axTX, [X], [T]; marker = :xcross, markersize = 18, color = c),
              text!(axTX, Point2f(1.0, T); text = rich("T", subscript(string(rank))), color = c,
                    align = (:right, :bottom), offset = (-4, 2), fontsize = 18, font = :bold))
    end

    "Re-draw all g–x panels (sorted by T: T₁ lowest at the bottom) and the isotherms."
    function redraw_slots!()
        foreach(p -> delete!(axTX, p), slotplts); empty!(slotplts)
        order = sortperm(points; by = p -> p.T)
        np    = length(points)
        for k in 1:NSLOTS
            rank = NSLOTS - k + 1                  # rank 1 = lowest T, drawn in the bottom panel
            if rank <= np
                p = points[order[rank]]
                draw_slot!(k, p.X, p.T, SLOT_COLORS[p.c], rank)
            else
                empty!(axGX[k])
                axGX[k].title = "– empty –"; axGX[k].titlecolor = :gray40
                axGX[k].subtitle = "pick a point in the T–X diagram, then press Compute g–x\n "
                limits!(axGX[k], -0.06, 1.06, 0, 1)
            end
        end
    end

    "Move the pending point into a g–x panel."
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

    "Programmatic access (used for testing / screenshots)."
    function add_point!(X, T)
        pending[] = (X, T); compute!()
    end

    "New W or number of pseudocompounds: new phases, new LP."
    function rebuild!()
        W[]      = sl_W.value[]
        phases[] = system_phases(W[], round(Int, sl_npc.value[]))
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
            status.text[] = @sprintf("Picked %.0f wt%% %s, T = %.0f °C  (%s)  – press Compute g–x", 100X, COMPONENTS[2],
                                     T, join(phase_label.(pl), " + "))
            tg_auto.active[] && compute!()
        end
        return Consume(false)
    end

    on(events(fig).mouseposition) do _
        if is_mouseinside(axTX.scene)
            X, T = mouseposition(axTX)
            if 0 <= X <= 1 && TLIM[1] <= T <= TLIM[2]
                pl = phase_list(phases[], equil(X, T))
                cursor.text[] = @sprintf("cursor: %.0f wt%% %s, T = %.0f °C  →  %s", 100X, COMPONENTS[2],
                                         T, join(phase_label.(pl), " + "))
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
    on(sl_P.value) do v                 # new pressure: everything is recomputed
        P[] = v
        init_grid!()
        redraw_slots!()
    end
    on(_ -> rebuild!(), sl_W.value)
    on(_ -> rebuild!(), sl_npc.value)

    # ------------------------------------------------------- initialisation
    init_grid!()
    redraw_slots!()
    return fig, (add_point = add_point!, init_grid = init_grid!,
                 set_P   = v -> set_close_to!(sl_P, v),
                 set_W   = v -> set_close_to!(sl_W, v),
                 set_npc = v -> set_close_to!(sl_npc, v))
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
