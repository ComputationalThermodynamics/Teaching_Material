#=
Example 3 – Interactive T–X diagram of a binary peritectic: forsterite – silica
================================================================================

Same as Example 2, but with an intermediate solid that melts incongruently and
with a liquid that unmixes at the silica-rich end.
Two components, forsterite Mg2SiO4 (Fo) and silica SiO2. Phases:

    Fo    forsterite Mg2SiO4,                              x = 0
    En    enstatite MgSiO3 = ½ Mg2SiO4 + ½ SiO2,           x = 0.5
          (protoenstatite, "Pr" in the diagram)
    Crs   cristobalite SiO2,                               x = 1
    L     liquid solution of Fo(l) and SiO2(l), any x

MOLES IN THE CALCULATION, WEIGHT IN THE PLOTS. The thermodynamics is written in
moles: x is the MOLE fraction of SiO2 in Mg2SiO4–SiO2, g is in kJ per mole of
these components, and the linear programme balances moles. The diagram of the
course (part 4: Morse 1980, Fig. 10.1, after Bowen & Andersen) is drawn in
WEIGHT per cent SiO2 of the whole system (MgO–SiO2 scale), so the plots use

    wt fraction SiO2  w = M_SiO2 / [(1-x) M_Fo + x M_SiO2]
                      (forsterite 42.7 wt%, enstatite 59.85 wt%, silica 100 wt%)
    phase amounts     mass fractions (lever rule in the wt% diagram)
    g, μ              per gram:  g [J/g] = g [kJ/mol] / M̄,  μ_i [J/g] = μ_i [kJ/mol] / M_i

A straight line in (x, g per mole) is also a straight line in (w, g per gram),
so the common tangent stays a straight line in the wt% plot.

Enstatite is stable at low T. On heating it does not melt to a liquid of its
own composition: at the peritectic temperature T_p it breaks down to Fo plus a
liquid that is *richer in SiO2 than En* (En → Fo + L). In g–x space this is the
temperature at which En lies exactly on the tangent from Fo to the liquid
curve; above T_p En lies above that tangent and is metastable. Between the
peritectic and the eutectic (En + Crs + L) the En + L field appears.
Silica-rich liquids unmix into two liquids (L + L) above 1695 °C: the g curve
of the liquid has a hump there (see Example 1 for a solvus).

The pure liquids follow from the melting data of the solids:
    g(l) - g(s) = ΔH (1 - T/T_m) + ΔV P

Liquid model (per mole of components; Redlich–Kister excess terms):
    g_L = (1-x) g_Fo(l) + x g_SiO2(l) + RT [x ln x + (1-x) ln(1-x)]
        + x (1-x) Σ_k (a_k - b_k T) (1-2x)^k,      k = 0, 1, 2
One term (k = 0) is the regular solution of Example 2; three terms are needed
here because the liquid mixes with a negative excess energy next to forsterite
and unmixes next to silica.

The "exact" boundaries are those of the continuous liquid (no pseudocompounds):
for every solid k of composition x_k the liquidus from
    (1-x_k) μ_Fo(L) + x_k μ_SiO2(L) = g_k
and the two-liquid field from the common tangent to the liquid curve.

CALIBRATION. The curves of the diagram were digitised (229 points, in wt% SiO2,
data/FoSiO2_liquidus_Morse_digitized.txt), converted to mole fractions and
fitted (calibration/fit_Example3_FoSiO2.jl). Melting points as in the diagram
(Fo 1890 °C, Crs 1713 °C); ΔH_Fo = 142 kJ/mol as in Example 4. Fitted:

    ΔH_SiO2 = 18.19 kJ/mol;   En: g = -2.414 kJ/mol relative to ½ Fo + ½ SiO2
    a = (129.53, 28.01, 111.27) kJ/mol,  b = (0.07848, 0.02667, 0.03521) kJ/mol/K

                           model                          diagram
    peritectic             1555.5 °C, L 61.3 wt% SiO2     1557 °C, ≈61.3 wt%
    eutectic               1545.3 °C, L 64.8 wt% SiO2     1543 °C, ≈65.0 wt%
    two liquids + Crs      1697.8 °C, L 69.3 and 98.7     1695 °C, ≈69.4 and 98.7 wt%
    liquidus misfit: rms 3.9 K (134 points); two-liquid field within 0.7 wt% (rms)

These are fitted parameters (only ΔH_Fo is a calorimetric value). The pressure
dependence (melting slopes dT_m/dP) is illustrative.

How to use
----------
As in Example 2: left-click in the T–X diagram, then "Compute g–x". The slider
"ΔW of L" adds a regular-solution term ΔW x(1-x) to the calibrated liquid
(ΔW < 0: the two-liquid field shrinks; ΔW > 0: it grows).

Suggested exercises: pick points at the En composition (59.85 wt% SiO2) just
below and just above T_p and compare the g–x diagrams; cool liquids with 55,
60.5, 63 and 67 wt% SiO2 and follow the phases (which solid appears first,
what happens at the peritectic and at the eutectic?); pick a point with
85 wt% SiO2 at 1750 °C (two liquids) and look at the hump of the liquid curve.

Run from a terminal:      julia Example3_peritectic_TX.jl
or from the Julia REPL:   include("Example3_peritectic_TX.jl")
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
const COMPONENTS = ["Fo", "SiO2"]                    # Mg2SiO4 and SiO2: the components of the calculation
const MOLAR_MASS = [140.69, 60.08]                   # Mg2SiO4, SiO2 [g/mol]: moles <-> weight
#                    name    amounts (mol)   g0       s     v
const SOLIDS = [Phase("Fo",  [1.0, 0.0],     0.0,     0.0,  0.0),
                Phase("En",  [0.5, 0.5],    -2.414,   0.0,  0.0),      # MgSiO3 = ½ Mg2SiO4 + ½ SiO2
                Phase("Crs", [0.0, 1.0],     0.0,     0.0,  0.0)]
# g of Fo and cristobalite is the reference (= 0 at every P,T). En: g = -2.414 kJ/mol relative
# to the mechanical mixture ½ Fo + ½ SiO2 (calibrated, see header). Units as in Example 0.

# melting data of the pure end-members: index of the solid in SOLIDS,
#   ΔH [kJ/mol]   T_m [°C] at P = 0   dT_m/dP [K/kbar]  ->  ΔV = ΔS dT_m/dP,  ΔS = ΔH/T_m
const MELTING = [(solid = 1, ΔH = 142.0,  Tm = 1890.0, dTdP = 4.7),    # Fo(s)  -> Fo(l)
                 (solid = 3, ΔH =  18.19, Tm = 1713.0, dTdP = 5.5)]    # Crs    -> SiO2(l)
# melting points of the diagram; ΔH_Fo as in Example 4, ΔH_SiO2 calibrated; dT_m/dP is illustrative

# excess Gibbs energy of the liquid: x(1-x) Σ_k (a_k - b_k T)(1-2x)^k  (calibrated, see header)
const RK_H = [129.53, 28.01, 111.27]                 # a_k [kJ/mol]
const RK_S = [0.07848, 0.02667, 0.03521]             # b_k [kJ/mol/K]

const LIQ_NAME  = "L"
const W_START   = 0.0                                # ΔW added to a_0 by the slider [kJ/mol]
const NPC       = 201                                # number of liquid pseudocompounds
const P_START   = 0.0                                # pressure [kbar]
const TLIM      = (1500.0, 1950.0)                   # temperature range [°C]
const NGRID     = 25                                 # initial grid nodes per axis
const NSLOTS    = 3                                  # number of g–x panels

const SLOT_COLORS  = [:royalblue3, :darkorange2, :forestgreen]
const LIQ_COLOR    = :firebrick
const FIELD_COLORS = ["#a6cee3", "#fdbf6f", "#b2df8a", "#fb9a99",
                      "#cab2d6", "#ffff99", "#8dd3c7", "#d9d9d9"]

# pure liquids Fo(l), SiO2(l) from the melting data (T_m converted to K; ΔS = ΔH/T_m, ΔV = ΔS dT_m/dP)
const LIQ_ENDS = [let s = SOLIDS[m.solid], Tm = m.Tm + 273.15
                      Phase(COMPONENTS[i] * "(l)", s.comp, s.g0 + m.ΔH, s.s + m.ΔH / Tm, s.v + m.ΔH / Tm * m.dTdP)
                  end for (i, m) in enumerate(MELTING)]

# =============================================================================
# THERMODYNAMICS (independent of the plotting)
# =============================================================================
K(T_C) = T_C + 273.15                                # temperatures are entered in °C

# ----- moles <-> weight. The calculation uses the mole fraction x of SiO2 in Mg2SiO4–SiO2,
#       the plots use the weight fraction w of SiO2 in the whole system (as in the diagram).
#       Each mole of either component contains one mole of SiO2, so w = M_SiO2 / M̄.
"mass [g] of one mole of a phase (or bulk) of mole fraction x"
molar_mass(x) = (1 - x) * MOLAR_MASS[1] + x * MOLAR_MASS[2]
"mole fraction of SiO2 -> weight fraction of SiO2"
x2w(x) = MOLAR_MASS[2] / molar_mass(x)
"weight fraction of SiO2 -> mole fraction of SiO2"
w2x(w) = clamp((MOLAR_MASS[1] - MOLAR_MASS[2] / w) / (MOLAR_MASS[1] - MOLAR_MASS[2]), 0.0, 1.0)

const X_MIN = x2w(0.0)                               # left end of the axis: forsterite, 42.7 wt% SiO2
const X_EN  = x2w(0.5)                               # enstatite, 59.85 wt% SiO2

# x-axis ticks: every 10 wt%; the solids are labelled by name
const XTICKS = ([X_MIN, 0.5, X_EN, 0.7, 0.8, 0.9, 1.0], ["Fo", "50", "En", "70", "80", "90", "100"])
const XLABEL = "wt% " * COMPONENTS[2]

"Partial excess functions of the term x(1-x)(1-2x)^k: its contributions to μ_Fo and μ_SiO2."
function rk_partial(x, k)
    u  = 1 - 2x
    G  = x * (1 - x) * u^k
    dG = u^(k + 1) - (k == 0 ? 0.0 : 2k * x * (1 - x) * u^(k - 1))
    return G - x * dG, G + (1 - x) * dG
end

"""
Liquid pseudocompounds: n liquids of fixed composition, evenly spaced on the wt% axis of
the diagram. Each one is an ordinary `Phase` (composition in moles, g = g0 - s T + v P per
mole of components); dW [kJ/mol] is added to the first excess term a_0.
"""
function liquid_pseudocompounds(dW, n)
    e1, e2 = LIQ_ENDS
    map(range(X_MIN, 1.0, length = n)) do w
        x  = w2x(w)
        hx = x * (1 - x) * (dW + sum(RK_H[k] * (1 - 2x)^(k - 1) for k in eachindex(RK_H)))   # excess enthalpy
        sx = x * (1 - x) * sum(RK_S[k] * (1 - 2x)^(k - 1) for k in eachindex(RK_S))          # excess entropy
        Phase(LIQ_NAME, [1 - x, x],
              (1 - x) * e1.g0 + x * e2.g0 + hx,
              (1 - x) * e1.s  + x * e2.s  - RGAS * (xlogx(x) + xlogx(1 - x)) + sx,
              (1 - x) * e1.v  + x * e2.v)
    end
end

"All phases: the solids followed by the n liquid pseudocompounds (compositions in mole fractions)."
system_phases(dW, n) = vcat(SOLIDS, liquid_pseudocompounds(dW, n))

"""
Solve at bulk composition w (WEIGHT fraction SiO2, as on the diagram), T [°C] and P [kbar].
The bulk is converted to moles; the result is in moles and kJ/mol.
"""
function equilibrium!(lp::GibbsLP, w, T, P)
    x = w2x(w)
    set_bulk!(lp, [1 - x, x])
    return solve!(lp, P, K(T))
end

"""
The physical phases of a solution r: name, amount (MASS fraction of the system),
weight fraction w and mole fraction x of SiO2.
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

phase_label(p) = p.name == LIQ_NAME ? @sprintf("%s(%.1f wt%% %s)", p.name, 100p.w, COMPONENTS[2]) : p.name

"""
Chemical potentials [kJ/mol] of Fo and SiO2 in the (continuous) liquid of MOLE
fraction x: the tangent to the g_L curve at x (T in °C).
    μ_Fo = g_Fo(l) + RT ln(1-x) + excess,   μ_SiO2 = g_SiO2(l) + RT ln x + excess
"""
function liquid_mu(x, T, P, dW)
    TK = K(T); RT = RGAS * TK
    gl = gibbs(LIQ_ENDS, P, TK)
    μ1, μ2 = gl[1] + RT * log(1 - x), gl[2] + RT * log(x)
    for k in eachindex(RK_H)
        L = RK_H[k] - TK * RK_S[k] + (k == 1 ? dW : 0.0)
        e1, e2 = rk_partial(x, k - 1)
        μ1 += L * e1; μ2 += L * e2
    end
    return [μ1, μ2]
end

"""
    exact_boundaries(dW, P, Ts)

Phase boundaries of the continuous liquid (no pseudocompounds), for every T [°C]:
* liquidus of every solid k (composition x_k): the liquid compositions x with
      (1-x_k) μ_Fo(L) + x_k μ_SiO2(L) = g_k.
  A root is kept only if it is a true equilibrium: no other solid is
  supersaturated and no other liquid lies below the tangent (no unmixing);
* the two-liquid field: the two ends of a common tangent to the liquid curve
  (a gap in its lower convex hull) that no solid undercuts.
Returns the points as vectors (x, T); x is a MOLE fraction (use x2w to plot it).
"""
function exact_boundaries(dW, P, Ts; ny = 1001)
    ys = range(1e-9, 1 - 1e-9, length = ny)
    xs = xcomp(SOLIDS)
    Xo, To = Float64[], Float64[]
    for T in Ts
        gs  = gibbs(SOLIDS, P, K(T))
        μs  = [liquid_mu(y, T, P, dW) for y in ys]                   # tangent of the liquid at every y
        gL  = [(1 - y) * m[1] + y * m[2] for (y, m) in zip(ys, μs)]  # g of the liquid
        stable(m) = all(gs[j] - ((1 - xs[j]) * m[1] + xs[j] * m[2]) >= -1e-6 for j in eachindex(SOLIDS)) &&
                    all(gL[q] - ((1 - ys[q]) * m[1] + ys[q] * m[2]) >= -1e-6 for q in 1:ny)
        for k in eachindex(SOLIDS)                                   # liquidus of solid k
            h(x) = (m = liquid_mu(x, T, P, dW); (1 - xs[k]) * m[1] + xs[k] * m[2] - gs[k])
            hv = [(1 - xs[k]) * m[1] + xs[k] * m[2] - gs[k] for m in μs]
            for i in 1:ny-1
                sign(hv[i]) == sign(hv[i+1]) && continue
                a, b = ys[i], ys[i+1]
                for _ in 1:40
                    m = (a + b) / 2
                    sign(h(m)) == sign(hv[i]) ? (a = m) : (b = m)
                end
                x = (a + b) / 2
                stable(liquid_mu(x, T, P, dW)) && (push!(Xo, x); push!(To, T))
            end
        end
        hx, hg = lower_hull(collect(ys), gL)                         # two liquids: gap in the hull of L
        for i in 1:length(hx)-1
            hx[i+1] - hx[i] > 2.5 / (ny - 1) || continue
            slope = (hg[i+1] - hg[i]) / (hx[i+1] - hx[i])
            line(x) = hg[i] + slope * (x - hx[i])                    # the common tangent
            all(gs[j] - line(xs[j]) >= -1e-6 for j in eachindex(SOLIDS)) || continue
            push!(Xo, hx[i], hx[i+1]); push!(To, T, T)
        end
    end
    return Xo, To
end

# =============================================================================
# APPLICATION
# =============================================================================
function build_app()
    # ---------------------------------------------------------------- state
    W      = Ref(W_START)                             # ΔW of the liquid (slider)
    P      = Ref(P_START)
    phases = Ref(system_phases(W[], NPC))
    lp     = Ref(GibbsLP(phases[], [0.5, 0.5]))

    points    = NamedTuple{(:X, :T, :c),Tuple{Float64,Float64,Int}}[]   # c = colour index
    pending   = Ref{Union{Nothing,Tuple{Float64,Float64}}}(nothing)
    color_of  = Dict{Any,Int}()                       # assemblage -> field colour (kept stable)
    fieldplts = Any[]                                 # plots of the current grid (re-created)
    slotplts  = Any[]                                 # isotherms / tie-lines   (re-created)

    # X is always the WEIGHT fraction of SiO2 (the coordinate of the diagram, X_MIN … 1)
    equil(X, T) = equilibrium!(lp[], X, T, P[])
    key(X, T)   = assemblage_key(phases[], equil(X, T).stable)

    # ---------------------------------------------------------------- layout
    fig = Figure(size = (1500, 950), fontsize = 15)
    Label(fig[0, 1:2], "Example 3 – binary peritectic (forsterite Fo, enstatite En, cristobalite Crs, liquid L): " *
          "Gibbs energy minimisation by linear programming", fontsize = 20, font = :bold, word_wrap = true, tellwidth = false)

    left  = fig[1, 1] = GridLayout()
    right = fig[1, 2] = GridLayout()
    colsize!(fig.layout, 2, Relative(0.5))

    axTX = Axis(right[1, 1]; xlabel = XLABEL, ylabel = "T [°C]",
                xticks = XTICKS,
                title = "T–X diagram  (left-click to pick a point)")
    limits!(axTX, X_MIN, 1, TLIM...)
    for name in (:rectanglezoom, :scrollzoom, :dragpan)   # keep the diagram fixed
        name in keys(interactions(axTX)) && deregister_interaction!(axTX, name)
    end

    # controls
    # one sub-grid per row, so that the widths of the rows do not add up (keeps the controls inside the window)
    ctrl = right[2, 1] = GridLayout(tellwidth = false)
    sg = SliderGrid(ctrl[1, 1],
        (label = "P [kbar]",        range = 0.0:0.5:30.0,  startvalue = P_START, format = "{:.1f}"),
        (label = "ΔW of L [kJ/mol]", range = -40.0:0.5:40.0, startvalue = W_START, format = "{:.1f}"),
        (label = "L pseudocomp.",   range = 3:2:401,       startvalue = NPC,     format = n -> "$n"),
        (label = "Grid nodes",      range = 5:1:100,       startvalue = NGRID,   format = n -> "$n × $n"))
    sl_P, sl_W, sl_npc, sl_grid = sg.sliders

    cb = ctrl[2, 1] = GridLayout(halign = :left)
    btn_init  = Button(cb[1, 1], label = "Initialise grid")
    btn_comp  = Button(cb[1, 2], label = "Compute g–x", buttoncolor = RGBf(0.80, 0.88, 1.0))
    btn_clear = Button(cb[1, 3], label = "Clear points")
    ct = ctrl[3, 1] = GridLayout(halign = :left)
    tg_fields = Toggle(ct[1, 1], active = true);  Label(ct[1, 2], "show fields",      halign = :left)
    tg_line   = Toggle(ct[1, 3], active = true);  Label(ct[1, 4], "exact boundaries", halign = :left)
    tg_auto   = Toggle(ct[1, 5], active = false); Label(ct[1, 6], "compute on click", halign = :left)

    status = Label(right[3, 1], "", halign = :left, tellwidth = false, color = :gray25)
    cursor = Label(right[4, 1], " ", halign = :left, tellwidth = false, color = :gray25)
    # rowsize!(right, 1, Relative(0.72))   (removed: the diagram now takes the height the controls leave)

    # g–x panels (top = highest T)
    axGX = [Axis(left[k, 1]; ylabel = "g [J/g]",
                 xticks = XTICKS,
                 xlabel = k == NSLOTS ? XLABEL : "")
            for k in 1:NSLOTS]
    linkxaxes!(axGX...)

    # pending point (cross)
    pend_pos = Observable(Point2f(NaN, NaN)); pend_vis = Observable(false)
    scatter!(axTX, pend_pos; marker = :xcross, markersize = 20, color = :black, visible = pend_vis)

    # ------------------------------------------------------------ functions
    "Initialise: scan the grid, refine the boundaries, draw the fields and the exact boundaries."
    function init_grid!()
        n   = round(Int, sl_grid.value[])
        Xs  = [X_MIN + (1 - X_MIN) * (i - 0.5) / n for i in 1:n]   # cell centres (wt fraction): the ends are not nodes
        Ts  = range(TLIM..., length = n)
        t0  = time()
        labels, asm = scan_grid(key, Xs, Ts)
        Xb, Tb      = refine_boundaries(key, Xs, Ts, labels)
        dt  = time() - t0
        Xl, Tl      = exact_boundaries(W[], P[], range(TLIM..., length = 600))
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
                        @sprintf("P = %.1f kbar, ΔW = %.1f kJ/mol, %d L pseudocompounds",
                                 P[], W[], length(phases[]) - length(SOLIDS))
    end

    "Draw the g–x diagram of panel k for bulk composition X (wt fraction SiO2) at temperature T."
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
        μ   = 1000 .* μmol ./ MOLAR_MASS               # intercepts of the tangent at Fo (X_MIN) and SiO2 (1)
        G   = 1000 * r.G / molar_mass(w2x(X))

        vlines!(ax, [X]; color = :gray50, linestyle = :dot)                              # bulk
        lines!(ax, xm[iL], g[iL]; color = LIQ_COLOR, linewidth = 2)                        # g of L
        length(iL) <= 41 && scatter!(ax, xm[iL], g[iL]; markersize = 6, color = LIQ_COLOR) # pseudocompounds
        hx, hg = lower_hull(xm, g)
        lines!(ax, hx, hg; color = c, linewidth = 3, linestyle = (:dot, :dense))           # G of system (all X)
        lines!(ax, [X_MIN, 1.0], μ; color = c, linestyle = :dash, linewidth = 1.5)        # tangent (μ)
        scatter!(ax, xm[iS], g[iS]; markersize = 16, strokewidth = 2, strokecolor = c,
                 color = [st[i] ? c : :white for i in iS])
        text!(ax, xm[iS], g[iS]; text = [p.name for p in SOLIDS], offset = (0, 12),
              align = (:center, :bottom), fontsize = 15)
        sL = filter(in(iL), r.stable)                                                    # stable liquid(s)
        isempty(sL) || scatter!(ax, xm[sL], g[sL]; markersize = 13, color = LIQ_COLOR,
                                strokewidth = 2, strokecolor = c)
        im = iL[argmax([minimum(abs.(xm[i] .- xm[iS])) for i in iL])]                     # far from the solids
        text!(ax, [xm[im]], [g[im]]; text = LIQ_NAME, offset = (0, 8), align = (:center, :bottom),
              fontsize = 16, font = :bold, color = LIQ_COLOR)
        scatter!(ax, [X], [G]; marker = :star5, markersize = 22, color = :gold, strokewidth = 1.5, strokecolor = :black)
        text!(ax, [X_MIN, 1.0], μ;
              text = [rich("μ", subscript(COMPONENTS[i]), @sprintf(" = %.1f J/g", round(μ[i], digits = 1) + 0.0)) for i in 1:2],
              align = [(:left, :top), (:right, :top)], offset = [(4, -6), (-4, -6)], fontsize = 13, color = c)

        lo, hi = extrema(vcat(g, μ)); pad = max(hi - lo, 1e-3)
        limits!(ax, X_MIN - 0.035, 1.035, lo - 0.25pad, hi + 0.25pad)

        amounts = join([@sprintf("%s %.0f wt%%", phase_label(p), 100p.amount) for p in pl], " + ")
        note = ""
        if liq_only
            note = "   |   μ: tangent to L at its composition"
        elseif length(r.stable) < length(COMPONENTS)
            note = "   |   1 phase: μ not unique (any line through it below the others)"
        else
            # a solid that is not stable but touches the tangent: we sit on the eutectic or peritectic
            on_line = [SOLIDS[i].name for i in iS
                       if !st[i] && abs(g[i] - (μ[1] * (1 - xm[i]) + μ[2] * (xm[i] - X_MIN)) / (1 - X_MIN)) < 1e-6 * (1 + abs(g[i]))]
            isempty(on_line) || (note = "   |   invariant (eutectic or peritectic): " * join(vcat(phase_label.(pl), on_line), ", ") * " coexist")
        end
        ax.title    = rich(rich("T", subscript(string(rank)), color = c),
                           @sprintf(" = %.0f °C,  bulk %.1f wt%% %s,  P = %.1f kbar", T, 100X, COMPONENTS[2], P[]))
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
                limits!(axGX[k], X_MIN - 0.035, 1.035, 0, 1)
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

    "Programmatic access (used for testing / screenshots): X = weight fraction SiO2."
    function add_point!(X, T)
        pending[] = (X, T); compute!()
    end

    "New ΔW or number of pseudocompounds: new phases, new LP."
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
            X = clamp(Float64(event.data[1]), X_MIN, 1.0)
            T = clamp(Float64(event.data[2]), TLIM...)
            pending[] = (X, T)
            pend_pos[] = Point2f(X, T); pend_vis[] = true
            pl = phase_list(phases[], equil(X, T))
            status.text[] = @sprintf("Picked %.1f wt%% %s, T = %.0f °C  (%s)  – press Compute g–x", 100X, COMPONENTS[2],
                                     T, join(phase_label.(pl), " + "))
            tg_auto.active[] && compute!()
        end
        return Consume(false)
    end

    on(events(fig).mouseposition) do _
        if is_mouseinside(axTX.scene)
            X, T = mouseposition(axTX)
            if X_MIN <= X <= 1 && TLIM[1] <= T <= TLIM[2]
                pl = phase_list(phases[], equil(X, T))
                cursor.text[] = @sprintf("cursor: %.1f wt%% %s, T = %.0f °C  →  %s", 100X, COMPONENTS[2],
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
