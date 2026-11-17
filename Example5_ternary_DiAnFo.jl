#=
Example 5 – Interactive ternary Di–An–Fo (diopside – anorthite – forsterite)
=============================================================================

Three components (Di = CaMgSi2O6, An = CaAl2Si2O8, Fo = Mg2SiO4) and four phases:

    Di(s), An(s), Fo(s)   pure solids
    L                     liquid solution of Di(l), An(l), Fo(l)

MOLES IN THE CALCULATION, WEIGHT IN THE PLOTS. The thermodynamics is written in
moles: phase compositions are MOLE fractions x of Di, An and Fo, g is in kJ per
mole of these components, and the linear programme (the same as in Examples
0–4, gibbs_lp.jl) balances moles; its duals are the chemical potentials in
kJ/mol. The Di–An–Fo diagram of the course (part 3) is in WEIGHT per cent, so
everything that is shown is converted for the display:

    compositions  wt fractions  w_i = x_i M_i / Σ_j x_j M_j
    phase amounts mass fractions (lever rule in the wt% triangle)
    g, μ          per gram:  g [J/g] = g [kJ/mol] / M̄,  μ_i [J/g] = μ_i [kJ/mol] / M_i

A plane in (x, g per mole) is also a plane in (w, g per gram), so the tangent
plane stays a plane over the wt% triangle. The liquid g surface is discretised
into pseudocompounds of fixed composition on a triangular grid (step 1/n in
wt fraction, i.e. evenly spaced in the triangle that is drawn).

Liquid model (per mole of "mixing units"; component i contributes q_i units per
formula unit, so the unit fractions are y_i = q_i x_i / Σ_j q_j x_j):

    G_L = Σ y_i [ΔH_i (1 - T/T_m,i) + ΔV_i P] / q_i              (pure liquids; solids = 0)
        + RT Σ y_i ln y_i                                        (ideal mixing of the units)
        + Σ_m (W_H,m - T W_S,m) Π_i y_i^e_m,i                    (subregular excess terms)

One mole of components (Di, An, Fo) holds Σ_j q_j x_j moles of units; g of the
liquid per mole of components is G_L times this number.

The parameters were CALIBRATED against the Di–An–Fo liquidus diagram (1 atm):
melting points (Di 1392, An 1553, Fo 1890 °C), the three binary eutectics, the
ternary eutectic (1270 °C, Di49 An44 Fo07), the crystallization path of the
bulk Di19 An30 Fo51 and the liquidus isotherms. With the parameters below:

                     model                     diagram
    Di–An eutectic   1274 °C, An 41.7 wt%      1274 °C, An 42 wt%
    Di–Fo eutectic   1382 °C, Fo 12.1 wt%      1387 °C, Fo 12.4 wt%
    An–Fo eutectic   1329 °C, Fo 28.3 wt%      1327 °C, Fo 28.3 wt%
    ternary E        1269 °C, Di49 An44 Fo07   1270 °C, Di49 An44 Fo07
    all 27 targets:  rms 6.9 K, max 27 K; the liquid is stable (no unmixing)
                     everywhere above the liquidus.

ΔH, q and W are *effective* (fitted) parameters, not measured enthalpies of
fusion. The pressure dependence uses the melting slopes dT_m/dP below; it is
illustrative (linear) and ignores the new phases that appear at high P.

How to use
----------
* Right: the ternary (Di bottom left, Fo bottom right, An top).
  - "isothermal section": the stable assemblage at the chosen T and P on a
    grid of bulk compositions, solved by LP (black dots: boundaries found by
    bisection between neighbouring cells).
  - "liquidus projection": primary phase fields, cotectics, eutectics and
    liquidus isotherms of the model (computed analytically).
  - "reference": the cotectics, eutectics and crystallization path of the
    diagram in part 3 (red), to compare with the model.
* Left-click in the triangle (or type "Di An Fo" in wt% into the box) to choose
  the bulk composition. The T and P sliders set the conditions.
* Left: the g surface of the liquid over the ternary base (3D, drag to rotate;
  g per gram over the wt% triangle), the solids (at g = 0, the reference), the
  tangent plane through the stable phases (its intercepts at the corners are
  μ_Di, μ_An, μ_Fo per gram; the molar values are printed below) and the bulk
  composition (star). In the triangle the coexisting phases are connected
  (tie-line or tie-triangle).

Run from a terminal:      julia Example5_ternary_DiAnFo.jl
or from the Julia REPL:   include("Example5_ternary_DiAnFo.jl")
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
const COMPONENTS = ["Di", "An", "Fo"]
const MOLAR_MASS = [216.55, 278.21, 140.69]         # CaMgSi2O6, CaAl2Si2O8, Mg2SiO4 [g/mol]
const TM_C       = [1392.0, 1553.0, 1890.0]         # melting points at 1 atm [°C]
const DTDP       = [9.0, 5.0, 4.7]                  # melting slopes dT_m/dP [K/kbar] (illustrative)

# calibrated liquid model (see header)
const Q   = [0.85, 1.80, 0.70]                      # mixing units per formula unit
const DH  = [328.4, 52.0, 118.9]                    # effective ΔH of fusion [kJ per mol formula]
#  excess terms Π y_i^e_i  (Di, An, Fo exponents)       W_H [kJ/mol]   W_S [kJ/mol/K]
const WEX = [([2, 1, 0], -15.75, 0.0320),           # y_Di² y_An
             ([1, 2, 0], 225.05, 0.1387),           # y_Di  y_An²
             ([2, 0, 1], -16.45, 0.0),              # y_Di² y_Fo
             ([1, 0, 2], -24.68, 0.0),              # y_Di  y_Fo²
             ([0, 2, 1], -22.53, 0.0),              # y_An² y_Fo
             ([0, 1, 2], -64.08, 0.0),              # y_An  y_Fo²
             ([1, 1, 1], -44.51, 0.0)]              # y_Di  y_An y_Fo

const LIQ_NAME  = "L"
const NPC_START = 40                                # liquid grid: step 1/n in wt fraction
const NSEC      = 24                                # isothermal section: cells per side
const T_START   = 1450.0                            # [°C]
const P_START   = 0.0                               # [kbar]
const X_START   = [0.19, 0.30, 0.51]                # bulk Di, An, Fo (crystallization path in part 3)
const TLIM_C    = (1200.0, 1900.0)
const ISOTHERMS = 1300:100:1800                     # liquidus isotherms of the projection [°C]

# reference: Di–An–Fo diagram of part 3 (wt fractions Di, An, Fo)
const REF_E       = [0.493, 0.439, 0.068]           # ternary eutectic, 1270 °C
const REF_BIN     = [([0.58, 0.42, 0.0], 1274), ([0.876, 0.0, 0.124], 1387), ([0.0, 0.717, 0.283], 1327)]
const REF_COTECT  = [[[0.58, 0.42, 0.0], REF_E],
                     [[0.876, 0.0, 0.124], [0.63, 0.29, 0.08], REF_E],
                     [[0.0, 0.717, 0.283], [0.34, 0.53, 0.13], [0.41, 0.49, 0.10], REF_E]]
const REF_PATH    = [([0.19, 0.31, 0.50], 1610), ([0.22, 0.35, 0.43], 1550), ([0.24, 0.38, 0.38], 1500),
                     ([0.27, 0.43, 0.30], 1450), ([0.32, 0.50, 0.18], 1350), ([0.34, 0.53, 0.13], 1290),
                     ([0.41, 0.49, 0.10], 1280)]

const PRIMARY_COLORS = ["#cfe0f3", "#fbd9bf", "#ebe99c"]             # Di, An, Fo fields
const FIELD_COLORS   = Dict("L" => "#d9cbe6", "Fo + L" => "#ebe99c", "An + L" => "#fbd9bf",
                            "Di + L" => "#cfe0f3", "An + Fo + L" => "#f4b183",
                            "Di + Fo + L" => "#9fd8c8", "Di + An + L" => "#f2a7b8",
                            "Di + An + Fo" => "#c8c8c8")
const EXTRA_COLORS   = ["#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3"]
const SEL_COLOR      = :royalblue3
const LIQ_COLOR      = :firebrick

# =============================================================================
# THERMODYNAMICS (independent of the plotting)
# =============================================================================
const TM_K = TM_C .+ 273.15
const DS   = DH ./ TM_K                             # kJ/mol/K
const DV   = DTDP .* DS                             # kJ/kbar/mol  (Clausius–Clapeyron)

K(T_C) = T_C + 273.15

# ----- moles <-> weight. The calculation uses the mole fractions x of Di, An, Fo,
#       the plots use their weight fractions w.
"weight fractions -> mole fractions"
w2x(w) = (n = w ./ MOLAR_MASS; n ./ sum(n))
"mole fractions -> weight fractions"
x2w(x) = (m = x .* MOLAR_MASS; m ./ sum(m))
"mass [g] of one mole of components with mole fractions x"
molar_mass(x) = dot(x, MOLAR_MASS)

"mole fractions x -> unit fractions y (component i contributes q_i mixing units per formula unit)"
unitfrac(x) = (n = x .* Q; n ./ sum(n))

mono(y, e) = y[1]^e[1] * y[2]^e[2] * y[3]^e[3]
dmono(y, e, k) = e[k] == 0 ? 0.0 : e[k] * mono(y, ntuple(i -> e[i] - (i == k), 3))
"∂(n G_ex)/∂n_k for one homogeneous monomial (per unit W)"
exbasis(y, e, k) = dmono(y, e, k) - (sum(e) - 1) * mono(y, e)

"Liquid chemical potentials [kJ/mol] of Di, An, Fo at mole fractions x (solids = 0)."
function liquid_mu(x, T_C, P)
    T = K(T_C); y = unitfrac(x)
    map(1:3) do k
        DH[k] * (1 - T / TM_K[k]) + DV[k] * P +
        Q[k] * (sum((wh - T * ws) * exbasis(y, e, k) for (e, wh, ws) in WEX) + RGAS * T * log(max(y[k], 1e-300)))
    end
end

"Liquidus temperature [°C] of solid k for a liquid of WEIGHT composition w (closed form)."
function tliq(w, k, P)
    y  = unitfrac(w2x(w))
    yk = max(y[k], 1e-300)
    H  = DH[k] / Q[k] + DV[k] * P / Q[k] + sum(wh * exbasis(y, e, k) for (e, wh, ws) in WEX)
    S  = DS[k] / Q[k] - RGAS * log(yk) + sum(ws * exbasis(y, e, k) for (e, wh, ws) in WEX)
    return H / S - 273.15
end
tliq(w, P) = maximum(tliq(w, k, P) for k in 1:3)

"""
Liquid pseudocompounds on the triangular grid of step 1/n in WEIGHT fraction (the nodes of
the triangle that is drawn). Each one is an ordinary `Phase` for the LP: composition in
mole fractions, g = g0 - sT + vP in kJ per mole of components (T in K, P in kbar).
"""
function liquid_pseudocompounds(n)
    nodes, _, _ = lattice(n)
    map(nodes) do w
        x  = w2x(w)
        y  = unitfrac(x)
        nu = dot(Q, x)                                              # mol of units per mol of components
        h  = dot(x, DH) + nu * sum(wh * mono(y, e) for (e, wh, ws) in WEX)
        s  = dot(x, DS) + nu * (-RGAS * sum(xlogx, y) + sum(ws * mono(y, e) for (e, wh, ws) in WEX))
        v  = dot(x, DV)
        Phase(LIQ_NAME, x, h, s, v)
    end
end

const SOLIDS = [Phase(COMPONENTS[k], [Float64(i == k) for i in 1:3], 0.0, 0.0, 0.0) for k in 1:3]

# ----------------------------------------------------------------------------
# triangular lattice: nodes (i, j) with w = (1 - (i+j)/N, i/N, j/N), i+j <= N
# ----------------------------------------------------------------------------
function lattice(N)
    idx   = Dict{Tuple{Int,Int},Int}()
    nodes = Vector{Vector{Float64}}()
    ij    = Tuple{Int,Int}[]
    for i in 0:N, j in 0:N-i
        push!(nodes, [1 - (i + j) / N, i / N, j / N]); push!(ij, (i, j)); idx[(i, j)] = length(nodes)
    end
    tris = NTuple{3,Int}[]
    for i in 0:N-1, j in 0:N-1-i
        push!(tris, (idx[(i, j)], idx[(i+1, j)], idx[(i, j+1)]))
        i + j <= N - 2 && push!(tris, (idx[(i+1, j)], idx[(i, j+1)], idx[(i+1, j+1)]))
    end
    return nodes, ij, tris
end

"cartesian position: Di (0,0), Fo (1,0), An (1/2, √3/2)"
tocart(w) = Point2f(w[3] + 0.5 * w[2], sqrt(3) / 2 * w[2])
function fromcart(x, y)
    a = y / (sqrt(3) / 2); f = x - 0.5 * a
    w = max.([1 - a - f, a, f], 0.0)
    return w ./ sum(w)
end

"g of the continuous liquid per gram [J/g] at WEIGHT composition w"
gliq(w, T_C, P) = (x = w2x(w); 1000 * dot(x, liquid_mu(x, T_C, P)) / molar_mass(x))

"""
Group the stable phases into physical phases. Two stable liquid pseudocompounds
belong to the same liquid if the continuous g surface between them lies below
their chord (no hump, i.e. no miscibility gap): on a discrete grid the LP may
pick pseudocompounds that are not direct neighbours. If there is a hump, they
are two different liquids (L + L). (The test is made per gram along a straight
line in wt fractions; this is equivalent to the test per mole in mole fractions.)
"""
function ternary_groups(phases, stable, T_C, P)
    nS     = length(SOLIDS)
    sol    = [[i] for i in stable if i <= nS]
    groups = [[i] for i in stable if i > nS]
    merged = true
    while merged
        merged = false
        for a in 1:length(groups)-1, b in a+1:length(groups)
            wa = sum(x2w(phases[i].comp) for i in groups[a]) ./ length(groups[a])
            wb = sum(x2w(phases[i].comp) for i in groups[b]) ./ length(groups[b])
            hump = maximum(gliq((1 - t) .* wa .+ t .* wb, T_C, P) -
                           ((1 - t) * gliq(wa, T_C, P) + t * gliq(wb, T_C, P)) for t in 0.1:0.1:0.9)
            if hump <= 1e-6
                append!(groups[a], groups[b]); deleteat!(groups, b); merged = true
                break
            end
        end
    end
    return vcat(sol, groups)
end

"Name, amount (MASS fraction of the system), wt composition w and mole composition x of every stable phase."
function phase_list(phases, r, T_C, P)
    pl = map(ternary_groups(phases, r.stable, T_C, P)) do g
        n = sum(r.α[g])                                        # moles of this phase
        x = sum(r.α[i] .* phases[i].comp for i in g) ./ n
        (name = phases[g[1]].name, mass = n * molar_mass(x), x = x)
    end
    mtot = sum(p.mass for p in pl)
    return [(name = p.name, amount = p.mass / mtot, w = x2w(p.x), x = p.x) for p in pl]
end
assemblage_key(pl) = join(sort([p.name for p in pl], by = n -> n == LIQ_NAME ? "~" : n), " + ")
pct(w) = @sprintf("Di%.0f An%.0f Fo%.0f", 100w[1], 100w[2], 100w[3])
phase_label(p) = p.name == LIQ_NAME ? "L(" * pct(p.w) * ")" : p.name

# ----------------------------------------------------------------------------
# analytic liquidus projection
# ----------------------------------------------------------------------------
"Segments where the field `f` (value per node) crosses `level` (marching triangles)."
function iso_segments(nodes, vals, tris, level)
    segs = Tuple{Vector{Float64},Vector{Float64}}[]
    for t in tris
        cr = Vector{Float64}[]
        for (a, b) in ((t[1], t[2]), (t[2], t[3]), (t[3], t[1]))
            va, vb = vals[a] - level, vals[b] - level
            if (va < 0) != (vb < 0)
                s = va / (va - vb); push!(cr, nodes[a] .+ s .* (nodes[b] .- nodes[a]))
            end
        end
        length(cr) == 2 && push!(segs, (cr[1], cr[2]))
    end
    return segs
end
seg_points(segs) = reduce(vcat, [[tocart(a), tocart(b), Point2f(NaN, NaN)] for (a, b) in segs]; init = Point2f[])

function binary_eutectic(i, j, P)
    w(x) = (v = zeros(3); v[i] = 1 - x; v[j] = x; v)
    f(x) = tliq(w(x), i, P) - tliq(w(x), j, P)
    a, b = 1e-6, 1 - 1e-6
    for _ in 1:60
        m = (a + b) / 2
        sign(f(m)) == sign(f(a)) ? (a = m) : (b = m)
    end
    return w(a), tliq(w(a), i, P)
end

function ternary_eutectic(P; n = 200)
    best, bw = Inf, X_START
    for i in 1:n-1, j in 1:n-1-i
        w  = [1 - (i + j) / n, i / n, j / n]
        Ts = (tliq(w, 1, P), tliq(w, 2, P), tliq(w, 3, P))
        s  = maximum(Ts) - minimum(Ts)
        s < best && (best = s; bw = w)
    end
    return bw, tliq(bw, P)
end

function liquidus_projection(P; N = 90)
    nodes, _, tris = lattice(N)
    Tk = [[tliq(w, k, P) for w in nodes] for k in 1:3]
    Tl = max.(Tk...)
    prim = [argmax([tliq(sum(nodes[collect(t)]) ./ 3, k, P) for k in 1:3]) for t in tris]
    iso = [(L, seg_points(iso_segments(nodes, Tl, tris, L))) for L in ISOTHERMS]
    cot = Point2f[]
    for (i, j) in ((1, 2), (1, 3), (2, 3))
        k = 6 - i - j
        for (a, b) in iso_segments(nodes, Tk[i] .- Tk[j], tris, 0.0)
            m = (a .+ b) ./ 2
            tliq(m, i, P) >= tliq(m, k, P) && push!(cot, tocart(a), tocart(b), Point2f(NaN, NaN))
        end
    end
    return (; nodes, tris, prim, iso, cot, Tl)
end

# =============================================================================
# APPLICATION
# =============================================================================
function build_app()
    # ---------------------------------------------------------------- state
    P      = Ref(P_START)
    T      = Ref(T_START)
    X      = Ref(copy(X_START))                       # bulk composition in WEIGHT fractions (as on the triangle)
    npc    = Ref(NPC_START)
    liq    = Ref(liquid_pseudocompounds(NPC_START))
    phases = Ref(vcat(SOLIDS, liq[]))
    lp     = Ref(GibbsLP(phases[], w2x(X[])))
    color_of = Dict{String,String}()
    secplts, projplts, selplts = Any[], Any[], Any[]
    debounce = Ref(0)

    # w = bulk in weight fractions; the LP is solved in moles (result in moles and kJ/mol)
    equil(w, T_C) = (set_bulk!(lp[], w2x(w)); solve!(lp[], P[], K(T_C)))
    fieldcolor(key) = get(FIELD_COLORS, key) do
        get!(color_of, key, EXTRA_COLORS[mod1(length(color_of) + 1, length(EXTRA_COLORS))])
    end

    # ---------------------------------------------------------------- layout
    fig = Figure(size = (1650, 980), fontsize = 15)
    Label(fig[0, 1:2], "Example 5 – ternary Di–An–Fo (Di, An, Fo, liquid L): " *
          "Gibbs energy minimisation by linear programming", fontsize = 20, font = :bold, word_wrap = true, tellwidth = false)
    left  = fig[1, 1] = GridLayout()
    right = fig[1, 2] = GridLayout()
    colsize!(fig.layout, 2, Relative(0.47))

    ax3 = Axis3(left[1, 1]; xlabel = "", ylabel = "", zlabel = "g [J/g]",
                aspect = (1, 0.866, 0.9), azimuth = 1.35π, elevation = 0.10π,
                xticksvisible = false, yticksvisible = false,
                xticklabelsvisible = false, yticklabelsvisible = false,
                xgridvisible = false, ygridvisible = false)
    info = Label(left[2, 1], " ", halign = :left, tellwidth = false, justification = :left)
    rowsize!(left, 1, Relative(0.88))

    axT = Axis(right[1, 1]; aspect = DataAspect(), title = "isothermal section  (left-click to pick a bulk composition)")
    hidedecorations!(axT); hidespines!(axT)
    limits!(axT, -0.24, 1.24, -0.12, 0.95)      # side margins hold the corner labels also in a small window
    for name in (:rectanglezoom, :scrollzoom, :dragpan)
        name in keys(interactions(axT)) && deregister_interaction!(axT, name)
    end

    # one sub-grid per row, so that the widths of the rows do not add up (keeps the controls inside the window)
    ctrl = right[2, 1] = GridLayout(tellwidth = false)
    sg = SliderGrid(ctrl[1, 1],
        (label = "T [°C]",            range = TLIM_C[1]:1.0:TLIM_C[2], startvalue = T_START, format = "{:.0f}"),
        (label = "P [kbar]",          range = 0.0:0.5:20.0,            startvalue = P_START, format = "{:.1f}"),
        (label = "L grid (1/n)",      range = 10:2:80,                 startvalue = NPC_START, format = n -> "n = $n"),
        (label = "section cells",     range = 6:2:60,                  startvalue = NSEC,    format = n -> "$n per side"))
    sl_T, sl_P, sl_npc, sl_sec = sg.sliders
    cb = ctrl[2, 1] = GridLayout(halign = :left)
    Label(cb[1, 1], "bulk wt%  Di An Fo:", halign = :right)
    tb = Textbox(cb[1, 2], placeholder = "19 30 51", width = 160)
    btn_sec = Button(cb[1, 3], label = "Update section")
    ct = ctrl[3, 1] = GridLayout(halign = :left)
    tg_sec  = Toggle(ct[1, 1], active = true);  Label(ct[1, 2], "isothermal section",  halign = :left)
    tg_proj = Toggle(ct[1, 3], active = true);  Label(ct[1, 4], "liquidus projection", halign = :left)
    tg_ref  = Toggle(ct[1, 5], active = false); Label(ct[1, 6], "reference (part 3)",  halign = :left)
    status = Label(right[3, 1], " ", halign = :left, tellwidth = false, color = :gray25)
    cursor = Label(right[4, 1], " ", halign = :left, tellwidth = false, color = :gray25)
    # rowsize!(right, 1, Relative(0.7))   (removed: the section now takes the height the controls leave)

    # ------------------------------------------------- static ternary decorations
    corners = [tocart(w) for w in ([1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0])]
    lines!(axT, [corners; corners[1:1]]; color = :black, linewidth = 2)
    for f in 0.1:0.1:0.9                                     # 10 % ticks on the edges
        for (a, b) in ((1, 3), (3, 2), (2, 1))
            w = zeros(3); w[a] = 1 - f; w[b] = f
            scatter!(axT, [tocart(w)]; marker = :circle, markersize = 3, color = :black)
        end
    end
    vlabel = text!(axT, corners; text = Observable(fill("", 3)),
                   align = [(:right, :top), (:center, :bottom), (:left, :top)],
                   offset = [(-4, -4), (0, 6), (4, -4)], fontsize = 16, font = :bold)

    # reference overlay (part 3 diagram)
    refpl = Any[]
    for c in REF_COTECT
        push!(refpl, lines!(axT, tocart.(c); color = :red, linewidth = 2, linestyle = :dash))
    end
    push!(refpl, scatter!(axT, vcat([tocart(REF_E)], [tocart(w) for (w, _) in REF_BIN]);
                          marker = :diamond, markersize = 12, color = :red))
    push!(refpl, scatterlines!(axT, [tocart(w) for (w, _) in REF_PATH]; color = :red, markersize = 7,
                               linestyle = :dot))
    push!(refpl, text!(axT, [tocart(w) for (w, _) in REF_PATH]; text = [string(t) for (_, t) in REF_PATH],
                       fontsize = 11, color = :red, offset = (6, 0), align = (:left, :center)))
    foreach(p -> translate!(p, 0, 0, 5), refpl)
    on(tg_ref.active; update = true) do v
        foreach(p -> p.visible = v, refpl)
    end

    # ------------------------------------------------------------ functions
    "Isothermal section at T: LP on the centroids of a triangular mesh of bulk compositions."
    function update_section!()
        foreach(p -> delete!(axT, p), secplts); empty!(secplts)
        N = round(Int, sl_sec.value[])
        nodes, _, tris = lattice(N)
        cents = [sum(nodes[collect(t)]) ./ 3 for t in tris]
        t0 = time()
        keyat(w) = assemblage_key(phase_list(phases[], equil(w, T[]), T[], P[]))
        keys_ = [keyat(w) for w in cents]
        # neighbours: triangles sharing an edge; bisect where the assemblage changes
        owner = Dict{Tuple{Int,Int},Int}(); bnd = Point2f[]
        nlp = length(cents)
        for (m, t) in enumerate(tris), (a, b) in ((t[1], t[2]), (t[2], t[3]), (t[3], t[1]))
            e = minmax(a, b)
            if haskey(owner, e)
                o = owner[e]
                if keys_[o] != keys_[m]
                    w1, w2, k1 = cents[o], cents[m], keys_[o]
                    for _ in 1:10
                        wm = (w1 .+ w2) ./ 2
                        keyat(wm) == k1 ? (w1 = wm) : (w2 = wm)
                    end
                    nlp += 10
                    push!(bnd, tocart((w1 .+ w2) ./ 2))
                end
            else
                owner[e] = m
            end
        end
        dt = time() - t0
        # coloured cells
        V = Vector{Point2f}(undef, 3length(tris)); F = zeros(Int, length(tris), 3); C = RGBAf[]
        for (m, t) in enumerate(tris)
            for (c, n) in enumerate(t)
                V[3m-3+c] = tocart(nodes[n]); F[m, c] = 3m - 3 + c
                push!(C, RGBAf(parse(RGBf, fieldcolor(keys_[m])), 0.85))
            end
        end
        msh = mesh!(axT, V, F; color = C, shading = NoShading, visible = tg_sec.active)
        translate!(msh, 0, 0, -10)
        bd = scatter!(axT, bnd; markersize = 4, color = :black, visible = tg_sec.active)
        push!(secplts, msh, bd)
        # field labels at the centre of every connected region of >= 3 cells
        for key in unique(keys_)
            ids = findall(==(key), keys_)
            length(ids) >= 3 || continue
            c = sum(tocart(cents[i]) for i in ids) / length(ids)
            push!(secplts, text!(axT, c; text = key, align = (:center, :center), fontsize = 15,
                                 font = :bold, color = :gray15, visible = tg_sec.active))
        end
        # analytic liquidus isotherm at T (the L / solid + L boundary must sit on it)
        nl, _, tl = lattice(60)
        iso = seg_points(iso_segments(nl, [tliq(w, P[]) for w in nl], tl, T[]))
        push!(secplts, lines!(axT, iso; color = LIQ_COLOR, linewidth = 2.5, linestyle = :dash,
                              visible = tg_sec.active))
        axT.title = @sprintf("isothermal section at %.0f °C, %.1f kbar  (left-click to pick a bulk composition)", T[], P[])
        status.text[] = @sprintf("Section: %d cells, %d boundary points, ≈ %d LPs (%d phases each) in %.2f s",
                                 length(tris), length(bnd), nlp, length(phases[]), dt)
    end

    "Liquidus projection of the model at pressure P."
    function update_projection!()
        foreach(p -> delete!(axT, p), projplts); empty!(projplts)
        pr = liquidus_projection(P[])
        V = Vector{Point2f}(undef, 3length(pr.tris)); F = zeros(Int, length(pr.tris), 3); C = RGBAf[]
        for (m, t) in enumerate(pr.tris)
            for (c, n) in enumerate(t)
                V[3m-3+c] = tocart(pr.nodes[n]); F[m, c] = 3m - 3 + c
                push!(C, RGBAf(parse(RGBf, PRIMARY_COLORS[pr.prim[m]]), 1.0))
            end
        end
        tint = mesh!(axT, V, F; color = C, shading = NoShading,
                     visible = lift((p, s) -> p && !s, tg_proj.active, tg_sec.active))
        translate!(tint, 0, 0, -11)
        push!(projplts, tint)
        for (L, pts) in pr.iso
            isempty(pts) && continue
            push!(projplts, lines!(axT, pts; color = (:black, 0.6), linestyle = :dash, linewidth = 1,
                                   visible = tg_proj.active))
            i = findfirst(p -> !isnan(p[1]), pts[div(length(pts), 2):end])
            i === nothing || push!(projplts, text!(axT, pts[div(length(pts), 2) - 1 + i]; text = string(L),
                                                   fontsize = 12, color = :gray20, visible = tg_proj.active))
        end
        push!(projplts, lines!(axT, pr.cot; color = :black, linewidth = 2.5, visible = tg_proj.active))
        eus = [binary_eutectic(1, 2, P[]), binary_eutectic(1, 3, P[]), binary_eutectic(2, 3, P[]), ternary_eutectic(P[])]
        push!(projplts, scatter!(axT, [tocart(w) for (w, _) in eus]; markersize = 10, color = :black,
                                 visible = tg_proj.active))
        push!(projplts, text!(axT, [tocart(w) for (w, _) in eus]; text = [@sprintf("%.0f", t) for (_, t) in eus],
                              fontsize = 13, font = :bold, offset = [(-6, 4), (0, -14), (8, 0), (8, -12)],
                              align = [(:right, :bottom), (:center, :top), (:left, :center), (:left, :top)],
                              visible = tg_proj.active))
        tm = [TM_C[k] + DTDP[k] * P[] for k in 1:3]
        vlabel.text[] = [@sprintf("Di %.0f", tm[1]), @sprintf("An %.0f", tm[2]), @sprintf("Fo %.0f", tm[3])]
        Ew, ET = eus[4]
        return Ew, ET
    end

    "g surface, tangent plane and phases for the current bulk composition, T and P."
    function update_selection!()
        foreach(p -> delete!(axT, p), selplts); empty!(selplts)
        empty!(ax3)
        w, T_C = X[], T[]
        r  = equil(w, T_C)                             # solved in moles: r.g, r.μ, r.G in kJ/mol
        pl = phase_list(phases[], r, T_C, P[])
        nS = length(SOLIDS)
        iL = [p for p in pl if p.name == LIQ_NAME]
        # μ: a single liquid has a unique tangent plane (continuous model); with solids the LP duals
        μmol = length(pl) == 1 && length(iL) == 1 && all(iL[1].w .> 1e-9) ? liquid_mu(iL[1].x, T_C, P[]) : r.μ   # [kJ/mol]
        # ---- display in weight: g and μ per gram [J/g]
        g  = 1000 .* r.g ./ [molar_mass(p.comp) for p in phases[]]
        μ  = 1000 .* μmol ./ MOLAR_MASS                # heights of the tangent plane at the corners
        G  = 1000 * r.G / molar_mass(w2x(w))

        # liquid surface
        nodes, _, tris = lattice(npc[])
        gl = g[nS+1:end]
        V  = [Point3f(tocart(n)..., gl[i]) for (i, n) in enumerate(nodes)]
        F  = permutedims(reduce(hcat, [collect(t) for t in tris]))
        mesh!(ax3, V, F; color = gl, colormap = :viridis, alpha = 0.8, transparency = true)
        # base triangle, prism edges and corner labels
        zlo = min(minimum(gl), minimum(μ), 0.0); zhi = max(maximum(gl), maximum(μ), 0.0)
        pad = 0.08 * (zhi - zlo + 1)
        zb  = zlo - pad
        c3(w, z) = Point3f(tocart(w)..., z)
        E   = ([1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0])
        lines!(ax3, [c3(e, zb) for e in (E..., E[1])]; color = :black, linewidth = 1.5)
        for (k, e) in enumerate(E)
            lines!(ax3, [c3(e, zb), c3(e, zhi + pad)]; color = (:gray40, 0.6), linestyle = :dot)
            text!(ax3, c3(e, zb); text = COMPONENTS[k], fontsize = 18, font = :bold,
                  align = (:center, :top), offset = (0, -6))
        end
        # tangent plane through μ
        Vp = [Point3f(tocart(E[k])..., μ[k]) for k in 1:3]
        mesh!(ax3, Vp, [1 2 3]; color = (SEL_COLOR, 0.35), transparency = true)
        lines!(ax3, [c3(E[k], μ[k]) for k in (1, 2, 3, 1)]; color = SEL_COLOR, linewidth = 1.5, linestyle = :dash)
        text!(ax3, [c3(E[k], μ[k]) for k in 1:3];
              text = [rich("μ", subscript(COMPONENTS[k]), @sprintf(" = %.1f", round(μ[k], digits = 1) + 0.0)) for k in 1:3],
              color = SEL_COLOR, fontsize = 13, offset = (6, 4))
        # solids (g = 0), filled when stable
        st = [any(p -> p.name == COMPONENTS[k], pl) for k in 1:3]
        scatter!(ax3, [c3(E[k], 0.0) for k in 1:3]; markersize = 16, strokewidth = 2, strokecolor = :black,
                 color = [s ? :black : :white for s in st])
        # stable phases on the plane, tie-line / tie-triangle, bulk
        pz = [c3(p.w, dot(p.w, μ)) for p in pl]
        length(pz) >= 2 && lines!(ax3, [pz; pz[1:1]]; color = SEL_COLOR, linewidth = 3)
        lq = [c3(p.w, dot(p.w, μ)) for p in iL]
        isempty(lq) || scatter!(ax3, lq; markersize = 14, color = LIQ_COLOR, strokewidth = 1.5, strokecolor = :black)
        lines!(ax3, [c3(w, zb), c3(w, G)]; color = :gray30, linestyle = :dot)
        scatter!(ax3, [c3(w, G)]; marker = :star5, markersize = 22, color = :gold, strokewidth = 1.5,
                 strokecolor = :black)
        limits!(ax3, -0.02, 1.02, -0.02, 0.89, zb, zhi + pad)

        amounts = join([@sprintf("%s %.0f wt%%", phase_label(p), 100p.amount) for p in pl], " + ")
        note = length(r.stable) < 3 && isempty(iL) ? "\n(fewer than 3 phases and no liquid: μ not unique, LP shows one choice)" : ""
        ax3.title = @sprintf("g at %.0f °C, %.1f kbar, bulk %s", T_C, P[], pct(w))
        info.text[] = "stable: " * amounts *
                      @sprintf("\nG = %.2f J/g  |  μ = (%.1f, %.1f, %.1f) J/g = (%.2f, %.2f, %.2f) kJ/mol", G,
                               (round.(μ, digits = 1) .+ 0.0)..., (round.(μmol, digits = 2) .+ 0.0)...) * note

        # ternary: bulk, phases, tie-lines
        pc = [tocart(p.w) for p in pl]
        length(pc) >= 2 && push!(selplts, lines!(axT, [pc; pc[1:1]]; color = SEL_COLOR, linewidth = 3))
        push!(selplts, scatter!(axT, pc; markersize = 11, color = SEL_COLOR),
              scatter!(axT, [tocart(w)]; marker = :xcross, markersize = 18, color = SEL_COLOR))
        foreach(p -> translate!(p, 0, 0, 10), selplts)
    end

    "Debounced section update while a slider is being dragged."
    function schedule_section!()
        debounce[] += 1; id = debounce[]
        @async try
            sleep(0.35)
            id == debounce[] && update_section!()
        catch err
            @error "section update failed" exception = (err, catch_backtrace())
        end
    end

    function rebuild!()
        npc[] = round(Int, sl_npc.value[])
        liq[] = liquid_pseudocompounds(npc[]); phases[] = vcat(SOLIDS, liq[])
        lp[]  = GibbsLP(phases[], w2x(X[]))
        update_selection!(); schedule_section!()
    end

    function set_bulk_point!(w)
        X[] = w ./ sum(w)
        update_selection!()
    end

    # ------------------------------------------------------------ callbacks
    register_interaction!(axT, :pick_point) do event::MouseEvent, axis
        if event.type === MouseEventTypes.leftclick
            set_bulk_point!(fromcart(event.data[1], event.data[2]))
        end
        return Consume(false)
    end
    on(events(fig).mouseposition) do _
        if is_mouseinside(axT.scene)
            x, y = mouseposition(axT)
            a = y / (sqrt(3) / 2); f = x - 0.5a
            if a >= 0 && f >= 0 && a + f <= 1
                w = fromcart(x, y)
                pl = phase_list(phases[], equil(w, T[]), T[], P[])
                cursor.text[] = @sprintf("cursor: %s  →  %s   |  liquidus %.0f °C", pct(w),
                                         join(phase_label.(pl), " + "), tliq(w, P[]))
            end
        end
        return Consume(false)
    end
    on(tb.stored_string) do s
        v = tryparse.(Float64, split(replace(s, "," => " ")))
        if length(v) == 3 && all(!isnothing, v) && all(>=(0), v) && sum(v) > 0
            set_bulk_point!(Float64.(v))
        else
            status.text[] = "Type three numbers: wt% Di An Fo, e.g. 19 30 51"
        end
    end
    on(sl_T.value) do v
        T[] = v; update_selection!(); schedule_section!()
    end
    on(sl_P.value) do v
        P[] = v; update_projection!(); update_selection!(); schedule_section!()
    end
    on(_ -> rebuild!(), sl_npc.value)
    on(_ -> schedule_section!(), sl_sec.value)
    on(_ -> update_section!(), btn_sec.clicks)

    # ------------------------------------------------------- initialisation
    update_projection!()
    update_section!()
    update_selection!()
    return fig, (set_bulk = set_bulk_point!, update_section = update_section!,
                 set_T = v -> set_close_to!(sl_T, v), set_P = v -> set_close_to!(sl_P, v),
                 ax3 = ax3, toggles = (section = tg_sec, projection = tg_proj, reference = tg_ref))
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
