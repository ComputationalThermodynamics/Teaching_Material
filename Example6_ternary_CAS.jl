#=
Example 6 – Interactive ternary CaO–Al2O3–SiO2 (CAS)
=====================================================

Three components (the oxides CaO, Al2O3, SiO2), a liquid and thirteen solids of
fixed composition (cement notation C = CaO, A = Al2O3, S = SiO2):

    Lime  CaO              Crn   corundum Al2O3       Crs / Trd  cristobalite / tridymite SiO2
    Pwo   CaSiO3 (CS)      C2S   Ca2SiO4              An    anorthite CaAl2Si2O8 (CAS2)
    Gh    gehlenite Ca2Al2SiO7 (C2AS)                 Mul   mullite Al6Si2O13 (A3S2)
    CA    CaAl2O4          CA2   CaAl4O7              CA6   CaAl12O19 (hibonite)
    C3A   Ca3Al2O6

The program is built like Example 5 (same linear programme, same layout), but the
end-members of the triangle are now oxides, and most solids are COMPOUNDS that lie
inside the triangle or on its edges. Some of them melt congruently (An, Gh, Pwo, C2S,
Mul, CA, CA2: each lies inside its own primary field), others incongruently (CA6, C3A:
they lie outside their field).

MOLES IN THE CALCULATION, WEIGHT IN THE PLOTS (as in Example 5). Phase compositions
are mole fractions of the oxides, g is in kJ per mole of oxides, and the duals of the
LP are the chemical potentials μ_CaO, μ_Al2O3, μ_SiO2 in kJ/mol. The diagram of the
course (data/Papers/CAS_Ternary.png) is in WEIGHT fractions, so everything shown is
converted: compositions in wt%, phase amounts in mass fractions, g and μ per gram.

Liquid model (per mole of "mixing units"; oxide i contributes q_i units, unit fractions
y_i = q_i x_i / Σ_j q_j x_j; the pure oxide LIQUIDS are the reference, g = 0):

    G_L = RT Σ y_i ln y_i + Σ_m (W_H,m - T W_S,m) Π_i y_i^e_m,i

One mole of oxides holds Σ_j q_j x_j moles of units; g of the liquid per mole of oxides
is G_L times this number.

Solids: every solid s (oxide mole fractions c_s) is defined from the liquid of the same
composition and a fusion term:

    g_s(T) = g_L(c_s, T) - ΔH_s (1 - T/T_m,s)

so a compound that melts congruently melts exactly at T_m,s whatever the liquid model
is. For CA6 and C3A, which melt incongruently, T_m,s is a fitted, metastable melting
point. Tridymite is cristobalite with g_Trd - g_Crs = -ΔH_t (1 - T/T_t), T_t = 1470 °C.

The parameters were CALIBRATED (calibration/fit_Example6_CAS.jl) against the liquidus
projection data/Papers/CAS_Ternary.png: 13 ternary invariant points, 11 points on the
edges, the labelled isotherms and one point inside every primary field (digitised in
data/CAS_liquidus_digitized.txt; the temperatures of the invariant points are from
Osborn & Muan 1960 and Jerebtsov & Mikhailov 2001). ΔH, q and W are *effective*
(fitted) parameters, not calorimetric values. The model reproduces every primary field
of the figure, except two very narrow fields that were left out (Ca3SiO5 between lime
and Ca2SiO4, rankinite Ca3Si2O7 between Ca2SiO4 and pseudowollastonite). With the
parameters below:

                                  model                    Osborn & Muan (1960)
    An + Pwo + Trd (eutectic)     1182 °C, C24 A14 S61     1170 °C, C23 A15 S62
    An + Gh + Pwo  (eutectic)     1279 °C, C37 A21 S42     1265 °C, C38 A20 S42
    An + Mul + Trd                1332 °C, C12 A21 S68     1345 °C, C10 A20 S70
    An + Mul + Crn                1483 °C, C15 A34 S52     1512 °C, C16 A37 S48
    C2S + CA + C3A                1336 °C, C49 A44 S7      1335 °C, C50 A44 S7
    all 13 ternary + 11 edge points: rms 27 K, max 86 K; liquid compositions rms 3 wt%
    every primary field of the figure is reproduced, and the liquid is stable (no
    unmixing) everywhere above the liquidus.

Known shortcoming: the real CaO–SiO2 liquid unmixes near SiO2 (two liquids above
≈ 1690 °C, compare Example 3); this makes the cristobalite liquidus next to the CaO–SiO2
edge very flat (1600–1700 °C in the figure). The model has no liquid immiscibility, so
its cristobalite liquidus there is 100–200 K too low. The lime liquidus at the Lime–C2S
boundary is ≈ 85 K too low.

1 atm only: there is no pressure slider (at high pressure other phases, e.g. garnet and
kyanite, appear, which are not in the model).

How to use
----------
* Right: the ternary (CaO bottom left, Al2O3 bottom right, SiO2 top), wt%.
  - "isothermal section": the stable assemblage at the chosen T on a grid of bulk
    compositions, solved by LP (black dots: boundaries found by bisection).
  - "liquidus projection": primary phase fields, field boundaries, liquidus isotherms
    and invariant points of the model (computed analytically).
  - "reference (figure)": the diagram of the course is drawn underneath, with the
    invariant points used for the calibration (red, literature temperatures).
* Left-click in the triangle (or type "CaO Al2O3 SiO2" in wt% into the box) to choose
  the bulk composition. The T slider sets the temperature.
* Left: the g surface of the liquid (3D, drag to rotate; g per gram over the wt%
  triangle), the solids (dots at their own composition, filled when stable), the
  tangent plane through the stable phases (its heights at the corners are μ_CaO,
  μ_Al2O3, μ_SiO2 per gram) and the bulk composition (star).

Run from a terminal:      julia Example6_ternary_CAS.jl
or from the Julia REPL:   include("Example6_ternary_CAS.jl")
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
const COMPONENTS = ["CaO", "Al2O3", "SiO2"]
const SHORT      = ["C", "A", "S"]                  # cement notation, used in the composition labels
const MOLAR_MASS = [56.077, 101.961, 60.084]        # g/mol

#>>> calibrated parameters (generated by calibration/fit_Example6_CAS.jl)
const Q   = [1.2487, 5.2993, 1.0000]                    # mixing units per oxide (CaO, Al2O3, SiO2)
#  excess terms Π y_i^e_i  (CaO, Al2O3, SiO2 exponents)   W_H [kJ/mol]   W_S [kJ/mol/K]
const WEX = [([2, 1, 0],  -389.289,  -0.115147),
             ([1, 2, 0],  -265.649,  -0.034331),
             ([2, 0, 1],  -984.780,  -0.617487),
             ([1, 0, 2],   276.224,   0.404416),
             ([0, 2, 1],   104.981,   0.086497),
             ([0, 1, 2],    69.689,   0.045273),
             ([1, 1, 1], -1142.876,  -0.221871),
             ([2, 2, 0],   143.181,   0.000650),
             ([2, 0, 2],   386.336,   0.001871),
             ([0, 2, 2],   144.527,   0.002655),
             ([2, 1, 1],  -100.215,  -0.000914),
             ([1, 2, 1],  -304.621,  -0.000142),
             ([1, 1, 2],   175.179,   0.000935),
             ([3, 0, 1], -1197.323,  -0.001348),
             ([1, 0, 3],  -447.052,  -0.002570),
             ([3, 1, 0],  -351.525,   0.000194),
             ([1, 3, 0],   -97.302,  -0.000537),
             ([0, 3, 1],    56.111,  -0.002550),
             ([0, 1, 3],    22.733,  -0.003439),
             ([4, 0, 1],    94.499,  -0.001936),
             ([3, 0, 2],  -464.791,   0.001387),
             ([2, 0, 3],    89.910,   0.000410),
             ([1, 0, 4],   636.543,  -0.003504),
             ([4, 1, 0],   285.617,   0.000824),
             ([1, 4, 0],    68.054,  -0.000912),
             ([3, 1, 1],  -672.263,  -0.000797)]
# effective enthalpies of fusion [kJ per mol of oxides]
const DH  = Dict("Lime" => 358.484, "Crs" => 38.050, "Crn" => 110.778, "Pwo" => 53.789, "C2S" => 28.687, "An" => 64.379, "Gh" => 50.036, "Mul" => 79.999, "CA" => 46.134, "CA2" => 71.524, "CA6" => 109.372, "C3A" => 89.238)
# fitted (metastable) melting points of the incongruently melting compounds [°C]
const TM_FIT = Dict("CA6" => 1893.31, "C3A" => 1538.23)
#<<<

# solids: key, label, oxides (CaO, Al2O3, SiO2) per formula unit, T_m [°C], ΔH [kJ per mol of oxides], colour.
# T_m of CA6 and C3A are fitted metastable melting points (they melt incongruently).
# Tridymite is not listed: it is made from cristobalite (DH_TRD, T_TRD below).
const SOLID_DATA = [
    ("Lime", "Lime CaO",             [1, 0, 0], 2572.0, DH["Lime"], "#c9e3a6"),
    ("Crs",  "Cristobalite SiO2",    [0, 0, 1], 1723.0, DH["Crs"],  "#e3e3e3"),
    ("Crn",  "Corundum Al2O3",       [0, 1, 0], 2054.0, DH["Crn"],  "#b9d3ee"),
    ("Pwo",  "Pseudowollastonite CS",[1, 0, 1], 1544.0, DH["Pwo"],  "#f6c89a"),
    ("C2S",  "Ca2SiO4 C2S",          [2, 0, 1], 2130.0, DH["C2S"],  "#f3e19b"),
    ("An",   "Anorthite CAS2",       [1, 1, 2], 1553.0, DH["An"],   "#e7a9c7"),
    ("Gh",   "Gehlenite C2AS",       [2, 1, 1], 1593.0, DH["Gh"],   "#bfe0d0"),
    ("Mul",  "Mullite A3S2",         [0, 3, 2], 1890.0, DH["Mul"],  "#d8c7ea"),
    ("CA",   "CaAl2O4 CA",           [1, 1, 0], 1604.0, DH["CA"],   "#a8d8e8"),
    ("CA2",  "CaAl4O7 CA2",          [1, 2, 0], 1765.0, DH["CA2"],  "#c6c6f0"),
    ("CA6",  "CaAl12O19 CA6",        [1, 6, 0], TM_FIT["CA6"], DH["CA6"], "#9ec5f4"),
    ("C3A",  "Ca3Al2O6 C3A",         [3, 1, 0], TM_FIT["C3A"], DH["C3A"], "#f0b9a0")]
const DH_TRD    = 1.0                               # kJ/mol, cristobalite -> tridymite
const T_TRD     = 1470.0                            # °C, inversion
const TRD_COLOR = "#f2f2f2"

const LIQ_NAME  = "L"
const NPC_START = 40                                # liquid grid: step 1/n in wt fraction
const NSEC      = 24                                # isothermal section: cells per side
const T_START   = 1400.0                            # [°C]
const X_START   = [0.28, 0.16, 0.56]                # bulk CaO, Al2O3, SiO2 (green dot of the course figure)
const TLIM_C    = (1100.0, 2650.0)
const ISOTHERMS = 1200:100:2500                     # liquidus isotherms of the projection [°C]

# reference: the course figure (PNG, wt%) and the invariant points used for the calibration
const REF_PNG     = joinpath(@__DIR__, "data", "Papers", "CAS_Ternary.png")
const REF_CORNERS = ((152.98, 845.88), (1032.70, 845.88), (594.38, 80.96))   # CaO, Al2O3, SiO2 in the PNG [px]
const REF_DATA    = joinpath(@__DIR__, "data", "CAS_liquidus_digitized.txt")

const LIQ_FIELD_COLOR = "#e6dcef"                   # L alone
const EXTRA_COLORS = ["#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5",
                      "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f", "#a6cee3", "#b2df8a", "#fb9a99", "#cab2d6"]
const SEL_COLOR = :royalblue3
const LIQ_COLOR = :firebrick

# =============================================================================
# THERMODYNAMICS (independent of the plotting)
# =============================================================================
K(T_C) = T_C + 273.15

# ----- moles <-> weight
"weight fractions -> mole fractions"
w2x(w) = (n = w ./ MOLAR_MASS; n ./ sum(n))
"mole fractions -> weight fractions"
x2w(x) = (m = x .* MOLAR_MASS; m ./ sum(m))
"mass [g] of one mole of oxides with mole fractions x"
molar_mass(x) = dot(x, MOLAR_MASS)

"mole fractions x -> unit fractions y (oxide i contributes q_i mixing units)"
unitfrac(x) = (n = x .* Q; n ./ sum(n))

mono(y, e) = y[1]^e[1] * y[2]^e[2] * y[3]^e[3]
dmono(y, e, k) = e[k] == 0 ? 0.0 : e[k] * mono(y, ntuple(i -> e[i] - (i == k), 3))
"∂(n G_ex)/∂n_k for one homogeneous monomial (per unit W)"
exbasis(y, e, k) = dmono(y, e, k) - (sum(e) - 1) * mono(y, e)

"h and s of the liquid per mole of oxides at mole fractions x (g = h - T s, T in K)."
function liquid_hs(x)
    y  = unitfrac(x)
    nu = dot(Q, x)                                  # mol of units per mol of oxides
    h  = nu * sum(wh * mono(y, e) for (e, wh, ws) in WEX)
    s  = nu * (-RGAS * sum(xlogx, y) + sum(ws * mono(y, e) for (e, wh, ws) in WEX))
    return h, s
end

"μ_k = A_k - T B_k of the oxides in the liquid of mole fractions x [kJ/mol]."
function liquid_AB(x)
    y = unitfrac(x)
    A = [Q[k] * sum(wh * exbasis(y, e, k) for (e, wh, ws) in WEX) for k in 1:3]
    B = [Q[k] * (-RGAS * log(max(y[k], 1e-300)) + sum(ws * exbasis(y, e, k) for (e, wh, ws) in WEX)) for k in 1:3]
    return A, B
end
"Liquid chemical potentials [kJ/mol] of CaO, Al2O3, SiO2 at mole fractions x."
liquid_mu(x, T_C) = ((A, B) = liquid_AB(x); A .- K(T_C) .* B)

# ----- solids: g = H - T S, from the liquid of the same composition minus the fusion term
const SOLID_NAMES = vcat([d[1] for d in SOLID_DATA], "Trd")
const SOLID_LABEL = Dict(vcat([d[1] => d[2] for d in SOLID_DATA], "Trd" => "Tridymite SiO2"))
const SOLID_COLOR = Dict(vcat([d[1] => d[6] for d in SOLID_DATA], "Trd" => TRD_COLOR))
const SOLID_X     = vcat([d[3] ./ sum(d[3]) for d in SOLID_DATA], [[0.0, 0.0, 1.0]])   # mole fractions
const SOLID_HS    = let hs = map(SOLID_DATA) do (_, _, ox, tm, dh, _)
        hL, sL = liquid_hs(ox ./ sum(ox))
        (hL - dh, sL - dh / K(tm))
    end
    hc, sc = hs[findfirst(d -> d[1] == "Crs", SOLID_DATA)]
    vcat(hs, [(hc - DH_TRD, sc - DH_TRD / K(T_TRD))])
end
const NS = length(SOLID_NAMES)
const SOLIDS = [Phase(SOLID_NAMES[k], SOLID_X[k], SOLID_HS[k][1], SOLID_HS[k][2], 0.0) for k in 1:NS]

"""
Liquidus temperatures [°C] of all solids for a liquid of WEIGHT composition w (closed form):
solid k is saturated when Σ_i c_k,i μ_i(L) = g_k, and both sides are linear in T.
"""
function tliq_all(w)
    A, B = liquid_AB(w2x(w))
    map(1:NS) do k
        num = dot(SOLID_X[k], A) - SOLID_HS[k][1]
        den = dot(SOLID_X[k], B) - SOLID_HS[k][2]
        den > 0 ? num / den - 273.15 : -Inf
    end
end
"Liquidus temperature [°C] of the liquid of WEIGHT composition w."
tliq(w) = maximum(tliq_all(w))

"""
Liquid pseudocompounds on the triangular grid of step 1/n in WEIGHT fraction (the nodes of
the triangle that is drawn). Each one is an ordinary `Phase` for the LP: composition in
mole fractions, g = g0 - sT in kJ per mole of oxides (T in K; 1 atm, so v = 0).
"""
function liquid_pseudocompounds(n)
    nodes, _, _ = lattice(n)
    map(nodes) do w
        x = w2x(w)
        h, s = liquid_hs(x)
        Phase(LIQ_NAME, x, h, s, 0.0)
    end
end

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

"cartesian position: CaO (0,0), Al2O3 (1,0), SiO2 (1/2, √3/2)"
tocart(w) = Point2f(w[2] + 0.5 * w[3], sqrt(3) / 2 * w[3])
function fromcart(x, y)
    s = y / (sqrt(3) / 2); a = x - 0.5 * s
    w = max.([1 - a - s, a, s], 0.0)
    return w ./ sum(w)
end

"g of the continuous liquid per gram [J/g] at WEIGHT composition w"
gliq(w, T_C) = (x = w2x(w); 1000 * dot(x, liquid_mu(x, T_C)) / molar_mass(x))

"""
Group the stable phases into physical phases (as in Example 5). Two stable liquid
pseudocompounds belong to the same liquid if the continuous g surface between them lies
below their chord (no hump); otherwise they are two liquids (L + L).
"""
function ternary_groups(phases, stable, T_C)
    sol    = [[i] for i in stable if i <= NS]
    groups = [[i] for i in stable if i > NS]
    merged = true
    while merged
        merged = false
        for a in 1:length(groups)-1, b in a+1:length(groups)
            wa = sum(x2w(phases[i].comp) for i in groups[a]) ./ length(groups[a])
            wb = sum(x2w(phases[i].comp) for i in groups[b]) ./ length(groups[b])
            hump = maximum(gliq((1 - t) .* wa .+ t .* wb, T_C) -
                           ((1 - t) * gliq(wa, T_C) + t * gliq(wb, T_C)) for t in 0.1:0.1:0.9)
            if hump <= 1e-6
                append!(groups[a], groups[b]); deleteat!(groups, b); merged = true
                break
            end
        end
    end
    return vcat(sol, groups)
end

"Name, amount (MASS fraction of the system), wt composition w and mole composition x of every stable phase."
function phase_list(phases, r, T_C)
    pl = map(ternary_groups(phases, r.stable, T_C)) do g
        n = sum(r.α[g])                                        # moles of this phase
        x = sum(r.α[i] .* phases[i].comp for i in g) ./ n
        (name = phases[g[1]].name, mass = n * molar_mass(x), x = x)
    end
    mtot = sum(p.mass for p in pl)
    return [(name = p.name, amount = p.mass / mtot, w = x2w(p.x), x = p.x) for p in pl]
end
assemblage_key(pl) = join(sort([p.name for p in pl], by = n -> n == LIQ_NAME ? "~" : n), " + ")
pct(w) = @sprintf("C%.0f A%.0f S%.0f", 100w[1], 100w[2], 100w[3])
phase_label(p) = p.name == LIQ_NAME ? "L(" * pct(p.w) * ")" : p.name

# ----------------------------------------------------------------------------
# analytic liquidus projection
# ----------------------------------------------------------------------------
"Segments where the field `vals` (value per node) crosses `level` (marching triangles)."
function iso_segments(nodes, vals, tris, level)
    segs = Tuple{Vector{Float64},Vector{Float64}}[]
    for t in tris
        cr = Vector{Float64}[]
        for (a, b) in ((t[1], t[2]), (t[2], t[3]), (t[3], t[1]))
            va, vb = vals[a] - level, vals[b] - level
            if isfinite(va) && isfinite(vb) && (va < 0) != (vb < 0)
                s = va / (va - vb); push!(cr, nodes[a] .+ s .* (nodes[b] .- nodes[a]))
            end
        end
        length(cr) == 2 && push!(segs, (cr[1], cr[2]))
    end
    return segs
end
seg_points(segs) = reduce(vcat, [[tocart(a), tocart(b), Point2f(NaN, NaN)] for (a, b) in segs]; init = Point2f[])

"""
Invariant point of the solids `ph` (2 on an edge, 3 inside) near the WEIGHT composition w0:
the liquid on which their liquidus temperatures are equal (Newton inside, bisection on an edge).
Returns (w, T) or `nothing` if there is none nearby.
"""
function invariant_point(ph, w0)
    Tw(w) = (v = max.(w, 1e-9); tliq_all(v ./ sum(v))[ph])
    if length(ph) == 3
        u = [w0[2], w0[3]]
        F(u) = (T = Tw([1 - u[1] - u[2], u[1], u[2]]); [T[1] - T[2], T[1] - T[3]])
        for _ in 1:40
            f = F(u); J = zeros(2, 2)
            for j in 1:2
                v = copy(u); v[j] += 1e-6; J[:, j] = (F(v) .- f) ./ 1e-6
            end
            δ = -(J \ f); s = 1.0
            while s > 1e-3 && minimum([1 - sum(u .+ s .* δ); u .+ s .* δ]) < 0; s /= 2; end
            u .+= s .* δ
            norm(δ) < 1e-10 && break
        end
        w = [1 - u[1] - u[2], u[1], u[2]]
    else
        z = findall(<(1e-9), w0); i, j = setdiff(1:3, z)
        wt(t) = (v = zeros(3); v[i] = 1 - t; v[j] = t; v)
        f(t) = (T = Tw(wt(t)); T[1] - T[2])
        a, b = max(w0[j] - 0.15, 1e-6), min(w0[j] + 0.15, 1 - 1e-6)
        sign(f(a)) == sign(f(b)) && return nothing
        for _ in 1:60
            m = (a + b) / 2; sign(f(m)) == sign(f(a)) ? (a = m) : (b = m)
        end
        w = wt(a)
    end
    (minimum(w) < -1e-9 || norm(w .- w0) > 0.15) && return nothing
    Ts = tliq_all(w)
    # it is an invariant point of the liquidus only if no other solid is higher there
    maximum(Ts) > Ts[ph[1]] + 1.0 && return nothing
    return w, Ts[ph[1]]
end

"Reference points (data/CAS_liquidus_digitized.txt): kind, phase indices, wt composition, T [°C]."
function read_reference(file)
    out = NamedTuple[]
    isfile(file) || return out
    for l in eachline(file)
        (startswith(l, "#") || isempty(strip(l))) && continue
        f = split(l)
        f[1] in ("inv", "bin") || continue
        ph = [findfirst(==(String(s)), SOLID_NAMES) for s in split(f[2], "+")]
        any(isnothing, ph) && continue
        w = parse.(Float64, f[3:5]); w ./= sum(w)
        push!(out, (kind = String(f[1]), ph = Int.(ph), w = w, T = parse(Float64, f[6])))
    end
    return out
end
const REFERENCE = read_reference(REF_DATA)

function liquidus_projection(; N = 90)
    nodes, _, tris = lattice(N)
    Tk   = [tliq_all(w) for w in nodes]                        # Tk[node][solid]
    Tl   = maximum.(Tk)
    prim = [argmax(tliq_all(sum(nodes[collect(t)]) ./ 3)) for t in tris]
    iso  = [(L, seg_points(iso_segments(nodes, Tl, tris, L))) for L in ISOTHERMS]
    # field boundaries: where the two highest liquidus surfaces cross
    present = sort(unique(prim))
    cot = Point2f[]
    for a in eachindex(present), b in a+1:length(present)
        i, j = present[a], present[b]
        for (u, v) in iso_segments(nodes, [T[i] - T[j] for T in Tk], tris, 0.0)
            T = tliq_all((u .+ v) ./ 2)
            max(T[i], T[j]) >= maximum(T) - 0.5 && push!(cot, tocart(u), tocart(v), Point2f(NaN, NaN))
        end
    end
    invs = filter(!isnothing, [invariant_point(r.ph, r.w) for r in REFERENCE])
    return (; nodes, tris, prim, iso, cot, Tl, invs)
end

# =============================================================================
# APPLICATION
# =============================================================================
function build_app()
    # ---------------------------------------------------------------- state
    T      = Ref(T_START)
    X      = Ref(copy(X_START))                       # bulk composition in WEIGHT fractions (as on the triangle)
    npc    = Ref(NPC_START)
    liq    = Ref(liquid_pseudocompounds(NPC_START))
    phases = Ref(vcat(SOLIDS, liq[]))
    lp     = Ref(GibbsLP(phases[], w2x(X[])))
    color_of = Dict{String,String}()
    secplts, selplts = Any[], Any[]
    debounce = Ref(0)

    # w = bulk in weight fractions; the LP is solved in moles (result in moles and kJ/mol), at 1 atm
    equil(w, T_C) = (set_bulk!(lp[], w2x(w)); solve!(lp[], 0.0, K(T_C)))
    function fieldcolor(key)
        key == LIQ_NAME && return LIQ_FIELD_COLOR
        parts = split(key, " + ")
        length(parts) == 2 && parts[2] == LIQ_NAME && return SOLID_COLOR[parts[1]]   # primary field: X + L
        get!(color_of, key, EXTRA_COLORS[mod1(length(color_of) + 1, length(EXTRA_COLORS))])
    end

    # ---------------------------------------------------------------- layout
    fig = Figure(size = (1650, 980), fontsize = 15)
    Label(fig[0, 1:2], "Example 6 – ternary CaO–Al2O3–SiO2 (13 solids and a liquid): " *
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
    limits!(axT, -0.24, 1.24, -0.12, 0.95)
    for name in (:rectanglezoom, :scrollzoom, :dragpan)
        name in keys(interactions(axT)) && deregister_interaction!(axT, name)
    end

    ctrl = right[2, 1] = GridLayout(tellwidth = false)
    sg = SliderGrid(ctrl[1, 1],
        (label = "T [°C]",            range = TLIM_C[1]:1.0:TLIM_C[2], startvalue = T_START, format = "{:.0f}"),
        (label = "L grid (1/n)",      range = 10:2:80,                 startvalue = NPC_START, format = n -> "n = $n"),
        (label = "section cells",     range = 6:2:60,                  startvalue = NSEC,    format = n -> "$n per side"))
    sl_T, sl_npc, sl_sec = sg.sliders
    cb = ctrl[2, 1] = GridLayout(halign = :left)
    Label(cb[1, 1], "bulk wt%  CaO Al2O3 SiO2:", halign = :right)
    tb = Textbox(cb[1, 2], placeholder = "28 16 56", width = 160)
    btn_sec = Button(cb[1, 3], label = "Update section")
    ct = ctrl[3, 1] = GridLayout(halign = :left)
    tg_sec  = Toggle(ct[1, 1], active = true);  Label(ct[1, 2], "isothermal section",  halign = :left)
    tg_proj = Toggle(ct[1, 3], active = true);  Label(ct[1, 4], "liquidus projection", halign = :left)
    tg_ref  = Toggle(ct[1, 5], active = false); Label(ct[1, 6], "reference (figure)",  halign = :left)
    status = Label(right[3, 1], " ", halign = :left, tellwidth = false, color = :gray25)
    cursor = Label(right[4, 1], " ", halign = :left, tellwidth = false, color = :gray25)

    # ------------------------------------------------- static ternary decorations
    corners = [tocart(w) for w in ([1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0])]
    lines!(axT, [corners; corners[1:1]]; color = :black, linewidth = 2)
    for f in 0.1:0.1:0.9                                     # 10 % ticks on the edges
        for (a, b) in ((1, 2), (2, 3), (3, 1))
            w = zeros(3); w[a] = 1 - f; w[b] = f
            scatter!(axT, [tocart(w)]; marker = :circle, markersize = 3, color = :black)
        end
    end
    text!(axT, corners; text = [@sprintf("CaO %.0f", SOLID_DATA[1][4]), @sprintf("Al2O3 %.0f", SOLID_DATA[3][4]),
                                @sprintf("SiO2 %.0f", SOLID_DATA[2][4])],
          align = [(:right, :top), (:left, :top), (:center, :bottom)],
          offset = [(-4, -4), (4, -4), (0, 6)], fontsize = 16, font = :bold)
    # compounds on the triangle (open squares) with their names
    for k in 1:NS
        SOLID_NAMES[k] in ("Lime", "Crn", "Crs", "Trd") && continue
        p = tocart(x2w(SOLID_X[k]))
        scatter!(axT, [p]; marker = :rect, markersize = 8, color = :white, strokecolor = :black, strokewidth = 1.2)
        onedge = minimum(x2w(SOLID_X[k])) < 1e-9
        text!(axT, p; text = SOLID_NAMES[k], fontsize = 12, color = :gray25,
              align = onedge ? (:center, :top) : (:left, :bottom), offset = onedge ? (0, -6) : (5, 3))
    end

    # reference overlay: the course figure underneath + the calibration points
    refpl = Any[]
    try
        img = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
            GLMakie.FileIO.load(REF_PNG)
        end
        (c0, a0, s0) = REF_CORNERS
        sx = 1 / (a0[1] - c0[1]); sy = (sqrt(3) / 2) / (c0[2] - s0[2])
        h, wd = size(img)
        im = image!(axT, (0 - c0[1]) * sx .. (wd - c0[1]) * sx, (c0[2] - h) * sy .. c0[2] * sy, rotr90(img); alpha = 0.55)
        translate!(im, 0, 0, -12)
        push!(refpl, im)
    catch err
        @warn "reference figure not loaded" exception = err
    end
    invref = [r for r in REFERENCE if r.kind == "inv"]
    push!(refpl, scatter!(axT, [tocart(r.w) for r in invref]; marker = :diamond, markersize = 11, color = :red))
    push!(refpl, text!(axT, [tocart(r.w) for r in invref]; text = [@sprintf("%.0f", r.T) for r in invref],
                       fontsize = 11, color = :red, offset = (6, -2), align = (:left, :top)))
    foreach(p -> translate!(p, 0, 0, 6), refpl[2:end])
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
        keyat(w) = assemblage_key(phase_list(phases[], equil(w, T[]), T[]))
        keys_ = [keyat(w) for w in cents]
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
        for key in unique(keys_)
            ids = findall(==(key), keys_)
            length(ids) >= 4 || continue
            c = sum(tocart(cents[i]) for i in ids) / length(ids)
            push!(secplts, text!(axT, c; text = key, align = (:center, :center), fontsize = 13,
                                 font = :bold, color = :gray15, visible = tg_sec.active))
        end
        # analytic liquidus isotherm at T (the L / solid + L boundary must sit on it)
        nl, _, tl = lattice(60)
        iso = seg_points(iso_segments(nl, [tliq(w) for w in nl], tl, T[]))
        push!(secplts, lines!(axT, iso; color = LIQ_COLOR, linewidth = 2.5, linestyle = :dash,
                              visible = tg_sec.active))
        axT.title = @sprintf("isothermal section at %.0f °C, 1 atm  (left-click to pick a bulk composition)", T[])
        status.text[] = @sprintf("Section: %d cells, %d boundary points, ≈ %d LPs (%d phases each) in %.2f s",
                                 length(tris), length(bnd), nlp, length(phases[]), dt)
    end

    "Liquidus projection of the model (computed once: there is no pressure dependence)."
    function draw_projection!()
        pr = liquidus_projection()
        V = Vector{Point2f}(undef, 3length(pr.tris)); F = zeros(Int, length(pr.tris), 3); C = RGBAf[]
        for (m, t) in enumerate(pr.tris)
            for (c, n) in enumerate(t)
                V[3m-3+c] = tocart(pr.nodes[n]); F[m, c] = 3m - 3 + c
                push!(C, RGBAf(parse(RGBf, SOLID_COLOR[SOLID_NAMES[pr.prim[m]]]), 1.0))
            end
        end
        tint = mesh!(axT, V, F; color = C, shading = NoShading,
                     visible = lift((p, s, r) -> p && !s && !r, tg_proj.active, tg_sec.active, tg_ref.active))
        translate!(tint, 0, 0, -11)
        for (L, pts) in pr.iso
            isempty(pts) && continue
            lines!(axT, pts; color = (:black, 0.55), linestyle = :dash, linewidth = 1, visible = tg_proj.active)
            i = findfirst(p -> !isnan(p[1]), pts[div(length(pts), 2):end])
            i === nothing || text!(axT, pts[div(length(pts), 2) - 1 + i]; text = string(L),
                                   fontsize = 11, color = :gray20, visible = tg_proj.active)
        end
        lines!(axT, pr.cot; color = :black, linewidth = 2.5, visible = tg_proj.active)
        # names of the primary fields (at the centre of the triangles of each field)
        for k in unique(pr.prim)
            ids = findall(==(k), pr.prim)
            c = sum(tocart(sum(pr.nodes[collect(pr.tris[i])]) ./ 3) for i in ids) / length(ids)
            text!(axT, c; text = SOLID_NAMES[k], fontsize = 14, font = :bold, color = (:gray10, 0.8),
                  align = (:center, :center),
                  visible = lift((p, s, r) -> p && !s && !r, tg_proj.active, tg_sec.active, tg_ref.active))
        end
        scatter!(axT, [tocart(w) for (w, _) in pr.invs]; markersize = 8, color = :black, visible = tg_proj.active)
        text!(axT, [tocart(w) for (w, _) in pr.invs]; text = [@sprintf("%.0f", t) for (_, t) in pr.invs],
              fontsize = 11, font = :bold, offset = (-5, 3), align = (:right, :bottom), visible = tg_proj.active)
    end

    "g surface, tangent plane and phases for the current bulk composition and T."
    function update_selection!()
        foreach(p -> delete!(axT, p), selplts); empty!(selplts)
        empty!(ax3)
        w, T_C = X[], T[]
        r  = equil(w, T_C)                             # solved in moles: r.g, r.μ, r.G in kJ/mol
        pl = phase_list(phases[], r, T_C)
        iL = [p for p in pl if p.name == LIQ_NAME]
        # μ: a single liquid has a unique tangent plane (continuous model); otherwise the LP duals
        μmol = length(pl) == 1 && length(iL) == 1 && all(iL[1].w .> 1e-9) ? liquid_mu(iL[1].x, T_C) : r.μ   # [kJ/mol]
        # ---- display in weight: g and μ per gram [J/g]
        g  = 1000 .* r.g ./ [molar_mass(p.comp) for p in phases[]]
        μ  = 1000 .* μmol ./ MOLAR_MASS                # heights of the tangent plane at the corners
        G  = 1000 * r.G / molar_mass(w2x(w))

        # liquid surface
        nodes, _, tris = lattice(npc[])
        gl = g[NS+1:end]
        V  = [Point3f(tocart(n)..., gl[i]) for (i, n) in enumerate(nodes)]
        F  = permutedims(reduce(hcat, [collect(t) for t in tris]))
        mesh!(ax3, V, F; color = gl, colormap = :viridis, alpha = 0.8, transparency = true)
        gs = g[1:NS]
        zlo = min(minimum(gl), minimum(μ), minimum(gs)); zhi = max(maximum(gl), maximum(μ), maximum(gs))
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
              text = [rich("μ", subscript(COMPONENTS[k]), @sprintf(" = %.0f", round(μ[k]) + 0.0)) for k in 1:3],
              color = SEL_COLOR, fontsize = 13, offset = (6, 4))
        # solids at their own composition, filled when stable
        st = [any(p -> p.name == SOLID_NAMES[k], pl) for k in 1:NS]
        ws = [x2w(SOLID_X[k]) for k in 1:NS]
        scatter!(ax3, [c3(ws[k], gs[k]) for k in 1:NS]; markersize = 12, strokewidth = 1.5, strokecolor = :black,
                 color = [s ? :black : :white for s in st])
        text!(ax3, [c3(ws[k], gs[k]) for k in 1:NS]; text = SOLID_NAMES, fontsize = 11, offset = (5, 3),
              color = [s ? :black : :gray40 for s in st])
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
        note = length(pl) < 3 && isempty(iL) ? "\n(fewer than 3 phases and no liquid: μ not unique, LP shows one choice)" : ""
        ax3.title = @sprintf("g at %.0f °C, 1 atm, bulk %s (wt%%)", T_C, pct(w))
        info.text[] = "stable: " * amounts *
                      @sprintf("\nG = %.1f J/g  |  μ = (%.0f, %.0f, %.0f) J/g = (%.1f, %.1f, %.1f) kJ/mol  (CaO, Al2O3, SiO2)", G,
                               (round.(μ) .+ 0.0)..., (round.(μmol, digits = 1) .+ 0.0)...) * note

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
            s = y / (sqrt(3) / 2); a = x - 0.5s
            if s >= 0 && a >= 0 && a + s <= 1
                w = fromcart(x, y)
                pl = phase_list(phases[], equil(w, T[]), T[])
                Ts = tliq_all(w)
                cursor.text[] = @sprintf("cursor: %s  →  %s   |  liquidus %.0f °C (%s)", pct(w),
                                         join(phase_label.(pl), " + "), maximum(Ts), SOLID_NAMES[argmax(Ts)])
            end
        end
        return Consume(false)
    end
    on(tb.stored_string) do s
        v = tryparse.(Float64, split(replace(s, "," => " ")))
        if length(v) == 3 && all(!isnothing, v) && all(>=(0), v) && sum(v) > 0
            set_bulk_point!(Float64.(v))
        else
            status.text[] = "Type three numbers: wt% CaO Al2O3 SiO2, e.g. 28 16 56"
        end
    end
    on(sl_T.value) do v
        T[] = v; update_selection!(); schedule_section!()
    end
    on(_ -> rebuild!(), sl_npc.value)
    on(_ -> schedule_section!(), sl_sec.value)
    on(_ -> update_section!(), btn_sec.clicks)

    # ------------------------------------------------------- initialisation
    draw_projection!()
    update_section!()
    update_selection!()
    return fig, (set_bulk = set_bulk_point!, update_section = update_section!,
                 set_T = v -> set_close_to!(sl_T, v), ax3 = ax3,
                 toggles = (section = tg_sec, projection = tg_proj, reference = tg_ref))
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
