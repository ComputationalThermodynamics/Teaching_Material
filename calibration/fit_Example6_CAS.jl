#=
fit_Example6_CAS.jl – calibration of the liquid model of Example6_ternary_CAS.jl

    julia calibration/fit_Example6_CAS.jl

The CaO–Al2O3–SiO2 liquidus projection (data/Papers/CAS_Ternary.png) is drawn in
WEIGHT fractions. Its invariant points, edge points, labelled isotherms and one point
inside every primary field were digitised (data/CAS_liquidus_digitized.txt); the
temperatures of the invariant points are literature values (see that file). The model
works in MOLES of the oxides CaO, Al2O3, SiO2; every point is converted before the fit.

Liquid (per mole of oxides; the pure oxide liquids are the reference, g = 0). Oxide i
contributes q_i mixing units, unit fractions y_i = q_i x_i / Σ q_j x_j, ν = Σ q_j x_j:

    g_L = ν [ R T Σ y_i ln y_i + Σ_m W_m Π_i y_i^e_m,i ]

Solids (stoichiometric compounds, composition c = oxide mole fractions):

    g_s(T) = g_L(c, T) - ΔH_s (1 - T/T_s)

i.e. the solid is the liquid of the same composition minus a fusion term. For a
compound that melts congruently T_s is its melting point (fixed); for CA6 and C3A,
which melt incongruently, T_s is a fitted (metastable) melting point. Tridymite is
cristobalite with g_Trd - g_Crs = -ΔH_t (1 - T/1743.15 K) (inversion at 1470 °C).

Every condition is linear in T, so the liquidus temperature of a solid s at liquid
composition x has a closed form:  T_s(x) = (Σ c_k A_k(x) - H_s) / (Σ c_k B_k(x) - S_s),
with μ_k(L) = A_k - T B_k. At an invariant point the liquidus temperatures of its
phases must equal the given temperature and all other solids must be lower. The misfit
in temperature is minimised (Levenberg–Marquardt). A penalty keeps the liquid stable
(convex g, no unmixing) on a grid of compositions, at the liquidus and at 2800 °C: with
a stable liquid no liquidus can rise above the melting point of its solid. Only the
seven third-degree excess terms get an entropy part W_S.

The result of every run is written to fit_Example6_CAS_start.jl and is the start of the
next run (delete that file to start afresh). The argument is the number of iterations:

    julia calibration/fit_Example6_CAS.jl 800

At the end the script prints the parameter block for the USER INPUT section of
Example6_ternary_CAS.jl (between the lines #>>> and #<<<). No packages are needed.
=#
using Printf, LinearAlgebra, Random
Random.seed!(1)

const R  = 8.31446261815324e-3                       # kJ/mol/K
const M  = [56.077, 101.961, 60.084]                 # CaO, Al2O3, SiO2 [g/mol]
K(t) = t + 273.15

# ------------------------------------------------------------------ solids
#              name    oxides (CaO, Al2O3, SiO2)   T_m [°C] (nothing = fitted)
const SOLIDS = [("Lime", [1, 0, 0], 2572.0),
                ("Crs",  [0, 0, 1], 1723.0),
                ("Trd",  [0, 0, 1], nothing),        # from Crs (inversion), not fitted
                ("Crn",  [0, 1, 0], 2054.0),
                ("Pwo",  [1, 0, 1], 1544.0),
                ("C2S",  [2, 0, 1], 2130.0),
                ("An",   [1, 1, 2], 1553.0),
                ("Gh",   [2, 1, 1], 1593.0),
                ("Mul",  [0, 3, 2], 1890.0),
                ("CA",   [1, 1, 0], 1604.0),
                ("CA2",  [1, 2, 0], 1765.0),
                ("CA6",  [1, 6, 0], nothing),
                ("C3A",  [3, 1, 0], nothing)]
const NS   = length(SOLIDS)
const NAME = first.(SOLIDS)
const C    = [s[2] ./ sum(s[2]) for s in SOLIDS]     # mole fractions of the oxides
const ITRD = findfirst(==("Trd"), NAME)
const ICRS = findfirst(==("Crs"), NAME)
const FIT_TM = [i for i in 1:NS if SOLIDS[i][3] === nothing && i != ITRD]   # CA6, C3A
const FIT_H  = [i for i in 1:NS if i != ITRD]
const DH_TRD = 1.0                                   # kJ/mol, Crs -> Trd
const T_TRD  = K(1470.0)

# ------------------------------------------------------------------ liquid excess terms
const EXPS = [[2, 1, 0], [1, 2, 0], [2, 0, 1], [1, 0, 2], [0, 2, 1], [0, 1, 2], [1, 1, 1],
              [2, 2, 0], [2, 0, 2], [0, 2, 2], [2, 1, 1], [1, 2, 1], [1, 1, 2],
              [3, 0, 1], [1, 0, 3], [3, 1, 0], [1, 3, 0], [0, 3, 1], [0, 1, 3],
              [4, 0, 1], [3, 0, 2], [2, 0, 3], [1, 0, 4], [4, 1, 0], [1, 4, 0], [3, 1, 1]]
const NW = length(EXPS)

mono(y, e) = y[1]^e[1] * y[2]^e[2] * y[3]^e[3]
dmono(y, e, k) = e[k] == 0 ? 0.0 : e[k] * mono(y, ntuple(i -> e[i] - (i == k), 3))
exbasis(y, e, k) = dmono(y, e, k) - (sum(e) - 1) * mono(y, e)
xlogx(v) = v > 0 ? v * log(v) : 0.0

# parameters θ = [log q_CaO, log q_Al2O3, W_H (NW), W_S (NW), log ΔH (FIT_H), T_m (FIT_TM, °C)]
const IH = 2 + 2NW                                   # offset of the log ΔH block
function par(θ)
    q  = [exp(θ[1]), exp(θ[2]), 1.0]
    W  = θ[3:2+NW]; WS = θ[3+NW:2+2NW]
    dh = zeros(NS); tm = zeros(NS)
    for (j, i) in enumerate(FIT_H); dh[i] = exp(θ[IH+j]); end
    for i in 1:NS; SOLIDS[i][3] === nothing || (tm[i] = K(SOLIDS[i][3])); end
    for (j, i) in enumerate(FIT_TM); tm[i] = K(θ[IH+length(FIT_H)+j]); end
    # h, s of every solid (g = h - T s)
    H = zeros(NS); S = zeros(NS)
    for i in 1:NS
        i == ITRD && continue
        hL, sL = hs_liquid(q, W, WS, C[i])
        H[i] = hL - dh[i]; S[i] = sL - dh[i] / tm[i]
    end
    H[ITRD] = H[ICRS] - DH_TRD; S[ITRD] = S[ICRS] - DH_TRD / T_TRD
    return (q = q, W = W, WS = WS, dh = dh, tm = tm, H = H, S = S)
end

"h and s of the liquid per mole of oxides at mole fractions x (g = h - T s)."
function hs_liquid(q, W, WS, x)
    n = q .* x; nu = sum(n); y = n ./ nu
    h = nu * sum(W[m] * mono(y, EXPS[m]) for m in 1:NW)
    s = nu * (-R * sum(xlogx, y) + sum(WS[m] * mono(y, EXPS[m]) for m in 1:NW))
    return h, s
end

"μ_k = A_k - T B_k of the three oxides in the liquid of mole fractions x."
function AB(p, x)
    n = p.q .* x; y = n ./ sum(n)
    A = [p.q[k] * sum(p.W[m] * exbasis(y, EXPS[m], k) for m in 1:NW) for k in 1:3]
    B = [p.q[k] * (-R * log(max(y[k], 1e-300)) + sum(p.WS[m] * exbasis(y, EXPS[m], k) for m in 1:NW)) for k in 1:3]
    return A, B
end

"Liquidus temperature [°C] of every solid at liquid mole fractions x."
function tliq_all(p, x)
    A, B = AB(p, x)
    map(1:NS) do i
        num = dot(C[i], A) - p.H[i]; den = dot(C[i], B) - p.S[i]
        den > 0 ? num / den - 273.15 : -1e4
    end
end

w2x(w) = (n = w ./ M; n ./ sum(n))
x2w(x) = (m = x .* M; m ./ sum(m))

# ------------------------------------------------------------------ data
struct Datum
    kind  :: String
    ph    :: Vector{Int}
    x     :: Vector{Float64}     # mole fractions
    w     :: Vector{Float64}     # wt fractions
    T     :: Float64
    wt    :: Float64
end
const DATA = let d = Datum[]
    for l in eachline(joinpath(@__DIR__, "..", "data", "CAS_liquidus_digitized.txt"))
        (startswith(l, "#") || isempty(strip(l))) && continue
        f = split(l)
        w = parse.(Float64, f[3:5]) ./ 100; w = w ./ sum(w)
        ph = [findfirst(==(String(s)), NAME) for s in split(f[2], "+")]
        push!(d, Datum(String(f[1]), ph, w2x(w), w, f[6] == "-" ? NaN : parse(Float64, f[6]), parse(Float64, f[7])))
    end
    d
end

const START = joinpath(@__DIR__, "fit_Example6_CAS_start.jl")   # last result, used as start of the next run (delete to start afresh)
const GAP_INV   = 5.0      # K: other solids must be this much below an invariant / edge / isotherm point
const GAP_FIELD = 15.0     # K: the primary phase must be this much above all others inside its field
const W_HINGE   = 1.0

"smooth hinge max(0, v) (softplus, width 2 K) – keeps the misfit differentiable"
hinge(v) = (s = 2.0; v > 20s ? v : s * log1p(exp(v / s)))
const W_STAB = 5.0         # weight of the liquid-stability penalty
const D2_MIN = 3.0         # kJ/mol: required curvature margin of the liquid
const T_TOP  = 2800.0      # °C, above every melting point: a stable liquid cannot have a liquidus above T_m

function resid(θ)
    p = par(θ); r = Float64[]
    for d in DATA
        Ts = tliq_all(p, d.x)
        if d.kind == "field"
            Tp = Ts[d.ph[1]]
            for i in 1:NS
                i in d.ph || (i == ITRD && d.ph[1] == ICRS) || (i == ICRS && d.ph[1] == ITRD) ||
                    push!(r, W_HINGE * hinge(Ts[i] - Tp + GAP_FIELD))
            end
        else
            for i in d.ph
                push!(r, d.wt * (Ts[i] - d.T))
            end
            for i in 1:NS
                i in d.ph && continue
                ((i == ITRD && ICRS in d.ph) || (i == ICRS && ITRD in d.ph)) && continue
                push!(r, W_HINGE * hinge(Ts[i] - d.T + GAP_INV))
            end
        end
    end
    # the liquid must not unmix above the liquidus
    append!(r, W_STAB .* hinge.(D2_MIN .- stab_d2(p)))
    # weak priors: keep the fusion enthalpies and the excess terms reasonable
    for (j, i) in enumerate(FIT_H); push!(r, 0.5 * (θ[IH+j] - log(40.0))); end
    # W_S only for the seven third-degree terms (as in Example 5); the higher terms are kept enthalpic
    for m in 1:NW; push!(r, 0.004 * θ[2+m]); push!(r, (m <= 7 ? 50.0 : 2000.0) * θ[2+NW+m]); end
    return r
end

"Levenberg–Marquardt with a numerical Jacobian."
function lm(resid, θ; iters = 400, verbose = false)
    θ = copy(θ); λ = 1e-2; r = resid(θ); f = sum(abs2, r)
    for it in 1:iters
        J = zeros(length(r), length(θ))
        for j in eachindex(θ)
            h = 1e-6 * max(1, abs(θ[j])); t = copy(θ); t[j] += h
            J[:, j] = (resid(t) .- r) ./ h
        end
        A, g = J' * J, J' * r
        ok = false
        for _ in 1:30
            δ = -(A + λ * Diagonal(diag(A) .+ 1e-9)) \ g
            rn = resid(θ .+ δ); fn = sum(abs2, rn)
            if isfinite(fn) && fn < f
                θ .+= δ; r = rn; f = fn; λ = max(λ / 3, 1e-9); ok = true; break
            end
            λ *= 4
        end
        verbose && it % 50 == 0 && @printf("  it %4d  f = %.4g\n", it, f)
        if !ok                                    # stuck: small random kick, keep it only if it helps later
            λ = 1e-2
            t = θ .+ 0.01 .* randn(length(θ)) .* max.(abs.(θ), 0.1)
            rt = resid(t); ft = sum(abs2, rt)
            ft < 1.5f && (θ = t; r = rt; f = ft)
        end
    end
    return θ
end

# ------------------------------------------------------------------ checks
"Model invariant point of the phases ph (2 or 3) near the wt composition w0 (Newton, in wt fractions)."
function invariant(p, ph, w0)
    Tw(w) = tliq_all(p, w2x(max.(w, 1e-9) ./ sum(max.(w, 1e-9))))[ph]
    if length(ph) == 3
        u = [w0[2], w0[3]]
        for _ in 1:40
            F(u) = (T = Tw([1 - u[1] - u[2], u[1], u[2]]); [T[1] - T[2], T[1] - T[3]])
            f = F(u); J = zeros(2, 2)
            for j in 1:2
                h = 1e-6; v = copy(u); v[j] += h; J[:, j] = (F(v) .- f) ./ h
            end
            δ = -(J \ f); s = 1.0
            while s > 1e-3 && (minimum([1 - sum(u .+ s .* δ); u .+ s .* δ]) < 0); s /= 2; end
            u .+= s .* δ
            norm(δ) < 1e-10 && break
        end
        w = [1 - u[1] - u[2], u[1], u[2]]
    else                                  # on an edge: bisection along the edge
        z = findall(<(1e-9), w0); i, j = setdiff(1:3, z)
        wt(t) = (v = zeros(3); v[i] = 1 - t; v[j] = t; v)
        f(t) = (T = Tw(wt(t)); T[1] - T[2])
        a, b = max(w0[j] - 0.15, 1e-6), min(w0[j] + 0.15, 1 - 1e-6)
        if sign(f(a)) == sign(f(b))
            return w0, NaN
        end
        for _ in 1:60
            m = (a + b) / 2; sign(f(m)) == sign(f(a)) ? (a = m) : (b = m)
        end
        w = wt(a)
    end
    Ts = tliq_all(p, w2x(w))
    return w, Ts[ph[1]]
end

"Primary phase map on a lattice of step 1/n (wt) and the list of field points with the wrong phase."
function field_check(p)
    bad = String[]
    for d in DATA
        d.kind in ("field", "iso") || continue
        Ts = tliq_all(p, d.x); k = argmax(Ts)
        ok = k == d.ph[1] || (NAME[k] in ("Crs", "Trd") && NAME[d.ph[1]] in ("Crs", "Trd"))
        ok || push!(bad, @sprintf("%s point %s: model primary phase %s", d.kind, NAME[d.ph[1]], NAME[k]))
    end
    return bad
end

# The CaO–SiO2 liquid really unmixes near SiO2 (two liquids above ≈ 1690 °C, Greig 1927); the flat
# 1600/1700 °C isotherms of the cristobalite field reflect this. The model has no such two-liquid
# field (it would need two-liquid data, as in Example 3): the liquid must be stable everywhere above
# the liquidus, and the cristobalite isotherms next to the CaO–SiO2 edge get a small weight. Setting
# a region here allows a miscibility gap in it (tested: the gap then grows far too large).
in_gap_region(w) = false

"Second differences of g_L (per mole, along the 3 lattice directions in mole fractions): (d2, w, T °C)."
function stab_list(p; n = 30)
    out = Tuple{Float64,Vector{Float64},Float64,Bool}[]
    h = 0.5 / n
    for i in 0:n, j in 0:n-i                              # the binary edges are included
        x = [1 - (i + j) / n, i / n, j / n]
        # d²g is linear in T: checking the liquidus and the top of the T range covers everything in between
        for T in (K(clamp(maximum(tliq_all(p, x)), 1000.0, T_TOP)), K(T_TOP))
            gl(z) = ((hh, ss) = hs_liquid(p.q, p.W, p.WS, z); hh - T * ss)
            for d in ([1, -1, 0], [1, 0, -1], [0, 1, -1])
                xa, xb = x .+ h .* d, x .- h .* d
                minimum([xa; xb]) < -1e-12 && continue
                push!(out, ((gl(xa) + gl(xb) - 2gl(x)) / h^2, x2w(x), T - 273.15, T == K(T_TOP)))
            end
        end
    end
    return out
end
# inside the region the gap must close at high T (stable at T_TOP), as a real miscibility gap does
stab_d2(p; n = 40) = [d for (d, w, T, top) in stab_list(p; n) if !in_gap_region(w) || top]

"Liquid stability: worst second difference outside and inside the region where unmixing is allowed."
function stability(p; n = 40)
    l = stab_list(p; n)
    out = [t for t in l if !in_gap_region(t[2]) || t[4]]; ins = [t for t in l if in_gap_region(t[2]) && !t[4]]
    return out[argmin(first.(out))], isempty(ins) ? (Inf, [NaN, NaN, NaN], NaN, false) : ins[argmin(first.(ins))]
end

pct(w) = @sprintf("%4.1f %4.1f %4.1f", 100w[1], 100w[2], 100w[3])

function report(θ, title)
    p = par(θ)
    println(title)
    @printf("  q = (%.3f, %.3f, %.3f)\n", p.q...)
    println("  W_H [kJ/mol]: ", join([@sprintf("%s %.2f", join(EXPS[m]), p.W[m]) for m in 1:NW], ",  "))
    println("  W_S [J/mol/K]: ", join([@sprintf("%s %.1f", join(EXPS[m]), 1000p.WS[m]) for m in 1:NW], ",  "))
    println("  ΔH [kJ per mol of oxides]: ", join([@sprintf("%s %.2f", NAME[i], p.dh[i]) for i in FIT_H], ",  "))
    println("  fitted T_m [°C]: ", join([@sprintf("%s %.1f", NAME[i], p.tm[i] - 273.15) for i in FIT_TM], ",  "))
    dT = Float64[]; dw = Float64[]
    println("  point                    figure (CaO Al2O3 SiO2, T)        model")
    for d in DATA
        d.kind in ("inv", "bin") || continue
        w, T = invariant(p, d.ph, d.w)
        push!(dT, T - d.T); push!(dw, 100 * norm(w .- d.w))
        @printf("  %-4s %-16s  %s  %6.0f     %s  %6.0f   (ΔT %+5.0f K, Δw %.1f wt%%)\n", d.kind,
                join(NAME[d.ph], "+"), pct(d.w), d.T, pct(w), T, T - d.T, 100 * norm(w .- d.w))
    end
    fin = isfinite.(dT)
    @printf("  invariant + edge points: T rms %.1f K, max %.1f K; composition rms %.2f wt%%, max %.2f wt%%\n",
            sqrt(sum(abs2, dT[fin]) / count(fin)), maximum(abs, dT[fin]), sqrt(sum(abs2, dw[fin]) / count(fin)), maximum(dw[fin]))
    for d in DATA
        d.kind == "iso" || continue
        Ts = tliq_all(p, d.x)
        @printf("  isotherm %-5s %s  %5.0f  model %5.0f (%s)\n", NAME[d.ph[1]], pct(d.w), d.T, maximum(Ts), NAME[argmax(Ts)])
    end
    bad = field_check(p)
    println(isempty(bad) ? "  all field and isotherm points lie in the right primary field" : "  WRONG FIELD: " * join(bad, "; "))
    so, si = stability(p)
    @printf("  liquid stability: min d²g/dx² = %.1f kJ/mol at %s, %.0f °C (%s)\n", so[1], pct(so[2]), so[3],
            so[1] > 0 ? "stable (corner: closed at the top of the T range)" : "UNMIXING")
    @printf("  high-SiO2 corner at the liquidus:  min d²g/dx² = %.1f kJ/mol at %s, %.0f °C (%s)\n", si[1], pct(si[2]), si[3],
            si[1] > 0 ? "no unmixing" : "two liquids possible")
end

θ0 = [0.0, 0.0,
      -40.0, -40.0, -150.0, -150.0, -40.0, -40.0, -100.0, zeros(NW - 7)..., zeros(NW)...,
      fill(log(40.0), length(FIT_H))..., 1800.0, 1450.0]
θ = isfile(START) ? include(START) : θ0
length(θ) == length(θ0) || (θ = θ0)
θ = lm(resid, θ; iters = isempty(ARGS) ? 400 : parse(Int, ARGS[1]), verbose = true)
open(START, "w") do io; println(io, repr(θ)); end
report(θ, "best fit:")

# ---------------------------------------------------------------- values for Example6_ternary_CAS.jl
"θ rebuilt from the rounded values that are printed for the example."
function rounded(θ)
    p = par(θ); t = copy(θ)
    t[1:2] = log.(round.(p.q[1:2], digits = 4))
    t[3:2+NW] = round.(p.W, digits = 3); t[3+NW:2+2NW] = round.(p.WS, digits = 6)
    for (j, i) in enumerate(FIT_H); t[IH+j] = log(round(p.dh[i], digits = 3)); end
    for (j, i) in enumerate(FIT_TM); t[IH+length(FIT_H)+j] = round(p.tm[i] - 273.15, digits = 2); end
    return t
end
θr = rounded(θ)
report(θr, "\nrounded values used in the example:")
p = par(θr)
println("\n# ---- paste into the USER INPUT block of Example6_ternary_CAS.jl ----")
@printf("const Q   = [%.4f, %.4f, %.4f]                    # mixing units per oxide (CaO, Al2O3, SiO2)\n", p.q...)
println("#  excess terms Π y_i^e_i  (CaO, Al2O3, SiO2 exponents)   W_H [kJ/mol]   W_S [kJ/mol/K]")
println("const WEX = [", join([@sprintf("(%s, %9.3f, %10.6f)", repr(EXPS[m]), p.W[m], p.WS[m]) for m in 1:NW], ",\n             "), "]")
println("# effective enthalpies of fusion [kJ per mol of oxides]")
println("const DH  = Dict(", join([@sprintf("\"%s\" => %.3f", NAME[i], p.dh[i]) for i in FIT_H], ", "), ")")
println("# fitted (metastable) melting points of the incongruently melting compounds [°C]")
println("const TM_FIT = Dict(", join([@sprintf("\"%s\" => %.2f", NAME[i], p.tm[i] - 273.15) for i in FIT_TM], ", "), ")")
