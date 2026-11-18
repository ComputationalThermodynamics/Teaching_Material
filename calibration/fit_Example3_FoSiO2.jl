#=
fit_Example3_FoSiO2.jl – calibration of the liquid model of Example3_peritectic_TX.jl

    julia calibration/fit_Example3_FoSiO2.jl

The Fo–SiO2 diagram (Morse 1980, Fig. 10.1, after Bowen & Andersen) is drawn in
WEIGHT per cent SiO2. Its curves were digitised
(data/FoSiO2_liquidus_Morse_digitized.txt). The thermodynamic model works in
MOLES, with the components Mg2SiO4 (Fo) and SiO2 and x = mole fraction of SiO2;
every point is converted (wt% SiO2 -> x) before the fit.

Liquid (per mole of components; the solids Fo and cristobalite are the reference, g = 0):

    g_L = (1-x) ΔH_Fo (1 - T/T_Fo) + x ΔH_SiO2 (1 - T/T_SiO2)
        + R T [x ln x + (1-x) ln(1-x)]
        + x (1-x) Σ_k (a_k - b_k T) (1-2x)^k ,   k = 0, 1, 2      (Redlich–Kister)

Enstatite MgSiO3 = ½ Mg2SiO4 + ½ SiO2 (x = 0.5):  g_En = constant.

Fixed:  T_Fo = 1890 °C, T_SiO2 = 1713 °C (figure), ΔH_Fo = 142 kJ/mol (as in Example 4).
Fitted: ΔH_SiO2, a_k, b_k (6 numbers) and g_En.

Conditions (all linear in T, so each has a closed-form temperature):
    Fo  liquidus   μ_Fo(L) = 0
    Crs liquidus   μ_SiO2(L) = 0
    En  liquidus   ½ μ_Fo(L) + ½ μ_SiO2(L) = g_En
    two liquids    μ_Fo and μ_SiO2 are the same in both liquids
The misfit in temperature is minimised (Levenberg–Marquardt); the labelled
invariant points of the figure are added with a weight. At the end the result is
checked with a true equilibrium calculation (lower convex hull of all phases).
No packages are needed.
=#
using Printf, LinearAlgebra

const R    = 8.31446261815324e-3                   # kJ/mol/K
const M    = (140.69, 60.08)                       # g/mol: Mg2SiO4, SiO2
const S_FO = 100 * M[2] / M[1]                     # wt% SiO2 of Mg2SiO4 (left end of the axis)
K(t) = t + 273.15
const T_FO, T_SI = K(1890.0), K(1713.0)
const H_FO = 142.0                                 # kJ/mol, fixed
const X_EN = 0.5
const NK   = 3                                     # Redlich–Kister terms

"wt% SiO2 of the figure -> mole fraction SiO2 in Mg2SiO4–SiO2"
s2x(s) = (w = (s - S_FO) / (100 - S_FO); (w / M[2]) / (w / M[2] + (1 - w) / M[1]))
x2s(x) = S_FO + (100 - S_FO) * x * M[2] / (x * M[2] + (1 - x) * M[1])

const DATA = [(String(f[1]), s2x(parse(Float64, f[2])), K(parse(Float64, f[3])))
              for f in (split(l) for l in eachline(joinpath(@__DIR__, "..", "data", "FoSiO2_liquidus_Morse_digitized.txt"))
                        if !startswith(l, "#") && !isempty(strip(l)))]
# invariant liquids of the figure: peritectic, eutectic, two liquids + cristobalite
const X_P, X_E, X_M1, X_M2 = s2x(61.25), s2x(64.95), s2x(69.35), s2x(98.65)
const W_INV = 3.0                                  # weight of the invariant points
const W_SOL = 25.0                                 # K per kJ/mol for the two-liquid conditions

# parameters θ = [ΔH_SiO2, a_0, a_1, a_2, b_0, b_1, b_2, g_En]
par(θ) = (HSi = θ[1], a = θ[2:1+NK], b = θ[2+NK:1+2NK], gEn = θ[2+2NK])

"Partial excess functions of the term x(1-x)(1-2x)^k for Fo and SiO2."
function exbasis(x, k)
    u  = 1 - 2x
    G  = x * (1 - x) * u^k
    dG = u^(k + 1) - (k == 0 ? 0.0 : 2k * x * (1 - x) * u^(k - 1))
    return G - x * dG, G + (1 - x) * dG
end
"μ = A - T B for Fo and SiO2 in the liquid of mole fraction x."
function AB(p, x)
    A1, B1 = H_FO, H_FO / T_FO - R * log(1 - x)
    A2, B2 = p.HSi, p.HSi / T_SI - R * log(x)
    for k in 0:NK-1
        e1, e2 = exbasis(x, k)
        A1 += p.a[k+1] * e1; B1 += p.b[k+1] * e1
        A2 += p.a[k+1] * e2; B2 += p.b[k+1] * e2
    end
    return A1, B1, A2, B2
end
mu(p, x, T) = ((A1, B1, A2, B2) = AB(p, x); (A1 - T * B1, A2 - T * B2))
T_Fo(p, x)  = ((A1, B1, _, _) = AB(p, x); A1 / B1)
T_Crs(p, x) = ((_, _, A2, B2) = AB(p, x); A2 / B2)
T_En(p, x)  = ((A1, B1, A2, B2) = AB(p, x); ((1 - X_EN) * A1 + X_EN * A2 - p.gEn) / ((1 - X_EN) * B1 + X_EN * B2))
gL(p, x, T) = ((m1, m2) = mu(p, x, T); (1 - x) * m1 + x * m2)

# the two limbs of the two-liquid field at the same temperature
const PAIRS = let L1 = [(x, T) for (b, x, T) in DATA if b == "L1"], L2 = [(x, T) for (b, x, T) in DATA if b == "L2"]
    [(x1, L2[j][1], T) for (x1, T) in L1 for j in (argmin([abs(T2 - T) for (_, T2) in L2]),) if abs(L2[j][2] - T) < 1]
end

liquidus_misfit(p) = [(b == "Fo" ? T_Fo(p, x) : b == "Pr" ? T_En(p, x) : T_Crs(p, x)) - T for (b, x, T) in DATA if b in ("Fo", "Pr", "Crs")]
function resid(θ)
    p = par(θ); r = liquidus_misfit(p)
    for (x1, x2, T) in PAIRS
        m1, m2 = mu(p, x1, T), mu(p, x2, T)
        push!(r, W_SOL * (m1[1] - m2[1]) / 3, W_SOL * (m1[2] - m2[2]))   # μ_Fo is poorly constrained next to pure SiO2
    end
    push!(r, W_INV * (T_Fo(p, X_P) - K(1557)), W_INV * (T_En(p, X_P) - K(1557)),
             W_INV * (T_En(p, X_E) - K(1543)), W_INV * (T_Crs(p, X_E) - K(1543)),
             W_INV * (T_Crs(p, X_M1) - K(1695)), W_INV * (T_Crs(p, X_M2) - K(1695)))
    return r
end

"Levenberg–Marquardt with a numerical Jacobian."
function lm(resid, θ; iters = 300)
    θ = copy(θ); λ = 1e-2; r = resid(θ); f = sum(abs2, r)
    for _ in 1:iters
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
        ok || break
    end
    return θ
end

# ---------------------------------------------------------------- true equilibrium (check)
"Phases on the lower convex hull at T [K]: list of two-phase gaps (x_a, x_b, name_a, name_b)."
function gaps(p, T; n = 4001)
    xs  = [1e-9; collect(range(0, 1, length = n))[2:end-1]; 1 - 1e-9]
    pts = [(x, gL(p, x, T), "L") for x in xs]
    push!(pts, (0.0, 0.0, "Fo"), (1.0, 0.0, "Crs"), (X_EN, p.gEn, "En"))
    sort!(pts, by = t -> (t[1], t[2]))
    h = eltype(pts)[]
    for q in pts
        !isempty(h) && q[1] == h[end][1] && continue
        while length(h) >= 2
            o, a = h[end-1], h[end]
            (a[1] - o[1]) * (q[2] - o[2]) - (a[2] - o[2]) * (q[1] - o[1]) <= 0 ? pop!(h) : break
        end
        push!(h, q)
    end
    return [(h[i][1], h[i+1][1], h[i][3], h[i+1][3]) for i in 1:length(h)-1
            if h[i+1][1] - h[i][1] > 2.5 / n || h[i][3] != "L" || h[i+1][3] != "L"]
end
has(p, T_C, a, b) = any(g -> g[3] == a && g[4] == b, gaps(p, K(T_C)))
function bisect(pred, lo, hi)
    for _ in 1:50
        m = (lo + hi) / 2
        pred(m) ? (lo = m) : (hi = m)
    end
    return (lo + hi) / 2
end
show_gaps(p, T_C) = join([@sprintf("%s(%.1f) + %s(%.1f)", a, x2s(xa), b, x2s(xb)) for (xa, xb, a, b) in gaps(p, K(T_C))], ",  ")

function report(θ, title)
    p = par(θ); m = liquidus_misfit(p)
    println(title)
    @printf("  ΔH_SiO2 = %.2f kJ/mol,  g_En = %.3f kJ/mol\n", p.HSi, p.gEn)
    @printf("  a = [%s] kJ/mol,  b = [%s] kJ/mol/K\n", join([@sprintf("%.2f", v) for v in p.a], ", "), join([@sprintf("%.5f", v) for v in p.b], ", "))
    @printf("  liquidus misfit (%d points): rms %.1f K, max %.1f K\n", length(m), sqrt(sum(abs2, m) / length(m)), maximum(abs, m))
    Tp = bisect(T -> has(p, T, "Fo", "En"), 1500.0, 1700.0)
    Te = bisect(T -> has(p, T, "En", "Crs"), 1450.0, 1650.0)
    Tm = bisect(T -> !has(p, T, "L", "L"), 1600.0, 1712.0)
    @printf("  peritectic %.1f °C (figure 1557):  %s\n", Tp, show_gaps(p, Tp + 0.2))
    @printf("  eutectic   %.1f °C (figure 1543):  %s\n", Te, show_gaps(p, Te + 0.2))
    @printf("  two liquids + Crs %.1f °C (figure 1695):  %s\n", Tm, show_gaps(p, Tm + 0.2))
    d1, d2 = Float64[], Float64[]
    for (x1, x2, T) in PAIRS, (xa, xb, a, b) in gaps(p, T)
        a == "L" && b == "L" && (push!(d1, x2s(xa) - x2s(x1)); push!(d2, x2s(xb) - x2s(x2)))
    end
    isempty(d1) || @printf("  two-liquid field, model - figure: Mg-rich limb rms %.2f wt%%, SiO2-rich limb rms %.2f wt%% SiO2 (%d temperatures)\n",
                           sqrt(sum(abs2, d1) / length(d1)), sqrt(sum(abs2, d2) / length(d2)), length(d1))
end

θ = lm(resid, [15.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0])
@printf("%d digitised points (wt%% -> mole fraction), %d two-liquid pairs\n", length(DATA), length(PAIRS))
report(θ, "best fit:")
# the values that are typed into Example3_peritectic_TX.jl
θr = [round(θ[1], digits = 2); round.(θ[2:1+NK], digits = 2); round.(θ[2+NK:1+2NK], digits = 5); round(θ[2+2NK], digits = 3)]
report(θr, "rounded values used in the example:")
