#=
fit_Example2_DiAn.jl – calibration of the liquid model of Example2_eutectic_TX.jl

    julia calibration/fit_Example2_DiAn.jl

The Di–An liquidus was read off the diagram in WEIGHT per cent
(data/DiAn_liquidus_digitized.txt). The thermodynamic model works with MOLE
fractions, so every point is first converted (wt% An -> x_An) and the model is
then fitted in mole fractions:

    Di liquidus:  ΔH_Di (1 - T/T_m,Di) + R T ln(1-x) + W x²     = 0
    An liquidus:  ΔH_An (1 - T/T_m,An) + R T ln x     + W (1-x)² = 0

Both are linear in T, so T_liquidus(x) has a closed form. ΔH_Di, ΔH_An and W are
found by least squares on the temperature misfit (Levenberg–Marquardt); the
melting points are those of the diagram. The eutectic of the diagram
(42 wt% An, 1274 °C) is added as an extra, weighted point on both curves.
No packages are needed.
=#
using Printf, LinearAlgebra

const R  = 8.31446261815324e-3                 # kJ/mol/K
const M  = (216.55, 278.21)                    # g/mol: CaMgSi2O6, CaAl2Si2O8
const TM = (1391.0 + 273.15, 1553.0 + 273.15)  # melting points of the diagram [K]
const EUT = (0.42, 1274.0 + 273.15)            # eutectic of the diagram: wt fraction An, T [K]
const W_EUT = 5.0                              # weight of the eutectic

w2x(w) = (w / M[2]) / (w / M[2] + (1 - w) / M[1])          # wt fraction An -> mole fraction An
x2w(x) = x * M[2] / (x * M[2] + (1 - x) * M[1])

const DATA = [(String(f[1]), w2x(parse(Float64, f[2]) / 100), parse(Float64, f[3]) + 273.15)
              for f in (split(l) for l in eachline(joinpath(@__DIR__, "..", "data", "DiAn_liquidus_digitized.txt"))
                        if !startswith(l, "#") && !isempty(strip(l)))]

TDi(θ, x) = (θ[1] + θ[3] * x^2) / (θ[1] / TM[1] - R * log(1 - x))
TAn(θ, x) = (θ[2] + θ[3] * (1 - x)^2) / (θ[2] / TM[2] - R * log(x))
misfit(θ) = [(b == "Di" ? TDi(θ, x) : TAn(θ, x)) - T for (b, x, T) in DATA]
resid(θ)  = vcat(misfit(θ), W_EUT * (TDi(θ, w2x(EUT[1])) - EUT[2]), W_EUT * (TAn(θ, w2x(EUT[1])) - EUT[2]))

"Levenberg–Marquardt with a numerical Jacobian."
function lm(resid, θ; iters = 200)
    θ = copy(θ); λ = 1e-3; r = resid(θ); f = sum(abs2, r)
    for _ in 1:iters
        J = reduce(hcat, [(resid(θ .+ (h = 1e-6 * max(1, abs(θ[j])); [i == j ? h : 0.0 for i in eachindex(θ)])) .- r) ./
                          (1e-6 * max(1, abs(θ[j]))) for j in eachindex(θ)])
        A, g = J' * J, J' * r
        ok = false
        for _ in 1:30
            δ = -(A + λ * Diagonal(diag(A))) \ g
            rn = resid(θ .+ δ); fn = sum(abs2, rn)
            if fn < f
                θ .+= δ; r = rn; f = fn; λ /= 3; ok = true; break
            end
            λ *= 4
        end
        ok || break
    end
    return θ
end

function eutectic(θ)
    a, b = 1e-4, 1 - 1e-4
    for _ in 1:80
        m = (a + b) / 2
        TDi(θ, m) > TAn(θ, m) ? (a = m) : (b = m)
    end
    return a, TDi(θ, a)
end

θ = lm(resid, [130.0, 130.0, -10.0])
xe, Te = eutectic(θ); m = misfit(θ)
@printf("%d points (wt%% -> mole fraction)\n", length(DATA))
@printf("ΔH_Di = %.1f kJ/mol,  ΔH_An = %.1f kJ/mol,  W = %.1f kJ/mol\n", θ...)
@printf("eutectic: %.1f mol%% An = %.1f wt%% An, %.1f °C\n", 100xe, 100x2w(xe), Te - 273.15)
@printf("liquidus misfit: rms %.2f K, max %.1f K\n", sqrt(sum(abs2, m) / length(m)), maximum(abs, m))
θr = round.(θ, digits = 1); xe, Te = eutectic(θr); m = misfit(θr)
@printf("with the rounded values used in the example (%.1f, %.1f, %.1f): eutectic %.1f wt%% An, %.1f °C; rms %.2f K\n",
        θr..., 100x2w(xe), Te - 273.15, sqrt(sum(abs2, m) / length(m)))
