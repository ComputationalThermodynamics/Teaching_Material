#=
gibbs_lp.jl  –  core of the interactive phase-equilibrium examples

Gibbs energy minimisation as a linear programme (same idea as LinProg1.ipynb):

    min   Σ_i g_i(P,T) α_i                (total Gibbs energy of the system)
    s.t.  Σ_i A[j,i] α_i = Xsys[j]        (mass balance for every oxide j)
          α_i ≥ 0                         (phase amounts cannot be negative)

The columns of A are normalised phase compositions (they sum to 1), so
Σ α_i = 1 follows from the mass balance and does not need its own constraint.
The Lagrange multipliers (duals) of the mass-balance constraints are the
chemical potentials μ_j of the oxides: the tangent line g = Σ_j μ_j X_j touches
every stable phase and lies below every metastable one.

Nothing in this file depends on the plotting package, so it can be used and
tested on its own.
=#

using JuMP, GLPK, LinearAlgebra, Printf

# ----------------------------------------------------------------------------
# Phases of fixed composition:  g(P,T) = g0 - s*T + v*P
# ----------------------------------------------------------------------------
struct Phase
    name :: String
    comp :: Vector{Float64}   # amount of each oxide (any scale, normalised below)
    g0   :: Float64           # Gibbs energy at T = 0, P = 0   [kJ/mol]
    s    :: Float64           # entropy                         [kJ/mol/K]
    v    :: Float64           # volume                          [kJ/kbar/mol = J/bar]
end

gibbs(p::Phase, P, T) = p.g0 - p.s * T + p.v * P
gibbs(phases::Vector{Phase}, P, T) = [gibbs(p, P, T) for p in phases]

"Constraint matrix: one column per phase, normalised so that each column sums to 1."
comp_matrix(phases::Vector{Phase}) = reduce(hcat, [p.comp ./ sum(p.comp) for p in phases])

"Composition of each phase on the binary X axis (fraction of the 2nd oxide)."
xcomp(phases::Vector{Phase}) = [p.comp[2] / sum(p.comp) for p in phases]

# ----------------------------------------------------------------------------
# Re-usable LP. The model is built once; only the objective (g changes with
# P,T) or the right-hand side (bulk composition) is updated between solves,
# which makes scanning a grid fast.
# ----------------------------------------------------------------------------
struct GibbsLP
    phases :: Vector{Phase}
    model  :: Model
    a      :: Vector{VariableRef}
    mass   :: Vector{ConstraintRef}
end

function GibbsLP(phases::Vector{Phase}, Xsys::AbstractVector)
    A      = comp_matrix(phases)
    nc, np = size(A)
    m      = Model(GLPK.Optimizer)
    set_silent(m)
    @variable(m, a[1:np] >= 0)
    @constraint(m, mass[j in 1:nc], sum(A[j, i] * a[i] for i in 1:np) == Xsys[j])
    @objective(m, Min, sum(a))          # coefficients are replaced by g in solve!
    return GibbsLP(phases, m, collect(a), collect(mass))
end

"Change the bulk composition of the system."
function set_bulk!(lp::GibbsLP, Xsys::AbstractVector)
    for j in eachindex(lp.mass)
        set_normalized_rhs(lp.mass[j], Xsys[j])
    end
    return lp
end

"""
    solve!(lp, P, T; tol=1e-6)

Minimise G at (P,T). Returns a NamedTuple with
  α       phase amounts (tiny values cleaned to 0)
  g       Gibbs energy of every phase
  G       Gibbs energy of the system (= α'g)
  μ       chemical potentials of the oxides (duals of the mass balance)
  stable  indices of the stable phases (sorted)
"""
function solve!(lp::GibbsLP, P, T; tol = 1e-6)
    g = gibbs(lp.phases, P, T)
    set_objective_coefficient(lp.model, lp.a, g)     # only the costs change with P,T
    optimize!(lp.model)
    is_solved_and_feasible(lp.model; dual = true) ||
        error("LP failed at P=$P, T=$T: $(termination_status(lp.model))")
    α = value.(lp.a)
    α[abs.(α) .< tol] .= 0.0
    μ = dual.(lp.mass)
    return (α = α, g = g, G = dot(α, g), μ = μ, stable = findall(>(0), α))
end

"Readable name of an assemblage, e.g. \"A + AB\"."
assemblage_name(phases, stable) = join((phases[i].name for i in stable), " + ")

# ----------------------------------------------------------------------------
# Solution phases as pseudocompounds (Example 2)
# ----------------------------------------------------------------------------
const RGAS = 8.31446261815324e-3          # gas constant [kJ/mol/K]

"x ln x with the limit 0 ln 0 = 0."
xlogx(x) = x > 0 ? x * log(x) : 0.0

"""
    solution_pseudocompounds(name, e1, e2, W, n)

Discretise a binary solution between the end-members e1 (X = 0) and e2 (X = 1)
into n phases of fixed composition ("pseudocompounds", X = 0, 1/(n-1), …, 1):

    g(X) = (1-X) g_e1 + X g_e2 + RT [X ln X + (1-X) ln(1-X)] + W X (1-X)

(ideal mixing + symmetric regular-solution / Margules term W [kJ/mol]).
Optionally W = W - T WS + P WV (keywords WS [kJ/mol/K], WV [kJ/kbar/mol]), and
`sites` mixing sites per formula unit (ideal entropy -sites·R[X lnX + (1-X)ln(1-X)],
e.g. 2 for Mg–Fe mixing in olivine (Mg,Fe)2SiO4).
This is still linear in P and T, so every pseudocompound is an ordinary
`Phase`: the configurational entropy simply goes into s. All pseudocompounds
get the same name; the LP then decides which of them are stable.
"""
function solution_pseudocompounds(name, e1::Phase, e2::Phase, W, n; WS = 0.0, WV = 0.0, sites = 1.0)
    c1, c2 = e1.comp ./ sum(e1.comp), e2.comp ./ sum(e2.comp)
    map(range(0.0, 1.0, length = n)) do x
        Phase(name, (1 - x) .* c1 .+ x .* c2,
              (1 - x) * e1.g0 + x * e2.g0 + W * x * (1 - x),
              (1 - x) * e1.s  + x * e2.s  - sites * RGAS * (xlogx(x) + xlogx(1 - x)) + WS * x * (1 - x),
              (1 - x) * e1.v  + x * e2.v  + WV * x * (1 - x))
    end
end

"""
    assemblage_groups(phases, stable)

Group the stable phases into physical phases: neighbouring pseudocompounds of the
same solution (consecutive indices, same name) are one phase whose composition
lies between them. Two separated pseudocompounds of the same solution are two
coexisting phases (e.g. L + L: immiscibility). Returns a vector of index groups.
"""
function assemblage_groups(phases, stable)
    groups = Vector{Int}[]
    for i in sort(stable)
        if !isempty(groups) && groups[end][end] == i - 1 && phases[i-1].name == phases[i].name
            push!(groups[end], i)
        else
            push!(groups, [i])
        end
    end
    return groups
end

"Names of the physical phases of an assemblage, e.g. [\"A\", \"L\"] (used as key)."
assemblage_key(phases, stable) = [phases[g[1]].name for g in assemblage_groups(phases, stable)]

# ----------------------------------------------------------------------------
# Initialisation: scan a grid and locate the field boundaries
# ----------------------------------------------------------------------------
"""
    scan_grid(key, xs, ys)

Evaluate `key(x, y)` (anything that identifies the stable assemblage) at every
node of the grid xs × ys. Returns
  labels      Matrix{Int} (length(xs) × length(ys)), an integer id per assemblage
  assemblages Vector of keys; assemblages[id] belongs to label id
"""
function scan_grid(key::Function, xs, ys)
    ids         = Dict{Any,Int}()
    assemblages = Any[]
    labels      = zeros(Int, length(xs), length(ys))
    for (j, y) in enumerate(ys), (i, x) in enumerate(xs)
        k = key(x, y)
        labels[i, j] = get!(ids, k) do
            push!(assemblages, k)
            length(assemblages)
        end
    end
    return labels, assemblages
end

"P–T grid (Example 0): the key is the vector of stable phase indices."
scan_grid(lp::GibbsLP, Ts, Ps) = scan_grid((T, P) -> solve!(lp, P, T).stable, Ts, Ps)

"""
    refine_boundaries(key, xs, ys, labels; nbis=14)

Wherever two neighbouring grid nodes have different assemblages, a boundary
lies between them: bisect along that grid edge to find it. Returns the
boundary points as vectors (xb, yb). The finer the grid, the more points; the
more bisection steps, the more accurate each point.
"""
function refine_boundaries(key::Function, xs, ys, labels; nbis = 14)
    xb, yb = Float64[], Float64[]
    function bisect(x1, y1, x2, y2)
        k1 = key(x1, y1)
        # A node lying exactly on a boundary is degenerate (either assemblage is
        # optimal) and may be re-solved differently: then there is nothing to bisect.
        k1 == key(x2, y2) && return
        for _ in 1:nbis
            xm, ym = (x1 + x2) / 2, (y1 + y2) / 2
            if key(xm, ym) == k1
                x1, y1 = xm, ym
            else
                x2, y2 = xm, ym
            end
        end
        push!(xb, (x1 + x2) / 2); push!(yb, (y1 + y2) / 2)
    end
    nx, ny = size(labels)
    for j in 1:ny, i in 1:nx
        i < nx && labels[i, j] != labels[i+1, j] && bisect(xs[i], ys[j], xs[i+1], ys[j])
        j < ny && labels[i, j] != labels[i, j+1] && bisect(xs[i], ys[j], xs[i], ys[j+1])
    end
    return xb, yb
end

"""
    field_regions(labels)

Split the grid into connected regions of equal label (4-neighbours). A field can
occur in several separate places (e.g. the liquid on both sides of an
immiscibility gap). Returns a vector of (label, node indices) per region.
"""
function field_regions(labels)
    nx, ny  = size(labels)
    seen    = falses(nx, ny)
    regions = Tuple{Int,Vector{CartesianIndex{2}}}[]
    for I in CartesianIndices(labels)
        seen[I] && continue
        l, stack, nodes = labels[I], [I], CartesianIndex{2}[]
        seen[I] = true
        while !isempty(stack)
            J = pop!(stack); push!(nodes, J)
            for d in ((1, 0), (-1, 0), (0, 1), (0, -1))
                K = J + CartesianIndex(d)
                checkbounds(Bool, labels, K) && !seen[K] && labels[K] == l || continue
                seen[K] = true; push!(stack, K)
            end
        end
        push!(regions, (l, nodes))
    end
    return regions
end

refine_boundaries(lp::GibbsLP, Ts, Ps, labels; nbis = 14) =
    refine_boundaries((T, P) -> solve!(lp, P, T).stable, Ts, Ps, labels; nbis)

"""
    reaction_line(phases, i, j, k, Ts)

Analytic univariant line on which phase k lies exactly on the tie-line i–j in
g–x space (e.g. A + B = AB). Because g is linear in P and T this is the
straight line  Δg0 - Δs T + Δv P = 0.  Returns P for every T in Ts.
"""
function reaction_line(phases, i, j, k, Ts)
    x  = xcomp(phases)
    wi = (x[j] - x[k]) / (x[j] - x[i])              # lever rule: x_k = wi x_i + wj x_j
    wj = 1 - wi
    Δ(f) = f(phases[k]) - wi * f(phases[i]) - wj * f(phases[j])
    Δg0, Δs, Δv = Δ(p -> p.g0), Δ(p -> p.s), Δ(p -> p.v)
    return [(Δs * T - Δg0) / Δv for T in Ts]
end

"Lower convex hull of points (x, g) – the Gibbs energy of the system vs X."
function lower_hull(x::AbstractVector, g::AbstractVector)
    idx = sortperm(collect(zip(x, g)))
    h   = Int[]
    for i in idx
        # same x as the previous point (e.g. solid and liquid end-member): keep only the lowest g
        !isempty(h) && x[i] == x[h[end]] && continue
        while length(h) >= 2
            o, a = h[end-1], h[end]
            cross = (x[a] - x[o]) * (g[i] - g[o]) - (g[a] - g[o]) * (x[i] - x[o])
            cross <= 0 ? pop!(h) : break
        end
        push!(h, i)
    end
    return x[h], g[h]
end
