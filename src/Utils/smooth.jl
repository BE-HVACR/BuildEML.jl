"""
    smooth_H(x, δ, t0)
    smooth_xH(x, δ, t0)

Basic smooth step helpers.

`smooth_H` is a smooth approximation of a Heaviside step centered at `t0`.
`smooth_xH` is a smooth approximation of `max(x - t0, 0)` with transition width `δ`.
"""
smooth_H(x, δ, t0) = 0.5 * (1 + (x - t0) / sqrt((x - t0)^2 + δ^2))
smooth_xH(x, δ, t0) = (x - t0) * smooth_H(x, δ, t0)

"""
    smooth_max(a, b; δ=1e-6)
    smooth_min(a, b; δ=1e-6)
    smooth_abs(x; δ=1e-6)
    smooth_sign(x; δ=1e-6)
    smooth_clamp(x, lo, hi; δ=1e-6)
    smooth_clamp01(x; δ=1e-6)

Common smooth algebraic helpers.

These functions provide smooth approximations of `max`, `min`, `abs`, `sign`, and interval clamping operations.
"""
smooth_max(a, b; δ = 1e-6) = b + smooth_xH(a - b, δ, 0.0)
smooth_min(a, b; δ = 1e-6) = a + b - smooth_max(a, b; δ = δ)
smooth_abs(x; δ = 1e-6) = smooth_max(x, -x; δ = δ)
smooth_sign(x; δ = 1e-6) = x / sqrt(x^2 + δ^2)
smooth_clamp(x, lo, hi; δ = 1e-6) = smooth_max(lo, smooth_min(x, hi; δ = δ); δ = δ)
smooth_clamp01(x; δ = 1e-6) = smooth_clamp(x, 0.0, 1.0; δ = δ)

"""
    smooth_on_eps(u; u_on=1.0, u_eps=1e-5, δ=1e-2, u_th=nothing)

Map `u` to a smooth activation level between `u_eps` and `u_on`.
"""
function smooth_on_eps(u; u_on = 1.0, u_eps = 1e-5, δ = 1e-2, u_th = nothing)
    s = smooth_clamp01(u; δ = δ)
    val = u_eps + (u_on - u_eps) * s
    return val
end

"""
    α_soft(sel; δ=1e-6)
    soft_blend(x_pos, x_neg, sel; δ=1e-6)

Smooth selector helpers.

`α_soft` returns a smooth selector weight based on the sign of `sel`.
`soft_blend` uses that weight to blend between `x_pos` and `x_neg`.
"""
α_soft(sel; δ = 1e-6) = 0.5 * (1 + sel / sqrt(sel^2 + δ^2))
soft_blend(x_pos, x_neg, sel; δ = 1e-6) = α_soft(sel; δ = δ) * x_pos + (1 - α_soft(sel; δ = δ)) * x_neg

"""
    snapconst(v; tol=1e-9)

Snap a nearly constant vector to an exactly constant vector for result
display. This helper is intended for post-processing rather than symbolic
model construction.
"""
snapconst(v; tol = 1e-9) = (maximum(v) - minimum(v) <= tol) ? fill(first(v), length(v)) : v

