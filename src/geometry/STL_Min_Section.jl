# =============================================================================
#  STL_Min_Section.jl  —  geometry (stages 1, 5)
#  Minimum cross-section along `axis`: n planes at cell centres (spacing h = L/n),
#  then n_ref planes within one cell of the coarse minimum. For a dogbone, A_min is
#  the neck (ligament) area used for the nominal tensile stress and for
#  G_F ≈ W_F / A_F (Cusatis et al. 2011, LDPM Part II).
#  A straight gauge gives a plateau of minima; its extent [x_lo, x_hi] (± h/2) is
#  returned and the position is the plateau centre.
#  L_eff = ∫ A_min / A(x) dx (midpoint rule on the n planes): the elastic elongation
#  at gauge stress σ is (σ/E) L_eff; L_eff = L for a prism.
#  Returns (A_min [mm^2], position [mm], A_max [mm^2], L_eff [mm], x_lo, x_hi [mm]).
#  FAITHFUL-NOTE: identical to v4 STL_Tools.jl stl_min_section.
# =============================================================================

function STL_Min_Section(STL, axis::Int; n::Int = 200, n_ref::Int = 50, rtol::Float64 = 1e-6)
    lo = STL.Lo[axis]; hi = STL.Hi[axis]; L = hi - lo; h = L/n
    cs = [lo + (k - 0.5)*h for k in 1:n]
    As = [STL_Section_Area(STL, axis, c) for c in cs]
    minimum(As) > 0 || error("STL_Min_Section: non-positive section area; the surface orientation is inconsistent.")
    k = argmin(As)
    cr = collect(range(max(lo + 1e-6*L, cs[k] - h), min(hi - 1e-6*L, cs[k] + h); length = n_ref))
    Ar = [STL_Section_Area(STL, axis, c) for c in cr]
    j = argmin(Ar)
    A_min = min(Ar[j], As[k])
    plat = findall(a -> a <= A_min*(1 + rtol), As)
    if length(plat) >= 2
        x_lo = cs[first(plat)] - h/2; x_hi = cs[last(plat)] + h/2
        x_c  = (x_lo + x_hi)/2
    else
        x_lo = cr[j] - h/2; x_hi = cr[j] + h/2
        x_c  = cr[j]
    end
    L_eff = h*sum(A_min ./ As)
    return A_min, x_c, maximum(As), L_eff, x_lo, x_hi
end
