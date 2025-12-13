export ghost_basis


"""
    NormalBasisRecipe

Container encoding the algebraic recipe used to generate a normal basis
along a given axis.
Only certain recipes are supported at this time.

Fields
- kind  :: Symbol
    One of :linear, :shifted_quad, :shifted_quad_sq

- shift :: Union{Nothing, NTuple}
    The translation τ used in shifted constructions, or `nothing` for :linear.
"""
struct NormalBasisRecipe
    kind::Symbol
    shift::Union{Nothing,Tuple}
end

# Test for non-vanishing Fourier coefficients at Float64 precision.
# This is a logical test for algebraic non-zeroness.
function _fourier_nonvanishing(x; threshold=1e-12)
    p = ifft(Float64.(x))
    return minimum(abs.(p)) > threshold
end


# Discover a valid normal-basis recipe for a single axis j.
function _discover_axis_recipe(K, j; verbose=false)
    ords = size(K)
    r = ndims(K)
    N = length(K)
    notj = Tuple(k for k in 1:r if k != j)

    # --- linear candidate ---
    x = dropdims(sum(K; dims=notj); dims=notj)
    if _fourier_nonvanishing(x)
        verbose && println("Axis $j: linear basis accepted.")
        return NormalBasisRecipe(:linear, nothing)
    end

    # --- shifted quadratic ---
    for m in 0:N-1
        τ = radix(m, ords)
        x = dropdims(sum(K .+ K .* circshift(K, τ); dims=notj); dims=notj)
        if _fourier_nonvanishing(x)
            verbose && println("Axis $j: shifted quadratic accepted with τ = $τ.")
            return NormalBasisRecipe(:shifted_quad, Tuple(τ))
        end
    end

    # --- shifted quadratic squared ---
    for m in 0:N-1
        τ = radix(m, ords)
        x = dropdims(sum((K .+ K .* circshift(K, τ)) .^ 2; dims=notj); dims=notj)
        if _fourier_nonvanishing(x)
            verbose && println("Axis $j: shifted quadratic squared accepted with τ = $τ.")
            return NormalBasisRecipe(:shifted_quad_sq, Tuple(τ))
        end
    end

    error("Failed to find a valid normal-basis recipe for axis $j.")
end


# Discover normal-basis recipes for all axes.
function _discover_normal_basis_recipes(K; verbose=false)
    r = ndims(K)
    recipes = Vector{NormalBasisRecipe}(undef, r)
    for j in 1:r
        recipes[j] = _discover_axis_recipe(K, j; verbose=verbose)
    end
    return recipes
end


# Evaluate a normal basis along axis j from a cached recipe.
function _eval_normal_basis_axis(K, j, recipe::Zauner.NormalBasisRecipe)
    ords = size(K)
    r = ndims(K)
    notj = Tuple(k for k in 1:r if k != j)

    if recipe.kind === :linear
        x = dropdims(sum(K; dims=notj); dims=notj)

    elseif recipe.kind === :shifted_quad
        τ = recipe.shift
        x = dropdims(sum(K .+ K .* circshift(K, τ); dims=notj); dims=notj)

    elseif recipe.kind === :shifted_quad_sq
        τ = recipe.shift
        x = dropdims(sum((K .+ K .* circshift(K, τ)) .^ 2; dims=notj); dims=notj)

    else
        error("Unknown normal-basis recipe kind $(recipe.kind).")
    end

    p = ifft(x)
    y = real.(ifft(1 ./ (ords[j] .* p)))
    return x, p, y
end


"""
    ghost_basis(K; recipes=nothing, verbose=false)

Compute normal bases and dual bases from the ghost overlap tensor `K`.

If `recipes` is not provided, a valid normal-basis recipe is discovered
(at low precision) and returned alongside the bases.

Returns
- x :: Vector{Vector{BigFloat}}          (normal bases)
- p :: Vector{Vector{Complex{BigFloat}}} (Fourier bases)
- y :: Vector{Vector{BigFloat}}          (dual bases)
- recipes :: Vector{NormalBasisRecipe}   (cached recipes)
"""
function ghost_basis(
    K::AbstractArray{BigFloat};
    recipes::Union{Nothing,Vector{NormalBasisRecipe}}=nothing,
    verbose=false
)
    if recipes === nothing
        verbose && println("Discovering normal-basis recipes.")
        recipes = _discover_normal_basis_recipes(K; verbose=verbose)
    end

    r = ndims(K)
    x = Vector{Vector{BigFloat}}(undef, r)
    p = Vector{Vector{Complex{BigFloat}}}(undef, r)
    y = Vector{Vector{BigFloat}}(undef, r)

    for j in 1:r
        x[j], p[j], y[j] = _eval_normal_basis_axis(K, j, recipes[j])
    end

    return x, p, y, recipes
end
