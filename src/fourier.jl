export fft_nd, ifft_nd

@doc """
    fft_nd(A::AbstractArray)
    fft_nd(A::AbstractArray [; dims])

Compute the `D`-dimensional discrete Fourier transform (DFT) of an array `A`,
where `D = ndims(A)`, by explicitly applying a one-dimensional FFT along each
dimension in sequence.
Conventions follow the FFTW and GenericFFT libraries.

If the optional keyword argument `dims` is provided (as any iterable collection),
then the list of dimensions in `dims` will be transformed.
If `dims` is not provided, all dimensions are transformed.

This function performs an in-place–style dimensional sweep:
each dimension is permuted to the front, reshaped into a matrix,
transformed with a 1D FFT along the first axis, and then reshaped back.
After all dimensions are processed, the original dimension order is restored.
"""
function fft_nd(A::AbstractArray;
    dims=nothing)
    D = ndims(A)
    s = size(A)
    N = length(A)

    # Normalize dims to a tuple of integers
    if dims === nothing
        dims_tuple = ntuple(identity, D)
    elseif dims isa Integer
        dims_tuple = (Int(dims),)
    else # Handle any iterable
        dims_tuple = Tuple(Int(d) for d in dims)
    end

    # Validate
    @assert all(d -> 1 ≤ d ≤ D, dims_tuple) "All dimensions must be in range 1:$D"
    dims_set = Set(dims_tuple)

    # p tracks the net permutation of dimensions
    p = ntuple(identity, D)
    for d in dims_set
        k = findfirst(==(d), p)
        q = (k, setdiff(1:D, k)...)
        A = permutedims(A, q)
        p = ntuple(i -> p[q[i]], D)
        rest = div(N, s[d])
        A = reshape(A, s[d], rest)
        A = fft(A, 1)
        s_perm = ntuple(i -> s[p[i]], D)
        A = reshape(A, s_perm...)
    end

    A = permutedims(A, invperm(p))
    return A
end


@doc """
    ifft_nd(A::AbstractArray)
    ifft_nd(A::AbstractArray [; dims])

Compute the `D`-dimensional inverse discrete Fourier transform (DFT) of an array `A`,
where `D = ndims(A)`, by explicitly applying a one-dimensional inverse FFT along each
dimension in sequence.
Conventions follow the FFTW and GenericFFT libraries.

If the optional keyword argument `dims` is provided (as any iterable collection),
then the list of dimensions in `dims` will be transformed.
If `dims` is not provided, all dimensions are transformed.

This function performs an in-place–style dimensional sweep:
each dimension is permuted to the front, reshaped into a matrix,
transformed with a 1D inverse FFT along the first axis, and then reshaped back.
After all dimensions are processed, the original dimension order is restored.
"""
function ifft_nd(A::AbstractArray;
    dims=nothing)
    D = ndims(A)
    s = size(A)
    N = length(A)

    # Normalize dims to a tuple of integers
    if dims === nothing
        dims_tuple = ntuple(identity, D)
    elseif dims isa Integer
        dims_tuple = (Int(dims),)
    else # Handle any iterable
        dims_tuple = Tuple(Int(d) for d in dims)
    end

    # Validate
    @assert all(d -> 1 ≤ d ≤ D, dims_tuple) "All dimensions must be in range 1:$D"
    dims_set = Set(dims_tuple)

    # p tracks the net permutation of dimensions
    p = ntuple(identity, D)
    for d in dims_set
        k = findfirst(==(d), p)
        q = (k, setdiff(1:D, k)...)
        A = permutedims(A, q)
        p = ntuple(i -> p[q[i]], D)
        rest = div(N, s[d])
        A = reshape(A, s[d], rest)
        A = ifft(A, 1)
        s_perm = ntuple(i -> s[p[i]], D)
        A = reshape(A, s_perm...)
    end

    A = permutedims(A, invperm(p))
    return A
end
