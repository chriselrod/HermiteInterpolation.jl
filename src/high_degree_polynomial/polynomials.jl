



struct HighDegreeInterpPoly{N,V}
    c::NTuple{N,V}
    w::V
end

struct ExponentPairs{N} end
ExponentPairs(::Val{N}) where {N} = ExponentPairs{N}()
Base.@pure ExponentPairs(N) = ExponentPairs{N}()

"""
@inline cmod(n, d) = n & (d - one(d))
function Base.getindex(::ExponentPairs{N}, i::T) where {N, T}
    d, r = divrem(i-one(T), T(N))
    r += one(T)
    r, one(T) + cmod(cmod(d, (T(N)-one(T))) + r, T(N))
end

function Base.getindex(::ExponentPairs{N}, i::T) where {N, T}
    d, r = fldmod1(i, T(N))
    r,  mod1( mod1(d, T(N)-one(T)) + r, T(N) )
end
"""
function Base.getindex(::ExponentPairs{N}, i::T) where {N, T}
    d, r = divrem(i-one(T), T(N))
    r += one(T)
    r, one(T) + ((d % (T(N)-one(T))) + r) % T(N)
end

Base.Vector(::ExponentPairs{N}, n::T) where {T,N} = getindex.(ExponentPairs{N}(), one(T):n)
like_size(::Type{Float64}, x) = Int64(x)
like_size(::Type{Float32}, x) = Cuint(x)

function coefs(seq::A, ::Val{N}) where {N, T, A <: AbstractArray{T,2}}
    n = size(seq, 2) #n * N = length(seq)
    #Need to iterative over pairs.
    
    l = cld( 2n, N ) #+ 1


    #cycles
    num_cycles = 2cld(n, binomial(N, 2))
    base_polynomial = fill( one(Complex{T}), num_cycles, N)

    #F it, lets stick with a vector!
    #Hmm, we actually want to guarantee last set does not feature overlap.
    #
    k = Vector(ExponentPairs{N}(), like_size(T, n))
    for (i,j) ∈ k
        for w ∈ 1:num_cycles
            base_polynomial[:,i]
            base_polynomial[:,j]
        end
    end

    cycle1 = binomial(N, 2)
    cycles1, remainder1 = divrem(n, cycle1)


    
    for i ∈ x
        
    end
end

struct ComplexViewPolynomial{T,A,D,N} where {T,A<:AbstractArray{T},D,N}
    ω::NTuple{D,T}
    c::A
end

@generated gen_omega(::Val{D}, ::Type{T} = Float64) where {D,T} = :(@ntuple $D i -> Complex{T}(exp(-2π*im*(i-one(T))/D)))

@generated function gen_omega_r(::Val{D}, ::Type{T} = Float64) where {D,T}
    l = (D-1) ÷ 2 #Truncated 1st element ( = 1 ), and last non-symmetric element if D is even ( = -1 )
    :(@ntuple $l i -> T(cos(2π*i/D)))
end
@generated function gen_omega_i(::Val{D}, ::Type{T} = Float64) where {D,T}
    l = (D-1) ÷ 2 #Truncated 1st elment ( = 0 ) and last non-symmetric element if D is even ( = 0 )
    :(@ntuple $l i -> T(sin(2π*(i-one(T))/D)))
end


function complex_polynomial_kernel(state, A::AbstractArray{T}, ω) where T
    d,p = @cartesianidx A #poly coeficient, and the variable in the multi-dimensional polynomial.
    return
end

struct EvaluationKernel{D,P,N} end

@generated function ComplexViewPolynomial(roots::AbstractArray{T}, ::EvaluationKernel{D,P,N})
                where {D,P,N,T,A <: AbstractArray{T,2}}
    ωᵣ = gen_omega_r(Val{D}(), T)
    ωᵢ = gen_omega_i(Val{D}(), T)
    complete_cycles, leftover = divrem( N, P )
    half_P = P ÷ 2
    first_set, repeats = divrem(leftover, half_P)
    Podd = isodd(P)
    if first_set == 0
        first_set = leftover
        repeats = 0
    else
        first_set = half_P
        if Podd
            repeats -= 1 ##Now we need Podd && repeats >= 0 to do the Podd branch.
        end
    end
    #Now, first set can assume it is writing into the matrix for the first time.
    #Then, repeats may modify earlier entires.
    #If not, we will have to fill out the rest of the matrix with appropriate 1s (r) and 0s (i).


    Deven = iseven(D)
    D_clipped = ( D - 1 ) ÷ 2
    quote

        base_poly = fill(one(T), $D, $P)
        @nexprs $half_P i -> begin
            i₂ = 2i
            i₁ = mod1( i₂ - 1, $P )
            i₂ = mod1( i₂, $P )
            # Index the array in the correct order.
            # Does this honestly matter? Is it just about the principle of the thing?
            # Note: the branches are known at compile time, and thus removed.
            if i₁ < i₂ 
                r₁ = roots[i₁, i]
                r₂ = roots[i₂, i]

                base_poly[1, i₁] *= one(T) - r₁
                base_poly[1, i₂] *= one(T) - r₂
                
                @nexprs $D_clipped d -> begin
                    dᵣ = 2d
                    dᵢ = dᵣ+1
                    
                    ωrτ = fma($ωᵣ[d], -r₁ , one(T))
                    bpᵣ = base_poly[dᵣ, i₁]
                    bpᵢ = base_poly[dᵢ, i₁]
                    base_poly[dᵣ, i₁] = bpᵣ*ωrτ + bpᵢ*r₁*ωᵢ[d]
                    base_poly[dᵢ, i₁] = bpᵢ*ωrτ - bpᵣ*r₁*ωᵢ[d]
                    ωrτ = fma($ωᵣ[d], -r₂ , one(T))
                    bpᵣ = base_poly[dᵣ, i₂]
                    bpᵢ = base_poly[dᵢ, i₂]
                    base_poly[dᵣ, i₂] = bpᵣ*ωrτ + bpᵢ*r₂*ωᵢ[d]
                    base_poly[dᵢ, i₂] = bpᵢ*ωrτ - bpᵣ*r₂*ωᵢ[d]
                end
                if $Deven
                    base_poly[$D, i₁] *= one(T) + r₁
                    base_poly[$D, i₂] *= one(T) + r₂
                end
            else
                #reverse order
            end
            
        end
        @nexprs $leftover i -> begin
            i₂ = 2i
            i₁ = mod1( i₂ - 1, $P )
            i₂ = mod1( i₂, $P )
            # Index the array in the correct order.
            # Does this honestly matter? Is it just about the principle of the thing?
            # Note: the branches are known at compile time, and thus removed.
            if i₁ < i₂ 
                r₁ = roots[i₁, i]
                r₂ = roots[i₂, i]

                base_poly[1, i₁] *= one(T) - r₁
                base_poly[1, i₂] *= one(T) - r₂
                
                @nexprs $D_clipped d -> begin
                    dᵣ = 2d
                    dᵢ = dᵣ+1
                    
                    ωrτ = fma($ωᵣ[d], -r₁ , one(T))
                    bpᵣ = base_poly[dᵣ, i₁]
                    bpᵢ = base_poly[dᵢ, i₁]
                    base_poly[dᵣ, i₁] = bpᵣ*ωrτ + bpᵢ*r₁*ωᵢ[d]
                    base_poly[dᵢ, i₁] = bpᵢ*ωrτ - bpᵣ*r₁*ωᵢ[d]
                    ωrτ = fma($ωᵣ[d], -r₂ , one(T))
                    bpᵣ = base_poly[dᵣ, i₂]
                    bpᵢ = base_poly[dᵢ, i₂]
                    base_poly[dᵣ, i₂] = bpᵣ*ωrτ + bpᵢ*r₂*ωᵢ[d]
                    base_poly[dᵢ, i₂] = bpᵢ*ωrτ - bpᵣ*r₂*ωᵢ[d]
                end
                if $Deven
                    base_poly[$D, i₁] *= one(T) + r₁
                    base_poly[$D, i₂] *= one(T) + r₂
                end
            else
                #reverse order
            end
            
        end
        ComplexViewPolynomial($ω, )
    end
end

function MatrixFFT(n::Int, ::Type{T} = Float64) where T
    out = fill(one(Complex{T}), n, n)
    for j ∈ 2:n
        out[j, 2] = exp( -2π*im*(j-1) / n )
    end
    for i ∈ 3:n, j ∈ i:n
        out[j,i] = out[j,i-1] * out[j,2]
    end
    Symmetric(out, :L) #./= √n
end