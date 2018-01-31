
@inline modv2(x, n) = rem(x - one(x), n) + one(x)


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
    :(@ntuple $l i -> T(sin(2π*i/D)))
end


function complex_polynomial_kernel(state, A::AbstractArray{T}, ω) where T
    d,p = @cartesianidx A #poly coeficient, and the variable in the multi-dimensional polynomial.
    return
end

"""
D really ought to be
D = cld(2N, P)+1
"""
struct EvaluationKernel{D,P,N,A} end

@generated function split_roots(roots::A, ::EvaluationKernel{D,P,N})
    where {D,P,N,T,A <: AbstractArray{T,2}}
    complete_cycles, leftover = divrem( N, P )

    quote
        root_poly = Array{T,2}(uninitialized, 2*$complete_cycles, $P)
        ind = 0
        @nexprs $P j -> begin
            @nexprs $complete_cycles i -> begin
                ind += 1
                k = modv2( j + i + cld( i, $P - 1) - 1 , $P )
                root_poly[2(i-1)+1,j] = roots[j,ind]
                root_poly[2i,k] = roots[k,ind]
            end
        end
        root_poly
    end
end
@generated function split_roots_and_inds(roots::A, ::EvaluationKernel{D,P,N})
                where {D,P,N,T,A <: AbstractArray{T,2}}
    complete_cycles, leftover = divrem( N, P )
    
    quote
        root_poly = Array{T,2}(uninitialized, 2*$complete_cycles, $P)
        root_inds = Array{T,2}(uninitialized, $complete_cycles, $P)
        ind = 0
        @nexprs $P j -> begin
            @nexprs $complete_cycles i -> begin
                ind += 1
                k = modv2( j + i + cld( i, $P - 1) - 1 , $P )
                root_inds[2(i-1)+1,j] = ind
                root_inds[2i,k] = ind
                root_poly[2(i-1)+1,j] = roots[j,ind]
                root_poly[2i,k] = roots[k,ind]
            end
        end
        root_poly,root_inds
    end
end
@generated function split_root_inds(::EvaluationKernel{D,P,N}, ::Type{T} = Cuint)
                where {D,P,N,T}
    complete_cycles, leftover = divrem( N, P )
    
    quote
        root_inds = Array{T,2}(uninitialized, $complete_cycles, $P)
        ind = 0
        @nexprs $complete_cycles i -> begin
            @nexprs $P j -> begin
                ind += 1
                k = modv2( j + i + cld( i, $P - 1) - 1 , $P )
                root_inds[2(i-1)+1,j] = ind
                root_inds[2i,k] = ind
            end
        end
        root_inds
    end
end
@inline function update_poly!(base_poly, ωᵣ, ωᵢ, d, i, r)
    dᵣ = 2d
    dᵢ = dᵣ+1    
    ωrτ = fma($ωᵣ[d], -r , one(eltype(base_poly)))
    bpᵣ = base_poly[dᵣ, i]
    bpᵢ = base_poly[dᵢ, i]
    base_poly[dᵣ, i] = bpᵣ*ωrτ + bpᵢ*r*ωᵢ[d]
    base_poly[dᵢ, i] = bpᵢ*ωrτ - bpᵣ*r*ωᵢ[d]
end
@inline function initial_write_poly!(base_poly, ωᵣ, ωᵢ, d, i, r)
    base_poly[2d, i] = fma($ωᵣ[d], -r , one(eltype(base_poly)))
    base_poly[2d+1, i] = - r*ωᵢ[d]
end
@inline function initial_write_poly!(base_poly, d, i)
    base_poly[2d, i] = one(eltype(base_poly))
    base_poly[2d+1, i] = zero(eltype(base_poly))
end

@generated function ComplexViewPolynomial(roots::AbstractArray{T,2}, ::EvaluationKernel{D,P,N,A})
                where {D,P,N,T,A <: AbstractArray}
    ωᵣ = gen_omega_r(Val{D}(), T)
    ωᵢ = gen_omega_i(Val{D}(), T)
    complete_cycles, leftover = divrem( N, P )
    offset = complete_cycles * P
    half_P, Podd = divrem(P, 2)
    first_set, repeats = divrem(leftover, half_P)
    if first_set == 0
        first_set = leftover
        finish_initialization = half_P - first_set
        repeats = 0
        Podd = 0
    else
        first_set = half_P
        finish_initialization = 0
        if repeats != 0
            repeats -= Podd ##Now we need Podd && repeats >= 0 to do the Podd branch.
        else
            Podd = 0
        end
    end
    #Now, first set can assume it is writing into the matrix for the first time.
    #Then, repeats may modify earlier entires.
    #If not, we will have to fill out the rest of the matrix with appropriate 1s (r) and 0s (i).


    D_clipped, Deven = divrem( D - 1, 2)#Deven asks whether we have a trailing real.
    quote
        base_poly = Array{T,2}(uninitialized, $D, $P)
        #first_set is set to however many need to be initialized.
        @nexprs $first_set i -> begin
            ind = i+$offset
            p₁ = mod1( 2ind - 1, $P )
            p₂ = mod1( 2ind, $P )
            if p₁ < p₂ 
                r₁ = roots[p₁, ind]
                r₂ = roots[p₂, ind]

                base_poly[1, p₁] = one(T) - r₁
                base_poly[1, p₂] = one(T) - r₂
                
                @nexprs $D_clipped d -> begin
                    #Originally, I expanded things here.
                    #But that was ugly and unclear.
                    #Inlining + constant propagation mean ugliness does not help performance.
                    #If ugliness does help, some argue we should avoid it anyway.
                    #The clean way is likely to catch up, be more maintainable, etc.
                    #On average, algorithmic advantages >>>> micro optimizations.
                    #Algorithmic advantages are enabled by code we understand.
                    initial_write_poly!(base_poly, ωᵣ, ωᵢ, d, p₁, r₁)
                    initial_write_poly!(base_poly, ωᵣ, ωᵢ, d, p₂, r₂)
                end
                if $Deven == 1
                    base_poly[$D, p₁] = one(T) + r₁
                    base_poly[$D, p₂] = one(T) + r₂
                end
            else
                r₂ = roots[p₂, ind]
                r₁ = roots[p₁, ind]

                base_poly[1, p₂] = one(T) - r₂
                base_poly[1, p₁] = one(T) - r₁
                @nexprs $D_clipped d -> begin
                    initial_write_poly!(base_poly, ωᵣ, ωᵢ, d, p₂, r₂)
                    initial_write_poly!(base_poly, ωᵣ, ωᵢ, d, p₁, r₁)
                end
                if $Deven == 1
                    base_poly[$D, p₂] = one(T) + r₂
                    base_poly[$D, p₁] = one(T) + r₁
                end
            end
        end
        # These are dimensions that have not been hit yet.
        @nexprs $finish_initialization i -> begin
            ind = i+$offset+$first_set
            p₁ = mod1( 2ind - 1, $P )
            p₂ = mod1( 2ind, $P )
            if p₁ < p₂#Need to set entries to one and zero, as appropriate.
                base_poly[1, p₁] = one(T)
                base_poly[1, p₂] = one(T)
                @nexprs $D_clipped d -> begin
                    initial_write_poly!(base_poly, d, p₁)
                    initial_write_poly!(base_poly, d, p₂)
                end
                if $Deven == 1
                    base_poly[$D, p₂] = one(T)
                    base_poly[$D, p₁] = one(T)
                end
            else
                base_poly[1, p₂] = one(T)
                base_poly[1, p₁] = one(T)
                @nexprs $D_clipped d -> begin
                    initial_write_poly!(base_poly, d, p₂)
                    initial_write_poly!(base_poly, d, p₁)
                end
                if $Deven == 1
                    base_poly[$D, p₂] = one(T)
                    base_poly[$D, p₁] = one(T)
                end
            end
        end
        
        if $Podd == 1
            ind = 1+$offset+$first_set
            p₁ = mod1( 2ind - 1, $P ) ###This has to be the one that was not hit before.
            p₂ = mod1( 2ind, $P )
            r₂ = roots[p₂, ind]
            if p₁ < p₂#Need to set entries to one and zero, as appropriate.
                @nexprs $D_clipped d -> begin
                    initial_write_poly!(base_poly, d, p₁)
                    update_poly!(base_poly, ωᵣ, ωᵢ, d, p₂, r₂)
                end
                if $Deven == 1
                    base_poly[$D, p₁] = one(T)
                    base_poly[$D, p₂] = one(T) + r₂
                end
            else
                @nexprs $D_clipped d -> begin
                    update_poly!(base_poly, ωᵣ, ωᵢ, d, p₂, r₂)
                    initial_write_poly!(base_poly, d, p₁)
                end
                if $Deven == 1
                    base_poly[$D, p₂] = one(T) + r₂
                    base_poly[$D, p₁] = one(T)
                end
            end

        end

        @nexprs $repeats i -> begin
            ind = i+$offset+$first_set+$Podd
            p₁ = mod1( 2ind - 1, $P )
            p₂ = mod1( 2ind, $P )
            if p₁ < p₂ 
                r₁ = roots[p₁, ind]
                r₂ = roots[p₂, ind]
                base_poly[1, p₁] *= one(T) - r₁
                base_poly[1, p₂] *= one(T) - r₂
                @nexprs $D_clipped d -> begin
                    update_poly!(base_poly, ωᵣ, ωᵢ, d, p₁, r₁)
                    update_poly!(base_poly, ωᵣ, ωᵢ, d, p₂, r₂)
                end
                if $Deven == 1
                    base_poly[$D, p₁] *= one(T) + r₁
                    base_poly[$D, p₂] *= one(T) + r₂
                end
            else
                r₂ = roots[p₂, ind]
                r₁ = roots[p₁, ind]
                base_poly[1, p₂] *= one(T) - r₂
                base_poly[1, p₁] *= one(T) - r₁
                @nexprs $D_clipped d -> begin
                    update_poly!(base_poly, ωᵣ, ωᵢ, d, p₂, r₂)
                    update_poly!(base_poly, ωᵣ, ωᵢ, d, p₁, r₁)
                end
                if $Deven == 1
                    base_poly[$D, p₂] *= one(T) + r₂
                    base_poly[$D, p₁] *= one(T) + r₁
                end
            end
        end

        roots_sorted = split_roots(roots, EvaluationKernel{$D,$P,$N}())
        for p ∈ 1:$P, i ∈ 1:2*$complete_cycles
            r = roots_sorted[i,p]
            base_poly[1, p] *= one(T) - r
            if $Deven == 1
                base_poly[$D, p] *= one(T) + r
            end
        end

        base_poly_a = A(base_poly)
        nodes = A(roots)
        roots_sorted_a = A(roots_sorted)

        ComplexViewPolynomial($ω, p)
    end
end
"""
roots are PxN
basis_poly is Dx2PxN
base_poly is DxP
"""
@generated function generate_basis_kernel(state, root_inds::GPUArray{Tuple{Cuint,Cuint}} roots::GPUArray, basis_poly::GPUArray{Float32,4}, base_poly::GPUArray, ωᵣ::NTuple{L,T}, ωᵢ::NTuple{L,T}, ::EvaluationKernel{D,P,N}) where {D,P,N, L, T}
    complete_cycles = N ÷ P
    L, Deven = divrem(D-1, 2)
    D_two_factor = Deven
    if Deven == 1
        for i ∈ 1:trunc(Int,log2(8))
            if D >> D_two_factor <<  D_two_factor != D
                D_two_factor -= 1
                break
            else
                D_two_factor += 1
            end
        end
    end
    quote
        n = @linearidx rootinds
        @nextract $P p i -> root_inds[Cuint(i),n] #We skip first two p_inds are factored out.

        @nexprs 2 j -> begin
            r = roots[p_j,n]
            basis_poly[one(Cuint),p_j,n] = base_poly[one(Cuint),p_j] / ( one(T) - r )
            @nexprs $L i -> begin
                realpart = one(T) - r * ωᵣ[Cuint(i)]
                imagpart = r * ωᵢ[Cuint(i)]
                bpreal = base_poly[Cuint(2)*Cuint(i),p_j]
                bpimag = base_poly[one(Cuint)+Cuint(2)*Cuint(i),p_j]

                denom = abs2(realpart) + abs2(imagpart)
                basis_poly[Cuint(2)*Cuint(i),p_j,n] = ( bpreal * realpart + bpimag * imagpart ) / denom
                basis_poly[one(Cuint)+Cuint(2)*Cuint(i),p_j,n] = ( bpimag * realpart - bpreal * imagpart ) / denom
            end
            if $Deven == 1
                basis_poly[$D,p_j,n] = base_poly[$D,p_j] / ( one(T) + r )
            end
        end
        
        @nexprs $P-2 j -> begin
            @nexprs $D d -> begin
                basis_poly[d,p_{j+2},n] = base_poly[d,p_{j+2}]
            end
        end
        @nexprs 2 j -> begin
            @nexprs $D d -> begin
                basis_poly[d,p_j+$P,n] = base_poly[d,p_j]
            end
        end
        @nexprs $P-2 j -> begin
            r = roots[p_{j+2},n]
            basis_poly[one(Cuint),$P+p_{j+2},n] = basis_poly[one(Cuint),p_{j+2},n] * ( one(T) - r )
            @nexprs $L i -> begin
                realpart = one(T) - r * ωᵣ[Cuint(i)]
                imagpart = r * ωᵢ[Cuint(i)]
                bpreal = basis_poly[Cuint(2)*Cuint(i),p_{j+2},n]
                bpimag = basis_poly[one(Cuint)+Cuint(2)*Cuint(i),p_{j+2},n]

                basis_poly[Cuint(2)*Cuint(i),$P+p_{j+2},n] = bpreal * realpart - bpimag * imagpart
                basis_poly[one(Cuint)+Cuint(2)*Cuint(i),$P+p_{j+2},n] = bpimag * realpart + bpreal * imagpart
            end
            if $Deven == 1
                basis_poly[$D,$P+p_{j+2},n] = basis_poly[$D,p_{j+2}] / ( one(T) + r )
            end
        end
        ### Now you have to ifft all of these.
        ### For now, this is just the naive algorithm.
        ### An FFT algorithm, like Cooley-Tukey, should be implemented later if this is actually a good idea.
        ### Should be easy to code an optimal one.
        ### May need to adjust D automatically via some heurstiics on number of calculations.
        @nexprs $P j -> begin
            @nexprs $D k -> begin
                @nexprs l -> begin
                    
                    @nexprs 2 i -> begin

                    end
                end
            end
        end
        nothing
    end
end

@inline function base_poly_ind(n::T, d::T) where T
    t, r = divrem(n-one(T), d)
    t + one(T), r + one(T)
end
"""
Generates base polynomial.
L = (D-1) ÷ 2
"""
@generated function generate_base_poly_kernel(state, base_poly::GPUArray, root_split::GPUArray, ωᵣ::NTuple{L,T}, ωᵢ::NTuple{L,T}, ::EvaluationKernel{D,P,N}) where {D,P,N,L,T}
    complete_cycles = N ÷ P
    quote
        p, d = base_poly_ind(GPUArrays.linear_index(state), Cuint($L))
        ωr = ωᵣ[p]
        ωi = ωᵢ[p]
        pᵣ = T(2)*p
        pᵢ = pᵣ + one(T)
        cᵣ = base_poly[pᵣ, d]
        cᵢ = base_poly[pᵢ, d]
        @nexprs 2*$complete_cycles rind -> begin
            r = root_split[Cuint(rind), d]
            rᵣ = one(Cuint) - r*ωr
            rᵢ = r*ωi
            cᵣ, cᵢ = ( rᵣ*cᵣ - rᵢ*cᵢ, rᵣ*cᵢ + rᵢ*cᵣ )
        end
        base_poly[pᵣ, d] = cᵣ
        base_poly[pᵢ, d] = cᵢ
        nothing
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

