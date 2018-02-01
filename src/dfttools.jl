
function MatrixDFT(n::Int, ::Type{T} = Float64) where T
    out = fill(one(Complex{T}), n, n)
    for j ∈ 2:n
        out[j, 2] = exp( -2π*im*(j-1) / n )
    end
    for i ∈ 3:n, j ∈ i:n
        out[j,i] = out[j,i-1] * out[j,2]
    end
    Symmetric(out, :L) #./= √n
end


### Broadly, the strategy used here is to try and move as much of the work to compile time as possible.
### In particular, calculating the values of cos and sin at different points on the unit circle.
### This is enabled by aggressive loop unrolling in the functions, so all terms are known.
### This approach probably does not scale well to larger problems.


ωr( ::Val{N}, ::Val{0} ) where N = 1.0
@generated function ωr( ::Val{N}, ::Val{i} ) where {N,i}
    r = i % N
    ###This is mostly an aesthetic choice, so that cosine evaluates to exactly zero when it should, instead of ≈ eps().
    if abs(4r) == N || abs(4r) == 3N 
        out = 0.0
    elseif abs(2r) == N
        out = -1.0
    elseif r == 0
        out = 1.0
    else
        out = cos(-2π*r/N)
    end
    out
end
ωi( ::Val{N}, ::Val{0} ) where N = 0.0
@generated function ωi( ::Val{N}, ::Val{i} ) where {N,i}
    r = i % N
    if 4r == N || 4r == -3N
        out = -1.0
    elseif 4r == 3N || 4r == -N
        out = 1.0
    elseif abs(2r) == N
        out = 0.0
    else
        out = sin(-2π*r/N)
    end
    out
end
@generated ωr( ::Val{N}, ::Val{i} ) where {N,i} = cos(-2π*i/N)
@generated ωi( ::Val{N}, ::Val{i} ) where {N,i} = sin(-2π*i/N)
ω( ::Val{N}, ::Val{0} ) where N = 1.0 + 0.0im 
@generated ω( ::Val{N}, ::Val{i} ) where {N,i} = ωr(Val{N}(), Val{i}()) + im*ωi(Val{N}(), Val{i}())


function dft_sum_expr(N, ::Type{T} = Float64) where T
    d = N ÷ 2
    out = :( (ffte[1] + ffto[1], ) )
    for i ∈ 2:N
        m = mod1(i,d)
        g = Complex{T}(ω(Val{N}(), Val{i-1}()))
        expr = :( ffte[$m] + $g * ffto[$m])
        push!(out.args, expr )
    end
    out
end

function dft_even_odd_split_expr(N)
    d = N ÷ 2
    e_arg = :( (x[1], ) )
    o_arg = :( (x[2], ) )
    for i ∈ 2:d
        inde = 2i-1
        indo = 2i
        push!(e_arg.args, :(x[$inde]))
        push!(o_arg.args, :(x[$indo]))
    end
    e_arg, o_arg
end


@inline dft(x::NTuple{1,T}) where T = x
@generated function dft(x::NTuple{n,T}) where {n,T}
    out = dft_sum_expr(n, T)
    e_arg, o_arg = dft_even_odd_split_expr(n)
    quote
        $(Expr(:meta, :inline))
        ffte = dft($e_arg)
        ffto = dft($o_arg)
        $out
    end
end


@inline function dft(x::NTuple{2,T}) where N
    ωr(Val{1}())
end



@generated function dft(x::NTuple{n,T}, ::Val{N} = Val{n}()) where {n,T,N}
    d, r = divrem(N, 2)

    quote
        $(Expr(:meta, :inline))
        @nexprs $d r -> begin

        end
    end
end
@inline function dft(x::NTuple{1,T}, ::Val{N} = Val{1}()) where {T,N} = x
@inline function dft(x::NTuple{2,T}, ::Val{N} = Val{2}()) where {T,N}
    ωr(Val{1}())
end
