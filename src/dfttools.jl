
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


ωr(N, i) = cospi(-2(i % N)/N)
ωi(N, i) = sinpi(-2(i % N)/N)
ω(N, i) = ωr(N, i) + im*ωi(N, i)

ωr( ::Val{N}, ::Val{0} ) where N = 1.0
@generated ωr( ::Val{N}, ::Val{i} ) where {N,i} = ωr(N, i)
ωi( ::Val{N}, ::Val{0} ) where N = 0.0
@generated ωi( ::Val{N}, ::Val{i} ) where {N,i} = ωi(N, i)
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
        push!(e_arg.args, :(x[$(2i-1)]))
        push!(o_arg.args, :(x[$(2i)]))
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


function recursion_level(d_array, recursion_levels)
    end_recursion_size = d_array[recursion_levels + 1] #Smallest component in which we can break this down.

    step = 1 << recursion_levels
    initial_eval = quote end
    for i ∈ 1:step
        base_step!(initial_eval.args, end_recursion_size, i, step, d_array[1])
    end
    initial_eval
end

function base_step(blength, offset, step, N, x = :x_, r = :r_)
    base = :( $(Symbol(r, offset)) = $(Symbol(x, offset)) + $(Symbol(x, offset + step)) * $(ω(N, step)))
    for i ∈ 2:blength-1
        push!(base.args[2].args, :( $(Symbol(x, offset + i*step)) * $(ω(N, i*step)) ))
    end
    base
end

function base_step!(initial, blength, offset, step, N, x = :x_, r = :r_)
    for j ∈ 1:blength        
        push!(initial.args, base_step(blength, offset, step, N, x, r))
    end
    initial
end

function compose( level, step )


end


qe = quote
        x = a + b
        y = c * d
        z = x * y + 4
    end

function recursion_meta(n)
    recursion_levels = 0
    d, r = divrem(n, 2)
    d_array = [n,d]
    while r != 1
        recursion_levels += 1
        d, r = divrem(d, 2)
        push!(d_array, d)
        d == 1 && break
    end
    d_array, recursion_levels
end

@generated function dfte(x::NTuple{n,T}) where {n,T}

    quote
        @nextract $n x x

    end
end

@nexprs $N i -> begin
    @nexprs $P-i j -> begin
        i * j
    end
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





function dft_call_expr(N)
    d = N ÷ 2
    e_arg = :( dftva(x[1]) )
    o_arg = :( dftva(x[2]) )
    for i ∈ 2:d
        push!(e_arg.args, :(x[$(2i-1)]))
        push!(o_arg.args, :(x[$(2i)]))
    end
    e_arg, o_arg
end
@inline dftva(x::T) where T <: Number = x
@generated function dftva(x::Vararg{T,n}) where {n,T}
    out = dft_sum_expr(n, T)
    e_arg, o_arg = dft_call_expr(n)
    quote
        $(Expr(:meta, :inline))
        ffte = $e_arg
        ffto = $o_arg
        $out
    end
end

