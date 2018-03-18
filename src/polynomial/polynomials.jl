#using Reduce, SpecialFunctions, Base.Cartesian, Sobol
#Reduce.Rational(false)

@generated function marginalize(x, ::Tuple{Val{N}}) where N
    quote
        out = zero(eltype(x))
        @nexprs $N i -> @inbounds out += x[i] *  marginal(Val{i-1}())
        out
    end
end

@generated function marginalize(x, ::Tuple{Val{N1},Val{N2}}) where {N1, N2}
    quote
        out = zero(eltype(x))
        @nexprs $N1 i -> begin
            @nexprs $N2 j -> begin
                out += x[j,i] * marginalize(Val{i-1}()) * marginalize(Val{j-1}())
            end
        end
        out
    end
end
@generated function marginalize(x, ::Tuple{Val{N1},Val{N2},Val{N3}}) where {N1, N2, N3}
    quote
        out = zero(eltype(x))
        @nexprs $N1 i -> begin
            @nexprs $N2 j -> begin
                @nexprs $N3 k -> begin
                    out += x[k,j,i] * marginalize(Val{i-1}()) * marginalize(Val{j-1}()) * marginalize(Val{k-1}())
                end
            end
        end
        out
    end
end
@generated function marginalize(x, ::Tuple{Val{N1},Val{N2},Val{N3},Val{N4}}) where {N1, N2, N3, N4}
    quote
        out = zero(eltype(x))
        @nexprs $N1 i -> begin
            @nexprs $N2 j -> begin
                @nexprs $N3 k -> begin
                    @nexprs $N4 l -> begin
                        out += x[l,k,j,i] * marginalize(Val{i-1}()) * marginalize(Val{j-1}()) * marginalize(Val{k-1}()) * marginalize(Val{l-1}())
                    end
                end
            end
        end
        out
    end
end

@generated function marginal_gen(::Val{N}, ::Type{T} = Float32) where {N,T}
    f = rcall(:(int( x^$N * exp(-x^2/2), x, -4.8, 4.8)))
    :(convert(T,$f))
end
@generated function marginal(::Val{N}, ::Type{T} = Float32) where {N,T}
    marginal_gen(Val{N}(), T)
end
function marginal(i::I, ::Type{T} = Float32) where {I <: Integer, T}
    out = @nif 16 d -> (i == 2(d-1)) d -> ( marginal(Val{2(d-1)}(), T) ) d -> ( zero(T) )
    out
end

@generated function uni_poly(x, coef_v, ::Val{N}) where N
    Nm1 = N - 1
    quote
        @inbounds out = coef_v[1]
        @nexprs $Nm1 i -> @inbounds out += x^i * coef_v[i+1]
        out * exp(-x^2/2)
    end
end


function coef(i::I, ::Type{T} = Float32) where {I <: Integer, T}
    out = @nif 160 d -> (i <= 2(d-1)) d -> ( 1/marginal(Val{2(d-1)}(), T) ) d -> ( zero(T) )
    out
end
function coef(p, T = Float32)
    c = similar(p, T)
    for i ∈ eachindex(c)
        c[i] = prod(x -> coef(x, T), p[i])
    end
    c
end

"""
a is the Array (nested VectorOfArray) form of the polynomial.
v is a vector of coefficients, typically created by solving a system of equations.
"""
poly(::Val{0}, n::NTuple{ns,I} = (32,), ::Type{S} = Float32) where {ns,S,I<:Integer} = fill(zero(S), n...)
@generated function poly(::Val{N}, n::NTuple{ns,I} = (32,), ::Type{S} = Float32) where {N,ns,S,I<:Integer}
    Nm1 = N - ns
    v_last = Symbol("v_", Nm1)
    quote
        v_0 = fill(zero($S), n...)#In preparation for the deprecation of zeros.
        @nextract $Nm1 v d -> [ v_{d-1} ]
        $v_last
    end
end

@generated function poly(v_coefs::Vector{T}, poly_exponents::Vector{NTuple{p,Cuint}}, n::NTuple{ns,I} = (32,), ::Type{S} = Float32) where {ns,p,T,S,I<:Integer}
    N = p - ns
    v_last = Symbol("v_", N)
    quote
        v_0 = poly(Val{$p}(), n, S)
        for i ∈ eachindex(v_coefs)
            ex = poly_exponents[i]
            @nextract $N v d -> begin
                for j ∈ length(v_{d-1}):ex[d]
                    push!(v_{d-1}, poly(Val{$N - d}(), n, S))
                end
                v_{d-1}[ ex[d] + 1 ]
            end
            (@nref $ns $v_last d -> ex[d + $N] + 1) = v_coefs[i]
        end
        v_0
    end
end
#const PolyDist{1,T} = Vector{T}
#const PolyDist{2,T} = Vector{Vector{T}}
struct PolyDist{N,T,V} <: AbstractArray{N,T}
    v::V
end
function PolyDist(β::Vector{T}, poly::Vector{NTuple{N,I}}) where {T,N,I<:Integer}

end

@generated function mss4(x, y, ::Val{S}, ::Val{N}) where {S,N}
    quote
        T = eltype(x)
        out = zero(T)
        @nexprs $N i -> @inbounds out += sum( vload( Vec{$S,T}, x, $S*(i-1)+1 ) * vload( Vec{$S,T}, y, $S*(i-1)+1 ) )
        out
    end
end
@generated function mssuneven(x,y,::Val{S},::Val{N},::Val{R}) where {S,N,R}
    rn = 1 + S*N
    quote
        T = eltype(x)
        out = zero(T)
        @nexprs $N i -> @inbounds out += sum( vload( Vec{$S,T}, x, $S*(i-1)+1 ) * vload( Vec{$S,T}, y, $S*(i-1)+1 ) )
        @inbounds out += sum( vload( Vec{$R,T}, x, $rn ) * vload( Vec{$R,T}, y, $rn ) )
        out
    end
end
@generated function mss(x, y, ::Val{N}, ::Val{S}) where {S,N}
    n, r = divrem(N,S)
    r == 0 ? :(mss4(x, y, Val{$S}(), Val{$n}())) : :(mssuneven(x,y, Val{$S}(), Val{$n}(), Val{$r}()))
end

Base.length(x::PolyDist) = length(x.v)
@generated function Base.size(x_0::PolyDist{N}) where N
    quote
        @nextract $N x i -> x_{i-1}[1]
        @ntuple $N i -> length(x_{i-1})
    end
end
Base.getindex(x::PolyDist, i) = x.v[i]
@generated function Base.getindex(x_0::PolyDist{N}, I::Vararg{Int, N}) where N
    quote
        @nextract $N x i -> x_{i-1}[I[i]]
    end
end
@generated function Base.setindex!(x_0::PolyDist{N}, v, I::Vararg{Int, N}) where N
    Nm1 = N - 1
    x_last = Symbol("x_", Nm1)
    quote
        @nextract $Nm1 x i -> x_{i-1}[I[i]]
        $x_last[I[$N]] = v
    end

end
function Base.setindex!(x_0::PolyDist{N,T,V}, v::Vs, i::I) where {N,T,Vs,V <: AbstractArray{Vs},I <: Integer}
    x_0.v[i] = v
end

@generated function marginalize_last(v_0::PolyDist{N,T}, s::NTuple{ns,Val{I}} = (Val{32}(),), nextshape = (32,)) where {N,T,ns,I}
    Nm1 = N - ns
    Nm2 = Nm1 - 1
    v_last = Symbol("v_", Nm1)
    m_out = Symbol("m_", Nm1)
    quote
        $m_out = poly(Val{$Nm1}(), nextshape )
        @nloops $Nm1 i d -> begin
            1:length(v_{$Nm1-d})
        end d -> begin
            d > 1 && begin
                m_{d-1} = poly(Val{d-1}(), nextshape )
                push!(m_d, m_{d-1})
            end
            v_{$N-d} = v_{$Nm1-d}[i_{d}]
        end begin
            m_1[i_1] = marginalize($v_last, s )
        end
        $m_out
    end
end
@generated function marginalize_last2(v_0::PolyDist{N,T}, s::NTuple{ns,Val{I}} = (Val{32}(),), nextshape = (32,)) where {N,T,ns,I}
    Nm1 = N - ns
    v_last = Symbol("v_", Nm1)
    m_out = Symbol("m_", Nm1)
    quote
        local $m_out
        @nloops $Nm1 i d -> begin
            1:length(v_{$Nm1-d})
        end d -> begin
            m_d = d == $Nm1 ? poly(Val{$Nm1}(), nextshape ) : m_{1+d}[i_{d}]
            v_{$N-d} = v_{$Nm1-d}[i_{d}]
        end begin
            m_1[i_1] = marginalize($v_last, s )
        end
        $m_out
    end
end

@generated function receive_and_marginalize_last(v_0::PolyDist{N,T}, message, s = (Val{32}(),) ) where {N,T}
    Nm1 = N - length(s)
    v_last = Symbol("v_", Nm1)
    m_out = Symbol("m_", Nm1)
    plan = plan_rfft(message)
    inv_plan = inv(plan)
    mess_fft = plan * message
    v_fft = similar(mess_fft)
    quote
        @nloops $Nm1 i d -> begin
            1:length(v_{$Nm1-d})
        end d -> begin
            m_{d} = d == $Nm1 ? poly(Val{$Nm1}()) : m_{1+d}[i_{d}]
            v_{$N-d} = v_{$Nm1-d}[i_{d}]
        end begin
            A_mul_B!($v_fft, $plan, $v_last)
            $v_fft .*= $mess_fft
            A_mul_B!($v_last, $inv_plan, $v_fft)
            $m_1[i_1] = marginalize($v_last, $(s...) )
        end
        $m_out
    end
end

@generated function max_degree_poly(::Val{N}, l::Int, T = Int) where N
    quote
        j_0 = l
        s_0 = 0
        out = Vector{NTuple{N,T}}(0)
        @nloops $N i p -> begin
            0:j_{$N-p}
        end p -> begin
            s_{$N-p+1} = s_{$N-p} + i_p
            j_{$N-p+1} = l - s_{$N-p+1}
        end begin
            push!(out, (@ntuple $N j -> T(i_{j})) )
        end
        out
    end
end
function tupGridd(x::NTuple{N,T}) where {N,T}
    ntuple( i -> begin
        x[i] > 0 ? :( $(Symbol('a' + i - 1, "_", x[i])) - $(Symbol('a' + i - 1, "_", x[i]-1)) ) : Symbol('a' + i - 1, "_", x[i])
     end, Val{N}())
end

@generated function lnorm_poly(::Val{N}, d::Int, l=0.5, T = Int) where N
    quote
        j_0 = d
        s_0 = 0
        reciprical = inv(l)
        dl = d^l
        out = Vector{NTuple{N,T}}(0)
        @nloops $N i p -> begin
            0:j_{$N-p}
        end p -> begin
            s_{$N-p+1} = s_{$N-p} + i_p^l
            j_{$N-p+1} = floor(Int, (dl - s_{$N-p+1})^reciprical)
        end begin
            push!(out, (@ntuple $N j -> T(i_{j})) )
        end
        out
    end
end

@generated function polycount(::Val{N}, l::Int) where N
    seq = ( 1:1, 1:3, 1:9, 1:14 )
    step = ( 1:1, 2:3, 4:9, 10:14 )
    quote
        l == 1 && return 1
        j_0 = l
        l_0 = 1
        s_0 = 0
        out = 0
        @nloops $N i p -> begin
            l_{$N-p}:j_{$N-p}
        end p -> begin
            s_{$N-p+1} = s_{$N-p} + i_p - 1
            j_{$N-p+1} = l - s_{$N-p+1}
            l_{$N-p+1} = max(1, (3-p)*l - s_{$N-p+1} )
        end begin
            @nloops $N k p -> begin
                p == $N ? $seq[i_p] : $step[i_p]
            end begin
                out += 1
            end
        end
        out
    end
end
@generated function polycount2(::Val{N}, seq = [ -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 ], l = length(seq)-1) where N
    quote
        j_0 = l
        l_0 = 1
        s_0 = 0
        out = 0
        @nloops $N i p -> begin
            l_{$N-p}:j_{$N-p}
        end p -> begin
            s_{$N-p+1} = s_{$N-p} + i_p - 1
            j_{$N-p+1} = l - s_{$N-p+1}
            l_{$N-p+1} = max(1, (3-p)*l - s_{$N-p+1} )
        end begin
            @nloops $N k p -> begin
                if p == $N
                    0:seq[i_p+1]
                else
                    1+seq[i_p]:seq[i_p+1]
                end
            end begin
                out += 1
            end
        end
        out
    end
end
@generated function polycount2(seq::NTuple{N,Vector{T}} = ([ -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 ],[ -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 ]), l = length(seq[1])-1) where {N,T <: Integer}
    quote
        j_0 = l
        l_0 = 1
        s_0 = 0
        out = 0
        @nloops $N i p -> begin
            l_{$N-p}:j_{$N-p}
        end p -> begin
            s_{$N-p+1} = s_{$N-p} + i_p - 1
            j_{$N-p+1} = l - s_{$N-p+1}
            l_{$N-p+1} = max(1, (3-p)*l - s_{$N-p+1} )
        end begin
            @nloops $N k p -> begin
                if p == $N
                    0:seq[p][i_p+1]
                else
                    1+seq[p][i_p]:seq[p][i_p+1]
                end
            end begin
                out += 1
            end
        end
        out
    end
end

@generated function pcvec(::Val{N}, seq = [ -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 ], l = length(seq)-2) where N
    quote
        l += 1
        j_0 = l
        l_0 = 1
        s_0 = 0
        out = Vector{NTuple{$N,Cuint}}(0)
        @nloops $N i p -> begin
            l_{$N-p}:j_{$N-p}
        end p -> begin
            s_{$N-p+1} = s_{$N-p} + i_p - 1
            j_{$N-p+1} = l - s_{$N-p+1}
            l_{$N-p+1} = max(1, (3-p)*l - s_{$N-p+1} )
        end begin
            @nloops $N k p -> begin
                if p == $N
                    0:seq[i_p+1]
                else
                    1+seq[i_p]:seq[i_p+1]
                end
            end begin
                push!(out, (@ntuple $N k -> Cuint(k)) )
            end
        end
        out
    end
end

([-1,0,1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],
[-1,0,1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18])
function poly_basis_gen(::Val{2})
    p = pcvec2(Val{2}(), (-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), 192)
    p, coef(p)
end
function poly_basis_gen(::Val{3})
    p = pcvec2(Val{2}(), ( -1, 0, 1, 2, 3, 4, 5, 6,  8, 10,  14 ), 256)
    p, coef(p)
end
function poly_basis_gen(::Val{4})
    p = pcvec2(Val{3}(), ( -1, 0, 1, 3,4, 10, 15 ), 320)
    p, coef(p)
end
function poly_basis_gen(::Val{5})
    p = pcvec2((
        [ -1, 0, 1, 3,11, 14 ],
        [ -1, 0, 1, 3,10, 15 ],
        [ -1, 0, 1, 3,10, 15 ],
        [ -1, 0, 1, 3,10, 15 ],
        [ -1, 0, 1, 3,10, 15]), 384)
    p, coef(p)
end
function poly_basis_gen(::Val{6})
    p = pcvec2((
        [ -1, 0, 1, 3,5, 15 ],
        [ -1, 0, 1, 3,5, 15 ],
        [ -1, 0, 1, 3,5, 15 ],
        [ -1, 0, 1, 3,5, 14 ],
        [ -1, 0, 1, 3,6, 14 ],
        [ -1, 0, 1, 3,6, 14 ]), 448)
    p, coef(p)
end
function poly_basis_gen(::Val{7})
    p = pcvec2((
        [ -1, 0, 1, 3,14, 14 ],
        [ -1, 0, 1, 3,14, 15 ],
        [ -1, 0, 1, 3,14, 15 ],
        [ -1, 0, 1, 3,14, 15 ],
        [ -1, 0, 1, 3,14, 15 ],
        [ -1, 0, 1, 3,13, 15 ],
        [ -1, 0, 1, 3,13, 15 ]), 1024)
    p, coef(p)
end
function poly_basis_gen(::Val{8})
    p = pcvec2((
        [ -1, 0, 1, 3,9, 15 ],
        [ -1, 0, 1, 3,9, 15 ],
        [ -1, 0, 1, 3,8, 15 ],
        [ -1, 0, 1, 3,8, 14 ],
        [ -1, 0, 1, 3,9, 14 ],
        [ -1, 0, 1, 3,9, 14 ],
        [ -1, 0, 1, 3,9, 14 ],
        [ -1, 0, 1, 3,9, 14 ]), 1052)
    p, coef(p)
end
function poly_basis_gen(::Val{9})
    p = pcvec2((
        [ -1, 0, 1, 3,4, 15 ],
        [ -1, 0, 1, 3,4, 15 ],
        [ -1, 0, 1, 3,4, 15 ],
        [ -1, 0, 1, 3,4, 14 ],
        [ -1, 0, 1, 3,5, 14 ],
        [ -1, 0, 1, 3,5, 14 ],
        [ -1, 0, 1, 3,5, 14 ],
        [ -1, 0, 1, 3,5, 14 ],
        [ -1, 0, 1, 3,5, 14 ]), 1280)
    p, coef(p)
end

function poly_basis_gen(::Val{10})
    p = pcvec2((
        [ -1, 0, 1, 3, 12, 15 ],
        [ -1, 0, 1, 3, 12, 15 ],
        [ -1, 0, 1, 3, 12, 15 ],
        [ -1, 0, 1, 2, 12, 15 ],
        [ -1, 0, 1, 2, 13, 15 ],
        [ -1, 0, 1, 2, 13, 15 ],
        [ -1, 0, 1, 2, 13, 15 ],
        [ -1, 0, 1, 2, 13, 15 ],
        [ -1, 0, 1, 2, 13, 15 ],
        [ -1, 0, 1, 2, 13, 14 ]), 2112)
    p, coef(p)
end
function poly_basis_gen(::Val{11})
    p = pcvec2((
        [ -1, 0, 1, 2, 11, 14 ],
        [ -1, 0, 1, 2, 11, 14 ],
        [ -1, 0, 1, 2, 11, 15 ],
        [ -1, 0, 1, 2, 11, 15 ],
        [ -1, 0, 1, 2, 11, 15 ],
        [ -1, 0, 1, 2, 10, 15 ],
        [ -1, 0, 1, 2, 10, 15 ],
        [ -1, 0, 1, 2, 10, 15 ],
        [ -1, 0, 1, 2, 10, 15 ],
        [ -1, 0, 1, 2, 10, 15 ],
        [ -1, 0, 1, 2, 10, 15 ]), 2304)
    p, coef(p)
end
function poly_basis_gen(::Val{12})
    p = pcvec2((
        [ -1, 0, 1, 2, 8, 14 ],
        [ -1, 0, 1, 2, 8, 14 ],
        [ -1, 0, 1, 2, 7, 15 ],
        [ -1, 0, 1, 2, 7, 15 ],
        [ -1, 0, 1, 2, 7, 15 ],
        [ -1, 0, 1, 2, 7, 15 ],
        [ -1, 0, 1, 2, 7, 15 ],
        [ -1, 0, 1, 2, 7, 15 ],
        [ -1, 0, 1, 2, 7, 14 ],
        [ -1, 0, 1, 2, 7, 14 ],
        [ -1, 0, 1, 2, 7, 14 ],
        [ -1, 0, 1, 2, 7, 14 ]), 2496)
    p, coef(p)
end
function poly_basis_gen(::Val{13})
    p = pcvec2((
        [ -1, 0, 1, 2, 9, 15 ],
        [ -1, 0, 1, 2, 9, 15 ],
        [ -1, 0, 1, 2, 10, 15 ],
        [ -1, 0, 1, 2, 10, 15 ],
        [ -1, 0, 1, 2, 10, 15 ],
        [ -1, 0, 1, 2, 10, 15 ],
        [ -1, 0, 1, 2, 10, 14 ],
        [ -1, 0, 1, 2, 10, 14 ],
        [ -1, 0, 1, 2, 10, 14 ],
        [ -1, 0, 1, 2, 10, 14 ],
        [ -1, 0, 1, 2, 10, 14 ],
        [ -1, 0, 1, 2, 10, 14 ],
        [ -1, 0, 1, 2, 10, 14 ]), 3584)
    p, coef(p)
end
function poly_basis_gen(::Val{14})
    p = pcvec2((
        [ -1, 0, 1, 2, 6, 15 ],
        [ -1, 0, 1, 2, 6, 15 ],
        [ -1, 0, 1, 2, 6, 15 ],
        [ -1, 0, 1, 2, 6, 14 ],
        [ -1, 0, 1, 2, 6, 14 ],
        [ -1, 0, 1, 2, 6, 14 ],
        [ -1, 0, 1, 2, 6, 14 ],
        [ -1, 0, 1, 2, 7, 14 ],
        [ -1, 0, 1, 2, 7, 14 ],
        [ -1, 0, 1, 2, 7, 14 ],
        [ -1, 0, 1, 2, 7, 14 ],
        [ -1, 0, 1, 2, 7, 14 ],
        [ -1, 0, 1, 2, 7, 14 ],
        [ -1, 0, 1, 2, 7, 14 ]), 3840)
    p, coef(p)
end
function poly_basis_gen(::Val{15})
    p = pcvec2((
        [ -1, 0, 1, 2, 4, 14 ],
        [ -1, 0, 1, 2, 4, 14 ],
        [ -1, 0, 1, 2, 4, 14 ],
        [ -1, 0, 1, 2, 4, 14 ],
        [ -1, 0, 1, 2, 4, 14 ],
        [ -1, 0, 1, 2, 3, 14 ],
        [ -1, 0, 1, 2, 3, 14 ],
        [ -1, 0, 1, 2, 3, 14 ],
        [ -1, 0, 1, 2, 3, 14 ],
        [ -1, 0, 1, 2, 3, 14 ],
        [ -1, 0, 1, 2, 3, 14 ],
        [ -1, 0, 1, 2, 3, 14 ],
        [ -1, 0, 1, 2, 3, 14 ],
        [ -1, 0, 1, 2, 3, 14 ],
        [ -1, 0, 1, 2, 3, 14 ]), 4096)
    p, coef(p)
end
poly_basis_gen(::Val{N}) where N = throw("Size $N not yet implemented.") 
@generated function pcvec2(::Val{N}, seq, total_size = (N+1)*64, ::Type{I} = Int) where {N,I<:Integer}
    quote
        l = length(seq) - 1
        j_0 = l
        l_0 = 1
        s_0 = 0
        out = Vector{NTuple{$N,I}}(undef, 0)
        sizehint!(out, total_size)
        @nloops $N i p -> begin
            l_{$N-p}:j_{$N-p}
        end p -> begin
            s_{$N-p+1} = s_{$N-p} + i_p - 1
            j_{$N-p+1} = l - s_{$N-p+1}
            l_{$N-p+1} = max(1, (3-p)*l - s_{$N-p+1} )
        end begin
            @nloops $N k p -> begin
                if p == $N
                    0:seq[i_p+1]
                else
                    1+seq[i_p]:seq[i_p+1]
                end
            end begin
                push!(out, (@ntuple $N p -> I(k_p)) )
            end
        end
        out
    end
end
@generated function pcvec2(seq::NTuple{N,Vector{T}}, total_size = (N+1)*64) where {N, T<: Integer}
    quote
        l = length(seq[1]) - 1
        j_0 = l
        l_0 = 1
        s_0 = 0
        out = Vector{NTuple{$N,T}}(undef, 0)
        sizehint!(out, total_size)
        @nloops $N i p -> begin
            l_{$N-p}:j_{$N-p}
        end p -> begin
            s_{$N-p+1} = s_{$N-p} + i_p - 1
            j_{$N-p+1} = l - s_{$N-p+1}
            l_{$N-p+1} = max(1, (3-p)*l - s_{$N-p+1} )
        end begin
            @nloops $N k p -> begin
                if p == $N
                    0:seq[p][i_p+1]
                else
                    1+seq[p][i_p]:seq[p][i_p+1]
                end
            end begin
                push!(out, (@ntuple $N k -> T(k)) )
            end
        end
        out
    end
end



@generated function pcvec3(::Val{N}) where N
    seq = ( (-1,0,1,2,3,4,5,6,7,8,9,10,11,12,12,13,13,13,13,13,13,14,14,15), )[N-1]
    l = length(seq) - 2
    quote
        l += 1
        j_0 = l
        l_0 = 1
        s_0 = 0
        out = Vector{NTuple{$N,Int}}(undef, 0)
        @nloops $N i p -> begin
            l_{$N-p}:j_{$N-p}
        end p -> begin
            s_{$N-p+1} = s_{$N-p} + i_p - 1
            j_{$N-p+1} = l - s_{$N-p+1}
            l_{$N-p+1} = max(1, (3-p)*l - s_{$N-p+1} )
        end begin
            @nloops $N k p -> begin
                if p == $N
                    0:seq[i_p+1]
                else
                    1+seq[i_p]:seq[i_p+1]
                end
            end begin
                push!(out, (@ntuple $N k) )
            end
        end
        out
    end
end
@generated function smolcount(::Val{N}, l::Int, λ = 1) where N
    quote
        j_0 = l
        lλ = l^λ
        λi = 1/λ
        l_0 = 1
        s_0 = 0
        @nloops $N i p -> begin
            l_{$N-p}:j_{$N-p}
        end p -> begin
            s_{$N-p+1} = s_{$N-p} + i_p^λ - 1
            j_{$N-p+1} = round(Int, (lλ - s_{$N-p+1})^λi)
            l_{$N-p+1} = max(1, (3-p)*l - s_{$N-p+1} )
        end begin
            println( ( @ntuple $N i ) )
        end
    end
end
