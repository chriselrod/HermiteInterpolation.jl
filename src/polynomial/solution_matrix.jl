

@generated function gen_polies(c::Float32, x::NTuple{N,T}) where {N,T,T2}
    Np1 = N+1
    quote
        out = Vector{Tuple{T2,NTuple{$N,T}}}($Np1)
        @nexprs $N (i -> out[i] = (T2(c*x[i]), ( @ntuple $N j -> i == j ? x[j] - one(T) : x[j] )) )
        out[$Np1] = ( T2(c), x )
        out
    end
end
@generated function gen_polies(x::NTuple{N,T}, ::Val{k}, c::T2 = one(Float32)) where {N,T,T2,k}
    N_level = binomial(N+k-1,N-1)
    quote
        buffer = Vector{T}(undef, $N)
        for i ∈ 1:$N
            buffer[i] = x[i]
        end
        N_0 = $N
        out = Vector{Tuple{T2,NTuple{$N,T}}}(undef, $N_level)
        itert = 0
        @nloops $k i l -> begin
            1:N_{$k-l}
        end l -> begin
            N_{$k-l+1} = i_l
            c_l = buffer[i_l]
            buffer[i_l] = max(zero(T), c_l - one(T))
        end l -> begin
        buffer[i_l] = c_l
        end begin
            itert += 1
            c_temp = c
            @nexprs $k j -> begin
                c_temp *= c_j
            end
            out[itert] = ( c_temp, ( @ntuple $N j -> buffer[j] )) 
        end
        out
    end
end

@generated function eval_poly_term(coef::T, node::NTuple{N,T}, exponent) where {N,T}
    quote
        iszero(coef) && return zero(T)
        @nexprs $N i -> coef *= node[i] ^ exponent[i]
        coef
    end
end


function create_linear_map( seq::Vector{NTuple{N,T}}, poly::Vector{NTuple{N,I}}, coefs::AbstractArray  ) where {N,T,I <: Integer}
    n = length(poly)
    np = length(seq)
    Np1 = N+1
    X = Matrix{T}(n,n)
    for i ∈ 1:n
        polies = gen_polies(coefs[i], poly[i], T)
        for j ∈ 1:np
            node = seq[j]
            l = Np1*(j-1)
            for (k,(c,ex)) ∈ enumerate(polies)
                X[l+k,i] = eval_poly_term(c, node, ex)
            end
        end
    end
    poly, X, inv(X)
end
function create_linear_map( seq::Vector{NTuple{N,T}}, poly::Vector{NTuple{N,I}}  ) where {N,T,I <: Integer}
    n = length(poly)
    np = length(seq)
    Np1 = N+1
    @assert np * Np1 == n
    X = Matrix{T}(n,n)
    for i ∈ 1:n
        polies = gen_polies(one(T), poly[i], T)
        for j ∈ 1:np
            node = seq[j]
            l = Np1*(j-1)
            for (k,(c,ex)) ∈ enumerate(polies)
                X[l+k,i] = eval_poly_term(c, node, ex)
            end
        end
    end
    poly, X, inv(X)
end



function linear_map_stable( seq::Vector{NTuple{N,T}}, poly::Vector{NTuple{N,I}} ) where {N,T,I <: Integer}
    n = length(poly)
    np = length(seq)
    Np1 = N+1
    X = Matrix{T}(undef, np*Np1, n)
    Xm = fill(zero(T), np*Np1, n)
    for i ∈ 1:np
        polies = gen_polies(poly[i], Val{1}(), one(T))
        for j ∈ 1:np
            node = seq[j]
            l = Np1*(j-1)
            for (k,(c,ex)) ∈ enumerate(polies)
                X[l+k,i] = eval_poly_term(c, node, ex)
            end
            X[l+Np1,i] = eval_poly_term(one(T), node, poly[i])
            if i == j
                Xm[l+Np1,i] = one(T)
            end
        end
    end
    for i ∈ np+1:n
        polies = gen_polies(poly[i], Val{1}(), one(T))#T(3)^( -prod(poly[i]) ) )
        for j ∈ 1:np
            node = seq[j]
            l = Np1*(j-1)
            for (k,(c,ex)) ∈ enumerate(polies)
                Xm[l+k,i] = X[l+k,i] = eval_poly_term(c, node, ex)
            end
            X[l+Np1,i] = eval_poly_term(one(T), node, poly[i])
        end
    end
    XtX = Xm' * X
    poly, X, Xm, XtX, inv(XtX)
end


function linear_map_stable_understandable( seq::Vector{NTuple{N,T}}, poly::Vector{NTuple{N,I}} ) where {N,T,I <: Integer}
    n = length(poly)
    np = length(seq)
    Np1 = N+1
    X = Matrix{T}(undef, np*Np1, n)
    Xm = fill(zero(T), np*Np1, n)
    for i ∈ 1:np
        polies = gen_polies(poly[i], Val{1}(), T(3)^( -prod(poly[i]) ) )
        for j ∈ 1:np
            node = seq[j]
            l = np + N*(j-1)
            for (k,(c,ex)) ∈ enumerate(polies)
                X[l+k,i] = eval_poly_term(c, node, ex)
            end
            X[j,i] = eval_poly_term(one(T), node, poly[i])
        end
        Xm[i,i] = one(T)
    end
    for i ∈ np+1:n
        polies = gen_polies(poly[i], Val{1}(), T(3)^( -prod(poly[i]) ) )
        for j ∈ 1:np
            node = seq[j]
            l = np + N*(j-1)
            for (k,(c,ex)) ∈ enumerate(polies)
                Xm[l+k,i] = X[l+k,i] = eval_poly_term(c, node, ex)
            end
            X[j,i] = eval_poly_term(one(T), node, poly[i])
        end
    end
    XtX = Xm' * X
    poly, X, Xm, XtX, inv(XtX)
end