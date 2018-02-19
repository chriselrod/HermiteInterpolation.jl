
struct Marginalize{T} end
Base.@pure Marginalize(marg...) = Marginalize{marg}()

struct Marginal{T} end
Base.@pure Marginal(marg...) = Marginal{marg}()

struct SparseGrid{T,N,L} <: AbstractArray{T,A<:AbstractArray{T},N,L}

    
    nodes::A
end

@generated function SparseGrid(::Type{T}, ::Val{N}, ::Val{L}) where {T <: Real, N, L}
    SparseGrid{T,N,L}(sparse_grid_mat(T, Val{N}(), l))
end
@generated function SparseGrid(::Type{A}, ::Val{N}, ::Val{L}) where {T, A <: AbstractArray{T}, N, L}
    SparseGrid{T,N,L}(A(sparse_grid_mat(T, Val{N}(), l)))
end
@generated function sized_setdiff(x::NTuple{Nold,I}, y::NTuple{Nnew,I}) where {Nold,Nnew,I<:Integer}
    Ndiff = Nold - Nnew
    quote
        @nextract $Nold o x
        @nextract $Nnew n y

        @nexprs $Ndiff i -> begin
            d_i = 1
        end
    end
end

@generated function SparseGrid(sg::SparseGrid{T,Nold,L}, ::Marginal{marg}) where {T,Nold,L,marg}
    Nnew = length(marg)
    Ndrop = Nold - Nnew
    quote
        increment = 0
        @nexprs $Nd i -> begin
            while (i + increment) in retain
                increment += 1
            end
            drop_i = i + increment
        end
        drop = @ntuple $Nd i -> drop_i

        out = Matrix{T}($Nn, polycount(Val{Nn}(), L))


    end
end
@generated function SparseGrid(sg::SparseGrid{T,No,L}, retain::NTuple{Nn,I}) where {T,No,L,Nn,I<:Integer}
    Nd = No - Nn
    quote
        increment = 0
        @nexprs $Nd i -> begin
            while (i + increment) in retain
                increment += 1
            end
            drop_i = i + increment
        end
        drop = @ntuple $Nd i -> drop_i

        out = Matrix{T}($Nn, polycount(Val{Nn}(), L))


    end
end

sparse_grid_mat(::Type{T}, ::Val{N}, l::Int) = sparse_grid_mat!(Matrix{T}($N, polycount(Val{N}(), l)), Val{N}(), l)

@generated function sparse_grid_mat!(out::Matrix{T}, ::Val{N}, l::Int) where {T,N}
    quote
        j_0 = l
        l_0 = 1
        s_0 = 0
        ind = 0
        @nloops $N i p -> begin
            l_{$N-p}:j_{$N-p}
        end p -> begin
            s_{$N-p+1} = s_{$N-p} + i_p - 1
            j_{$N-p+1} = l - s_{$N-p+1}
            l_{$N-p+1} = max(1, (3-p)*l - s_{$N-p+1} )
        end begin
            @nloops $N j p -> begin
                p == $N ? $hermite_nodes[i_p] : $hermite_nodes_diffs[i_p]
            end begin
                ind += 1
                @nexprs $N k -> begin
                    out[k,ind] = j_k
                end

            end
        end
        out
    end
end

#Ridiculously ad hoc.
#It is possible to skip a LOT of these iterations.
#Likely also to not have to iterate at all.
#This is absolutely terrible.
@generated function sg_sub2ind(sg::SparseGrid{T,N,L}, ind::NTuple{N,I}) where {N,T,L,I <: Integer}
    Np1 = N+1
    quote
        j_0 = $L
        l_0 = 1
        s_0 = 0
        out = 0
        @nloops $N i p -> begin
            l_{$N-p}:j_{$N-p}
        end p -> begin
            s_{$N-p+1} = s_{$N-p} + i_p - 1
            j_{$N-p+1} = $L - s_{$N-p+1}
            l_{$N-p+1} = max(1, (3-p)*$L - s_{$N-p+1} )
        end begin
            @nloops $N k p -> begin
                max_l = length(hermite_nodes[i_p])
                if $No == 1
                    min_l = i_1 == 1 ? 1 : length(hermite_nodes[i_p-1])+1
                elseif p == $N
                    min_l = 1
                elseif i_p == 1
                    min_l = 1
                else
                    min_l = length(hermite_nodes[i_p-1])+1
                end
                min_l:max_l
            end begin
                out += 1
                @nif $Np1 d -> k_d != ind[d] d -> () d -> return out
            end
        end
        out
    end
end

@generated function sg_sub2ind(::Val{L}, ind::NTuple{N, I}) where {L,N,I<:Integer}
    Np1 = N+1
    quote
        j_0 = $L
        l_0 = 1
        s_0 = 0
        out = 0
        @nloops $N i p -> begin
            l_{$N-p}:j_{$N-p}
        end p -> begin
            s_{$N-p+1} = s_{$N-p} + i_p - 1
            j_{$N-p+1} = $L - s_{$N-p+1}
            l_{$N-p+1} = max(1, (3-p)*$L - s_{$N-p+1} )
        end begin
            @nloops $N k p -> begin
                max_l = length(hermite_nodes[i_p])
                if $No == 1
                    min_l = i_1 == 1 ? 1 : length(hermite_nodes[i_p-1])+1
                elseif p == $N
                    min_l = 1
                elseif i_p == 1
                    min_l = 1
                else
                    min_l = length(hermite_nodes[i_p-1])+1
                end
                min_l:max_l
            end begin
                out += 1
                @nif $Np1 d -> k_d != ind[d] d -> () d -> return out
            end
        end
        throw("Out of bounds.")
    end
end

function gen_sub(marg::NTuple{N,I}) where {N,I<:Integer}
    ntuple(i -> Symbol() Val{N}())
end
function symbolify(x::NTuple{N,T}, sym::Symbol = :j_) where {N,T}
    ntuple( p -> Symbol(sym, x[p]) , Val{N}())
end
@generated function tm(::Marginal{marg}) where marg
    Nin = length(marg) #@ntuple $Nin p -> j_{ $marg[p] }
    tnin = 2Nin
    tupj = Expr(:tuple, symbolify(marg)... )
    quote
        @nexprs $tnin p -> begin
            j_p = 4p
        end
        $tupj
    end
end

#@generated function _marginalize(sg::SparseGrid{T,No,L}, marg::NTuple{Nin}, margout::NTuple{Nout}) where {T,No,L, Nin, Nout}
#end
@generated function marginalize(sg::SparseGrid{T,No,L}, ::Marginal{marg}, ::Marginalize{margout}) where {T,No,L, marg, margout}
    _marginalize(sg::SparseGrid{T,No,L}, Marginal{marg}(), Marginalize{margout}())
end

## Need to replace Marginal{marg} and Marginalize{margout} with marg::NTuple and margout::NTuple
## add out and in_interp
@generated function _marginalize(sg::SparseGrid{T,No,L}, ::Marginal{marg}, ::Marginalize{margout}) where {T,No,L, marg, margout}
    Nin = length(marg) #@ntuple $Nin p -> j_{ $marg[p] }
    Nout = length(margout)
    tupj_marg = Expr(:tuple, symbolify(marg)... ) #For sub2ind-ing
    tupj_margout = Expr(:tuple, symbolify(margout)... ) #For marginalizing out
    tupj_margout_course = Expr(:tuple, symbolify(margout, :i_)...)
    calc_l = :(sum($tupj_margout_course) - $(length(margout)))
    quote
        ex = quote @fastmath begin end end
        exa = ex.args[2].args[2].args


        j_0 = $L
        l_0 = 1
        s_0 = 0
        ind = 0
        @nloops $No i p -> begin
            l_{$N-p}:j_{$N-p}
        end p -> begin
            s_{$N-p+1} = s_{$N-p} + i_p - 1
            j_{$N-p+1} = $L - s_{$N-p+1}
            l_{$N-p+1} = max(1, (3-p)*$L - s_{$N-p+1} )
        end begin
            @nloops $N j p -> begin
#                max_l = p == $N ? length(hermite_nodes[i_p]) : length(hermite_nodes_diffs[i_p])
#                1:max_l

                max_l = length(hermite_nodes[i_p])
                if $No == 1
                    min_l = i_1 == 1 ? 1 : length(hermite_nodes[i_p-1])+1
                elseif p == $N
                    min_l = 1
                elseif i_p == 1
                    min_l = 1
                else
                    min_l = length(hermite_nodes[i_p-1])+1
                end
                min_l:max_l
            end begin
                ind += 1
           #     @nexprs $N k -> begin

                l = $L + $No - sum( ( @ntuple $No i ) ) 

                marg_ind = sg_sub2ind( Val{L}(), $tupj_marg )
                weight_val = calc_weight( Val{L}(), l, $tupj_margout, $tupj_margout_course)
                expr = :( temp_out = in_interp[1,$ind] * $weight_val )
                push!(exa, expr)
                for j ∈ 1:$Nout #jth derivative
                    deriv_ind = $margout[j]+1
                    weight = calc_weight( Val{L}(), l, $tupj_margout, $tupj_margout_course, j)
                    expr = :( temp_out += in_interp[$deriv_ind,$ind] * $weight )
                    push!(exa, expr)
                end
                push!(exa, :(out[1,$marg_ind] += temp_out))
                for j ∈ 2:$Nin+1
                    deriv_ind = $marg[j-1] + 1
                    push!(exa, :(out[$j,$marg_ind] += in_interp[$deriv_ind,$ind] * $weight_val ))
                end

         #           out[k,ind] = j_k
         #       end

            end
        end
        ex


    end
end

macro nnexprs(N, itersym, rangeexpr, args...)
    _nnexprs(N, itersym, rangeexpr, args...)
end
import Base.Cartesian: inlineanonymous
function _nnexprs(N, itersym, rangeexpr, args...)

#    body = args[end]
#    ex = Expr(:escape, body)
    ex = args
    for dim ∈ 1:N
        itervar = inlineanonymous(itersym, dim) #for :i, dim produces i_dim
        rng = inlineanonymous(rangeexpr, dim) #subs the dim appropriately into expression and simplifies
        ex = quote
            @nexprs $(esc(rng)) $(esc(itervar)) -> $ex
        end

    end
    ex
end
nnt = quote
        @nexprs 4 i_4 -> begin
            @nexprs 5 i_3 -> begin
                @nexprs 6 i_2 -> begin
                    @nexprs 7 i_1 -> begin

                        @ntuple 4 i
                    end
                end
            end
        end
    end

    
function Base.start(sg::SparseGrid{T,N}) where {T,N}

end
function Base.next(sg::SparseGrid{T,N}, state) where {T,N}

end
function Base.done(sg::SparseGrid{T,N}, state) where {T,N}

end
Base.IteratorSize(sg::SparseGrid) = HaseLength()
Base.IteratorEltype(sg::SparseGrid{T,N}) = HasEltype()
Base.eltype(sg::SparseGrid{T,N}) = Tuple{Tuple{NTuple{N,T},T},SparseGrid{T,N}}

@generated function Base.size(sg::SparseGrid{T,N,L}) where {T,N,L}
    ntuple( length(hermite_nodes) , Val{L}())
end
@generated function Base.length(sg::SparseGrid{T,N,L}) where {T,N,L}
    polycount(Val{N}(), L)
end


@inline function Base.getindex(sg::SparseGrid{T,N}, i::Int) where {T,N}
    ntuple( j -> sg.nodes[j,i], Val{N}())
end
@inline function Base.getindex(sg::SparseGrid{T,N}, sub::Vararg{I, N}) where {T,N,I<:Integer}
    
end

function Base.setindex!(sg::SparseGrid, v, i)
    
end

function Base.setindex!(sg::SparseGrid{T,N}, v, i::Vararg{I, N}) where {T,N,I<:Integer}
    
end

function Base.firstindex(sg::SparseGrid)

end

function Base.lastindex(sg::SparseGrid)

end

Base.IndexStyle(::SparseGrid) = IndexLinear()



@generated function SparseIterator(f, sg::SparseGrid{T,N}, args) where {T,N}


    quote
        j_0 = l
        s_0 = 0
        @nloops $N i p -> begin
            0:j_{$N-p}
        end p -> begin
            s_{$N-p+1} = s_{$N-p} + i_p
            j_{$N-p+1} = l - s_{$N-p+1}
        end begin
            f(args, (@ntuple $N j -> T(i_{j}))... )
        end
        out

    end
end


@generated function polycount(::Val{N}, l::Int) where N
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
                p == $N ? $hermite_nodes[i_p] : $hermite_nodes_diffs[i_p]
            end begin
                out += 1
            end
        end
        out
    end
end
@generated function pc2(::Val{N}, l::Int) where N
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
                max_l = length(hermite_nodes[i_p])
                if $N == 1
                    min_l = i_1 == 1 ? 1 : length(hermite_nodes[i_p-1])+1
                elseif p == $N
                    min_l = 1
                elseif i_p == 1
                    min_l = 1
                else
                    min_l = length(hermite_nodes[i_p-1])+1
                end
                min_l:max_l
            end begin
                location = @ntuple $N d -> k_d
                @show  location
                out += 1
            end
        end
        out
    end
end