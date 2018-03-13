
struct Marginalize{T} end
Base.@pure Marginalize(marg...) = Marginalize{marg}()

struct Marginal{T} end
Base.@pure Marginal(marg...) = Marginal{marg}()

struct SparseGrid{T,N,L,A} <: AbstractArray{T,N}
    nodes::A
end
@generated grid_length(::Val{N}, ::Val{L}) where {N,L} = _grid_length(Val{N}(), L)
@generated function _grid_length(::Val{N}, L::Int) where {N}
    quote
        ilim_0 = L
        s_0 = 0
        ind = 0
        out = 0
        @nloops $N i p -> begin
            0:ilim_{$N-p}
        end p -> begin
            s_{$N-p+1} = s_{$N-p} + i_p
            ilim_{$N-p+1} = L - s_{$N-p+1}
        end begin
            @nloops $N j p -> begin
                1:rule_lengths[i_p+1]
            end begin
                out += 1
            end
        end
        out
    end
end

@generated function SparseGrid(f, data, ::Val{N}, ::Val{L}, ::Type{T} = Float64) where {T,N,L}

    quote

        ilim_0 = $L
        s_0 = 0
        ind = 0

        out = Vector{$T}(grid_length(Val{$N}(), Val{$L}()))

        @nloops $N i p -> begin
            0:ilim_{$N-p}
        end p -> begin
            s_{$N-p+1} = s_{$N-p} + i_p
            ilim_{$N-p+1} = $L - s_{$N-p+1}
        end begin

            #Maxmimum length left over to the marg-outs
            # max_l = $L - sum( $tupj_marg_course )
            # l_remaining = max_l - sum( $tupj_margout_course ) 
            l_remaining = $L - sum( ( @ntuple $N i ) )

            @nloops $N j p -> begin
                HermiteQuadratureRules.gk_nodes[i_p+1]
            end begin
                ind += 1
                out[ind] = f(data, (@ntuple $N j)...)
            end


        end

    end
end

@generated function SparseGrid(::Type{T}, ::Val{N}, ::Val{L}) where {T <: Real, N, L}
    SparseGrid{T,N,L}(sparse_grid_mat(T, Val{N}(), L))
end
@generated function SparseGrid(::Type{A}, ::Val{N}, ::Val{L}) where {T, A <: AbstractArray{T}, N, L}
    SparseGrid{T,N,L}(A(sparse_grid_mat(T, Val{N}(), L)))
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

sparse_grid_mat(::Type{T}, ::Val{N}, l::Int) where {T,N} = sparse_grid_mat!(Matrix{T}(N, polycount(Val{N}(), l)), Val{N}(), l)

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

function gen_sub(marg::NTuple{N,I}) where {N,I<:Integer}
    ntuple(i -> Symbol(), Val{N}())
end
function symbolify(x::NTuple{N,T}, sym::Symbol = :j_) where {N,T}
    ntuple( p -> Symbol(sym, x[p]) , Val{N}())
end

#@generated function _marginalize(sg::SparseGrid{T,No,L}, marg::NTuple{Nin}, margout::NTuple{Nout}) where {T,No,L, Nin, Nout}
#end
@generated function marginalize(sg::SparseGrid{T,No,L}, ::Marginal{marg}, ::Marginalize{margout}) where {T,No,L, marg, margout}
    _marginalize(sg::SparseGrid{T,No,L}, marg, margout)
end

function marginal_symbol(course_inds::NTuple{N}, fine_inds::NTuple{N}, symbolstring = "m") where N
    for (i,j) in zip(course_inds, fine_inds)
        symbolstring *= "_" * string(i) * "_" * string(j)
    end
    Symbol(symbolstring)
end

function set_symbol_to_zero(course_inds::NTuple{N}, fine_inds::NTuple{N}, ::Type{T} = Float64, symbolstring = "m") where {N,T}
    symbol = marginal_symbol(course_inds, fine_inds, symbolstring)
    z = zero(T)
    :($symbol = $z)
end


## Need to replace Marginal{marg} and Marginalize{margout} with marg::NTuple and margout::NTuple
## add out and in_interp
# @generated function _marginalize(sg::SparseGrid{T,N,L}, marg::NTuple{Nin,Int}, margout::NTuple{Nout}) where {T,N,L,Nin,Nout}
@generated function _marginalize(sg::SparseGrid{T,N,L}, ::Marginal{marg}, ::Marginalize{margout}) where {T,N,L,marg,margout}
    Nin = length(marg)
    Nout = length(margout)
    # calc_l = :(sum($tupj_margout_course) - $(length(margout)))
    tupj_marg = Expr(:tuple, symbolify(marg, :j_)... ) #This tuple is for the fine indices we keep.
    tupj_margout = Expr(:tuple, symbolify(margout, :j_)... ) #This tuple gives the fine indices we're marginalizing out.
    tupj_marg_course = Expr(:tuple, symbolify(marg, :i_)... ) #This tuple gives the course indices we keep.
    tupj_margout_course = Expr(:tuple, symbolify(margout, :i_)...) #This tuple gives the course indices we're marginalizing out.
    quote
        ex = quote @fastmath begin in_interp = sg.nodes end end
        exa = ex.args[2].args[3].args# 0.7
        # exa = ex.args[2].args[2].args# 0.6

        ilim_0 = $L
        s_0 = 0

        marg_length = 0
        #Here we set all the marginal values to zero, so that later we can "+=" them without worry.
        @nloops $Nin i p -> begin
            0:ilim_{$Nin-p}
        end p -> begin
            s_{$Nin-p+1} = s_{$Nin-p} + i_p
            ilim_{$Nin-p+1} = $L - s_{$Nin-p+1}
        end begin
            @nloops $Nin j p -> begin
                1:rule_lengths[i_p+1]
            end begin
                marg_length += 1
                push!(exa,  set_symbol_to_zero( ( @ntuple $Nin i ) , ( @ntuple $Nin j ), T ) )
            end
        end

        ilim_0 = $L
        # s_0 = 0
        ind = 0

        @nloops $N i p -> begin
            0:ilim_{$N-p}
        end p -> begin
            s_{$N-p+1} = s_{$N-p} + i_p
            ilim_{$N-p+1} = $L - s_{$N-p+1}
        end begin

            #Maxmimum length left over to the marg-outs
            # max_l = $L - sum( $tupj_marg_course )
            # l_remaining = max_l - sum( $tupj_margout_course ) 
            l_remaining = $L - sum( ( @ntuple $N i ) )

            @nloops $N j p -> begin
                1:rule_lengths[i_p+1]
            end begin
                ind += 1
                marg_symbol = marginal_symbol( $tupj_marg_course, $tupj_marg )
                weight_val = calc_weights( l_remaining, $tupj_margout_course, $tupj_margout )
                expr = :( $marg_symbol += in_interp[$ind] * $weight_val )
                push!(exa, expr)
            end
        end

        push!(exa, :(marginal_interp = Vector{$T}(undef, $marg_length)))        
        marg_length = 0
        #Here we fill in the marginal vector.
        @nloops $Nin i p -> begin
            0:ilim_{$Nin-p}
        end p -> begin
            s_{$Nin-p+1} = s_{$Nin-p} + i_p
            ilim_{$Nin-p+1} = $L - s_{$Nin-p+1}
        end begin
            @nloops $Nin j p -> begin
                1:rule_lengths[i_p+1]
            end begin
                marg_length += 1
                marg_symbol = marginal_symbol( ( @ntuple $Nin i ), ( @ntuple $Nin j ) )
                push!(exa,  :(marginal_interp[$marg_length] = $marg_symbol ) )
            end
        end
        Nin = $Nin
        output = :(SparseGrid{$T,$Nin,$L,Vector{$T}}(marginal_interp))
        push!(exa, output )
        ex
    end
end

## Need to replace Marginal{marg} and Marginalize{margout} with marg::NTuple and margout::NTuple
## add out and in_interp
@generated function _marginalize_deriv(sg::SparseGrid{T,No,L}, ::Marginal{marg}, ::Marginalize{margout}) where {T,No,L, marg, margout}
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
Base.IteratorEltype(sg::SparseGrid) = HasEltype()
Base.eltype(sg::SparseGrid{T,N}) where {T,N} = Tuple{Tuple{NTuple{N,T},T},SparseGrid{T,N}}

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