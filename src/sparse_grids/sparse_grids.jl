


struct SparseGrid{T,N} <: AbstractArray{T,N}

    
    location::Vector{Int}
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

function Base.size(sg::SparseGrid{T,N}) where {T,N}

end
function Basel.length(sg::SparseGrid{T,N}) where {T,N}

end


function Base.getindex(sg::SparseGrid{T,N}, i::Int) where {T,N}

end

function Base.getindex(sg::SparseGrid{T,N}, i::Vararg{I, N}) where {T,N,I<:Integer}

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