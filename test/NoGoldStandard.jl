


# Distributions does not work yet on Julia 0.7
beta_lpdf(x, alpha, beta) = (alpha-1)*log(x) + (beta-1)*log(1-x)

using StaticArrays
const TypedLength{N,T} = Union{NTuple{N,T}, StaticVector{N,T}}

function generate_call(fname, n)
    out = :( $fname( data, params[1] ) )
    for i ∈ 2:n
        push!(out.args, :(params[$i]))
    end
    out
end
# Expr(:meta, :inline)
@generated function NGS(data::Tuple, params::TypedLength{N,T}) where {N,T}
    generate_call(:NoGoldStandard, N )
end

"""Log Density of the No Gold-Standard Model"""
function NoGoldStandard(priordata::Tuple, S₁::T, S₂::T, C₁::T, C₂::T, pi::Vararg{T,K}) where {T,K}
    n, alpha_pi, beta_pi, alpha_S1, beta_S1, alpha_S2, beta_S2,
        alpha_C1, beta_C1, alpha_C2, beta_C2 = priordata
    target = beta_lpdf(S₁, alpha_S1, beta_S1) + beta_lpdf(S₂, alpha_S2, beta_S2)
    target += beta_lpdf(C₁, alpha_C1, beta_C1) + beta_lpdf(C₂, alpha_C2, beta_C2)
    # target += sum(pdf.(Beta( alpha_pi, beta_pi), pi))
    # target += pdf(Beta( alpha_S1, beta_S1), S₁)
    # target += pdf(Beta( alpha_S2, beta_S2), S₂)
    # target += pdf(Beta( alpha_C1, beta_C1), C₁)
    # target += pdf(Beta( alpha_C2, beta_C2), C₂)
    @inbounds for k in 1:K
        target += beta_lpdf( pi[k], alpha_pi, beta_pi )

        target += n[k,1] * log( pi[k] * (   S₁ *   S₂ ) + (1-pi[k]) * ((1-C₁)*(1-C₂)) )
        target += n[k,2] * log( pi[k] * (   S₁ *(1-S₂)) + (1-pi[k]) * ((1-C₁)*   C₂ ) )
        target += n[k,3] * log( pi[k] * ((1-S₁)*   S₂ ) + (1-pi[k]) * (   C₁ *(1-C₂)) )
        target += n[k,4] * log( pi[k] * ((1-S₁)*(1-S₂)) + (1-pi[k]) * (   C₁ *   C₂ ) )

        target += n[k,5] * log( pi[k]*   S₁  + (1-pi[k])*(1-C₁) )
        target += n[k,6] * log( pi[k]*(1-S₁) + (1-pi[k])*   C₁  )

        target += n[k,7] * log( pi[k]*   S₂  + (1-pi[k])*(1-C₂) )
        target += n[k,8] * log( pi[k]*(1-S₂) + (1-pi[k])*   C₂  )
    end
    target
end

logit(x) = log(x/(1-x))
invlogit(x) = 1/(1+exp(-x))
function dinvlogit(x)
    expnx = exp(-x)
    expnx / abs2(1+expnx)
end


# using LineSearches, Optim, BenchmarkTools
# rosenbrock(x) =  (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
# initial_x = zeros(2);
# hage = BFGS();
# back = BFGS(;linesearch = LineSearches.BackTracking());
# od = TwiceDifferentiable(rosenbrock, initial_x; autodiff = :forward);
# @benchmark optimize($od, $initial_x, $hage)
# @benchmark optimize($od, $initial_x, $back)


# Optim.minimizer(ans)
# td = TwiceDifferentiable(f, initial_x; autodiff = :forward)
# Optim.minimizer(optimize(td, initial_x, Newton()))

using Compat
using Base.Cartesian, Compat.Random#, LinearAlgebra
struct Multinomial{N,T}
    n::Int
    p::NTuple{N,T}
    cumulative_p::NTuple{N,T}
end
function Multinomial(n, p::NTuple{N,T}) where {N,T}
    total = Ref{T}(zero(T)) #Cumsum not defined for ntuples; this gets around closure bug.
    cumulative_p = ntuple( i -> begin
        total[] += p[i]
    end, Val{N}())
    Multinomial(n, p, cumulative_p)
end

function Random.rand(m::Multinomial{N,T}, n::Int...) where {N,T}
    rand!(Array{NTuple{N,Int}}(undef, n...), m)
end
@generated function Random.rand!(out::AbstractArray{NTuple{N,I}}, m::Multinomial{N,T}) where {N,T,I<:Integer}
    quote
        for i ∈ eachindex(out)
            @nexprs $N d -> x_d = zero(I)
            for j ∈ 1:m.n
                r = rand()
                @nif $N d -> r < m.cumulative_p[d] d -> x_d += one(I)
            end
            out[i] = rand(m)
        end
        out
    end
end
@generated function Random.rand(m::Multinomial{N,T}) where {N,T}
    quote
        @nexprs $N d -> x_d = 0
        for j ∈ 1:m.n
            r = rand()
            @nif $N d -> r < m.cumulative_p[d] d -> x_d += 1
        end
        @ntuple $N d -> x_d
    end
end
# struct CorrErrors{p, T<: Real}
#     π::SVector{p,T}
#     S::SVector{2,T}
#     C::SVector{2,T}
#     covsc::SVector{2,T}
# end
# function CorrErrors(π::NTuple{p,T}, S::Tuple{T,T}, C::Tuple{T,T}, covsc::Tuple{T,T} = (0.0,0.0)) where {p,T}
#     CorrErrors{p,T}(SVector(π), SVector(S), SVector(C), SVector(covsc))
# end
struct CorrErrors{P, T<: Real}
    π::NTuple{P,T}
    S::NTuple{2,T}
    C::NTuple{2,T}
    covsc::NTuple{2,T}
end
function CorrErrors(π::NTuple{p,T}, S::Tuple{T,T}, C::Tuple{T,T}, covsc::Tuple{T,T} = (0.0,0.0)) where {p,T}
    CorrErrors{p,T}(π, S, C, covsc)
end
function common_p(Θ::CorrErrors{P,Float64}) where P
  [(Θ.π[i]*(Θ.S[1]*Θ.S[2]+Θ.covsc[1]) + (1-Θ.π[i])*((1-Θ.C[1])*(1-Θ.C[2])+Θ.covsc[2]), Θ.π[i]*(Θ.S[1]*(1-Θ.S[2])-Θ.covsc[1]) + (1-Θ.π[i])*((1-Θ.C[1])*Θ.C[2]-Θ.covsc[2]), Θ.π[i]*(Θ.S[2]*(1-Θ.S[1])-Θ.covsc[1]) + (1-Θ.π[i])*(Θ.C[1]*(1-Θ.C[2])-Θ.covsc[2]), Θ.π[i]*((1-Θ.S[1])*(1-Θ.S[2])+Θ.covsc[1]) + (1-Θ.π[i])*(Θ.C[1]*Θ.C[2]+Θ.covsc[2])) for i ∈ 1:P]
end
function p_i(Θ::CorrErrors{P,Float64}, i::Int) where P
  [(Θ.π[j]*Θ.S[i] + (1-Θ.π[j])*(1-Θ.C[i]), Θ.π[j]*(1-Θ.S[i]) + (1-Θ.π[j])*Θ.C[i]) for j ∈ 1:P]
end
function gen_data(Θ::CorrErrors{P,T}, n_common::Int, n_1_only::Int, n_2_only::Int) where {P,T}
    double_test = common_p(Θ)
    p_1_only = p_i(Θ, 1)
    p_2_only = p_i(Θ, 2)
    out = Matrix{Int}(8, P)
    for i ∈ 1:P
        out[1:4,i] .= rand(Multinomial(n_common, double_test[i]))
        out[5:6,i] .= rand(Multinomial(n_1_only, p_1_only[i]))
        out[7:8,i] .= rand(Multinomial(n_2_only, p_2_only[i]))
    end
    out
end

corr_errors_independent2 = CorrErrors((0.1,        0.4), (0.9, 0.95), (0.85,0.97));
corr_errors_independent3 = CorrErrors((0.1,  0.25, 0.4), (0.9, 0.95), (0.85,0.97));
corr_errors_independent4 = CorrErrors((0.1,0.2,0.3,0.4), (0.9, 0.95), (0.85,0.97));

n_small_2  = gen_data(corr_errors_independent2,  100,  200,  16);
n_medium_2 = gen_data(corr_errors_independent2,  400,  800,  64);
n_big_2    = gen_data(corr_errors_independent2, 1600, 3200, 252);

n_small_3  = gen_data(corr_errors_independent3,  100,  200,  16);
n_medium_3 = gen_data(corr_errors_independent3,  400,  800,  64);
n_big_3    = gen_data(corr_errors_independent3, 1600, 3200, 252);

n_small_4  = gen_data(corr_errors_independent4,  100,  200,  16);
n_medium_4 = gen_data(corr_errors_independent4,  400,  800,  64);
n_big_4    = gen_data(corr_errors_independent4, 1600, 3200, 252);



const cds2 = (n_small_4, 1,1,1,1,1,1,1,1,1,1);
x = ntuple(i -> rand(), Val(8));
using BenchmarkTools
@benchmark NGS($ds2, $x)


# julia> ForwardDiff.gradient(x -> NGS(cds2,x), svx)
# 8-element SArray{Tuple{8},Float64,1,8}:
#   8.671839440018804e12 
#   9.607806107425758e13 
#  -1.514460021428038e14 
#  -4.5768269817453856e14
#  -1.966459100582506e14 
#  -4.929286840442453e13 
#   2.1047957073464297e13
#  -1.0733791260531581e14

# julia> @benchmark ForwardDiff.gradient(x -> NGS(cds2,x), $svx)
# BenchmarkTools.Trial: 
#   memory estimate:  0 bytes
#   allocs estimate:  0
#   --------------
#   minimum time:     945.061 ns (0.00% GC)
#   median time:      949.848 ns (0.00% GC)
#   mean time:        952.880 ns (0.00% GC)
#   maximum time:     1.670 μs (0.00% GC)
#   --------------
#   samples:          10000
#   evals/sample:     33

# gradconf = ForwardDiff.GradientConfig(x -> NGS(cds2, x), svx)


# #Slower?

# julia> @benchmark ForwardDiff.vector_mode_gradient(x -> NGS(cds2,x), $svx, $gc)
# BenchmarkTools.Trial: 
#   memory estimate:  320 bytes
#   allocs estimate:  4
#   --------------
#   minimum time:     1.329 μs (0.00% GC)
#   median time:      1.434 μs (0.00% GC)
#   mean time:        1.542 μs (8.11% GC)
#   maximum time:     1.256 ms (99.60% GC)
#   --------------
#   samples:          10000
#   evals/sample:     10


