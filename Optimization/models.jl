


function ode02!(du,u,p,t)
    (S, R) = u
    (rₛ, rᵣ, K) = p
    N = S+R
    growthₛ = rₛ*S*(1-N/K) - α*u*S - dₛ*u*S
    growthᵣ = rᵣ*R*(1-N/K) + α*u*S - dᵣ*u*R
    @inbounds begin
        du[1] = growthₛ
        du[2] = growthᵣ
    end
    nothing
end;

function ode03!(du,u,p,t)
    (S, R) = u
    (rₛ, rᵣ, K, ϵ, γ) = p
    N = S+R
    growthₛ = rₛ*S*(1-N/K) - (ϵ+α*u)*S - dₛ*u*S + γ*R
    growthᵣ = rᵣ*R*(1-N/K) + (ϵ+α*u)*S - dᵣ*u*R - γ*R
    @inbounds begin
        du[1] = growthₛ
        du[2] = growthᵣ
    end
    nothing
end;
