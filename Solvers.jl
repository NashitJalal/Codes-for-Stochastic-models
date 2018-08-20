M=readdlm("m.txt")
s=size(M)[1]
det(M)              ## Size of M matrix
###### Solving the ODE using RK4 method ####
@time begin
    function fokker(du,u,p,t)
        for i=1:s
            sum=0
            for j=1:s
                sum+=M[i,j]*u[j]
            end
            du[i]=sum
        end
    end
    using DifferentialEquations
    u0 =zeros(s)
    u0[1]=1.0
    tspan = (0.0001,0.01)
    prob = ODEProblem(fokker,u0,tspan)
    sol = solve(prob,RK4())
    a=sol(0.0002)
       end
