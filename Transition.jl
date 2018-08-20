include("Scriptreader.jl")
x=readdir(Parameters.pot)
Bienergy=zeros(Parameters.n)
cd(Parameters.pot)
l=length(x)
for i=1:l
    if i==1
    Bienergy=readdlm("Potential1.txt")
else
    Bienergy=hcat(Bienergy,readdlm("Potential$i.txt"))   # Making Matrix containing the potentials at different states
  end
end
cd(Parameters.kineticfile)
x=readdir(Parameters.kineticfile)
mat=Array{Float64,2}                            # Matrix for Storing Kinetic rate constants
mat=readdlm("Kinetic.txt")
cd(Parameters.h)

######### Forward and Backward Jump Rates############
#if Parameters.choice="M"

    L=zeros(Parameters.n,Parameters.n,Parameters.chemical)
    for p=1:Parameters.chemical  # Loop for M matrix for each chemical state
        du=zeros(Parameters.n-1); dp=zeros(Parameters.n); dl=zeros(Parameters.n-1)
        Bi=Bienergy[:,p]
        for i=1:Parameters.n-1  # For Lower Diagonal
            alpha=-(Bi[i+1]-Bi[i])/Parameters.kT
            if alpha==0
                F=Parameters.D/(Parameters.grid)^2
            else
                F= Parameters.D/(Parameters.grid)^2 * -alpha/(exp(-alpha)-1)
            end
            dl[i]=F
        end    # .. For Lower Diagonal


        for i=1:Parameters.n    # For principal diagonal
            if i==Parameters.n
                alpha= -(Bi[1]-Bi[Parameters.n])/Parameters.kT
                if alpha==0
                    F=Parameters.D/(Parameters.grid)^2
                else
                    F= Parameters.D/(Parameters.grid)^2 * -alpha/(exp(-alpha)-1)
                end

            else
                alpha= -(Bi[i+1]-Bi[i])/Parameters.kT
                if alpha==0
                    F=Parameters.D/(Parameters.grid)^2
                else
                    F=Parameters.D/(Parameters.grid)^2 * -alpha/(exp(-alpha)-1)
                end
            end

            if i==1
                alpha= -(Bi[1]-Bi[Parameters.n])/Parameters.kT
                if alpha==0
                    B=Parameters.D/(Parameters.grid)^2
                else
                    B= Parameters.D/(Parameters.grid)^2 * alpha/(exp(alpha)-1)
                end
            else
                alpha= -(Bi[i]-Bi[i-1])/Parameters.kT;
                if alpha==0
                   B=Parameters.D/(Parameters.grid)^2
                 else
                   B= Parameters.D/(Parameters.grid)^2 * alpha/(exp(alpha)-1)
                 end
            end
         dp[i]=-(F+B)
end # .. For Principal Diagonal

for i=1:Parameters.n-1    # For Upper Diagonal
    alpha= -(Bi[i+1]-Bi[i])/Parameters.kT
    if alpha==0
        B=Parameters.D/(Parameters.grid)^2
    else
       B=Parameters.D/(Parameters.grid)^2 * alpha/(exp(alpha)-1)
    end
   du[i]=B
end    #.. For Upper Diagonal
L[:,:,p]=Tridiagonal(dl,dp,du)
                         # Tridiagonal Matrix

###########Applying Periodic Boundary Conditions###############
alpha= -(Bi[1]-Bi[Parameters.n])/Parameters.kT
if alpha==0
     L[Parameters.n,1,p]= Parameters.D/(Parameters.grid)^2
     L[1,Parameters.n,p]= Parameters.D/(Parameters.grid)^2
 else
      L[Parameters.n,1,p]= Parameters.D/(Parameters.grid)^2 * alpha/(exp(alpha)-1)  # B_1\2
      L[1,Parameters.n,p]= Parameters.D/(Parameters.grid)^2* -alpha/(exp(-alpha)-1) # F_N+1/2
end
end

##### K (Matrix of Kinetic Rate Constant)  #######
cd(Parameters.kineticfile)
kin =readdlm("Kinetic.txt")
K=zeros(Parameters.chemical,Parameters.chemical)   # Matrix for Storing Rate Constants
b=0
for i=1:Parameters.chemical         # Matrix stored according to theoritical indexing
    if Parameters.chemical ==2
        K[1,2]=kin[1,1]
        K[2,1]=kin[1,2]
    else
    b=i+1
    if b>Parameters.chemical
        b=1
    end
    K[i,b]=kin[i,1]
    K[b,i]=kin[i,2]
  end
end
  cd(Parameters.h)   # Root directory
##### M (Transition matrix with Kinetic Rate Constant) Matrix #######
M=zeros(1,Parameters.n*Parameters.chemical)  # Storing the M Matrix
zer=zeros(Parameters.n,Parameters.n)
one=eye(Parameters.n)
for i=1:Parameters.chemical
    a=0
    temp=zeros(4,4)      # Temporary file for storing the rows
    for j=1:Parameters.chemical
        a=a+1
        if j==i
            if j==1
                temp=L[:,:,1]-(K[1,2]*one)-(K[1,Parameters.chemical]*one)
            elseif j==Parameters.chemical
                temp=hcat(temp,L[:,:,Parameters.chemical]-(K[Parameters.chemical,1]*one)-(K[Parameters.chemical,Parameters.chemical-1]*one))
            else
                temp=hcat(temp,L[:,:,i]-(K[i,i+1]*one)-(K[i,i-1]*one))

            end

    else
        if j==1
            if K[a,i]==0
                temp=zer
            else
                temp=K[a,i]*one
            end
        else
            if K[a,i]==0
                temp=hcat(temp,zer)
            else
                temp=hcat(temp,K[a,i]*one)
            end
        end
    end
end
M=vcat(M,temp)                      # Concatinating to build the M matrix
end
row=1
M=M[setdiff(1:end,row),:]                # Removing the extra Ist row in M matrix
writedlm("m.txt",M)

########## Finding the Stationary Solution#############
J=zeros(Parameters.n,1)
M[1,:]=1
 B=zeros(Parameters.n*Parameters.chemical)
  B[1]=1
 Z=M\B;
 for i=1:Parameters.n
     J[i]=(Z[i]+Z[i+Parameters.n])/(Parameters.grid*10^9) ##### Probability Distribution
 end
 writedlm("StationarySolution.txt",J)
 ###########################
