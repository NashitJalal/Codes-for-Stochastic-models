# Potential File for sawtooth potential
# Linear potential with asymmetry,peak,lowest free energy and state as input
include("Scriptreader.jl")
importall Parameters
function potentialsaw(a::Float64,peak::Float64,low::Float64,state::Int)  # Arguments as Asymmetry,Peak,Low and Chemical State
node=zeros(Parameters.n)

for i=1:Parameters.n
    node[i]=(i-1)*Parameters.grid
end
    BE=zeros(Parameters.n)
    for  i=1:Parameters.n
        d=(i-1)*Parameters.grid
        if d<=a*Parameters.period
            BE[i]=peak*Parameters.kT+(((low*Parameters.kT-peak*Parameters.kT)*(node[i]))/(a*Parameters.period))          # Before the asymmetry point
        else

            BE[i]=low*Parameters.kT+(((peak*Parameters.kT-low*Parameters.kT)*(node[i]-a*Parameters.period))/(Parameters.period-a*Parameters.period)) # After Asymmetry point
        end
    end
    cd(Parameters.pot)
        writedlm("Potential$state.txt",BE)    # Writing the potenial files
        cd(Parameters.h)    # Changing Directory to home directory
    end
