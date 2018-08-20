module Parameters
f=open("kinscript.txt")   # Script file
entry=String[]
entry=readlines(f)
for i=1:11
    x=findfirst(isequal('%'),entry[i])
    entry[i]=entry[i][1:x-1]
end
choice=entry[1]          # Monomeric or Dimeric (M/D)
chemical=entry[2]
chemical=parse(Int,chemical)   # No of chemical states
period=entry[3]
period=parse(Float64,period)  # Period
n=entry[4]
n=parse(Int,n)                # Nodes
D=entry[5]
D=parse(Float64,D)            # Diffusion Constant
T=entry[6]
T=parse(Float64,T)            # Temperature in Kelvin
pot=entry[7]                  # Path for potential file
k=entry[8]
k=parse(Float64,k)            # Spring Constant
I=entry[9]                    # Interaction between heads
kineticfile=entry[10]         # Path for kinetic rate Constants
grid=period/(n-1)
kb=1.38064852e-23
kT=kb*T
h=entry[11]                   #Root Directory
#export all
end
