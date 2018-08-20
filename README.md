
# Title : Generalized code for Stochastic Models ~ {KAUSTUBH RANE RESEARCH GROUP}
  
  About Kinesin Motors........
Kinesins are biological protein motors in animal cells.They are associated with the cargo movements such as vesicles, organelles and chromosomes along a track called Microtubule.
Broadly the Kinesin can be classified as Single headed and Double headed Kinesin.These molecular motors gets energy by hydrolysis of ATP. Physically the phenomenon can be studied as a Brownian Motion while mathematically we can model it through Fokker Planck equation as the system follows Langevin Dynamics.

  About the Package....
 The package provides a generalized code for studying the biomolecular transport phenomenon of both the monomeric as well as Dimeric Kinesin with any number of chemical states. It include potential generators as well as different mathematical methods to solve the Stochastic Partial Differnetial Equation. 
  It basically uses the WPE method to generate the Transition Matrix which are then solved by various numerical techniques. The entire code assumes a 1D motion with spring interaction between heads in case of dimeric Kinesin.
  Coding has been done in such a way that modification of the package is easy for future work. The codes are written in Julia Language.
  <https://julialang.org/>

Instructions...
 1. "kinscript.txt" is the input data file which can be edited by the user. Place it in the same directory with other .jl files. Place potentials with name "Potential$State.txt" in "Potential Folder".
 2. You can generate different potentials using generator functions ("Sawtoothpotential.jl",       ). The files will automatically be stored in the "Potential" Folder
 3. Enter the kinetic rate constants in the "Kinetic.txt" file in "Kineticrates" folder . It contains 2 column of Forward and Backward rates respectively.
 4. Execute the "Scriptreader.jl "file followed by the "Transition.jl"
 5. Execute any of the Solver Files (                                             )
 6. The Results are stored in the "Output" folder.
 
 
 "Sawtoothpotential.jl" --- call potentialsaw(Asymmetry, Peak value of free energy,Lowest free energy,State number as Input)
 Note :  All entry in decimal except the State Number.......  

#### Note : Currently the code can only be used for finding stationary distributions.Rest of the codes will be updated soon

Contributors: Nashit Jalal , Atul Sharma , Kaustubh Rane.

