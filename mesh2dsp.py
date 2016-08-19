from sys import argv
import os
import subprocess
import string
import re
import math

script, stlFile, modelName, young, poisson, density, t = argv

fileName = stlFile[0:len(stlFile)-4]
print fileName

#subprocess.call(["gmsh", "-2", "-o", stlFile])
#mshFile = fileName + ".msh"
#subprocess.call(["ElmerGrid", "14", "2", mshFile)

# Writing of the sif file#
##########################

sifFile = open("case.sif", "w")
sifFile.write("""Header
Mesh DB "." " """,fileName,""" "
Include Path ""
Results Directory ""
End

Simulation
Coordinate System = "Cartesian 3D"
Coordinate Mapping(3) = 1 2 3
Simulation Type = "Steady State"
Steady State Max Iterations = 1
Solver Input File = "eigen_values.sif"
Output File = "eigen_values.dat"
Post File = "eigen_values.ep"
End

Body 1
Equation = 1
Material = 1
End
Material 1
""")
sifFile.write("Youngs Modulus = ")
sifFile.write(young)
sifFile.write("\nPoisson Ratio = ")
sifFile.write(poisson)
sifFile.write("\nDensity = ")
sifFile.write(density)
sifFile.write("""
End

Equation 1
Stress Analysis = True
End

Solver 1
Equation = "Stress Analysis"
Eigen Analysis = Logical True
Eigen System Values = Integer 20
Linear System Solver = "direct"
Variable = "Displacement"
Variable Dofs = 3
Linear System Iterative Method = "BiCGStab"
Linear System Max Iterations = 1000
Linear System Convergence Tolerance = 1.0e-08
Linear System Abort Not Converged = True
Linear System Preconditioning = "ILU0"
Linear System Residual Output = 1
Steady State Convergence Tolerance = 1.0e-05
Nonlinear System Convergence Tolerance = 1.0e-05
Nonlinear System Max Iterations = 1
Nonlinear System Newton After Iterations = 3
Nonlinear System Newton After Tolerance = 1.0e-02
Nonlinear System Relaxation Factor = 1
Linear System Precondition Recompute = 1
End

Boundary Condition 1
Target Boundaries(1) = 1
Displacement 1 = 0
Displacement 2 = 0
Displacement 3 = 0
End

Boundary Condition 2
Target Boundaries(1) = 2
Displacement 1 = 0
Displacement 2 = 0
Displacement 3 = 0
End

""")
sifFile.close()
           
output = subprocess.check_output(["ElmerSolver case.sif"], shell=True)
index1 = string.rfind(output, "Computed Eigen Values")
index2 = string.rfind(output, "EigenSolve")
output = output[index1:index2]

#Extracting the list of eigen values
eigenVS = re.findall("\d+\.\d+", output)
eigenV = [float('%.3f'%float((x))) for x in eigenVS]
eigenV = filter(lambda s: s != 0.0, eigenV)

frequencies = [math.sqrt(x)/(2*math.pi) for x in eigenV]
print frequencies
nModes = len(frequencies)
masses = [1] * nModes
print masses
t60 = t;



# Writing of the dsp file #
###########################
file = open("modalModel.dsp", "w")

file.write("import(\"architecture/pm.lib\");\n")
file.write("import(\"music.lib\");\n\n")
file.write("pi = 4*atan(1.0);\n")
file.write("nModes = ")
file.write(str(nModes))
file.write(";\n")
file.write("eigenValues = ("); #writing the frequencies list
k = 0
for i in frequencies :
    file.write(str(i))
    if(k+1 < nModes):
        file.write(", ")
    k += 1
file.write(");\n");
file.write("massEigenValues = ("); #writing the masses list
k = 0
for i in masses :
    file.write(str(i))
    if(k+1 < nModes):
        file.write(", ")
    k += 1
file.write(");\n");
file.write("modeFreqs = par(i,nModes,sqrt(take(i+1,eigenValues))*2*pi);\n")
file.write("modeGains = par(i,nModes,take(i+1,massEigenValues));\n")
file.write("t60 = ")
file.write(str(t60))
file.write(";\n\n")
file.write(modelName)
file.write(" = modalModel(nModes,modeFreqs,modeGains,t60);");

file.close();
