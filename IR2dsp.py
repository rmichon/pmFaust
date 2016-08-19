from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from sys import argv
import subprocess
from scipy.io.wavfile import read
import peakutils

script, sound, model, Nwindow, t = argv

N = int(Nwindow)
(fs, x2) = read(sound)
t60 = t
modelName = model
x = x2[0:N] #windowing brutal
x = x/max(x)

#FFT
fft_size = N
X = np.abs(np.fft.fft(x, fft_size))

#computing corresponding frequencies
time_step = 1 / 44100
freqs = np.fft.fftfreq(x.size, time_step) 
idx = np.argsort(freqs) 

#detecting peaks
indexes = peakutils.indexes(X[idx], thres=0.01, min_dist=fft_size/100)
    
#plt.plot(freqs[idx], X[idx])
#plt.show()

#Storing frequencies and modes for each bp filters
peaksFreq = []
peaksMode = []
for i in indexes:
    if freqs[idx][i] > 0:
        peaksFreq.append(freqs[idx][i])
        peaksMode.append(X[idx][i])
peaksMode = peaksMode/(max(peaksMode))

print "peaks frequencies :"
print peaksFreq
print "corresponding gains :"
print peaksMode


# Writing the dsp file #
###########################
file = open(modelName + ".dsp", "w")

file.write("import(\"architecture/pm.lib\");\n")
file.write("import(\"music.lib\");\n\n")
file.write("pi = 4*atan(1.0);\n")
file.write("nModes = ")
file.write(str(len(peaksMode)))
file.write(";\n")
file.write("eigenValues = ("); #writing the frequencies list
k = 0
for i in peaksFreq :
    file.write(str(i))
    if(k+1 < len(peaksMode)):
        file.write(", ")
    k += 1
file.write(");\n");
file.write("massEigenValues = ("); #writing the masses list
k = 0
for i in peaksMode :
    file.write(str(i))
    if(k+1 < len(peaksMode)):
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
file.write('\n gate = button("gate");')
file.write('\n process = impulseExcitation(gate) : ' + modelName + ' <: _,_;')

file.close();



