import("stdfaust.lib");

modeFilter(f,t60) = fi.tf2(b0,b1,b2,a1,a2)
with{
    b0 = 1;
	b1 = 0;
	b2 = -1;
	w = 2*ma.PI*f/ma.SR;
	r = pow(0.001,1/float(t60*ma.SR));
	a1 = -2*r*cos(w);
	a2 = r^2;
};

mode(f,t60,gain) = modeFilter(f,t60)*gain;

freq(0) = 100;
freq(1) = 200;
freq(2) = 300;
freq(3) = 400;

t60(0) = 0.9;
t60(1) = 0.8;
t60(2) = 0.6;
t60(3) = 0.5;

gain(0) = 0.9;
gain(1) = 0.9;
gain(2) = 0.5;
gain(3) = 0.6;

nModes = 4;

model = _ <: par(i,nModes,mode(freq(i),t60(i),gain(i))) :> _;

process = model;