import("stdfaust.lib");
import("../pm.lib"); // TODO remove

// TODO: ready to go in demo.lib
ksString_demo = gate : impulseExcitation : ksString( (freq:f2l), damping )
with{
	f = hslider("v:karplus/[0]freq[0]",440,50,1000,0.01);
	bend = hslider("v:karplus/[1]bend[1][hide:1][midi:pitchwheel]",1,0,10,0.01) : si.polySmooth(gate,0.999,1);
	gain = hslider("v:karplus/[2]gain",0.8,0,1,0.01);
	damping = hslider("v:karplus/[3]damping[midi:ctrl 1]",0.5,0,1,0.01);
	t = button("v:karplus/[4]gate");
	s = hslider("sustain[hide:1][midi:ctrl 64]",0,0,1,1);

	gate = t+s : min(1);
	freq = f*bend;
};

process = ksString_demo;