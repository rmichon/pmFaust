import("stdfaust.lib");
import("../pm.lib");




////////////////////////////////////////////////////////
// CLARINET MODEL STUFF -> SHOULD MOVE THE PM.LIB
////////////////////////////////////////////////////////

reed(offset,slope) = reedTable : min(1) : max(-1)
with {
	reedTable = *(slope) + offset;
};

clarinetReed(stiffness) = *(-1) <: *(reed(0.7,tableSlope))
with{
	tableSlope = -0.44 + 0.26*stiffness;
};

// wood brass bell
wbBell(opening) = si.smooth(opening);

// TODO may be add lowpass on noise
blower(pressure,breathGain,breathCutoff) = pressure + breathNoise
with{
	breathNoise = no.noise : fi.lowpass(2,breathCutoff) : *(pressure*breathGain);
};

blower_ui = blower(pressure,breathGain,breathCutoff)
with{
	pressure = hslider("v:blower/[0]pressure",0,0,1,0.01) : si.smoo;
	breathGain = hslider("v:blower/[1]breathGain",0.1,0,1,0.01)*0.05;
	breathCutoff = hslider("v:blower/[2]breathCutoff",2000,20,20000,0.1);
};

singleReedMod(length,pressure,reedStiffness,bellOpening) = endChain(modelChain)
with{
	delTuning = 3;
	delLength = length*ma.SR/speedOfSound/2 - delTuning;
	bore = waveguide(maxDel,delLength);
	modelChain = terminations( clarinetReed(reedStiffness), chain( in(pressure) : bore : out ), wbBell(bellOpening) );
};

singleReedMod_ui(pressure) =  singleReedMod(tubeLength,pressure,reedStiffness,bellOpening)*0.9
with{
	tubeLength = hslider("v:singleReed/[0]tubeLength",0.8,0.01,3,0.01) : si.smoo;
	reedStiffness = hslider("v:singleReed/[1]reedStiffness",0.5,0,1,0.01);
	bellOpening = hslider("v:singleReed/[2]bellOpening",0.5,0,1,0.01);
};

singleReedMIDI =  singleReedMod(tubeLength,blow,reedStiffness,bellOpening)*0.9
with{
	f = hslider("v:singleReedMidi/[0]freq",440,50,1000,0.01);
	bend = hslider("v:singleReedMidi/[1]bend[hide:1][midi:pitchwheel]",1,0,10,0.01) : si.polySmooth(gate,0.999,1);
	gain = hslider("v:singleReedMidi/[2]gain",0.8,0,1,0.01);
	reedStiffness = hslider("v:singleReedMidi/[3]reedStiffness[midi:ctrl 1]",0.5,0,1,0.01);
	bellOpening = hslider("v:singleReedMidi/[4]bellOpening[midi:ctrl 3]",0.5,0,1,0.01);
	vibratoFreq = hslider("v:singleReedMidi/[5]vibratoFreq",6,1,10,0.01);
	vibratoGain = hslider("v:singleReedMidi/[6]vibratoGain",0.25,0,1,0.01)*0.01;
	envAttack = hslider("v:singleReedMidi/[7]envAttack",1,0,30,0.01)*0.001;
	t = button("v:singleReedMidi/[8]gate");
	s = hslider("v:singleReedMidi/sustain[hide:1][midi:ctrl 64]",0,0,1,1);

	gate = t+s : min(1);
	vibrato = 1+os.osc(vibratoFreq)*vibratoGain*envelope;
	freq = f*bend*vibrato;
	envelope = gate*0.6 : si.smooth(ba.tau2pole(envAttack));

	tubeLength = freq : f2l;
	pressure = envelope*vibrato;
	blow = blower(pressure,0.05,2000);
};

clarinetInstr_demo = hgroup("clarinet",blower_ui : singleReedMod_ui);

// clarinet MIDI instr goes here

process = singleReedMIDI <: _,_;

//////////////////////////////////////////////////////////////
// END OF CLARINET MODEL STUFF -> SHOULD MOVE THE PM.LIB
//////////////////////////////////////////////////////////////


/*
ksString(length,damping,excitation) = endChain(ksChain)
with{
	delTuning = 6;
	delLength = length*ma.SR/speedOfSound/2 - delTuning;
	refCoef = (1-damping)*0.2+0.8;
	refFilter = _ <: (_+_')/2*refCoef;
	ksChain = terminations(_,chain(in(excitation) : waveguide(maxDel,delLength) : out),refFilter);
};

f = hslider("freq",400,50,1000,0.01);
g = button("gate");
damping = hslider("damping",0.5,0,1,0.01);

process = g : impulseExcitation : ksString( (f:f2l), damping );
*/











//up(upProc) = _,(_ <: upProc,_),_ : +,_,_;
//down(downProc) = (_ <: _,downProc),_,_ : _,+,_;
//foo = _,_,_;

//process = chain( down(*(-0.99)) : up(*(-0.99)) );
//process = chain( down(*(-0.99)) : foo );

/*
f = hslider("freq",400,50,1000,0.01);
g = button("gate");

d = ma.SR/f/2;

myStr = 0,_,0 : chain( down(*(-0.99)) : waveguide(512,d) : up(*(-0.99)) ) : !,_,!;

process = g : ba.impulsify : myStr <: _,_;
*/

//process = chain( waveguide(512,d) : up(*(0.9)) );

//process = chain( waveguide(512,3) : down(_) : up(_) : waveguide(512,50));

//process = chain(waveguide(512,d) : waveguide(256,d) : waveguide(254,d));

//process = fullTerminations(*(0.9),chain(waveguide(512,2) : waveguide(512,2)),_);





// EXAMPLE 1

//process = feedback(_,+(1));

// EXAMPLE 2

/*
myString(freq,gate,pos) = acousticExcitation(gate,n) : resonator
with{
	nMax = 512;
	nAdjust = 4; // because 2 waveguides are used
	n = ma.SR/freq/2 - nAdjust;
	nUp = n*(pos);
	nDown = n*(1-pos);
	terminations = _ <: (_+_')/2 : *(-1);
	resonator(x) = fullTerminations(terminations,chain(waveguide(nMax,nUp) : input(x) : waveguide(nMax,nDown) : output),terminations);
};

f = hslider("freq",400,50,1000,0.01);
p = hslider("pos",0.5,0,1,0.01);
g = button("gate");

process = myString(f,g,p);
*/