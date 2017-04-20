import("stdfaust.lib");
import("../pm.lib");

/*
fluteEmbouchure(pressure) = (_ <: _,_),_,_ : _,*(0.5)+(pressure-_*0.5 : fluteJetTable),_;

fluteHead = lTermination(*(absorption),basicBlock)
with{
  absorption = 0.95; // TODO: may be should be controllable
};

fluteFoot = rTermination(basicBlock,*(absorption) : dispersion)
with{
  dispersion = si.smooth(0.7);
  absorption = 0.95;
};

// TODO: problem with tuning!!!
basicFluteModel(tubeLength,mouthPosition,pressure) = endChain(fluteChain) : fi.dcblocker
with{
  maxTubeLength = 3; // meters
  tubeTuning = 0.27;
  tLength = tubeLength+tubeTuning;
  embouchurePos = 0.27 + (mouthPosition-0.5)*0.4;
  tted = tLength*embouchurePos;
  eted = tLength*(1-embouchurePos);
  fluteChain = chain(fluteHead : openTube(maxTubeLength,tted) : fluteEmbouchure(pressure) : openTube(maxTubeLength,eted) : fluteFoot : out);
};

basicFluteModel_ui(pressure) =
basicFluteModel(tubeLength,mouthPosition,pressure)*outGain
with{
  tubeLength = hslider("v:fluteModel/[0]tubeLength",0.8,0.01,3,0.01) : si.smoo;
  mouthPosition = hslider("v:fluteModel/[1]mouthPosition",0.5,0,1,0.01) : si.smoo;
  outGain = hslider("v:fluteModel/[2]outGain",0.5,0,1,0.01);
};

basicFlute_ui = hgroup("flute",blower_ui : basicFluteModel_ui);

basicFlute_ui_MIDI = basicFluteModel(tubeLength,mouthPosition,blow)*outGain
with{
  f = hslider("v:fluteMidi/v:[0]midi/[0]freq",440,50,1000,0.01);
	bend = hslider("v:fluteMidi/v:[0]midi/[1]bend[hide:1][midi:pitchwheel]",1,0,10,0.01) : si.polySmooth(gate,0.999,1);
	gain = hslider("v:fluteMidi/v:[0]midi/[2]gain",0.9,0,1,0.01);
  s = hslider("v:fluteMidi/v:[0]midi/[3]sustain[hide:1][midi:ctrl 64]",0,0,1,1);
	t = button("v:fluteMidi/v:[0]midi/[4]gate");
  mouthPosition = hslider("v:fluteMidi/v:[1]model/[0]mouthPosition",0.5,0,1,0.01) : si.smoo;
  envAttack = hslider("v:fluteMidi/v:[2]other/[0]envAttack",10,0,30,0.01)*0.001;
  vibratoFreq = hslider("v:fluteMidi/v:[2]other/[1]vibratoFreq",5,1,10,0.01);
	vibratoGain = hslider("v:fluteMidi/v:[2]other/[2]vibratoGain",0.5,0,1,0.01)*0.04;
  outGain = hslider("v:fluteMidi/v:[2]other/[3]outGain",0.5,0,1,0.01);

  gate = t+s : min(1);
	freq = f*bend;
	envelope = gate*gain : si.smooth(ba.tau2pole(envAttack));

	tubeLength = freq : f2l;
	pressure = envelope;
	blow = blower(pressure,0.05,2000,vibratoFreq,vibratoGain);
};

process = basicFlute_ui_MIDI <: _,_;
*/

//////////////////////////////////
// BRASS
//////////////////////////////////

/*
brassLipsMouthPieceTable(tubeLength,lipsTension) = *(0.03) : lipFilter <: * : clipping
with{
  clipping = min(1) : max(-1);
  freq = (tubeLength : l2f)*pow(4,(2*lipsTension)-1);
  filterR = 0.994009;
  a1 = -2*filterR*cos(ma.PI*2*freq/ma.SR);
  lipFilter = fi.tf2(1,0,0,a1,filterR);
};

brassLipsMouthPiece(tubeLength,lipsTension,pressure) = lTermination(mpInteraction,basicBlock)
with{
  reflexion = *(0.85);
  p = pressure*0.3;
  mpInteraction = reflexion <:
  (p-_ : brassLipsMouthPieceTable(tubeLength,lipsTension) <: *(p),1-_),_ : _,* : + : fi.dcblocker;
};

basicBrassModel(tubeLength,lipsTension,mute,pressure) = endChain(brassChain)
with{
  maxTubeLength = 3; // meters
  lengthTuning = 0; // TODO
  tunedLength = tubeLength + lengthTuning;
  brassChain = chain(brassLipsMouthPiece(tunedLength,lipsTension,pressure) : openTube(maxTubeLength,tunedLength) : wBell(mute) : out);
};

brassModel_ui = basicBrassModel(tubeLength,lipsTension,mute,pressure)
with{
  pressure = hslider("v:brassModel/[0]pressure",0,0,1,0.01) : si.smoo;
  tubeLength = hslider("v:brassModel/[1]tubeLength",0.5,0.01,2.5,0.01) : si.smoo;
  lipsTension = hslider("v:brassModel/[2]lipsTension",0.5,0,1,0.01) : si.smoo;
  mute = hslider("v:brassModel/[3]mute",0.5,0,1,0.01) : si.smoo;
};

brassModel_ui_MIDI = basicBrassModel(tubeLength,lipsTension,mute,pressure)*outGain
with{
  f = hslider("v:brassMidi/v:[0]midi/[0]freq",440,50,1000,0.01);
	bend = hslider("v:brassMidi/v:[0]midi/[1]bend[hide:1][midi:pitchwheel]",1,0,10,0.01) : si.polySmooth(gate,0.999,1);
	gain = hslider("v:brassMidi/v:[0]midi/[2]gain",0.5,0,1,0.01);
  s = hslider("v:brassMidi/v:[0]midi/[3]sustain[hide:1][midi:ctrl 64]",0,0,1,1);
	t = button("v:brassMidi/v:[0]midi/[4]gate");
  lTension = hslider("v:brassMidi/v:[1]model/[0]lipsTension",0.5,0,1,0.01) : si.smoo;
  mute = hslider("v:brassMidi/v:[1]model/[1]mute",0.5,0,1,0.01) : si.smoo;
  envAttack = hslider("v:brassMidi/v:[2]other/[0]envAttack",10,0,30,0.01)*0.001;
  vibratoFreq = hslider("v:brassMidi/v:[2]other/[1]vibratoFreq",5,1,10,0.01);
  vibratoGain = hslider("v:brassMidi/v:[2]other/[2]vibratoGain",0.5,0,1,0.01)*0.04;
  outGain = hslider("v:brassMidi/v:[2]other/[3]outGain",0.5,0,1,0.01);

  gate = t+s : min(1);
	vibrato = 1+os.osc(vibratoFreq)*vibratoGain*envelope;
  lipsTension = lTension*vibrato;
  freq = f*bend;
	envelope = gate*gain : si.smooth(ba.tau2pole(envAttack));

	tubeLength = freq : f2l;
	pressure = envelope*vibrato;
};

process = brassModel_ui_MIDI <: _,_;
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
