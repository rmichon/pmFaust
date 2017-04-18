import("stdfaust.lib");
import("../pm.lib");

fluteJetTable = _ <: *(* : -(1)) : clipping
with{
  clipping = min(1) : max(-1);
};

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

basicFluteModel(tubeLength,mouthPosition,pressure) = endChain(fluteChain)
with{
  maxTubeLength = 2; // meters
  tubeTuning = 0.27;
  tLength = tubeLength+tubeTuning;
  embouchurePos = 0.27 + (mouthPosition-0.5)*0.2;
  tted = tLength*embouchurePos;
  eted = tLength*(1-embouchurePos);
  fluteChain = chain(fluteHead : openTube(maxTubeLength,tted) : fluteEmbouchure(pressure) : openTube(maxTubeLength,eted) : fluteFoot : out);
};

pressure = hslider("pressure",0,0,1,0.01) : si.smoo;
length = hslider("freq",440,50,2000,0.01) : l2f : si.smoo;
pos = hslider("pos",0.3,0,1,0.01) : si.smoo;
outGain = hslider("outGain",0.5,0,1,0.01);

process = blower_ui : basicFluteModel(length,pos) : fi.dcblocker*outGain <: _,_;


//////////////////////////////////
// BRASS
//////////////////////////////////

/*
brassLipsMouthPiece(freq,lipTension) = *(0.03) : fi.fi.resonbp(lipFilterFrequency,2,1) <: * : clipping
with{
  clipping = min(1);
  lipFilterFrequency = freq*pow(4,(2*lipTension)-1);
};

basicBrassModel(tubeLength) = endChain(brassChain)
with{
  maxTubeLength = 1; // meters
  brassChain = chain(openTube(maxTubeLength,tubeLength) : out);
};

process = basicBrassModel(0.5);
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
