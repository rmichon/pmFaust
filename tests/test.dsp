import("stdfaust.lib");
import("../pm.lib");

//process = waveguide(512,8);

// TODO: Now "real" plucked string goes here


/*
bowTable(offset,slope) = pow(abs(sample) + 0.75, -4) : min(1)
with{
	sample = +(offset)*slope;
};

bow(bowPressure) = bowTable(0,tableSlope)
with{
	tableSlope = 5 - (4*bowPressure);
};

rStringRigidTermination = rTermination(basicBlock,*(-1));

lStringRigidTermination = lTermination(*(-1),basicBlock);

nut = lStringRigidTermination;

multStringNut(nStrings) = par(i,nStrings,nut);

bridge(reflexion,absorption) = rTermination(basicBlock,reflectance) : _,transmittance,_
with{
	// absorption is typically 0.6
	reflectance = *(-reflexion) : si.smooth(absorption);
	transmittance = _; // TODO: should be the inverse of refelctance
};


// TODO misssing sympathetic resonances
multStringBridge(nStrings,reflexion,absorption) = par(i,nStrings,singleBridge(reflexion,absorption)) :> _,_,_;

violinBody = reflectance,transmittance,_
with{
	transmittance = fi.resonbp(500,2,1);
	reflectance = _;
};

interact(b) = (_,_ <: b,_,_ :> _,_),_;

bowInteraction(bowPressure,bowVelocity) = interact(bowSystem)
with{
	bowSystem = + : bowVelocity-_ <: *(bow(bowPressure)) <: _,_;
};

// TODO: freq for other models for high level
bowedStringModel(bowPressure,bowVelocity,bridgeReflexion,bridgeAbsorption,bowPosition,stringLength) = endChain(modelChain)
with{
	stringTuning = 0.95;
	stringL = stringLength*stringTuning;
	ntbd = stringL*bowPosition;
	btbd = stringL*(1-bowPosition);
	modelChain = chain(
			   nut :
			   openString(ntbd) :
			   bowInteraction(bowPressure,bowVelocity) :
			   openString(btbd) :
			   bridge(bridgeReflexion,bridgeAbsorption) :
			   violinBody :
			   out
	);
};

bowVel = hslider("bowVel",0,0,1,0.01) : si.smoo;
bowPress = hslider("bowPress",0.5,0,1,0.01);
bowPos = hslider("bowPos",0.7,0,1,0.01);
bridgeReflexion = hslider("bridgeReflexion",0.95,0,1,0.01);
bridgeAbsorption = hslider("bridgeAbsorption",0.6,0,1,0.01);
length = hslider("length",0.75,0,1,0.01);

process = bowedStringModel(bowPress,bowVel,bridgeReflexion,bridgeAbsorption,bowPos,length) <: _,_;
*/

////////////////////////////////////////////////////////
// CLARINET MODEL STUFF -> SHOULD MOVE THE PM.LIB
////////////////////////////////////////////////////////


// clarinet MIDI instr goes here

process = basicClarinet_ui_MIDI <: _,_;
//l = hslider("l",0.5,0,1,0.01);
//process = openString(1,0.5);

//process = de.delay((1 : l2s/2 : ma.np2), l);

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
