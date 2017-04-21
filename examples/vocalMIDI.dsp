declare name "VocalMIDI";
declare description "Simple MIDI-controllable source-filter vocal synthesizer.";
declare license "MIT";
declare copyright "(c)Romain Michon, CCRMA (Stanford University), GRAME";

import("stdfaust.lib");
import("../pm.lib");

process = SFFormantModel_ui_MIDI <: _,_;
