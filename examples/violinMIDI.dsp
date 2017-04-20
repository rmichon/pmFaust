declare name "ViolinMidi";
declare description "Simple MIDI-controllable violin physical model.";
declare license "MIT";
declare copyright "(c)Romain Michon, CCRMA (Stanford University), GRAME";

import("stdfaust.lib");
import("../pm.lib");

process = violin_ui_MIDI <: _,_;
