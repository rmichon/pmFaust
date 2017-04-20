declare name "NylonGuitarMidi";
declare description "Simple acoustic guitar model with nylon strings.";
declare license "MIT";
declare copyright "(c)Romain Michon, CCRMA (Stanford University), GRAME";

import("stdfaust.lib");
import("../pm.lib");

process = nylonGuitar_ui_MIDI <: _,_;
