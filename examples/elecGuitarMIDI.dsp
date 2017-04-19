declare name "ElecGuitarMidi";
declare description "Simple electric guitar model without effect chain.";
declare license "MIT";
declare copyright "(c)Romain Michon, CCRMA (Stanford University), GRAME";

import("stdfaust.lib");
import("../pm.lib");

// TODO: We could potentially add an audio effect chain here

process = elecGuitar_ui_MIDI <: _,_;
