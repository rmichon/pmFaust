declare name "KarplusStrong";
declare description "Simple call of the Karplus-Strong model for the Faust physical modeling library";
declare license "MIT";
declare copyright "(c)Romain Michon, CCRMA (Stanford Univeristy), GRAME";

import("stdfaust.lib");
import("../pm.lib");

process = ks_ui_MIDI <: _,_;
