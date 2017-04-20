declare name "Violin";
declare description "Simple violin physical model controlled with continuous parameters.";
declare license "MIT";
declare copyright "(c)Romain Michon, CCRMA (Stanford University), GRAME";

import("stdfaust.lib");
import("../pm.lib");

process = violin_ui <: _,_;
