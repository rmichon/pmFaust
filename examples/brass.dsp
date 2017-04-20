declare name "Brass";
declare description "Simple brass instrument physical model with physical parameters.";
declare license "MIT";
declare copyright "(c)Romain Michon, CCRMA (Stanford University), GRAME";

import("stdfaust.lib");
import("../pm.lib");

process = brass_ui <: _,_;
