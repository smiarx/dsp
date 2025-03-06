OnePoleLP : Filter {

	*ar { arg in = 0.0, freq = 800, mul = 1.0, add = 0.0;
		^this.multiNew('audio', in, freq).madd(mul, add)
	}
}
OnePoleHP : Filter {

	*ar { arg in = 0.0, freq = 800, mul = 1.0, add = 0.0;
		^this.multiNew('audio', in, freq).madd(mul, add)
	}
}
OnePoleAP : Filter {

	*ar { arg in = 0.0, freq = 800, mul = 1.0, add = 0.0;
		^this.multiNew('audio', in, freq).madd(mul, add)
	}
}

SVFLP : Filter {

	*ar { arg in = 0.0, freq = 800, res=0.70710678, mul = 1.0, add = 0.0;
		^this.multiNew('audio', in, freq, res).madd(mul, add)
	}
}

SVFHP : Filter {

	*ar { arg in = 0.0, freq = 800, res=0.70710678, mul = 1.0, add = 0.0;
		^this.multiNew('audio', in, freq, res).madd(mul, add)
	}
}

SVFAP : Filter {

	*ar { arg in = 0.0, freq = 800, res=0.70710678, mul = 1.0, add = 0.0;
		^this.multiNew('audio', in, freq, res).madd(mul, add)
	}
}

SVFBP : Filter {

	*ar { arg in = 0.0, freq = 800, res=0.70710678, mul = 1.0, add = 0.0;
		^this.multiNew('audio', in, freq, res).madd(mul, add)
	}
}

SVFNotch : Filter {

	*ar { arg in = 0.0, freq = 800, res=0.70710678, mul = 1.0, add = 0.0;
		^this.multiNew('audio', in, freq, res).madd(mul, add)
	}
}

LadderLP : Filter {

	*ar { arg in = 0.0, freq = 800, res=0, mul = 1.0, add = 0.0;
		^this.multiNew('audio', in, freq, res).madd(mul, add)
	}
}

LadderHP : Filter {

	*ar { arg in = 0.0, freq = 800, res=0, mul = 1.0, add = 0.0;
		^this.multiNew('audio', in, freq, res).madd(mul, add)
	}
}

LadderAP : Filter {

	*ar { arg in = 0.0, freq = 800, res=0, mul = 1.0, add = 0.0;
		^this.multiNew('audio', in, freq, res).madd(mul, add)
	}
}
