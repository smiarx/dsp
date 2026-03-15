VoiceBox : UGen {
	*ar { arg in = 0.0, lambda=0.995, warpin=0, warpout=0, ptop=0.99, ptrho=0.993, ptlambda=0.99, shift=1.0, mul=1.0, add=0.0;
		^this.multiNew('audio', in, lambda, warpin,warpout, ptop, ptrho, ptlambda, shift).madd(mul, add)
	}
}

ZC : UGen {
	*ar { arg in = 0.0, dur=0.05, mul=1.0, add=0.0;
		^this.multiNew('audio', in, dur).madd(mul, add)
	}
}

PitchTracker : UGen {
	*ar { arg in = 0.0;
		^this.multiNew('audio', in)
	}
}
