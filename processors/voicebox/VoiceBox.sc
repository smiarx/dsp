VoiceBox : UGen {
	*ar { arg in = 0.0, shift=1.0, formantShift=1.0, mul=1.0, add=0.0;
		^this.multiNew('audio', in, shift, formantShift).madd(mul, add)
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
