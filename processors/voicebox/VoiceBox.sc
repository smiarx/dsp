VoiceBox : UGen {
	*ar { arg in = 0.0, lambda=0.995, warpin=0, warpout=0, mul=1.0, add=0.0;
		^this.multiNew('audio', in, lambda, warpin,warpout).madd(mul, add)
	}
}
