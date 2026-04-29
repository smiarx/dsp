VoiceBox : UGen {
	*ar { arg in = 0.0, shift=1.0, formantShift=1.0, quantizeSpeed=0.02, quantizeConfirmTime=0.002, a4=440, scale=(0..11), mul=1.0, add=0.0;
		^this.multiNew('audio', in, shift, formantShift, quantizeSpeed, quantizeConfirmTime, a4, (1<<(scale%12)).reduce('|')).madd(mul, add)
	}
}

ZC : UGen {
	*ar { arg in = 0.0, dur=0.05, mul=1.0, add=0.0;
		^this.multiNew('audio', in, dur).madd(mul, add)
	}
}

PitchTracker : MultiOutUGen
{
    *ar
    {
		arg in;
        ^this.multiNew('audio', in);
    }

    init
    {
        arg... theInputs;
        inputs = theInputs;
        ^this.initOutputs(2, rate);
    }
}

HMM : UGen {
	*ar { arg in = 0.0, transition=0.5, mul=1.0, add=0.0;
		^this.multiNew('audio', in, transition).madd(mul, add)
	}
}
