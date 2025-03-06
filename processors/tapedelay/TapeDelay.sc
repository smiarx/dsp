TapeDelay : MultiOutUGen
{
    *ar
    {
		arg inputArray, delay=0.1, feedback=0.9, lpf=2000, hpf=120, saturation=(-10), drift=0.0, mode=0, drywet=1.0;
        ^this.multiNewList(['audio', delay, feedback, lpf, hpf, saturation, drift, mode, drywet]++ inputArray.asArray);
    }

    checkInputs
    {
        var numArgs        = 8;
        var numAudioInputs = this.numInputs - numArgs;

        if (numAudioInputs != 2, { ^("input is not stereo"); })
            ;
        numAudioInputs.do({
            | i | var index = numArgs + i;
            if (inputs.at(index).rate != 'audio', {
                    ^("input is not audio rate:" + inputs.at(index) +
                      inputs.at(index).rate);
                })
                ;
        });
        ^this.checkValidInputs;
    }

    init
    {
        arg... theInputs;
        inputs = theInputs;
        ^this.initOutputs(2, rate);
    }
}
