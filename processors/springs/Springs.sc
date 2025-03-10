Springs : MultiOutUGen
{
    *ar
    {
        arg inputArray, r=0.75, freq=2500, td=0.05, t60=3.5, diffusion=0, chaos=0.0, spread=1, width=1, drywet=1;
        ^this.multiNewList(['audio',r,freq,td,t60,diffusion,chaos,spread,width,drywet]++ inputArray.asArray);
    }

    checkInputs
    {
        var numArgs        = 9;
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
