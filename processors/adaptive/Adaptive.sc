Adaptive {
    var <buffer, <order;

    *new { arg order, server;
        var buffer, numFrames;
        order = order ? 8;
        numFrames = (order + 4 - 1);
        numFrames = numFrames - (numFrames%4);
        buffer = Buffer.alloc(server, numFrames);
		buffer.zero;
        ^super.newCopyArgs(buffer, order);
    }

    free {
        buffer.free();
    }
}

RLS : Filter {
	*ar { arg in = 0.0, adaptive, mul=1.0, add=0.0;
		^this.multiNew('audio', in, adaptive.order, adaptive.buffer).madd(mul, add)
	}
}

AdaptiveReconstruct : Filter {
	*ar { arg in = 0.0, adaptive, mul=1.0, add=0.0;
		^this.multiNew('audio', in, adaptive.order, adaptive.buffer).madd(mul, add)
	}
}
