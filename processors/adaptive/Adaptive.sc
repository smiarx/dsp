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
	*ar { arg in = 0.0, adaptive, lambda=0.995, mul=1.0, add=0.0;
		^this.multiNew('audio', in, adaptive.order, adaptive.buffer, lambda).madd(mul, add)
	}
}

RLSWarped : Filter {
	*ar { arg in = 0.0, adaptive, lambda=0.995, warp=0, mul=1.0, add=0.0;
		^this.multiNew('audio', in, adaptive.order, adaptive.buffer, lambda, warp).madd(mul, add)
	}
}

AdaptiveReconstruct : Filter {
	*ar { arg in = 0.0, adaptive, mul=1.0, add=0.0;
		^this.multiNew('audio', in, adaptive.order, adaptive.buffer).madd(mul, add)
	}
}

AdaptiveReconstructWarped : Filter {
	*ar { arg in = 0.0, adaptive, warp=0, mul=1.0, add=0.0;
		^this.multiNew('audio', in, adaptive.order, adaptive.buffer, warp).madd(mul, add)
	}
}
