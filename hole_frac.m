function waveform = hole(r,r0,w,s,depthfrac) 

waveform = ...
    1 * ...
    ( ...
        1 ...
        + depthfrac * tanh( s * ( r - r0 - w/2 ) ) ...
        - depthfrac * tanh( s * ( r - r0 + w/2 ) ) ...
    );