function wave = pulse(t,freq,duty)
    tperiod = 1/freq;
    if ( (t/tperiod - floor(t/tperiod)) < duty )
        wave = 1;
    else
        wave = 0;
    end
end