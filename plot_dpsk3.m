clear *;
pkg load signal;
subplot_count = 4;
[loopfilter, dsr] = audioread("LoopFilter.wav", "native");
[datafilter, dsr] = audioread("DataFilter.wav", "native");
[thirdmixer, dsr] = audioread("ThirdMixer.wav", "native");
[secondmixer, dsr] = audioread("SecondMixer.wav", "native");
[firstmixer, dsr] = audioread("FirstMixer.wav", "native");
[phaseaccumulator, dsr] = audioread("PhaseAccumulator.wav", "native");
[samplepulse, dsr] = audioread("SamplePulse.wav", "native");
[filteredsignal, dsr] = audioread("FilteredSignal.wav", "native");

subplot(subplot_count,1,1);
plot(loopfilter);
hold on;
plot(firstmixer);
hold off;
subplot(subplot_count,1,2);
plot(datafilter);
hold on;
plot(samplepulse, '-o');
hold off;
subplot(subplot_count, 1,3);
plot(thirdmixer);

subplot(subplot_count, 1,4);
plot(filteredsignal);
hold on;
plot(samplepulse, '-o')
hold off;


