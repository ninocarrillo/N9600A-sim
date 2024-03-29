clear *;
pkg load signal;
subplot_count = 3;
[wavetable, dsr] = audioread("Wavetable.wav", "native");
[audiosine, dsr] = audioread("AudioSine.wav", "native");
[audiocosine, dsr] = audioread("AudioCosine.wav", "native");
[audiodither, dsr] = audioread("AudioDither.wav", "native");
subplot(subplot_count,1,1);
plot(wavetable);
subplot(subplot_count,1,2);
plot(audiosine);
hold on;
plot(audiocosine);
hold off;
subplot(subplot_count, 1,3);
plot(audiodither);