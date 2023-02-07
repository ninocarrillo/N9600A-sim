clear *;
pkg load signal;
subplot_count = 2;
[wavetable, dsr] = audioread("Wavetable.wav", "native");
[audio, dsr] = audioread("Audio.wav", "native");
subplot(subplot_count,1,1);
plot(wavetable);
subplot(subplot_count,1,2);
plot(audio);
