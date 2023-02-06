clear *;
pkg load signal;
subplot_count = 2;
[demod_sig, demod_fs] = audioread("DemodSignal.wav");
[filtered_sig, filtered_fs] = audioread("FilteredSignal.wav");
[phase_acc, phase_fs] = audioread("PhaseAccumulator.wav");
[PLL_control, pll_fs] = audioread("PLLControl.wav","native");
%subplot(subplot_count,1,1);
%plot(filtered_sig);
subplot(subplot_count, 1,1);
plot(demod_sig);
hold on;
%subplot(subplot_count, 1,3);
plot(phase_acc);
hold off;
subplot(subplot_count, 1,2);
plot(PLL_control);

