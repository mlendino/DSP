% Michael Lendino ECE310 DSP PSET 1 Problem 2 MATLAB 
clc;
clear all;
%% Superimpose graphs of the phase response of Hmin(z) H(z) and A(z) all unwrapped in degrees
% zeros in the numerator and poles in the denominator of H(z)
Hzero = [-4/3 -5 inf]';
Hpole = [1/2 -3 -3]';
% zeros in the numerator and poles in the denominator of A(z)
Azero = [-1/5 -3 -3 -3/4]';
Apole = [-5 -1/3 -1/3 -4/3]';
% zeros in the numerator and poles in the denominator of Hmin(z)
Hminzero = [-1/5 -3/4 inf 0]';
Hminpole = [1/2 -1/3 -1/3]';
%Thank you to Ali Rahman for helping me figure out what sampling rate to
%choose for phasez() below
%create the transfer function for each set of poles and zeros, and create
%the unwrapped phase response 
[Hnum, Hdenom] = zp2tf(Hzero, Hpole, 3/2);
[Hphase, wH] = phasez(Hnum, Hdenom, 10000);

[Anum, Adenom] = zp2tf(Azero, Apole, 20/27);
[Aphase, wA] = phasez(Anum, Adenom, 10000);

[Hminnum, Hmindenom] = zp2tf(Hminzero, Hminpole, 10/9);
[Hminphase, wHm] = phasez(Hminnum, Hmindenom, 10000);

figure; 
grid on;
plot(wH, Hphase * 180/pi);
hold on;
plot(wA, Aphase * 180/pi);
hold on;
plot(wHm, Hminphase * 180/pi);
legend('H', 'A', 'H_(min)');
title('Phase Response (Unwrapped in Degrees) of a filter H(z), an allpass filter A(z), and minimum-phase H_(min)(z)');
xlabel('Frequency (Hz)');
ylabel('Phase Response (Unwrapped in Degrees)'); 