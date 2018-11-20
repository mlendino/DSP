% Michael Lendino ECE310 DSP PSET 6  
clc;
clear all;
%% Comparing a Chebyshev window and a Kaiser window
%Construct a length 31 Chebyshev window with 30dB peak sidelobe level. 
% Normalize the window so its DC value is 1.
chebywindow = chebwin(31,30);

[f,w] = freqz(chebywindow, 1000);
fcheby = f/f(1);
%Use zplane to look at the zeros of the window. By calculating the zeros, cmpute
%the (null-to-null) mainlobe width relative to that of the rectangular window (i.e.,
%divide by 4pi/N).
[z, p, k] = tf2zpk(chebywindow);
figure('Name','Zeros of the Chebyshev Window','NumberTitle','off');
zplane(z,p)
title('Zeros of the Chebyshev Window');
MaincalcWidth = abs(2*(z(1))/((4*pi)/31));
%null to null mainlobe width is 4.9338
%Graph the magnitude response of the Chebyshev window, with magnitude on a
%decibel scale.
magcheb = 20*log10(abs(fcheby))
figure('Name','Magnitude Response of the Chebyshev Window','NumberTitle','off');
plot(w/pi, magcheb)
title('Superimposed Kaiser Window on Chebyshev Window Magnitude Response');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
%By trial and error find the beta so a Kaiser window of the same length has
%the same mainlobe width ((superimpose a graph of jH (!)j for the kaiser
%window until the first null lines up)
kaiserr = kaiser(31,3.02)
[h,w1] = freqz(kaiserr, 1000);
fkaiser = h/h(1);

[z1,p1,k1] = tf2zpk(kaiserr)
magkaiser = 20*log10(abs(fkaiser));
hold on
plot(w1/pi,magkaiser)
legend('Magnitude Response of the Chebyshev Window', 'Magnitude Response of the Kaiser Window');
hold off
%zdifference = abs(z(1)) - abs(z1(1))
figure('Name','Superimposed Zeros of Chebyshev and Kaiser Windows','NumberTitle','off');
zplane(z,z1)
legend('Chebyshev zeros', 'Kaiser zeros');

%beta is 3.02

%wvtool(kaiserr);
%peak sidelobe level is 24.8dB via wvtool in the above comment, its
%commented because it's annoying when it pops up

%Compute the fraction of the energy in the sidelobes for each window
%(beyond the first null)
EnergyTotalCheby = sum(abs(fcheby).^2);
zeroCheb = islocalmin(abs(fcheby));
EnergyLobeCheby = sum(abs(fcheby(49:end)).^2);
chebFrac = EnergyLobeCheby/EnergyTotalCheby
%do similar thing for kaiser
EnergyTotalKaiser = sum(abs(fkaiser).^2);
zeroKaiser = islocalmin(abs(fkaiser));
EnergyLobeKaiser = sum(abs(fkaiser(49:end)).^2);
Kaiserrac = EnergyLobeKaiser/EnergyTotalCheby

%% We wish to design a digital bandpass filter meeting certain specifications
%Compute the deviations for the bands
deltaStop = abs(10^(-1.5))
deltaPass = abs((10^(-0.1) - 1)/(10^(-0.1) + 1))

%Use MATLAB tools to estimate the filter order, and design the filters.
%From the documentation, we have that:
fsamp = 10e6;
fcuts = [1.5e6 2e6 3e6  3.5e6];
mags = [0 1 0];
devs = [0.0316 0.1146 0.0316];
%Obtaining Kaiser
[n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fsamp);
%The filter order is 31
hh = fir1(n+1,Wn,ftype,kaiser(n+2,beta),'noscale');

[H,f] = freqz(hh,1,1024,fsamp);

%For Parks McClellan
[n,f,a,wf] = firpmord(fcuts,mags,devs, fsamp)
%The filter order is 19
b = firpm(n+4,f,a,wf);
b = [b zeros(1,abs(length(b) - length(hh)))]
[h,wff] = freqz(b,1,1024, fsamp);
%However Armaan asked you in class why it was underdesigned and you
%suggested to try n+1 and n+2 but they still didn't work if I did this
%correctly, so the new filter orders are 23 and 32

figure('Name','Chebyshev and Kaiser Stem Plots','NumberTitle','off');
subplot(2,1,1)
stem(1:length(b), b)
title('Stem Plot for Parks McClellan Filter');
xlabel('n');
ylabel('h[n]');
subplot(2,1,2)
stem(1:length(hh), hh)
title('Stem Plot for Kaiser Filter');
xlabel('n');
ylabel('h[n]');

StopLowerBound = linspace(0, 1.5e6);
passBounds = linspace(2e6, 3e6);
StopUpperBound = linspace(3.5e6, 5e6);

figure('Name','Magnitude Response of Parks McClellan and Kaiser Filters','NumberTitle','off');
subplot(2,1,1)
plot(wff, 20*log10(abs(H)));
title('Magnitude Response for Kaiser Filter');
xlabel('Frequency (Hz)');
ylabel('dB');
hold on
ylim([-50 2]);
plot(StopLowerBound, -30*ones(size(StopLowerBound)), 'k--');
plot(passBounds, -2*ones(size(passBounds)), 'k--');
plot(passBounds, 0*ones(size(passBounds)), 'k--');
plot(StopUpperBound, -30*ones(size(StopUpperBound)), 'k--');
hold off 
subplot(2,1,2)
plot(wff, 20*log10(abs(h)));
title('Magnitude Response for Parks McClellan Filter');
xlabel('Frequency (Hz)');
ylabel('dB');
hold on
ylim([-50 2]);
plot(StopLowerBound, -30*ones(size(StopLowerBound)), 'k--');
plot(passBounds, -2*ones(size(passBounds)), 'k--');
plot(passBounds, 0*ones(size(passBounds)), 'k--');
plot(StopUpperBound, -30*ones(size(StopUpperBound)), 'k--');
hold off

%Theyre underdesigned, but the order was modified slightly to have them better fit the specs 

%This quantity matches w, the weights, thus we confirm the property in part
%e.
doublecheck = deltaPass./deltaStop
