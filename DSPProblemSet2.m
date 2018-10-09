% Michael Lendino ECE310 DSP PSET 2 MATLAB 
clc;
clear all;

%% Compute B and w0 and find the specifications of the equivalent lowpass filter
Apass = 2;
Astop = 30;
Fpass = [100e3, 130e3];
Fstop = [90e3, 140e3];
wpass = 2*pi*Fpass;
wstop = 2*pi*Fstop;

B = wpass(2)- wpass(1);
w0 = sqrt(wpass(2)*wpass(1));

wproto90 = ((wstop(1)).^2 - (w0).^2)./(B*wstop(1));
wproto140 = ((wstop(2)).^2 - (w0).^2)./(B*wstop(2));
%Since the magnitude of wproto140 is less than the magnitude of wproto90,
%pick wproto140
%% Use the explicit formulas for the Butterworth and Chebyshev filter orders for the lowpass filters
%For Butterworth, formula given in notes
nButter = ceil((log10((10^(Astop/10)-1)/(10^(Apass/10)-1))/log10(wproto140/1))/2);
%For Chebyshev, formula given in notes
nCheby = ceil(acosh(sqrt(10^(Astop/10)-1)/(10^(Apass/10)-1))/acosh(wproto140/1));

%% Call the appropriate functions to design the bandpass filters and check that the orders match
[n,Wn] = buttord([wpass],[wstop],Apass,Astop,'s');
[z,p,k] = butter(n, Wn, 's');
[n2,Wn2] = cheb1ord([wpass],[wstop],Apass,Astop,'s');
[z2,p2,k2] = cheby1(n2,Apass,Wn2,'s');
[n3,Wn3] = cheb2ord([wpass],[wstop],Apass,Astop,'s');
[z3,p3,k3] = cheby2(n3,Astop,Wn3,'s');
[n4,Wn4] = ellipord([wpass],[wstop],Apass,Astop,'s');
[z4,p4,k4] = ellip(n4,Apass,Astop,Wn4,'s'); 
%since nButter=9=n and nCheby=5=n2=n3, the orders are confirmed

%% Obtain zplane plots for each of the filters
zplane(z,p);
title('Pole Zero Plot for the Butterworth Bandpass Filter');
grid on;
figure;
zplane(z2,p2);
title('Pole Zero Plot for the Chebyshev I Bandpass Filter');
grid on;
figure;
zplane(z3,p3);
title('Pole Zero Plot for the Chebyshev II Bandpass Filter');
grid on;
figure;
zplane(z4,p3);
title('Pole Zero Plot for the Elliptic Bandpass Filter');
grid on;

%% Design the elliptic lowpass prototype filter
wProtoP = 1;
[n5,Wn5] = ellipord(wProtoP,wproto140,Apass,Astop,'s');
[z5,p5,k5] = ellip(n5,Apass,Astop,Wn5,'s');
figure;
zplane(z5,p5);
title('Poles Zero Plot for Elliptic Lowpass Prototype Filter');
grid on;

%% Obtain magnitude and response plots for the four bandpass filters; on the graphs, superimpose dashed line segments to indicate the magnitude specifications
%need to make these to superimpose dashed line segments in correct region 
StopLowerBound = linspace(0, 90);
passBounds = linspace(100, 130);
StopUpperBound = linspace(140, 200);

[b,a] = zp2tf(z,p,k);
[b2,a2] = zp2tf(z2,p2,k2);
[b3,a3] = zp2tf(z3,p3,k3);
[b4,a4] = zp2tf(z4,p4,k4);
%graphs were jagged before, thank you Ali for telling me to add more
%sampling points in the third coordinate
[h,w] = freqs(b,a,100000);
[h2,w2] = freqs(b2,a2,100000);
[h3,w3] = freqs(b3,a3,100000);
[h4,w4] = freqs(b4,a4,100000);

%need to graph the magnitude in dB
Butterdb = 20*log10(abs(h));
%phase in degrees unwrapped
ButterPhase = unwrap(angle(h)) * (180/(pi));
figure;
grid on;
subplot(2,1,1);
hold on;
plot(w/(2e3*pi),Butterdb);
%placing the dashed lines
plot(StopLowerBound, -30*ones(size(StopLowerBound)), 'k--');
plot(passBounds, -2*ones(size(passBounds)), 'k--');
plot(passBounds, 0*ones(size(passBounds)), 'k--');
plot(StopUpperBound, -30*ones(size(StopUpperBound)), 'k--');
hold off;
ylabel('Magnitude (dB)');
xlabel('Frequency (kHz)');
title('Magnitude Frequency Response for Butterworth Bandpass');
grid on;
%frequency on the x axis from 0 to 200kHz on a linear scale
xlim([0 200]);
%magnitude in dB
ylim([-40 5]);
subplot(2,1,2);
plot(w/(2e3*pi), ButterPhase);
ylabel('Phase (Degrees)');
xlabel('Frequency (kHz)');
xlim([0 200]);
title('Phase Frequency Response for Butterworth Bandpass');
grid on;

%need to graph the magnitude in dB
ChebyIdb = 20*log10(abs(h2));
%phase in degrees unwrapped
ChebyIPhase = unwrap(angle(h2)) * (180/(pi));
figure;
grid on;
subplot(2,1,1);
hold on;
plot(w2/(2e3*pi),ChebyIdb);
%placing the dashed lines
plot(StopLowerBound, -30*ones(size(StopLowerBound)), 'k--');
plot(passBounds, -2*ones(size(passBounds)), 'k--');
plot(passBounds, 0*ones(size(passBounds)), 'k--');
plot(StopUpperBound, -30*ones(size(StopUpperBound)), 'k--');
hold off;
ylabel('Magnitude (dB)');
xlabel('Frequency (kHz)');
title('Magnitude Frequency Response for Chebyshev I Bandpass');
grid on;
%frequency on the x axis from 0 to 200kHz on a linear scale
xlim([0 200]);
%magnitude in dB
ylim([-40 5]);
subplot(2,1,2);
plot(w2/(2e3*pi), ChebyIPhase);
ylabel('Phase (Degrees)');
xlabel('Frequency (kHz)');
xlim([0 200]);
title('Phase Frequency Response for Chebyshev I Bandpass');
grid on;

%need to graph the magnitude in dB
ChebyIIdb = 20*log10(abs(h3));
%phase in degrees unwrapped
ChebyIIPhase = unwrap(angle(h3)) * (180/(pi));
figure;
grid on;
subplot(2,1,1);
hold on;
plot(w3/(2e3*pi),ChebyIIdb);
%placing the dashed lines
plot(StopLowerBound, -30*ones(size(StopLowerBound)), 'k--');
plot(passBounds, -2*ones(size(passBounds)), 'k--');
plot(passBounds, 0*ones(size(passBounds)), 'k--');
plot(StopUpperBound, -30*ones(size(StopUpperBound)), 'k--');
hold off;
ylabel('Magnitude (dB)');
xlabel('Frequency (kHz)');
title('Magnitude Frequency Response for Chebyshev II Bandpass');
grid on;
%frequency on the x axis from 0 to 200kHz on a linear scale
xlim([0 200]);
%magnitude in dB
ylim([-40 5]);
subplot(2,1,2);
plot(w3/(2e3*pi), ChebyIIPhase);
ylabel('Phase (Degrees)');
xlabel('Frequency (kHz)');
xlim([0 200]);
title('Phase Frequency Response for Chebyshev II Bandpass');
grid on;

%need to graph the magnitude in dB
Ellipticdb = 20*log10(abs(h4));
%phase in degrees unwrapped
EllipticPhase = unwrap(angle(h4)) * (180/(pi));
figure;
grid on;
subplot(2,1,1);
hold on;
plot(w4/(2e3*pi),Ellipticdb);
%placing the dashed lines
plot(StopLowerBound, -30*ones(size(StopLowerBound)), 'k--');
plot(passBounds, -2*ones(size(passBounds)), 'k--');
plot(passBounds, 0*ones(size(passBounds)), 'k--');
plot(StopUpperBound, -30*ones(size(StopUpperBound)), 'k--');
hold off;
ylabel('Magnitude (dB)');
xlabel('Frequency (kHz)');
title('Magnitude Frequency Response of Elliptic Bandpass');
grid on;
%frequency on the x axis from 0 to 200kHz on a linear scale
xlim([0 200]);
%magnitude in dB
ylim([-40 5]);
subplot(2,1,2);
plot(w4/(2e3*pi), EllipticPhase);
ylabel('Phase (Degrees)');
xlabel('Frequency (kHz)');
xlim([0 200]);
title('Phase Frequency Response of Elliptic Bandpass');
grid on;

%% Obtain zoomed in graphs  of the magnitude response of the four filters in the passband all superimposed on one plot
figure;
hold on;
plot(w/(2e3*pi), Butterdb);
plot(w2/(2e3*pi), ChebyIdb);
plot(w3/(2e3*pi), ChebyIIdb);
plot(w4/(2e3*pi), Ellipticdb);
%superimposing the lines at the correct attenuation
plot(linspace(95, 135), -2*ones(100), 'k--');
plot(linspace(95, 135), 0*ones(100), 'k--');
hold off;
ylabel('Magnitude (dB)');
xlabel('Frequency (kHz)');
legend({'Butterworth', 'Chebyshev I', 'Chebyshev II', 'Elliptic'}, 'Location', 'South');
title('Zoomed in Magnitude Response of the Four Filters in the Passband');
grid on;
xlim([95 135]);
ylim([-3 1]);

%% Compute the attenuation at each of the stopband edges for each of the filters listed in a table
%Thank you to Leart and Ali for telling me about the existence and functionality of
%interpl()
ButterStopLower = interp1(w/(2e3*pi), Butterdb, 90);
ButterStopUpper = interp1(w/(2e3*pi), Butterdb, 140);

ChebyIStopLower = interp1(w2/(2e3*pi), ChebyIdb, 90);
ChebyIStopUpper = interp1(w2/(2e3*pi), ChebyIdb, 140);

ChebyIIStopLower = interp1(w3/(2e3*pi), ChebyIIdb, 90);
ChebyIIStopUpper = interp1(w3/(2e3*pi), ChebyIIdb, 140);

EllipticStopLower = interp1(w4/(2e3*pi), Ellipticdb, 90);
EllipticStopUpper = interp1(w4/(2e3*pi), Ellipticdb, 140);

attenuationTable = [ButterStopLower, ButterStopUpper; ChebyIStopLower, ChebyIStopUpper; ChebyIIStopLower, ChebyIIStopUpper; EllipticStopLower, EllipticStopUpper];
