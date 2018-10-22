% Michael Lendino ECE310 DSP PSET 3
%sorry for not doing this in a very efficient way, I did this late at night
%and was tired
clc;
clear all;

%% Call the appropriate functions to design digital bandpass filters with a sampling rate of 400kHz
Apass = 2;
Astop = 30;
Fpass = [100e3, 130e3];
Fstop = [90e3, 140e3];
%normalize frequency to sampling rate
Fsamp = 400e3;
FpassBL = [100e3, 130e3]/(Fsamp/2);
FstopBL = [90e3, 140e3]/(Fsamp/2);
wpass = 2*pi*Fpass;
wstop = 2*pi*Fstop;

%% Call the appropriate functions to design the bandpass filters and check that the orders match
%First we have to design the analog filters
[n,Wn] = buttord([wpass],[wstop],Apass,Astop,'s');
[b1,a1] = butter(n, Wn, 's');
[zb,pb,kb] = butter(n, Wn, 's');
[n2,Wn2] = cheb1ord([wpass],[wstop],Apass,Astop,'s');
[z2,p2,k2] = cheby1(n2,Apass,Wn2,'s');
[b2,a2] = cheby1(n2,Apass,Wn2,'s');
[n3,Wn3] = cheb2ord([wpass],[wstop],Apass,Astop,'s');
[z3,p3,k3] = cheby2(n3,Astop,Wn3,'s');
[b3,a3] = cheby2(n3,Astop,Wn3,'s');
[n4,Wn4] = ellipord([wpass],[wstop],Apass,Astop,'s');
[z4,p4,k4] = ellip(n4,Apass,Astop,Wn4,'s'); 
[b4,a4] = ellip(n4,Apass,Astop,Wn4,'s'); 
%Order for each analog filter
orderButterA = 2*n;
orderChebyIA = 2*n2;
orderChebyIIA = 2*n3;
orderEllipA = 2*n4;
%since nButter=9=n and nCheby=5=n2=n3, the orders are confirmed
%Now we have to do digital via bilinear transform
[n1,Wn1] = buttord([FpassBL],[FstopBL],Apass,Astop);
[z1,p1,k1] = butter(n1, Wn1); 
[b5,a5] = butter(n1, Wn1);
[n22,Wn22] = cheb1ord([FpassBL],[FstopBL],Apass,Astop);
[z22,p22,k22] = cheby1(n22,Apass,FpassBL);
[b6,a6] = cheby1(n22,Apass,FpassBL);
[n33,Wn33] = cheb2ord([FpassBL],[FstopBL],Apass,Astop);
[z33,p33,k33] = cheby2(n33,Astop,FstopBL);
[b7,a7] = cheby2(n33,Astop,FstopBL);
[n44,Wn44] = ellipord([FpassBL],[FstopBL],Apass,Astop);
[z44,p44,k44] = ellip(n44,Apass,Astop,Wn44); 
[b77,a77] = ellip(n44,Apass,Astop,Wn44); 
%order for each digital via bilinear transform filter
orderButterBL = 2*n1;
orderChebyIBL = 2*n22;
orderChebyIIBL = 2*n33;
orderEllipBL = 2*n44;

%Now we have to do design via impulse invariance
[bz1,az1] = impinvar(b1,a1,Fsamp);
[zz1,pp1,kk1] = tf2zp(bz1,az1);
[bz2,az2] = impinvar(b2,a2,Fsamp);
[zz2,pp2,kk2] = tf2zp(bz2,az2);
[bz3,az3] = impinvar(b3,a3,Fsamp);
[zz3,pp3,kk3] = tf2zp(bz3,az3);
[bz4,az4] = impinvar(b4,a4,Fsamp);
[zz4,pp4,kk4] = tf2zp(bz4,az4);
%order for the impulse invariance filters
Ni1 = filtord(bz1,az1);
Ni2 = filtord(bz2,az2);
Ni3 = filtord(bz3,az3);
Ni4 = filtord(bz4,az4);

%get pole zero plots for the analog, digital design via bilinear transform
%and for the digital via impulse invariance
figure('Name','Analog Filters','NumberTitle','off'); 
subplot(2,2,1);
zplane(zb,pb);
grid on;
title('Pole Zero Plot for the Butterworth Bandpass Filter');
subplot(2,2,2);
zplane(z2,p2);
grid on;
title('Pole Zero Plot for the Chebyshev I Bandpass Filter');
subplot(2,2,3);
zplane(z3,p3);
grid on;
title('Pole Zero Plot for the Chebyshev II Bandpass Filter');
subplot(2,2,4);
zplane(z4,p4);
grid on;
title('Pole Zero Plot for the Elliptic Bandpass Filter');

figure('Name','Bilinear Transforms','NumberTitle','off');
subplot(2,2,1);
zplane(z1,p1);
grid on;
title('Digital: Butter via Bilinear Transform');
subplot(2,2,2);
zplane(z22,p22);
grid on;
title('Digital: Chebyshev I via Bilinear Transform');
subplot(2,2,3);
zplane(z33,p33);
grid on;
title('Digital: Chebyshev II via Bilinear Transform');
subplot(2,2,4);
zplane(z44,p44);
grid on;
title('Digital: Elliptic via Bilinear Transform');

figure('Name','Impulse Invariance Design','NumberTitle','off');
subplot(2,2,1);
zplane(zz1,pp1);
grid on;
title('Digital: Butter via Impulse Invariance');
subplot(2,2,2);
zplane(zz2,pp2);
grid on;
title('Digital: Chebyshev I via Impulse Invariance');
subplot(2,2,3);
zplane(zz3,pp3);
grid on;
title('Digital: Chebyshev II via Impulse Invariance');
subplot(2,2,4);
zplane(zz4,pp4);
grid on;
title('Digital: Elliptic via Impulse Invariance');

%get magnitude response plots
[b111,a111] = zp2tf(zb,pb,kb);
[b21,a21] = zp2tf(z2,p2,k2);
[b31,a31] = zp2tf(z3,p3,k3);
[b41,a41] = zp2tf(z4,p4,k4);

[b51,a51] = zp2tf(z1,p1,k1);
[b61,a61] = zp2tf(z22,p22,k22);
[b71,a71] = zp2tf(z33,p33,k33);
[b771,a771] = zp2tf(z44,p44,k44);

[b81,a81] = zp2tf(zz1,pp1,kk1);
[b91,a91] = zp2tf(zz2,pp2,kk2);
[b101,a101] = zp2tf(zz3,pp3,kk3);
[b1111,a1111] = zp2tf(zz4,pp4,kk4);

[h,w] = freqs(b111,a111,100000);
[h2,w2] = freqs(b21,a21,100000);
[h3,w3] = freqs(b31,a31,100000);
[h4,w4] = freqs(b41,a41,100000);

[h5,w5] = freqz(b51,a51,100000, Fsamp);
[h6,w6] = freqz(b61,a61,100000, Fsamp);
[h7,w7] = freqz(b71,a71,100000, Fsamp);
[h77,w77] = freqz(b771,a771,100000, Fsamp);

[h8,w8] = freqz(b81,a81,100000, Fsamp);
[h9,w9] = freqz(b91,a91,100000, Fsamp);
[h10,w10] = freqz(b101,a101,100000, Fsamp);
[h11,w11] = freqz(b1111,a1111,100000, Fsamp);

%need to graph the magnitude in dB
Butterdb = 20*log10(abs(h));
ChebyIdb = 20*log10(abs(h2));
ChebyIIdb = 20*log10(abs(h3));
Ellipticdb = 20*log10(abs(h4));

ButterBLdb = 20*log10(abs(h5));
ChebyIBLdb = 20*log10(abs(h6));
ChebyIIBLdb = 20*log10(abs(h7));
EllipticBLdb = 20*log10(abs(h77));

ButterIMdb = 20*log10(abs(h8));
ChebyIdbIM = 20*log10(abs(h9));
ChebyIIdbIM = 20*log10(abs(h10));
EllipticIMdb = 20*log10(abs(h11));

figure('Name','Analog Magnitude Response','NumberTitle','off');
subplot(2,2,1);
plot(w/(2*pi),Butterdb);
grid on;
title('Analog: Butterworth Magnitude Response');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
%frequency on the x axis from 0 to Nyquist Bandwidth on a linear scale
xlim([0 (Fsamp/2)]);
%magnitude in dB
ylim([-50 1]);
subplot(2,2,2);
plot(w2/(2*pi),ChebyIdb);
grid on;
title('Analog: Chebyshev I Magnitude Response');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
%frequency on the x axis from 0 to Nyquist Bandwidth on a linear scale
xlim([0 (Fsamp/2)]);
%magnitude in dB
ylim([-50 1]);
subplot(2,2,3);
plot(w3/(2*pi),ChebyIIdb);
grid on;
title('Analog: Chebyshev II Magnitude Response');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
%frequency on the x axis from 0 to Nyquist Bandwidth on a linear scale
xlim([0 (Fsamp/2)]);
%magnitude in dB
ylim([-50 1]);
subplot(2,2,4);
plot(w4/(2*pi),Ellipticdb);
grid on;
title('Analog: Chebyshev II Magnitude Response');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
%frequency on the x axis from 0 to Nyquist Bandwidth on a linear scale
xlim([0 (Fsamp/2)]);
%magnitude in dB
ylim([-50 1]);

figure('Name','Digital Design via Bilinear Transform Magnitude Response','NumberTitle','off');
subplot(2,2,1);
plot(w5,ButterBLdb);
grid on;
title('Digital via Bilinear: Butterworth Magnitude Response');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
%frequency on the x axis from 0 to Nyquist Bandwidth on a linear scale
xlim([0 (Fsamp/2)]);
%magnitude in dB
ylim([-50 1]);
subplot(2,2,2);
plot(w6,ChebyIBLdb);
grid on;
title('Digital via Bilinear: Chebyshev I Magnitude Response');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
%frequency on the x axis from 0 to Nyquist Bandwidth on a linear scale
xlim([0 (Fsamp/2)]);
%magnitude in dB
ylim([-50 1]);
subplot(2,2,3);
plot(w7,ChebyIIBLdb);
grid on;
title('Digital via Bilinear: Chebyshev II Magnitude Response');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
%frequency on the x axis from 0 to Nyquist Bandwidth on a linear scale
xlim([0 (Fsamp/2)]);
%magnitude in dB
ylim([-50 1]);
subplot(2,2,4);
plot(w77,EllipticBLdb);
grid on;
title('Digital via Bilinear: Elliptic Magnitude Response');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
%frequency on the x axis from 0 to Nyquist Bandwidth on a linear scale
xlim([0 (Fsamp/2)]);
%magnitude in dB
ylim([-50 1]);

figure('Name','Digital Design via Impulse Invariance Magnitude Response','NumberTitle','off');
subplot(2,2,1);
plot(w8,ButterIMdb);
grid on;
title('Digital via Impulse Invariance: Butterworth Magnitude Response');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
%frequency on the x axis from 0 to Nyquist Bandwidth on a linear scale
xlim([0 (Fsamp/2)]);
%magnitude in dB
ylim([-50 1]);
subplot(2,2,2);
plot(w9,ChebyIdbIM);
grid on;
title('Digital via Impulse Invariance: Chebyshev I Magnitude Response');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
%frequency on the x axis from 0 to Nyquist Bandwidth on a linear scale
xlim([0 (Fsamp/2)]);
%magnitude in dB
ylim([-50 1]);
subplot(2,2,3);
plot(w10,ChebyIIdbIM);
grid on;
title('Digital via Impulse Invariance: Chebyshev II Magnitude Response');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
%frequency on the x axis from 0 to Nyquist Bandwidth on a linear scale
xlim([0 (Fsamp/2)]);
%magnitude in dB
ylim([-50 1]);
subplot(2,2,4);
plot(w11,EllipticIMdb);
grid on;
title('Digital via Impulse Invariance: Elliptic Magnitude Response');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
%frequency on the x axis from 0 to Nyquist Bandwidth on a linear scale
xlim([0 (Fsamp/2)]);
%magnitude in dB
ylim([-50 1]);

%now we wish to find the time constants for each of the 12 filters
%for the four analog filters, butterworth, ChebyI, Cheby II, and elliptic
tau1 = max(1./abs(real(pb)));
tau2 = max(1./abs(real(p2)));
tau3 = max(1./abs(real(p3)));
tau4 = max(1./abs(real(p4)));
%for the eight digital filters, four bilinear, four impulse invariant,
%butterworth, chebyI, chebyII, and elliptic
tau5 = max(abs(1./log(abs(p1))))/Fsamp;
tau6 = max(abs(1./log(abs(p22))))/Fsamp;
tau7 = max(abs(1./log(abs(p33))))/Fsamp;
tau8 = max(abs(1./log(abs(p44))))/Fsamp;

tau9 = max(abs(1./log(abs(pp1))))/Fsamp;
tau10 = max(abs(1./log(abs(pp2))))/Fsamp;
tau11 = max(abs(1./log(abs(pp3))))/Fsamp;
tau12 = max(abs(1./log(abs(pp4))))/Fsamp;
%Now we wish to calculate the specifications of the analog bandpass filter
%obtained via prewarping
Omega1 = tan(0.5*(90e3*2*pi/Fsamp));
Omega2 = tan(0.5*(100e3*2*pi/Fsamp));
Omega3 = tan(0.5*(130e3*2*pi/Fsamp));
Omega4 = tan(0.5*(140e3*2*pi/Fsamp));





