% Michael Lendino ECE310 DSP PSET 4 MATLAB 
clc;
clear all;
%% 2: Compute the straddle loss by evaluating the Dirichlet sinc at the frequency offset
%part c: see section title
fOrig = 10e3; 
fs = 50e3; 
%256 point DFT with bin frequencies as given, 
N = 256; 
k = [0:N-1]; 
wbin = 2*pi*k/N; 
kmax = round(fOrig/fs*N);
%need to take the smallest difference between the theoretical and actual
%bin frequency
fOff = min(abs(wbin - 2*pi*fOrig/fs));
%evaluating straddle loss at the dirichlet sinc
dirichlet = diric(fOff, 250);  
%expressing this as a decibal loss from the peak value of the dirichlet
%sinc
Loss = 20*log10(diric(0,250)) -  20*log10(dirichlet); 

%part d: employ hamming window and compute straddle loss
%returns 250 point symmetric hamming window
H = hamming(250); 
Hw = freqz(H, 1, [0 fOff]);
%compare the magnitude of the DTFT of the hamming window for some w' to the
%peak of the mainlobe at DC
HammingOneTwoTwoEl = 20*log10(abs(Hw(1))) - 20*log10(abs(Hw(2)));  

%part e: 250 samples replaced with 512 DFT recalculate both straddle losses
N2 = 512;
k2 = [0:N2-1]; 
wbin2 = 2*pi*k2/N2;
%same process and reasoning as above
fOff2 = min(abs(wbin2 - 2*pi*fOrig/fs)); 
dirichlet = diric(fOff2, 250); 
Loss2 = 20*log10(diric(0,250)) -  20*log10(dirichlet); %loss in dB

Hw2 = freqz(H, 1, [0 fOff2]);
%same loss even when increasing the size of the DFT
HammingOneTwoTwoEl2 = 20*log10(abs(Hw2(1))) - 20*log10(abs(Hw2(2))); 
%% 6: Create a sinewave corrupted by noise as follows
f = 20e6; 
A = 2; 
Fs = 100e6;
T = 1/Fs;
L = 1000;
t = (0:L-1)*T;
n = 1024;
f1 = (-n/2:n/2-1)*(Fs/n);
fx = (0:n-1)*(Fs/n);
%The sinewave is real, at 20MHz, with amplitude 2, reference phase 0 degrees 
%sampled at a rate 100MHz. It is corrupted by AWGN with variance 0.2.
y = A*sin(2*pi*f*t) + sqrt(0.2)*randn(size(t)); 
dB = 30;

%Chebyshev window of length 1000 with 30dB peak sidelobe level
window = chebwin(L,dB);
yChebWind = y.*window';
%block of 100 samples colected and 1024 DFT computed
yDFT = fft(yChebWind, n);
fshift = 20*log10(fftshift(yDFT));
figure('Name','Magnitude Spectrum','NumberTitle','off');
%Graph magnitude spectrum (in dB) with frequency from -fs/2 to fs/2 on a Hz
%scale determined from the DFT coefficients
plot(f1,fshift);
grid on;
title('Magnitude Spectrum Via DFT');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
xlim([(-Fs/2) (Fs/2)]);

%% 7: x[n] is the original unwindowed datablock from problem 6, which in our case is y
%w0[n] is a subsampling function with support from 0 to 999 and within that
%range it is 1 for every 4th sample with an offset of 3 and 0 otherwise,
%i.e. w0[4m+3]=1
y = A*sin(2*pi*f*t) + sqrt(0.2)*randn(size(t)); 
s = 1000;
W = zeros(1000,1)';

for z = 0:249
        W(4*z+3) = 1;
end

xHat = W.*y;
xHatChebWind = xHat.*window';
%block of 100 samples colected and 1024 DFT computed
xHatDFT = fft(xHatChebWind, n);
fshift2 = 20*log10(abs(fftshift(xHatDFT)));
figure('Name','Magnitude Spectrum','NumberTitle','off');
%Graph magnitude spectrum (in dB) with frequency from -fs/2 to fs/2 on a Hz
%scale determined from the DFT coefficients
plot(f1,fshift2);
grid on;
title('Magnitude Spectrum Via DFT');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
xlim([(-Fs/2) (Fs/2)]);
%I notice there are more peaks and the magnitude decreased about 12dB
%The peaks for this new DFT occur at -4.502e7, -2.998e7,
%-2.002e7, -4.98e6, 4.98e6, 2.002e7, 2.998e7, 4.502e7
%The old peaks were -2.002e7 and 2.002e7, so the pattern (?) is the new
%peaks are the old peaks convolved with the DFT of the window function,
%also you know the frequencies because the sampling frequency is 100e6Hz
%and by implementing this window function we are downsampling by 4, so we
%expect the new sampling frequency at 25kHz so aliasing would occur at
%+\-20e6 + n25e6, and if you were to write all of these out, they
%approximately match the values we got

%Compute and graph a stem plot of the magnitude of the 1000 point DFT of
%w0[n] on a linear scale, also plotted is a stem plot of w0[n] 

w0DFT = abs(fftshift(fft(W,1000)))';
%figure;
%stem(W)
%stem plot above is the stem plot of w0[n], uncomment for visual aid and
%thought triggering
figure('Name','Magnitude of 1000 point DFT of w0[n]','NumberTitle','off');
stem(w0DFT)  
grid on;
title('Magnitude Spectrum Via 1000 Point DFT');
ylabel('Magnitude');
xlabel('DFT index k');
%The DFT takes impulse trains to impulse trains, in this case the DFT maps
%the impulse train of w0[n] to the 1000 point DFT of the subsampling
%function, which is a train of deltas; The subsampling function is a train
%at 25MHz and when you're multiplying in time you're convolving with those
%frequencies so it makes sense you're only getting 0,250,500,750 as those
%are all less than 999














