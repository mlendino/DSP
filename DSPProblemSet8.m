%Michael Lendino DSP PSET 8
clc;
clear all;
%% Design an 8th order bandpass digital elliptic filter with 1.5 dB ripple
%%in the passband, 30dB ripple in the stopband, and passband from 0.3 to
%%0.6 on a scale where 1 is the Nyquist bandwidth

Apass = 1.5;
Astop = 30;
wpass = [0.3, 0.6];
n = 8;

[z,p,k] = ellip(n/2,Apass,Astop,wpass,'bandpass');

%Convert to tf form and obtain a graph of the magnitude response
[b,a] = zp2tf(z,p,k);
[h,w] = freqz(b,a,100000);
ellipdB = 20*log10(abs(h));
figure('Name','Magnitude Response for 8th Order Bandpass Digital Elliptic Filter','NumberTitle','off');
plot(w/pi,ellipdB);
grid on;
title('Digital: Elliptic Magnitude Response');
ylabel('Magnitude (dB)');
xlabel('Normalized Digital Radian Frequency');

%Convert the zpk form to sos form with up ordering and again with down
%ordering with L infinity scaling

[sosu, gu] = zp2sos(z,p,k,'up',inf);
[sosd, gd] = zp2sos(z,p,k,'down',inf);

[z1,p1,k1] = sos2zp(sosu,gu);
[z2,p2,k2] = sos2zp(sosd,gd);
%P1 is increasing, and P2 is decreasing, as desired
P1 = abs(p1);
P2 = abs(p2);
%look at sosu matrix to figure out which entries correspond to the
%coefficients, 4 rows so going to have 4 transfer functions for up and down
A1u = freqz(1, sosu(1, 4:6)', 100000);
A2u = freqz(1, sosu(2, 4:6)', 100000);
A3u = freqz(1, sosu(3, 4:6)', 100000);
A4u = freqz(1, sosu(4, 4:6)', 100000);

B1u = freqz(sosu(1, 1:3)', 1, 100000);
B2u = freqz(sosu(2, 1:3)', 1, 100000);
B3u = freqz(sosu(3, 1:3)', 1, 100000);
B4u = freqz(sosu(4, 1:3)', 1, 100000);

A1d = freqz(1, sosd(1, 4:6)', 100000);
A2d = freqz(1, sosd(2, 4:6)', 100000);
A3d = freqz(1, sosd(3, 4:6)', 100000);
A4d = freqz(1, sosd(4, 4:6)', 100000);

B1d = freqz(sosd(1, 1:3)', 1, 100000);
B2d = freqz(sosd(2, 1:3)', 1, 100000);
B3d = freqz(sosd(3, 1:3)', 1, 100000);
B4d = freqz(sosd(4, 1:3)', 1, 100000);

magtf1u = 20*log10(abs(gu*A1u));
magtf2u = 20*log10(abs(gu*A1u.*A2u.*B1u));
magtf3u = 20*log10(abs(gu*A1u.*A2u.*A3u.*B1u.*B2u));
magtf4u = 20*log10(abs(gu*A1u.*A2u.*A3u.*A4u.*B1u.*B2u.*B3u));

magtf1d = 20*log10(abs(gu*A1d));
magtf2d = 20*log10(abs(gu*A1d.*A2d.*B1d));
magtf3d = 20*log10(abs(gu*A1d.*A2d.*A3d.*B1d.*B2d));
magtf4d = 20*log10(abs(gu*A1d.*A2d.*A3d.*A4d.*B1d.*B2d.*B3d));

figure('Name','Magnitude Response for Successive Cumulative Transfer Functions (UP)','NumberTitle','off');
hold on;
plot(w/pi,magtf1u);
plot(w/pi,magtf2u);
plot(w/pi,magtf3u);
plot(w/pi,magtf4u);
grid on;
title('Magnitude Response for Successive Cumulative Transfer Functions (UP)');
legend('Delay Chain 1', 'Delay Chain 2', 'Delay Chain 3', 'Delay Chain 4');
ylabel('Magnitude (dB)');
xlabel('Normalized Digital Radian Frequency');
hold off;

figure('Name','Magnitude Response for Successive Cumulative Transfer Functions (DOWN)','NumberTitle','off');
hold on;
plot(w/pi,magtf1d);
plot(w/pi,magtf2d);
plot(w/pi,magtf3d);
plot(w/pi,magtf4d);
grid on;
title('Magnitude Response for Successive Cumulative Transfer Functions (DOWN)');
legend('Delay Chain 1', 'Delay Chain 2', 'Delay Chain 3', 'Delay Chain 4');
ylabel('Magnitude (dB)');
xlabel('Normalized Digital Radian Frequency');
hold off;

B1uz = sosu(1, 1:3);
B4dz = flipud(sosd(4, 1:3));
check = B1uz/B4dz;
%They match up to the scaling factor of 0.2927, so it checks out

%% Sensitivity Properties of Parallel Allpass Realizations
%Use fi to compute fixed point representation of the coefficients for
%b=5bits. Use these values as standard double precision to examine the
%resulting quantized transfer functions for the quantized transfer
%functions for the original form and for the parallel allpass form
numcoefforig = [0.1336 0.0568 0.0563 0.1336];
denomcoefforig = [1 -1.5055 1.2630 -0.3778];
[Horig, wo] = freqz(numcoefforig, denomcoefforig, 1e4);

ydo = fi(denomcoefforig, 1, 5, 3);
ado = ydo.data;

yno = fi(numcoefforig, 1, 5, 6);
bno = yno.data;

%Check that the ideal gain is 1 at w=0 and 0 at w=pi, then compute the
%actual gains for the filters with quantized coefficients and report the
%error in dB
Hquant = freqz(bno, ado, 1e4); %transfer function with quantized coefficients
err0 = 20*log10(abs(abs(Horig(1)) - abs(Hquant(1)))); % the error is -21.7499 dB when w = 0
errpi = 20*log10(abs(abs(Horig(end)) - abs(Hquant(end)))); %the error is -80.3538 dB when w = pi

%Compute the frequency response for all three filters for 1e4 points from 0
%to pi. 
bAllpass1 = [-0.4954 1];
aAllpass1 = [1 -0.4954];

yAa1 = fi(aAllpass1, 1, 5, 3);
aA1 = yAa1.data;

yAb1 = fi(bAllpass1, 1, 5, 3);
bA1 = yAb1.data;
tfA1 = tf(bA1, aA1);

bAllpass2 = [0.7626 -1.0101 1];
aAllpass2 = [1 -1.0101 0.7626];

yAa2 = fi(aAllpass2, 1, 5, 3);
aA2 = yAa2.data;

yAb2 = fi(bAllpass2, 1, 5, 3);
bA2 = yAb2.data;

H02 = freqz(bA2, aA2, 1e4);
tfA2 = tf(bA2, aA2);

HsumA = minreal(0.5*(tfA1 + tfA2));
[bsumA, asumA] = tfdata(HsumA);
HsumA = freqz(bsumA{1,1}, asumA{1,1}, 1e4);

%Compute the maximum of |H(w)-Hq0(w)| and |H(w) - HQA(w)|;
maxdiff1 = max(abs(Horig - Hquant)); 
maxdiff2 = max(abs(Horig - HsumA));

magHorig = 20*log10(abs(Horig));
magHquant = 20*log10(abs(Hquant));
magHsum = 20*log10(abs(HsumA));

%Superimpose plots of the magnitude responses of the three filters on a
%decibal scale 
figure('Name','Magnitude Responses of Each Filter','NumberTitle','off');
plot(wo, magHorig);
hold on;
plot(wo, magHquant);
plot(wo, magHsum);
grid on;
legend ('H', 'HQ0', 'HQA');
xlabel('Normalized Digital Radian Frequency');
ylabel('Magnitude (dB)');
title('Magnitude Responses of Each Filter');
ylim([-40 1]);
hold off;

%Computing the maximum deviation of the filter gains from the gain of the
%infinite precision filter in the passband, thank you aziza for helping me
%figure out how to do this part of this question especially dealing with
%the indices

[gaindiff0, i0] = max(abs(abs(Horig(1:3001)) - abs(Hquant(1:3001))));  
[gaindiff1, i1] = max(abs(abs(Horig(1:3001)) - abs(Hquant(1:3001)))); 
[gaindiff2, i2] = max(abs(abs(Horig(1:3001)) - abs(Hquant(1:3001)))); 

%Compute the deviation from the equiripple characteristics (check for local
%maxima and miinima in the passband and stopband and observe to what extent
%they are not uniform)

%first local min of passband at index ~1825, found by looking at the graph
%at around where the local min was, then going into variable and confirming
%it, similar process for the other indices

[localmin1orig, iminp1orig] = min(abs(Horig(1:3001)));  
[localmin1quant, iminp1quant] = min(abs(Hquant(1:2500)));  
diffmin1p = abs(localmin1orig - localmin1quant);  

%2nd local min at edge of passband 
localmin2orig = abs(Horig(3001));
localmin2quant = abs(Hquant(3001));
diffmin2p = abs(localmin2orig - localmin2quant); 

%1st local max in passband 
localmax1orig = abs(Horig(1));
localmax1quant = abs(Hquant(1));
diffmax2p = abs(localmax1orig - localmax1quant); 

%2nd local max in passband
[localmax2orig, imaxp2c] = max(abs(Horig(1825:3001)));  
[localmax2qnat, imaxp20] = max(abs(Hquant(1825:3001)));  
diffmax22p = abs(localmax2orig - localmax2qnat);

%1st local min of stopband 
[localmin1sorig, imins1c] = min(abs(Horig(3001:5000)));   
[localmin1squant, imins10] = min(abs(Hquant(3001:end)));  
imins10 = imins10 + 3001;
imins1c = imins1c + 3001;
diffmin1s = abs(localmin1sorig - localmin1squant);  

%2nd local min in stopband 
localmin2sorig = abs(Horig(end));
localmin2squant = abs(Hquant(end));
diffmin2s = abs(localmin2sorig - localmin2squant); 


%1st local max in stopband
localmax1sorig = abs(Horig(3002));
localmax1squant = abs(Hquant(3002));
diffmax2s = abs(localmax1sorig - localmax1squant); 

%2nd local max in stopband 
[localmax2sorig, imax2cs] = max(abs(Horig(4100:end)));   
[localmaxsquant, imax20s] = max(abs(Hquant(4100:end)));  
imax20s = imax20s + 4100;
imax2cs = imax2cs + 4100;
diffmax22s = abs(localmax2sorig - localmaxsquant); 

%Superimpose graphs of the group delay in the passband 
grpdelayorig = grpdelay(numcoefforig, denomcoefforig, 3000, [0 0.3*pi]);
grpdelayquant = grpdelay(bno, ado, 3000, [0 .3*pi]);
grpdelaysum = grpdelay(bsumA{1,1}, asumA{1,1}, 3000, [0 0.3*pi]);

figure('Name','Group Delays in the Passband','NumberTitle','off');
plot(wo(1:3000)*pi, grpdelayorig)
hold on
plot(wo(1:3000)*pi, grpdelayquant)
plot(wo(1:3000)*pi, grpdelaysum)
legend ('Group Delay of H', 'Group Delay of HQ0', 'Group Delay of HQA')
xlabel('Normalized Digital Radian Frequency')
ylabel('Seconds')
xlim([0 0.3*pi])
title('Group Delays in the Passband')
hold off
%The parallel all-pass realization is more accurate 

%Compute the poles and zeros of the original filter and each of the reduced
%precision realizations, not the individual allpass factors
[zorig, porig, korig] = tf2zp(numcoefforig, denomcoefforig);
[zQ0, pQ0, kQ0] = tf2zp(bno, ado);
[zQA, pQA, kQA] = tf2zp(bsumA{1,1}, asumA{1,1});

figure('Name','Poles and Zeros for Each Filter','NumberTitle','off');
subplot(3,1,1)
zplane(zorig, porig)

subplot(3,1,2)
zplane(zQ0, pQ0)

subplot(3,1,3)
zplane(zQA, pQA)

%The zero didn't move? Maybe I plotted the wrong things or calculated
%something wrong, but they all look the same to me on the graph.

