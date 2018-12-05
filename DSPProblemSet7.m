%Michael Lendino DSP PSET 7 
clc;
clear all
%% Question 5
syms k;
syms z;
N = 12;

h = [0.15774243 0.69950381 1.06226376 0.44583132 ...
-0.3199866 -0.18351806 0.13788809 0.03892321 ...
-0.04466375 0.000783251 0.006756062 -0.001523534]
%Thank you aziza for helping me decipher the homework and get the coefficients for each FIR filter below 
h1 = zeros(12,1);
f0 = zeros(12,1);
f1 = zeros(12,1);
for k = 0:11
    h1(k+1) = ((-1)^k)*h(12-k);
    f0(k+1) = h(12-k);
    f1(k+1) = ((-1)^k)*h(k+1);
end

[H0w, w] = freqz(h, 1, 1000);
[H1w, w] = freqz(h1, 1, 1000);
[F0w, w] = freqz(f0, 1, 1000);
[F1w, w] = freqz(f1, 1, 1000);

%Compute and graph the magnitude of H0(w) and H1(w) superimposed on the
%same axes both on a linear scale from 0 to pi
figure('Name','Magnitude Responses of H0(w) and H1(w) ','NumberTitle','off');
plot(w,abs(H0w))
xlim([0 pi]);
xlabel('Frequency (linear scale)');
ylabel('Magnitude');
grid on;
hold on;
plot(w,abs(H1w))
hold off;

% figure('Name','Magnitude Responses of F0(w) and F1(w) ','NumberTitle','off');
% plot(w,abs(F0w))
% hold on;
% plot(w,abs(F1w))
% hold off;
% 

%Confirm their magnitude responses are the same, their computation is lines
%24 and 25.
err1 = max(abs(H0w) - abs(F0w));
err2 = max(abs(H1w) - abs(F1w));

%The value of the sum of the squares of the magnitudes of H0(w) and H1(w) should be a constant. The constant is 4.000
constt = (abs(H0w).^2) + (abs(H1w).^2);

%H0 and H1 are a QMF iff pair of H0(w) = H1(w+pi) i.e. H1 ought to be the
%reflection about pi/2. Since this value is pretty much zero, theyre
%confirmed QMF.

QMF = max(abs(H1w) - flipud(abs(H0w)));

%Construct E(z), i.e. the polyphase matrix corresponding to the filter bank, and verify that it is paraunitary
%effectively upsampling by 2, for k=l=0, so 12 turns to 24
e00 = zeros(1,24)
e00 = h(1:2:end);

e01 = zeros(1,24)
e01 = h(2:2:end);

e10 = zeros(1,24)
e10 = h1(1:2:end);

e11 = zeros(1,24)
e11 = h1(2:2:end);

E00z = dsppset(e00) 
E01z = dsppset(e01)
E10z = dsppset(e10)
E11z = dsppset(e11)

E = [E00z E01z ; E10z E11z];
Ehp = [para(E00z) para(E10z) ; para(E01z) para(E11z)];
paraunitary = vpa(expand(simplify(Ehp*E)),2);
%There is no No? but c=2.

%R(z) is the polyphase matrix corresponding to the synthesis bank
r00 = zeros(1,24)
r00 = f0(2:2:end)

r01 = zeros(1,24)
r01 = f0(1:2:end)

r10 = zeros(1,24)
r10 = f1(2:2:end)

r11 = zeros(1,24)
r11 = f1(1:2:end)

R00z = dsppset(r00) 
R01z = dsppset(r01)
R10z = dsppset(r10)
R11z = dsppset(r11)

R = [R00z R01z; R10z R11z];
RE = vpa(expand(simplify(R*E)),2)
%i am not sure what n0 is :(