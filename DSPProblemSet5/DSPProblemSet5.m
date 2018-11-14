% Michael Lendino ECE310 DSP PSET 5 MATLAB 
clc;
clear all;
%% Use freq2z to compute H(w), find the maximum absolute value of imaginary part and also the maximum numerical difference between the real part returned by freqz2, zero out imaginary part after that

h = 1/6*[1 4 1; 4 -20 4; 1 4 1];
[H, f1, f2] = freqz2(h, [100 100]);

Herror = max(max(abs(imag(H))));
%I think what I did is correct, so the error is 0

%Contour Plot
figure('Name','Contour Plot for H(w)','NumberTitle','off');
contour(pi*f1, pi*f2, H)
grid on;
xlabel('Horizontal Digital Radian Frequency');
ylabel('Vertical Digital Radian Frequency');
title('Contour Plot for H(w)');
%Surface Plot
figure('Name','Surface Plot for H(w)','NumberTitle','off');
surface(pi*f1, pi*f2, H)
xlabel('Horizontal Digital Radian Frequency (wx)');
ylabel('Vertical Digital Radian Frequency (wy)');
zlabel('Magnitude Response for H(w)');
grid on;
title('Surface Plot for H(w)');
%I'm pretty sure this is a lowpass filter 

%Examining the slice H(wx,0) for wx between -pi and pi and superimposing
%the ideal curve if this were to exactly match the frequency response of
%the laplacian
Hslice = H(:,50);
figure('Name','Slice and Ideal Curve','NumberTitle','off');
plot(pi*f1, Hslice)
hold on
plot(pi*f1,-abs(pi*f1).^2)
legend('The Slice H(wx,0)', 'Frequency Response of Laplacian');
xlabel('Digital Radian Frequency (wx)');
ylabel('Magnitude Response');
grid on;
xlim([-pi pi])
hold off
%Superimpose plots of H and the spectrum of the Laplacian with axes f1 and
%f2
[ff1,ff2] = meshgrid(f1, f2);
laplacianSpectrum = -(ff1.^2 + ff2.^2) * pi;
figure('Name','Plots of H and the Spectrum of the Laplacian','NumberTitle','off');
surf(f1, f2, laplacianSpectrum)
hold on
surf(f1, f2, H)
xlabel('Normalized Digital Radian Frequency (f1)');
ylabel('Normalized Digital Radian Frequency (f2)');
zlabel('Magnitude Response');
grid on;
hold off

%% Upsampling by 2: g(2n) = h(n) and at all points n=(n1,n2) where both coordinates are not even, g(n)=0, we insert zeros in every other row and every other column of h

m = floor(rand*30)+ 1;
n = floor(rand*30)+ 1;

if mod(m,2) == 0
    m = m + 1;
end

if mod(n,2) == 0
    n = n + 1;
end

M = floor(rand(m,n)*30) + 1;
N = interp1(M);

%Applying interpolator function to specific h from previous problem to
%obtain g
g = interp1(h);
G = freqz2(g, [100 100]);

figure('Name','Contour Plot for G(w1,w2)','NumberTitle','off');
contour(pi*f1, pi*f2, G)
xlabel('Horizontal Digital Radian Frequency');
ylabel('Veritcal Digital Radian Frequency');
grid on;
title('Contour Plot for G(w1,w2)');

figure('Name','Surface plot for G(w1,w2)','NumberTitle','off');
surface(pi*f1, pi*f2, G)
xlabel('Horizontal Digital Radian Frequency');
ylabel('Vertical Digital Radian Frequency');
zlabel('Magnitude Response for G(w1,w2)');
grid on;
title('Surface Plot for G(w1,w2)');
%This is not a plausible approximation to the Laplacian operator
%This could be fixed by lowpass filtering and then scaling the output.

%% Sobel Edge Detection
%following matrix is an approximation of the partial derivative with
%respect to y
hy = 1/8*[-1 -2 -1; 0 0 0; 1 2 1];
hx = hy'; 
%obtain 2D graphical representations for the magnitude responses of hx and
%hy and confirm that these correspond to horizontal and vertical highpass
%filters
Hy = freqz2(hy, [100 100]);
Hx = freqz2(hx, [100 100]);

figure('Name','Contour Plot for Hy(w)','NumberTitle','off');
contour(pi*f1, pi*f2, abs(Hy))
xlabel('Horizontal Digital Radian Frequency');
ylabel('Vertical Digital Radian Frequency');
grid on;
title('Contour Plot for Hy(w)');
%Based on my understanding, this is the vertical highpass filter

figure('Name','Contour Plot for Hx(w)','NumberTitle','off');
contour(pi*f1, pi*f2, abs(Hx))
xlabel('Digital Radian Frequency (wx)');
ylabel('Digital Radian Frequency (wy)');
grid on;
title('Contour Plot for Hx(w)')
%Based on my understanding, this is the horizontal highpass filter

%Import the image circuit.tif using imread, convert it to dobule and apply
%sobel edge detection process, use image to display the original image and
%imtool to display the 0-1 edge identification matrix
%Thank you Armaan for helping me with this part of the question! 
ckt = double(imread('circuit.tif'));
figure('Name','Original Image','NumberTitle','off');
image(ckt)
M1 = Sobel(ckt, 10, 2);
%Compared with median value of the image
med = median(ckt(:));
M2 = Sobel(ckt, med, 2);

%Sobel edge detection modified with simplified gradient approximation
M3 = Sobel(ckt, 10, 1);


