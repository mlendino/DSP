function N = Sobel(g, T, L)
hy = 1/8*[-1 -2 -1; 0 0 0; 1 2 1];
hx = hy'; 

delx = conv2(hx,g);
dely = conv2(hy,g);
%removing appropriate rows and columns after convolution to obtain an
%output image the same size as the original and then computing the contrast
%and comparing to the threshold
delx = delx(2:end-1, 2:end-1);
dely = dely(2:end-1, 2:end-1);

%Gradient vector as defined for L-2 Norm
if L == 2
deltag = sqrt(delx.^2 + dely.^2);
end

%Simplified approximation to the gradient for L-1 Norm
if L ==1
deltag = abs(delx) + abs(dely);
end
N = deltag > T;

imtool(N)
end