function M = interp1(A)
m = length(A(:,1));
n = length(A(1,:));
%produces matrix so that outermost rows and columns is not an appended row
%or column of all 0's
M = zeros(2*m - 1, 2*n - 1);

for i = 1:(2*m -1) 
    for j = 1:(2*n -1)
        if mod(i,2) ~=0 && mod(j,2)~=0
            M(i,j) = A(ceil(i/2), ceil(j/2));
        end
    end
end
end
