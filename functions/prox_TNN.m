function [X] = prox_TNN(Y,rho)

[n1,n2,n3] = size(Y);
n12 = min(n1,n2);
Y = fft(Y,[],3);

U = zeros(n1,n12,n3);
V = zeros(n2,n12,n3);
S = zeros(n12,n12,n3);

for i = 1 : n3
    [U(:,:,i),s,V(:,:,i)] = svd(Y(:,:,i),'econ');
    s = diag(s);
    s = max(s-rho,0);    
    S(:,:,i) = diag(s);

    Y(:,:,i) = U(:,:,i)*S(:,:,i)*V(:,:,i)';
end

X = ifft(Y,[],3);
end


