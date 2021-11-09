% Funkcja cconvDFT(x,y)
function splot = cconvDFT(x,y)

if(length(x) ~= length(y))
    return;
end

%DF - dyskretne przekszta³cenie Fouriera
%Z twierdzenia o DFT splotu mamy: DF(splot x i y) = N * DF(x) * DF(y),
%gdzie N jest d³ugoœci¹ wektora x i y.

N = length(x);
dfsplotu = zeros(N,1);
dfx = fft(x,N);
dfy = fft(y,N);
for n = 1:N
    dfsplotu(n) = 4 * dfx(n) * dfy(n); 
end
splot = ifft(1/4 * dfsplotu);