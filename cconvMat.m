% Funkcja cconvMat(x,y)
function splot = cconvMat(x,y)
if (length(x) ~= length(y))
    return;
end
N = length(x);
Y = zeros(N,1);
for k = 1:N
    Y(k) = y(mod(1 - k, N) + 1);
end
X = zeros(N);
for k = 1:N
    for l = 1:N
        X(k,l) = x(mod(k + l - 2,N) + 1);
    end
end
splot = X * Y;