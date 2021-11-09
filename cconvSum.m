% FunkcjacconvSum(x,y)
function splot = cconvSum(x,y)
if (length(x) ~= length(y))
    return;
end
N = length(x);
splot = zeros(1, N);
for n = 1:N
    for k = 1:N
    if(n + 1 - k < 1)      
        splot(n) = splot(n) + x(k) * y(n + 1 - k + N);
    else
        splot(n) = splot(n) + x(k) * y(n + 1 - k);
    end
    end
end