%% Zadanie 1
IM1 = imread('http://pages.mini.pw.edu.pl/~blaszczykl/dydaktyka/ASiSP/lab5/zagadka1.jpg');
IM2 = imread('http://pages.mini.pw.edu.pl/~blaszczykl/dydaktyka/ASiSP/lab5/zagadka2.jpg');

FT1 = fftshift(fft2(IM1));
FT2 = fftshift(fft2(IM2));

n1 = size(IM1,2);   n2 = size(IM1,1);
x = (-floor(n1/2):1:ceil(n1/2)-1)/(n1/2);
y = (-floor(n2/2):1:ceil(n2/2)-1)/(n2/2);

% Widma amplitudowe.
subplot(1,2,1)
    imshow(log(abs(FT1) + 1), [],'XData',[min(x(:)) max(x(:))],...
        'YData',[min(y(:)) max(y(:))])
subplot(1,2,2)
    imshow(log(abs(FT2) + 1), [],'XData',[min(x(:)) max(x(:))],...
        'YData',[min(y(:)) max(y(:))])

%% Widma fazowe
subplot(1,2,1)
    imshow(angle(FT1), [],'XData',[min(x(:)) max(x(:))],...
        'YData',[min(y(:)) max(y(:))])
subplot(1,2,2)
    imshow(angle(FT2), [],'XData',[min(x(:)) max(x(:))],...
        'YData',[min(y(:)) max(y(:))])
    
% Widzimy, dla obrazu zagadka1.jpg czêstotliwoœæ jednego piku wynosi f1_{x,0} = 0.05, f1_{y,0} = -0.1351,
% a dla zagadka2.jph f2_{x,0} = 0.125, f2_{x,y} = -0.03378.

    
%% Zadanie 2
    
% Pamiêtamy, ¿e dla obrazu zagadka1.jpg czêstotliwoœæ jednego piku wynosi f1_{x,0} = 0.05, f1_{y,0} = -0.1351,
% a dla zagadka2.jph f2_{x,0} = 0.125, f2_{x,y} = -0.03378.

% Sprawdzimy dzia³anie naszego filtru na obrazku zagadka1.jpg z zadania 1
% dla sigma = 0.1, n = 2.
sigma = 0.1;
n = 2;
IM = imread('http://pages.mini.pw.edu.pl/~blaszczykl/dydaktyka/ASiSP/lab5/zagadka1.jpg');

FT = fftshift(fft2(IM));

n1 = size(IM,2);
x = (-floor(n1/2):1:ceil(n1/2)-1)/(n1/2);
y = (-floor(n2/2):1:ceil(n2/2)-1)/(n2/2);
[X,Y] = meshgrid(x,y);
H = 1 - 1 ./ (1 + sqrt((X - 0.05).^2 + (Y + 0.1351).^2)./sigma).^(2 * n) - 1 ./ (1 + sqrt((X + 0.05).^2 + (Y - 0.1351).^2)./sigma).^(2 * n);
FT_f = FT.*H;
IM_f = ifft2(ifftshift(FT_f));
imshow(abs(IM_f), [])
% Widzimy wiêc, ¿e filtr dzia³a poprawnie.
%% Sprawdzimy dzia³anie naszego filtru na obrazku zagadka2.jpg z zadania 1 dla sigma = 0.935, n = 2.55.
sigma = 0.935;
n = 2.55;
IM = imread('http://pages.mini.pw.edu.pl/~blaszczykl/dydaktyka/ASiSP/lab5/zagadka2.jpg');

FT = fftshift(fft2(IM));

n1 = size(IM,2);
x = (-floor(n1/2):1:ceil(n1/2)-1)/(n1/2);
y = (-floor(n2/2):1:ceil(n2/2)-1)/(n2/2);
[X,Y] = meshgrid(x,y);
H = 1 - 1 ./ (1 + sqrt((X - 0.125).^2 + (Y + 0.03378).^2)./sigma).^(2 * n) - 1 ./ (1 + sqrt((X + 0.125).^2 + (Y - 0.03378).^2)./sigma).^(2 * n);
FT_f = FT.*H;
IM_f = ifft2(ifftshift(FT_f));
imshow(abs(IM_f), [])
% Widzimy wiêc, ¿e taki filtr dzia³a prawid³owo.

%% Zadanie 3

% Odszyfrowujemy obraz zagadka1.jpg.
% Bierzemy sigma = 0.1, n = 2.
sigma = 0.1;
n = 2;
IM = imread('http://pages.mini.pw.edu.pl/~blaszczykl/dydaktyka/ASiSP/lab5/zagadka1.jpg');

FT = fftshift(fft2(IM));

n1 = size(IM,2);
x = (-floor(n1/2):1:ceil(n1/2)-1)/(n1/2);
y = (-floor(n2/2):1:ceil(n2/2)-1)/(n2/2);
[X,Y] = meshgrid(x,y);
H = 1 - 1 ./ (1 + sqrt((X - 0.05).^2 + (Y + 0.1351).^2)./sigma).^(2 * n) - 1 ./ (1 + sqrt((X + 0.05).^2 + (Y - 0.1351).^2)./sigma).^(2 * n);
FT_f = FT.*H;
IM_f = ifft2(ifftshift(FT_f));
imshow(abs(IM_f), [])
% Zatem obraz zagadka1.jpg to obraz B.

%% Odszyfrowujemy obraz zagadka2.jpg.
% Bierzemy sigma = 0.935, n = 2.55.
sigma = 0.935;
n = 2.55;
IM = imread('http://pages.mini.pw.edu.pl/~blaszczykl/dydaktyka/ASiSP/lab5/zagadka2.jpg');

FT = fftshift(fft2(IM));

n1 = size(IM,2);
x = (-floor(n1/2):1:ceil(n1/2)-1)/(n1/2);
y = (-floor(n2/2):1:ceil(n2/2)-1)/(n2/2);
[X,Y] = meshgrid(x,y);
H = 1 - 1 ./ (1 + sqrt((X - 0.125).^2 + (Y + 0.03378).^2)./sigma).^(2 * n) - 1 ./ (1 + sqrt((X + 0.125).^2 + (Y - 0.03378).^2)./sigma).^(2 * n);
FT_f = FT.*H;
IM_f = ifft2(ifftshift(FT_f));
imshow(abs(IM_f), [])
% St¹d obraz zagadka2.jpg to obraz A.