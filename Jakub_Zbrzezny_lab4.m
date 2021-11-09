%% Zadanie 1

N = 512; Fs = 1000; T = N/Fs;
t = 0:1/Fs:(N-1)/Fs;
f0 = 62.5;  f1 = 125;   f2 = 300;

% Sygna� na wej�ciu systemu.
x = sin(2*pi*f0 * t) + 0.2 * cos(2* pi * f1 * t) - 0.5 * sin(2*pi*f2*t) + 1.2;

f = -(N/2)/T:1/T:(N/2-1)/T;
X = fftshift(fft(x));

k1 = 0.01;
k2 = 0.35;
tau1 = 1e-2;
tau2 = 1e-4;
H = (abs(f) <= 300) .* (k1 ./(2 * pi * 1i * tau1 * f + 1)) + (300 < abs(f)) .* (abs(f) <= 1000) .* ...
   (k2./(2 * pi * 1i * tau2 * f + 1)); 
% Filtr H jest symetryczny.

Y = X .* H;

% Sygna� na wyj�ciu systemu.
y = ifft(ifftshift(Y))
%% Zadanie 2

disp('Sprawdzamy poprawno�� dzia�ania skryptu z zadania 1.');
[audio, Fs] = audioread('dzwiek.wav');
% Odtwarzanie d�wi�ku pierwotnego.
soundsc(audio, Fs)

N = length(audio); T = N/Fs;
f = (-floor(N/2):1:(ceil(N/2) - 1))/T;
X = fftshift(fft(audio));
k1 = 0.01;
k2 = 0.35;
tau1 = 1e-2;
tau2 = 1e-4;
H = (abs(f) <= 300) .* (k1 ./(2 * pi * 1i * tau1 * f + 1)) + (300 < abs(f)) .* (abs(f) <= 1000) .* ...
   (k2./(2 * pi * 1i * tau2 * f + 1));
% Filtr H jest symetryczny.

Y = X .* H.';
y = ifft(ifftshift(Y));
pause(8);
% Odtwarzanie d�wi�ku otrzymanego po filtracji.
soundsc(y, Fs)
% S�yszymy, �e otrzymany d�wi�k brzmi tak, jakby by� odtwarzany ze starego
% radia.
%% Zadanie 3

N = 1536; Fs = 3000; T = N/Fs;
t = 0:1/Fs:(N-1)/Fs;
f0 = 62.5;  f1 = 125;   f2 = 300;

% Sygna� na wej�ciu systemu.
x = sin(2*pi*f0 * t) + 0.2 * cos(2* pi * f1 * t) - 0.5 * sin(2*pi*f2*t) + 1.2;

% Widmo sygna�u na wej�ciu systemu.
X = fftshift(fft(x));

% Definiujemy filtr H.
k1 = 0.01;
k2 = 0.35;
tau1 = 1e-2;
tau2 = 1e-4;
f = -(N/2)/T:1/T:(N/2-1)/T;
H = (abs(f) <= 300) .* (k1 ./(2 * pi * 1i * tau1 * f + 1)) + (300 < abs(f)) .* (abs(f) <= 1000) .* ...
   (k2./(2 * pi * 1i * tau2 * f + 1)); 
% Filtr H jest symetryczny.

% Widmo sygna�u na wyj�ciu systemu.
Y = X .* H;

% Por�wnujemy widmo amplitudowe sygna�u na wej�ciu i wyj�ciu systemu.
figure(1);
plot(f, abs(X), 'b.', f, abs(Y), 'g.', 'MarkerSize', 4)
xlabel('f')
title('Widmo amplitudowe sygna�u na wej�ciu i wyj�ciu systemu')
legend('|X(f)|', '|Y(f)|')
axis([-1700 1700 0 3.5])
% Na wykresie widzimy, �e amplituda widma sygna�u sfiltrowanego jest
% mniejsza ni� widma sygna�u na wej�ciu systemu.
% Dla |f| <= 1500 |X(f)| jest niezerowe. Dla pozosta�ych f
% |X(f)| = 0.
% |X(f)| jest kawa�kami ci�g�e w przedzia�ach [-1500, -300), [-300, 300],
% [300, 1500]. Punkty f = -300, f = 300 s� punktami nieci�g�o�ci |X(f)|.
% Dla |f| <= 1000 |Y(f)| jest niezerowe. 
% Dla |f| > 1000 |Y(f)| jest r�wne 0, co wynika wprost z definicji filtru
% H, gdy� dla |f| > 1000 H(f) = 0.
% |Y(f)| jest kawa�kami ci�g�e w przedzia�ach [-1500, -1000), [-1000, -300), [-300, 300],
% (300, 1000], (1000, 1500].
% Cztery punkty: f = -1000, f = -300 f = 300, f = 1000 s� punktami nieci�g�o�ci |Y(f)|.
% W przedziale [-300, 300] widmo amplitudowe sygna�u y jest bardzo ma�e w por�wnaniu z widmem amplitudowym
% na sumie przedzia��w {-1000, -300), (300, 1000].


%% Teraz por�wnujemy widmo fazowe sygna�u na wej�ciu i wyj�ciu systemu.
figure(1);
plot(f, angle(X), 'b.', f, angle(Y), 'g.', 'MarkerSize', 4)
xlabel('f')
title('Widmo fazowe sygna�u na wej�ciu i wyj�ciu systemu')
legend('arg(X(f))', 'arg(Y(f))')
axis([-1700 1700 -10 10])
% Widzimy na wykresie, �e:
% Dla |f| <= 1500 arg(X(f)) jest niezerowe. Dla pozosta�ych f
% arg(X(f)) = 0.
% Dla |f| <= 1000 arg(Y(f)) jest niezerowe. 
% Dla |f| > 1000 arg(Y(f)) jest r�wne 0, co wynika wprost z definicji filtru
% H, gdy� dla |f| > 1000 H(f) = 0.
% arg(X(f)) jest kawa�kami liniowe oraz ci�g�e w przedzia�ach [-1500, -300), [-300,0],
% [0, 300], (300, 1500]. Punkty f = -1500, f = -300, f = 0, f = 300, f = 1500 s�
% punktami nieci�g�o�ci arg(X(f)).
% arg(Y(f)) jest te� kawa�kami ci�g�e w przedzia�ach [-1500, -1000), [-1000, -300), [-300,
% 0], [0, 300], (300, 1000], (1000, 1500]. Punkty f = -1000, f = -300, f = 0, f = 300, f = 1000 s�
% punktami nieci�g�o�ci arg(X(f)).
% arg(Y(f)) jest kawa�kami liniowe, ale tylko w przedzia�ach [-1000, -300),
% (300, 1000]. Dla pozosta�ych f arg(X(f)) jest nieliniowe.

% Dla |f| <= 1000 |H(f)|, arg(H(f)) s� r�ne od 0.
% Dla pozosta�ych f |H(f)| = 0, arg(H(f)) = 0.
% Zatem filtr H jest dolnoprzepustowy.