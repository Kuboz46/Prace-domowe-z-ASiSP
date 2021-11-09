%% Zadanie 1

N = 512; Fs = 1000; T = N/Fs;
t = 0:1/Fs:(N-1)/Fs;
f0 = 62.5;  f1 = 125;   f2 = 300;

% Sygna³ na wejœciu systemu.
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

% Sygna³ na wyjœciu systemu.
y = ifft(ifftshift(Y))
%% Zadanie 2

disp('Sprawdzamy poprawnoœæ dzia³ania skryptu z zadania 1.');
[audio, Fs] = audioread('dzwiek.wav');
% Odtwarzanie dŸwiêku pierwotnego.
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
% Odtwarzanie dŸwiêku otrzymanego po filtracji.
soundsc(y, Fs)
% S³yszymy, ¿e otrzymany dŸwiêk brzmi tak, jakby by³ odtwarzany ze starego
% radia.
%% Zadanie 3

N = 1536; Fs = 3000; T = N/Fs;
t = 0:1/Fs:(N-1)/Fs;
f0 = 62.5;  f1 = 125;   f2 = 300;

% Sygna³ na wejœciu systemu.
x = sin(2*pi*f0 * t) + 0.2 * cos(2* pi * f1 * t) - 0.5 * sin(2*pi*f2*t) + 1.2;

% Widmo sygna³u na wejœciu systemu.
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

% Widmo sygna³u na wyjœciu systemu.
Y = X .* H;

% Porównujemy widmo amplitudowe sygna³u na wejœciu i wyjœciu systemu.
figure(1);
plot(f, abs(X), 'b.', f, abs(Y), 'g.', 'MarkerSize', 4)
xlabel('f')
title('Widmo amplitudowe sygna³u na wejœciu i wyjœciu systemu')
legend('|X(f)|', '|Y(f)|')
axis([-1700 1700 0 3.5])
% Na wykresie widzimy, ¿e amplituda widma sygna³u sfiltrowanego jest
% mniejsza ni¿ widma sygna³u na wejœciu systemu.
% Dla |f| <= 1500 |X(f)| jest niezerowe. Dla pozosta³ych f
% |X(f)| = 0.
% |X(f)| jest kawa³kami ci¹g³e w przedzia³ach [-1500, -300), [-300, 300],
% [300, 1500]. Punkty f = -300, f = 300 s¹ punktami nieci¹g³oœci |X(f)|.
% Dla |f| <= 1000 |Y(f)| jest niezerowe. 
% Dla |f| > 1000 |Y(f)| jest równe 0, co wynika wprost z definicji filtru
% H, gdy¿ dla |f| > 1000 H(f) = 0.
% |Y(f)| jest kawa³kami ci¹g³e w przedzia³ach [-1500, -1000), [-1000, -300), [-300, 300],
% (300, 1000], (1000, 1500].
% Cztery punkty: f = -1000, f = -300 f = 300, f = 1000 s¹ punktami nieci¹g³oœci |Y(f)|.
% W przedziale [-300, 300] widmo amplitudowe sygna³u y jest bardzo ma³e w porównaniu z widmem amplitudowym
% na sumie przedzia³ów {-1000, -300), (300, 1000].


%% Teraz porównujemy widmo fazowe sygna³u na wejœciu i wyjœciu systemu.
figure(1);
plot(f, angle(X), 'b.', f, angle(Y), 'g.', 'MarkerSize', 4)
xlabel('f')
title('Widmo fazowe sygna³u na wejœciu i wyjœciu systemu')
legend('arg(X(f))', 'arg(Y(f))')
axis([-1700 1700 -10 10])
% Widzimy na wykresie, ¿e:
% Dla |f| <= 1500 arg(X(f)) jest niezerowe. Dla pozosta³ych f
% arg(X(f)) = 0.
% Dla |f| <= 1000 arg(Y(f)) jest niezerowe. 
% Dla |f| > 1000 arg(Y(f)) jest równe 0, co wynika wprost z definicji filtru
% H, gdy¿ dla |f| > 1000 H(f) = 0.
% arg(X(f)) jest kawa³kami liniowe oraz ci¹g³e w przedzia³ach [-1500, -300), [-300,0],
% [0, 300], (300, 1500]. Punkty f = -1500, f = -300, f = 0, f = 300, f = 1500 s¹
% punktami nieci¹g³oœci arg(X(f)).
% arg(Y(f)) jest te¿ kawa³kami ci¹g³e w przedzia³ach [-1500, -1000), [-1000, -300), [-300,
% 0], [0, 300], (300, 1000], (1000, 1500]. Punkty f = -1000, f = -300, f = 0, f = 300, f = 1000 s¹
% punktami nieci¹g³oœci arg(X(f)).
% arg(Y(f)) jest kawa³kami liniowe, ale tylko w przedzia³ach [-1000, -300),
% (300, 1000]. Dla pozosta³ych f arg(X(f)) jest nieliniowe.

% Dla |f| <= 1000 |H(f)|, arg(H(f)) s¹ ró¿ne od 0.
% Dla pozosta³ych f |H(f)| = 0, arg(H(f)) = 0.
% Zatem filtr H jest dolnoprzepustowy.