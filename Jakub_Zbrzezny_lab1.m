%% Jakub Zbrzezny, numer indeksu: 286689

%% Zadanie 1
%% (1)

% Przyjmujemy, ze T = pi, t0 = pi/3.

% Tworzymy widmo amplitudowe i fazowe sygnalu x(t).
% x(t) jest suma po wszystkich n calkowitych z cn * exp(-j*n*omega0*t).
% x(t) jest rowne a0 + suma od 1 do nieskonczonosci z (an * cos(n omega0 t)
% + bn * sin(n omega0 t).
% Z pracy domowej otrzymalismy, ze:
% omega0 = pi/T;
% a0 = 0;
% dla k >= 1 ak = ((4*T) / (pi * k^2 * (T^2 - t0^2))) * (cos(k*pi*t0/T) - 
% (-1)^k);
% dla k >= 1 bk = ((4*T)/(pi * k^2 * (T^2 - t0^2))) * sin(k*pi*t0/T).
% Z wykladu wiemy, ze:
% dla n >= 0 an = 2 Re(cn);
% dla n >= 1 bn = -2 Im(cn).
% Ale w Matlabie nie mozemy wygenerowac widma dla nieskonczenie wielu c,
% wiec wygenerujemy dla n takich, ze -10 <= n <= 10.


%% Widmo amplitudowe
% Widmem amplitudowym sygnalu x(t) jest {|cn|: n jest liczba calkowita}.
clear;
T = pi;
t0 = pi/3;
omega0 = pi/T;
a0 = 0;

for k=-10:1:10
ak = ((4*T) / (pi * k^2 * (T^2 - t0^2))) * (cos(k*pi*t0/T) - (-1)^k);
bk = ((4*T)/(pi * k^2 * (T^2 - t0^2))) * sin(k*pi*t0/T);
end
k = -10:1:10;
ck(find(k>=0)) = (ak - 1i*bk)/2;
ck(find(k==0)) = ak(1);
ck(find(k<0)) = conj(fliplr(ck(find(k>0))));

plot(k,abs(ck),'bo')
hold on
for i=1:length(k), line([k(i) k(i)], [0 abs(ck(i))], 'Color', 'blue'); end
hold off
axis([k(1) k(end) -pi pi])
xlabel('k')
ylabel('abs(ck)')
legend('Widmo amplitudowe')
% Widzimy, ¿e elementy widma amplitudowego s¹ bliskie zera.

%% Widmo fazowe
% Widmem fazowym x(t) jest {arg(cn): n jest liczba calkowita}.
clear;
T = pi;
t0 = pi/3;
omega0 = pi/T;

a0 = 0;
for k=-10:1:10
    ak = ((4*T) / (pi * k^2 * (T^2 - t0^2))) * (cos(k*pi*t0/T) - (-1)^k);
    bk = ((4*T)/(pi * k^2 * (T^2 - t0^2))) * sin(k*pi*t0/T);
end
k = -10:1:10;
ck(find(k>=0)) = (ak - 1i*bk)/2;
ck(find(k==0)) = ak(1);
ck(find(k<0)) = conj(fliplr(ck(find(k>0))));

plot(k,angle(ck),'bo')
hold on
for i=1:length(k), line([k(i) k(i)], [0 angle(ck(i))], 'Color', 'red'); end
hold off
axis([k(1) k(end) -pi pi])
xlabel('k')
ylabel('angle(ck)')
legend('Widmo fazowe')
% Zauwazmy, ze dla dn = angle(cn): dn = -d(-n), poniewaz ak = a(-k) oraz 
% bk = b(-k) (arg(ak/2 - i bk/2) = -arg(a(-k)/2 + i b(-k)/2)).

%% (2)
clear;
% Odtwarzanie sygna³u na podstawie widma.
T = pi;
t0 = pi/3;
omega0 = pi/T;

a0 = 0;
k = 1:1:5000;
ak = ((4*T) ./ (pi.*k.^2 .* (T.^2 - t0.^2))) .* (cos(k.*pi.*t0./T) - (-1).^k);
bk = ((4.*T)./(pi .* k.^2 .* (T.^2 - t0.^2))) .* sin(k.*pi.*t0./T);

t = -3.5:0.001:3.5;
x = ak(1)/2 *ones(size(t));
for K = 2:1:length(k)
   x = x + ak(K) * cos(k(K)*omega0*t) + bk(K) * sin(k(K)*omega0*t); 
end
plot(t,x,'b-')
xlabel('t')
ylabel('x(t)')
legend('Sygna³ odtworzony')

%% Zadanie 2

% W pracy domowej otrzymalismy, ze transformata Fouriera sygnalu y(t)
% wynosi (2*T/(omega^2(T^2-t0^2)))*(exp(-j*omega*t0) - exp(-j*omega*T)).
clear;

syms t k T t0 n
assume(in(k,'integer'))
assume(in(t0,'real'))
assume(T>0)
assumeAlso(t0<T & t0>-T)
a = @(x,t,k,T) int(x*cos(k*pi*t/T)/T, t, -T, T);
b = @(x,t,k,T) int(x*sin(k*pi*t/T)/T, t, -T, T);
xs = @(x,t,n,T) a(x,t,0,T)/2 + ...
    symsum(a(x,t,k,T)*cos(k*pi*t/T) + b(x,t,k,T)*sin(k*pi*t/T), k, 1, n);

% Dla t z przedzialu [-T, t0):
% y(t) = t/(T + t0) + T/(T + t0)
% Dla t z przedzialu [t0, T):
% y(t) = t/(t0 - T) + T/(T-t0)
y = @(t0,T) piecewise(t<t0, t/(T + t0) + T/(T + t0), ...
    t>t0, t/(t0 - T) + T/(T-t0));

% Przy wykorzystaniu obliczen symbolicznych w MATLABie otrzymalismy, ze 
% tranformata Fouriera syganlu y(t) wynosi
% (643923891931320*T^2*cos((pi*t)/T) + 643923891931320*T^2*cos((pi*(t - t0))/T) + 1518580440096720*T^2*cos((pi*t)/T)^3 - 7761928852928250*T^2*cos((pi*t)/T)^4 - 27259582672498176*T^2*cos((pi*t)/T)^5 + 107632080093938400*T^2*cos((pi*t)/T)^6 + 239252401892275200*T^2*cos((pi*t)/T)^7 - 761112566378564400*T^2*cos((pi*t)/T)^8 - 1128760409839462400*T^2*cos((pi*t)/T)^9 + 3075202288398240000*T^2*cos((pi*t)/T)^10 + 3114026852197847040*T^2*cos((pi*t)/T)^11 - 7533313727094355200*T^2*cos((pi*t)/T)^12 - 5171565657277071360*T^2*cos((pi*t)/T)^13 + 11381136269422233600*T^2*cos((pi*t)/T)^14 + 5086953140239466496*T^2*cos((pi*t)/T)^15 - 10359173739502748160*T^2*cos((pi*t)/T)^16 - 2728666566033408000*T^2*cos((pi*t)/T)^17 + 5206749961914941440*T^2*cos((pi*t)/T)^18 + 614880809032089600*T^2*cos((pi*t)/T)^19 - 1109859860302921728*T^2*cos((pi*t)/T)^20 + 211688968716225*T^2*pi^2 + 1518580440096720*T^2*cos((pi*(t - t0))/T)^3 + 7761928852928250*T^2*cos((pi*(t - t0))/T)^4 - 27259582672498176*T^2*cos((pi*(t - t0))/T)^5 - 107632080093938400*T^2*cos((pi*(t - t0))/T)^6 + 239252401892275200*T^2*cos((pi*(t - t0))/T)^7 + 761112566378564400*T^2*cos((pi*(t - t0))/T)^8 - 1128760409839462400*T^2*cos((pi*(t - t0))/T)^9 - 3075202288398240000*T^2*cos((pi*(t - t0))/T)^10 + 3114026852197847040*T^2*cos((pi*(t - t0))/T)^11 + 7533313727094355200*T^2*cos((pi*(t - t0))/T)^12 - 5171565657277071360*T^2*cos((pi*(t - t0))/T)^13 - 11381136269422233600*T^2*cos((pi*(t - t0))/T)^14 + 5086953140239466496*T^2*cos((pi*(t - t0))/T)^15 + 10359173739502748160*T^2*cos((pi*(t - t0))/T)^16 - 2728666566033408000*T^2*cos((pi*(t - t0))/T)^17 - 5206749961914941440*T^2*cos((pi*(t - t0))/T)^18 + 614880809032089600*T^2*cos((pi*(t - t0))/T)^19 + 1109859860302921728*T^2*cos((pi*(t - t0))/T)^20 - 211688968716225*t0^2*pi^2)/(423377937432450*pi^2*(T^2 - t0^2))
% Jest to wynik zwracany przez kod simplify(xs(y(t0,T), t, 20,T)). 
% Znacz¹co odbiega od wyniku otrzymanego w pracy domowej.

%% Widmo amplitudowe
% Widmem amplitudowym sygnalu y(t) jest {|Y(omega)|: omega jest
% rzeczywiste}.

% Narysujemy wykres widma w przedziale [-10, 10].
% Przyjmujemy, ¿e T = 2 oraz t0 = -1/2.

T = 2;
t0 = -1/2;

Y = @(omega) (2.*T./(omega.^2.*(T.^2-t0.^2))).*(exp(-i.*omega*t0) - exp(-i.*omega*T));
omega = 0.1:0.1:5;
plot(omega, abs(Y(omega)))
title('Widmo amplitudowe sygna³u y(t)')
xlabel('omega')
ylabel('|Y(omega)|')
% Z wykresu wnioskujemy, ze amplituda sygnalu y(t) maleje wykladniczo do
% zera.
%% Widmo fazowe
% Widmem fazowym sygnalu y(t) jest {arg(Y(omega)): omega jest
% rzeczywiste}.

% Narysujemy wykres widma w przedziale [-10, 10].
% Przyjmujemy, ¿e T = 2 oraz t0 = -1/2.

T = 2;
t0 = -1/2;

Y = @(omega) (2.*T./(omega.^2.*(T.^2-t0.^2))).*(exp(-i.*omega*t0) - exp(-i.*omega*T));
omega = 0.1:0.1:40;
plot(omega,angle(Y(omega)))
title('Widmo fazowe sygna³u y(t)')
xlabel('omega')
ylabel('arg(Y(omega))')
% Z wykresu widzimy, ze faza sygnalu y(t) jest okresowa o okresie okolo 12.5.

%% Zadanie 3

% Wygenerujemy melodie "Wlazl kotek na plotek".
% Jej zapis nutowy to: sol, mi, mi, fa, re, re, do, mi, sol, sol, mi, mi,
% fa, re, re, do, mi, do, do, mi, mi, fa, re, re, do, mi, sol, sol, mi, 
% mi, fa, re, re, do, mi, do.

clear;

fs = 44100;
t = 0:1/fs:0.5;
step = nthroot(2,12);
% Kolejnoœæ nut: do, re, mi, fa, sol, la, si.
N1 = 440;              n1 = sin(2*pi*N1*t); %do
N2 = N1*(step^2);     n2 = sin(2*pi*N2*t); %re
N3 = N1*(step^3);    n3 = sin(2*pi*N3*t); %mi
N4 = N1*(step^4);     n4 = sin(2*pi*N4*t); %fa
N5 = N1*(step^5);     n5 = sin(2*pi*N5*t); %sol
N6 = N1*(step^6);    n6 = sin(2*pi*N6*t); %la
N7 = N1*(step^7);    n7 = sin(2*pi*N7*t); %si
Pause = zeros(1,floor(0.05*fs));

line1 = [n5,n3,n3,n4,n2,n2,n1,n3,n5,Pause];
line2 = [n5,n3,n3,n4,n2,n2,n1,n3,n1,Pause];
line3 = [n1,n3,n3,n4,n2,n2,n1,n3,n5,Pause];
line4 = [n5,n3,n3,n4,n2,n2,n1,n3,n1,Pause];
song = [line1, line2, line3, line4];
sound(song,fs)

% Teraz wygenerujemy przyblizone widmo tego sygnalu.
ck = fft(ifftshift(song))/length(song);

%% Pokazemy jego widmo amplitudowe.
k = (0:1:20000)*fs/length(song);
plot(k, abs(ck(1:length(k))), 'r.')
title('Przybli¿one widmo amplitudowe melodii')
hold on
for i = 1:length(k), line([k(i) k(i)], [0 abs(ck(i))], 'Color', 'red'); end
hold off
axis([k(1) k(end) -0.1 0.5])
xlabel('k')
ylabel('|ck|')

%% Pokazemy jego widmo fazowe.
k = (0:1:150)*fs/length(song);
plot(k, angle(ck(1:length(k))), 'r.')
title('Przybli¿one widmo fazowe melodii')
hold on
for i = 1:length(k), line([k(i) k(i)], [0 angle(ck(i))], 'Color', 'red'); end
hold off
axis([k(1) k(end) -0.1 0.5])
xlabel('k')
ylabel('arg(ck)')