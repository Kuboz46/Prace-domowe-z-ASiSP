% Jakub Zbrzezny, numer indeksu: 286689

%% Zadanie 1
%% (1)
% Przyjmujemy, ze f0 = 1/3, A = 3.
f0 = 1/3;
A = 3;
fs = 1e2;
t = 0:1/fs:10;
phi = rand(1,1)*2*pi - pi;
% ksi(n) = A*sin(2*pi*f0*tn+phi), tn = n/fs.
% Teraz wyswietlamy przykladowa realizacje sygnalu ksi(n).
plot(t,A*sin(2*pi*f0*t+phi))
%% (2)
% Tworzymy przyblizona funkcje autokorelacji sygnalu ksi(n) poprzez
% usrednienie naszego wyniku
It = 1e5;
autocor = zeros(1, 2*length(t)-1);
for iter = 1:It
    phi = rand(1,length(t))*2*pi - pi;
    ksi = A*sin(2*pi*f0*t+phi);
    autocor = autocor + xcorr(ksi,'unbiased');
end
[~,lags] = xcorr(ksi,'unbiased');
stem(lags,autocor/It)
% W pracy domowej otrzymalismy, ze jednowymiarowa funkcja autokorelacji w
% punkcie n wynosi A^2 * pi *cos(2*pi*f0*n/fs)
% Zatem wynik otrzymany poprzez usrednianie w MATLABie jest zblizony do
% wyniku otrzymanego w pracy domowej.
%% Zadanie 2
clear;
file = load('286689_delay.mat');
fs = 500;
syg = file.delayed;
delay = randi([50, 250])/1000;
td = (0:1:length(syg)-1)/fs;

delayed = [zeros(1,length(td)-1), file.delayed];
signal  = [file.signal, zeros(1,length(td)-1)];

[corr, lags] = xcorr(signal, delayed, 'unbiased');
plot(lags/fs,corr)
% Z wykresu funkcji korelacji widzimy, ze sygnal delayed jest opozniony o
% 40 sekund wzgledem sygnalu oryginalnego.

%% Zadanie 3
%% (1)
clear;
file = load('286689_ekg.mat');
ecg = file.ecg';
fs = 360;
t = (0:1:length(ecg) -1)/fs;
f01 = 0.5;
x1 = sin(2*pi*f01*t); 
x1N = x1 + 2*randn(1,length(t));

f02 = 50;
x2 = sin(2*pi*f02*t); 
x2N = x2 + 2*randn(1,length(t));

autocor = zeros(1,2*length(t)-1);
It = 200;
for iter = 1:It
    x1N = x1 + 2*randn(1,length(t)); % Skladowa o czestotliwosci 0.5 Hz.
    x2N = x2 + 2*randn(1,length(t)); % Skladowa o czestotliwosci 50 Hz.
    autocor = autocor + xcorr(x1N,'unbiased') + xcorr(x2N,'unbiased');
    
    subplot(3,1,1)
    plot(t,x1N + x2N)
    
    subplot(3,1,2)
    plot(t,ecg/iter)
    
    subplot(3,1,3)
    plot(t,ecg/iter - x1 - x2) % Usuwanie skladowej o czestotliwosci  
    % 0.5 Hz i skladowej o czestotliwosci 50 Hz.
    pause(1e-6)
end
% Ta petla dlugo trwa, okolo 9 minut.
%%
figure(2)
[~,lags] = xcorr(x1N + x2N,'biased');
stem(lags,autocor/It)

%% (2)
clear;
file = load('286689_ekg.mat');
ecg = file.ecg';
fs = 360;
t = (0:1:length(ecg)-1)/fs;
k = (-floor(length(ecg)/2):1:ceil(length(ecg)/2)-1)*fs/length(ecg);
ecgFFT = fftshift(fft(ecg));
plot(k,abs(ecgFFT));
axis([-75 75 0 300])
ecgFFT(abs(k)<=0.5) = 0;
ecgFFT(abs(k)==50) = 0;
ecg2 = ifft(ifftshift(ecgFFT));
%%
[autocor, lags] = xcorr(ecg2','biased');
plot(lags/fs,autocor)
%%
[~,i] = max(autocor.*(lags/fs >= 0.1));
T = lags(i)/fs;
F = 60/T; % Jest to srednia czestotliwosc tetna badanego pacjenta.
disp(F)
% Widzimy, ze wynosi ona 52.9412 1/min.