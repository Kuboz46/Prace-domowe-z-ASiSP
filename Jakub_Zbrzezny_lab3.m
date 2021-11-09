%% Zadanie 2
disp('Por�wnamy wyniki liczenia splotu sygna��w {2, 8, 6, 6} i {6, 6, 8, 9} dla ka�dej z trzech metod z zadania 1.');
x = [2 8 6 6];
y = [6 6 8 9];
disp('W pracy domowej otrzymali�my, �e splot cykliczny wynosi {168,162,154,154}.');

disp('Testujemy funkcj� cconvSum.');
splot = cconvSum(x,y);
disp('Wynik wyszed�:');
disp(splot);
disp('Jest to wynik taki sam jak w pracy domowej.');
disp('Zatem ta funkcja dzia�a poprawnie.');

disp('Teraz sprawdzimy poprawno�� funkcji cconvMat.');
splot = cconvMat(x,y);
disp('Wynik wyszed�:');
disp(splot);
disp('Jest to wi�c znowu wynik taki sam jak w pracy domowej.');
disp('St�d ta funkcja r�wnie� dzia�a poprawnie.');

disp('A teraz przetestujemy funkcj� cconvDFT.');
splot = cconvDFT(x,y);
disp('Wynik wyszed�:');
disp(splot);
disp('Tu r�wnie� mamy taki sam wynik jak w pracy domowej.');
disp('Tak wi�c ta funkcja tak�e dzia�a poprawnie.');
disp('St�d wszystkie te zaimplementowane metody dzia�aj� prawid�owo.')
%% Zadanie 3
% Zadanie 3
disp('Por�wnamy z�o�ono�� obliczeniow� implementowanych metod i wbudowanej funkcji cconv w zale�no�ci od wymiaru N splatanych wektor�w.');
t = zeros(4,255);
N = 2:256;
for n = 2:256
    x = randi([0,9],1,n);
    y = randi([0,9],1,n);  
    tic
    cconvSum(x,y);
    t(1,n-1) = toc;
    
    tic
    cconvMat(x,y);
    t(2,n-1) = toc;
    
    tic
    cconvMat(x,y);
    t(3,n-1) = toc;
    tic
    
    cconv(x,y,n);
    t(4,n-1) = toc;
    tic
end
% Teraz zamieszczamy nasze wyniki na wykresie.
figure(1);
plot(N,t(1,:),'b.',N,t(2,:),'g.',N,t(3,:),'r.', N, t(4,:), 'c.')
xlabel('n')
ylabel('t')
legend('cconvSum','cconvMat','cconvDFT','cconv')
% Czas jest w sekundach.
axis([0 256 0 1e-3])

disp('Na wykresie widzimy, �e dla bardzo ma�ych n wbudowana metoda cconv dzia�a najwolniej. Znacznie szybsze s� zaimplementowane metody cconvSum, cconvMat, cconvDFT.');
disp('Natomiast dla wi�kszych n najwolniejszymi metodami s� cconvMat i cconvDFT.'); 
disp('Wyra�nie wida�, �e metoda cconvSum jest szybsza od tamtych metod. A najszybciej liczy wbudowana metoda cconv.');
disp('Wynika to st�d, �e z�o�ono�� obliczeniowa meotdy cconv jest liniowa, cconvSum kwadratowa, a cconvMat i cconvDFT rz�du n^3.');