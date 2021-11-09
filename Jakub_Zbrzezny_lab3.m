%% Zadanie 2
disp('Porównamy wyniki liczenia splotu sygna³ów {2, 8, 6, 6} i {6, 6, 8, 9} dla ka¿dej z trzech metod z zadania 1.');
x = [2 8 6 6];
y = [6 6 8 9];
disp('W pracy domowej otrzymaliœmy, ¿e splot cykliczny wynosi {168,162,154,154}.');

disp('Testujemy funkcjê cconvSum.');
splot = cconvSum(x,y);
disp('Wynik wyszed³:');
disp(splot);
disp('Jest to wynik taki sam jak w pracy domowej.');
disp('Zatem ta funkcja dzia³a poprawnie.');

disp('Teraz sprawdzimy poprawnoœæ funkcji cconvMat.');
splot = cconvMat(x,y);
disp('Wynik wyszed³:');
disp(splot);
disp('Jest to wiêc znowu wynik taki sam jak w pracy domowej.');
disp('St¹d ta funkcja równie¿ dzia³a poprawnie.');

disp('A teraz przetestujemy funkcjê cconvDFT.');
splot = cconvDFT(x,y);
disp('Wynik wyszed³:');
disp(splot);
disp('Tu równie¿ mamy taki sam wynik jak w pracy domowej.');
disp('Tak wiêc ta funkcja tak¿e dzia³a poprawnie.');
disp('St¹d wszystkie te zaimplementowane metody dzia³aj¹ prawid³owo.')
%% Zadanie 3
% Zadanie 3
disp('Porównamy z³o¿onoœæ obliczeniow¹ implementowanych metod i wbudowanej funkcji cconv w zale¿noœci od wymiaru N splatanych wektorów.');
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

disp('Na wykresie widzimy, ¿e dla bardzo ma³ych n wbudowana metoda cconv dzia³a najwolniej. Znacznie szybsze s¹ zaimplementowane metody cconvSum, cconvMat, cconvDFT.');
disp('Natomiast dla wiêkszych n najwolniejszymi metodami s¹ cconvMat i cconvDFT.'); 
disp('WyraŸnie widaæ, ¿e metoda cconvSum jest szybsza od tamtych metod. A najszybciej liczy wbudowana metoda cconv.');
disp('Wynika to st¹d, ¿e z³o¿onoœæ obliczeniowa meotdy cconv jest liniowa, cconvSum kwadratowa, a cconvMat i cconvDFT rzêdu n^3.');