I = 1:30;
[X, Y] = meshgrid(I, I);

load('60x60_aborted.mat');

conds = conds(I, I);
condests = condests(I, I);
times = times(I, I);
density = density(I, I);


figure()
mesh(X, Y, conds);
title('Konditionszahl der Steifigkeitmatrix, SVD-Berechnung');
xlabel('Anzahl Sinus-Basisfunktionen im Ort');
ylabel('Anzahl Legendre-Basisfunktionen in der Zeit');

figure()
mesh(X, Y, condests);
title('Konditionszahl der Steifigkeitmatrix, condest-Schätzung');
xlabel('Anzahl Sinus-Basisfunktionen im Ort');
ylabel('Anzahl Legendre-Basisfunktionen in der Zeit');

figure()
mesh(X, Y, times);
title('Dauer der Berechnung der Steifigkeitmatrix');
xlabel('Anzahl Sinus-Basisfunktionen im Ort');
ylabel('Anzahl Legendre-Basisfunktionen in der Zeit');

figure()
mesh(X, Y, density);
title('Dichte der Steifigkeitmatrix');
xlabel('Anzahl Sinus-Basisfunktionen im Ort');
ylabel('Anzahl Legendre-Basisfunktionen in der Zeit');
