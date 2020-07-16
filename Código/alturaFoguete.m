% Autor: Felipe Bogaerts de Mattos
% Equipe Supernova Rocketry UFJF
% Codigo compativel com a planilha SRM do Nakka.
% Para utilizar, salve a planilha com o nome 'Calculos(Nakka).xlsx' na
% mesma pasta do codigo.
% As EDOs assumem gravidade constante, diminuicao linear da densidade do ar
% com a altura (preciso ate 4000 metros), tempo de voo maximo de 60
% segundos, coeficiente de arrasto constante, projecao vertical circular, 
% voo vertical, voo sem paraquedas

%% Leitura dos dados da planilha

clear all; close all; clc;

% tf = tempo final de queima
tBurnout = xlsread('Calculos(Nakka).xlsx','Performance','L910:L910');
initial = [0;0;0;0];

% T a leitura do tempo na planilha
T = xlsread('Calculos(Nakka).xlsx','Performance','L28:L910');
% FT a leitura da forca fornecida pelo motor
FT = xlsread('Calculos(Nakka).xlsx','Performance','J28:J910');
[rows1, columns1] = size(FT);
[rows2, columns2] = size(T);
% Mp a leitura da massa do propelente
Mp = xlsread('Calculos(Nakka).xlsx','Pressure','U28:U862');

dominioT = 60;

i = 862;
while i <= 882
    Mp(i+1) = 0;
    i = i + 1;
end

intervaloT = 0.001;
% Para que as matrizes tenham mesma dimensao & estender o dominio t 
% do grafico:
while T < dominioT -2*intervaloT
    FT(end+1) = 0;
    Mp(end+1) = 0;
    T(end+1) = T(end) + intervaloT;
    rows1=rows1+1;
    rows2=rows2+1;
end

clear i rows1 rows2 intervaloT columns1 columns 2 dominioT;


%% Variaveis de input:

% Coeficiente de arrasto, raio do foguete (m) e altura inicial (m):
Cf = 0.95; r = (67e-3)/2; h0 = 0;
% Massa do foguete sem o motor (Kg) e massa da carga util (Kg):
Mf = 2;
Mc = 0.8;
Md = Mf + Mc;
% Massa do motor vazio (Kg):
Mm = 0.85;


%% Resolucao das EDOs:

[t,f] = ode45(@(t,f) odes(t,f,T,FT,Mp,Cf,r,h0,Md,Mm),[0 T(end)],initial);

[Vmax, i2] = max(f(:,1));
%i1 = indice da altura maxima na matriz H
[h, i1] = max(f(:,2));
%tA = tempo ate apogeu
tempoApogeu = t(i1);

clear i1 i2;


%% Indice da matriz f

j = 1;
while f(j,2) > -0.1
    j = j+1;
end

tempoVoo = t(j);
lenghtf = size(f);
f(j:lenghtf(1,1),:) = [];
t(j:lenghtf(1,1),:) = [];
clear lenghtf;


%% Distancia da base de lancamento

b = j-1;
Distancia = sqrt((f(b,3)^2)+(f(b,4)^2));


%% Forca de arrasto instantanea

A = pi*r^2; 
a = 0.2/1500; p0 = 1.2;

lenghtf = size(f);
j = 1;
while j <= lenghtf(1,1)
    moduloV = (f(j,1)/abs(f(j,1)));
    Drag(j,1) = moduloV*(A*Cf*(p0-a*(f(j,2)+h0)))*(f(j,1)^2)/2;
    j = j + 1;
end

clear moduloV lenghtf;


%% Graficos:

figure(1)
plot(t,f(:,1),'r');
text(tBurnout,max(f(:,1)),sprintf('     Vmax: %.1f m/s', Vmax))
grid on
hold on
plot(t,f(:,2),'k-.');
hold on
text(tempoApogeu,max(f(:,2)),sprintf('         Apogeu: %.0f metros', h))
xlabel('Tempo (s)');
legend('Veloc. (m/s)','Altura (m)','Location','Northwest')
grid on

saveas(gcf,'apogeu.png')

figure(2)
plot3(f(:,3),f(:,4),f(:,2),'b--o')
text(Distancia,f(b,2),sprintf('  Pouso: %.0f m', Distancia))
grid on
title('Trajetoria')
xlabel('Distancia ao Norte da base de lancamento (m)')
zlabel('Altura (m)')
ylabel('Distancia ao Leste da base de lancamento (m)')

saveas(gcf,'trajetoria.png')

figure(3)
grid on
yyaxis left
plot(t,Drag)
ylabel('Forca de arrasto (N)')
yyaxis right
plot(t,f(:,1))
ylabel('Velocidade (m/s)')
xlabel('Tempo (s)')

saveas(gcf,'forca_arrasto.png')

clc


%% Prints dos resultados:
fprintf('Altura final: %.2f metros.\n\n',h);
fprintf('Velocidade maxima: %.2f metros/segundo.\n\n',Vmax);
fprintf('Velocidade maxima: %.2f km/h.\n\n',Vmax*3.6);
fprintf('Velocidade maxima: %.2f Mach.\n\n',Vmax/340);
fprintf('Tempo ate apogeu: %.2f segundos.\n\n',tempoApogeu);
fprintf('Tempo final de queima: %.2f segundos.\n\n',tBurnout);
fprintf('Massa inicial do propelente: %.3f Kg.\n\n',Mp(1,1));
fprintf('Massa inicial do foguete: %.3f Kg.\n\n',Md + Mm + Mp(1,1));
fprintf('Tempo de voo: %.1f segundos.\n\n',tempoVoo);
fprintf('Distancia da base no pouso: %.0f metros.\n\n',Distancia);

%clear tBurnout tempoApogeu tempoVoo T FT;
%clear all;


function rk = odes(t,f,T,FT,Mp,Cf,r,h0,Md,Mm);

format long

FT = interp1(T,FT,t);
Mp = interp1(T,Mp,t);
% Massa total instantanea:
M = Md + Mm + Mp;
% Funcao da densidade do ar:
a = 0.2/1500;
p0 = 1.2;
% Valores constantes:
g = 9.80665; angulo_N = degtorad(89); angulo_L = degtorad(89);
A = pi*r^2; D = (A*Cf*(p0 - a*(f(2) + h0)))/2;

if f(1) < 0
    x = -1; k = 0;
else
    x = 1; k = pi;
end

rk(1) = (FT -x*D*f(1)^2)/(M) -g*sin(angulo_N);

rk(2) = f(1)*sin(angulo_N)*sin(angulo_L);

rk(3) = f(2)*cos(angulo_N);

rk(4) = f(2)*cos(angulo_L);

rk = rk(:);

% f(1) corresponde a velocidade; 
% f(2) corresponde a altura;
% f(3) corresponde a distancia ao Norte;
% f(4) corresponde a distancia ao Leste;

end


