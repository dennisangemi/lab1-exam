%% Compund pendulum
% Physics lab 1 exam
% 
% manca da calcolare rappresentare g con la propagazione degli errori.

% cleaning
clc
clear

% importing data
df1=readtable("..\data\exp-data-1.csv");
df2=readtable("..\data\exp-data-2.csv");
tools=readtable("..\data\tools.csv");


% count configurations
uc=unique(df2.configuration); % unique configuration
nc=length(uc);                % number of configs
%% Prima esplorazione e distribuzione delle misure indirette

% line chart
% linechart=figure;
% for i=1:nc
%     subplot(3,5,i)
%     plot(table2array(df2(df2.configuration==uc(i),"event")),table2array(df2(df2.configuration==uc(i),"time_ms")))        
%     title(strcat("d_{CM} = ", string(table2array(df1(df1.configuration==uc(i),"distance_cm"))-50), " cm"))
% end
% histogram

hst=figure;
for i=1:nc
    subplot(3,5,i)
    histogram(table2array(df2(df2.configuration==uc(i),"time_ms")))     
    % title(strcat("Configuration",string(uc(i)), "distanza dal CM = ", string(table2array(df1(df1.configuration==uc(i),"distance_cm"))-50), " cm"))
    title(strcat("d_{CM} = ", string(table2array(df1(df1.configuration==uc(i),"distance_cm"))-50), " cm"))
end
%% Qual è l'incertezza sui tempi?
% Per stimarla, calcolo la deviazione standard sui tempi per ogni configurazione

% calcolo la deviazione standard dei tempi misurati
for i=1:nc
    % calcolo la media
    tm(i)=round(mean(table2array(df2(df2.configuration==uc(i),"time_ms"))),0);

    % calcolo la deviazione
    deviation(i)=round(std(table2array(df2(df2.configuration==uc(i),"time_ms"))),0);
end

% creo tabella
dev = table(uc, tm', deviation','VariableNames',{'configuration','mean_ms','deviation'});

% calcolo rapporto tra deviazione/media)
dev.ratio = round((dev.deviation./dev.mean_ms).*100,2);
%% Rigetto di dati

% calcolo distanza dalla media in valore assoluto per ogni valore (scarti
% normalizzati t)

% join con dev per ottenere tempi medi (mean_ms)
df2 = join(df2,dev,"Keys","configuration");
df2 = df2(:,["configuration","event","time_ms","mean_ms","deviation"]);
df2.difference = abs(df2.time_ms - df2.mean_ms);

% calcola di quante deviazioni standard il valore è fuori
df2.out = round(df2.difference ./ df2.deviation,1);

% elimino i valori che sono fuori di più di 3.5 sigma
df2 = df2(df2.out <= 3.5,:)
% aggiorno variabili
uc=unique(df2.configuration); % unique configuration
nc=length(uc);                % number of configs

% histogram
hst=figure;
for i=1:nc
    subplot(3,5,i)
    histogram(table2array(df2(df2.configuration==uc(i),"time_ms")),15)     
    % title(strcat("Configuration",string(uc(i)), "distanza dal CM = ", string(table2array(df1(df1.configuration==uc(i),"distance_cm"))-50), " cm"))
    title(strcat("d_{CM} = ", string(table2array(df1(df1.configuration==uc(i),"distance_cm"))-50), " cm"))
end
%% Rigetto di dati

toDelete = df2.configuration == 14;
df2(toDelete,:) = [];
%% 
% Alla luce dell'applicazione del criterio di Chauvenet, ricalcolo media e deviazione

% delete old mean and deviation
df2 = df2(:,["configuration", "event", "time_ms"]);

% calcolo la deviazione standard dei tempi misurati
for i=1:nc
    % calcolo la media
    tm(i)=round(mean(table2array(df2(df2.configuration==uc(i),"time_ms"))),0);
    % calcolo la deviazione
    deviation(i)=round(std(table2array(df2(df2.configuration==uc(i),"time_ms"))),0);
end

% creo tabella
% dev = table(uc, tm', deviation','VariableNames',{'configuration','mean_ms','deviation'});

% calcolo rapporto tra deviazione/media)
% dev.ratio = round((dev.deviation./dev.mean_ms).*100,2)
%%
% aggiorno variabili
uc=unique(df2.configuration); % unique configuration
nc=length(uc);                % number of configs
l=1;                          % pendulum length
g=9.81;                       % gravitational acceleration
dg=0.01;                      % error gravitational acceleration
dt=tools.uncertainty(1);      % error t
dr=tools.uncertainty(2);      % errror distance

% overwrite distance error
dr = 0.002; % meters

% creating empty array
tm=zeros(nc,1);     % tempi medi
dtm=zeros(nc,1);    % error t
d=zeros(nc,1);      % distance from CM
dd=zeros(nc,1);     % error distance
r=[0:0.0001:0.5];   % theoretical distance
gc=zeros(nc,1);     % gravitational acceleration calculated
dgc=zeros(nc,1);    % error gc
regc=zeros(nc,1);   % relative error gravitational acceleration
cfrg=zeros(nc,1);   % position first significant digit gravitational acceleration
uomd=string(zeros(nc,1));   % uom distance
uomg=string(zeros(nc,1));   % uom g
uomt=string(zeros(nc,1));   % uom t

% multiply *2, converting ms2s and rounding
% df2.period_s=round(df2.time_ms.*2./1000,1)
df2.period_s=df2.time_ms*2/1000
dt = 2*dt
%% Analisi statistica

% test del chi quadro su una distribuzione
% scelgo arbitrariamente la configurazione 6, grafico distribuzione ed
% eseguo test chi quadro

test = table2array(df2(df2.configuration == 6,"time_ms"))
sigmatest = std(test)
meantest = mean(test)
% set number of bins
nb = 4;
bin = (1:4)'
% number of values for each bin
% nv = zeros(nb,1);
% nv(1) = lengthtest(test < meantest - sigmatest,:));
% nv(2) = lengthtest(test(test > meantest - sigmatest & test < meantest,:));
% nv(3) = lengthtest(test(test > meantest  & test < meantest + sigmatest,:));
% nv(4) = lengthtest(test(test > meantest + sigmatest,:));

% si poteva fare la stessa con histcounts 

% specifico edges
edges = [min(test) ,meantest - sigmatest, meantest, meantest + sigmatest, max(test) ];

% compute the count
N = histcounts(test,edges)
% probability
p = [0.16 0.34 0.34 0.16]

% valori attesi
va = Np
% inserire 
histfit(test)
%%


% calcolo gc2 con formula più corta (g calculated 2)
% for i = 1:nc
%     gc2(i) = (2.*pi./tm(i)).^2  .*  ((l.^2)/(12.*d(i)) + d(i));
% end


% per ogni configurazione calcolo la media del periodo e la sua deviazione
% standard



sigma_t = zeros(nc,1);

for i = 1:nc
    tm(i) = mean(table2array(df2(df2.configuration == uc(i),"period_s")));
    sigma_t(i) = std(table2array(df2(df2.configuration == uc(i),"period_s")));
end

o4 = table(uc,tm,sigma_t,'VariableNames',{'configuration','mean_period','sigma_t'})
% rounding
o4.mean_period = round(o4.mean_period,2);
% o4.sigma_t = round(o4.sigma_t,2);

% join o4 with df1 in order to obtain distance
o4 = join(o4,df1,"Keys","configuration");
o4.distance_m = (o4.distance_cm-50)/100;

% filter out distance_cm
o4 = o4(:,["configuration","mean_period","sigma_t","distance_m"]);

% compute teoretical period
o4.teo_period = teot(1,o4.distance_m);

tt=(2.*pi./sqrt(g)).*sqrt(((l.^2)./(12.*r))+r); % theoretical curve

% preview
o4
% plotting
plt1=figure;
plot(o4.distance_m,o4.mean_period,'o')
xlabel('Distanza dal CM (m)')
ylabel('Periodo T (s)')
xlim([0,0.5])
ylim([0 5])
hold on
plot(r,tt)
hold off
legend('data','theoretical curve')
plt2=figure;
errorbar(o4.distance_m,o4.mean_period,o4.sigma_t,o4.sigma_t,'.')
xlabel('Distanza dal CM (m)')
ylabel('Periodo T (s)')
xlim([0,0.5])
hold on
plot(r,tt)
hold off
ylim([1.4 2.3])
xlim([0.05 0.5])
legend('data','theoretical curve')
% calcolare chi quadro su questo grafico
% compute chi square
chio4 = sum(((o4.mean_period-o4.teo_period)./o4.sigma_t).^2)
% chi ridotto
chio4./(height(o4)-1)
%% Metodo dei minimi quadrati

d = o4.distance_m;
y = tm;
x = sqrt(  (repelem(l,length(d))'.^2)./(12.*d) + d   );

% plot y and x
% plot(x,y,'o')
% xlabel("x")
% ylabel("y")
% determino parametri per applicare metodo minimi quadrati

% numero di punti
n = length(tm);

delta = n.*sum(x.^2) - (sum(x)).^2;

% intercetta
a = (   (sum(x.^2).*sum(y))  - (sum(x).*sum(x.*y))   )./delta

% coefficiente angolare
b = (   (n.*sum(x.*y)) - (sum(x).*sum(y))   )./delta
plt3=figure;
errorbar(x,o4.mean_period,o4.sigma_t,o4.sigma_t,'.')
hold on
plot(0:1,a+b.*(0:1)')
hold off
legend("data","fit",'Location','southeast')
xlim([0.7 1])
ylim([1.45 2.1])
plt4=figure;
plot(x,y,'o')
hold on
plot(0:1,a+b.*(0:1)')
hold off
legend("data","fit",'Location','southeast')
xlim([0.7 1])
ylim([1.45 2.1])
% stimo g dal coefficiente angolare
% gmq sta per "g minimi quadrati"
gmq = (2.*pi./b).^2
% add coefficiente di pearson


%% Determino g (media)

% join df1 e df2
df2 = join(df2,df1,"Keys","configuration");
df2.distance_cm = df2.distance_cm - 50;
df2.distance_m = df2.distance_cm/100;
df2 = df2(:,["configuration","event","period_s","distance_m"]);

% calculate g
for i = 1:height(df2)
    df2.g_calc(i) = (2.*pi./df2.period_s(i)).^2  .*  ((l.^2)/(12.*df2.distance_m(i)) + df2.distance_m(i));
end

% preview 
df2

% calculate mean and standard deviation of g
for i = 1:nc
    g_calc_mean(i) = mean(table2array(df2(df2.configuration == uc(i),"g_calc")));
    g_sigma(i) = std(table2array(df2(df2.configuration == uc(i),"g_calc")));
end

o5 = table(uc,g_calc_mean',g_sigma','VariableNames',{'configuration','g','sigma'})
% rounding
o5.g = round(o5.g,2);
o5.sigma = round(o5.sigma,2);

% preview
o5
% chi quadro
cg = sum(((o5.g-repelem(9.81,height(o5))')./o5.sigma).^2)
% mean g
round(mean(o5.g),1)
mean(o5.sigma)
%%
% plotting
plt5=figure;
errorbar(o5.configuration,o5.g,o5.sigma,'.')
hold on
plot(0:length(o5.configuration)+1,repelem(9.81,length(0:length(o5.configuration))+1,1))
hold off
xlim([0 height(o5)+1])
ylim([7.5 11])
legend("Experimental data","Accepted value","Location","southeast")
xlabel("Configuration")
ylabel("g (m/s^2)")
%% Exporting

% converting array to string
% d=string(d);
% tm=string(tm);
% gc=string(gc);
% regc=string(regc);
% 
% for i=1:nc
%     d(i)=sprintf('%.3f',d(i));
%     tm(i)=sprintf('%.2f',tm(i));
%     gc(i)=sprintf('%.1f',gc(i));
%     regc(i)=sprintf('%.2f',regc(i));
% end

% generating output 1
% output1=array2table(cat(2,string(uc),d,gc,dgc,regc,uomg),"VariableNames",{'configuration','distance_CM','gravitational_acceleration','uncertainty','relative_error','uom'});
% output1=sortrows(output1,"distance_CM");
% output1=output1(:,{'configuration','gravitational_acceleration','uncertainty','relative_error','uom'}) %remove distance_CM

% generating output2
% output2=sortrows(array2table(cat(2,string(uc),d,dd,uomd,tm,dtm,uomt),"VariableNames",{'configuration','distance_CM','uncertainty_distance','uom_distance','time','uncertainty_time','uom_time'}),"distance_CM")
% generating output3
% output3=array2table(cat(2,gm,dgm,uomg(1),regm),"VariableNames",{'gravitational_acceleration','uncertainty','uom','relative_error'})
%%
% exporting csv
% writetable(output1,'..\data\output-data-1.csv','Delimiter',',','Encoding','UTF-8')
% writetable(output2,'..\data\output-data-2.csv','Delimiter',',','Encoding','UTF-8')
% writetable(output3,'..\data\output-data-3.csv','Delimiter',',','Encoding','UTF-8')

% exporting img
saveas(hst,'..\img\histogram.png');
saveas(plt1,'..\img\plot1.png');
saveas(plt2,'..\img\plot2.png');
saveas(plt3,'..\img\plot3.png');
saveas(plt4,'..\img\plot4.png');
saveas(plt5,'..\img\plot5.png');

% exporting mlx2m
mlxloc = fullfile(pwd,'livescript.mlx');
fileout = 'script.m';
matlab.internal.liveeditor.openAndConvert(mlxloc,fileout);
%% Functions

function chi = chi(x,y,sigma) 
    % x valore atteso
    % y valore osservato
    % sigma deviazione
    for i = 1:length(x)
        c(i) = ((y(i)-x(i))./sigma(i)).^2;
    end
    chi = sum(c);
end

function tt = teot(l,r)
    % l lunghezza pendolo
    % r distanza asse di rotazione - centro di massa
    tt=(2.*pi./sqrt(9.81)).*sqrt(((l.^2)./(12.*r))+r);
end