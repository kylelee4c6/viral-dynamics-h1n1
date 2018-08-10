gamma = 1/6; %recovery rate
mu = 0;
%mu = 1/(365*65); %birth rate = death rate
beta = 5.3306e-05; %transmission rate
b = 0.1354; %protection rate
c = 2.5963e-06; %transmission rate (protected)
S0 = 18234;
%S0 = 1000;
N = 4000;

rv_protected = @(v) ((beta + ((b*c)./(mu + v)))./(mu + gamma)) .*...
    (N/(mu + b));

vax_rates = linspace(0, 1);
fontlabs = 'Times New Roman';


figure(2)
plot(vax_rates, rv_protected(vax_rates))
plot_1_title='Protected Model';
xlabel('Vaccination Rate','FontSize',16,'FontName',fontlabs, ...
    'interpreter','latex');  
ylabel('$R_v$','FontSize',16,'FontName',fontlabs, ...
    'interpreter','latex'); 
title(plot_1_title,'FontSize',16,'FontName', ...
    'Times New Roman','interpreter','latex');