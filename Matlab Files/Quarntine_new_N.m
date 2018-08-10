gamma = 1/6; %recovery rate
mu = 0;
%mu = 1/(365*65); %birth rate = death rate
beta = 5.3306e-05; %transmission rate
b = 0.1354; %protection rate
c = 2.5963e-06; %transmission rate (protected)
N = 4000;

rv_quarantine = @(v) (beta/(mu + b + gamma)).*(N./(v + mu));


quarantine_rate = 9.337812662124634e-05;

vax_rates = linspace(0, 1);
fontlabs = 'Times New Roman';

figure(1)
hold on
plot(vax_rates, rv_quarantine(vax_rates))
plot(quarantine_rate, 1, 'r*')
plot_3_title='Quarantine Model';
xlabel('Vaccination Rate','FontSize',16,'FontName',fontlabs, ...
    'interpreter','latex');  
ylabel('$R_v$','FontSize',16,'FontName',fontlabs, ...
    'interpreter','latex'); 
title(plot_3_title,'FontSize',16,'FontName', ...
    'Times New Roman','interpreter','latex');
legend({'$R_v$', '$p = 9.34*10^{-5}$'}, 'Location',...
    'northeast', 'interpreter','latex')
