gamma = 1/6; %recovery rate
mu = 1/(365*65); %birth rate = death rate
beta = 5.3306e-05; %transmission rate
b = 0.1354; %protection rate
c = 2.5963e-06; %transmission rate (protected)
S0 = 18234;

rv_sivr = @(v) S0.*beta.*(1-v)./(gamma + mu);
rv_protected = @(v) ((beta + ((b*c)./(mu + v)))./(mu + gamma)) .*...
    (mu/(mu + b)) .* S0;
rv_quarantine = @(v) (beta/(mu + b + gamma)).*(mu./(v + mu)).*S0;

%sivr_rate = bisection(@(p) rv_sivr(p) - 1, 0, 1, 1e-8);
sivr_rate = 0.828382097184658;
%quarantine_rate = bisection(@(p) rv_quarantine(p) - 1, 0, 1, 1e-8);
quarantine_rate = 9.337812662124634e-05;

vax_rates = linspace(0, 1);
fontlabs = 'Times New Roman';

figure(1)
hold on
plot(vax_rates, rv_sivr(vax_rates))
plot(sivr_rate, 1, 'r*')
plot_0_title='Kermack-McKendrick Model w/ Vaccination';
xlabel('Vaccination Rate','FontSize',16,'FontName',fontlabs, ...
    'interpreter','latex');  
ylabel('$R_v$','FontSize',16,'FontName',fontlabs, ...
    'interpreter','latex'); 
title(plot_0_title,'FontSize',16,'FontName', ...
    'Times New Roman','interpreter','latex');
legend({'$R_v$', '$p = 0.8284$'}, 'Location',...
    'northeast', 'interpreter','latex')

figure(2)
plot(vax_rates, rv_protected(vax_rates))
plot_1_title='Protected Model';
xlabel('Vaccination Rate','FontSize',16,'FontName',fontlabs, ...
    'interpreter','latex');  
ylabel('$R_v$','FontSize',16,'FontName',fontlabs, ...
    'interpreter','latex'); 
title(plot_1_title,'FontSize',16,'FontName', ...
    'Times New Roman','interpreter','latex');

figure(4)
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
