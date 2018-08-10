N = 18234; %total population
[t_s, y_s] = ode45(@sir_sol,[0 100],[18223 11 0 0]);
[t_q, y_q] = ode45(@quarantine_sol,[0 100],[18223 11 0 0 0]);
[t_p, y_p] = ode45(@protected_sol,[0 100],[18223 11 0 0 0]);
[t_pa, y_pa] = ode45(@protected_alt_sol,[0 100],[18223 11 0 0 0]);
data = xlsread('h1n1InfectionPrevalence.csv');

fontlabs = 'Times New Roman';

%SIVR Model Plot
figure(1)
hold on
plot(t_s, y_s(:,2)./N, 'm-')
% plot(data(:,1), data(:,2), 'ko')
plot_0_title='Kermack-McKendrick Model w/ Vaccination';
xlabel('Days Since Outbreak','FontSize',16,'FontName',fontlabs, ...
    'interpreter','latex');  
ylabel('Infection Prevalence','FontSize',16,'FontName',fontlabs, ...
    'interpreter','latex'); 
title(plot_0_title,'FontSize',16,'FontName', ...
    'Times New Roman','interpreter','latex');
% legend({'Model', 'Data'}, 'Location',...
%     'northeast', 'interpreter','latex')
% axis([0 100 0 0.5]);

%Quarantine Model Plot
figure(2)
hold on
plot(t_q, y_q(:,2)./N, 'r-')
% plot(data(:,1), data(:,2), 'ko')
plot_1_title='Quarantine Model';
xlabel('Days Since Outbreak','FontSize',16,'FontName',fontlabs, ...
    'interpreter','latex');  
ylabel('Infection Prevalence','FontSize',16,'FontName',fontlabs, ...
    'interpreter','latex'); 
title(plot_1_title,'FontSize',16,'FontName', ...
    'Times New Roman','interpreter','latex');
% legend({'Model', 'Data'}, 'Location',...
%     'northeast', 'interpreter','latex')
% axis([0 100 0 0.3]);

%Protected Model Plot
figure(3)
hold on
plot(t_p, y_p(:,2)./N, 'b-')
plot(data(:,1), data(:,2), 'ko')
plot_2_title='Protected Model';
xlabel('Days Since Outbreak','FontSize',16,'FontName',fontlabs, ...
    'interpreter','latex');  
ylabel('Infection Prevalence','FontSize',16,'FontName',fontlabs, ...
    'interpreter','latex'); 
title(plot_2_title,'FontSize',16,'FontName', ...
    'Times New Roman','interpreter','latex');
axis([0 100 0 0.04]);
legend({'Model', 'Data'}, 'Location',...
    'northeast', 'interpreter','latex')

%Alt. Protected Model Plot
figure(4)
hold on
plot(t_pa, y_pa(:,2)./N, 'g-')
% plot(data(:,1), data(:,2), 'ko')
plot_3_title='Alternative Protected Model';
xlabel('Days Since Outbreak','FontSize',16,'FontName',fontlabs, ...
    'interpreter','latex');  
ylabel('Infection Prevalence','FontSize',16,'FontName',fontlabs, ...
    'interpreter','latex'); 
title(plot_3_title,'FontSize',16,'FontName', ...
    'Times New Roman','interpreter','latex');
% legend({'Model', 'Data'}, 'Location',...
%     'northeast', 'interpreter','latex')
% axis([0 100 0 0.04]);

%Combined Plots
figure(5)
hold on
plot(t_s, y_s(:,2)./N, 'm-')
plot(t_q, y_q(:,2)./N, 'r-')
plot(t_p, y_p(:,2)./N, 'b-')
plot(t_pa, y_pa(:,2)./N, 'g-')
plot(data(:,1), data(:,2), 'ko')
plot_4_title= 'All Models';
xlabel('Days Since Outbreak','FontSize',16,'FontName',fontlabs, ...
    'interpreter','latex');  
ylabel('Infection Prevalence','FontSize',16,'FontName',fontlabs, ...
    'interpreter','latex'); 
title(plot_4_title,'FontSize',16,'FontName', ...
    'Times New Roman','interpreter','latex');
legend({'SIVR', 'Quarantine', 'Protected', 'Alt. Protected', 'Data'}, ...
    'Location', 'northeast', 'interpreter','latex')
axis([0 100 0 0.5]);

function dydt = sir_sol(t, y)
gamma = 1/6; %recovery rate
mu = 0.0042; %birth rate = death rate
beta = 5.45e-5; %transmission rate
p = 0;
%y(1) = S, y(2) = I, y(3) = V, y(4) = R
dydt(1) = (1-p)*mu - mu*y(1) - beta*y(2)*y(1);
dydt(2) = beta*y(2)*y(1) - mu*y(2) - gamma*y(2);
dydt(3) = (p-y(3))*mu;
dydt(4) = gamma*y(2) - mu*y(4);
dydt = dydt';
end

function dydt = quarantine_sol(t, y)
gamma = 1/6; %recovery rate
mu = 0.0042; %birth rate = death rate
beta = 5.45e-5; %transmission rate
b = 0.15; %quarantine rate
p = 0;
%y(1) = S, y(2) = I, y(3) = Q, y(4) = R, y(5) = V
dydt(1) = mu - beta*y(2)*y(1) - mu*y(1) - p*y(1);
dydt(2) = beta*y(2)*y(1) - b*y(2) - gamma*y(2)- mu*y(2);
dydt(3) = b*y(2) - mu*y(3) - gamma*y(3);
dydt(4) = gamma*y(2) + gamma*y(3) - mu*y(4);
dydt(5) = p*y(1) - mu*y(5);
dydt = dydt';
end

function dydt = protected_sol(t, y)
gamma = 1/6; %recovery rate
mu = 0.0042; %birth rate = death rate
beta = 5.45e-5; %transmission rate
b = 0.15; %protection rate
c = 3.74e-6; %transmission rate (protected)
p = 0;
%y(1) = S, y(2) = I, y(3) = P, y(4) = V, y(5) = R
dydt(1) = mu - b*y(1) - beta*y(2)*y(1) - mu*y(1) - p*y(1);
dydt(2) = beta*y(2)*y(1) + c*y(3)*y(2) - gamma*y(2)-mu*y(2);
dydt(3) = b*y(1) - c*y(3)*y(2)- mu*y(3);
dydt(4) = mu*p - mu*y(4);
dydt(5) = gamma*y(2) - mu*y(5);
dydt = dydt';
end

function dydt = protected_alt_sol(t, y)
gamma = 1/6; %recovery rate
mu = 0.0042; %birth rate = death rate
beta = 5.45e-5; %transmission rate
b = 0.15; %protection rate
c = 0; %transmission rate (protected)
p = 0;
%y(1) = S, y(2) = I, y(3) = P, y(4) = V, y(5) = R
dydt(1) = mu - b*y(1) - beta*y(2)*y(1) - mu*y(1) - p*y(1);
dydt(2) = beta*y(2)*y(1) + c*y(3)*y(2) - gamma*y(2)-mu*y(2);
dydt(3) = b*y(1) - c*y(3)*y(2)- mu*y(3);
dydt(4) = mu*p - mu*y(4);
dydt(5) = gamma*y(2) - mu*y(5);
dydt = dydt';
end