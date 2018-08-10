function f = kermMcKend(b,S0,y)
% f = kermMcKend(a,b,d,S0,P0,y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Kermack-McKendrick model
%R0 = (b*S0)/y without vaccination
f = (b*S0)/y;
f1 = []

% R0 = f(1-p) with vaccination
% Let p be the percet of the population that is vaccinated, ranging from
% 0 percent to 100 percent, increasing the percentange vaccinated by 1/1000
% at each step.
p1 = 0:1000
p2 = p1/1000
r0 = [];
r0 = f*(1-p2)

hold on
plot(p2,r0,'b')
hold on

plot(0,f,'r*')
xlabel('Percentage of People Vaccinated')
ylabel('Basic Reproduction Number R')
axis([0 1 0 7])
legend('With Vaccination','Without Vaccination R0 = 5.9589','location','southwest')
end

