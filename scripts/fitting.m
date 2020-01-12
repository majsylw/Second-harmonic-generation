clc
clear all
close all

dane = importdata('dane.txt');

M = 2^10;
time = 100;

tau = linspace(0, time, M);
dtau = tau(2) - tau(1);

parameters_Mashinsky = [2.6019, 0.04377];
parameters_Keiser = [2.1322, 0.04831];

AESHG = diff_fitting_AESHG(parameters_Mashinsky);
ESHG = diff_fitting_ESHG(parameters_Mashinsky);
AWSHG = diff_fitting_AWSHG(parameters);
WSHG = diff_fitting_WSHG(parameters);

figure(1)
set(gca, 'fontsize', 18)
plot(dane(:,1),dane(:,2), '.k')
hold on
plot(tau*5, ESHG, 'b', 'LineWidth',2)
plot(tau*5, AESHG, 'r', 'LineWidth',2)
%plot(tau*5, WSHG, 'g', 'LineWidth',2)
%plot(tau*5, AWSHG, 'c', 'LineWidth',2)
grid on
legend('Experimental data','ESHG', 'AESHG', 'Location', 'northwest')
%legend('Experimental data', 'ESHG', 'AESHG', 'WSHG', 'AWSHG', 'Location', 'northwest')
set(gca, 'xlim', [0 240])
xlabel('t [min]'); ylabel('\eta [%]');
print('-f1', 'Dopasowanie_two_Keiser_short_with_pomp_Mashinsky','-dpng')

figure(2)
set(gca, 'fontsize', 18)
plot(dane(:,1),dane(:,2), '.k')
hold on
plot(tau*5, ESHG, 'b', 'LineWidth',2)
plot(tau*5, AESHG, 'r', 'LineWidth',2)
plot(tau*5, WSHG, 'g', 'LineWidth',2)
plot(tau*5, AWSHG, 'c', 'LineWidth',2)
grid on
legend('Data', 'ESHG', 'AESHG', 'WSHG', 'AWSHG', 'Location', 'best')
set(gca, 'xlim', [0 500])
xlabel('t [min]'); ylabel('\eta [%]');
print('-f2', 'Dopasowanie_all_Mashinsky','-dpng')
close all

figure(1)
set(gca, 'fontsize', 18)
plot(dane(:,1),dane(:,2), '.k')
hold on
plot(tau*5, ESHG, 'b', 'LineWidth',2)
plot(tau*5, AESHG, 'r', 'LineWidth',2)
grid on
legend('', 'ESHG', 'AESHG', 'Location', 'northwest')
set(gca, 'xlim', [0 240])
xlabel('t [min]'); ylabel('\eta [%]');
print('-f1', 'Dopasowanie_Mashinsky_2_short','-dpng')

figure(2)
set(gca, 'fontsize', 18)
plot(dane(:,1),dane(:,2), '.k')
hold on
plot(tau*5, ESHG, 'b', 'LineWidth',2)
plot(tau*5, AESHG, 'r', 'LineWidth',2)
grid on
legend('Data', 'ESHG', 'AESHG', 'Location', 'best')
set(gca, 'xlim', [0 500])
xlabel('t [min]'); ylabel('\eta [%]');
print('-f2', 'Dopasowanie_Mashinsky_2','-dpng')
