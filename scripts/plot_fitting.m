% dane = importdata('dane.txt');
%x = lsqcurvefit(@diff_fitting, [2.2, 0.05], dane(:, 1), dane(:, 2));
dane = importdata('dane.txt');
x = [2.1322, 0.04831];
zz = diff_fitting(x, dane(:,1));

figure(1)
set(gca, 'fontsize', 18)
plot(dane(:,1),dane(:,2), '.k')
hold on
plot(dane(:, 1), zz, 'r', 'LineWidth',2)
grid on
legend('Experimental data', 'AESHG', 'Location', 'best')
set(gca, 'xlim', [0 240])
xlabel('t [min]'); ylabel('\eta [%]');
%text(10,5.5,strcat('E_1 = ',num2str(x(1))),'Color','k','FontSize',18);
%text(10,4.9,strcat('E_2 = ',num2str(x(2))),'Color','k','FontSize',18);
print('-f1', 'Fit_Keiser_1','-dpng')