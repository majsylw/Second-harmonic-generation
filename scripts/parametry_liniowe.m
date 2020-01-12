clc
clear all
dane = importdata('dane.txt'); % wydajnosc eksperymentalna
n_eff = importdata('n_eff.txt'); % dlugosc fali, n_p, n_sh, alpha_p, alpha_sh

% wykresy dla wydajnosci
figure(1)
set(gca, 'fontsize', 18)
plot(dane(:,1),ones(1,length(dane(:,1)))* dane(401,2),'--','LineWidth',2,'Color','k');
text(10,5.,'5.37 %','Color','k','FontSize',20);
hold on
plot(dane(:,1),ones(1,length(dane(:,1)))* 0.07,'--','LineWidth',2,'Color','k');
text(10,0.4,'0.07 %','Color','k','FontSize',20);
plot(dane(:,1),ones(1,length(dane(:,1)))* 0.05,'--','LineWidth',2,'Color',[.5 .5 .5]);
plot(dane(:,1),ones(1,length(dane(:,1)))* 0.08,'--','LineWidth',2,'Color',[.5 .5 .5]);
vline(200,'--r','')
text(202,3.2,'200 min','Color','r','FontSize',14);
plot(dane(:,1),dane(:,2), 'xb')
grid on
xlabel('t [min]'); ylabel('\eta [%]');
set(gca, 'xlim', [0 240])
hold off
print('-f1','eksperymentalna_wydajnosc_linear','-dpng')

figure(2)
set(gca,'FontName','Times','FontSize',18)
plot(dane(:,1),ones(1,length(dane(:,1)))* dane(401,2),'--','LineWidth',2,'Color','k');
text(10,4.,'5.37 %','Color','k','FontSize',20);
hold on
plot(dane(:,1),ones(1,length(dane(:,1)))* 0.07,'--','LineWidth',2,'Color','k');
text(190,0.13,'0.07 %','Color','k','FontSize',20);
plot(dane(:,1),ones(1,length(dane(:,1)))* 0.05,'--','LineWidth',2,'Color',[.5 .5 .5]);
text(190,0.035,'0.05 %','Color',[.5 .5 .5],'FontSize',20);
plot(dane(:,1),ones(1,length(dane(:,1)))* 0.08,'--','LineWidth',2,'Color',[.5 .5 .5]);
text(10,0.15,'0.08 %','Color',[.5 .5 .5],'FontSize',20);
plot(dane(:,1),dane(:,2), 'xb')
xlabel('t [min]'); ylabel('\eta [%]');
set(gca, 'xlim', [0 240])
set(gca,'XTick',0:60:240);
set(gca,'Ytick',[1e-1 1e0]);
grid on
hold off
set(gca,'YScale','log')
print('-f2','eksperymentalna_wydajnosc_log','-dpng')

% wspolczynniki zalamania
f = figure(3);
set(gca, 'fontsize', 18)
plot(n_eff(find(n_eff(:,1)==1064),1),n_eff(find(n_eff(:,1)==1064),3),'or','LineWidth',2);
hold on
text(1070,1.48,'1.469254','Color','r','FontSize',16);
plot(n_eff(find(n_eff(:,1)==1064),1),n_eff(find(n_eff(:,1)==1064),2),'ob','LineWidth',2);
text(1070,1.44,'1.456047','Color','b','FontSize',16);
grid on
xlabel('\lambda [nm]'); ylabel('n_{eff}');
set(gca, 'ylim', [1.39 1.56])
legend('n_p', 'n_{sh}', 'Location', 'best')
line(n_eff(:,1),n_eff(:,3),'LineWidth',2,'Color','r');
plot(n_eff(:,1),n_eff(:,2),'LineWidth',2,'Color','b');
ax1 = gca; 
set(ax1,'XColor','b','YColor','k') 
ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','left',...
           'Color','none',...
           'XColor','r','YColor',get(f,'color'));
set(ax2, 'fontsize', 18)
set(ax2,'Xtick',[250 400:200:1000]);
set(ax2,'Ytick',[1.4:0.02:1.56]);
xlabel(ax2,'\lambda [nm]');
%set(ax2,'Xlabel','\lambda [nm]');
%xlabel(ax2,'\lambda [nm]')
xlim([250 1000])
hold off
print('-f3','wspolczynnik_zalamania','-dpng')

% tlumienie
lambda = 0.5:0.01:2.7;
figure(4)
set(gca, 'fontsize', 18)
alphaUV = (201.3*lambda.^-4)*1e-3;
alphaIR = (3.947e10*exp(-42.73./lambda))*1e-3;
plot(lambda*10^3,alphaUV,'b--','LineWidth',2);
hold on
plot(lambda(find(lambda==0.53))*10^3,alphaUV(find(lambda==0.53))+alphaIR(find(lambda==0.53)),'ok','LineWidth',4);
text(600,2.5,'2.5512','Color','k','FontSize',16);
plot(lambda(find(lambda==1.06))*10^3,alphaUV(find(lambda==1.06))+alphaIR(find(lambda==1.06)),'ok','LineWidth',4);
text(1070,0.22,'0.1594','Color','k','FontSize',16);
plot(lambda*10^3,alphaIR,'r--','LineWidth',2);
plot(lambda*10^3,alphaUV+alphaIR,'k','LineWidth',3);
hold off
set(gca,'YScale','log')
xlabel('\lambda [nm]')
ylabel('\alpha [dB/m]')
grid on
set(gca,'XLim',[500 2500]);
set(gca,'YLim',[1e-2 1e1]);
set(gca,'YTick',[1e-1 1e0]);
print('-f4','tlumienie','-dpng')