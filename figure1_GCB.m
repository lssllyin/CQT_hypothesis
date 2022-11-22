clear;clc;close all;

figure;
% set(gcf,'position',[-1509,85,1371,871]);
set(gcf,'position',[50,5,1371,871]);

subplot(2,3,2);
set(gca,'box','on');
x = 0.05:0.1:0.95;
y3= 0.5.*x+0.25;
y4= 1.2.*x-0.1;
plot(x,y3,'-','color',[0.5 0.5 0.5],'Linewidth',4);
hold on;
plot(x,y4,':','color',[0.5 0.5 0.5],'Linewidth',4);
axis([0 1 0 1]);
line([0.2 0.2],[0 0.35],'Linestyle','--','color','k');
line([0.5 0.5],[0 0.50],'Linestyle','--','color','k');
line([0.8 0.8],[0 0.85],'Linestyle','--','color','k');
plot(0.2,0.5.*0.2+0.25,'ok','Markersize',18,'Markerfacecolor','k');
plot(0.2,1.2.*0.2-0.1,'ok','Markersize',18,'Markerfacecolor','w');
plot(0.5,1.2.*0.5-0.1,'sk','Markersize',18,'Markerfacecolor','k');
plot(0.8,1.2.*0.8-0.1,'^k','Markersize',18,'Markerfacecolor','w');
plot(0.8,0.5.*0.8+0.25,'^k','Markersize',18,'Markerfacecolor','k');

text(0.22,0.05,'$$T_1$$','fontsize',18,'interpreter','latex');
text(0.52,0.05,'$$T_2$$','fontsize',18,'interpreter','latex');
text(0.82,0.05,'$$T_3$$','fontsize',18,'interpreter','latex');
ylabel('ln(R{\fontsize{15}s})','fontsize',20);
xlabel({'Temperature'},'fontsize',20);
set(gca,'xtick',[],'ytick',[]);
axis([0 1 0 1]);
set(gca,'position',[0.60    0.60    0.3040    0.3412]);
text(-0.08, 1.1,'(b)','fontsize',18,'fontweight','bold');

h=subplot(2,3,1);
set(gca,'box','on');
x = 0.1:0.1:0.9;
y = 0.9:-0.1:0.1;
plot(x,y,'-','color',[0.5 0.5 0.5],'Linewidth',4);
axis([0 1 0 1]);
hold on;
plot(0.2,0.8,'ok','Markersize',18,'Markerfacecolor','w');
plot(0.8,0.2,'ok','Markerfacecolor','k','Markersize',18);
text(0.23, 0.82, 'High Q\fontsize{15}10','fontsize',18);
text(0.62, 0.14, 'Low Q\fontsize{15}10','fontsize',18);

p1 = [0 -0.15];
p2 = [1 -0.15];
dp = p2-p1;
quiver(p1(1),p1(2),dp(1),dp(2),0,'color','k');
set(get(h,'children'),'clipping','off')% turn off clippings
text(0, -0.20, 'Recalcitrant','fontsize',18);
text(0.85, -0.20, 'Labile','fontsize',18);

xlabel('Carbon Quality [ln(R{\fontsize{15}s})]','fontsize',20);
ylabel({'Temperature Sensitivity';'Q\fontsize{15}10'},'fontsize',20,'position',[-0.07 0.5 0]);
set(gca,'xtick',[],'ytick',[]);
text(-0.05,0,'Low','fontsize',18,'rotation',90);
text(-0.05,0.85,'High','fontsize',18,'rotation',90);
text(-0.2,1.1,'(a)','fontsize',18,'fontweight','bold');
set(gca,'position',[0.1300    0.60    0.3040    0.3412]);


%% for data
dat = xlsread('Figure1.xlsx','Soil_Respiration');
% % save 'soil_respiration.mat' 'dat';
% % load soil_respiration.mat;
a   = dat(:,12);
lna = log(a);
Q10 = dat(:,13);
b   = log(Q10)/10;
expf= @(p,x) p(1)+p(2).*x;
for i=1:length(b)
    SR1(i,1) = expf([lna(i,1),b(i,1)],0);
    SR0(i,1) = expf([lna(i,1),b(i,1)],43);%43.5
    SR2(i,1) = expf([lna(i,1),b(i,1)],80);
end

subplot(2,3,4);
h1 = plot(SR1,Q10,'ok','Markerfacecolor',[.5 .5 .5]);
hold on;
[a1,b1,r1,p1]=linearregression(SR1,Q10,[.5 .5 .5],'-');
set(gca,'fontsize',20);
text1 = strcat('{\itr^2}=',num2str(r1,'%4.2f'),', {\itp}<0.001');
text(-1,4.8,text1,'fontsize',18);
axis([-2 3 2 5]);
ylabel('Q\fontsize{15}10','fontsize',20);
xlabel({'ln(R\fontsize{15}0)'},'fontsize',20);
text(-3, 5.3,'(c)','fontsize',18,'fontweight','bold');

subplot(2,3,5);
plot(SR0,Q10,'ok');
hold on;
[a1,b1,r1,p1]=linearregression(SR0,Q10,[.5 .5 .5],'-');
set(gca,'fontsize',20);
text3 = strcat('{\itr^2}=',num2str(r1,'%4.2f'),', {\itp}=',num2str(p1,'%4.2f'));
text(4,4.85,text3,'fontsize',18);
axis([3 8 2 5]);
xlabel({'ln(R{\fontsize{15}43})'},'fontsize',20);
ylabel('Q\fontsize{15}10','fontsize',20);
text(2.0, 5.3,'(d)','fontsize',18,'fontweight','bold');
set(gca,'position',[0.4042    0.1100    0.2039    0.3412]);

subplot(2,3,6);
tem = 0:1:60;tem=tem';
for k=1:length(tem)
    for i=1:length(b)
        SR(i,1) = expf([lna(i,1),b(i,1)],tem(k));
    end
    [r1,p1] = corrcoef(SR,Q10);
    r(k,1) = r1(1,2);
    p(k,1) = p1(1,2);
end
idx = find(p>0.05);
plot(tem,r,'k','linewidth',2);
hold on;
idx0 = find(tem ==0);
idx25 = find(tem ==25);
idx43 = find(tem ==43);
plot(tem(idx0),r(idx0),'ok','Markerfacecolor',[0.5 0.5 0.5],'Markersize',10);
plot(tem(idx25),r(idx25),'ok','Markersize',10,'Markerfacecolor','k');
plot(tem(idx43),r(idx43),'ok','Markersize',10,'Markerfacecolor','w');

set(gca,'fontsize',20);
axis([0 60 -0.8 0.4]);
set(gca,'ytick',-0.8:0.4:0.4);
text(0.5,-0.75,'$$ln(R_{0})$$','fontsize',18,'interpreter','latex');
text(27,-0.4,'$$ln(R_{25})$$','fontsize',18,'interpreter','latex');
text(41,-0.12,'$$ln(R_{43})$$','fontsize',18,'interpreter','latex');
xlabel('T{\fontsize{15}ref} [\circC]','fontsize',20);
ylabel('Correlation Coefficient','fontsize',20);
text(-18, 0.52,'(e)','fontsize',18,'fontweight','bold');
tem(idx(1))
tem(idx(end))
% h=gcf;
% set(h,'Units','Centimeters');
% set(h,'PaperPositionMode','auto','PaperSize',[13,9]);
% print(h, '-dpdf', 'Figure1.pdf');