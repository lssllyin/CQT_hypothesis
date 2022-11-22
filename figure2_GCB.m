clear;clc;close all;
% load the data
[dat,txt] = xlsread('Figure2.xlsx','uncat_cat');
sub   = txt(2:end,1);[idx1,ids]=unique(sub);
lnk25 = dat(:,2);
nEa   = dat(:,6)*4.184;%kJ/mol
cEa   = dat(:,12)*4.184;%kJ/mol

set(gcf,'position',[66,325,1175,511]);
% panel a
subplot(1,2,1);
plot(lnk25,nEa,'ok','Markerfacecolor','w','Markersize',12);
hold on;
plot(lnk25,cEa,'ok','Markerfacecolor','k','Markersize',12);
[na,nb,nr,np]=linearregression(lnk25(ids),nEa(ids),'k','--');
[ca,cb,cr,cp]=linearregression(lnk25,cEa,'k','none');
line([-50 -3],[nanmean(cEa)+1 nanmean(cEa)+1],'Linestyle','-','color',[0.5 0.5 0.5 0.3],'Linewidth',52);
set(gca,'fontsize',20);
% calculate Q10
Qc20 = [exp(nanmin(cEa)*1000*(1/293.15-1/303.15)/8.314)...
       exp(nanmax(cEa)*1000*(1/293.15-1/303.15)/8.314)]%Q10 between 20-30 for catalyst-reactions;;
Qn20 = [exp(nanmin(nEa)*1000*(1/273.15-1/283.15)/8.314)...
       exp(nanmax(nEa)*1000*(1/273.15-1/283.15)/8.314)]%Q10 between 20-30 for uncatalyst-reactions;
% create the legend text
nstr = strcat('r^2=',num2str(nr,'%4.2f'),', {\itp}<0.0001');
cstr = strcat('r^2=',num2str(cr,'%4.3f'), ', {\itp}=',num2str(cp,'%4.2f'));
text(-23,220,nstr,'fontsize',14);
text(-23,190,cstr,'fontsize',14);

axis([-50 -3 30 240]);
line([-30 -30],[180 240],'Linestyle','-','color','k');
line([-30 -3],[180 180],'Linestyle','-','color','k');

hl = legend('Uncatalysed (n=56)','Catalysed (n=21)');
set(hl,'box','off','fontsize',14);
set(hl,'position',[0.2521 0.7378 0.2213 0.2143]);

xlabel('ln(k{\fontsize{12}non}) at  25 \circC [s^{-1}]','fontsize',20);
ylabel('E{\fontsize{16}a} [kJ mol^{-1}]','fontsize',20);

text(-60, 255,'(a)','fontsize',20,'fontweight','bold');

x1 = [-2.5 -2.5 -2.5 -2.5];
y1 = 50:50:200;
Q10 = [2 4 8 15];
text(x1', y1' ,num2str(Q10'),'fontsize',20);
text(3, 120 ,'Q{\fontsize{12}10}','fontsize',20,'rotation',90);
set(gca,'position',[0.12 0.1797 0.3347 0.7453]);

%% panel b
subplot(1,2,2);
kb     = 1.38064852e-23;% Boltzmann constant,J K-1
h      = 6.626070040e-34;%Planck’s constant,J s
R      = 8.314;
eyring = @(p,x) log(kb/h)+log(x)...
                -(p(1)-x*p(2))./(R.*x);
tem    = 0:1:200;
tem    = tem+273.15;
lnk25  = dat(:,2);
dh     = dat(:,4)*1000*4.184;
ds     = dat(:,5)/298.15*1000*4.184;
% calculate rates based on Eyring Equation
for i = 1:length(dh)
    lnk(i,:)   = eyring([dh(i) ds(i)],tem);%lnk from -50 to 200 C
end
% do jackknife
for k = 1:size(lnk,2)
    %uncatalysed
    ujack = jackknife(@corr,lnk(ids,k),nEa(ids));
    onr(k,1) = mean(ujack);
    nse = sqrt(55/56*sum((ujack-onr(k,1)).^2));%standard error, Page 279
    nrL(k,1) = onr(k,1) - nse;
    nrU(k,1) = onr(k,1) + nse;
    % catalysted
    idx = ~isnan(cEa);
    cjack = jackknife(@corr,lnk(idx,k),cEa(idx));
    ocr(k,1) = mean(cjack);
    ose = sqrt(20/21*sum((cjack-ocr(k,1)).^2));
    crL(k,1) = ocr(k,1) - ose;
    crU(k,1) = ocr(k,1) + ose;
end
tem1 = tem'-273.15;
plot(tem1,onr,'--k','linewidth',2);
hold on;
plot(tem1,ocr,'-k','linewidth',2);
idx = tem>=263.15 & tem<=333.15;
fill([tem1;flipud(tem1)],[nrU;flipud(nrL)],...
    [0.5 .5 0.5],'Edgecolor','none','FaceAlpha',0.2);
fill([tem1(idx);flipud(tem1(idx))],[nrU(idx);flipud(nrL(idx))],...
    [0.5 0.5 0.5],'Edgecolor','none','FaceAlpha',0.5);

fill([tem1;flipud(tem1)],[crU;flipud(crL)],...
    [0.5 .5 0.5],'Edgecolor','none','FaceAlpha',0.2);
fill([tem1(idx);flipud(tem1(idx))],[crU(idx);flipud(crL(idx))],...
    [0.5 0.5 0.5],'Edgecolor','none','FaceAlpha',0.5);

xlabel('T{\fontsize{15}ref} [\circC]','fontsize',20);
ylabel('Correlation Coefficient','fontsize',20);

set(gca,'position',[0.63 0.1797 0.3347 0.7453]);
set(gca,'fontsize',20,'ytick',-0.9:0.3:0.3);
legend('Uncatalysed','Catalysed');
set(legend,'box','off','position',[0.6328    0.5367    0.1566    0.1301],'fontsize',16);
axis([0 200 -1 0.38]);
text(-45, 0.46,'(b)','fontsize',20,'fontweight','bold');

% h=gcf;
% set(h,'Units','Centimeters');
% set(h,'PaperPositionMode','auto','PaperSize',[12,5.5]);
% print(h, '-dpdf', 'Figure2.pdf');


