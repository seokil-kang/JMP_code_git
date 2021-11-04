clear; clc; close all;
load data_LTW
load data_Cochrane

lw = 1.5;
Cindx = 2;
Lindx = 3;
figure('name','Data comparison between Cochrane(2021) and Leeper, Traum, Walker(2017)','color','w')
subplot(2,2,1)
plot(data_Cochrane.data_date,data_Cochrane.data(:,Cindx)*100,'x-','linewidth',lw)
hold on
plot(data_LTW.data_date,data_LTW.data(:,Lindx),'d-','linewidth',lw)
%legend(data_Cochrane.var_name(Cindx),data_LTW.var_name(Lindx))
%legend boxoff
title('output growth')

subplot(2,2,2)
Cindx = 3;
Lindx = 10;
plot(data_Cochrane.data_date,data_Cochrane.data(:,Cindx)*100,'x-','linewidth',lw)
hold on
plot(data_LTW.data_date,data_LTW.data(:,Lindx),'d-','linewidth',lw)
%legend(data_Cochrane.var_name(Cindx),data_LTW.var_name(Lindx))
%legend boxoff
title('inflation')

subplot(2,2,3)
Cindx = 6;
Lindx = 12;
plot(data_Cochrane.data_date,data_Cochrane.data(:,Cindx)*25,'x-','linewidth',lw)
hold on
plot(data_LTW.data_date,data_LTW.data(:,Lindx),'d-','linewidth',lw)
%legend(data_Cochrane.var_name(Cindx),data_LTW.var_name(Lindx))
%legend boxoff
title('interest rate')

subplot(2,2,4)
Cindx = 5;
Lindx = 9;
plot(data_Cochrane.data_date(2:end),diff(data_Cochrane.data(:,Cindx)*100),'x-','linewidth',lw)
hold on
plot(data_LTW.data_date,data_LTW.data(:,Lindx),'d-','linewidth',lw)
%legend(data_Cochrane.var_name(Cindx),data_LTW.var_name(Lindx))
%legend boxoff
legend('Cochrane(2021)','Leeper et al(2017)','orientation','horizontal','location','north')
legend boxoff
title('debt value')