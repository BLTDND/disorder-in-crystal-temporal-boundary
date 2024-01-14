close all
clear all
clc
tic
%% input disorder
%% before
V_before=0.25;
V_after=-0.4;
disorder_C_before=0;               
disorder_C_after=0;
gamma_c0=0.3;
disorder_gamma=0; 
level=100;
sta_num=30;
Disorder_C_before=zeros(level,sta_num);
E_2=zeros(100,1000/2,sta_num,level);
E_t_stor_before=zeros(2,100,1000,sta_num,level);
%random disorder
progressBar=waitbar(0,"disorder_all");
for disorder_level=1:level
    waitbar(disorder_level/level, progressBar, sprintf('进度 %d/%d', disorder_level, level));
    disorder_C_before=1/level*disorder_level;
    %disorder_C_after=1/level*disorder_level;
    for count=1:sta_num
        [t_ratio,E_t_stor0,E_omiga_2]=cal_disorder(V_before,V_after,disorder_C_before,disorder_C_after,gamma_c0,disorder_gamma);
        Disorder_C_before(disorder_level,count)=t_ratio;
        E_t_stor_before(:,:,:,count,disorder_level)=E_t_stor0;
        E_2(:,:,count,disorder_level)=E_omiga_2;
    end 
end
close(progressBar);
save Disorder_C_before.mat
save('E_t_stor_before.mat', 'E_t_stor_before', '-v7.3');
save E_2.mat
load Disorder_C_before.mat
x=reshape(1:level,[level,1]);
y_mean=mean(Disorder_C_before,2);
yneg=max(Disorder_C_before,[],2);
ypos=min(Disorder_C_before,[],2);
T=y_mean;
R=1-T;
figure(1)
% errorbar(x,R,yneg/30,ypos/30,'o')
% xlabel('disorder level');ylabel('R&t'),title("disorder before")
% hold on
% errorbar(x,T,yneg/30,ypos/30,'o')
% legend('R', 'T');
% 绘制拟合曲线
x_level=1:level;
order =4;
coefficients_disorder_before = polyfit(x_level, T, order);


xFit_level= linspace(min(x_level), max(x_level), 100);
yFit_disorder_order_before = polyval(coefficients_disorder_before, xFit_level);
plot(x,T,"r*");
hold on 
plot(xFit_level,yFit_disorder_order_before,"r-");
xlabel("level");ylabel("T");
title("Disorder C before")
legend("T","T FIT");




%% after
V_before=0.25;
V_after=-0.4;
disorder_C_before=0;               
disorder_C_after=0;
gamma_c0=0.3;
disorder_gamma=0; 
level=100;
sta_num=30;
Disorder_C_after=zeros(level,sta_num);
E_2_=zeros(100,1000/2,sta_num,level);
E_t_stor_after=zeros(2,100,1000,sta_num,level);
%random disorder
progressbar=waitbar(0,"disorder order");
 for disorder_level=1:level
     %disorder_C_before=1/level*disorder_level;
     waitbar(disorder_level/level, progressbar, sprintf('进度 %d/%d', disorder_level, level));
     disorder_C_after=1/level*disorder_level;
     for count=1:sta_num
         [t_ratio,E_t_stor0,E_omiga_2]=cal_disorder(V_before,V_after,disorder_C_before,disorder_C_after,gamma_c0,disorder_gamma);
         Disorder_C_after(disorder_level,count)=t_ratio;
         E_t_stor_after(:,:,:,count,disorder_level)=E_t_stor0;
         E_2_(:,:,count,disorder_level)=E_omiga_2;
     end 
 end
close(progressbar);
save Disorder_C_after.mat
save('E_t_stor_after.mat', 'E_t_stor_after', '-v7.3');
% writematrix(E_t_stor_after,'E_t_stor_after.csv');
save E_2_.mat
load Disorder_C_after.mat
x=reshape(1:level,[level,1]);
y_mean=mean(Disorder_C_after,2);
yneg=max(Disorder_C_after,[],2);
ypos=min(Disorder_C_after,[],2);
T=y_mean;
R=1-T;
figure(1)
% errorbar(x,R,yneg/30,ypos/30,'o')
% xlabel('disorder level');ylabel('R&T'),title("disorder all")
% hold on
% errorbar(x,T,yneg/30,ypos/30,'o')
title("Disorder C after")
% legend('R', 'T');
% 绘制拟合曲线
x_level=1:level;
order =4;
coefficients_disorder_after= polyfit(x_level, T, order);


xFit_level= linspace(min(x_level), max(x_level), 100);
yFit_disorder_after = polyval(coefficients_disorder_after, xFit_level);
plot(x,T,"r*");
hold on 
plot(xFit_level,yFit_disorder_after,"r-");
xlabel("level");ylabel("T");
legend("T","T FIT");
title("Dicsorder C after");
toc