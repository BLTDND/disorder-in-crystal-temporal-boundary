close all
clear all
clc
tic
% debug
% V_before=0.25
% V_after=-0.4
% disorder_C_before=0
% disorder_C_after=0
% gamma_c0=0.1
% disorder_gamma=0
% [t_ratio_r]=cal_disorder(V_before,V_after,disorder_C_before,disorder_C_after,gamma_c0,disorder_gamma)
%% input disorder
V_before=0.25;
V_after=-0.4;
disorder_C_before=0;               
disorder_C_after=0;
gamma_c0=0.3;
disorder_gamma=0; 
level=100;
sta_num=5;
Disorder_C=zeros(level,sta_num);
E_2=zeros(100,1000/2,sta_num,level);
 n=0;
E_t_stor0_all=zeros(2,100,1000,sta_num,level);
%random disorder
progressBar = waitbar(0, 'disorder_C');
 for disorder_level=1:level
     waitbar(disorder_level/level, progressBar, sprintf('进度 %d/%d', disorder_level, level));
     disorder_C_before=1/level*disorder_level;
     disorder_C_after=1/level*disorder_level;
     for count=1:sta_num
         [t_ratio,E_t_stor0,E_omiga_2]=cal_disorder(V_before,V_after,disorder_C_before,disorder_C_after,gamma_c0,disorder_gamma);
         Disorder_C(disorder_level,count)=t_ratio;
         E_t_stor0_all(:,:,:,count,disorder_level)=E_t_stor0;
         E_2(:,:,count,disorder_level)=E_omiga_2;
     end 
 end
close(progressBar)
save Disorder_C
save('E_t_stor0_all.mat', 'E_t_stor0_all', '-v7.3');
save E_2
load("Disorder_C.mat")
x=reshape(1:level,[level,1]);
y_mean=mean(Disorder_C,2);
yneg=max(Disorder_C,[],2);
ypos=min(Disorder_C,[],2);
T=y_mean;
R=1-T;
figure(1)
% errorbar(x,R,yneg/30,ypos/30,'o')
% xlabel('disorder phase');ylabel('Reflactiion'),title("disorder all")
% hold on
% errorbar(x,T,yneg,ypos,'o')
% xlabel('disorder phase');ylabel('transmission'),title("相位均匀分布的无序")
title("Disorder C temporal")
% legend('R', 'T');
% 绘制拟合曲线
x_level=1:level;
order =4;
coefficients_disorder_temporal= polyfit(x_level, T, order);


xFit_level= linspace(min(x_level), max(x_level), 1000);
yFit_disorder_temporal = polyval(coefficients_disorder_temporal, xFit_level);
plot(x,T,"r*");
hold on 
plot(xFit_level,yFit_disorder_temporal,"r-");
xlabel("level");ylabel("T");
legend("T","T FIT");
title("Disorder C temporal");
%% sure disorder
V_before=0.25;
V_after=-0.4;
disorder_C_before=0;               
disorder_C_after=0;
gamma_c0=0.3;
disorder_gamma=0; 
level=100;
sta_num=5;
Disorder_C_order=zeros(level,sta_num);
E_2_=zeros(100,1000/2,sta_num,level);
E_t_stor_all=zeros(2,100,1000,sta_num,level);
progressbar=waitbar(0,"disorder_C_order");
for disorder_level=1:level
    waitbar(disorder_level/level, progressbar, sprintf('进度 %d/%d', disorder_level, level));
    disorder_C_before=1/level*disorder_level;
    disorder_C_after=1/level*disorder_level;
    for count=1:sta_num
        rng("shuffle");
        random_value=rand(1,1000)-0.5;
        [t_ratio,E_t_stor1,E_omiga_2]=cal_disorder_order(V_before,V_after,disorder_C_before,disorder_C_after,gamma_c0,disorder_gamma,random_value);
        Disorder_C_order(disorder_level,count)=t_ratio;
        E_t_stor_all(:,:,:,count,disorder_level)=E_t_stor1;
        E_2_(:,:,count,disorder_level)=E_omiga_2;
    end 
end
close(progressbar);
save Disorder_C_order.mat 
save('E_t_stor_all.mat', 'E_t_stor_all', '-v7.3');
save E_2_.mat
x=reshape(1:level,[level,1]);
y_mean=mean(Disorder_C_order,2);
yneg=max(Disorder_C_order,2);
ypos=min(Disorder_C_order,2);
T=y_mean;
R=1-y_mean;

figure(2)
title("Disorder C order")
% legend('R', 'T');
% 绘制拟合曲线
x_level=1:level;
order =4;
coefficients_disorder_order= polyfit(x_level, T, order);


xFit_level= linspace(min(x_level), max(x_level), 1000);
yFit_disorder_order= polyval(coefficients_disorder_order, xFit_level);
plot(x,T,"r*");
hold on 
plot(xFit_level,yFit_disorder_order,"r-");
xlabel("level");ylabel("T");
legend("T","T FIT");
title("Disorder C order");
% xlabel('disorder phase');ylabel('Reflactiion'),title("disorder partcial")
toc