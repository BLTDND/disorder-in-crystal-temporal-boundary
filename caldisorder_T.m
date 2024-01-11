load Disorder_C_before.mat
load Disorder_C.mat
load Disorder_C_after.mat

y_mean1=mean(Disorder_C_before,2);
y_mean2=mean(Disorder_C_after,2);
y_mean3=mean(Disorder_C,2);
figure
plot(y_mean1,"*")
hold on
plot(y_mean2,"*")
hold on
plot(y_mean3,"ob")
legend("before","after","all")
