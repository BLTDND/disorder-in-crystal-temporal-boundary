
clear all

close all
clc
load E_t_stor_before.mat
[a,b,c,d,e]=size(E_t_stor_before);
E_t_stor=mean(E_t_stor_before,4);
figure


for level_=1:e
    E_cal_co=abs(reshape(E_t_stor(1,:,:,:,level_)+E_t_stor(2,:,:,:,level_),[],1));
%     E_cal_co=E_cal_co(1:1000:end);
    x_l=(1:2*length(E_cal_co)-1)-length(E_cal_co)/2;
%     for tau= 1:length(E_cal_co)%xindex=y_index+tau
%         numerator_sum=0;
%         denominator_sum=0;
%         tau_=tau-length(E_cal_co)/2;
%         for yindex=1:length(E_cal_co)
%             if 0<tau_+yindex && tau_+yindex<length(E_cal_co)
%                 numerator_sum=numerator_sum+E_cal_co(tau_+yindex)'*E_cal_co(yindex)*1000*0.075;    
%             end
%             denominator_sum=denominator_sum+sqrt(abs(E_cal_co(yindex)))*1000*0.075;
%         end    
%         E_co_time(tau)=numerator_sum/denominator_sum; 
%     end
    E_co_time=xcorr(E_cal_co, 'coeff');
    subplot(1,2,1),plot(x_l,E_co_time,"*"); title("g^1")
    hold on
end
for level_=1:e
     E_cal_co=abs(reshape(E_t_stor(1,:,:,:,level_)+E_t_stor(2,:,:,:,level_),[],1));
    E_cal_co=E_cal_co(1:100:end);
    x_l=(1:2*length(E_cal_co)-1)-length(E_cal_co)/2;
    N=length(E_cal_co);
    E_co_time2=zeros(2*N-1,1);
    for k=1:N
        E_co_time2(k:N+k-1) = E_co_time2(k:N+k-1) + E_cal_co(k)*E_cal_co(1:N);
    end
    E_co_time2=E_co_time2/N;
%     for tau= 1:length(E_cal_co)%xindex=y_index+tau
%         numerator_sum=0;
%         denominator_sum=0;
%         tau_=tau-length(E_cal_co)/2;
%         times=0;
%         for yindex=1:length(E_cal_co)
%             if 0<tau_+yindex && tau_+yindex<length(E_cal_co)
%                 numerator_sum=numerator_sum+E_cal_co(tau_+yindex)'*E_cal_co(tau_+yindex)*E_cal_co(yindex)'*E_cal_co(yindex)*1000*0.075;
%                 
%             end
%             denominator_sum=denominator_sum+(abs(E_cal_co(yindex))^2)*1000*0.075;
%         end
%         E_co_time2=numerator_sum/(denominator_sum)^2; 
%     end
    subplot(1,2,2),plot(x_l,E_co_time2,"*"); title("g^2")
    hold on
end