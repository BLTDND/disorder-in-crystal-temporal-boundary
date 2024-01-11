function [t_ratio,E_t_stor,E_omiga_2]=cal_disorder(V_before,V_after,disorder_C_before,disorder_C_after,gamma_c0,disorder_gamma)
% V_before - 相位偏移平均值-50round前
% V_after - 相位偏移平均值-50round后
% disorder_C_* - 突变前后相位无序宽度
% disorder   - 耦合无序宽度
% gamma_c0 - 耦合平均值
format long
L= 15; % m
c= 0.3; %m/ns 
T_RT = 75; %ns
delta_t=0.075; %ns corresponds to a switching rate of 13.33 GHz
n_g=1.5; %typical value for silica optical fibers at telecom wavelengths
v_g=c/n_g; % the group velocity
 
k0=pi/3;%the central carrier wave frequency
Omega=2*pi*v_g/L; % FSR e9
omega0=k0*v_g;
N_z=1000;
m=linspace(-50,50,N_z);%for standing wave
omega_m=m.*Omega+omega0;
N_omega_m=length(omega_m);
d_omega=omega_m(2)-omega_m(1);
sigma_f = 2*Omega;
sigma_T = v_g/sigma_f;


kai=gamma_c0/delta_t; 
N_in_trip=T_RT/delta_t;

%% initial state 
C=@(V) V/(2*delta_t);
phi_plus=@(k,C) [2*C.*cos(k)+sqrt((2*C.*cos(k)).^2+kai^2); kai]./sqrt((2*C.*cos(k)+sqrt((2*C.*cos(k)).^2+kai^2)).^2+kai^2);
phi_initial=phi_plus(k0,C(V_before));


%% 
N_t=100;

delta_z=L/N_z;
E_initial=zeros(2,N_z/2);
for n_z=1:N_z
    z=n_z*delta_z;
    if mod(n_z,2)==1
        E_initial(1,floor(n_z/2)+1)=sqrt(sigma_f)/power(pi,1/4).*exp(-(v_g*k0-z).^2./(2*sigma_T^2)).*double(phi_initial(1)); % Eqn.21 
    else
        E_initial(2,floor(n_z/2))=sqrt(sigma_f)/power(pi,1/4).*exp(-(v_g*k0-z).^2./(2*sigma_T^2)).*double(phi_initial(2)); % Eqn.21 
    end
end

s=v_g*delta_t/delta_z/2;
%Lax-Friendrichs mathod
D0=(diag((s+1/2)*ones(N_z/2-1,1),1))+(diag((-s+1/2)*ones(N_z/2-1,1),-1));
E_t=zeros(2,N_z/2,N_in_trip);
E_omiga_A=zeros(N_t,N_omega_m/2);
E_omiga_B=zeros(N_t,N_omega_m/2);
% couple_m=[sqrt(1-gamma_c^2),    -1i*gamma_c; ...
%              -1i*gamma_c,       sqrt(1-gamma_c^2)];
N_z0=floor(N_z/4); % the location in the ring where we sample the field over time
E_t_omiga_single=zeros(2,N_in_trip,N_omega_m/2);
V_t=zeros(N_in_trip*N_t,1);
E_t(:,:,1)=E_initial;
E_omiga=zeros(2,N_t,N_omega_m/2);
E_omiga_2=zeros(N_t,N_omega_m/2);
z_pm=floor(N_z/4);
z_0=floor(N_z/8);
E_t_stor=zeros(2,N_t,N_omega_m);
for n_trip=1:N_t 
    if n_trip<=N_t/2
%             V=V_before;  
         V=V_before+disorder_C_before*(rand-1);
    else
%             V=V_after;
        V=V_after+disorder_C_after*(rand-1);
    end
    for tt=1:N_in_trip
        n_t=(n_trip-1)*N_in_trip+tt;
        t=n_t*delta_t;

        V_t(n_t)=V*cos(Omega*t);
        % for ring A&B z+pm
        if tt>1
            E_t(:,:,tt)=E_t(:,:,tt-1)*D0; %eqn(17)
            E_t(:,1,tt)=(s+1/2)*E_t(:,end,tt-1)+(1/2-s)*E_t(:,1,tt-1); 
            E_t(:,end,tt)=(s+1/2)*E_t(:,end-1,tt-1)+(1/2-s)*E_t(:,end,tt-1);
        end
        E_t(1,z_pm,tt)=exp(1i*V_t(n_t)).*E_t(1,z_pm,tt); %eqn(18) L/2:the location of the phase modulator z_pm
        E_t(2,z_pm,tt)=exp(-1i*V_t(n_t)).*E_t(2,z_pm,tt);
        gamma_c=gamma_c0+(rand-1)*disorder_gamma;
        couple_m=[sqrt(1-gamma_c^2),    -1i*gamma_c; ...
             -1i*gamma_c,       sqrt(1-gamma_c^2)];
        E_t(:,z_0,tt)=(couple_m)*E_t(:,z_0,tt); %eqn(19),eqn(20)
        E_t_stor(:,n_trip,tt)=E_t(:,N_z0,tt);%zeros(2,N_in_trip,N_omiga/2)
        %  Fourier transform the fields in time
        for n_omega=1:N_omega_m
            if mod(n_omega,2)==1
                E_t_omiga_single(1,tt,floor(n_omega/2)+1) = E_t(1,N_z0,tt).*exp(-1i*omega_m(n_omega)*t)*delta_t;   %eqn (22)  
            else
                E_t_omiga_single(2,tt,floor(n_omega/2)) = E_t(2,N_z0,tt).*exp(-1i*omega_m(n_omega)*t)*delta_t;   %eqn (22)  
            end
        end
    end
    E_t(:,:,1)=E_t(:,:,end)*D0; %eqn(17)
    E_t(:,1,1)=(s+1/2)*E_t(:,end,end)+(1/2-s)*E_t(:,2,end);
    E_t(:,end,1)=(s+1/2)*E_t(:,end-1,end)+(1/2-s)*E_t(:,1,end);

     for omega_index=1:N_omega_m/2
        E_omiga_A(n_trip,omega_index)=sum(E_t_omiga_single(1,:,omega_index));
        E_omiga_B(n_trip,omega_index)=sum(E_t_omiga_single(2,:,omega_index));
        E_omiga(:,n_trip,omega_index)=[E_omiga_A(n_trip,omega_index);E_omiga_B(n_trip,omega_index)];
        E_omiga_2(n_trip,omega_index)=E_omiga(:,n_trip,omega_index)'*E_omiga(:,n_trip,omega_index);
    end
end
[~,saperate_point]=find(E_omiga_2(49,:)==max(E_omiga_2(49,:)));

t_ratio=sum(E_omiga_2(end,1:saperate_point))/sum(E_omiga_2(end,:));
end



















