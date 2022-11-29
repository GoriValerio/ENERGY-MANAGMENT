clc;
clearvars;
close all;
clear all;
% set(gca,'TickLabelInterpreter','latex');

%% load data
load('Assignment1_data.mat');

%% variabels

omega = ice.Points(:,1);
tau = ice.Points(:,2);

omega_em = em.Points(:,1);
tau_em = em.Points(:,2);

h0 = 42500; %[j/g]   values [g/s]

%% part 1 ICE
omega_min = min(omega);
omega_max = max(omega);
tau_min = min(tau);
tau_max = max(tau);

[omega2d, tau2d] = ndgrid(omega_min:1:omega_max, tau_min:1:tau_max);
% z1 = griddata(omega, tau, ice.Values, omega2d, tau2d, 'linear');
mass_flow = ice(omega2d, tau2d);

%% part 1 EM
omega_min_em = min(omega_em);
omega_max_em = max(omega_em);
tau_min_em = min(tau_em);
tau_max_em = max(tau_em);

[omega2d_em, tau2d_em] = ndgrid(omega_min_em:1:omega_max_em, tau_min_em:1:tau_max_em);
z1_em = griddata(omega_em,tau_em, em.Values, omega2d_em, tau2d_em, 'linear');

%% efficiency ICE

P_in_ice = h0*ice(omega2d, tau2d); % J/s 
P_out_ice = tau2d.*omega2d; % [w]

eta = P_out_ice./P_in_ice;

dimension = size(eta);

for i = 1:dimension(1)
    for k = 1:dimension(2)
        if(eta(i, k) < 0 || eta(i, k)>1)
            eta(i,k) = NaN;
        end
    end
end


%% 

%% Willans lines (ICE)
Pp = omega.*tau;

for i = 1 : 1: length(Pp)
    if(Pp(i) < 0)
        Pp(i) = NaN;
    end
end


%% eta em eff

 P_in_em = tau2d_em.*omega2d_em;
P_out_em = em(omega2d_em, tau2d_em);

eta_em = P_out_em./P_in_em;

dimension_em = size(eta_em);

for i = 1:dimension_em(1)
    for k = 1:dimension_em(2)
        if(eta_em(i, k) < 0 || eta_em(i, k)>1)
            eta_em(i,k) = NaN;
        end
    end
end

%% mirror em eff

P_in_em = tau2d_em.*omega2d_em;
P_out_em = em(omega2d_em, tau2d_em);

eta_em_mirror = P_in_em./P_out_em;

dimension_em = size(eta_em);

for i = 1:dimension_em(1)
    for k = 1:dimension_em(2)
        if(eta_em_mirror(i, k) < 0 || eta_em_mirror(i, k)>1)
            eta_em_mirror(i,k) = NaN;
        end
    end
end



%% OOL
% P_out_ice_max = max(P_out_ice);
% P_out_ice_min = min(P_out_ice);
% [min_fuel, index_min_ice_omega] = min(mass_flow); 
% Power_line = linspace(P_out_ice_min, P_out_ice_max);

P_out_ice_lin = linspace(0, 136e3, 1001);
for i= 1:length(P_out_ice_lin)
    tau_ice_lin_2d(i,:) = linspace(0, 200, 1001);
    P_out_ice_lin_2d(i,1:1001) = P_out_ice_lin(i);
    omega_ool(i,:) = P_out_ice_lin_2d(i,:)./tau_ice_lin_2d(i,:); 
    Fuel_ool(i,:) = ice(omega_ool(i,:), tau_ice_lin_2d(i,:));
    [Fuel_min(i,1), Fuel_min_index(i,1)] = min(Fuel_ool(i,:));
    omega_ool_min(i,1) = omega_ool(i, Fuel_min_index(i));
    tau_ice_min(i,1) = tau_ice_lin_2d(i, Fuel_min_index(i));
end

%%
tau2d_pos = tau2d(tau2d > 0);
%% plotting
figure(1)
contourf(omega2d, tau2d, mass_flow, 30);
colormap(hot);
hold on;
plot(omega_ool_min, tau_ice_min);
ylabel('Torque $\tau$ [Nm]','interpreter','latex');
xlabel('Angular velocity $\omega$ [rad/s]', 'interpreter','latex');
title('ICE', 'interpreter','latex');
hold off;

%%
figure(2)
contourf(omega2d_em, tau2d_em, z1_em, 30);
hold on;
colormap(hot);
ylabel('Torque $\tau$ [Nm]','interpreter','latex');
xlabel('Angular velocity $\omega$ [rad/s]', 'interpreter','latex');
title('EM', 'interpreter','latex');
hold off;
%% eta

figure(3)
contourf(omega2d, tau2d, eta, 30);
hold on;
ylabel('Torque $\tau$ [Nm]','interpreter','latex');
xlabel('Angular velocity $\omega$ [rad/s]', 'interpreter','latex');
title('ICE', 'interpreter','latex');

%% willans lines plot

figure(4)
for i = 1:1:15
plot(Pp(i*25-24:i*25), ice.Values(i*25-24:i*25));
hold on
end
ylabel('Fuel consumption [g/s]','interpreter','latex');
xlabel('Power [w]', 'interpreter','latex');
title('ICE', 'interpreter','latex');
plot(P_out_ice_lin_2d, Fuel_min, '-or')
hold off;

%%
figure(5)
contourf(omega2d_em, tau2d_em, eta_em, 30);
hold on;
contourf(omega2d_em, tau2d_em, eta_em_mirror, 30);
% colormap(hot);
ylabel('Torque $\tau$ [Nm]','interpreter','latex');
xlabel('Angular velocity $\omega$ [rad/s]', 'interpreter','latex');
title('EM', 'interpreter','latex');

%% 1.2

t = linspace(0,1800, 1801); 

WLTP_v_kph = WLTP.Values;
WLTP_v = WLTP_v_kph./3.6;
%%
figure(6);
plot(t, WLTP_v_kph);
grid on;
ylabel('Velocity $v$ [km/h]','interpreter','latex');
xlabel('time $t$ [s]', 'interpreter','latex');
title('WLTP cycle', 'interpreter','latex');

%% Pw
% $m_vh ̇v = Ftract − m_vh g sin α − croll m_vh g cos α ·signv −cdragv2$;
% α = 0
% 

m_vh = 1500; % kg
g = 9.81; 
Cdrag = 0.455; % [kg/m] 
Croll = 0.001;

%%
dV = gradient(WLTP_v);
dt = 1; % s
acc = dV/dt;
%%
P_requested = (m_vh.*acc + Croll.*m_vh.*g.*sign(WLTP_v) + Cdrag.*WLTP_v.^2).*WLTP_v; % [N]
F_requested = (m_vh.*acc + Croll.*m_vh.*g.*sign(WLTP_v) + Cdrag.*WLTP_v.^2);

v_from_F = 0;

%%
figure(7);
plot(t,P_requested);
% plot(t, acc);

%% PI controller
Kp = 1;
Ki = 1;
% c_pi = pid(Kp, Ki);

% v_driver = 0;

% for i:length(WLTP_v)
%     v_driver(i) = WLTP_v(i);
%     acc_driver(i) = ;
% %     F_requested_driver = (m_vh.*acc_driver + Croll.*m_vh.*g.*sign(WLTP_v) + Cdrag.*WLTP_v.^2);
% end


%????

%% energy of cycle, cumulative positive propulsion power

P_requested_pos = P_requested( P_requested>=0 );
E_req = mean(P_requested_pos)*max(t); %[w*s = J]
E_req_kwh = E_req/3.6e6;

kwh_p_100km = E_req_kwh/((mean(WLTP_v_kph)*max(t)/3600))*100 %[kwh/100km]

% check v2 from matlab? -  GIT works woho

