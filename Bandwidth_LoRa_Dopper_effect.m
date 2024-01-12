close all, clear all
%% Paper title: Understanding the Limits of LoRa Direct-to-Satellite: The Doppler Perspectives
%% Authors: M. Asad Ullah, G. Pasolini, K. Mikhaylov and H. Alves
%% Journal: IEEE Open Journal of the Communications Society, vol. 5, pp. 51-63, 2024

%% Key variable and parameters

% BW = Bandwidth
% Application_payload = payload
% LDRO = Low Data Rate Optimization
% SF = Spreading factor
% E  = Elevation angle
% Rp = Average inter arrival time
% L_static = loss due to static Doppler shift
% L_dynamic = loss due to dynamic Doppler shift

BW=31.25e3;
Fc = 433e6;  % carrier frequency in Hz
Application_payload = 55;
LDRO = true;
Rp = 5000;
LoRa_ToA

%% Norbay's pass duration

startTime = datetime(2022,9,29,20,54,02,00);
stopTime = datetime(2022,9,29,21,05,50,00); 
sampleTime = 1;  % simulation step 1s

%% Create a satellite scenarios
sc = satelliteScenario(startTime,stopTime,sampleTime);
sat = satellite(sc,"Norby.tle");     % Read TLE

%% Create a ground station scenarios
lat = [55.01880];
lon = [82.933952];
gs = groundStation(sc,lat, lon,'MinElevationAngle',1);

%% Access analysis
ac = access(sat,gs);
intvls = accessIntervals(ac);
acStatus = accessStatus(ac);

%% Orbital analysis
numHours = stopTime - startTime;
numSeconds = seconds(numHours);
[az,el,r] = deal(zeros(numSeconds+1,1));
[satV,dir] = deal(zeros(3,numSeconds+1));
[latency,dopV] = deal(NaN(numSeconds+1,1));
c = physconst("Lightspeed");

for iMiliSec = 0:1:numSeconds
  time(iMiliSec+1) = startTime + seconds(iMiliSec);
  idx = iMiliSec+1;
  pos(idx,:) = states(sat,time(iMiliSec+1),"CoordinateFrame","geographic");
  
  if acStatus(idx)
    % Calculate azimuth, elevation, and range from the satellite to the
    % ground station
    [az(idx),el(idx),r(idx)] = aer(sat,gs,time(iMiliSec+1));

    % Calculate latency
    latency(idx) = r(idx) / c;
    % Calculate slant range
    Slant_range(idx) = r(idx);
    
    % Calculate satellite velocity in North/East/Down (NED) frame of
    % satellite. Physically, this is the ECEF velocity, represented in NED
    % frame. Therefore, the relative velocity with respect to the ground
    % station is also the same.

    [~,satV(:,idx)] = states(sat,time(iMiliSec+1),"CoordinateFrame","geographic");

    % Calculate the direction of gs with respect to sat in sat NED frame
    dir(:,idx) = [cosd(el(idx))*cosd(az(idx)); ...
      cosd(el(idx))*sind(az(idx)); ...
      -sind(el(idx))];

    % Calculate the velocity along the line between the gs and the sat.
    % This velocity determines the doppler frequency.
    dopV(idx) = dot(satV(:,idx),dir(:,idx));
  end
end

%% Calculations
Index = find(Slant_range~=0);
FD = ((((c+dopV(Index)) ./ (c)) * Fc) - Fc);  % doppler shift in Hz

Length_Doppler_Array = length(FD);
Doppler_Array_Indexes = 1:1:Length_Doppler_Array;
Simulation_time_samples = 1:1/1e3:table2array(intvls(1,6));

%% Convert Slant range to Elevation angle
E = height2el(min(Slant_range)-1,1,Slant_range);
digits(64)
%% Interpolation to create more samples
Interpolated_FD = interp1(Doppler_Array_Indexes,FD,Simulation_time_samples);
Interpolated_Dynamic_Doppler = interp1(Doppler_Array_Indexes(1:1:end-1),diff(FD),Simulation_time_samples(1:1:end-1));
Interpolated_time = interp1(Doppler_Array_Indexes,time,Simulation_time_samples);
Interpolated_Slant_distance = interp1(Doppler_Array_Indexes,Slant_range,Simulation_time_samples);
Interpolated_el = interp1(Doppler_Array_Indexes,E,Simulation_time_samples);

Index_simulated_packets = 1:Rp:length(Interpolated_Dynamic_Doppler)-ToA(6);
%% LoRa tolerance to static Doppler shift
F_static = (25/100).*BW;
%% packet losses due to static doppler shift
L_static = [];
for Doppler_loop=1:Rp:length(Interpolated_FD)-ToA(6)
    if abs(Interpolated_FD(Doppler_loop))>=F_static
        L_static(Doppler_loop) = 1;
    end
end

%% experienced dynamic Doppler shift
for lop=1:1:length(Index_simulated_packets)-1
        Delta_FE_SF7(lop,:) = Interpolated_FD(Index_simulated_packets(lop))-Interpolated_FD(Index_simulated_packets(lop)+ToA(1,:)); % SF7
        Delta_FE_SF10(lop,:) = Interpolated_FD(Index_simulated_packets(lop))-Interpolated_FD(Index_simulated_packets(lop)+ToA(4,:)); %SF10
        Delta_FE_SF12(lop,:) = Interpolated_FD(Index_simulated_packets(lop))-Interpolated_FD(Index_simulated_packets(lop)+ToA(6,:));  %SF12
end

%% LoRa tolerance to dynamic Doppler shift
if(LDRO==true)
        F_Dynamic_SF7=(16*BW)/(3*2^SF(1));
        F_Dynamic_SF10=(16*BW)/(3*2^SF(4));
        F_Dynamic_SF12=(16*BW)/(3*2^SF(6));
else
        F_Dynamic_SF7=(BW)/(3*2^SF(1));
        F_Dynamic_SF10=(BW)/(3*2^SF(4));
        F_Dynamic_SF12=(16*BW)/(3*2^SF(6));
end
%PDR

Simulation_time_t=-find(Interpolated_el==max(Interpolated_el))+1:1:length(Interpolated_el)-find(Interpolated_el==max(Interpolated_el))-ToA(6);

f1=figure(11)
set(gcf, 'Units', 'centimeters'); % set units 
%set(gcf,'TickLabelInterpreter', 'latex')
LineWidth = 2;
axesFontSize = 18;
legendFontSize = 14;
% setting size & position
afFigurePosition = [2 2 19 19]; % [pos_x pos_y width_x width_y] 
set(gcf, 'Position', afFigurePosition,'PaperSize',[19 11],'PaperPositionMode','auto'); % [left bottom width height], setting printing properties 

subplot(311)
scatter(Simulation_time_t(1:5e3:end)/1e3, Interpolated_el(Index_simulated_packets),'b','filled');
L_static_ones = find(L_static);
hold on
scatter(Simulation_time_t(find(L_static))/1e3, Interpolated_el(find(L_static)),'red','filled');

L_dynamic_index = Index_simulated_packets((Delta_FE_SF7>= F_Dynamic_SF7));
L_dynamic_angle = Interpolated_el(Index_simulated_packets((Delta_FE_SF7>= F_Dynamic_SF7)));

if(nnz(Delta_FE_SF7>= F_Dynamic_SF7)>0)
    hold on
    scatter(Simulation_time_t(Index_simulated_packets((Delta_FE_SF7>= F_Dynamic_SF7)))/1e3,Interpolated_el(Index_simulated_packets((Delta_FE_SF7>= F_Dynamic_SF7))),'ko','filled')
    L_joint_index= ismember(L_dynamic_index,L_static_ones);
    if(nnz(L_joint_index)>0)
        scatter(Simulation_time_t(L_dynamic_index(L_joint_index))/1e3,L_dynamic_angle(L_joint_index),[],[0.9290 0.6940 0.1250],'filled');
    end
    lost_packets_SF7=unique([find(L_static) Index_simulated_packets((Delta_FE_SF7>= F_Dynamic_SF7)) L_joint_index]);
else
lost_packets_SF7=unique([find(L_static) Index_simulated_packets((Delta_FE_SF7>= F_Dynamic_SF7))]);

end
% 
grid on
text(-350,70,'SF=7','Color','black','FontSize',12)
text(-350,90,'BW=31.25 kHz','Color','black','FontSize',12)
axis([Simulation_time_t(1)/1e3 Simulation_time_t(end)/1e3 0 100])
%  
subplot(312)
scatter(Simulation_time_t(1:5e3:end)/1e3, Interpolated_el(Index_simulated_packets),'b','filled');
hold on
scatter(Simulation_time_t(find(L_static))/1e3, Interpolated_el(find(L_static)),'red','filled');
L_dynamic_index = Index_simulated_packets((Delta_FE_SF10>= F_Dynamic_SF10));
L_dynamic_angle = Interpolated_el(Index_simulated_packets((Delta_FE_SF10>= F_Dynamic_SF10)));

if(nnz(Delta_FE_SF10>= F_Dynamic_SF10)>0)
    hold on
    scatter(Simulation_time_t(Index_simulated_packets((Delta_FE_SF10>= F_Dynamic_SF10)))/1e3,Interpolated_el(Index_simulated_packets((Delta_FE_SF10>= F_Dynamic_SF10))),'ko','filled')
    L_joint_index= ismember(L_dynamic_index,L_static_ones);
    if(nnz(L_joint_index)>0)
        scatter(Simulation_time_t(L_dynamic_index(L_joint_index))/1e3,L_dynamic_angle(L_joint_index),[],[0.9290 0.6940 0.1250],'filled')
    end
    
end
text(-350,70,'SF=10','Color','black','FontSize',12)
text(-350,90,'BW=31.25 kHz','Color','black','FontSize',12)
ylabel('Elevation Angle','Color','k','Interpreter','Latex','FontSize', axesFontSize)
axis([Simulation_time_t(1)/1e3 Simulation_time_t(end)/1e3 0 100])
grid on
subplot(313)
% 
scatter(Simulation_time_t(1:5e3:end)/1e3, Interpolated_el(Index_simulated_packets),'b','filled');
hold on
scatter(Simulation_time_t(find(L_static))/1e3, Interpolated_el(find(L_static)),'red','filled');
L_dynamic_index = Index_simulated_packets((Delta_FE_SF12>= F_Dynamic_SF12));
L_dynamic_angle = Interpolated_el(Index_simulated_packets((Delta_FE_SF12>= F_Dynamic_SF12)));

if(nnz(Delta_FE_SF12>= F_Dynamic_SF12)>0)
    hold on
    scatter(Simulation_time_t(Index_simulated_packets((Delta_FE_SF12>= F_Dynamic_SF12)))/1e3,Interpolated_el(Index_simulated_packets((Delta_FE_SF12>= F_Dynamic_SF12))),'ko','filled')
    L_joint_index= ismember(L_dynamic_index,L_static_ones);
    if(nnz(L_joint_index)>0)
        scatter(Simulation_time_t(L_dynamic_index(L_joint_index))/1e3,L_dynamic_angle(L_joint_index),[],[0.9290 0.6940 0.1250],'filled');
    else

    end
    
end
grid on
 text(-350,70,'SF=12','Color','black','FontSize',12)
text(-350,90,'BW=31.25 kHz','Color','black','FontSize',12)
 axis([Simulation_time_t(1)/1e3 Simulation_time_t(end)/1e3 0 100])
 xlabel('Simulation time [s]','Interpreter','Latex','FontSize', axesFontSize)
 legend('successful','failed: static','failed: dynamic','failed: joint','Location','southeast','NumColumns',4, 'Interpreter', 'Latex','FontSize',12)