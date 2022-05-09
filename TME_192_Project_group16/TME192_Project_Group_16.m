%% TME192 Active Safety
% Project 2019/20
% 
% Group Number 
% Student 1: Fikri Farhan, Witjaksono; fikrif@student.chalmers.se
% Student 2: Nithin , Revadal; nithinr@student.chalmers.se
% Student 3: Subhash, Ravi Kumar; ravisu@student.chalmers.se
%
% Student template
% Add your code to answer the tasks and questions below
% Please comment your code, this will help us to better understand what you are doing and to speed up the correction
%
%% -------------------------------------------------------------------- %%
clc
clear
close all

% Load provided Dataset

load RadarData.mat

% Create new Dataset to store results in
RadarDataCorrected_group16 =[] ;


% Parameters for identifying braking maneuver
%BrakingThresholds = ;
% for i= 1:1:size(RadarDataCorrected,2)
%         find RadarDataCorrected.VehicleAcceleration(i) < -0.5
%         RadarDataCorrected.VehicleAcceleration(i+1) < -0.5
%         RadarDataCorrected.VehicleAcceleration(i+2)< -0.5
%         RadarDataCorrected.VehicleAcceleration(i+3)< -0.5       
% end        
% Start
BrakingThresholds.start.acceleration = [] ;% Acceleration in m/s^2
BrakingThresholds.start.speed = []; % Speed in m/s
BrakingThresholds.start.yawrate = []; % Yaw Rate in rad/s

% End  
BrakingThresholds.end.acceleration =[] ; % Acceleration in m/s^2
BrakingThresholds.end.speed = []; % Speed in m/s
BrakingThresholds.end.yawrate = []; % Yaw Rate in rad/s

% Should plots be created?
plotcreation = 1; % 0 for no, 1 for yes
%% -------------------------------------------------------------------- %
%% Task 1 - A hypothetical rear-end conflict situation
% Given parameters
initial_velocity = 90*5/18; % Initial velocity [m/s]
reaction_time = 1.5; % Reaction time [s]
a = -5; % Deceleration [m/s^2]

% Calculations
% Sub-task-1
final_velocity = 0; % Final velocity [m/s]
time = (final_velocity-initial_velocity)/a; % Time for collision with brakes applied [s]
minimum_range = (initial_velocity*time)+((1/2)*a*time^2); % Minimum range at which brake needs to be applied [m]

%Sub-task-2
ttc = (minimum_range/initial_velocity);%Time to collision [s]

%Sub-task-3
Distance1 = initial_velocity*reaction_time; %Distance covered before reacting [m]
Distance_ac = Distance1+minimum_range; %Distance to avoid collision[m]

% Since the distance needed to avoid collision is greater than the available distance (90m),
%  its not possible to avoid the collision.


%% -------------------------------------------------------------------- %
% Task 2 - Study critical braking behavior with experimental data
load('RadarData');
RadarDataCorrected_group16 = struct('RadarRange',[],'RadarRangeRate',[],'RadarAccel',[],'VehicleSpeed',[],'VehicleYawRate',[],'TestPerson',[],'RadarTime',[],'VehicleTime',[],'FileName',[]); 
% data cleaning

NaN_radarr = arrayfun(@(x) sum(isnan(x.RadarRange)), RadarData);
total_array_radarr = arrayfun(@(x) size(RadarData(x).RadarRange,1), 1:numel(RadarData));

NaN_radarrangerate = arrayfun(@(x) sum(isnan(x.RadarRangeRate)), RadarData);
total_array_radarrangerate = arrayfun(@(x) size(RadarData(x).RadarRangeRate,1), 1:numel(RadarData));

NaN_radaraccel = arrayfun(@(x) sum(isnan(x.RadarAccel)), RadarData);
total_array_radaraccel = arrayfun(@(x) size(RadarData(x).RadarAccel,1), 1:numel(RadarData));

NaN_speed = arrayfun(@(x) sum(isnan(x.VehicleSpeed)), RadarData);
total_array_speed = arrayfun(@(x) size(RadarData(x).VehicleSpeed,1), 1:numel(RadarData));

NaN_yawrate = arrayfun(@(x) sum(isnan(x.VehicleYawRate)), RadarData);
total_array_yawrate = arrayfun(@(x) size(RadarData(x).VehicleYawRate,1), 1:numel(RadarData));

NaN_radartime = arrayfun(@(x) sum(isnan(x.RadarTime)), RadarData);
total_array_radartime = arrayfun(@(x) size(RadarData(x).RadarTime,1), 1:numel(RadarData));

NaN_vehicletime = arrayfun(@(x) sum(isnan(x.VehicleTime)), RadarData);
total_array_vehicletime = arrayfun(@(x) size(RadarData(x).VehicleTime,1), 1:numel(RadarData));

NaN_testperson = arrayfun(@(x) sum(isnan(x.TestPerson)), RadarData);
total_array_testperson = arrayfun(@(x) size(RadarData(x).TestPerson,1), 1:numel(RadarData));

NaN_filename = arrayfun(@(x) sum(isnan(x.FileName)), RadarData);
total_array_filename = arrayfun(@(x) size(RadarData(x).FileName,1), 1:numel(RadarData));

for i = size(RadarData,2):-1:1
    if  NaN_radarr(i) < 0.8*total_array_radarr(i)       
    
        RadarDataCorrected_group16(i).RadarRange =  RadarData(i).RadarRange;
        
    end



    if NaN_radarrangerate(i) <  0.8*total_array_radarrangerate(i)
       RadarDataCorrected_group16(i).RadarRangeRate = RadarData(i).RadarRangeRate;
             
    end
    if NaN_radaraccel(i) < 0.8*total_array_radaraccel(i)
       RadarDataCorrected_group16(i).RadarAccel = RadarData(i).RadarAccel;
    
    end
    
    if NaN_speed(i) < 0.8*total_array_speed(i)
       RadarDataCorrected_group16(i).VehicleSpeed = RadarData(i).VehicleSpeed;
    end
    
    if  NaN_yawrate(i)< 0.8*total_array_yawrate(i)
        RadarDataCorrected_group16(i).VehicleYawRate = RadarData(i).VehicleYawRate;
    end
    
    if NaN_radartime(i)< 0.8*total_array_radartime(i)
         RadarDataCorrected_group16(i).RadarTime = RadarData(i).RadarTime;
    end
    
    if NaN_vehicletime(i) < 0.8*total_array_vehicletime(i)
        RadarDataCorrected_group16(i).VehicleTime =RadarData(i).VehicleTime;
    end
    
     if NaN_testperson(i)< 0.8*total_array_testperson(i)
         RadarDataCorrected_group16(i).TestPerson = RadarData(i).TestPerson;
    end
    
    if NaN_filename(i) < 0.8*total_array_filename(i)
        RadarDataCorrected_group16(i).FileName =RadarData(i).FileName;
    end
end     
% for i = 1:length(RadarData) 
%  if all(isempty(RadarData(i).RadarRange))
%    index (i)= true;
%  else
%      index(i) = false;
%      
%  end
%  
% end
% RadarDataCorrected.RadarRange=RadarData.RadarRange(index);

RadarRange = arrayfun(@(x) any(ismissing(x.RadarRange)), RadarData);

%Remove the whole run if one of the field is empty
for n = size(RadarDataCorrected_group16,2):-1:1
    if isempty(RadarDataCorrected_group16(n).RadarRange) 
        RadarDataCorrected_group16(n) = [];
    end
end
% The kinematic data with the delay in Radar Time (+200 ms)
for n = size(RadarDataCorrected_group16,2):-1:1
RadarDataCorrected_group16(n).RadarTime(:,1) = RadarDataCorrected_group16(n).RadarTime(:,1)-0.2;
% RadarDataCorrected(n).RadarTime(:,1) = RadarDataCorrected(n).RadarTime(:,1)-RadarDataCorrected(n).RadarTime(1,1);

end
% Normalized Vehicle Time
for n = size(RadarDataCorrected_group16,2):-1:1
RadarDataCorrected_group16(n).VehicleTime(:,1) = RadarDataCorrected_group16(n).VehicleTime(:,1)-RadarDataCorrected_group16(n).VehicleTime(1,1);
end


%%%%Interpolate with fillmissing
 for i = length(RadarDataCorrected_group16):-1:1 
     RadarDataCorrected_group16(i).RadarRange = fillmissing(RadarDataCorrected_group16(i).RadarRange,'linear');
     RadarDataCorrected_group16(i).RadarAccel = fillmissing(RadarDataCorrected_group16(i).RadarAccel,'linear');
     RadarDataCorrected_group16(i).RadarRangeRate = fillmissing(RadarDataCorrected_group16(i).RadarRangeRate,'linear');
     RadarDataCorrected_group16(i).VehicleYawRate = fillmissing(RadarDataCorrected_group16(i).VehicleYawRate,'linear');
     RadarDataCorrected_group16(i).VehicleSpeed = fillmissing(RadarDataCorrected_group16(i).VehicleSpeed,'linear');
 end
 
 % Create New Field Acceleration
for n = size(RadarDataCorrected_group16,2):-1:1
    n1 = gradient(RadarDataCorrected_group16(n).VehicleSpeed);
    n2 = gradient(RadarDataCorrected_group16(n).VehicleTime);
    RadarDataCorrected_group16(n).VehicleAcceleration = n1./n2;
end

 %Find slope that is negative as a start of the braking and when it becomes
 %zero again as the end
 
for i= size(RadarDataCorrected_group16,2):-1:1
   if RadarDataCorrected_group16(i).VehicleAcceleration(:,1)<0
     BrakingThresholds.start.acceleration(:,1) = RadarDataCorrected_group16(i).VehicleSpeed;
   end
end

%% Braking Maneuver
deceleration_threshold= 0.5; 
deceleration_duration = 2;

 %Find the start and end of the braking occurence and clear if it is less
 %than 2 second long
for n = size(RadarDataCorrected_group16,2):-1:1
   
    start{n}= find(RadarDataCorrected_group16(n).VehicleAcceleration(:,1) < -deceleration_threshold );
%     start_sorted{n} = diff([0,diff(start{1,28}==1,0])
    if numel(start{n})<20
        start{n} = [];
    end
end
% %% Plot creation to help you (and us) visualize results
% % Plot speed, acceleration and radar range per test run
% % Add where your identified braking maneuver starts and where it ends (e.g. by adding a line)
%  if plotcreation ~= 0
%  for i= size(RadarDataCorrected,2):-1:1
%  figure(i)
% plot(RadarDataCorrected(i).VehicleTime,RadarDataCorrected(i).VehicleYawRate)
% hold on
% plot(RadarDataCorrected(i).VehicleTime,RadarDataCorrected(i).VehicleSpeed)
%     plot(RadarDataCorrected(i).VehicleTime,RadarDataCorrected(i).VehicleAcceleration)
% %     plot(RadarDataCorrected(i).RadarTime,RadarDataCorrected(i).RadarRange)
%   legend('Yaw Rate','Speed','Acceleration','RadarRange')
%   hold off
%  end
%  
%  end
 % %
% %% 
% % Plot Radar Range raw data and your corrected data set
% for i=size(RadarData,2):-1:1
% figure(i)
% hold on
%  plot(RadarData(i).RadarTime,RadarData(i).RadarRange)
% legend('rawdata', 'Corrected Raw Data')
% hold off
% end
 
%  end
% %%
% % Create Table for Task 2
% % table with data per lap
% Table_Task2= ;
% 
% for i=1:length(RadarData)
%     Table_Task2.Mean_Acc_TR(i) = ;
%     Table_Task2.Max_Acc_TR(i) = ;
%     Table_Task2.Speed_at_onset_TR(i) = ;
%     Table_Task2.RadarRange_at_onset_TR(i) = ;
%     Table_Task2.TTC_TR(i) = ;
%     Table_Task2.Testperson(i) = ;
%     
% end
% 
% % Print mean results
% Table_Task2.Mean_Acc = median(Table_Task2.Mean_Acc_TR)
% Table_Task2.Max_Accn = median(Table_Task2.Max_Acc_TR)
% Table_Task2.Speed_at_onset = median(Table_Task2.Speed_at_onset_TR)
% Table_Task2.TTC = median(Table_Task2.TTC_TR)
% Table_Task2.Number_TR = length(RadarDataCorrected)
% 
% % Save table
% 
% 
% %% -------------------------------------------------------------------- %
% % Task 3 - Driver Behavior Analysis
load RadarDataCorrected

for i = size(RadarDataCorrected,2):-1:1
  
     brake_start_indx(i) = RadarDataCorrected(i).brake_start_indx;
     RadarDataCorrected(i).RadarRange_brake_onset =  RadarDataCorrected(i).RadarRange(brake_start_indx(i));
     RadarDataCorrected(i).Speed_brake_onset = RadarDataCorrected(i).VehicleSpeed(brake_start_indx(i));
     RadarDataCorrected(i).TTC = RadarDataCorrected(i).RadarRange_brake_onset./RadarDataCorrected(i).Speed_brake_onset;
   
     Speed_brake_onset(i) = RadarDataCorrected(i).Speed_brake_onset;
     TTC(i) =    RadarDataCorrected(i).TTC;
     
end
%% PLOT THE SPEED ON BRAKE ONSET VS TTC (x-y axis)
      

%% Plot Whole Participant TTC
      
      figure(1)
      
      X = [ Speed_brake_onset(1,1:12),TTC(1,1:12),Speed_brake_onset(1,13:18),TTC(1,13:18),...
          Speed_brake_onset(1,19:28),TTC(1,19:28),Speed_brake_onset(1,29:38),TTC(1,29:38),...
          Speed_brake_onset(1,39:45),TTC(1,39:45),Speed_brake_onset(1,46:56),TTC(1,46:56),...
         Speed_brake_onset(1,57:67),TTC(1,57:67),Speed_brake_onset(1,68:77),TTC(1,68:77)];
      
      plot(X(1,1:12),X(1,13:24),'go',X(1,25:30),X(1,31:36),'bo'...
           ,X(1,37:46),X(1,47:56),'ro',X(1,57:66),X(1,67:76),'co'...
           ,X(1,77:83),X(1,84:90),'mo',X(1,91:101),X(1,102:112),'yo'...
           ,X(1,113:123),X(1,124:134),'ko',X(1,135:144),X(1,145:154),'rx');
       
       
       title('All Participant Scatter Plot')
       xlabel('Speed at Brake Start [m/s]'); ylabel('TTC at Brake Start [s]');
     legend( 'TTC_{Driver 1}', 'TTC_{Driver 2}','TTC_{Driver 4}','TTC_{Driver 5}'...
            ,'TTC_{Driver 6}', 'TTC_{Driver 7}','TTC_{Driver 8}','TTC_{Driver 9}','Location','bestoutside')
            
        %% Spectral Clustering Analysis
          Y1 = [Speed_brake_onset(1,1:12); TTC(1,1:12)];
          Y1 = Y1.';
          Y2 = [Speed_brake_onset(1,13:18);TTC(1,13:18)];
          Y2 = Y2.';
          Y4 = [Speed_brake_onset(1,19:28);TTC(1,19:28)];
          Y4 = Y4.';
          Y5 = [Speed_brake_onset(1,29:38);TTC(1,29:38)];
          Y5 = Y5.';
          Y6 = [Speed_brake_onset(1,39:45);TTC(1,39:45)];
          Y6 = Y6.';
          Y7 = [Speed_brake_onset(1,46:56);TTC(1,46:56)];
          Y7 = Y7.';
          Y8 = [Speed_brake_onset(1,57:67);TTC(1,57:67)];
          Y8 = Y8.';
          Y9 = [Speed_brake_onset(1,68:77);TTC(1,68:77)];
          Y9 = Y9.';
          
          figure(2)
        Y = [Y1;Y2;Y4;Y5;Y6;Y7;Y8;Y9];
        
        
%         
%         opts = statset('Display','final');
%         [idxx,C] = kmeans(Y,3,'Distance','cityblock',...
%                   'Replicates',5,'Options',opts);

%         plot(Y(idxx==1,1),Y(idxx==1,2),'r.','MarkerSize',12)
%         hold on
%         plot(Y(idxx==2,1),Y(idxx==2,2),'b.','MarkerSize',12)
%         plot(C(:,1),C(:,2),'kx',...
%         'MarkerSize',15,'LineWidth',3) 
%         legend('Cluster 1','Cluster 2','Centroids',...
%          'Location','NW')
%         title 'Cluster Assignments and Centroids'
%         xlabel('Speed [m/s]'); ylabel('Speed [m/s]');
%         hold off
        
       idxx = spectralcluster(Y,3,'SimilarityGraph','epsilon','Radius',2);
       gscatter(Y(:,1),Y(:,2),idxx);
        xlabel('TTC at Brake Start [s]'); ylabel('Speed at Brake Start [m/s]');
        legend( '1st cluster', '2nd Cluster','3rd Cluster','Location','bestoutside')
      
% 

%% Plot Individual Participant TTC with Clustering
%Plot the first 4 participant
figure(3)      
t = tiledlayout(2,2);
%       
%       % First Participant 
       nexttile
        idx = spectralcluster(Y1,3,'SimilarityGraph','epsilon','Radius',2);
       gscatter(Y1(:,1),Y1(:,2),idx);
        xlabel('Speed [m/s]'); ylabel('TTC [s]');
        legend( '1st cluster', '2nd Cluster','3rd Cluster','Location','bestoutside')
         title('Participant 1')
       
        % Second Participant
        nexttile
        idx2 = spectralcluster(Y2,3,'SimilarityGraph','epsilon','Radius',2);
       gscatter(Y2(:,1),Y2(:,2),idx2);
        xlabel('Speed [m/s]'); ylabel('TTC [s]');
        legend( '1st cluster', '2nd Cluster','3rd Cluster','Location','bestoutside')
         title('Participant 2')
       
       % Third Participant
        nexttile
        idx4 = spectralcluster(Y4,3,'SimilarityGraph','epsilon','Radius',2);
       gscatter(Y4(:,1),Y4(:,2),idx4);
       xlabel('Speed [m/s]'); ylabel('TTC [s]');
        legend( '1st cluster', '2nd Cluster','3rd Cluster','Location','bestoutside')
         title('Participant 4')
       
        % Fourth Participant
        nexttile
        idx5 = spectralcluster(Y5,3,'SimilarityGraph','epsilon','Radius',2);
       gscatter(Y5(:,1),Y5(:,2),idx5);
        xlabel('Speed [m/s]'); ylabel('TTC [s]');
        legend( '1st cluster', '2nd Cluster','3rd Cluster','Location','bestoutside')
         title('Participant 5')
       
       
       t.Padding = 'none';
       t.TileSpacing = 'none';
       hold on
       
        %Plot the last 4 participant
        figure(4)
      t2 = tiledlayout(2,2);
      % Sixth Participant 
       nexttile
        idx6 = spectralcluster(Y6,3,'SimilarityGraph','epsilon','Radius',2);
       gscatter(Y6(:,1),Y6(:,2),idx6);
      xlabel('Speed [m/s]'); ylabel('TTC [s]');
        legend( '1st cluster', '2nd Cluster','3rd Cluster','Location','bestoutside')
         title('Participant 6')
       
        % Seventh Participant
      nexttile
        idx7 = spectralcluster(Y7,3,'SimilarityGraph','epsilon','Radius',2);
       gscatter(Y7(:,1),Y7(:,2),idx7);
        xlabel('Speed [m/s]'); ylabel('TTC [s]');
        legend( '1st cluster', '2nd Cluster','3rd Cluster','Location','bestoutside')
         title('Participant 7')
       
       % Eighth Participant
      nexttile
        idx8 = spectralcluster(Y8,3,'SimilarityGraph','epsilon','Radius',2);
       gscatter(Y8(:,1),Y8(:,2),idx8);
        xlabel('Speed [m/s]'); ylabel('TTC [s]');
        legend( '1st cluster', '2nd Cluster','3rd Cluster','Location','bestoutside')
         title('Participant 8')
       
        % Ninth Participant
       nexttile
        idx9 = spectralcluster(Y9,3,'SimilarityGraph','epsilon','Radius',2);
       gscatter(Y9(:,1),Y9(:,2),idx9);
        xlabel('Speed [m/s]'); ylabel('TTC [s]');
        legend( '1st cluster', '2nd Cluster','3rd Cluster','Location','bestoutside')
         title('Participant 9')
       
       t2.Padding = 'none';
       t2.TileSpacing = 'none';
       
       %% Plot the TTC distribution
       
       % Calculate the 5%-ile and 95%-ile
       NinetyFive_Percentile = prctile(TTC,95)
       Five_Percentile = prctile(TTC,5)
       
         TTC_r = reshape(TTC,[],1);
        pdTTC = fitdist(TTC_r,'Normal')
        x3 = [-3:.1:3];
        y3 = normpdf(x3,0,1);
        
        mu = 2.02632;
        sigma = 0.650244;
%       
       figure (5)
       hold on
       plot(x3,y3)
       title('TTC Distribution Plot');
       ylabel('Probability Density');
       z3 = xline(0,'--','Mean');
       z3.LabelVerticalAlignment = 'middle';
       z3.LabelHorizontalAlignment = 'center';
       
       sigma_95percentile = (NinetyFive_Percentile-mu)/sigma;
       v3 = xline(sigma_95percentile,'--','95%-ile');
       v3.LabelVerticalAlignment = 'top';
       v3.LabelHorizontalAlignment = 'center';
       
       sigma_5percentile = (Five_Percentile-mu)/sigma;
       r3 = xline(sigma_5percentile,'--','5%-ile');
       r3.LabelVerticalAlignment = 'top';
       r3.LabelHorizontalAlignment = 'center';
       
       p = normspec([sigma_5percentile,sigma_95percentile]);
      
       hold off
       

       %% Calculate the Mean Acceleration and Maximum Acceleration
       for i = size(RadarDataCorrected,2):-1:1
               brake_end_indx(i) = RadarDataCorrected(i).brake_end_indx;
       end
        for i = size(RadarDataCorrected,2):-1:1
                   
               RadarDataCorrected(i).SelectedAcceleration(:,1) = RadarDataCorrected(i).VehicleAccel((brake_start_indx(i):brake_end_indx(i)),1);
%                
               RadarDataCorrected(i).MaximumAcceleration = min(RadarDataCorrected(i).SelectedAcceleration(:,1));
               MaximumAcceleration(i) = RadarDataCorrected(i).MaximumAcceleration;
               RadarDataCorrected(i).MeanAcceleration = mean((RadarDataCorrected(i).SelectedAcceleration(:,1)));
               MeanAcceleration(i) = RadarDataCorrected(i).MeanAcceleration;
                
        end
       
       %% Plot the Maximum Acceleration and Mean Acceleration Normal
          %Distribution
        
          %% Plot Maximum Acceleration Distribution and Analyze
        MaximumAcceleration_r = reshape(MaximumAcceleration,[],1);
        pd1 = fitdist(MaximumAcceleration_r,'Normal')
        x1 = [-3:.1:3];
        y1 = normpdf(x1,0,1);
        
        NinetyFive_Percentile_max = prctile(MaximumAcceleration,95)
        Five_Percentile_max = prctile(MaximumAcceleration,5)
        
         mu_max = -4.97971;
        sigma_max = 2.66164  ;
        
        figure (7)
         hold on
       plot(x1,y1)
       title('Maximum Deceleration Distribution Plot');
        ylabel('Probability Density');
        
        sigma_95percentile_max = (NinetyFive_Percentile_max-mu_max)/sigma_max;
       a3 = xline(sigma_95percentile_max,'--','95%-ile');
       a3.LabelVerticalAlignment = 'top';
       a3.LabelHorizontalAlignment = 'center';
       
       sigma_5percentile_max = (Five_Percentile_max-mu_max)/sigma_max;
       b3 = xline(sigma_5percentile_max,'--','5%-ile');
       b3.LabelVerticalAlignment = 'top';
       b3.LabelHorizontalAlignment = 'center';
       
       p = normspec([sigma_5percentile_max,sigma_95percentile_max]);
     
       hold off
      
        %% Plot Mean Acceleration Distribution and Analyze
        
        MeanAcceleration_r = reshape(MeanAcceleration,[],1);
        pd2 = fitdist(MeanAcceleration_r,'Normal')
        x2 = [-3:.1:3];
        y2 = normpdf(x2,0,1);
        
         mu_mean =  -3.03818  ;
        sigma_mean = 1.42419 ;
        
%         cdf = makedist('Normal')
%         z2 = cdf(x2
%         cum_df = makedist('Normal','mu',mu_mean,'sigma',sigma_mean);
%         p = cdf(cum_df,x2);
%         

        NinetyFive_Percentile_mean = prctile(MeanAcceleration,95)
        Five_Percentile_mean = prctile(MeanAcceleration,5)
        
        
        
       figure (9)
        hold on
        plot(x2,y2)
       title('Mean Deceleration Distribution Plot');
         ylabel('Probability Density');
         
        
         
       sigma_95percentile_mean = (NinetyFive_Percentile_mean-mu_mean)/sigma_mean;
       c3 = xline(sigma_95percentile_mean,'--','95%-ile');
       c3.LabelVerticalAlignment = 'top';
       c3.LabelHorizontalAlignment = 'center';
       
       sigma_5percentile_mean = (Five_Percentile_mean-mu_mean)/sigma_mean;
       d3 = xline(sigma_5percentile_mean,'--','5%-ile');
       d3.LabelVerticalAlignment = 'top';
       d3.LabelHorizontalAlignment = 'center';
       
       p = normspec([sigma_5percentile_mean,sigma_95percentile_mean]);

       hold off
       figure(11)
       histfit(TTC)
       xlabel('TTC [s]')
       ylabel('Probability (%)')
       hold off 
       
       figure(12)
       histfit(MaximumAcceleration,6)
       xlabel('Maximum Deceleration [m/s^2]')
       ylabel('Probability (%)')
       hold off
       figure(13)
       histfit(MeanAcceleration,8)
       xlabel('Mean Deceleration [m/s^2]')
       ylabel('Probability (%)')
       
% 
% % Plot TTC at brake onset vs speed at brake onset?  
% figure()
% 
% 
% % Do the different drivers have both high and low TTC values or are there clusters in the data? 
% figure()
% 
% 
% % What does the TTC distribution look like? What are the 5th and 95th percentile?  
% 
% 
% % Plot distributions for mean accelerations and maximum accelerations (results from Table 2). What
% % information about driver behavior can you obtain from these?
% figure ()
% 
% %% -------------------------------------------------------------------- %
% % Task 4 - Active Safety System Design 
% 
% % FCW system (conservative)
% 
% 
% % FCW system (aggressive)
% 
% 
% % AEB system
% 
% 
% 
