% Exercise: Field data acquisition and analysis
% Author:   Alberto Morando - morando@chalmers.se
%           Christian-Nils Boda - boda@chalmers.se
% Date: Oct 2018
% Institution: Chalmers, course: TME 192 Active Safety

close all
clear all
clc

%% WARNING
% Before using this script, convert the BAG file into a CSV file with the
% function `bag2csv`. For example `bag2csv('2018-10-24-11-11-25_0.bag')`

%% DATA LOADING
% `import_lidar_data` Fill a structure containing the recorded data. The structure is as follow:
% +-------------+-------+---------------------------------------------------------------------+
% | Timestamp_s | [s]   | n x 1 vector, where n is the number of timestamps (or frames)       |
% | Angle_rad   | [rad] | 1 x n vector, where n is the number of scanned angles.              |
% |             |       | The range of angles is [-1.658063, 1.658063] radians (-95, +95 deg) |
% | Distance_m  | [m]   | m x n matrix, where m is the number of timestamps (or frames),      |
% |             |       | and n is the number the scanned angle.                              |
% |             |       | Each value represents the distance at a given time stamp and angle  |
% +-------------+-------+---------------------------------------------------------------------+

filename_csv = '2019-10-17-09-48-14_0.csv' % <--- !!! Change this one with your own file !!!

data = import_lidar_data( filename_csv );

% Sampling settings
sample_period_s = median(diff(data.Timestamp_s));

%% SELECT THE EVENT OF INTEREST
% Find the relevant events for your dataset using the fuction playLidar().
% With `playLidar` determine what is the optimal FOV to detect the
% obstacles, and write down the start and stop time index (frame) for each event into `index_of_events`.
% The matrix `index_of_events` should follow the following template:
%
% +-------+------+
% | start | stop |
% +-------+------+
% |     1 |  123 |
% |   156 |  236 |
% |   ... |  ... |
% +-------+------+

% Uncomment the line to use `playLidar`
% playLidar(fullfile('C:\Users\Fikri Farhan\Desktop\Active Safety Project\Exercise 6', filename_csv))

% The following event indexes and field of view contraint are valid only for
% the sample dataset: 2018-10-24-11-11-25_0.csv included with the exercise.
index_of_events = [
    924, 966;
    2030, 2084;
    3067, 3128;
    3129, 3183;
    3184, 3275;
    4320, 4387;
    4440, 4504];

FOV_limits_deg = 65;

FOV_limits_rad =  deg2rad([-FOV_limits_deg/2, FOV_limits_deg/2]);

% Extract the index for the columns within the selected field of view (FOV)
FOV_span_idx = [find(data.Angle_rad >= FOV_limits_rad(1), 1, 'first') : find(data.Angle_rad <= FOV_limits_rad(end), 1, 'last')];

%% PLOT A SNAPSHOT OF THE LIDAR MEASUREMENT
% Select a _single_ frame from one of your events and plot the LIDAR
% measurement
frameNumber = 924;  % <--- !!! Change this with respect to your data!!!

figure('name', 'Example of FOV')
hold on

% PLot the LIDAR origin
plot(0, 0, 'r^', 'markersize', 10, 'MarkerFaceColor', 'r');

% Plot the distance measurements in the full field of view
Angle = data.Angle_rad(1,:);
Dis = data.Distance_m(frameNumber,:);
[X, Y] = pol2cart(Angle,Dis);
plot(Y, X, '-k');

% Extract the data in the chosen FOV
[X_fov,Y_fov] = pol2cart(data.Angle_rad(1,FOV_span_idx), data.Distance_m(frameNumber,FOV_span_idx));
plot(Y_fov, X_fov, '-k', 'linewidth', 2);

% M the closest point in the field of view
n=6;
[val_1,idx]=min(abs(X_fov-n));
x_closest = X_fov(idx);

n_2=-0.1359;
[val_2,idx_2]= min(abs(Y_fov-n_2));
y_closest = Y_fov(idx_2);

plot(y_closest, x_closest, 'xb', 'linewidth', 1.2, 'markersize', 8);

% Plot the boundaries of the selected field of view
xlim_1 = linspace(0,-9.142,10);
xlim_2 = linspace(0,7.5,10);
ylim_1 = linspace(0,14,10);
ylim_2 = linspace(0,12,10);
plot(xlim_1,ylim_1 , ':r', 'linewidth', 1.5) % left limit
plot(xlim_2, ylim_2, ':r', 'linewidth', 1.5) % right limit
  
% For a better visualization you should limit the axis in the figure
set(gca, 'xlim', ([-10 10]), 'ylim', ([0 16]))
axis equal

xlabel('Distance [m]'); 
ylabel('Distance [m]');

legend('LIDAR', 'Data_{all}', 'Data_{FOV}', 'Closest point', 'Boundaries of ','the field of view')

%% PLOT DISTANCE OF TARGET ACROSS TIME

% Fetch the time for each Event 
Time_1= data.Timestamp_s((index_of_events(1,1):index_of_events(1,2)),1)-data.Timestamp_s(index_of_events(1,1),1);
Time_2= data.Timestamp_s((index_of_events(2,1):index_of_events(2,2)),1)-data.Timestamp_s(index_of_events(2,1),1);
Time_3= data.Timestamp_s((index_of_events(3,1):index_of_events(3,2)),1)-data.Timestamp_s(index_of_events(3,1),1);
Time_4= data.Timestamp_s((index_of_events(4,1):index_of_events(4,2)),1)-data.Timestamp_s(index_of_events(4,1),1);
Time_5= data.Timestamp_s((index_of_events(5,1):index_of_events(5,2)),1)-data.Timestamp_s(index_of_events(5,1),1);
Time_6= data.Timestamp_s((index_of_events(6,1):index_of_events(6,2)),1)-data.Timestamp_s(index_of_events(6,1),1);
Time_7= data.Timestamp_s((index_of_events(7,1):index_of_events(7,2)),1)-data.Timestamp_s(index_of_events(7,1),1);
%

% Find Angle using FOV
Angle_FOV = data.Angle_rad(1,FOV_span_idx);

% Convert the polar coordinate to cartesian
[X_whole_1, Y_whole_1] = pol2cart(Angle_FOV,data.Distance_m(index_of_events(1,1):index_of_events(1,2),FOV_span_idx));
[X_whole_2, Y_whole_2] = pol2cart(Angle_FOV,data.Distance_m(index_of_events(2,1):index_of_events(2,2),FOV_span_idx));
[X_whole_3, Y_whole_3] = pol2cart(Angle_FOV,data.Distance_m(index_of_events(3,1):index_of_events(3,2),FOV_span_idx));
[X_whole_4, Y_whole_4] = pol2cart(Angle_FOV,data.Distance_m(index_of_events(4,1):index_of_events(4,2),FOV_span_idx));
[X_whole_5, Y_whole_5] = pol2cart(Angle_FOV,data.Distance_m(index_of_events(5,1):index_of_events(5,2),FOV_span_idx));
[X_whole_6, Y_whole_6] = pol2cart(Angle_FOV,data.Distance_m(index_of_events(6,1):index_of_events(6,2),FOV_span_idx));
[X_whole_7, Y_whole_7] = pol2cart(Angle_FOV,data.Distance_m(index_of_events(7,1):index_of_events(7,2),FOV_span_idx));
% 
% Create Approximation Matrix using starting point and ending point in 
% play Lidar
x_true_1= linspace(6,1,43);
x_true_2= linspace(6,4,55);
x_true_3= linspace(10,2,62);
x_true_4= linspace(7,1.7,55);
x_true_5= linspace(10,0.4,92);
x_true_6= linspace(7,0.6,68);
x_true_7= linspace(3.5,0.3,65);
%
% Transpose the approximation matrix
x_true_t_1 = x_true_1';
x_true_t_2 = x_true_2';
x_true_t_3 = x_true_3';
x_true_t_4 = x_true_4';
x_true_t_5 = x_true_5';
x_true_t_6 = x_true_6';
x_true_t_7 = x_true_7';

% Replicate the distance approximation to make it easier to find diff
x_true_r_1 = repmat(x_true_t_1,1,521);
x_true_r_2 = repmat(x_true_t_2,1,521);
x_true_r_3 = repmat(x_true_t_3,1,521);
x_true_r_4 = repmat(x_true_t_4,1,521);
x_true_r_5= repmat(x_true_t_5,1,521);
x_true_r_6 = repmat(x_true_t_6,1,521);
x_true_r_7 = repmat(x_true_t_7,1,521);
 
% Find difference between approximation and whole real data
diffs_1 = X_whole_1-x_true_r_1;
diffs_2 = X_whole_2-x_true_r_2;
diffs_3 = X_whole_3-x_true_r_3;
diffs_4 = X_whole_4-x_true_r_4;
diffs_5 = X_whole_5-x_true_r_5;
diffs_6 = X_whole_6-x_true_r_6;
diffs_7 = X_whole_7-x_true_r_7;

% Find minimum difference between approximation and whole real data and
% Save the Index
 [val_3_1, idx_3_1]= min(abs(diffs_1),[],2);
 [val_3_2, idx_3_2]= min(abs(diffs_2),[],2);
 [val_3_3, idx_3_3]= min(abs(diffs_3),[],2);
 [val_3_4, idx_3_4]= min(abs(diffs_4),[],2);
 [val_3_5, idx_3_5]= min(abs(diffs_5),[],2);
 [val_3_6, idx_3_6]= min(abs(diffs_6),[],2);
 [val_3_7, idx_3_7]= min(abs(diffs_7),[],2);
 %
 % Find closest point to rough approximation
 for i = 1:length(idx_3_1)
X_whole_closest_1(i) = X_whole_1(i,idx_3_1(i,1));
 end
 for i = 1:length(idx_3_2)
X_whole_closest_2(i) = X_whole_2(i,idx_3_2(i,1));
 end
 for i = 1:length(idx_3_3)
X_whole_closest_3(i) = X_whole_3(i,idx_3_3(i,1));
 end
 for i = 1:length(idx_3_4)
X_whole_closest_4(i) = X_whole_4(i,idx_3_4(i,1));
 end
 for i = 1:length(idx_3_5)
X_whole_closest_5(i) = X_whole_5(i,idx_3_5(i,1));
 end
 for i = 1:length(idx_3_6)
X_whole_closest_6(i) = X_whole_6(i,idx_3_6(i,1));
 end
  for i = 1:length(idx_3_7)
X_whole_closest_7(i) = X_whole_7(i,idx_3_7(i,1));
  end
 
 %Transpose Matrix closest
X_whole_closest_t_1 = X_whole_closest_1.';
X_whole_closest_t_2 = X_whole_closest_2.';
X_whole_closest_t_3 = X_whole_closest_3.';
X_whole_closest_t_4 = X_whole_closest_4.';
X_whole_closest_t_5 = X_whole_closest_5.';
X_whole_closest_t_6 = X_whole_closest_6.';
X_whole_closest_t_7 = X_whole_closest_7.';

% Use smoothdata function to filter
 X_whole_smooth_1 = smoothdata(X_whole_closest_t_1,2);
 X_whole_smooth_2 = smoothdata(X_whole_closest_t_2,2);
 X_whole_smooth_3 = smoothdata(X_whole_closest_t_3,2);
 X_whole_smooth_4 = smoothdata(X_whole_closest_t_4,2);
 X_whole_smooth_5 = smoothdata(X_whole_closest_t_5,2);
 X_whole_smooth_6 = smoothdata(X_whole_closest_t_6,2);
 X_whole_smooth_7 = smoothdata(X_whole_closest_t_7,2);
 %
 % Remove the artifacts in 2nd and 4th Event since the graph strange for
 % both even after smoothing
    d = [0; diff( X_whole_smooth_2 )] ;
    isValid = ~logical( cumsum( -sign( d ) .* (abs( d ) > 1) )) ;
    t = 1 : numel( X_whole_smooth_2 ) ;
    X_whole_smooth_2_clean = interp1( t(isValid), X_whole_smooth_2(isValid), t ) ;
    X_whole_smooth_2_clean_t = X_whole_smooth_2_clean.';
    
    d = [0; diff( X_whole_smooth_4 )] ;
    isValid = ~logical( cumsum( -sign( d ) .* (abs( d ) > 1) )) ;
    t = 1 : numel( X_whole_smooth_4 ) ;
    X_whole_smooth_4_clean = interp1( t(isValid), X_whole_smooth_4(isValid), t ) ;
    X_whole_smooth_4_clean_t = X_whole_smooth_4_clean.';
 
    % Plot the Distance Events
 figure(2)
 plot(Time_1(:,1),X_whole_smooth_1, Time_2(:,1),X_whole_smooth_2_clean,Time_3(:,1),X_whole_smooth_3,...
     Time_4(:,1),X_whole_smooth_4_clean,Time_5(:,1),X_whole_smooth_5,Time_6(:,1),X_whole_smooth_6,...
     Time_7(:,1),X_whole_smooth_7,'linewidth', 1.5)
 hold on 
 xlabel('Time [s]'); ylabel('Distance [m]')
 legend('1st Event','2nd Event','3rd Event','4th Event','5th Event','6th Event','7th Event')

%% PLOT RELATIVE SPEED TO TARGET ACROSS TIME
figure (3);
hold on
   % 1st Event Relative Velocity
   Relative_distance_change_1 = gradient(X_whole_smooth_1);
   Relative_time_change_1 = gradient(Time_1);
   a = Relative_distance_change_1./Relative_time_change_1;
   a_smooth = smoothdata(a,2);
   % 2nd Event Relative Velocity
   Relative_distance_change_2 = gradient(X_whole_smooth_2_clean_t);
   Relative_time_change_2 = gradient(Time_2);
   b = Relative_distance_change_2./Relative_time_change_2;
   b_smooth = smoothdata(b,2);
   % 3rd Event Relative Velocity
   Relative_distance_change_3 = gradient(X_whole_smooth_3);
   Relative_time_change_3 = gradient(Time_3);
   c = Relative_distance_change_3./Relative_time_change_3;
   c_smooth = smoothdata(c,2);
   % 4th Event Relative Velocity
   Relative_distance_change_4 = gradient(X_whole_smooth_4_clean_t);
   Relative_time_change_4 = gradient(Time_4);
   d = Relative_distance_change_4./Relative_time_change_4;
   d_smooth = smoothdata(d,2);
   % 5th Event Relative Velocity
   Relative_distance_change_5 = gradient(X_whole_smooth_5);
   Relative_time_change_5 = gradient(Time_5);
   e = Relative_distance_change_5./Relative_time_change_5;
   e_smooth = smoothdata(e,2);
   % 6th Event Relative Velocity
   Relative_distance_change_6 = gradient(X_whole_smooth_6);
   Relative_time_change_6 = gradient(Time_6);
   f = Relative_distance_change_6./Relative_time_change_6;
   f_smooth = smoothdata(f,2);
   % 7th Event Relative Velocity
   Relative_distance_change_7 = gradient(X_whole_smooth_7);
   Relative_time_change_7 = gradient(Time_7);
   g = Relative_distance_change_7./Relative_time_change_7;
   g_smooth = smoothdata(g,2);
   
   %Plot the Relative Velocity Events
   plot(Time_1(:,1),a_smooth,Time_2(:,1),b_smooth,Time_3(:,1),c_smooth,...
       Time_4(:,1),d_smooth,Time_5(:,1),e_smooth,Time_6(:,1),f_smooth,...
       Time_7(:,1),g_smooth,'linewidth', 1.5)
xlabel('Time [s]'); ylabel('Relative speed [m/s]')
 legend('1st Event','2nd Event','3rd Event','4th Event','5th Event','6th Event','7th Event')

%% PLOT TIME TO COLLISION TO TARGET ACROSS TIME
figure (4);
hold on
% Calculate the TTC = distance/relative speed
TTC_1 = abs(X_whole_smooth_1./a_smooth);
TTC_2 = abs(X_whole_smooth_2./b_smooth);
TTC_3 = abs(X_whole_smooth_3./c_smooth);
TTC_4 = abs(X_whole_smooth_4./d_smooth);
TTC_5 = abs(X_whole_smooth_5./e_smooth);
TTC_6 = abs(X_whole_smooth_6./f_smooth);
TTC_7 = abs(X_whole_smooth_7./g_smooth);

% Smoothen the TTC
TTC_1_smooth = smoothdata(TTC_1,2);
TTC_2_smooth = smoothdata(TTC_2,2);
TTC_3_smooth = smoothdata(TTC_3,2);
TTC_4_smooth = smoothdata(TTC_4,2);
TTC_5_smooth = smoothdata(TTC_5,2);
TTC_6_smooth = smoothdata(TTC_6,2);
TTC_7_smooth = smoothdata(TTC_7,2);

% Plot the TTC
   plot(Time_1(:,1),TTC_1_smooth,Time_2(:,1),TTC_2_smooth,Time_3(:,1),TTC_3_smooth,...
       Time_4(:,1),TTC_4_smooth,Time_5(:,1),TTC_5_smooth,Time_6(:,1),TTC_6_smooth,...
       Time_7(:,1),TTC_7_smooth,'linewidth', 1.5)

xlabel('Time [s]'); ylabel('TTC [s]')
legend('1st Event','2nd Event','3rd Event','4th Event','5th Event','6th Event','7th Event')
%% --- Question ---
% What safety measure would you use to design a warning that alert the user
% that is about to collide with an obstacle? You  may want to use one of the safety measure 
% you computed in this script of find a more effective one.
% What value of such measure would you use to trigger a warning?

% Answer: The user can be alerted by creating beep sound and flashing
% light on the headsup display or vibrate the steering.Form the safety measure 
% that was computed TTC is more effective to trigger a warning to the user when compared
% to relative speed and distance.The warning can be triggred at higher value
% for TTC to keep the user safe as TTC give the time for two vehicles before
% collision.
