41%% ----- Exercise 01 - Decode CAN data -----
% Version: 2018
% Course: TME 192 Active Safety
%         Chalmers
% Author: Alberto Morando (morando@chalmers.se)
%         Marco Dozza (dozza@chalmers.se)

clear all
close all
clc

%% ### Auxiliary functions
% Function to translate hex to binary
hex2bin = @(data) dec2bin(hex2dec(data), 8);

%% ### Import the data
% There are multiple ways to import data from a txt file. One way is to use 'readtable' (read the doc for more
% information). Note that 'format' is a list of `%x` where `x` represent
% the different data types (e.g. `s` for string, `d` for decimal). For the
% exercise, we already defined the format for you.
FORMAT_SPEC = '%f%d%s%s%c%d%s%s%s%s%s%s%s%s';
CANdata = readtable('CANlog.txt', 'format', FORMAT_SPEC);

% Add variable names. This would be helpful for analysing the data
CANdata.Properties.VariableNames{1} = 'Timestamp';
CANdata.Properties.VariableNames{3} = 'Identifier';
%CANdata.Properties.VariableNames{7} = 'Message(hex)';

%% ### Decode the signals
% Select the frames based on the identifier. You could use `strcmpi` on
% the `Identifier` column
str1 = strcmpi(CANdata.Identifier,'401');
str2 = strcmpi(CANdata.Identifier,'1CC');
%% * Find vehicle speed signal to decode *
speed_index = str1;

% Extract the bytes of interest
speed_hex = CANdata{speed_index,end};

% Convert hex -> bin 
speed_bin = hex2bin(speed_hex);

% convert bin -> dec (use 'bin2dec'). Apply the 'factor' and `offset`
speed = bin2dec(speed_bin);
speed_final = speed.*1;
%% * Find the accelerator pedal position signal to decode *

accelerator_index = str2;

% Extract the bytes of interest
accelerator_hex_reshaped = CANdata{accelerator_index,end-3:end-2};
acc1 = accelerator_hex_reshaped(:,1);
acc2 = accelerator_hex_reshaped(:,2);
accelerator_hex = strcat(acc1,acc2);

% Convert hex -> bin 
accelerator_bin_sep = hex2bin(accelerator_hex);
accelerator_bin= accelerator_bin_sep(:,7:16);
    
% convert bin -> dec (use 'bin2dec'). Apply the 'factor' and `offset`
accelerator = bin2dec(accelerator_bin).*0.098;
%accelerator_final = accelerator.*0.098;

%% Plot the signals with respect to time. Include labels and legend
speed_time = CANdata{speed_index,1};
accelerator_time = CANdata{accelerator_index,1};

figure(1)
plot(speed_time,speed_final)
hold on
plot(accelerator_time,accelerator)
xlabel('Timestamp [s]')
ylabel('Vehicle Speed [km/h]')
legend('Vehicle Speed','Accelerator pedal position')
hold off

%% Downsample the signals
% Remove the delay between the two signals. That is, make sure both signals
% starts from zero.
speed_time_from_0 = speed_time;
accelerator_time(accelerator_time<0.05) = 0;
accelerator_time_from_0 = accelerator_time; %manually inputed from accelerator_time, excl. the first element

% Create a master clock, at 1Hz, common to both signals
sync_time = linspace(0,40,41);

% Downsample by interpolation (use 'interp1')
speed_sync = interp1(speed_time_from_0,speed_final,sync_time);
accelerator_sync = interp1(accelerator_time_from_0,accelerator,sync_time);

%% Plot the signals with respect to time. Include labels and legend. Limit the axis.

figure(2)
plot(sync_time,speed_sync,'bo-',speed_time,speed_final,'b--')
hold on
plot(sync_time,accelerator_sync,'ro-',accelerator_time,accelerator,'r--')
xlim([28 36])
ylim([15 30])
xlabel('Timestamp [s]')
ylabel('Vehicle Speed [km/h]')
legend('Vehicle Speed_ SYNC','Vehicle Speed','Accelerator pedal position','Accelerator pedal position_ SYNC')
hold off
%% Issue a warning when the driver is about to exceed the speed limit
% Speed limit: 50 km/h, upper tolerance: 50 + 10 km/h
% Driver reaction time: 1.0
% SPEED_LIMIT = 50;
% TOLERANCE = 10;
% REACTION_TIME = 1.0;

R_t= 1;
v_tolerance=60/3.6;
a=zeros(1,length(speed_sync));
a(2:end)=diff(speed_sync./3.6);


predicted_speed=speed/3.6 + a*R_t

predicted_speed=speed_sync(1)

k=1
while predicted_speed < v_tolerance
predicted_speed=speed_sync(k)/3.6 +a(k)*R_t;
warning_time=sync_time(k)
pause(0.30);
k=k+1
end
fprintf(['Warning activated at time ' num2str(warning_time) ' s.']);

    
return