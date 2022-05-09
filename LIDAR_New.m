% Student template script for TME192 LIDAR exercise
% 190901 Jonas B?rgman


% Load the data if not available. you may have to set specific path
if ~exist('oData')
    load('oData')
end


% Initiate a plot
fig=figure(1);

% Set the coordinates for what to show
fPlotCoordsX=[6407050,6407120];
fPlotCoordsY=[1276550,1276650];

% Initiate an AVI. 
% STUDENT: YOU HAVE TO CHANGE THE PATH!
aviobj = VideoWriter(['D:\lid' datestr(now,30) '.avi']);
open(aviobj);

% Loop through all times in the Sensor Fused data
for iIndex=1:length(oData.iTimeSF)

   % Get the specific time for this index from Sensor Fusion data		
   time=oData.iTimeSF(iIndex);
   
   % Find the closest LIDAR time corresponding to the Sensor Fusion time 
   iLIDARIndex=find(oData.iLidarTime>time,1);

   
   % FROM HERE ON STUDENT CODE
   
   % Do the translations and coordinate transformations to extract the
   %   LIDAR reflections in the coordinate system of RT90 (GPS antenna
   %   mounting position)
    X1 =((oData.fLIDAR_X{iIndex}) + oData.fLIDARposX - oData.fGPSposX);
    Y1 =((oData.fLIDAR_Y{iIndex}) + oData.fLIDARposY -oData.fGPSposY);
   % Add the RT90 position (global coordinates from GPS), but in order to 
   %  be able to add them vehicle data it will have to be projected on the 
   %  RT90 coordinate system using the heading.
      Positioning_X = (X1 *cos(oData.fHeadingSF(iIndex))-Y1*sin(oData.fHeadingSF(iIndex)));
      Positioning_Y = (X1 *sin(oData.fHeadingSF(iIndex))+Y1*cos(oData.fHeadingSF(iIndex)));
           % Add to the RT90 cartesian coordinate system. When the code for
           % the two global coordinates is ready, the output should be in 
           % these two variables. They should be the final output pasted 
           % into the function in the Matlab Grader. If you have multiple
           % lines of code to create the two variables, all should be
           % pasted into Matlab Grader.
    
           fXechoGlobal= ( Positioning_X)+ oData.fXRT90SF(iIndex);
           fYechoGlobal= ( Positioning_Y)+ oData.fYRT90SF(iIndex); 

   % END OF STUDENT CODE (if you want, more can be added)

   % Plot single LIDAR
  %plot(oData./LIDARX{1},oData./LIDARY{1},'.')
   % Plot the lidar in RT90 coodrinate system	   
   plot(fXechoGlobal,fYechoGlobal,'.b')
  hold on 
   % Plot the vehicle position (the GPS antenna) too
   plot(oData.fXRT90SF(iIndex),oData.fYRT90SF(iIndex),'.r','MarkerSize',30)
 hold off
   
   % Add your name to the plot
   %%% STUDENT: You should change this to be a list of all the names in your group
   text(6407050,1276560,'Fikri Farhan Witjaksono, Adipta Laha, Aashish Udaykumar Kodgi')
   
   % Set the axis of the plot
   axis([fPlotCoordsX fPlotCoordsY])
   hold on;
   
   % Get it as an avi-frame
   F = getframe(fig);
   % Add the frame to the avi
   writeVideo(aviobj,F);
   %aviobj = addframe(aviobj,F);
   
end

% Close the AVI from Matlab
close(aviobj);
