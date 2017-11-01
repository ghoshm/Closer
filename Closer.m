% Closer 

%% Data Import 

% Data should be delimited in Excel to consist of 
 % Rows - frames 
 % Columns - No. of fish * 2 (lights then fish) 

% Data should then be imported as three matricies 
 % Lights_Binary - the light status of each fish @ every frame as either: 
  % 0 - off 
  % 1 - on 
 Lights_Binary = Untitled; 
 clear Untitled 
 
 % Fish -  delta pixels for each fish @ every frame
 Data = Untitled; 
 clear Untitled 

 % Defining Day vs Night Periods (Day = 1, Night = 0) 
 DN = Untitled;
 clear Untitled
 DN(DN >= 23) = 0;
 DN(DN < 9) = 0;
 DN(DN > 1) = 1;

% DN(1:size(Lights_Binary,1),1) = 0; 
% DN(1:(find(sum(Lights_Binary,2) == 0,1,'first'))-1,1) = 1; 
% DN(find(sum(Lights_Binary,2) == 0,1,'last')+1:end,1) = 1; 

% Hard Coded Variables 
% Frame Rate % Inactive period length that triggers the closed loop  
fps = 40;
ip = 45; % In seconds)  

%% Variable Building 

% V1 Light Ticks - Roughly 0.5s per fish 
 % A 1 at dark/light transitions and the 
 % Response latency at the response point  
 % Note that Day Light is not tagged 
 
tic
Light_Ticks = Lights_Binary; 
%Light_Ticks(end,:) = 0; % Prevents counting past the end of the experiment 

% Removing the day Lights on
for f = 1:size(Light_Ticks,2) % For each fish  
    scrap = Light_Ticks(:,f);
    scrap(DN == 1) = 0; % Remove the Day Light Ticks 
    Light_Ticks(:,f) = scrap; 
    clear scrap 
end 

progressbar('Fish') %Initialise progress bars 
for f = 1:size(Light_Ticks,2) % For each fish 
    scrap = diff(Light_Ticks(:,f)); % Diff to find Transitions   
    locs_start = find(scrap == 1); % Find Lights On 
    locs_start = locs_start + 1; % Convert Time 
    
    for t = locs_start' % For each Light on Stimulus  
            end_point = find(Light_Ticks(t:end,f) == 0,1,'first') -1; % Find lights off  
            Light_Ticks(((t+end_point)-1),f) = end_point; % Tag end point
            Light_Ticks(t+1:((t+end_point)-2),f) = 0; % Clear middle points 
            clear end_point
    end
    
    clear scrap locs_start 
    progressbar(f/size(Light_Ticks,2)); %Update the Fish progressbar
    
end 

clear f t  
toc

%% V2 Inactive Bouts - Roughly 5 mins per fish  
% Sleep - a 1 @ active to inactive transitions and a number at the end of 
    % each inactive bout counting it's length 
Sleep = Data;
for f = 1:size(Light_Ticks,2) % For each fish 
    scrap = Sleep(:,f); 
    scrap(scrap == 0) = 1; % Inactivity == 1 
    scrap(scrap > 1) = 0; % Activity == 0 
    Sleep(:,f) = scrap; 
    clear scrap 
end 

% Inactive Bout Length Loop  
tic 
progressbar('Fish','Inactive Bout') %Initialise progress bars 
for f = 1:size(Light_Ticks,2) % For each fish 
  
    scrap = diff(Sleep(:,f)); % Diff to find transitions   
    locs_start = find(scrap == 1); % Find Inactive Bout Starts  
    locs_start = locs_start + 1; % Convert Time 
    
    a = 1; % Start inactive bout counter 
    for t = locs_start' % For each inactive bout 
            end_point = find(Sleep(t:end,f) == 0,1,'first') -1; % Find it's end   
            Sleep(((t+end_point)-1),f) = end_point; % Tag end point
            Sleep(t+1:((t+end_point)-2),f) = 0; % Clear middle points 
            clear end_point
            progressbar([],a/size(locs_start,1)); %Update the Inactive Bout progressbar
            a = a+1; % Add to counter  
    end
    
    clear a locs_start scrap  
    progressbar(f/size(Light_Ticks,2)); %Update the Fish progressbar
end 

clear f t 
toc 

%% Re-scale the Data 

for f = 1:size(Data,2) % For each fish  
    top = max(Data(:,f)); % Find the max value 
    % Note that this assumes that zero will be the lowest value 
    % This should be the case! 
    Data(:,f) = Data(:,f) / top; % Rescale 
end 

%% V3 Response Latencies, Windows & Targets  
  % response_latencies - a cell for each fish storing the time points (:,1) 
   % & latencies of each repsonse (:,2)
  % response_latencies_bin - a single vector containing all response
   % latencies (in frames) 
  % response_windows - 1 second of data either side of each response with
   % the response centered @ the fps point 
  % targets - a cell: 
    % Rows: 
        % 1 - target fish response windows 
        % 2 - non-target fish response windows 
    % Columns - each fish 

tic
response_latencies_bin = []; 
progressbar('Fish') %Initialise progress bars 
for f = 1:size(Light_Ticks,2) % For each fish 
        
    % Response Latencies & Windows 
    locs = find(Light_Ticks(:,f) > 1); % Find the responses 
    locs = locs(locs > (fps - 1)); % Remove response too close to the start of the experiment 
    locs = locs(locs < (size(Data,1) - fps)); %Remove responses too close to the end of the experiment 

    response_latencies{1,f}(1:size(locs,1),1) = locs; % Store time points
    response_latencies{1,f}(1:size(locs,1),2) = Light_Ticks(locs,f); % Store latencies
    response_latencies_bin = [response_latencies_bin' Light_Ticks(locs,f)']'; % Store latencies 
        
    for r = 1:size(locs,1) % For each response
        response_windows{1,f}(r,1:fps*2) = ...
            Data((locs(r)-(fps-1)):(locs(r)+fps),f); % Store 1s of Data either side 
    end
    
    % Targets 
    target_window = fps*1; % Look xs either side of the stimulus
    lights = find(Light_Ticks(:,f) == 1); % Find when the lights flash 
    lights = lights(lights > (target_window - 1)); % Remove stimuli too close to the start of the experiment 
    lights = lights(lights < (size(Data,1) - target_window)); % Remove stimuli too close to the end of the experiment 
    
    total = [];
    non_target_fish = setdiff(1:size(Light_Ticks,2),f); % Use fish other than the "target" fish 
    
    for l = 1:size(lights,1) % For each light flash      
        targets{1,f}(l,1:target_window*2) = ...
                    Data((lights(l)-(target_window-1)):(lights(l)+target_window),f); % Store xs of Data either side
        for f2 = non_target_fish % For each non target fish 
            if sum(Lights_Binary((lights(l)-(target_window-1)):(lights(l)+target_window),f2))...
                    == 0 % & The Lights are off for this fish 
                if isempty(total) == 1 % If its the first fish (conditions)
                    total = Data((lights(l)-(target_window-1)):(lights(l)+target_window),f2)'; % Store 1s of Data either side
                else
                    total = [total ; Data((lights(l)-(target_window-1)):(lights(l)+target_window),f2)']; % Add it to this vector
                end
            end
        end         
    end
    
    targets{2,f} = total; % Save the data for non-target fish under the target fish 
    
    clear locs lights total non_target_fish 
    
    progressbar(f/size(Light_Ticks,2)); %Update the Fish progressbar
end 

clear f r l f2
toc 

%% Figure Working

% Light Status Figures 
for f = 1:size(Data,2) % For each fish 
    
    figure; hold on; 
        smoothed = smooth(Data(:,f),fps*60); % Smooth the data into minutes 
        % Note that smoothing into longer windows looks nicer but can
        % easily smooth into the inactive periods 
        top = max(smoothed); % Find the maximum delta pixels 
        % Plot a rectangle over the night period 
        rectangle('Position',[find(DN==0,1,'first') 0 ...
            ((find(DN==0,1,'last')-find(DN==0,1,'first'))+1) (top+(top*0.1))],...
            'FaceColor',[0.9608    0.9608    0.9608],'Edgecolor',[1 1 1]); 
    
    for r = 1:size(response_latencies{1,f},1) % For each response 
        rectangle('Position',[(response_latencies{1,f}(r,1) - ...
            (response_latencies{1,f}(r,2) - 1)) 0 response_latencies{1,f}(r,2) ...
                (top+(top*0.1))],... 
                'FaceColor',[1 1 1],'Edgecolor',[1 1 1]); % Place a white rectangle 
    end 
   
    plot(smoothed); % Plot the Data  
    axis([(fps*2), size(Data,1)-(fps*2),...
        0, max(smoothed(fps*2:size(Data,1)-(fps*2)))]); % Set the axis 
    clear smoothed top 
    
end 

% Response Latencies - frames 

    % Sequential 
    figure; hold on; 
    for f = 1:size(response_latencies,2) % For each fish 
        plot(response_latencies{1,f}(:,2)); %Plot each latency 
    end;

    % Over Time 
    figure; hold on; 
    for f = 1:size(response_latencies,2) % For each fish 
        scatter(response_latencies{1,f}(:,1),response_latencies{1,f}(:,2),'filled')
    end;

    % Over Time - Scatter Hist 
    for f = 1:size(response_latencies,2) % For each fish 
        if f == 1
            response_latencies_hist = response_latencies{1,f}; 
        else 
            response_latencies_hist = [response_latencies_hist ; response_latencies{1,f}]; 
        end 
    end 
    
    figure; hold on; 
    scatterhist(response_latencies_hist(:,1),response_latencies_hist(:,2),'Kernel','on')
    
    % Mean 
    notBoxPlot(response_latencies_bin)
    nanmedian(response_latencies_bin)/fps 
    
% Response Windows 
figure; hold on;
for f = 1:size(response_windows,2) % For each fish
    plot(response_windows{1,f}','color',[0.4392    0.5020    0.5647]); 
     % Plot each response 
end

% Sleep 
t_sleep = Sleep; % Thresholded Sleep 
t_sleep(t_sleep < fps * 60) = NaN; % Apply a Hard Coded Definition of sleep 

% Totals in each period  
for f = 1:size(Data,2) % For each fish 
    % Day 1 
    sleep_total(1,f) = nansum(t_sleep(1:find(DN==0,1,'first')-1,f)); % Sum the sleep
    sleep_bout_counts(1,f) = size(find(isnan(t_sleep(1:find(DN==0,1,'first')-1,f))==0),1); % Count the 
        % Number of bouts 
    
    % Night 1 
    sleep_total(2,f) = nansum(t_sleep(find(DN==0,1,'first'):find(DN==0,1,'last'),f));
    sleep_bouts_counts(2,f) = size(find(isnan(t_sleep(find(DN==0,1,'first'):find(DN==0,1,'last'),f))==0),1); 
    
    % Day 2 
    sleep_total(3,f) = nansum(t_sleep(find(DN==0,1,'last')+1:end,f));
    sleep_bouts_counts(3,f) = size(find(isnan(t_sleep(find(DN==0,1,'last')+1:end,f))==0),1); 
    
end 

plotSpread((t_sleep(find(DN==0,1,'first'):find(DN==0,1,'last'),:))/(fps*60))
plot(sleep_total)
plot(sleep_bouts_counts) % Plot sleep bouts across periods 

% Across time 
for f = 1:size(Sleep,2) % For each fish 
    
    locs = find(isnan(t_sleep(:,f))==0); % Find Sleep Bouts (ends)
    sleep_bouts{1,f}(:,1) = (locs - t_sleep(locs,f)) + 1; % Store Sleep Bout Starts 
    sleep_bouts{1,f}(:,2) = t_sleep(locs,f); % Store Sleep Bout Lengths 
    
    clear locs 
end 

% Sequential 
    figure; hold on; 
    for f = 1:size(sleep_bouts,2) % For each fish 
        if f <= 3
            plot(sleep_bouts{1,f}(:,2),'r'); %Plot each sleep bout length  
        else 
            plot(sleep_bouts{1,f}(:,2),'b'); %Plot each sleep bout length  
        end 
    end;

    % Over Time 
    figure; hold on; 
    for f = 1:size(sleep_bouts,2) % For each fish 
        if f <= 3
            scatter(sleep_bouts{1,f}(:,1),sleep_bouts{1,f}(:,2),'r','filled')
        else 
            scatter(sleep_bouts{1,f}(:,1),sleep_bouts{1,f}(:,2),'b','filled')
        end 
    end;

    % Over Time - Scatter Hist 
    clear groups sleep_bouts_hist
    for f = 1:size(sleep_bouts,2) % For each fish 
        if f == 1
            sleep_bouts_hist = sleep_bouts{1,f}; 
            groups(1:size(sleep_bouts{1,f},1),1) = 0; 
        else 
            sleep_bouts_hist = [sleep_bouts_hist ; sleep_bouts{1,f}]; 
            
            if f <= 3 
                groups = [groups ; zeros(size(sleep_bouts{1,f},1),1)];
            else 
                groups = [groups ; ones(size(sleep_bouts{1,f},1),1)];
            end 
        end 
    end 
    
    
    figure; hold on; 
    scatterhist(sleep_bouts_hist(:,1),sleep_bouts_hist(:,2),'Group',groups,'Kernel','on')
    
%% Figures for Show 

%% Light Isolation 
for t = 1:size(targets,2) % For each fish  
    if t == 1 % If it's the first fish 
        yt = targets{1,1}; % Yes Target (yt) = targeted fish response windows 
        nt = targets{2,1}; % Not Target (nt) = the other fish (with lights off) 
    else 
        yt = [yt' targets{1,t}']'; % Add to Yes Target 
        nt = [nt' targets{2,t}']'; % Add to Not Target 
    end 
end

figure; hold on;
a = shadedErrorBar(1:size(nt,2),nanmean(nt),nanstd(nt)/sqrt(size(nt,1)),'b');
b = shadedErrorBar(1:size(yt,2),nanmean(yt),nanstd(yt)/sqrt(size(yt,1)),'r');

%Nice Figure; 
box off; 
set(gca,'Fontsize',44) 
xlabel('Time (Frames)')
ylabel('Activity (delta px)')
x = [size(nt,2)/2 size(nt,2)/2] 
scrap = ylim; 
y = [0 scrap(2)]; clear scrap 
plot(x,y,'k--','linewidth',3)
set(gca, 'Layer','top')
[h,icons,plots,str] = legend([a.mainLine,b.mainLine], ...
    horzcat('No Stimulus (e = ',num2str(size(nt,1)),')') ...
    ,horzcat('Stimulus (e = ',num2str(size(yt,1)),')')...
    ,'location','northwest');
set(h,'Fontsize',30)
set(icons(1:2),'Fontsize',34)
set(icons,'Linewidth',20)
legend('boxoff')

%% Example Traces - cutting 5mins either side 
figure; 
for f = 1:size(Data,2) % For each fish 
    smoothed(:,f) = smooth(Data((fps*60*5):size(Data,1)-(fps*60*5),f),fps*60);
     % Note that smoothing into longer windows looks nicer but can
        % easily smooth into the inactive periods 
end 
smoothed_top = max(max(smoothed)/4); % Find the highest value & Quarter (Hard Coded   

for f = 1:size(Data,2) % For each fish 
    subplot(2,3,f)% Hard coded
    hold on; 
        % Plot a rectangle over the night period 
        rectangle('Position',[find(DN==0,1,'first') 0 ...
            ((find(DN==0,1,'last')-find(DN==0,1,'first')) + 1) (smoothed_top+(smoothed_top*0.1))],...
            'FaceColor',[0.9608    0.9608    0.9608],'Edgecolor',[1 1 1]); 
    
    for r = 1:size(response_latencies{1,f},1) % For each response 
        rectangle('Position',[(response_latencies{1,f}(r,1) - ...
            (response_latencies{1,f}(r,2) - 1)) 0 response_latencies{1,f}(r,2) ...
                (smoothed_top+(smoothed_top*0.1))],... 
                'FaceColor',[1 1 1],'Edgecolor',[1 1 1]); % Place a white rectangle 
    end
    
    if f <= 3
       plot(smoothed(:,f),'r'); % Plot the Data  
    else 
       plot(smoothed(:,f),'b')
    end 
    
    axis([(fps*60*5), size(Data,1)-(fps*60*5),...
        0, smoothed_top]); % Set the axis 
    set(gca,'fontsize',22) 
    set(gca, 'Layer','top')
    xlabel('Time (Frames)') 
    ylabel ('Activity (delta px)');
    
end  

%% Sleep Deprivation Figure
figure; 
subplot(1,2,1); hold on; 
rectangle('Position',[0 0 7 ...
    max(max((t_sleep(find(DN==0,1,'first'):find(DN==0,1,'last'),:))/(fps*60)))],...
            'FaceColor',[0.9608    0.9608    0.9608],'Edgecolor',[1 1 1]);
plotSpread((t_sleep(find(DN==0,1,'first'):find(DN==0,1,'last'),:))/(fps*60),...
    'distributionColors',...
    {'r','r','r',[0.3922    0.5843    0.9294],[0.3922    0.5843    0.9294]...
    [0.3922    0.5843    0.9294]}) % Plot Sleep Bouts During the Night 
title('Night Sleep Bout Lengths') 
xlabel('Fish'); ylabel('Sleep Bout Length (Mins)'); 
set(gca,'Fontsize',34); set(gca,'Layer','top'); 
axis([0 7 0 ...
    max(max((t_sleep(find(DN==0,1,'first'):find(DN==0,1,'last'),:))/(fps*60)))])

subplot(1,2,2); hold on; 
rectangle('Position',[1.5 0 2 (max(max(sleep_total)))/(fps*60)],...
            'FaceColor',[0.9608    0.9608    0.9608],'Edgecolor',[1 1 1]);
plot((sleep_total(:,4:6))/(fps*60),'color',[0.3922    0.5843    0.9294],'linewidth',3)
plot((sleep_total(:,1:3))/(fps*60),'r','linewidth',3)
title('Total Sleep') 
xlabel('Time Window'); ylabel('Total Sleep (Mins)'); 
set(gca,'Fontsize',34); set(gca,'Layer','top'); 
axis([0.5 3.5 0 (max(max(sleep_total)))/(fps*60)])


%% Merging Data 

     
%% Video Analysis 

% Load Video

% Select a video file
[filename, pathname] = uigetfile('*.avi', 'Select a Video file'); %Select a geno file
if isequal(filename,0) %If no file is selected
    error('No File Selected') %Show Error
else %If selected
    disp(['User selected ', fullfile(pathname, filename)]) %Show selected filename
end
v = VideoReader(strcat(pathname,filename));
frames = 1; 
    while hasFrame(v) % Read all video frames 
        video(:,:,:,frames) = readFrame(v);
        frames = frames + 1; % Count all video frames 
    end
clear frames 

%% Making the Video 

%Positions Vector 
sc = [160 , 160]; % Starting coordinate  
    x = 350; % X distance between wells  
    y = 390; % Y distance between wells  
 
p = 1; % For each position 
for r = 1:2 % For each row (Hard coded)
    for c = 1:3 % For each column (Hard coded)
        well_positions(p,:) = [sc(1)+((c*x)-x),sc(2)+((r*y)-y)];
        p = p+1;
    end
end;

% Axis (0.15 by 0.15 size) - Hard Coded Locations 
axis_positions(1,:) = [0.15 0.50]; 
    axis_positions(2,:) = [0.43 0.50]; 
    axis_positions(3,:) = [0.71 0.50]; 
    axis_positions(4,:) = [0.15 0.13]; 
    axis_positions(5,:) = [0.43 0.13]; 
    axis_positions(6,:) = [0.71 0.13]; 
    
% Making the Video 
    % Pre-allocate a strucutre to store video frames 
s(size(video,4)) = struct('cdata',[],'colormap',[]); 

data_frames = 22140:22539; 

set(0,'DefaultFigureVisible','off') % Suppress figures from popping up 

for f = 1:size(video,4) %For each frame 
    
    h = axes('position',[0 0 1 1]); % Initial axis  
    image(video(:,:,:,f)); hold on; %Draw the frame 
    set(gca,'XTick',[]); % Remove x-axis ticks 
    set(gca,'YTick',[]); % Remove y-axis ticks 

    % Light Status 
    
    for p = 1:size(well_positions,1) % For each well  
        if Lights_Binary(data_frames(f),p) == 0 % If the lights are off 
            scatter(well_positions(p,1),well_positions(p,2),300,'MarkerFaceColor',...
            [0.9608    0.9608    0.9608],'MarkerEdgeColor','none'); % Draw a White Dot 
        else % If the lights are on 
            scatter(well_positions(p,1),well_positions(p,2),300,'MarkerFaceColor',...
            'y','MarkerEdgeColor','none'); % Draw a yellow dot  
        end
    end
    
    % Fish Behaviour 
    
    for a = 1:size(axis_positions,1) % For each axis  
        axis_structure(a) = axes('position',[axis_positions(a,1),axis_positions(a,2)...
            0.15 0.15],'ylim',[0 inf]); % Position an axis 
        hold on; % Set hold on to this axis 
        plot(axis_structure(a),1:fps,Data...
            ((data_frames(f)-fps)+1:data_frames(f),a),'color',...
            [0.3922    0.5843    0.9294],'linewidth',3); % Plot 1s of data before this frame
         % To this frame 
        scatter(axis_structure(a),fps,Data...
            (data_frames(f),a),45,'MarkerFaceColor',...
            'k','MarkerEdgeColor','none'); % Draw a black Dot at the current value 
        axis('off');
        
    end 
    
    % Storing frames with drawings on top 
    drawnow;
    s(f) = getframe(h); % Grab the frame of the largest axis 
    f % Report back the frames 
    close all % ! % Note this two commands incrase the writing speed 
    hold off % !  % By decreasing the amount of data held in RAM 
    
end; 

set(0,'DefaultFigureVisible','on') % Allow figures to pop up again 

% Open a new figure and play the movie from the structure
movie(s,1,v.FrameRate) % Play the movie at the native framerate  

%Write the Video 
vOut = VideoWriter('C:\Users\Marcus\Desktop\Mpeg4.avi','MPEG-4');
vOut.FrameRate = v.FrameRate/4;
open(vOut)
for k = 1:numel(s)
    writeVideo(vOut,s(k)); 
end
close(vOut)



    

 


