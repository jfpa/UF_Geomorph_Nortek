
%% Load Nortek data files and store according to waves and currents data
% This function uses the both ways that the files are converted by the
%   program used to extract the data directly from the instruments. Vary
%   according to the update of the sofware.
% The program converts the AWAC data. (the same as Aquadopp but including
%   STrk data).
% Path and filename
% % % [pathname,filename]=filepath_data(a,b,c,d,e);
function [curr,wave] = ...
    load_rawdata_Nortek_AWAC(pathname,filename,coord,ins);
clc, close all
fprintf('======================================\n')
fprintf('LOAD RAW DATA FROM NORTEK AWAC\n')
fprintf('Geomorphology Laboratory\n')
fprintf('Department of Geological Sciences\n')
fprintf('University of Florida\n')
fprintf('Gainesville, FL, USA\n')
fprintf('Summer of 2016\n')
fprintf('======================================\n')

tic;

fprintf('1. LOAD TEXT FILES...\n')

% Load files for Nortek formatting
sen = load([pathname filename '.sen']);
fprintf('.sen file loaded\n')

v1 = load([pathname filename '.v1']);
v2 = load([pathname filename '.v2']);
v3 = load([pathname filename '.v3']);
fprintf('.v1, .v2, and .v3 files loaded\n')

a1 = load([pathname filename '.a1']);
a2 = load([pathname filename '.a2']);
a3 = load([pathname filename '.a3']);
fprintf('.a1, .a2, and .a3 files loaded\n')

wad = load([pathname filename '.wad']);

fprintf('.wad file loaded\n')

whd = load([pathname filename '.whd']);
fprintf('.whd file loaded\n')


% Deployment information
% Retrieve data from .hdr file
fprintf('\nDeployment information from .hdr file\n')
fprintf('----------------------------------------------\n')

fid = fopen([pathname filename '.hdr']);

for k = 1:22
    % Read current line within text
    tline = fgets(fid);
    
    if k == 11 % Number of cells
        disp(tline)
        Ncells = cell2mat(textscan(tline,'%*s %*s %*s %f %*s'));
        
    elseif k == 12 % Cell size, m
        disp(tline)
        cellsize = cell2mat(textscan(tline,'%*s %*s %f %*s'))/100;
        
    elseif k == 15 % Transmit pulse length, m
        disp(tline)
        Lp = cell2mat(textscan(tline,'%*s %*s %*s %f %*s'));

    elseif k == 16 % Blanking distance, m
        disp(tline)
        blankdist = cell2mat(textscan(tline,'%*s %*s %f %*s'));
        
    elseif k == 21 % Number of wave samples
        disp(tline)
        Nw = cell2mat(textscan(tline,'%*s %*s %*s %*s %*s %f %*s'));

    elseif k == 22 % Frequency of data collection (burst mode)
        disp(tline)
        freq = cell2mat(textscan(tline,'%*s %*s %*s %*s %f %*s'));
        
    end
end

fclose(fid);

% Vector with the position of each cell center measured from the instrument
% Remember, first cell center is located at 'blankdist + 1 cellsize'
cells=(blankdist+(1:Ncells)*cellsize)';

fprintf('--------------------------------------------------------------\n')

fprintf('Deployment information retrieved...\n')

%% Currents data
fprintf('2. CURRENTS DATA...\n')

% z_p: upward distance from the instrument) which are averaged to obtain
%   velocity values.
curr.z_p = cells;

% Velocities
curr.u = v1;
curr.v = v2;
curr.w = v3;

% Echo intensity
curr.RSSI1 = a1;
curr.RSSI2 = a2;
curr.RSSI3 = a3;

% Time
yearsc = sen(:,3);
monthsc = sen(:,1);
daysc = sen(:,2);
hoursc = sen(:,4);
minsc = sen(:,5);
secsc = sen(:,6);
curr.t = datenum(yearsc,monthsc,daysc,hoursc,minsc,secsc);
curr.tY = curr.t-datenum(yearsc(1),0,0); % Time in yeardays

% Pressure
curr.P = sen(:,14);

% Temperature
curr.T = sen(:,15);

% Heading
curr.heading = sen(:,11);

% Pitch
curr.pitch = sen(:,12);

% Roll
curr.roll = sen(:,13);

% Voltage
curr.volt = sen(:,9);

% Message
fprintf('Currents data retrieved\n')

%% Wave data
% Message
fprintf('3. WAVES DATA...\n')

% Time (complete wave data)
% Assuming that the .wad file includes the time data
yearswad = wad(:,3);
monthswad = wad(:,1);
dayswad = wad(:,2);
hourswad = wad(:,4);
minswad = wad(:,5);
secswad = wad(:,6);
twad = datenum(yearswad,monthswad,dayswad,hourswad,minswad,secswad); % Initial time for each burst

% Time (initial datum for each burst; .whd 'always' has the time in it, but
%   DOUBLE CHECK THOSE TEXT FILES!)
yearswhd = whd(:,3);
monthswhd = whd(:,1);
dayswhd = whd(:,2);
hourswhd = whd(:,4);
minswhd = whd(:,5);
secswhd = whd(:,6);
twhd = (datenum(yearswhd,monthswhd,dayswhd,hourswhd,minswhd,secswhd))'; % Initial time for each burst

% Data collection parameters
Nb = length(whd(:,1)); % Number of bursts

twvec = (0:1:Nw-1)/freq/(3600*24); % Time vector for each burst (in days)

% Preallocation
wave.t = nan(Nw,Nb); % Time for each wave datum
wave.P = nan(Nw,Nb); % Pressure recorded at bursts
wave.STrk = nan(Nw,Nb); % Acoustic Surface Track (AST) distance recorded at bursts
wave.u = nan(Nw,Nb); % Velocity (East)
wave.v = nan(Nw,Nb); % Velocity (North)
wave.w = nan(Nw,Nb); % Velocity (Up)
wave.RSSI1 = nan(Nw,Nb); % RSSI Amplitude (Beam 1)
wave.RSSI2 = nan(Nw,Nb); % RSSI Amplitude (Beam 2)
wave.RSSI3 = nan(Nw,Nb); % RSSI Amplitude (Beam 3)

% Temperature
wave.T = whd(:,17)';

% Heading
wave.heading = whd(:,12)';

% Pitch
wave.pitch = whd(:,13)';

% Roll
wave.roll = whd(:,14)';

% Voltage
wave.volt = whd(:,10)';


% Indices for each burst may be located in first column of within .wad
%   file; otherwise, there will be values of time (month, day, year,...)
% If time was recorded in both .wad and .whd files, first computed time
% value (in vectors twad and twhd) must be the same; on the contrary, they
%   will be VERY different.
fprintf('Indices for wave structure array...\n')
if twad(1) == twhd(1) % If time was recorded in .wad file...
    
    for i = 1:Nb
        % Counting progress...
        if i == round(Nb/4) || i == round(Nb/50) || ...
                i == round(Nb/75) || i == Nb
            fprintf(['Creating wave structure array indices: ' ...
                num2str(round(i/Nb*100)) '%% \n'])
        end
        
        % Find time in .wad file that is greater than the initial time of
        %   burst 'i' that is given by 'twhd(i)' value.
        % For the last burst (assuming last burst was not completed):
        if i == Nb
            indexburst{i} = find(twad>=twhd(i));
        else
            
            % For every burst except the last one:
            indexburst{i} = find(twad>=twhd(i)&twad<twhd(i+1));
        end
    end
    
else % If time was NOT recorded in .wad file...
    for i = 1:Nb
        % Counting progress...
        if i == round(Nb/4) || i == round(Nb/2) || ...
                i == round(3*Nb/4) || i == Nb
            fprintf(['Creating wave structure array indices: ' ...
                num2str(round(i/Nb*100)) '%% \n'])
        end
        
        % If time was NOT recorded in .wad file, each burst is marked in column
        %   1 within .wad file. So, look for 'i' values of this column of .wad
        %   and create the indices for each burst.        
        indexburst{i} = find(wad(:,1)==i);
    end
end

fprintf('Wave structure array...\n')
% Storage of wave data using indices for each burst in 'indexburst'
if twad(1) == twhd(1) % .wad and .whd files include the time data
    for i = 1:Nb
        % Counting...
        if i == round(Nb/4) || i == round(Nb/2) || ...
                i == round(3*Nb/4) || i == Nb
            fprintf(['Creating wave structure array: ' ...
                num2str(round(i/Nb*100)) '%% \n'])
        end
        
        % Extract the variables associated with the index (concatenate NaNs if
        %   there is an incomplete burst)
        wave.t(:,i) = [twad(indexburst{i});nan(Nw-length(indexburst{i}),1)]; % Time in Matlab days
        wave.P(:,i) = [wad(indexburst{i},7);nan(Nw-length(indexburst{i}),1)]; % Pressure (dbar)
        wave.STrk(:,i) = ...
            ([wad(indexburst{i},8);nan(Nw-length(indexburst{i}),1)]+...
            [wad(indexburst{i},9);nan(Nw-length(indexburst{i}),1)])/2; % Acoustic Surface Tracking distance, m
        wave.u(:,i) = [wad(indexburst{i},12);nan(Nw-length(indexburst{i}),1)]; % u (East), m/s
        wave.v(:,i) = [wad(indexburst{i},13);nan(Nw-length(indexburst{i}),1)]; % v (North), m/s
        wave.w(:,i) = [wad(indexburst{i},14);nan(Nw-length(indexburst{i}),1)]; % w (Up), m/s
        wave.RSSI1(:,i) = [wad(indexburst{i},15);nan(Nw-length(indexburst{i}),1)]; % echo amplitude (East), in counts
        wave.RSSI2(:,i) = [wad(indexburst{i},16);nan(Nw-length(indexburst{i}),1)]; % echo amplitude (North), in counts
        wave.RSSI3(:,i) = [wad(indexburst{i},17);nan(Nw-length(indexburst{i}),1)]; % echo amplitude (Up), in counts
    end
    
else % .wad file DOES NOT include the time for each datum of wave data
    
    for i = 1:Nb
        % Counting...
        if i == round(Nb/4) || i == round(Nb/2) || ...
                i == round(3*Nb/4) || i == Nb
            fprintf(['Creating wave structure array: ' ...
                num2str(round(i/Nb*100)) '%% \n'])
        end
        
        % Extract the variables associated with the index (concatenate NaNs
        %   if there is an incomplete burst)
        wave.t(:,i) = twhd(1,i)+twvec'; % Time in Matlab days
        wave.P(:,i) = [wad(indexburst{i},3);nan(Nw-length(indexburst{i}),1)]; % Pressure (dbar)
        wave.STrk(:,i) = ...
            ([wad(indexburst{i},4);nan(Nw-length(indexburst{i}),1)]+...
            [wad(indexburst{i},5);nan(Nw-length(indexburst{i}),1)])/2; % Acoustic Surface Tracking distance, m
        wave.u(:,i) = [wad(indexburst{i},8);nan(Nw-length(indexburst{i}),1)]; % u (East), m/s
        wave.v(:,i) = [wad(indexburst{i},9);nan(Nw-length(indexburst{i}),1)]; % v (North), m/s
        wave.w(:,i) = [wad(indexburst{i},10);nan(Nw-length(indexburst{i}),1)]; % w (Up), m/s
        wave.RSSI1(:,i) = [wad(indexburst{i},11);nan(Nw-length(indexburst{i}),1)]; % echo amplitude (Beam 1), in counts
        wave.RSSI2(:,i) = [wad(indexburst{i},12);nan(Nw-length(indexburst{i}),1)]; % echo amplitude (Beam 2), in counts
        wave.RSSI3(:,i) = [wad(indexburst{i},13);nan(Nw-length(indexburst{i}),1)]; %  echo amplitude (Beam 3), in counts
    end
end

wave.tY = wave.t-datenum(min(yearswhd),0,0); % Time in Yeardays

fprintf('Waves data retrieved\n')

%% Save file with raw data
save([pathname filename '-raw.mat'],'curr','wave','coord');

% Show messages about finished job and elapsed time
endt = toc;

fprintf(['Raw data from ' ins ' loaded\n'])
fprintf([filename '-raw.mat file created\n'])
fprintf(['Elapsed time: ' num2str(endt/60) ' min\n'])
end