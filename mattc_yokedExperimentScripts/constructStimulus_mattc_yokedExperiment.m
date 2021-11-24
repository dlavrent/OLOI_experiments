function [stimTimes, stim, duration] = constructStimulus_mattc_yokedExperiment(chargeTime,fname, flyNumber)
% function [stimTimes, stim, duration] = constructStimulus_mattc_yokedExperiment
% fname is the name of the processed flyTracks file whose behavioral-based
% odor experience you want to deliver
% flynumber is the specific fly whose behavioral odor experience you want
% to deliver


if nargin < 1
    chargeTime = 0;
end

load(fname, 'flyTracks');                   % load flyTracks data structure
firstOdorFrame = find(floor(flyTracks.relTimes) == (flyTracks.stim{2}(1) + flyTracks.chargeTime),1);
lastOdorFrame = find(ceil(flyTracks.relTimes) == flyTracks.stim{2}(end),1,'last');

% create position vector
% assuming that MCH is the reference odor
a = (flyTracks.headLocal(firstOdorFrame:lastOdorFrame,2,flyNumber) > flyTracks.corridorPos(end));
b = (flyTracks.headLocal(firstOdorFrame:lastOdorFrame,2,flyNumber) < flyTracks.corridorPos(1));

% find fraction of tunnel occupied as a function of odor period
% -1 = MCH, 0 = corridor (no odor), 1 = OCT
tunnelPosition=a-b; % can be either -1, 0, or 1, which determines which solenoid is open

newOdorFrames=[1; find(diff(tunnelPosition))+1];

odorBlockTimes=[flyTracks.relTimes(newOdorFrames)];
newOdorFrames(odorBlockTimes>(flyTracks.relTimes(lastOdorFrame)-flyTracks.relTimes(firstOdorFrame)))=[];
odorBlockTimes(odorBlockTimes>(flyTracks.relTimes(lastOdorFrame)-flyTracks.relTimes(firstOdorFrame)))=[];

nBlocks = length(odorBlockTimes);   % number of odor blocks
    
odorDurations=diff([odorBlockTimes flyTracks.relTimes(lastOdorFrame)-flyTracks.relTimes(firstOdorFrame)]); % seconds

% tooShort=odorDurations<1;
% newOdorFrames(tooShort)=[];
% odorDurations(tooShort)=[];
% nBlocks=nBlocks-sum(tooShort);

% remove odor epochs shorter than a threshold
%odorBinThreshold=3; % seconds


odors = {'MCH' 'OCT' 'Air'};
%conc = [0.1 0.15];%  default -aMW ctrl               % proportion saturated vapor         
conc = [0.10 0.16];

OdorAonRight = false;                                                        % KH141217 - hardcode this for consistency

% if odor A is on right, need to make sure concentration swaps as well
if OdorAonRight == true
   concT=conc;
   conc(2)=concT(1);
   conc(1)=concT(2);
end

odorDur = 180;                        % in sec
isi = 0;                            % in sec

preTime = 185;                       % wait time before first odor block
postTime = 10;                      % wait time after last odor block


% Read valve assignments from csv file
fid = fopen('C:\Users\khonegger\Documents\MATLAB\TunnelSoftware\Acquisition\odors.csv');
v = textscan(fid, '%s %d %d', 'delimiter', ','); % Format: {Odor, SideA, SideB}
fclose(fid);

% Alternate sides, with random starting side
%OdorAonTop = repmat(randperm(2)-1, [1, ceil(nBlocks/2)]); % Logical vector
% indicating
% blocks with
% odors{1} on top
% (Side B)


%Odor matrix: [top/bottom concentration; top/bottom valves]
for qq = 1:nBlocks
   
    if tunnelPosition(newOdorFrames(qq))==(1)
        valve(1) = v{3}(strmatch(odors{2}, v{1}));
        valve(2) = v{2}(strmatch(odors{2}, v{1}));
        stim(qq).odor = [conc(2), conc(2); double(valve)];
    elseif tunnelPosition(newOdorFrames(qq))==(-1)
        valve(1) = v{2}(strmatch(odors{1}, v{1}));
        valve(2) = v{3}(strmatch(odors{1}, v{1}));
        stim(qq).odor = [conc(1), conc(1); double(valve)];
    else
        valve(1) = v{2}(strmatch(odors{3}, v{1}));
        valve(2) = v{3}(strmatch(odors{3}, v{1}));
        stim(qq).odor = [0.0, 0.0; double(valve)];
    end
    

    
    
    % Sort odor labels [Side A, Side B]
    if sum((v{2}==valve(1)))
        stim(qq).labels = [v{1}(v{2}==valve(1)) v{1}(v{3}==valve(2))];
    else
        stim(qq).labels = [v{1}(v{2}==valve(2)) v{1}(v{3}==valve(1))];
    end
end

% 
% 
% if OdorAonRight
%     if v{1}{find(v{3} == valve(1))} == odors{1}
%         conc(1) = conc(1);
%     else
%         conc(1) = conc(2);
%     end
%     if v{1}{find(v{2} == valve(2))} == odors{1}
%         conc(2) = conc(1);
%     else
%         conc(2) = conc(2);
%     end
% else
%     if v{1}{find(v{2} == valve(1))} == odors{1}
%         conc(1) = conc(1);
%     else
%         conc(1) = conc(2);
%     end
%     if v{1}{find(v{3} == valve(2))} == odors{1}
%         conc(2) = conc(1);
%     else
%         conc(2) = conc(2);
%     end
% end


%Build stimulus epochs
% multiply times by 10 to get accuracy down to 100 ms
lastTime = 10*(preTime);                 % Runs on first block only

for i=1:nBlocks
    startTime = (lastTime);
    stim(i).times = startTime:(startTime + 10*chargeTime + 10*odorDurations(i));
    lastTime = max(stim(i).times) - 10*chargeTime;
end


duration = lastTime + 10*postTime + 10*chargeTime;

stimTimes = zeros(1,duration + 1);

for qqq=1:nBlocks
    stimTimes(stim(qqq).times) = qqq;
end