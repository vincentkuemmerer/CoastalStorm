function [icStorm,icCalm,info] = StormIdentification(Hs,TempRes,ST,ID,IT,MSD)

% -------------------------------------------------------------------------
%
% This function identifies storms from a significant wave
% height (Hs) timeseries. The following parameters must be given to the
% function (example values based on the Outer Hebrides):
%
% - Hs: The Hs timeseries [double,1*n]
% 
% - TempRes: The temporal resolution is the time interval in between 
%                observations (e.g. 1 hour) [double,1], in hours
%
% - ST: The storm threshold ST determines if wave conditions are 
%       considered extreme (e.g. 6.31 meter - 95th percentile). [double,1]
%
% - ID: The independence duration ID is the minimum time between two 
%       independent consecutive storms
%       (e.g. 48 hour - extremal index). [double,1]
%
% - IT: The independence threshold (IT) is another Hs treshold that
%       represents storm dissipation and must be smaller than ST but 
%       above the mean Hs for the coastal region 
%       (e.g. 3.8 meter - mean winter Hs) [double,1]
%
% - MSD: The minimum storm duration (MSD) is the minimum time of a storm
%        to cause relevant coastal impacts. Note that this value is
%        automatically adjusted in case r is not 1 hour 
%        (e.g. 6 hour) [double,1]
%
% The function creates 3 outputs:
% 
% - icStorm: The storm value indices
% - icCalm: The calm period indices
% - info: Cell array containing information about the effect of storm 
%         criteria on the number of extracted storms    
%
% -------------------------------------------------------------------------

% ------ Variable initialization -------

if nargin<7 % Check if MSD is provided
  MSD = 6;
end

if TempRes > 1 % Adjust MSD and ID in case TemporalRes > 1
    MSD = MSD / TempRes;
    ID = round(ID / TempRes);
end

if ID < 1 % ID cannot be smaller than 1
    ID = 1;
end

zE = 0;
zS = 0;

W = logical(Hs>ST); % Logical parameter of ST exceedance

info = cell(4,2); % Cell array with information
info{1,1} = "Number of storms using ID";
info{2,1} = "Storms not considered because of MSD";
info{3,1} = "Extra storms created using IT";
info{4,1} = "Number of storms using all criteria";

% ----- Get indices of Storm Start and Storm End (EoS_ind,SoS_ind) ------

for k = ID+1:length(W)-ID % Loop over timeseries with ID
    if (W(k-1) == 1) && (W(k) == 0) % Storm End index
        if sum(W(k:k+ID-1)) == 0
            zE = zE+1;
            EoS_ind(zE) = k-1; % Last value with Hs > ST
        end  
    end
    if (W(k-1) == 0) && (W(k) == 1) % Storm Start index
        if sum(W(k-ID:k-1)) == 0
            zS = zS+1;
            SoS_ind(zS) = k; % First value with Hs > ST
        end
    end    
end
EoS_ind = nonzeros(EoS_ind)'; % remove zeros
SoS_ind = nonzeros(SoS_ind)';

% ------ Check if storms occur at the start or end of the ts -----
% (ID is restricting the storm identification in this period)

SoS_start = 0;
SoS_end = 0;
EoS_start = 0;
EoS_end = 0;

if any(W(1:ID) == 1) % Check for storms at the start of ts 
    SoS_start = find((W(1:ID) == 1),1);
    EoS_start = SoS_start;
    while (W(EoS_start) == 1) && (EoS_start ~= length(W))...
            && (sum(W(k:k+ID-1)) == 0) % replaced k with EoS_start
        EoS_start = EoS_start+1;
    end
    EoS_start = EoS_start-1;
    if EoS_start > ID
        EoS_start = 0;
    end
end 

if any(W((length(W)-ID):end) == 1) % Check for storms at the end of ts
    SoS_end = length(W)-ID+find((W((length(W)-ID):end) == 1),1)-1;
    EoS_end = SoS_end;
    while (W(EoS_end) == 1) && (EoS_end ~= length(W))
        EoS_end = EoS_end+1;
    end
    EoS_end = EoS_end-1;
    if EoS_end == length(W)-1
        EoS_end = EoS_end+1;
    end
    if sum(W((length(W)-2*ID):(length(W)-ID))) == 0
        SoS_end = length(W)-ID+find((W((length(W)-ID):end) == 1),1)-1;
    else
        SoS_end = 0;
    end
end

% Add to the EoS_ind and SoS_ind in case storms occured
if (EoS_start ~= 0) && (EoS_start == EoS_end)
    EoS_ind = EoS_end;
    SoS_ind = SoS_start;
    SoS_start = 0;
    SoS_end = 0;
end
if SoS_start ~= 0
    SoS_ind = horzcat(SoS_start,SoS_ind);
end
if EoS_start ~= 0
    EoS_ind = horzcat(EoS_start,EoS_ind);
end
if SoS_end ~= 0
    SoS_ind = horzcat(SoS_ind,SoS_end);
end
if EoS_end ~= 0
    EoS_ind = horzcat(EoS_ind,EoS_end);
end
EoS_ind = unique(EoS_ind); % remove duplicates
SoS_ind = unique(SoS_ind);
EoS_ind = nonzeros(EoS_ind)'; % remove zeros
SoS_ind = nonzeros(SoS_ind)';

% ----- Pair start and end Storm values (icStorm) ------

icS = cell(1,length(SoS_ind)); % Preallocation
icStorm = cell(1,length(SoS_ind));
StormHs_check = zeros(1,length(SoS_ind));
nIT_c = zeros(1,length(SoS_ind));
m = 1;

for k = 1:(length(SoS_ind)) % Loop over number of storm ends
    icS{k} = SoS_ind(k):EoS_ind(k); % Get all storm indices
        
    % ----- Apply IT -----
    StormHs = Hs(icS{k})';
    StormHs_check(k) = any((StormHs <= IT)); % Check for storm values <= IT
    
    if StormHs_check(k) == true % Check if any Storm Hs <= IT
        % Check how many IT interruptions occur
        nIT = 20; % Possible storm interruptions  - arbitrarily 
        for j = 1:nIT % Loop over number of IT interruptions
            if j == 1
                StormHs2 = Hs(icS{k}(1:end))';
            else
                StormHs2 = Hs(icStorm{m-1})';
            end
            % Cut storm in 2
            FirstIdx = find((StormHs2 <= IT),1); % Find the first value < IT
            if isempty(FirstIdx) == true % If FirstIdx is empty it continues with the next iteration
                continue
            end
            nIT_c(k) = nIT_c(k) + 1;
            zE1(j) = FirstIdx; % Assign index
            zS2(j) = FirstIdx; % Assign index
            while StormHs2(zE1(j)) <= ST % Get the end of the 1st storm
                zE1(j) = zE1(j) - 1;
            end
            while StormHs2(zS2(j)) <= ST % Get the start of the 2nd storm
                zS2(j) = zS2(j) + 1;
            end
            % Create new storm index varibale 
            if j>1
                iSvals2_nz = icStorm{m-1};
                icStorm{m-1} = iSvals2_nz(1:zE1(j)); % Define boundaries of the 1nd storm
                icStorm{m} = iSvals2_nz(zS2(j):end); % Define boundaries of the 2st storm
                m = m+1;
            else
                icStorm{m} = icS{k}(1:zE1(j)); % Define boundaries of the 1nd storm
                icStorm{m+1} = icS{k}(zS2(j):end); % Define boundaries of the 2st storm
                m = m+2;
            end
        end
    else
        icStorm{m} = icS{k}(1:end);
        m = m+1;
    end
end
info{1,2} = length(icS);
info{3,2} = sum(nIT_c);

% ----- Apply Minimum Storm Duration (MSD) correction -----

MSD_count = 0;
for k = length(icStorm):-1:1
    if length(icStorm{k}) < MSD
        icStorm(k) = [];
        MSD_count = MSD_count+1;
    end
end 
info{2,2} = MSD_count;
info{4,2} = length(icStorm);

% ----- Get storm value indices (icStorm) -----

% Get new start and end values of storms
SoS_ind2 = zeros(1,length(icStorm));
EoS_ind2 = zeros(1,length(icStorm));

for k = 1:length(icStorm)
    SoS_ind2(k) = icStorm{k}(1);
    EoS_ind2(k) = icStorm{k}(end);
end

% ----- Get calm period indices (icCalm) -----

icCalm = cell(1,length(icStorm)-1);

FirstiNSvals = 1:SoS_ind2(1);
LastiNSvals = (EoS_ind2(end)+1):length(Hs);
for k = 1:(length(icStorm)-1)
    icCalm{k} = (EoS_ind2(k)+1):(SoS_ind2(k+1)-1);
end
icCalm = horzcat(FirstiNSvals,icCalm,LastiNSvals);

end
