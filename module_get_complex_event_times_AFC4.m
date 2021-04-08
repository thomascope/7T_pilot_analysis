function tempDesign = module_get_complex_event_times_AFC4(subj_initials,testing_date,nRuns,nVolumes)


%clear all
[workingdir,file,ext] = fileparts(which('module_get_event_times.m'));
behaviour_folder = [workingdir '/behavioural_data'];

addpath(genpath('/imaging/mlr/users/tc02/toolboxes'));

%% control parameters

nSessions=1;
nWords  = 16;  % there were 16 different words presented in the experiment
nTrialtypes = 8; % there were 6 conditions - Match low, Match high, Mismatch low, Mismatch high, Writtenonly, Response. For legacy purposes 2 extra space provided for neutral conditions, one of which is now used to model the written word and one the button press.
tr           = 2.5;   % the EPI volumes were acquired every 2.5 seconds
stimdur      = 0.5;   % the audio and visual stimuli lasted approx 0.5 seconds
stimdelay    = 1.7; % the audio stimulus was presented 1700ms after the scanner pulse
cuedelay = 1; % the written stimulus was presented 1000ms after the scanner pulse
responsedelay = 1.05; % the response cue was presented 1050ms after the spoken cue
condition_order = {'Match low','Match high','Mismatch low','Mismatch high','Neutral low','Neutral high','Writtenonly','Response'}; %NB: In the new AFC there were no neutral trials, so these are zeros

% pre allocate cell arrays
tempOnsets   = cell(1,nRuns);
design       = cell(1,nRuns);
tpatterns    = cell(1,nSessions);

% the idendity of the subject will be subject initials
thisSubject = subj_initials;

% prepare design matrices
tempOnsets = cell(1,nTrialtypes*nWords+2);
tempDesign = cell(1,nRuns);

for runI=1:nRuns
    
    cd(behaviour_folder)
    this_file = dir(['AFC*' subj_initials '_Run_' num2str(runI) '_' testing_date '.mat']);
    cd(workingdir)
    
    fileName = this_file.name;
    %read the scanning log file
 
    [startpulses,stimType,stimNumber,stimName,writtenNumber,all_rts] = extract_pulsenumbers_from_AFC4_paradigm_wordidentity([behaviour_folder '/' fileName],runI);
    
    starttime = ((startpulses-1) * tr) + stimdelay; % Remember first pulse occurs at time zero
    writtentime = ((startpulses-1) * tr) + cuedelay;
    
    if size(stimNumber,2) ~= size(starttime,2)
        warning('Something has gone wrong in recording the stimulus types. Padding with NaNs but you MUST check your data')
        stimNumber = [stimNumber NaN(1,(size(starttime,2)-size(stimNumber,2)))];
    end
    
    for condJ = 1:nTrialtypes
        for condI=1:nWords
            tempOnsets{condI+((condJ-1)*nWords)}=starttime(stimNumber==condI & stimType==condJ);
            if condJ == 5 %Written only onsets
                tempOnsets{condI+((condJ-1)*nWords)}=writtentime(writtenNumber==condI);
            elseif condJ == 6 && condI == 1 %Model button presses
                tempOnsets{condI+((condJ-1)*nWords)}=starttime(stimType==8)+1.05+(all_rts(all_rts~=0)/1000);
            end
        end
    end
    
    tempOnsets{nTrialtypes*nWords+1}=starttime(stimNumber==0 & stimType==7);
    tempOnsets{nTrialtypes*nWords+2}=starttime(stimNumber==0 & stimType==8);
    
    
    
    
    tempDesign{runI} = tempOnsets;
end
