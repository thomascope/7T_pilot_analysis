function Vocode_report_PINFA_VESPASTIM_updated

% A new script written for the VESPA (Vocoder evaluation of speech
% perception in aphasia) study, by T Cope, 17/2/14.
% Presents a target, vocoded word
% After an interval of 1050ms (+/- up to 50ms), presents 4AFC
% Correct option is in ratio 2:1:1:2:
% 1) A base word
% 2) An onset phonemic switch from (1)
% 3) An offset phonemic switch from (1)
% 4) An unrelated word, drawn from a permuted list of (1), sharing no
% phonemic elements

global lisName age res_box
if isempty(lisName) ~= 1 || isempty(age) ~= 1 || isempty(res_box) ~= 1
    useprev = input('Would you like to use previous subject demographics? Please enter y or n: ','s');
    while strcmp(useprev,'y') == 0 & strcmp(useprev,'n') == 0
        useprev = input('Please enter y or n only: ','s');
    end
    if strcmp(useprev,'n')==1
        clear all
        global lisName age res_box
    end
end

Screen('Preference', 'SkipSyncTests', 1)

% simulation mode?
simulate = 0;

randn('state',sum(100*clock));
rand('state',sum(100*clock));

dateTime = clock;
dn = datenum(dateTime(1:3));
dato = [datestr(dn,11) datestr(dn,5) datestr(dn,7)];

if isempty(lisName) == 1 || isempty(age) == 1 || isempty(res_box) == 1
    lisName = input('Enter listener initials: ','s');
    age = input('Enter listener age: ');
    res_box = input('\nEnter response device, i.e. k for keyboard or r for response box: ','s');
    while strcmp(res_box,'k') == 0 & strcmp(res_box,'r') == 0
        res_box = input('Please enter k or r only: ','s');
    end
    if strcmp(res_box,'r')==1
        res_box = 1;
    else
        res_box = 0;
    end    
end
% feedback = input('Would you like feedback? Please enter y or n: ','s');
% while strcmp(feedback,'y') == 0 & strcmp(feedback,'n') == 0
%     feedback = input('Please enter y or n only: ','s');
% end

preload = input('\nwould you like to preload sounds? Please enter y or n: ','s');
while strcmp(preload,'y') == 0 & strcmp(preload,'n') == 0
    preload = input('Please enter y or n only: ','s');
end
if strcmp(preload,'y')==1
    preload = 1;
else
    preload = 0;
end

feedback = 'n'; % never give feedback during this version of task - not implemented anyway yet!

% define directories
DIRS.base = pwd;
DIRS.stimuli = strcat(DIRS.base,{'\stimuli\normalised\'});
DIRS.stimuli = char(DIRS.stimuli);
DIRS.data = strcat(DIRS.base,{'\pinfa_data\'});
DIRS.data = char(DIRS.data);
DIRS.dataFigs = strcat(DIRS.base,{'\pinfa_datafigs\'});
DIRS.dataFigs = char(DIRS.dataFigs);

% define save filename
mfname = mfilename;
underscores = strfind(mfname,'_');
dfname = [mfname(1:underscores(1)-1) '_VESPA'];
dfilename = [dfname '_' num2str(age) '_' lisName '_' dato];


% load mat file containing: fixed stimulus_list, target_list, difficulty_list (number of vocode channels) 
load vocode_report_params;
NTrials = length(target_list);              % number of trials per run

% number of full sets (6 trials, 2:1:1:2) of practice
NPTrialsets = 2;
NPTrials = NPTrialsets*6;
clear_practice = 1; % use clear speech in practice?

% stimuli
srate = 44100;                              % sampling rate in Hz
sfactor = 0.2;                              % crude amplitude adjustment
wdms = 20;                                  % window in ms raised cosine

% define SOAs in ms as 1050 jittered +/- by up to 50ms
all_SOAs_ms = 1050 + 50.*rand(NTrials,1).*sign(rand(NTrials,1)-0.5);
all_iti_ms = 850 + 50.*rand(NTrials,1).*sign(rand(NTrials,1)-0.5);

% set experiment parameters
vocode_channels = [3,6,15]; % Number of vocoder channels in presented speech
target_types = [{'onset_switch'},{'base'},{'offset_switch'},{'mismatch'}];

% mat save params
save([DIRS.data dfilename],'preload','lisName','age','NTrials','NPTrials','wdms','srate','all_SOAs_ms','all_iti_ms',...
    'stimulus_list', 'target_list','res_box', 'feedback', 'vocode_channels', 'target_types');

all_practice_mixes = [];
all_mixes = [];
all_practicetrial_targets = [];
all_trial_targets = [];
all_practiceoptions = [];
all_options = [];

%generate random sequence of targets for practice trials
practice_target_list = [1,1,2,3,4,4];
practice_targets = practice_target_list(randperm(length(practice_target_list)));
if NPTrialsets>1
    for i=1:(NPTrialsets-1)
        practice_targets = [practice_targets practice_target_list(randperm(length(practice_target_list)))];
    end
end
practice_difficulty_list = [1,2,3,1,2,3];
practice_difficulties = practice_difficulty_list(randperm(length(practice_difficulty_list)));
if NPTrialsets>1
    for i=1:(NPTrialsets-1)
        practice_difficulties = [practice_difficulties practice_difficulty_list(randperm(length(practice_difficulty_list)))];
    end
end

% Randomly draw practice stimuli from experimental
if NTrials >= NPTrials
    practicestims = randperm(NTrials);
    practicestims = practicestims(1:NPTrials);
else
    error('You are trying to do more practice runs than real runs!');
end
practice_stimulus_list = stimulus_list(:,practicestims);

% Preload practice sounds if selected
if preload == 1
    disp('Loading practice sounds');
    if clear_practice == 1
        for i = 1:NPTrials
            soundfilename = [DIRS.stimuli '\Originals\' char(practice_stimulus_list(practice_targets(i),i)) '.wav'];
            [psig{i}.y,fs_check] = audioread(soundfilename);
            psig{i}.y(:,2) = psig{i}.y;
            psig{i}.y = psig{i}.y';
            if srate ~= fs_check,
                error('Sampling rate not as expected')
                
            end
        end
    else
        for i = 1:NPTrials
            soundfilename = [DIRS.stimuli char(practice_stimulus_list(practice_targets(i),i)) '_noise_n_greenwood_half_30_' num2str(vocode_channels(practice_difficulties(i))) '.wav'];
            [psig{i}.y,fs_check] = audioread(soundfilename);
            psig{i}.y(:,2) = psig{i}.y;
            psig{i}.y = psig{i}.y';
            if srate ~= fs_check,
                error('Sampling rate not as expected')
                i
            end
        end
    end
end

if preload == 1
    disp('Loading main trial sounds');
    for i = 1:NTrials
        soundfilename = [DIRS.stimuli char(stimulus_list(target_list(i),i)) '_noise_n_greenwood_half_30_' num2str(vocode_channels(difficulty_list(i))) '.wav'];
        [sig{i}.y,fs_check] = audioread(soundfilename);
        sig{i}.y(:,2) = sig{i}.y;
        sig{i}.y = sig{i}.y';
        if srate ~= fs_check,
            error('Sampling rate not as expected')
            i
        end
    end
end

% start experiment

try
    ListenChar(2);
    
    % check for opengl compatability
    AssertOpenGL;
    
    % get screen
    screens = Screen('Screens');
    screenNumber = max(screens);
    HideCursor;
   
    
    % set window
    pixdepth = 32;
    buffermode = 2; % double buffer
    [w, wRect]=Screen('OpenWindow', screenNumber, 0, [], pixdepth, buffermode);
    [width, height] = Screen('WindowSize', w);
    Screen(w,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
% 
% % For debug only
%     pixdepth = 32;
%     buffermode = 2; % double buffer
%     [w, wRect]=Screen('OpenWindow', screenNumber, 0, [100,100,500,500], pixdepth, buffermode);
%     [width, height] = Screen('WindowSize', w);
%     Screen(w,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    % set colours
    gray = GrayIndex(screenNumber);
    white = WhiteIndex(screenNumber);
    black = BlackIndex(screenNumber);
    
    % set priority
    priorityLevel = MaxPriority(w);
    Priority(priorityLevel);
    
    % clear screen
    Screen('FillRect',w, gray);
    Screen('Flip', w);

    % set text parameters
    yoffset = 0;
    textsize = 32;
    Screen(w,'TextFont','Courier New');
    Screen('TextStyle', w, [1]); %set text to bold
    
    % instructions
    if (res_box ==1)        
        drawText('You will hear a word and then see some choices' , w, yoffset - 60, white, textsize);
        drawText('Press the button that matches what you heard', w, yoffset, white, textsize);
        drawText('Press any button to begin a practice', w, yoffset + 60, white, textsize);
    else
        drawText('You will hear a word and then see some choices' , w, yoffset - 60, white, textsize);
        drawText('Press the number that matches what you heard', w, yoffset, white, textsize);
        drawText('Press the mouse to begin a practice', w, yoffset + 60, white, textsize);
    end
    Screen('Flip', w);
    
    % wait for input
    buttons = 0;
    if (res_box == 1)
        blink_wait_usb
    else
        while ~any(buttons) % wait for press
            [x,y,buttons] = GetMouse;
            % Wait 10 ms before checking the mouse again to prevent
            % overload of the machine at elevated Priority()
            WaitSecs(0.01);
        end
    end
    
    % clear screen
    Screen('Flip', w);
    
    % wait a bit before starting trial
    WaitSecs(1.000);
    

    
    % loop through practice trials
    for trial=1:(NPTrials),
        if simulate
            dat.correct(trial) = round(rand()+(roughscore-0.5));
        else
            Screen('Flip', w); %clear screen
            disp(['Started trial: ' int2str(trial)]);
            % wait a bit between trials (ITI)
            WaitSecs(all_iti_ms(trial)/1000);
            
            trialstart = GetSecs;
            dat.pT.start(trial) = 0;
            
            % initialize key
            [KeyIsDown, endrt, KeyCode] = KbCheck;
            
            % load sound for trial if not preloaded
            if preload == 0
                if clear_practice == 1
                    soundfilename = [DIRS.stimuli '\Originals\' char(practice_stimulus_list(practice_targets(trial),trial)) '.wav'];
                    [psig.y,fs_check] = audioread(soundfilename);
                    psig.y(:,2) = psig.y;
                    psig.y = psig.y';
                    if srate ~= fs_check,
                        error('Sampling rate not as expected');
                    end
                else
                    soundfilename = [DIRS.stimuli char(practice_stimulus_list(practice_targets(trial),trial)) '_noise_n_greenwood_half_30_' num2str(vocode_channels(practice_difficulties(trial))) '.wav'];
                    [psig.y,fs_check] = audioread(soundfilename);
                    psig.y(:,2) = psig.y;
                    psig.y = psig.y';
                    if srate ~= fs_check,
                        error('Sampling rate not as expected');
                    end
                end
            end
            
            %start trial
            tic
            trialstart = GetSecs;
            if preload == 1
                Snd('Play',psig{trial}.y,srate);
                Snd('Wait');
            elseif preload == 0
                Snd('Play',psig.y,srate);
                Snd('Wait');
            end
            
            options = practice_stimulus_list(:,trial);
            mixitup = randperm(size(practice_stimulus_list,1));
            options = options(mixitup);
            practicetrial_target = find(mixitup==practice_targets(trial));
            all_practice_mixes = [all_practice_mixes; mixitup];
            all_practicetrial_targets = [all_practicetrial_targets; practicetrial_target];
            all_practiceoptions = [all_practiceoptions, options]; 
            sound_time = toc;
            
            if all_SOAs_ms(trial)/1000 > sound_time
                pause(all_SOAs_ms(trial)/1000 - sound_time)
                drawnow
            else
                disp('Error, the sound took longer to play than the SOA, proceeding immediately!')
            end
            
            listoptions=['1: ', options{1},'   2: ', options{2},'   3: ', options{3},'   4: ', options{4}];
            drawText(listoptions, w, yoffset, white, textsize);
            Screen('Flip', w);
         
            startrt = GetSecs-trialstart;
            if (res_box == 1)
                response_buttons = read_usbbox;
                response = find(response_buttons == 1);
                drawnow
                endrt = GetSecs-trialstart;
            else
                while (1),
                    % quit experiment, if desired
                    if ( KeyIsDown==1 & KeyCode(27)==1 ) %was "esc" pressed?
                        Screen('FillRect', w, gray);
                        Screen('Flip', w);
                        WaitSecs(0.500);
                        Screen('CloseAll');
                        ShowCursor;
                        ListenChar;
                        Priority(0);
                        Snd('Quiet');
                        disp('Experimenter exit ...');
                        error('exit');
                    end
                    if ( KeyCode(KbName('1'))==1 | KeyCode(KbName('2'))==1 | KeyCode(KbName('3'))==1 | KeyCode(KbName('4'))==1 )
                        break;
                    end
                    [KeyIsDown, endrt, KeyCode] = KbCheck;
                    endrt = GetSecs-trialstart;
                end
            end
            
            rt = round(1000*(endrt-startrt));   % get rt
            
            if (res_box==1)
                resp={['button_' num2str(response)]};
            else
                resp = KbName(KeyCode==1);           % get key pressed
            end
            
            %record reaction times
            dat.pT.startrt(trial) = startrt;
            dat.pT.endrt(trial) = endrt;
            dat.pT.respRT_ms(trial)  = rt;
            
            %record response as a string
            dat.pT.resp_str{trial}   = resp;
            
            %record response
            if (res_box == 1)
                dat.pT.resp(trial) = response;
            else
                if KeyCode(KbName('1'))==1
                    dat.pT.resp(trial) = 1;
                elseif KeyCode(KbName('2'))==1
                    dat.pT.resp(trial) = 2;
                elseif KeyCode(KbName('3'))==1
                    dat.pT.resp(trial) = 3;
                elseif KeyCode(KbName('4'))==1
                    dat.pT.resp(trial) = 4;
                end
            end
            
            %code if response correct            
            if dat.pT.resp(trial) == practicetrial_target
                dat.pT.correct(trial) = 1;
            elseif dat.pT.resp(trial) > 4
                dat.pT.correct(trial) = -1; %miss or incorrect response
            else
                dat.pT.correct(trial) = 0;
            end
            
            dat.pT.target(trial) = practicetrial_target;

            save([DIRS.data dfilename],'dat','all_practice_mixes','all_practiceoptions','all_practicetrial_targets','all_SOAs_ms','all_iti_ms','practice_difficulties','practice_stimulus_list','practice_targets','practicestims', '-append');
            
        end
    end
    
    percentcorrect = num2str(ceil(mean(dat.pT.correct)*100));
    
    % practice feedback
    if (res_box ==1)
        practicefeedback = ['Well done, you got: ' percentcorrect ' percent correct!'];
        drawText(practicefeedback , w, yoffset - 60, white, textsize);
        drawText('Press any button to start the real experiment!', w, yoffset + 60, white, textsize);
    else
        practicefeedback = ['Well done, you got: ' percentcorrect ' percent correct!'];
        drawText(practicefeedback , w, yoffset - 60, white, textsize);
        drawText('Press the mouse to start the real experiment!', w, yoffset + 60, white, textsize);
    end
    Screen('Flip', w);
    
    % wait for input
    buttons = 0;
    if (res_box == 1)
        blink_wait_usb
    else
        while ~any(buttons) % wait for press
            [x,y,buttons] = GetMouse;
            % Wait 10 ms before checking the mouse again to prevent
            % overload of the machine at elevated Priority()
            WaitSecs(0.01);
        end
    end 
    
    % now do real experiment
    for trial=1:(NTrials),
        if simulate
            dat.correct(trial) = round(rand()+(roughscore-0.5));
        else
            Screen('Flip', w); %clear screen
            disp(['Started trial: ' int2str(trial)]);
            % wait a bit between trials (ITI)
            WaitSecs(all_iti_ms(trial)/1000);
            
            trialstart = GetSecs;
            dat.pT.start(trial) = 0;
            
            % initialize key
            [KeyIsDown, endrt, KeyCode] = KbCheck;
                     
            % load sound for trial
            if preload == 0
                soundfilename = [DIRS.stimuli char(stimulus_list(target_list(trial),trial)) '_noise_n_greenwood_half_30_' num2str(vocode_channels(difficulty_list(trial))) '.wav'];
                [sig.y,fs_check] = audioread(soundfilename);
                sig.y(:,2) = sig.y;
                sig.y = sig.y';
                if srate ~= fs_check,
                    error('Sampling rate not as expected');
                end
            end
            
            %start trial
            tic
            trialstart = GetSecs;
            if preload == 1
                Snd('Play',sig{trial}.y,srate);
                Snd('Wait');
            elseif preload == 0
                Snd('Play',sig.y,srate);
                Snd('Wait');
            end
            
            options = stimulus_list(:,trial);
            mixitup = randperm(size(stimulus_list,1));
            options = options(mixitup);
            trial_target = find(mixitup==target_list(trial));
            all_mixes = [all_mixes; mixitup];
            all_trial_targets = [all_trial_targets; trial_target];
            all_options = [all_options; options]; 
            sound_time = toc;
            
            if all_SOAs_ms(trial)/1000 > sound_time
                pause(all_SOAs_ms(trial)/1000 - sound_time)
                drawnow
            else
                disp('Error, the sound took longer to play than the SOA, proceeding immediately!')
            end
            
            listoptions=['1: ', options{1},'   2: ', options{2},'   3: ', options{3},'   4: ', options{4}];
            drawText(listoptions, w, yoffset, white, textsize);
            Screen('Flip', w);
            
            startrt = GetSecs-trialstart;
            if (res_box == 1)
                response_buttons = read_usbbox;
                response = find(response_buttons == 1);
                drawnow
                endrt = GetSecs-trialstart;
            else
                while (1),
                    % quit experiment, if desired
                    if ( KeyIsDown==1 & KeyCode(27)==1 ) %was "esc" pressed?
                        Screen('FillRect', w, gray);
                        Screen('Flip', w);
                        WaitSecs(0.500);
                        Screen('CloseAll');
                        ShowCursor;
                        ListenChar;
                        Priority(0);
                        Snd('Quiet');
                        disp('Experimenter exit ...');
                        error('exit');
                    end
                    if ( KeyCode(KbName('1'))==1 | KeyCode(KbName('2'))==1 | KeyCode(KbName('3'))==1 | KeyCode(KbName('4'))==1 )
                        break;
                    end
                    [KeyIsDown, endrt, KeyCode] = KbCheck;
                    endrt = GetSecs-trialstart;
                end
            end
            
            rt = round(1000*(endrt-startrt));   % get rt
            
            if (res_box==1)
                resp={['button_' num2str(response)]};
            else
                resp = KbName(KeyCode==1);           % get key pressed
            end
            
            %record reaction times
            dat.T.startrt(trial) = startrt;
            dat.T.endrt(trial) = endrt;
            dat.T.respRT_ms(trial)  = rt;
            
            %record response as a string
            dat.T.resp_str{trial}   = resp;
            
            %record response
            if (res_box == 1)
                dat.T.resp(trial) = response;
            else
                if KeyCode(KbName('1'))==1
                    dat.T.resp(trial) = 1;
                elseif KeyCode(KbName('2'))==1
                    dat.T.resp(trial) = 2;
                elseif KeyCode(KbName('3'))==1
                    dat.T.resp(trial) = 3;
                elseif KeyCode(KbName('4'))==1
                    dat.T.resp(trial) = 4;
                end
            end
            
            %code if response correct            
            if dat.T.resp(trial) == trial_target
                dat.T.correct(trial) = 1;
            elseif dat.T.resp(trial) > 4
                dat.T.correct(trial) = -1; %miss or incorrect response
            else
                dat.T.correct(trial) = 0;
            end
            
            save([DIRS.data dfilename],'dat','all_mixes','all_options','all_trial_targets','all_SOAs_ms','all_iti_ms','difficulty_list','stimulus_list','target_list', '-append');
            
        end
        
    end
    % cleanup at end of experiment
    disp('Normal exit')
    Screen('CloseAll');
    Snd('Quiet');
    ShowCursor;
    fclose('all');
    ListenChar;
    Priority(0);
catch ME
    % catch error
    disp('Catch error');
    Screen('CloseAll');
    Snd('Quiet');
    ShowCursor;
    fclose('all');
    ListenChar;
    Priority(0);
    ME
    save(['lasterror'], 'ME');
end
                





%%%%%%%%%%% draw text function %%%%%%%%%%%
function drawText(text, win, yoffset, color, textsize)

[width, height] = Screen('WindowSize', win);
Screen('TextSize', win, textsize);
[textbox, textbox2] = Screen('TextBounds', win, text);
Screen('DrawText', win, text, (width-textbox(3))/2, (height-textbox(4))/2 + yoffset, color);

return;

clear all
            