function [] = CreateRegressors(GLMAnalDir,BehavDataDir,PhysioDataFile,MotionDataDir,SubjID)

NumSess = length(BehavDataDir);

for sess=1:NumSess
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create Stimulus onset regressors
    
    DataFile = dir(sprintf('%s%s_run%d_*.txt',BehavDataDir{sess},SubjID,sess));
    DataFile = fullfile(BehavDataDir{sess},DataFile(1).name);
    fid = fopen(DataFile);
    tmp = textscan(fid,'%f%f%s%f%f%f%f%f%f%f%f%f%f%f%f%s%f%f%f%f%s%s%f','Delimiter','\t','HeaderLines',1);
    
    % Make event times relative to time of first scan
    t0 = tmp{7}(1);
    tmp{6} = tmp{6}-t0; % Sound onset
    tmp{7} = tmp{7}-t0; % Volume onset
    tmp{10} = tmp{10}-t0; % Button press
    
    counter = 0;
    
    for itemset=1:8
        counter = counter + 1;
        conIdx{counter} = find( ismember(tmp{:,21},'clearSpeech')==1 & tmp{:,18}==2 & tmp{:,20}==1 & tmp{:,17}==itemset & (tmp{:,19}==1 | tmp{:,19}==2) & tmp{:,5}==0 ); % Non-targets
        names{counter} = sprintf('Strong+M_Set%d_Item%d',itemset,1);
        onsets{counter} = double(tmp{6}(conIdx{counter}));
        durations{counter} = zeros(length(onsets{counter}),1);
        
        counter = counter + 1;
        conIdx{counter} = find( ismember(tmp{:,21},'clearSpeech')==1 & tmp{:,18}==2 & tmp{:,20}==1 & tmp{:,17}==itemset & (tmp{:,19}==3 | tmp{:,19}==4) & tmp{:,5}==0 ); % Non-targets
        names{counter} = sprintf('Strong+M_Set%d_Item%d',itemset,3);
        onsets{counter} = double(tmp{6}(conIdx{counter}));
        durations{counter} = zeros(length(onsets{counter}),1);
    end
    for itemset=1:8
        counter = counter + 1;
        conIdx{counter} = find( ismember(tmp{:,21},'clearSpeech')==1 & tmp{:,18}==1 & tmp{:,20}==1 & tmp{:,17}==itemset & (tmp{:,19}==1 | tmp{:,19}==2) & tmp{:,5}==0 ); % Non-targets
        names{counter} = sprintf('Weak+M_Set%d_Item%d',itemset,1);
        onsets{counter} = double(tmp{6}(conIdx{counter}));
        durations{counter} = zeros(length(onsets{counter}),1);
        
        counter = counter + 1;
        conIdx{counter} = find( ismember(tmp{:,21},'clearSpeech')==1 & tmp{:,18}==1 & tmp{:,20}==1 & tmp{:,17}==itemset & (tmp{:,19}==3 | tmp{:,19}==4) & tmp{:,5}==0 ); % Non-targets
        names{counter} = sprintf('Weak+M_Set%d_Item%d',itemset,3);
        onsets{counter} = double(tmp{6}(conIdx{counter}));
        durations{counter} = zeros(length(onsets{counter}),1);
    end
    for itemset=1:8
        counter = counter + 1;
        conIdx{counter} = find( ismember(tmp{:,21},'clearSpeech')==1 & tmp{:,18}==2 & tmp{:,20}==2 & tmp{:,17}==itemset & (tmp{:,19}==1 | tmp{:,19}==2) & tmp{:,5}==0 ); % Non-targets
        names{counter} = sprintf('Strong+MM_Set%d_Item%d',itemset,1);
        onsets{counter} = double(tmp{6}(conIdx{counter}));
        durations{counter} = zeros(length(onsets{counter}),1);
        
        counter = counter + 1;
        conIdx{counter} = find( ismember(tmp{:,21},'clearSpeech')==1 & tmp{:,18}==2 & tmp{:,20}==2 & tmp{:,17}==itemset & (tmp{:,19}==3 | tmp{:,19}==4) & tmp{:,5}==0 ); % Non-targets
        names{counter} = sprintf('Strong+MM_Set%d_Item%d',itemset,3);
        onsets{counter} = double(tmp{6}(conIdx{counter}));
        durations{counter} = zeros(length(onsets{counter}),1);
    end
    for itemset=1:8
        counter = counter + 1;
        conIdx{counter} = find( ismember(tmp{:,21},'clearSpeech')==1 & tmp{:,18}==1 & tmp{:,20}==2 & tmp{:,17}==itemset & (tmp{:,19}==1 | tmp{:,19}==2) & tmp{:,5}==0 ); % Non-targets
        names{counter} = sprintf('Weak+MM_Set%d_Item%d',itemset,1);
        onsets{counter} = double(tmp{6}(conIdx{counter}));
        durations{counter} = zeros(length(onsets{counter}),1);
        
        counter = counter + 1;
        conIdx{counter} = find( ismember(tmp{:,21},'clearSpeech')==1 & tmp{:,18}==1 & tmp{:,20}==2 & tmp{:,17}==itemset & (tmp{:,19}==3 | tmp{:,19}==4) & tmp{:,5}==0 ); % Non-targets
        names{counter} = sprintf('Weak+MM_Set%d_Item%d',itemset,3);
        onsets{counter} = double(tmp{6}(conIdx{counter}));
        durations{counter} = zeros(length(onsets{counter}),1);
    end
%     for itemset=1:8
%         counter = counter + 1;
%         conIdx{counter} = find( ismember(tmp{:,21},'noiseFinal')==1 & tmp{:,18}==2 & tmp{:,17}==itemset & (tmp{:,19}==1 | tmp{:,19}==2) & tmp{:,5}==0 ); % Non-targets
%         names{counter} = sprintf('Strong+Noise_Set%d_Item%d',itemset,1);
%         onsets{counter} = double(tmp{6}(conIdx{counter}));
%         durations{counter} = zeros(length(onsets{counter}),1);
%         
%         counter = counter + 1;
%         conIdx{counter} = find( ismember(tmp{:,21},'noiseFinal')==1 & tmp{:,18}==2 & tmp{:,17}==itemset & (tmp{:,19}==3 | tmp{:,19}==4) & tmp{:,5}==0 ); % Non-targets
%         names{counter} = sprintf('Strong+Noise_Set%d_Item%d',itemset,3);
%         onsets{counter} = double(tmp{6}(conIdx{counter}));
%         durations{counter} = zeros(length(onsets{counter}),1);
%     end
%     for itemset=1:8
%         counter = counter + 1;
%         conIdx{counter} = find( ismember(tmp{:,21},'noiseFinal')==1 & tmp{:,18}==1 & tmp{:,17}==itemset & (tmp{:,19}==1 | tmp{:,19}==2) & tmp{:,5}==0 ); % Non-targets
%         names{counter} = sprintf('Weak+Noise_Set%d_Item%d',itemset,1);
%         onsets{counter} = double(tmp{6}(conIdx{counter}));
%         durations{counter} = zeros(length(onsets{counter}),1);
%         
%         counter = counter + 1;
%         conIdx{counter} = find( ismember(tmp{:,21},'noiseFinal')==1 & tmp{:,18}==1 & tmp{:,17}==itemset & (tmp{:,19}==3 | tmp{:,19}==4) & tmp{:,5}==0 ); % Non-targets
%         names{counter} = sprintf('Weak+Noise_Set%d_Item%d',itemset,3);
%         onsets{counter} = double(tmp{6}(conIdx{counter}));
%         durations{counter} = zeros(length(onsets{counter}),1);
%     end
    for itemset=1:8
        counter = counter + 1;
        conIdx{counter} = find( ismember(tmp{:,21},'noiseInitial')==1 & tmp{:,17}==itemset & (tmp{:,19}==1 | tmp{:,19}==2) & tmp{:,5}==0 ); % Non-targets
        names{counter} = sprintf('Noise+Speech_Set%d_Item%d',itemset,1);
        onsets{counter} = double(tmp{6}(conIdx{counter}));
        durations{counter} = zeros(length(onsets{counter}),1);
        
        counter = counter + 1;
        conIdx{counter} = find( ismember(tmp{:,21},'noiseInitial')==1 & tmp{:,17}==itemset & (tmp{:,19}==3 | tmp{:,19}==4) & tmp{:,5}==0 ); % Non-targets
        names{counter} = sprintf('Noise+Speech_Set%d_Item%d',itemset,3);
        onsets{counter} = double(tmp{6}(conIdx{counter}));
        durations{counter} = zeros(length(onsets{counter}),1);    
    end
    
    counter = counter + 1;
    conIdx{counter} = find( ismember(tmp{:,21},'noiseBoth')==1 & tmp{:,5}==0 ); % Non-targets
    names{counter} = sprintf('Noise+Noise');
    onsets{counter} = double(tmp{6}(conIdx{counter}));
    durations{counter} = zeros(length(onsets{counter}),1);
    
    for itemset=1:8
        for i=1:4
            if i<=2
                counter = counter + 1;
                conIdx{counter} = [find( ismember(tmp{:,21},'clearSpeech')==1 & (tmp{:,18}==2 & tmp{:,17}==itemset & tmp{:,5}==0 & tmp{:,20}==1 & tmp{:,19}==i) );...
                                   find( ismember(tmp{:,21},'clearSpeech')==1 & (tmp{:,18}==2 & tmp{:,17}==itemset & tmp{:,5}==0 & tmp{:,20}==2 & tmp{:,19}==i+2) );...
                                   find( ismember(tmp{:,21},'noiseFinal')==1 & (tmp{:,18}==2 & tmp{:,17}==itemset & tmp{:,5}==0 & tmp{:,20}==1 & tmp{:,19}==i) ) ]; % Non-targets
                names{counter} = sprintf('StrongSyl1+Clear_Set%d_Item%d',itemset,i);
                onsets{counter} = double(tmp{6}(conIdx{counter}));
                durations{counter} = zeros(length(onsets{counter}),1);
            else
                counter = counter + 1;
                conIdx{counter} = [find( ismember(tmp{:,21},'clearSpeech')==1 & (tmp{:,18}==2 & tmp{:,17}==itemset & tmp{:,5}==0 & tmp{:,20}==1 & tmp{:,19}==i) );...
                                   find( ismember(tmp{:,21},'clearSpeech')==1 & (tmp{:,18}==2 & tmp{:,17}==itemset & tmp{:,5}==0 & tmp{:,20}==2 & tmp{:,19}==i-2) );...
                                   find( ismember(tmp{:,21},'noiseFinal')==1 & (tmp{:,18}==2 & tmp{:,17}==itemset & tmp{:,5}==0 & tmp{:,20}==1 & tmp{:,19}==i) ) ]; % Non-targets                
                names{counter} = sprintf('StrongSyl1+Clear_Set%d_Item%d',itemset,i);
                onsets{counter} = double(tmp{6}(conIdx{counter}));
                durations{counter} = zeros(length(onsets{counter}),1);
            end
        end
    end
    for itemset=1:8
        for i=1:4
            if i<=2            
                counter = counter + 1;
                conIdx{counter} = [find( ismember(tmp{:,21},'clearSpeech')==1 & (tmp{:,18}==1 & tmp{:,17}==itemset & tmp{:,5}==0 & tmp{:,20}==1 & tmp{:,19}==i) );...
                                   find( ismember(tmp{:,21},'clearSpeech')==1 & (tmp{:,18}==1 & tmp{:,17}==itemset & tmp{:,5}==0 & tmp{:,20}==2 & tmp{:,19}==i+2) );...
                                   find( ismember(tmp{:,21},'noiseFinal')==1 & (tmp{:,18}==1 & tmp{:,17}==itemset & tmp{:,5}==0 & tmp{:,20}==1 & tmp{:,19}==i) ) ]; % Non-targets                
                names{counter} = sprintf('WeakSyl1+Clear_Set%d_Item%d',itemset,i);
                onsets{counter} = double(tmp{6}(conIdx{counter}));
                durations{counter} = zeros(length(onsets{counter}),1);
            else              
                counter = counter + 1;
                conIdx{counter} = [find( ismember(tmp{:,21},'clearSpeech')==1 & (tmp{:,18}==1 & tmp{:,17}==itemset & tmp{:,5}==0 & tmp{:,20}==1 & tmp{:,19}==i) );...
                                   find( ismember(tmp{:,21},'clearSpeech')==1 & (tmp{:,18}==1 & tmp{:,17}==itemset & tmp{:,5}==0 & tmp{:,20}==2 & tmp{:,19}==i-2) );...
                                   find( ismember(tmp{:,21},'noiseFinal')==1 & (tmp{:,18}==1 & tmp{:,17}==itemset & tmp{:,5}==0 & tmp{:,20}==1 & tmp{:,19}==i) ) ]; % Non-targets
                names{counter} = sprintf('WeakSyl1+Clear_Set%d_Item%d',itemset,i);
                onsets{counter} = double(tmp{6}(conIdx{counter}));
                durations{counter} = zeros(length(onsets{counter}),1);
            end
        end
    end
    
    counter = counter + 1;
    conIdx{counter} = find( tmp{:,12}==1|tmp{:,13}==1 );
    names{counter} = 'HitsAndFalseAlarms';
    onsets{counter} = double(tmp{6}(conIdx{counter}));
    durations{counter} = zeros(length(onsets{counter}),1);
    
%     counter = counter + 1;
%     conIdx{counter} = find( tmp{:,14}==1&tmp{:,12}==0 );
%     names{counter} = 'Misses';
%     onsets{counter} = double(tmp{6}(conIdx{counter}));
%     durations{counter} = zeros(length(onsets{counter}),1);
%     if isempty(onsets{counter})
%         onsets{counter} = tmp{7}(end);
%         durations{counter} = 0;
%     end
    
    %     counter = counter + 1;
    %     conIdx{counter} = find( tmp{:,5}==1 );
    %     names{counter} = 'Pause';
    %     onsets{counter} = double(tmp{6}(conIdx{counter}))+double(tmp{9}(conIdx{counter}));
    %     durations{counter} = zeros(length(onsets{counter}),1);
    
    save([GLMAnalDir '/Regressors_Sess_' num2str(sess) '_' SubjID],'names', 'onsets', 'durations');
    save([GLMAnalDir '/DataFile_Session_' num2str(sess)],'DataFile');
    clear names onsets durations DataFile
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Motion regressors
    MotionRegFile=spm_select('FPList',MotionDataDir{sess},'^rp.*\.txt$');
    MotionReg=textread(MotionRegFile);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Physio Regressors
    if ~isempty(PhysioDataFile{sess});
        PhysioRegFile=PhysioDataFile{sess};
        load(PhysioRegFile);
        PhysioReg =R;
        clear R
    else
        PhysioReg =[];
    end
    
    % combine Motion and Physio
    try
        R =[MotionReg PhysioReg];
    catch % if missing physio reg for part of session
        PhysioReg = [PhysioReg; zeros(size(MotionReg,1)-size(PhysioReg,1),size(PhysioReg,2))];
        R =[MotionReg PhysioReg];
    end
    save([GLMAnalDir '/MotionPhysio_Sess_' num2str(sess) '_' SubjID],'R');
    save([GLMAnalDir '/MotionRegFile_Session_' num2str(sess)],'MotionRegFile');
    if ~isempty(PhysioReg)
        save([GLMAnalDir '/PhysioRegFile_Session_' num2str(sess)],'PhysioRegFile');
    end
    
    clear R PhysioReg MotionReg
end