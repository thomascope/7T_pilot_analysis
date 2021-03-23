function [] = CreateRegressors(GLMAnalDir,BehavDataDir,PhysioDataFile,MotionDataDir,SubjID)

NumSess = length(BehavDataDir);

for sess=1:NumSess
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create Stimulus onset regressors
    
    DataFile = dir(sprintf('%s%s_run%d_*.txt',BehavDataDir{sess},SubjID,sess));
    DataFile = fullfile(BehavDataDir{sess},DataFile(1).name);
    fid = fopen(DataFile);
    tmp = textscan(fid,'%f%f%s%f%f%f%f%f%f%f%f%f%f%f%f%s%f%f%f%f%s%s%f','Delimiter','\t','HeaderLines',1);
    
    %
    fid2 = fopen('items_selected_full_postHoc.txt');
    tmp2 = textscan(fid2,'%s%s%s%s%f%f%f%f%f%f%f','Delimiter','\t','HeaderLines',0);
    probM = tmp2{7};
    probM(find(isnan(probM))) = [];
    probM = tiedrank(probM);
    probM = probM-mean(probM);
    probItemsM = strrep(tmp2{3},'-','');
    probItemsM(find(cellfun(@isempty,probItemsM))) = [];
    clear probMM
    for i=1:4:length(probM)
        probMM(i:i+3,1) = [probM(i+2:i+3); probM(i:i+1)];
    end
    probItemsMM = strrep(tmp2{4},'-','');
    probItemsMM(find(cellfun(@isempty,probItemsMM))) = [];
    
    probSyl1M = tmp2{11};
    probSyl1M(find(isnan(probSyl1M))) = [];
    probSyl1M = tiedrank(probSyl1M);
    probSyl1M = probSyl1M-mean(probSyl1M);
    clear probSyl1MM
    for i=1:4:length(probSyl1M)
        probSyl1MM(i:i+3,1) = [probSyl1M(i+2:i+3); probSyl1M(i:i+1)];
    end
    fclose(fid2);
    
    % Make event times relative to time of first scan
    t0 = tmp{7}(1);
    tmp{6} = tmp{6}-t0; % Sound onset
    tmp{7} = tmp{7}-t0; % Volume onset
    tmp{10} = tmp{10}-t0; % Button press
    
    counter = 0;
    
    counter = counter + 1;
    conIdx{counter} = find( ismember(tmp{:,21},'clearSpeech')==1 & tmp{:,20}==1 & tmp{:,5}==0 ); % Non-targets
    names{counter} = 'M';
    onsets{counter} = double(tmp{6}(conIdx{counter}));
    durations{counter} = zeros(length(onsets{counter}),1);
    orth{counter} = 0;
    pmod(counter).name{1} = 'M_Prob';
    pmod(counter).poly{1} = 1;
    for i=1:length(conIdx{counter})
        pmod(counter).param{1}(i) = probM(find(ismember(probItemsM,tmp{16}(conIdx{counter}(i)))));
    end
    pmod(counter).name{2} = 'M_ProbSyl1';
    pmod(counter).poly{2} = 1;
    for i=1:length(conIdx{counter})
        pmod(counter).param{2}(i) = probSyl1M(find(ismember(probItemsM,tmp{16}(conIdx{counter}(i)))));
    end

    counter = counter + 1;
    conIdx{counter} = find( ismember(tmp{:,21},'clearSpeech')==1 & tmp{:,20}==2 & tmp{:,5}==0 ); % Non-targets
    names{counter} = 'MM';
    onsets{counter} = double(tmp{6}(conIdx{counter}));
    durations{counter} = zeros(length(onsets{counter}),1);
    orth{counter} = 0;
    pmod(counter).name{1} = 'MM_Prob';
    pmod(counter).poly{1} = 1;
    for i=1:length(conIdx{counter})
        pmod(counter).param{1}(i) = probMM(find(ismember(probItemsMM,tmp{16}(conIdx{counter}(i)))));
    end
    pmod(counter).name{2} = 'MM_ProbSyl1';
    pmod(counter).poly{2} = 1;
    for i=1:length(conIdx{counter})
        pmod(counter).param{2}(i) = probSyl1MM(find(ismember(probItemsMM,tmp{16}(conIdx{counter}(i)))));
    end
 
    counter = counter + 1;
    conIdx{counter} = find( ismember(tmp{:,21},'noiseFinal')==1 & tmp{:,20}==1 & tmp{:,5}==0 ); % Non-targets
    names{counter} = 'Clear+Noise';
    onsets{counter} = double(tmp{6}(conIdx{counter}));
    durations{counter} = zeros(length(onsets{counter}),1);
    orth{counter} = 0;
    pmod(counter).name{1} = 'Clear+Noise_Prob';
    pmod(counter).poly{1} = 1;
    for i=1:length(conIdx{counter})
        pmod(counter).param{1}(i) = probM(find(ismember(probItemsM,tmp{16}(conIdx{counter}(i)))));
    end
    pmod(counter).name{2} = 'Clear+Noise_ProbSyl1';
    pmod(counter).poly{2} = 1;
    for i=1:length(conIdx{counter})
        pmod(counter).param{2}(i) = probSyl1M(find(ismember(probItemsM,tmp{16}(conIdx{counter}(i)))));
    end
    
    counter = counter + 1;
    conIdx{counter} = find( ismember(tmp{:,21},'noiseInitial')==1 & tmp{:,20}==1 & tmp{:,5}==0 ); % Non-targets
    names{counter} = 'Noise+Clear';
    onsets{counter} = double(tmp{6}(conIdx{counter}));
    durations{counter} = zeros(length(onsets{counter}),1);

    counter = counter + 1;
    conIdx{counter} = find( ismember(tmp{:,21},'noiseBoth')==1 & tmp{:,5}==0 ); % Non-targets
    names{counter} = sprintf('Noise+Noise');
    onsets{counter} = double(tmp{6}(conIdx{counter}));
    durations{counter} = zeros(length(onsets{counter}),1);
    
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
    
    save([GLMAnalDir '/Regressors_Sess_' num2str(sess) '_' SubjID],'names', 'onsets', 'durations','pmod','orth');
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