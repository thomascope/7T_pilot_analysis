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
    clear X
    fid2 = fopen('items_selected_full_postHoc.txt');
    tmp2 = textscan(fid2,'%s%s%s%s%f%f%f%f%f%f%f','Delimiter','\t','HeaderLines',0);
    itemsM = [strrep(tmp2{3},'-','');];
    itemsM(find(cellfun(@isempty,itemsM))) = [];
    itemsMM = [strrep(tmp2{4},'-','');];
    itemsMM(find(cellfun(@isempty,itemsMM))) = [];
    items = [itemsM(1:64); itemsMM(1:64)];
    fclose(fid2);
    load('/group/language/data/ediz.sohoglu/projects/fMRI_2017/scripts_fMRI/PredictorsFromSimSegVersion/predictor_sim_univariate_peAbs_syl1.mat','predictorForfMRI');
    predictorForfMRI(predictorForfMRI==0) = 1e-3;
    predictorForfMRI = log(predictorForfMRI);
    predictorForfMRI = zscore(predictorForfMRI);
    X(:,1) = predictorForfMRI;
    load('/group/language/data/ediz.sohoglu/projects/fMRI_2017/scripts_fMRI/PredictorsFromSimSegVersion/predictor_sim_univariate_peAbs_syl2.mat','predictorForfMRI');
    predictorForfMRI(predictorForfMRI==0) = 1e-3;
    predictorForfMRI = log(predictorForfMRI);
    predictorForfMRI = zscore(predictorForfMRI);
    X(:,2) = predictorForfMRI;
    
    % Make event times relative to time of first scan
    t0 = tmp{7}(1);
    tmp{6} = tmp{6}-t0; % Sound onset
    tmp{7} = tmp{7}-t0; % Volume onset
    tmp{10} = tmp{10}-t0; % Button press
    tmp{16}; % Items presented
    
    counter = 0;
    
    counter = counter + 1;
    conIdx{counter} = find( ismember(tmp{:,21},'clearSpeech')==1 & tmp{:,5}==0 ); % Non-targets
    names{counter} = 'Clear+Clear';
    onsets{counter} = double(tmp{6}(conIdx{counter}));
    durations{counter} = zeros(length(onsets{counter}),1);
    orth{counter} = 0;
    pmod(counter).name{1} = 'Syl1';
    pmod(counter).poly{1} = 1;
    for i=1:length(conIdx{counter})
        pmod(counter).param{1}(i) = X(find(ismember(items,tmp{16}(conIdx{counter}(i)))),1);
    end
    pmod(counter).name{2} = 'Syl2';
    pmod(counter).poly{2} = 1;
    for i=1:length(conIdx{counter})
        pmod(counter).param{2}(i) = X(find(ismember(items,tmp{16}(conIdx{counter}(i)))),2);
    end
 
    counter = counter + 1;
    conIdx{counter} = find( ismember(tmp{:,21},'noiseFinal')==1 & tmp{:,20}==1 & tmp{:,5}==0 ); % Non-targets
    names{counter} = 'Clear+Noise';
    onsets{counter} = double(tmp{6}(conIdx{counter}));
    durations{counter} = zeros(length(onsets{counter}),1);
    
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