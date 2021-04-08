function am_test_mvpa_search(search_this_cond)
%{
in this tarball you will find the 2 sessions from the same subject (i.e. 101295 is session 1 and
101363 is session 2).

The two sessions are aligned, and the structural of the 1st session is also
aligned within a subject.

The 18 first images (after applying reorderBrainBetas are
the subject's own images from his photo-album, the following 18 images are
the other paired subject's own images (we will not analyse that subject's fMRI data today
and then the last 36 images are the general images that every subject from every
pair saw.

the similarity judgments RDM are contained in the
pair1_subj1_extra.mat in the structure behav_and_icons.RDM

Here are some more notes on the stimulus indices:
 
OwnBodyParts = 1:3; OwnFaces = 4:8; OwnPet = 9; OwnPlaces = 10:15; OwnObjects=16:18;
OtherBodyParts = 19:21; OtherFaces = 22:26; OtherPet = 27; OtherPlaces = 28:33; OtherObjects=34:36;
GeneralBodyParts = 37:44; GeneralFaces = 45:52; GeneralPets = [53 54]; GeneralPlaces=55:66; GeneralObjects=67:72;

note that own refer to the subjects' own images, and other to the paired
subject's own images, and general to objects that were seen by all subjects.

bodies = [OwnBodyParts OtherBodyParts GeneralBodyParts];
nbodies = length(bodies);
faces = [OwnFaces OwnPet OtherFaces OtherPet GeneralFaces GeneralPets];
nfaces = length(faces);
places = [OwnPlaces OtherPlaces GeneralPlaces];
nplaces = length(places);
objects = [OwnObjects OtherObjects GeneralObjects];
nobjects = length(objects);

vowel1 = [faces bodies];
vowel2 = [places objects];


% ALREADY DONE
% 1 - preprocess and split the data
% 2 - estimate single-subject activity patterns

% LINEAR SVM
% 2 - review activity-pattern estimation
        inspect the design matrices for the runs in both sessions
        load the single-subject activity patterns (t-patterns)
% 3 - select voxels
        constrain the t-patterns to the hIT masks (actually called  VOT)
% 4 - train the classifier
%       train a linear SVM classifier to distinguish animate from inanimate
%       objects
% 5 - test the classifier
% 6 - statistical inference
%       use a condition-label permutation test to determine whether the
%       classifier performs above chance in this particular subject

% RSA
% 7 - generate an RDM for each session
% 8 - compare RDMs across sessions
% 9 - relate the session-averaged hIT RDM to the similarity judgments RDM

%}

%clear all
[workingdir,file,ext] = fileparts(which('am_test_mvpa.m'));

analyse.lSVM=0;
analyse.RSA=0;
analyse.searchlight=1;

addpath(genpath('/imaging/mlr/users/tc02/toolboxes'));

%% control parameters

subj_initials = 'am';
testing_date = '20170711';
nSessions    = 1;
nRuns        = 4;   % there was a total of 6 runs per subject and session
nWords  = 9;  % there were 9 different vowels presented in the experiment
nTrialtypes = 8; % there were 8 conditions - Match low, Match high, Mismatch low, Mismatch high, Neutral low, Neutral high, Writtenonly, Response
nVolumes     = 210; % there were at least 200 EPI volumes acquired per run
tr           = 2.5;   % the EPI volumes were acquired every 3 seconds
stimdur      = 0.6;   % the audio stimuli lasted approx 0.5 seconds
stimdelay    = 2.20; % the audio stimulus was presented 2200ms after the scanner pulse
% XXX STIMDELAY NOT YET IMPLEMENTED. INSTEAD MODEL WRITTEN AND SPOKEN
% TOGETHER AS A 3 SECOND STIMULUS PRESENTED AT SCANNER PULSE
stimdur      = 2.5;
condition_order = {'Match low','Match high','Mismatch low','Mismatch high','Neutral low','Neutral high','Writtenonly','Response'};

% % load the subject's image ordering (the same ordering is valid for both sessions)

% pre allocate cell arrays
tempOnsets   = cell(1,nRuns);
design       = cell(1,nRuns);
tpatterns    = cell(1,nSessions);

% the idendity of the subject will be subject initials
thisSubject = subj_initials;

% % load hIT mask
% mask = logical(spm_read_vols(spm_vol(fullfile(workingdir,thisSubject,'masks','VOT','bilateral.VOT.nii'))));

structural_fname = 'mwc1DATA_0011.nii';
%wholespace_mask = 'wholespace_mask.nii';
wholespace_mask = 'Wholespace_test.nii';
%%% NEED A BETTER MASK
%mask = logical(spm_read_vols(spm_vol('Univariate_Mask.nii')));
mask = logical(spm_read_vols(spm_vol('Univariate_cluster.nii')));

% % vectorize the mask
vmask = squish(mask,3);


% loop over the trial types


% prepare design matrices
tempOnsets = cell(1,nTrialtypes*nWords+2);
tempDesign = cell(1,nRuns);

for runI=1:nRuns
    
    this_file = dir(['*' subj_initials '_Run_' num2str(runI) '_' testing_date '.mat']);
    fileName = this_file.name;
    %read the scanning log file
 
    [startpulses,stimType,stimNumber,stimName] = extract_pulsenumbers_from_second_fullparadigm(fileName,runI);
    
    starttime = ((startpulses-1) * tr) + stimdelay; % Remember first pulse occurs at time zero
    
    if size(stimNumber,2) ~= size(starttime,2)
        warning('Something has gone wrong in recording the stimulus types. Padding with NaNs but you MUST check your data')
        stimNumber = [stimNumber NaN(1,(size(starttime,2)-size(stimNumber,2)))];
    end
    
    for condJ = 1:nTrialtypes
        for condI=1:nWords
            tempOnsets{condI+((condJ-1)*nWords)}=starttime(stimNumber==condI & stimType==condJ);
        end
    end
    
    tempOnsets{nTrialtypes*nWords+1}=starttime(stimNumber==0 & stimType==7);
    tempOnsets{nTrialtypes*nWords+2}=starttime(stimNumber==0 & stimType==8);
    
    tempDesign{runI} = tempOnsets;
end

%Create a design matrix for each run 
design = cell(1,nRuns);
for runI=1:nRuns
    design{runI} = zeros(nVolumes,length(tempOnsets));
    for condI=1:length(tempOnsets)
        design{runI}(round((tempDesign{runI}{condI}-stimdelay)/tr)+1,condI) = 1;
    end
    design{runI} = design{runI}(1:nVolumes,:);
end

% plot one run's design matrix
figure;
set(gcf,'Position',[100 100 400 800],'Color','w')
imagesc(design{1});colormap(gray)
ylabel('\bf{number of volumes}')
xlabel('\bf{number of conditions}')
title(sprintf('design matrix for run 1'),'FontSize',14,'Fontweight','bold')

if ~(exist(['tpatterns_masked_' subj_initials '.mat'],'file'))
    %Create whole brain masked input files
    
    if isempty(dir(['mPNFA*' subj_initials '_Run_1_' testing_date '.nii']));
        input_files = dir(['PNFA*' subj_initials '_Run_*_' testing_date '.nii']);
        cbupool(nRuns)
        parfor i = 1:nRuns
            spm_mask(structural_fname, input_files(i).name, 0.5)
        end
        matlabpool close
    end
    
    % Split the data by run
    tpatterns = cell(1,2);
    wholebrain_tpatterns = cell(1,2);
    for session = 1:2
        data = cell(1,session:2:nRuns);
        % load the nifti data;
        for runI=1:(nRuns/2) %XXX Relies on an even number of runs at present
            fprintf('***\t importing EPI time-series for run %d \t***\n',session+(2*(runI-1)))
            this_file = dir(['mP*' subj_initials '_Run_' num2str(session+(2*(runI-1))) '_' testing_date '.nii']);
            fileName = this_file.name;
            thisRun = fullfile(workingdir,fileName);
            data{runI} = single(spm_read_vols(spm_vol(thisRun)));
            data{runI} = data{runI}(:,:,:,1:nVolumes); %Trunkate data at nVolumes to prevent errors if this doens't match between runs.
        end
        
        % run GLMdenoise
        %     results = GLMdenoisedata(design,data,stimdur,tr, ...
        %         'optimize',[],struct('numboots',100,'numpcstotry',20,'wantparametric',1), ...
        %         []);
        % XXX Suboptimal - bootstrapping turned off to prevent out-of-memory crash
        clear results %Necessary for memory management
        results = GLMdenoisedata(design(session:2:nRuns),data,stimdur,tr, ...
            'optimize',[],struct('numboots',50,'numpcstotry',8,'wantparametric',1), ...
            ['testfigures_session' num2str(session)]);
        
        % limit the betas to the valid conditions
        modelmd = results.modelmd{2}(:,:,:,1:(nWords*nTrialtypes));
        % limit the standard errors to the valid conditions
        modelse = results.modelse{2}(:,:,:,1:(nWords*nTrialtypes));
        % get the pooled error
        poolse  = sqrt(mean(modelse.^2,4));
        % normalise the betas by the pooled error to get t-patterns
        modelmd = bsxfun(@rdivide,modelmd,poolse);
        
        % show the unmasked patterns
        figure(1);
        set(gcf,'Position',[100 100 600 800],'Color','w')
        subplot(3,1,1)
        imagesc(makeimagestack(mean(modelmd,4),-1))
        title('\bf{unmasked mean t-pattern}')
        subplot(3,1,2)
        imagesc(makeimagestack(mask,-1))
        title('\bf{hIT mask}')
        subplot(3,1,3)
        imagesc(makeimagestack(mean(modelmd,4).*mask,-1)) %XXX Commented out at present as mask larger than image
        title('\bf{mean t-pattern masked}')
        
        vmodelmd     = squish(modelmd,3);
        
        
        % masked t-pattenrs will be nvoxels x conditions
        tpatterns{session}    = vmodelmd(vmask,:);
        % With no mask for later searchlight
        wholebrain_tpatterns{session} = vmodelmd;
    end
end

% save the data if not existing already otherwise load it
if ~(exist(['tpatterns_masked_' subj_initials '.mat'],'file'))
    % save the t-patterns
    save(['tpatterns_masked_' subj_initials '.mat'], 'tpatterns','wholebrain_tpatterns')
else
    % if the t-patterns already exist --> load it
    %load('tpatterns.mat')
    load(['tpatterns_masked_' subj_initials '.mat'])
end


%% vowel classification using linear SVM

if analyse.lSVM
    
    
    for condJ = 1:nTrialtypes
      % define vowel vectors  
%   
%     vowel1=[1:3]+(9*(condJ-1));
%     vowel2=[4:6]+(9*(condJ-1));
%     vowel3=[7:9]+(9*(condJ-1));
    
    vowel1=[1:3];
    vowel2=[4:6];
    vowel3=[7:9];
    
    % control variables
    %libSVMsettings='-s 1 -t 0'; % nu-SVM, linear
    libSVMsettings='-s 1 -t 2'; % nu-SVM, radial basis
    nRandomisations=100; %XXX CHANGE TO 1000 WHEN CODE IS FINISHED
        
    % linear SVM
    cvFolds=[1 2; 2 1]; % columns = folds, row 1 = session used for training, row 2 = session used for testing
    for foldI=1:size(cvFolds,2)
        % define training and test data sets
        tpatternsTrain=tpatterns{cvFolds(1,foldI)}; tpatternsTrain=double(tpatternsTrain');
        tpatternsTest=tpatterns{cvFolds(2,foldI)}; tpatternsTest=double(tpatternsTest');
        
        tpatternsTrain = tpatternsTrain(1+((condJ-1)*9):9+((condJ-1)*9),:);
        tpatternsTest = tpatternsTest(1+((condJ-1)*9):9+((condJ-1)*9),:);
        % define class lables
        labels=zeros((size(vowel1,2)+size(vowel2,2)+size(vowel3,2)),1);
        labels(vowel1)=1;
        labels(vowel2)=2; 
        labels(vowel3)=3; 
        % train and test the classifier
        model=svmtrain(labels,tpatternsTrain,libSVMsettings);
        %model=svmtrain(labels,tpatternsTrain);
        [labelsPredicted,accuracy,decVals]=svmpredict(labels,tpatternsTest,model);
        accuracy_fold(foldI)=accuracy(1);
        
        % create null distribution for statistical inference
        for randomisationI=1:nRandomisations
            % randomise labels (for training)
            labelsRand=labels(randperm(length(labels)));
            % train and test the classifier using the randomised training labels
            modelRand=svmtrain(labelsRand,tpatternsTrain,libSVMsettings);
            [labelsPredictedRand,accuracyRand,decValsRand]=svmpredict(labels,tpatternsTest,modelRand);
            accuracy_randomisation_fold(randomisationI,foldI)=accuracyRand(1);
        end % randomisationI
    end % foldI
    
    % statistical inference
    accuracy=mean(accuracy_fold);
    accuracyH0=mean(accuracy_randomisation_fold,2);
    p=1-relRankIn_includeValue_lowerBound(accuracyH0,accuracy);
    
    % visualise results
    figure(100+condJ); clf;
    % plot null distribution
    hist(accuracyH0); hold on;
    % plot accuracy (mean across folds) found in the data
    xlim([5 95]); xlims=xlim;
    plot(accuracy,0,'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',8);
    ylims=ylim;
    text(accuracy,0.04*ylims(2),'\bfdata','Color','r');
    % plot statistical result
    text(0.85*xlims(2),0.9*ylims(2),['p = ',sprintf('%1.4f',p)]);
    % label axes
    xlabel('accuracy');
    ylabel('frequency');
    title({'\fontsize{11}null distribution of classification accuracy',['\fontsize{8}',num2str(nRandomisations),' stimulus-label randomisations']})
    end
end


%% RSA

if analyse.RSA
    for condJ = 1:nTrialtypes
    %Define matrix based on shared features
    judgmentRDM.RDM = zeros(9,9);
    judgmentRDM.RDM(1:10:end) = 1;
    judgmentRDM.RDM(2:30:end) = 2/3;
    judgmentRDM.RDM(3:30:end) = 1/3;
    judgmentRDM.RDM(10:30:end) = 2/3;
    judgmentRDM.RDM(12:30:end) = 2/3;
    judgmentRDM.RDM(20:30:end) = 2/3;
    judgmentRDM.RDM(19:30:end) = 1/3;
    judgmentRDM.RDM = 1-judgmentRDM.RDM;
    judgmentRDM.name = 'feature_vowel';
    
    % now we make an RDM per session
    RDMs=struct();
    empty_conditions = {};
    for sessionI=1:2
        empty_conditions{sessionI} = logical(sum(abs(tpatterns{sessionI}(:,1+((condJ-1)*9):9+((condJ-1)*9))))==0);
    end
    if ~isempty(empty_conditions)
        all_empty_conditions = logical(empty_conditions{1}+empty_conditions{2});
        display(['Warning, condition(s) ' num2str(find(all_empty_conditions)) ' are empty in at least one run. They will be ignored'])
        %pause
    end
    
    for sessionI=1:2
        % the correlation distance patterns are computed using the pdist function
        thesepatterns = tpatterns{sessionI}(:,1+((condJ-1)*9):9+((condJ-1)*9))';
        thesepatterns = thesepatterns(logical(1-all_empty_conditions),:);
        RDMs(1,sessionI).RDM   = squareform(pdist(thesepatterns,'correlation'));
        RDMs(1,sessionI).name  = sprintf('Univariate Mask | session %d | condition %s excluding %s',sessionI,condition_order{condJ},num2str(find(all_empty_conditions)));
        RDMs(1,sessionI).color = [];
    end
    
    % show the 2 session RDMs
    figI=condJ;
    figure(figI);clf;set(gcf,'Position',[100 100 800 800],'Color','w')
    showRDMs(RDMs,figI);
    
    % ---------------------------------------------------------------------
    % compare RDMs across sessions
    r12=corr([vectorizeRDMs(RDMs(1).RDM)]',[vectorizeRDMs(RDMs(2).RDM)]');
    % ---------------------------------------------------------------------
    
    avgRDM = averageRDMs_subjectSession(RDMs,'subject');
    avgRDM.name=sprintf('RDM across sessions | condition %s',condition_order{condJ});
    
    figure(figI);set(gcf,'Position',[100 100 800 800],'Color','w')
    showRDMs(avgRDM,figI);
    
    
    % define the labels and indices for familiar and unfamiliar images
    reductionLabels = {'vowel1','vowel2','vowel3'};
    reductionvectors = logical(1-all_empty_conditions);
    nobjects = sum(logical(1-all_empty_conditions));
    
    %reduction = reductionvectors;

    Vowel1 = 1:3; Vowel2 = 4:6; Vowel3 = 7:9; 
    
    nCols=3;
    cmap=RDMcolormap;
    colors=cmap([1 111 222],:);
    options.categoryColors=zeros(9,3);
    options.categoryColors(Vowel1,:)=repmat(colors(1,:),length(Vowel1),1);
    options.categoryColors(Vowel2,:)=repmat(colors(2,:),length(Vowel2),1);
    options.categoryColors(Vowel3,:)=repmat(colors(3,:),length(Vowel3),1);
    options.categoryColors(all_empty_conditions,:) = [];
    
    options.spheres=2;
    options.cols=options.categoryColors;
    options.replicability=0;
    options.view=1;
    
    %D=avgRDM.RDM(reduction,reduction);
    D=avgRDM.RDM;
    %[pats_mds_2D,stress,disparities]=extractMDS(D,2,options);
    [pats_mds_2D,stress,disparities]=mdscale(D,2);
    
    % draw the mds
    nImages=size(pats_mds_2D,1);
    
    % compute image size
    imageAreaProportion=.5;
    boundingBoxArea=max(prod(range(pats_mds_2D)),max(range(pats_mds_2D))^2/10);
    totalImageArea=boundingBoxArea*imageAreaProportion;
    imageWidth=sqrt(totalImageArea/nImages);
    
    % smooth alpha channel
    transparentCol=[128 128 128 2];
    hsize=5*transparentCol(4);
    sigma=1*transparentCol(4);
    kernel=fspecial('gaussian', hsize, sigma);
    
    markerSize=85;
    figure(figI+10);clf;
    set(gcf,'Position',[ 100 100 800 800],'Color',[1 1 1],'Renderer','OpenGL','BackingStore','on'); % much better
    axes('Position',[0.05 0.2 0.9 0.75])
 
    
    hold on
%     for imageI=1:nImages
%         
%         plot(pats_mds_2D(imageI,1),pats_mds_2D(imageI,2),'o','MarkerFaceColor',[128 128 128]./255,'MarkerEdgeColor',[128 128 128]./255,'MarkerSize',markerSize+20);
%         
%         
%         %[xs,ys,rgb3]=size(imageStruct(imageI).image);
%         
%         if reductionI==1
%             transparent=imageIcons(imageI).image(:,:,1)==transparentCol(1) & imageIcons(imageI).image(:,:,2)==transparentCol(2) & imageIcons(imageI).image(:,:,3)==transparentCol(3);
%         else
%             transparent=imageIcons(imageI+36).image(:,:,1)==transparentCol(1) & imageIcons(imageI+36).image(:,:,2)==transparentCol(2) & imageIcons(imageI+36).image(:,:,3)==transparentCol(3);
%         end
%         if numel(transparentCol)==4
%             % smooth alpha channel
%             opacity=imfilter(double(1-transparent),kernel);
%         else
%             opacity=~transparent;
%         end
%         if reductionI==1
%             if imageI<=18
%                 plot(pats_mds_2D(imageI,1),pats_mds_2D(imageI,2),...
%                     'o','MarkerFaceColor',[128 128 128]./255,'MarkerEdgeColor',[128 128 128]./255,'MarkerSize',markerSize+20);
%             end
%         end
%         
%         %     plot(pats_mds_2D(imageI,1),pats_mds_2D(imageI,2),...
%         %         'o','MarkerFaceColor',options.categoryColors(imageI,:),'MarkerEdgeColor',options.categoryColors(imageI,:),'MarkerSize',markerSize);
%         %
%         %imagesc(npats_mds_2D(imageI,1),npats_mds_2D(imageI,2),imageIcons(imageI+36).image);
%         if reductionI==1
%             image('CData',imageIcons(imageI).image,'XData',[pats_mds_2D(imageI,1)-imageWidth/2, pats_mds_2D(imageI,1)+imageWidth/2],'YData',[pats_mds_2D(imageI,2)+imageWidth/2, pats_mds_2D(imageI,2)-imageWidth/2],'AlphaData',opacity);
%         else
%             image('CData',imageIcons(imageI+36).image,'XData',[pats_mds_2D(imageI,1)-imageWidth/2, pats_mds_2D(imageI,1)+imageWidth/2],'YData',[pats_mds_2D(imageI,2)+imageWidth/2, pats_mds_2D(imageI,2)-imageWidth/2],'AlphaData',opacity);
%         end
%     end
%     axis tight equal off
%     annotation('textbox',[0 .90 1 0.1],'EdgeColor','none','String','MDS plot for unfamiliar images',...
%         'HorizontalAlignment','center','FontWeight','bold','FontSize',18);
%     
%     axes('Position',[0.7 0 0.25 0.25])
%     hold on;
    % add a micro mds plot
    markerSize=18;
    for imageI=1:nImages
        
        plot(pats_mds_2D(imageI,1),pats_mds_2D(imageI,2),...
            'o','MarkerFaceColor',options.categoryColors(imageI,:),'MarkerEdgeColor',options.categoryColors(imageI,:),'MarkerSize',markerSize);
    end
    axis tight equal off
    
    % relate activation pattern and candidate models
    userOptions.candRDMdifferencesTest='conditionRFXbootstrap';
    userOptions.nBootstrap=100; % XXX CHange to 10000 when code finalised
    judgments=cell(1,2);
    judgments{1}=judgmentRDM;
    judgments{1}.RDM = judgments{1}.RDM(reductionvectors,reductionvectors);
    judgments{2}.name = 'vowels_only';
    judgments{2}.RDM = ones(9,9);
    judgments{2}.RDM(1:10:end) = 1;
    judgments{2}.RDM(2:30:end) = 0.5;
    judgments{2}.RDM(3:30:end) = 0.5;
    judgments{2}.RDM(10:30:end) = 0.5;
    judgments{2}.RDM(12:30:end) = 0.5;
    judgments{2}.RDM(20:30:end) = 0.5;
    judgments{2}.RDM(19:30:end) = 0.5;
    judgments{2}.RDM = judgments{2}.RDM(reductionvectors,reductionvectors);
    
    
    userOptions.figureIndex = [20+condJ, 30+condJ];
    stats_p_r=compareRefRDM2candRDMs(avgRDM, judgments, userOptions);
    
    % ---------------------------------------------------------------------
    % add a second model
    
    
    
    
    % ---------------------------------------------------------------------
    end
end

if analyse.searchlight
    
    if ~(exist(['tpatterns_masked_' subj_initials '.mat'],'file'))
        input('Must run the SVM first, to create the whole brain t-patterns for searchlight analysis')
    else
        % if the t-patterns already exist --> load it
        %load('tpatterns.mat')
        load(['tpatterns_masked_' subj_initials '.mat'])
    end
    
    brainmask = logical(spm_read_vols(spm_vol(wholespace_mask)));
    binaryMasks_nS.subj1.mask = brainmask;
    
    fullBrainVols.subj1 = zeros(size(wholebrain_tpatterns{1}));
    fullBrainVols.subj1 = repmat(fullBrainVols.subj1,[1,1,size(wholebrain_tpatterns,2)]);
    for session = 1:size(wholebrain_tpatterns,2)
        fullBrainVols.subj1(:,:,session) = wholebrain_tpatterns{session};
    end
    clear wholebrain_tpatterns %For memory

    addpath('/imaging/mlr/users/tc02/7T_full_paradigm_pilot_second_am/Preprocessed_Images')
    for condJ = search_this_cond:search_this_cond
        %Define matrix based on shared features
        theseBrainVols.subj1 = fullBrainVols.subj1(:,1+((condJ-1)*9):9+((condJ-1)*9),:);
        
        judgmentRDM.RDM = zeros(9,9);
        judgmentRDM.RDM(1:10:end) = 1;
        judgmentRDM.RDM(2:30:end) = 2/3;
        judgmentRDM.RDM(3:30:end) = 1/3;
        judgmentRDM.RDM(10:30:end) = 2/3;
        judgmentRDM.RDM(12:30:end) = 2/3;
        judgmentRDM.RDM(20:30:end) = 2/3;
        judgmentRDM.RDM(19:30:end) = 1/3;
        judgmentRDM.RDM = 1-judgmentRDM.RDM;
        judgmentRDM.name = 'feature_vowel';
        
        % now we make an RDM per session
        RDMs=struct();
        empty_conditions = {};
        for sessionI=1:2
            empty_conditions{sessionI} = logical(sum(abs(tpatterns{sessionI}(:,1+((condJ-1)*9):9+((condJ-1)*9))))==0);
        end
        if ~isempty(empty_conditions)
            all_empty_conditions = logical(empty_conditions{1}+empty_conditions{2});
            display(['Warning, condition(s) ' num2str(find(all_empty_conditions)) ' are empty in at least one run. They will be ignored'])
            %pause
        end
        reductionvectors = logical(1-all_empty_conditions);
        nobjects = sum(logical(1-all_empty_conditions));
        
        % relate activation pattern and candidate models
        clear judgments
        judgments(1)=judgmentRDM;
        judgments(1).RDM = judgments(1).RDM(reductionvectors,reductionvectors);
        judgments(2).name = 'vowels_only';
        judgments(2).RDM = ones(9,9);
        judgments(2).RDM(1:10:end) = 1;
        judgments(2).RDM(2:30:end) = 0.5;
        judgments(2).RDM(3:30:end) = 0.5;
        judgments(2).RDM(10:30:end) = 0.5;
        judgments(2).RDM(12:30:end) = 0.5;
        judgments(2).RDM(20:30:end) = 0.5;
        judgments(2).RDM(19:30:end) = 0.5;
        judgments(2).RDM = judgments(2).RDM(reductionvectors,reductionvectors);
        
        theseBrainVols.subj1 = theseBrainVols.subj1(:,reductionvectors,:);
        
        userOptions.figureIndex = [20+condJ, 30+condJ];
        
        betaCorrespondence = [];
        betaCorrespondence.identifier = '0001';
        betaCorrespondence = repmat(betaCorrespondence,[size(fullBrainVols.subj1,3),nobjects]);
        
        userOptions.betaPath = '/imaging/mlr/users/tc02/7T_full_paradigm_pilot_second_am/Preprocessed_Images/am_univariate_gmmasked/beta_[[betaIdentifier]].nii';
        %userOptions.subjectNames = {subj_initials};
        userOptions.analysisName = ['Condition_' num2str(condJ)];
        userOptions.rootPath = '/imaging/mlr/users/tc02/7T_full_paradigm_pilot_second_am/Preprocessed_Images/Searchlight';
        if ~exist(userOptions.rootPath,'dir')
            mkdir(userOptions.rootPath)
        end
        userOptions.voxelSize = [1.4 1.4 1.4];
        userOptions.searchlightRadius = 7;
        
        %[rMaps_sS, maskedSmoothedRMaps_sS, searchlightRDMs, rMaps_nS, nMaps_nS] = fMRISearchlight(theseBrainVols,binaryMasks_nS,judgments,betaCorrespondence,userOptions);
        
        fMRISearchlight(theseBrainVols,binaryMasks_nS,judgments,betaCorrespondence,userOptions);
    end
    

end
