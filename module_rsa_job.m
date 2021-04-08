function [avgRDM, stats_p_r] = module_rsa_job(tpattern_numbers,mask_path,data_path,condition_number,condition_name)

condJ = condition_number+100; %To avoid figure numbers over-writing open SPM windows

addpath(genpath('/imaging/mlr/users/tc02/toolboxes')); %Where is the RSA toolbox?


%Define matrices based on shared features

%Squares based on vowels
judgmentRDM.RDM = zeros(16,16);
judgmentRDM.RDM(1:17:end) = 1;
judgmentRDM.RDM(2:68:end) = 1/3;
judgmentRDM.RDM(3:68:end) = 1/3;
judgmentRDM.RDM(4:68:end) = 1/3;
judgmentRDM.RDM(17:68:end) = 1/3;
judgmentRDM.RDM(19:68:end) = 1/3;
judgmentRDM.RDM(20:68:end) = 1/3;
judgmentRDM.RDM(33:68:end) = 1/3;
judgmentRDM.RDM(34:68:end) = 1/3;
judgmentRDM.RDM(36:68:end) = 1/3;
judgmentRDM.RDM(49:68:end) = 1/3;
judgmentRDM.RDM(50:68:end) = 1/3;
judgmentRDM.RDM(51:68:end) = 1/3;

judgmentRDM.RDM = 1-judgmentRDM.RDM;
judgmentRDM.name = 'vowels only';

judgments=cell(1,4);
judgments{1}=judgmentRDM;

% figure
% imagesc(judgments{1}.RDM)

%Squares based on vowels and onset consonants
judgments{2}.RDM = zeros(16,16);
judgments{2}.RDM(1:17:end) = 1;
judgments{2}.RDM(2:68:end) = 2/3;
judgments{2}.RDM(3:68:end) = 1/3;
judgments{2}.RDM(4:68:end) = 1/3;
judgments{2}.RDM(17:68:end) = 2/3;
judgments{2}.RDM(19:68:end) = 1/3;
judgments{2}.RDM(20:68:end) = 1/3;
judgments{2}.RDM(33:68:end) = 1/3;
judgments{2}.RDM(34:68:end) = 1/3;
judgments{2}.RDM(36:68:end) = 2/3;
judgments{2}.RDM(49:68:end) = 1/3;
judgments{2}.RDM(50:68:end) = 1/3;
judgments{2}.RDM(51:68:end) = 2/3;

judgments{2}.RDM(15,3) = 1/3;
judgments{2}.RDM(15,4) = 1/3;
judgments{2}.RDM(16,4) = 1/3;
judgments{2}.RDM(16,3) = 1/3;
judgments{2}.RDM(3,16) = 1/3;
judgments{2}.RDM(4,16) = 1/3;
judgments{2}.RDM(4,15) = 1/3;
judgments{2}.RDM(3,15) = 1/3;

judgments{2}.RDM = 1-judgments{2}.RDM;
judgments{2}.name = 'onset_vowel';

% figure
% imagesc(judgments{2}.RDM)

%Squares based on vowels and offset consonants

judgments{3}.RDM = zeros(16,16);
judgments{3}.RDM(1:17:end) = 1;
judgments{3}.RDM(2:68:end) = 1/3;
judgments{3}.RDM(3:68:end) = 2/3;
judgments{3}.RDM(4:68:end) = 1/3;
judgments{3}.RDM(17:68:end) = 1/3;
judgments{3}.RDM(19:68:end) = 1/3;
judgments{3}.RDM(20:68:end) = 2/3;
judgments{3}.RDM(33:68:end) = 2/3;
judgments{3}.RDM(34:68:end) = 1/3;
judgments{3}.RDM(36:68:end) = 1/3;
judgments{3}.RDM(49:68:end) = 1/3;
judgments{3}.RDM(50:68:end) = 2/3;
judgments{3}.RDM(51:68:end) = 1/3;

judgments{3}.RDM(1,16) = 1/3;
judgments{3}.RDM(1,14) = 1/3;
judgments{3}.RDM(16,1) = 1/3;
judgments{3}.RDM(14,1) = 1/3;
judgments{3}.RDM(3,16) = 1/3;
judgments{3}.RDM(3,14) = 1/3;
judgments{3}.RDM(16,3) = 1/3;
judgments{3}.RDM(14,3) = 1/3;

judgments{3}.RDM(5,9) = 1/3;
judgments{3}.RDM(7,9) = 1/3;
judgments{3}.RDM(9,5) = 1/3;
judgments{3}.RDM(9,7) = 1/3;
judgments{3}.RDM(5,11) = 1/3;
judgments{3}.RDM(7,11) = 1/3;
judgments{3}.RDM(11,5) = 1/3;
judgments{3}.RDM(11,7) = 1/3;

judgments{3}.RDM(6,10) = 1/3;
judgments{3}.RDM(8,10) = 1/3;
judgments{3}.RDM(10,6) = 1/3;
judgments{3}.RDM(10,8) = 1/3;
judgments{3}.RDM(6,12) = 1/3;
judgments{3}.RDM(8,12) = 1/3;
judgments{3}.RDM(12,6) = 1/3;
judgments{3}.RDM(12,8) = 1/3;

judgments{3}.RDM = 1-judgments{3}.RDM;
judgments{3}.name = 'offset_vowel';

% figure
% imagesc(judgments{3}.RDM)

%Squares based on all shared features

judgments{4}.RDM = zeros(16,16);
judgments{4}.RDM(1:17:end) = 1;
judgments{4}.RDM(2:68:end) = 2/3;
judgments{4}.RDM(3:68:end) = 2/3;
judgments{4}.RDM(4:68:end) = 1/3;
judgments{4}.RDM(17:68:end) = 2/3;
judgments{4}.RDM(19:68:end) = 1/3;
judgments{4}.RDM(20:68:end) = 2/3;
judgments{4}.RDM(33:68:end) = 2/3;
judgments{4}.RDM(34:68:end) = 1/3;
judgments{4}.RDM(36:68:end) = 2/3;
judgments{4}.RDM(49:68:end) = 1/3;
judgments{4}.RDM(50:68:end) = 2/3;
judgments{4}.RDM(51:68:end) = 2/3;

judgments{4}.RDM(1,16) = 1/3;
judgments{4}.RDM(1,14) = 1/3;
judgments{4}.RDM(16,1) = 1/3;
judgments{4}.RDM(14,1) = 1/3;
judgments{4}.RDM(3,16) = 1/3;
judgments{4}.RDM(3,14) = 1/3;
judgments{4}.RDM(16,3) = 1/3;
judgments{4}.RDM(14,3) = 1/3;

judgments{4}.RDM(5,9) = 1/3;
judgments{4}.RDM(7,9) = 1/3;
judgments{4}.RDM(9,5) = 1/3;
judgments{4}.RDM(9,7) = 1/3;
judgments{4}.RDM(5,11) = 1/3;
judgments{4}.RDM(7,11) = 1/3;
judgments{4}.RDM(11,5) = 1/3;
judgments{4}.RDM(11,7) = 1/3;

judgments{4}.RDM(6,10) = 1/3;
judgments{4}.RDM(8,10) = 1/3;
judgments{4}.RDM(10,6) = 1/3;
judgments{4}.RDM(10,8) = 1/3;
judgments{4}.RDM(6,12) = 1/3;
judgments{4}.RDM(8,12) = 1/3;
judgments{4}.RDM(12,6) = 1/3;
judgments{4}.RDM(12,8) = 1/3;

judgments{4}.RDM(15,3) = 1/3;
judgments{4}.RDM(15,4) = 1/3;
judgments{4}.RDM(16,4) = 1/3;
judgments{4}.RDM(16,3) = 1/3;
judgments{4}.RDM(3,16) = 1/3;
judgments{4}.RDM(4,16) = 1/3;
judgments{4}.RDM(4,15) = 1/3;
judgments{4}.RDM(3,15) = 1/3;

judgments{4}.RDM = 1-judgments{4}.RDM;
judgments{4}.name = 'all_features';

% figure
% imagesc(judgments{4}.RDM)

% now we make an RDM per session
RDMs=struct();

for sessionI=1 %At the moment, only implements average
    % the correlation distance patterns are computed using the pdist function
    %First create t patterns from mask
    mask = spm_read_vols(spm_vol(mask_path));
    for pat = 1:length(tpattern_numbers)
        raw_pattern = spm_read_vols(spm_vol([data_path 'spmT_' sprintf('%04d',tpattern_numbers(pat)) '.nii']));
        if min(size(mask) == size(raw_pattern)) ~= 1
            error('Your mask and t-patterns are not in the same space, try again')
        end
        thesepatterns(pat,:) = raw_pattern(mask~=0)';
    end
        
    RDMs(1,sessionI).RDM   = squareform(pdist(thesepatterns,'correlation'));
    RDMs(1,sessionI).name = condition_name; 
    %RDMs(1,sessionI).name  = sprintf('Univariate Mask | session %d | condition %s excluding %s',sessionI,condition_order{condJ},num2str(find(all_empty_conditions)));
    RDMs(1,sessionI).color = [];
end

% % show the 2 session RDMs
% figI=condj;
% figure(figI);clf;set(gcf,'Position',[100 100 800 800],'Color','w')
% showRDMs(RDMs,figI);

% % ---------------------------------------------------------------------
% % compare RDMs across sessions
% r12=corr([vectorizeRDMs(RDMs(1).RDM)]',[vectorizeRDMs(RDMs(2).RDM)]');
% % ---------------------------------------------------------------------

% avgRDM = averageRDMs_subjectSession(RDMs,'subject');
avgRDM = RDMs(1,1); %At the moment only average is implemented
%avgRDM.name=sprintf('RDM across sessions | condition %s',condition_order{condJ});

figI = condJ;
figure(figI);set(gcf,'Position',[100 100 800 800],'Color','w')
showRDMs(avgRDM,figI);


% define the labels and indices for familiar and unfamiliar images
reductionLabels = {'Vowel1','Vowel2','Vowel3','Vowel4'};
nobjects = length(tpattern_numbers);

%reduction = reductionvectors;

Vowel1 = 1:4; Vowel2 = 5:8; Vowel3 = 9:12; Vowel4 = 13:16;

nCols=4;
cmap=RDMcolormap;
colors=cmap([1 80 155 222],:);
options.categoryColors=zeros(length(tpattern_numbers),3);
options.categoryColors(Vowel1,:)=repmat(colors(1,:),length(Vowel1),1);
options.categoryColors(Vowel2,:)=repmat(colors(2,:),length(Vowel2),1);
options.categoryColors(Vowel3,:)=repmat(colors(3,:),length(Vowel3),1);
options.categoryColors(Vowel4,:)=repmat(colors(4,:),length(Vowel4),1);

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

userOptions.figureIndex = [20+condJ, 30+condJ];
stats_p_r=compareRefRDM2candRDMs(avgRDM, judgments, userOptions);

% ---------------------------------------------------------------------
% add a second model



% ---------------------------------------------------------------------