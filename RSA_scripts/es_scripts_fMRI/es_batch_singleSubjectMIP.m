addpath /group/language/data/ediz.sohoglu/matlab/MIP-C

clearvars
close all
clc

spm fmri

GLMAnalPD = '/imaging/es03/fMRI_2017/GLMAnalysisMNISmooth6mm';

Subj = {'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24' 'subj25' 'subj26' 'subj27'};

SubjToAnalyze = [3,4,5,6,8,10,11,14,15,16,17,18,19,20,21,22,23,24,26,27];
ContrastsToAnalyse = [9 16];

thresh = .001; % Uncorrected

for k=1:length(SubjToAnalyze)
    
    for c=1:length(ContrastsToAnalyse)
                
        SubjCurrent=Subj{SubjToAnalyze(k)};
                
        GLMAnalDir =[GLMAnalPD '/' SubjCurrent];
                
        fmri_mip_color(fullfile(GLMAnalDir,'SPM.mat'), ContrastsToAnalyse(c), thresh, 'none');
        
        load(fullfile(GLMAnalDir,'SPM.mat'));
        contrastName = strrep(SPM.xCon(ContrastsToAnalyse(c)).name,' - All Sessions','');
        contrastName = strrep(contrastName,' ','');
        
        title([SubjCurrent ' ' contrastName]);
        saveas(gcf,[ SubjCurrent '_' contrastName],'pdf');
        close(gcf);
        
    end
    
end

