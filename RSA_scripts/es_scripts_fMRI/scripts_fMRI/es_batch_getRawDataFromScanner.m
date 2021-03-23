clearvars

AnalysisFolder = '/imaging/es03/fMRI_2017';
Subj = {'subj4' 'subj5' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24' 'subj25' 'subj26' 'subj27'};
RawDataFolders = {
    '/mridata/cbu/CBU170700_MR17011/20171025_102118'... % 'subj4'
    '/mridata/cbu/CBU170713_MR17011/20171027_134934'... % 'subj5'
    '/mridata/cbu/CBU170721_MR17011/20171101_101145'... % 'subj7'
    '/mridata/cbu/CBU170727_MR17011/20171102_134813'... % 'subj8'
    '/mridata/cbu/CBU170729_MR17011/20171103_114323'... % 'subj9'
    '/mridata/cbu/CBU170737_MR17011/20171106_133716'... % 'subj10'
    '/mridata/cbu/CBU170746_MR17011/20171108_102003'... % 'subj11'
    '/mridata/cbu/CBU170751_MR17011/20171109_134345'... % 'subj12'
    '/mridata/cbu/CBU170755_MR17011/20171110_121958'... % 'subj13'
    '/mridata/cbu/CBU170761_MR17011/20171113_133519'... % 'subj14'
    '/mridata/cbu/CBU170764_MR17011/20171114_130021'... % 'subj15'
    '/mridata/cbu/CBU170769_MR17011/20171115_121829'... % 'subj16'
    '/mridata/cbu/CBU170772_MR17011/20171116_115635'... % 'subj17'
    '/mridata/cbu/CBU170777_MR17011/20171117_134332'... % 'subj18'
    '/mridata/cbu/CBU170780_MR17011/20171120_120728'... % 'subj19'
    '/mridata/cbu/CBU170785_MR17011/20171121_115917'... % 'subj20'
    '/mridata/cbu/CBU170797_MR17011/20171124_115017'... % 'subj21'
    '/mridata/cbu/CBU170801_MR17011/20171127_115951'... % 'subj22'
    '/mridata/cbu/CBU170823_MR17011/20171201_121255'... % 'subj23'
    '/mridata/cbu/CBU170827_MR17011/20171204_122259'... % 'subj24'
    '/mridata/cbu/CBU170835_MR17011/20171206_122644'... % 'subj25'
    '/mridata/cbu/CBU170857_MR17011/20171211_142317'... % 'subj26'
    '/mridata/cbu/CBU170862_MR17011/20171212_134545'    % 'subj27'  
    };
Ndummies = 4;

SubjToAnalyze = [21];

for k=1:length(SubjToAnalyze)
    RawDataSubjectFolder = fullfile(RawDataFolders{SubjToAnalyze(k)});
    SubFolders = dir(RawDataSubjectFolder);
    SubFolders = SubFolders(3:end); % Dont include first '.' and '..' folders from dir()
    for ii=1:length(SubFolders)
        Files = dir(fullfile(RawDataSubjectFolder,SubFolders (ii).name,'*.dcm'));
        Hdr = {};
        for iii=1:length(Files)
            Hdr(iii) = spm_dicom_headers(fullfile(RawDataSubjectFolder,SubFolders (ii).name,Files(iii).name));
        end
        if strfind(SubFolders (ii).name,'run1') & length(Files) > 280
            Hdr = Hdr(1+Ndummies:end);
            OutputFolder = fullfile(AnalysisFolder,'PreprocessAnalysis',Subj{SubjToAnalyze(k)},'Functional','Sess_1');
            if ~exist(OutputFolder,'dir');
                mkdir(OutputFolder);
            else
                error('Folder %s for %s already exists!',OutputFolder,Subj{SubjToAnalyze(k)});
            end
            spm_dicom_convert(Hdr,'all','flat','nii',OutputFolder);
        elseif strfind(SubFolders (ii).name,'run2') & length(Files) > 280
            Hdr = Hdr(1+Ndummies:end);
            OutputFolder = fullfile(AnalysisFolder,'PreprocessAnalysis',Subj{SubjToAnalyze(k)},'Functional','Sess_2');
            if ~exist(OutputFolder,'dir');
                mkdir(OutputFolder);
            else
                error('Folder %s for %s already exists!',OutputFolder,Subj{SubjToAnalyze(k)});
            end
            spm_dicom_convert(Hdr,'all','flat','nii',OutputFolder);
        elseif strfind(SubFolders (ii).name,'run3') & length(Files) > 280
            Hdr = Hdr(1+Ndummies:end);
            OutputFolder = fullfile(AnalysisFolder,'PreprocessAnalysis',Subj{SubjToAnalyze(k)},'Functional','Sess_3');
            if ~exist(OutputFolder,'dir');
                mkdir(OutputFolder);
            else
                error('Folder %s for %s already exists!',OutputFolder,Subj{SubjToAnalyze(k)});
            end
            spm_dicom_convert(Hdr,'all','flat','nii',OutputFolder);
        elseif strfind(SubFolders (ii).name,'run4') & length(Files) > 280
            Hdr = Hdr(1+Ndummies:end);
            OutputFolder = fullfile(AnalysisFolder,'PreprocessAnalysis',Subj{SubjToAnalyze(k)},'Functional','Sess_4');
            if ~exist(OutputFolder,'dir');
                mkdir(OutputFolder);
            else
                error('Folder %s for %s already exists!',OutputFolder,Subj{SubjToAnalyze(k)});
            end
            spm_dicom_convert(Hdr,'all','flat','nii',OutputFolder);
        elseif strfind(SubFolders (ii).name,'run5') & length(Files) > 280
            Hdr = Hdr(1+Ndummies:end);
            OutputFolder = fullfile(AnalysisFolder,'PreprocessAnalysis',Subj{SubjToAnalyze(k)},'Functional','Sess_5');
            if ~exist(OutputFolder,'dir');
                mkdir(OutputFolder);
            else
                error('Folder %s for %s already exists!',OutputFolder,Subj{SubjToAnalyze(k)});
            end
            spm_dicom_convert(Hdr,'all','flat','nii',OutputFolder);
        elseif strfind(SubFolders (ii).name,'MPRAGE')
            OutputFolder = fullfile(AnalysisFolder,'PreprocessAnalysis',Subj{SubjToAnalyze(k)},'Structural');
            if ~exist(OutputFolder,'dir');
                mkdir(OutputFolder);
            else
                error('Folder %s for %s already exists!',OutputFolder,Subj{SubjToAnalyze(k)});
            end
            spm_dicom_convert(Hdr,'all','flat','nii',OutputFolder);
        elseif strfind(SubFolders (ii).name,'fieldmap')
            if length(Files)==96
                OutputFolder = fullfile(AnalysisFolder,'PreprocessAnalysis',Subj{SubjToAnalyze(k)},'FieldMaps','Magnitude');
                if ~exist(OutputFolder,'dir');
                    mkdir(OutputFolder);
                else
                    error('Folder %s for %s already exists!',OutputFolder,Subj{SubjToAnalyze(k)});
                end
                spm_dicom_convert(Hdr,'all','flat','nii',OutputFolder);
            elseif length(Files)==48
                OutputFolder = fullfile(AnalysisFolder,'PreprocessAnalysis',Subj{SubjToAnalyze(k)},'FieldMaps','Phase');
                if ~exist(OutputFolder,'dir');
                    mkdir(OutputFolder);
                else
                    error('Folder %s for %s already exists!',OutputFolder,Subj{SubjToAnalyze(k)});
                end
                spm_dicom_convert(Hdr,'all','flat','nii',OutputFolder);
            end
        end
    end
end
