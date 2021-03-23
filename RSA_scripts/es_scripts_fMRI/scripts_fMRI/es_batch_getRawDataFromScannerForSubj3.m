clearvars

AnalysisFolder = '/imaging/es03/fMRI_2017';
Subj = {'subj3'};
RawDataFolders = {'/mridata/cbu/CBU170676_MR17011/20171016_115459'};
Ndummies = 5;

SubjToAnalyze=[1];

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
        if strfind(SubFolders (ii).name,'run1')
            Hdr = Hdr(1+Ndummies:end);
            OutputFolder = fullfile(AnalysisFolder,'PreprocessAnalysis',Subj{SubjToAnalyze(k)},'Functional','Sess_1');
            if ~exist(OutputFolder,'dir');
                mkdir(OutputFolder);
            else
                error('Folder %s for %s already exists!',OutputFolder,Subj{SubjToAnalyze(k)});
            end
            spm_dicom_convert(Hdr,'all','flat','nii',OutputFolder);
        elseif strfind(SubFolders (ii).name,'run2')
            Hdr = Hdr(1+Ndummies:end);
            OutputFolder = fullfile(AnalysisFolder,'PreprocessAnalysis',Subj{SubjToAnalyze(k)},'Functional','Sess_2');
            if ~exist(OutputFolder,'dir');
                mkdir(OutputFolder);
            else
                error('Folder %s for %s already exists!',OutputFolder,Subj{SubjToAnalyze(k)});
            end
            spm_dicom_convert(Hdr,'all','flat','nii',OutputFolder);
        elseif strfind(SubFolders (ii).name,'run3')
            Hdr = Hdr(1+Ndummies:end);
            OutputFolder = fullfile(AnalysisFolder,'PreprocessAnalysis',Subj{SubjToAnalyze(k)},'Functional','Sess_3');
            if ~exist(OutputFolder,'dir');
                mkdir(OutputFolder);
            else
                error('Folder %s for %s already exists!',OutputFolder,Subj{SubjToAnalyze(k)});
            end
            spm_dicom_convert(Hdr,'all','flat','nii',OutputFolder);
        elseif strfind(SubFolders (ii).name,'run4')
            Hdr = Hdr(1+Ndummies:end);
            OutputFolder = fullfile(AnalysisFolder,'PreprocessAnalysis',Subj{SubjToAnalyze(k)},'Functional','Sess_4');
            if ~exist(OutputFolder,'dir');
                mkdir(OutputFolder);
            else
                error('Folder %s for %s already exists!',OutputFolder,Subj{SubjToAnalyze(k)});
            end
            spm_dicom_convert(Hdr,'all','flat','nii',OutputFolder);
        elseif strfind(SubFolders (ii).name,'run5')
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
