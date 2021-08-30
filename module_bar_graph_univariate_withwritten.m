addpath ./plotting
masknames = {
       'Left_Precentral_Univariate_Interaction_combined';
     'Left_STG_Univariate3mm_15>3'
    'Left_Frontal_Univariate_MM>M'
    'Left_Temporal_Univariate_MM>M'
    'Left_Angular_Univariate_Interaction_combined'
    };
secondlevelnames = {
    'Clarity Congruency Interaction Positive';
    };
if ~exist('./univariate_bars')
    mkdir('./univariate_bars')
end

for this_mask = 1:length(masknames)
    for this_second_level = 1:length(secondlevelnames)
        Y = spm_read_vols(spm_vol(['./atlas_Neuromorphometrics/' masknames{this_mask} '.nii']),1);
        indx = find(Y>0);
        [x,y,z] = ind2sub(size(Y),indx);
        XYZ = [x y z]';
        
        load(['/imaging/mlr/users/tc02/PINFA_preprocessed_2021/stats6_3/' secondlevelnames{this_second_level} '/SPM.mat'])
        written_SPM = load(['/imaging/mlr/users/tc02/PINFA_preprocessed_2021/stats6_3/Written > Normal/SPM.mat']);
        all_data = [];
        all_subj_weighted_activations = [];
        all_subj_unweighted_activations = [];
        for crun = 1:size(subjects,2)
            disp(['Working on subject ' subjects{crun}]);
            this_split = strsplit(SPM.xY.P{crun},'con_');
            this_SPM_dir = this_split{1};
            this_contrast_num = str2num(this_split{2}(1:end-4));
            this_SPM = load([this_SPM_dir 'SPM.mat']);
            num_sess = length(this_SPM.SPM.Sess);
            num_conditions = (length(this_SPM.SPM.xCon(this_contrast_num).c) - num_sess)/num_sess;
            num_active_conditions = sum((this_SPM.SPM.xCon(this_contrast_num).c(1:num_conditions)~=0));
            if crun == 1
                template_active_conditions = num_active_conditions;
            else
                assert(num_active_conditions==template_active_conditions,'The contrasts do not seem to be specified the same way in all files')
            end
            for i = 1:num_sess
                k = 0;
                for j = 1:num_conditions
                    if this_SPM.SPM.xCon(this_contrast_num).c((i-1)*num_conditions+j)~=0
                        k = k+1;
                        all_data(crun,i,k) = mean(spm_get_data([this_SPM_dir this_SPM.SPM.Vbeta((i-1)*num_conditions+j).fname], XYZ),2);
                    end
                end
            end
            
            this_split = strsplit(written_SPM.SPM.xY.P{crun},'con_');
            this_SPM_dir = this_split{1};
            this_contrast_num = str2num(this_split{2}(1:end-4));
            this_SPM = load([this_SPM_dir 'SPM.mat']);
            num_sess = length(this_SPM.SPM.Sess);
            num_conditions = (length(this_SPM.SPM.xCon(this_contrast_num).c) - num_sess)/num_sess;
            for i = 1:num_sess
                for j = 1:num_conditions
                    if this_SPM.SPM.xCon(this_contrast_num).c((i-1)*num_conditions+j)>0
                        all_data(crun,i,k+1) = mean(spm_get_data([this_SPM_dir this_SPM.SPM.Vbeta((i-1)*num_conditions+j).fname], XYZ),2);
                    end
                end
            end
            
            all_subj_weighted_activations(crun,:) = mean(normr(squeeze(all_data(crun,1:num_sess,:))));
            all_subj_unweighted_activations(crun,:) = mean(squeeze(all_data(crun,1:num_sess,:)));
        end
        
        figure
        %         barweb([mean(all_subj_weighted_activations(group==1,:));mean(all_subj_weighted_activations(group==2,:))],[std(all_subj_weighted_activations(group==1,:))/sqrt(sum(group==1));std(all_subj_weighted_activations(group==1,:))/sqrt(sum(group==1))],[],{'Controls','Patients'},['Extracted normalised betas in ' masknames{this_mask}],[],'Weighted beta value',[],[],{'Match 3';'Match 15';'MisMatch 3';'MisMatch 15';'Written'}) ;
        %Reorder to match behaviour
        all_subj_weighted_activations = all_subj_weighted_activations(:,[1,3,2,4,5]);
        barweb([mean(all_subj_weighted_activations(group==1,:));mean(all_subj_weighted_activations(group==2,:))],[std(all_subj_weighted_activations(group==1,:))/sqrt(sum(group==1));std(all_subj_weighted_activations(group==2,:))/sqrt(sum(group==2))],[],{'Controls','Patients'},['Extracted normalised betas in ' masknames{this_mask}],[],'Weighted beta value',[],[],{'Match 3';'MisMatch 3';'Match 15';'MisMatch 15';'Written'}) ;
        ylim([(min(mean(all_subj_weighted_activations))-4*max(std(all_subj_weighted_activations)/sqrt(sum(group==2)))),(max(mean(all_subj_weighted_activations))+4*max(std(all_subj_weighted_activations)/sqrt(sum(group==2))))])
        saveas(gcf, ['./univariate_bars/Extracted normalised betas in ' masknames{this_mask} '.pdf']);
        saveas(gcf, ['./univariate_bars/Extracted normalised betas in ' masknames{this_mask} '.png']);
        
        figure
        %         barweb([mean(all_subj_unweighted_activations(group==1,:));mean(all_subj_unweighted_activations(group==2,:))],[std(all_subj_unweighted_activations(group==1,:))/sqrt(sum(group==1));std(all_subj_unweighted_activations(group==1,:))/sqrt(sum(group==1))],[],{'Controls','Patients'},['Extracted betas in ' masknames{this_mask}],[],'Unweighted beta value',[],[],{'Match 3';'Match 15';'MisMatch 3';'MisMatch 15';'Written'}) ;
        %Reorder to match behaviour
        all_subj_unweighted_activations = all_subj_unweighted_activations(:,[1,3,2,4,5]);
        barweb([mean(all_subj_unweighted_activations(group==1,:));mean(all_subj_unweighted_activations(group==2,:))],[std(all_subj_unweighted_activations(group==1,:))/sqrt(sum(group==1));std(all_subj_unweighted_activations(group==2,:))/sqrt(sum(group==2))],[],{'Controls','Patients'},['Extracted betas in ' masknames{this_mask}],[],'Unweighted beta value',[],[],{'Match 3';'MisMatch 3';'Match 15';'MisMatch 15';'Written'}) ;
        ylim([(min(mean(all_subj_unweighted_activations))-4*max(std(all_subj_unweighted_activations)/sqrt(sum(group==2)))),(max(mean(all_subj_unweighted_activations))+4*max(std(all_subj_unweighted_activations)/sqrt(sum(group==2))))])
        saveas(gcf, ['./univariate_bars/Extracted betas in ' masknames{this_mask} '.pdf']);
        saveas(gcf, ['./univariate_bars/Extracted betas in ' masknames{this_mask} '.png']);
    end
end