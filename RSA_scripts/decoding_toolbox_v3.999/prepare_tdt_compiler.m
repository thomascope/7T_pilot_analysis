function prepare_tdt_compiler

% This function can be added to your script that you want to compile if you
% want to make sure the Matlab compiler compiles all relevant files in TDT.
% It also contains details on how to compile. It is a workaround calling
% all relevant files in TDT once manually, which should invoke their
% inclusion. However, it seems to be failsafe as compared to other methods
% of compilation.
%
% Martin Hebart 2018/07/31

% first and only check: is decoding_defaults on the path
try
    decoding_defaults;
catch 
    error('Please add TDT to the Matlab path before running this preparation script')
end

%% First of all, since you may need SPM, you need to make sure your version
% of it on the server contains a 'Contents.txt' file (SPM needs this to run
% on the server in compiled mode). While on the server, execute the
% following code (commented here):
% sts = copyfile(fullfile(spm('Dir'),'Contents.m'),...
%                fullfile(spm('Dir'),'Contents.txt'));

%% Then the usual approach would be to call mcc from the command line and
% add all relevant subfolders to compilation, possibly both for SPM and for
% TDT. An example would be the following (can be made a shell script by
% replacing the pairs in the input, e.g. replace '-a','tdt' by -a tdt .

% scriptname = 'myscript.m'; % this is the name of your script to be compiled, perhaps the full path
% tdt_dir = 'location_of_tdt'; % if on path, use tdt_dir = fileparts(decoding_defaults);
% spm_dir = 'location_of_spm'; % if on path, use spm_dir = spm('dir');
% mcc('mcc -C -mv -R -singleCompThread -R -nodisplay',...
%  '-I',tdt_dir,...
%  '-I',fullfile(tdt_dir,'checks'),...
%  '-I',fullfile(tdt_dir,'feature_transformation'),...
%  '-I',fullfile(tdt_dir,'utils'),...
%  '-I',fullfile(tdt_dir,'general'),...
%  '-I',fullfile(tdt_dir,'decoding_software'),...
%  '-I',fullfile(tdt_dir,'decoding_software/libsvm3.17'),...
%  '-I',fullfile(tdt_dir,'decoding_software/libsvm3.17/matlab'),...
%  '-I',fullfile(tdt_dir,'decoding_software/newton_svm'),...
%  '-I',fullfile(tdt_dir,'decoding_software/lda'),...
%  '-I',fullfile(tdt_dir,'decoding_software/correlation_classifier'),...
%  '-I',fullfile(tdt_dir,'decoding_software/liblinear1.94'),...
%  '-I',fullfile(tdt_dir,'decoding_software/liblinear1.94/blas'),...
%  '-I',fullfile(tdt_dir,'decoding_software/liblinear1.94/matlab'),...
%  '-I',fullfile(tdt_dir,'decoding_software/ensemble'),...
%  '-I',fullfile(tdt_dir,'decoding_software/pattern_similarity'),...
%  '-I',fullfile(tdt_dir,'transform_results'),...
%  '-I',fullfile(tdt_dir,'transform_results/subfunctions'),...
%  '-I',fullfile(tdt_dir,'parameter_selection'),...
%  '-I',fullfile(tdt_dir,'statistics'),...
%  '-I',fullfile(tdt_dir,'statistics/prevalence_inference'),...
%  '-I',fullfile(tdt_dir,'demos'),...
%  '-I',fullfile(tdt_dir,'demos/demo_haxby2001'),...
%  '-I',fullfile(tdt_dir,'templates'),...
%  '-I',fullfile(tdt_dir,'design'),...
%  '-I',fullfile(tdt_dir,'access_images'),...
%  '-I',fullfile(tdt_dir,'access_images/spm12'),...
%  '-I',fullfile(tdt_dir,'access_images/afni'),...
%  '-I',fullfile(tdt_dir,'access_images/spm8'),...
%  '-I',fullfile(tdt_dir,'access_images/spm2'),...
%  '-I',fullfile(tdt_dir,'access_images/spm8b'),...
%  '-I',fullfile(tdt_dir,'access_images/spm5'),...
%  '-I',fullfile(tdt_dir,'access_images/spm12b'),...
%  '-I',fullfile(tdt_dir,'feature_selection'),...
%  '-I',fullfile(tdt_dir,'feature_selection/general'),...
%  '-I',fullfile(tdt_dir,'feature_selection/embedded_subfunctions'),...
%  '-I',fullfile(tdt_dir,'feature_selection/filter_subfunctions'),...
%  '-I',fullfile(tdt_dir,'visualization'),...
%  '-I',spm_dir,...
%  '-I',fullfile(spm_dir,'spm_orthviews'),...
%  '-I',fullfile(spm_dir,'config'),...
%  '-I',fullfile(spm_dir,'matlabbatch'),...
%  '-I',fullfile(spm_dir,'matlabbatch/cfg_confgui'),...
%  '-I',fullfile(spm_dir,'matlabbatch/cfg_basicio'),...
%  '-I',fullfile(spm_dir,'matlabbatch/cfg_basicio/src'),...
%  '-I',fullfile(spm_dir,'matlabbatch/cfg_basicio/src/02_var_ops'),...
%  '-I',fullfile(spm_dir,'matlabbatch/cfg_basicio/src/03_run_ops'),...
%  '-I',fullfile(spm_dir,'matlabbatch/cfg_basicio/src/01_file_dir_ops'),...
%  '-I',fullfile(spm_dir,'matlabbatch/cfg_basicio/src/01_file_dir_ops/01_dir_ops'),...
%  '-I',fullfile(spm_dir,'matlabbatch/cfg_basicio/src/01_file_dir_ops/02_file_ops'),...
%  '-I',fullfile(spm_dir,'matlabbatch/examples'),...
%  '-I',fullfile(spm_dir,'src'),...
%  '-I',fullfile(spm_dir,'toolbox'),...
%  '-I',fullfile(spm_dir,'compat'),...
%  '-I',fullfile(spm_dir,'help'),...
%  '-I',fullfile(spm_dir,'help/images'),...
%  '-I',fullfile(spm_dir,'EEGtemplates'),...
%  '-I',fullfile(spm_dir,'canonical'),...
%  '-I',fullfile(spm_dir,'batches'),...
%  '-I',fullfile(spm_dir,'rend'),...
%  '-I',fullfile(spm_dir,'external'),...
%  '-I',fullfile(spm_dir,'tpm'),...
%  scriptname);

%% Now run actual call

cfg.verbose = 0;
try decoding(cfg); end %#ok<*TRYNC>
try decoding_create_maskfile(cfg); end
try decoding_defaults(cfg); end
try decoding_describe_data(cfg); end
try decoding_example(cfg); end
try decoding_example_afni; end
try decoding_feature_selection(cfg); end
try decoding_feature_transformation(cfg); end
try decoding_generate_output(cfg); end
try decoding_load_data(cfg); end
try decoding_load_misc(cfg); end
try decoding_parameter_selection(cfg); end
try decoding_prepare_feature_selection(cfg); end
try decoding_prepare_searchlight(cfg); end
try decoding_scale_data(cfg); end
try decoding_statistics(cfg); end
try decoding_tutorial(cfg); end
try decoding_write_results(cfg); end
try decoding_write_similarity(cfg); end
try dispv(cfg); end
try fill_passed_data(cfg); end
try genpath_nosvn(cfg); end
try get_ind(cfg); end
try get_n_decodings(cfg); end
try inherit_settings(cfg); end
try load_mask(cfg); end
try report_files(cfg); end
try residuals_from_spm(cfg); end
try warningv(cfg); end
try check_software(cfg); end
try compare_volumes(cfg); end
try get_decoding_out_from_passed_data(cfg); end
try get_filenames(cfg); end
try read_header(cfg); end
try read_image(cfg); end
try read_resultdata(cfg); end
try read_voxels(cfg); end
try write_image(cfg); end
try check_software_afni(cfg); end
try get_filenames_afni(cfg); end
try read_header_afni(cfg); end
try read_image_afni(cfg); end
try read_voxels_afni(cfg); end
try write_image_afni(cfg); end
try check_software_spm12(cfg); end
try get_filenames_spm12(cfg); end
try read_header_spm12(cfg); end
try read_image_spm12(cfg); end
try read_voxels_spm12(cfg); end
try write_image_spm12(cfg); end
try check_software_spm12b(cfg); end
try get_filenames_spm12b(cfg); end
try read_header_spm12b(cfg); end
try read_image_spm12b(cfg); end
try read_voxels_spm12b(cfg); end
try write_image_spm12b(cfg); end
try check_software_spm2(cfg); end
try get_filenames_spm2(cfg); end
try read_header_spm2(cfg); end
try read_image_spm2(cfg); end
try read_voxels_spm2(cfg); end
try write_image_spm2(cfg); end
try check_software_spm5(cfg); end
try get_filenames_spm5(cfg); end
try read_header_spm5(cfg); end
try read_image_spm5(cfg); end
try read_voxels_spm5(cfg); end
try write_image_spm5(cfg); end
try check_software_spm8(cfg); end
try get_filenames_spm8(cfg); end
try read_header_spm8(cfg); end
try read_image_spm8(cfg); end
try read_voxels_spm8(cfg); end
try write_image_spm8(cfg); end
try check_software_spm8b(cfg); end
try get_filenames_spm8b(cfg); end
try read_header_spm8b(cfg); end
try read_image_spm8b(cfg); end
try read_voxels_spm8b(cfg); end
try write_image_spm8b(cfg); end
try check_datatrans(cfg); end
try check_libsvm(cfg); end
try check_result_fits_cfg(cfg); end
try check_verbosity(cfg); end
try decoding_basic_checks(cfg); end
try correlation_classifier_test(cfg); end
try correlation_classifier_train(cfg); end
try distance_test(cfg); end
try distance_train(cfg); end
try ensemble_balance_test(cfg); end
try ensemble_balance_train(cfg); end
try lda_test(cfg); end
try lda_train(cfg); end
try liblinear_test(cfg); end
try liblinear_train(cfg); end
try libsvm_test(cfg); end
try libsvm_train(cfg); end
try newton_test(cfg); end
try newton_train(cfg); end
try passdata_test(cfg); end
try passdata_train(cfg); end
try similarity_test(cfg); end
try similarity_train(cfg); end
try correlation_classifier(cfg); end
try ens_bal_te(cfg); end
try ens_bal_tr(cfg); end
try ldapredict(cfg); end
try ldatrain(cfg); end
try nsvm_test(cfg); end
try nsvm_train(cfg); end
try pattern_similarity(cfg); end
try pattern_similarity_fast(cfg); end
try combine_designs(cfg); end
try design_from_afni(cfg); end
try design_from_spm(cfg); end
try make_design_alldata(cfg); end
try make_design_boot(cfg); end
try make_design_boot_cv(cfg); end
try make_design_cv(cfg); end
try make_design_permutation(cfg); end
try make_design_rsa(cfg); end
try make_design_rsa_cv(cfg); end
try make_design_separate(cfg); end
try make_design_similarity(cfg); end
try make_design_similarity_cv(cfg); end
try make_design_xclass(cfg); end
try make_design_xclass_cv(cfg); end
try sort_design(cfg); end
try RFEselect_test(cfg); end
try RFEselect_train(cfg); end
try feature_selection_embedded(cfg); end
try feature_selection_filter(cfg); end
try RFE(cfg); end
try eget(cfg); end
try fget(cfg); end
try fvalue(cfg); end
try uget(cfg); end
try wget(cfg); end
try transfeat_PCA(cfg); end
try correl(cfg); end
try correlmat(cfg); end
try cov2(cfg); end
try covq(cfg); end
try covshrink_lw(cfg); end
try covshrink_lw2(cfg); end
try covshrink_oas(cfg); end
try cveuclidean2(cfg); end
try euclidean(cfg); end
try euclidean2(cfg); end
try interp1dec(cfg); end
try select_peak(cfg); end
try squareformq(cfg); end
try tdt_fileparts(cfg); end
try transform_vol(cfg); end
try uniqueq(cfg); end
try uniqueqs(cfg); end
try wildcard2regexp(cfg); end
try param_string_number(cfg); end
try f_cdf(cfg); end
try stats_binomial(cfg); end
try stats_permutation(cfg); end
try t_cdf(cfg); end
% try prevalenceCore(cfg); end
% try prevalenceTDT(cfg); end
% try prevalence_loaddata(cfg); end
% try prevalence_savedata_to_images(cfg); end
try saveMRImage(cfg); end
try decoding_transform_results(cfg); end
try transres_MSS_diff(cfg); end
try transres_SVM_pattern(cfg); end
try transres_SVM_pattern_alldata(cfg); end
try transres_SVM_weights(cfg); end
try transres_SVM_weights_plusbias(cfg); end
try transres_accuracy_matrix(cfg); end
try transres_accuracy_matrix_minus_chance(cfg); end
try transres_accuracy_pairwise(cfg); end
try transres_accuracy_pairwise_minus_chance(cfg); end
try transres_binomial_probability(cfg); end
try transres_confusion_matrix(cfg); end
try transres_model_parameters(cfg); end
try transres_ninputdim(cfg); end
try transres_other(cfg); end
try transres_other_average(cfg); end
try transres_other_meandist(cfg); end
try transres_rsa_beta(cfg); end
try transres_rsa_corr_kendall(cfg); end
try transres_rsa_corr_pearson(cfg); end
try transres_rsa_corr_spearman(cfg); end
try transres_signed_decision_values(cfg); end
try AUCstats(cfg); end
try AUCstats_matrix(cfg); end
try dprimestats(cfg); end
try tdtutil_change_voxelsize_or_bb(cfg); end
try tdtutil_change_voxelsize_or_bb_of_betas_and_reduce_SPM(cfg); end
try display_design(cfg); end
try display_progress(cfg); end
try display_regressor_names(cfg); end
try display_volume(cfg); end
try plot_design(cfg); end
try plot_selected_voxels(cfg); end
try save_fig(cfg); end