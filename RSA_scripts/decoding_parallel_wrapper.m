function results = decoding_parallel_wrapper(cfg,misc,searchlights_per_worker,worker_number)
cfg.searchlight.subset = ((worker_number-1)*searchlights_per_worker)+1:worker_number*searchlights_per_worker;
cfg.results.resultsname = cellstr(['parallel_loop_' num2str(worker_number)]);
addpath /group/language/data/thomascope/spm12_fil_r6906/ % Your SPM path for the workers
spm('ver'); % Needed or sometimes the decoding toolbox complains in parallel that SPM is not initialised.
results = decoding(cfg,[],misc);