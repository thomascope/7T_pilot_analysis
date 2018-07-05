function module_topup_job(base_image_path, reversed_image_path, outpath, minvols, filestocorrect)

fsldir = '/imaging/local/software/fsl/fsl64/fsl-5.0.3/fsl';

setenv('FSLDIR',fsldir);  % this to tell where FSL folder is
setenv('FSLOUTPUTTYPE', 'NIFTI'); % this to tell what the output type 

tsize=minvols;

%First merge the two directions for topup deformation field creation
cmd = [fsldir '/bin/fslmerge -t ' outpath 'both_directions_fortopup ' base_image_path ' ' reversed_image_path ];
system(cmd)

%Now write txt file containing acquisition parameters (assumes y-direction
%flip for phase encoding)
base_head = spm_vol(base_image_path);
reversed_head = spm_vol(reversed_image_path);
if ~(abs(base_head.private.diminfo.slice_time.duration-reversed_head.private.diminfo.slice_time.duration) < 0.0001)
    error('Something went wrong - acquisition parameters are not the same in the two images')
end
epi_params = [0 -1 0 base_head.private.diminfo.slice_time.duration; 0 1 0 reversed_head.private.diminfo.slice_time.duration]; 
dlmwrite([outpath 'epi_acq_param.txt'],epi_params,'delimiter',' ')

%Then calculate topup (can take a longish time)
cmd = [fsldir '/bin/topup --imain=' outpath 'both_directions_fortopup --datain=' outpath 'epi_acq_param.txt --config=b02b0.cnf --out=' outpath 'epi_topup_results'];
disp('Calculating topup, this can take some time')
system(cmd)

%Now apply topup
disp('Calculations complete, now applying topup')
for i = 1:length(filestocorrect)
    splitpath = strsplit(filestocorrect{i}, '/');
    
    cmd = [fsldir '/bin/fslsplit ' filestocorrect{i} ' ' filestocorrect{i}(1:end-4) '_split -t'];
    system(cmd);
    
    for TR = 1:tsize
        vol = TR-1;
        
        cmd = [fsldir '/bin/applytopup --imain=' filestocorrect{i}(1:end-4) '_split' sprintf('%04d',vol) ' --inindex=1 --topup=' outpath 'epi_topup_results --datain=' outpath 'epi_acq_param.txt --method=jac --out=' filestocorrect{i}(1:end-length(splitpath{end})) 'topup_' splitpath{end}(1:end-4) '_split' sprintf('%04d',vol) ' --verbose'];
        system(cmd);
    end
    
    cmd = [fsldir '/bin/fslmerge -t ' outpath 'topup_' splitpath{end} ' ' filestocorrect{i}(1:end-length(splitpath{end})) 'topup_' splitpath{end}(1:end-4) '_split*'];
    system(cmd);
end
