function crop_images(vols, no_images_to_crop)
% Coregister source images vols and apply to other images
% if present
% crop neck of source image
% Crop the first no_images_to_crop images and apply
% transformations to all.
%If cropping multiple images they should all have matching voxel size.
%But it may crop close enough.
%If the images do not match, e.g. T1 and T2 and you want to crop both
%then this function should be performed twice.

do_gunzip=0; %Gunzip the nifti files?
do_delete_sub=0;

if nargin < 2
    other='';
end

spm('defaults','pet');
% spm_jobman('initcfg');

%Change paths to match your setup
addpath(fileparts(mfilename('fullpath')));

%tp = fileparts(mfilename('fullpath'));

% subs={
% 'sub-25382' 'ses-20160111'
%     };

% source_image='mp2rage.nii';

%Only first image is used as a coregistration source
%Amend for partitipants with SWI. Don't crop the SWI.
%imagesanat={
%    '_acq-mprage_T1w.nii'
%    };

%imagesswi={
%    '_acq-mag_part-mag_GRE.nii'
%    '_acq-phase_part-phase_GRE.nii' 
%    };

if(1==exist('vols','var'))
    if ischar(vols)
        vols=cellstr(vols);
    end
    
    for vj=1:size(vols,1)
    
        [~,~,ext]=spm_fileparts(vols{vj,:});
        if(strcmp(ext,'.gz'))
%             volsgz=vols(vj,:);
			disp(['Gunzipping ' vols{vj,:}])
            vols(vj,:)=gunzip(vols{vj,:});
            do_delete_sub=1; %Mark sub for deletion
        end
    end
end

if 1 ~= exist('no_images_to_crop','var')
    no_images_to_crop=1; 
end

disp(['Cropping ', num2str(no_images_to_crop), ' images'])
% vols=strvcat(vols,other);
% vols=sub;

%Make sure to set to use elderly target (4) and crop first 4 images only
%     nii_setOrigin12x(vols, 4, true, no_images_to_crop,[-76 -112 -68; 76 76 104]);
%changed -z to -67 to force even number of slices in cropped brain.
%pout=nii_setOrigin12x(vols, 4, true, no_images_to_crop,[-76 -112 -67; 76 76 104]); %DEFAULT LIBERAL
pout=nii_setOrigin12x(vols, 4, true, no_images_to_crop,[-76 -112 -47; 76 76 104]);

if isempty(pout)
	disp(['EMPTY POUT is ' pout]) 
    for j=1:size(vols,1)
        [op,of,oe]=fileparts(char(vols{j,:}));
        pout{j,1}=fullfile(op,strcat('p',of,oe));
        movefile(char(vols{j,:}),char(pout{j})); 
    end
end

if do_delete_sub
	%Only delete the sub if it has been cropped. Otherwise rename it with a p
	%to be like the skull stripped images as it has been coregistered. 
	%if isempty(pout)
	%	disp(['EMPTY POUT is ' pout]) 
	%	[op,of,oe]=fileparts(char(sub));
     %   pout={fullfile(op,strcat('p',of,oe))};
	%	movefile(char(sub),char(pout));
	%else
	%	disp(['POUT is ' pout])
	%	delete(char(sub)); 
	%end
    for j=1:size(vols,1)
        disp(['Deleting sub ' vols(j,:)])
        delete(char(vols(j,:)));
    end
end

%zip the output
if do_gunzip
	for j=1:size(vols,1)
    	gzip(pout{j});
    	delete(pout{j});    
	end
end

