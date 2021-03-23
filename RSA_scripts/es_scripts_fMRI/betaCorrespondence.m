%  betaCorrespondence.m is a simple nullary function which should return three things:
%  	preBeta:	a string which is at the start of each file containing a beta image
%  	betas:	a struct indexed by (session, condition) containing a sting unique to each beta image
%  	postBeta:	a string which is at the end of each file containing a beta image, not containing the file .suffix
%
%	use "[[subjectName]]" as a placeholder for the subject's name as found in userOptions.subjectNames if necessary
%  
%  Cai Wingfield 1-2010

function betas = betaCorrespondence_H(userOptions)

%pre and post beta empty so as to fit with existing expectations of toolbox

mkdir([pwd,'/BetaList'])

preBeta = '';

if strcmp(userOptions.sourceImage,'average_con')
    %find path using toolbox format
    readPath = replaceWildcards(userOptions.betaSelector,'[[subjectName]]', userOptions.subjectNames);
    %select all con images
    all_con=spm_select('FPList',readPath{1},'^con.*\.nii$');
    for c=1:size(all_con,1)
        %read in SPM files
        temp=spm_vol(all_con(c,:));
       
        %create list of the averaged cons
        list_average{c,1}=temp.descrip;
        beta(1,c).identifier=sprintf('con_%04d',c);
    end
    clear readPath c
        
elseif strcmp(userOptions.sourceImage,'average_spmt')
    %find path using toolbox format
    readPath = replaceWildcards(userOptions.betaSelector,'[[subjectName]]', userOptions.subjectNames);
    %select all con images
    all_con=spm_select('FPList',readPath{1},'^spmT.*\.nii$');
    c=0;
    for i=1:size(all_con,1)
        %read in SPM files
        temp=spm_vol(all_con(i,:));
       
        %create list of the averaged cons
        %if regexp(temp.descrip,'Sn\(A\)')>0
        c=c+1;
        list_average{c,1}=temp.descrip;
        %end
    end
    clear i c temp
    
    for i=1:size(userOptions.conditionLabels,2)
        %identify associated condition label files
        ind=find(double(~cellfun('isempty',strfind(list_average,userOptions.conditionLabels{i}))));         
        
        %way to find con number directly from file data
        dash=strfind(list_average{ind},'-');
        eND=strfind(list_average{ind},':');    
        
        %remove space
        con_number=list_average{ind}(dash+11:eND-1);        
        
        if str2num(con_number) < 10
        beta(1,i).identifier=['spmT_000',con_number];
        elseif str2num(con_number) < 100 & str2num(con_number) > 9 
        beta(1,i).identifier=['spmT_00',con_number];
        else str2num(con_number) > 99
        beta(1,i).identifier=['spmT_0',con_number];
        end
        clear con_number dash eND ind
    end

    clear i readPath
    
elseif strcmp(userOptions.sourceImage,'single_con')
    %find path using toolbox format
    readPath = replaceWildcards(userOptions.betaSelector,'[[subjectName]]', userOptions.subjectNames);
    %select all con images
    all_con=spm_select('FPList',readPath{1},'^con.*\.nii$');
    c=0;
    for i=1:size(all_con,1)
        %read in SPM files
        temp=spm_vol(all_con(i,:));
       
        %create list of the averaged cons
        if regexp(temp.descrip,'Sn\([1-4]\)')>0
        c=c+1;
        list_average{c,1}=temp.descrip;
        end
    end
    clear i c temp

    for i=1:size(userOptions.conditionLabels,2)
        %identify associated condition label files
        ind=find(double(~cellfun('isempty',strfind(list_average,userOptions.conditionLabels{i}))));         
        
        %way to find con number directly from file data
        for dd=1:size(ind,1)
        dash(dd)=strfind(list_average{ind(dd)},'-');
        eND(dd)=strfind(list_average{ind(dd)},':');    
        con_number{dd}=list_average{ind(dd)}(dash(dd)+2:eND(dd)-1); 
        start_num_run{dd}=str2num(list_average{ind(dd)}(regexp(list_average{ind(dd)},'\([1-4]\)')+1));
         
        if str2num(con_number{dd}) < 10
        beta(start_num_run{dd},i).identifier=['con_000',con_number{dd}];
        elseif str2num(con_number{dd}) > 9 & str2num(con_number{dd}) < 100
        beta(start_num_run{dd},i).identifier=['con_00',con_number{dd}];
        else str2num(con_number{dd}) > 99
        beta(start_num_run{dd},i).identifier=['con_0',con_number{dd}];
        end
        end
        %remove space
              
        clear con_number dash eND ind start_num_run dd 
    end

    clear i readPath
elseif strcmp(userOptions.sourceImage,'single_spmt')
     %find path using toolbox format
    readPath = replaceWildcards(userOptions.betaSelector,'[[subjectName]]', userOptions.subjectNames);
    %select all con images
    all_con=spm_select('FPList',readPath{1},'^spmT.*\.nii$');
    c=0;
    for i=1:size(all_con,1)
        %read in SPM files
        temp=spm_vol(all_con(i,:));
       
        %create list of the averaged cons
        if regexp(temp.descrip,'Sn\([1-4]\)')>0
        c=c+1;
        list_average{c,1}=temp.descrip;
        end
    end
    clear i c temp
    
    for i=1:size(userOptions.conditionLabels,2)
        %identify associated condition label files
        ind=find(double(~cellfun('isempty',strfind(list_average,userOptions.conditionLabels{i}))));         
        
        %way to find con number directly from file data
        for dd=1:size(ind,1)
        dash(dd)=strfind(list_average{ind(dd)},'-');
        eND(dd)=strfind(list_average{ind(dd)},':');    
        con_number{dd}=list_average{ind(dd)}(dash(dd)+11:eND(dd)-1); 
        start_num_run{dd}=str2num(list_average{ind(dd)}(regexp(list_average{ind(dd)},'\([1-4]\)')+1))
         
        if str2num(con_number{dd}) < 10
        beta(start_num_run{dd},i).identifier=['spmT_000',con_number{dd}];
        elseif str2num(con_number{dd}) > 9 & str2num(con_number{dd}) < 100
        beta(start_num_run{dd},i).identifier=['spmT_00',con_number{dd}];
        else str2num(con_number{dd}) > 99
        beta(start_num_run{dd},i).identifier=['spmT_0',con_number{dd}];
        end
        end
        %remove space 
       
        clear con_number dash eND ind start_num_run dd
    end
    clear i readPath
else
    display('==Specifiy the correct name for input images==');
    return
end

postBeta = '.nii';

for session = 1:size(beta,1)
	for condition = 1:size(beta,2)
		betas(session,condition).identifier = [preBeta beta(session,condition).identifier postBeta];
	end%for
end%for
clear i c condition session 

pathSave=[cell2mat([pwd,'/BetaList/beta_',userOptions.subjectNames,'.mat'])]
save (pathSave,'betas', 'all_con', 'list_average','userOptions')
clear list_average all_con beta