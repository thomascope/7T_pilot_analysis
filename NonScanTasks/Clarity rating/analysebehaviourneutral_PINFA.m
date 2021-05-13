function dat = analysebehaviourneutral_PINFA

DIRS = getdirs_psych;
fpath = DIRS.data;

if ~exist(DIRS.dataFigs,'dir')
    mkdir(DIRS.dataFigs)
end

multi = 1;
Mfilename = '_All'; %for concatenated runs the name will be: ['SubjInitials' 'Mfilename']

Timepoints = {'Pre','Post'};

Patfiles = ls([fpath 'Patients\P3E1*' Timepoints{1} '.txt']);
Confiles = ls([fpath 'Controls\P3E1*' Timepoints{1} '.txt']);
startnumber = strfind(Patfiles(1,:),'P7P');
Patinits = Patfiles(:,startnumber:startnumber+4); %Assumes all 5 letter ID codes and 2 digit ages.
for i = 1:size(Patinits,1)
    if Patinits(i,4) == '_'
        Patinits(i,4) = [' '];
    end
end
disp('The following patients found: ');
Patinits
startnumber = strfind(Confiles(1,:),'P7C');
Coninits = unique(Confiles(:,startnumber:startnumber+4),'rows');
for i = 1:size(Coninits,1)
    if Coninits(i,4) == '_'
        Coninits(i,4) = [' '];
    end
end
disp('The following controls found: ');
Coninits
for I = 1:size(Patfiles,1)
    PatSes{I} = ['Patients\' Patfiles(I,:)];
end
for I = 1:size(Confiles,1)
    ConSes{I} = ['Controls\' Confiles(I,:)];
end

man = input('Process these? [1] or get files manually [0] : ');
if ~man,
    S = 1; Ses = [];
    while S
        F = uigetfile('*.mat','Open Each File - Cancel to Quit');
        if F == 0, break; end
        Ses{S}  = F;
        S = S+1;
    end
end

for this_time = 1:length(Timepoints)
    clear patdatA condatA
    for S = 1:size(PatSes,2),
        if this_time ~= 1
            PatSes{S} = strrep(PatSes{S},Timepoints{this_time-1},Timepoints{this_time});
        end
        disp(['Working on ' Timepoints{this_time} ' patient file: ' fpath PatSes{S}]);
        fname = [fpath PatSes{S}];
        if ~exist(fname,'file')
            continue
        end
        
        fid=fopen(fname,'r','n','Unicode');  % open file
        all_resps = [];
        all_congruencies = {};
        all_bands = [];
        while 1
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            if strfind(tline,'CollectResponse.RESP: ')
                this_resp = strsplit(tline, 'CollectResponse.RESP: ');
                all_resps = [all_resps, str2num(this_resp{end})];
            elseif strfind(tline,'Congruency: ')
                this_congruency = strsplit(tline, 'Congruency: ');
                all_congruencies = [all_congruencies, this_congruency{end}];
            elseif strfind(tline,'Bands: ')
                this_bands = strsplit(tline, 'Bands: ');
                all_bands = [all_bands, str2num(this_bands{end})];
            end
        end
        fclose(fid);
        
        all_bands = all_bands(1:length(all_resps)); %If run trunkated, then response is collected last, get surviving data.
        all_congruencies = all_congruencies(1:length(all_resps));
        
        clear patdat
        patdat.filename = fname;
        
        patdat.trialarray = zeros(length(all_resps),3);
        patdat.trialarray(:,3) = all_resps;
        patdat.trialarray(:,2) = all_bands;
        patdat.trialarray(strcmp(all_congruencies,'Match'),1) = 1;
        patdat.trialarray(strcmp(all_congruencies,'Mismatch'),1) = 2;
        patdat.trialarray(strcmp(all_congruencies,'Neutral'),1) = 3;
        
        if patdat.trialarray(1,3) > 4
            patdat.trialarray(:,3) = patdat.trialarray(:,3)-4;
        end
        
        patdat.match4=[];
        patdat.match8 = [];
        patdat.match16 = [];
        patdat.mismatch4 = [];
        patdat.mismatch8 = [];
        patdat.mismatch16 = [];
        patdat.neutral4 = [];
        patdat.neutral8 = [];
        patdat.neutral16 = [];
        for i = 1:size(patdat.trialarray,1)
            if patdat.trialarray(i,1)==1 && patdat.trialarray(i,2)==4
                patdat.match4=[patdat.match4, patdat.trialarray(i,3)];
            elseif patdat.trialarray(i,1)==2 && patdat.trialarray(i,2)==4
                patdat.mismatch4=[patdat.mismatch4, patdat.trialarray(i,3)];
            elseif patdat.trialarray(i,1)==3 && patdat.trialarray(i,2)==4
                patdat.neutral4=[patdat.neutral4, patdat.trialarray(i,3)];
            elseif patdat.trialarray(i,1)==1 && patdat.trialarray(i,2)==8
                patdat.match8=[patdat.match8, patdat.trialarray(i,3)];
            elseif patdat.trialarray(i,1)==2 && patdat.trialarray(i,2)==8
                patdat.mismatch8=[patdat.mismatch8, patdat.trialarray(i,3)];
            elseif patdat.trialarray(i,1)==3 && patdat.trialarray(i,2)==8
                patdat.neutral8=[patdat.neutral8, patdat.trialarray(i,3)];
            elseif patdat.trialarray(i,1)==1 && patdat.trialarray(i,2)==16
                patdat.match16=[patdat.match16, patdat.trialarray(i,3)];
            elseif patdat.trialarray(i,1)==2 && patdat.trialarray(i,2)==16
                patdat.mismatch16=[patdat.mismatch16, patdat.trialarray(i,3)];
            elseif patdat.trialarray(i,1)==3 && patdat.trialarray(i,2)==16
                patdat.neutral16=[patdat.neutral16, patdat.trialarray(i,3)];
            end
        end
        
        patdat.meansarray = [mean(patdat.match4),mean(patdat.neutral4),mean(patdat.mismatch4);mean(patdat.match8),mean(patdat.neutral8),mean(patdat.mismatch8);mean(patdat.match16),mean(patdat.neutral16),mean(patdat.mismatch16)];
        stdsarray = [std(patdat.match4),std(patdat.mismatch4),std(patdat.neutral4);std(patdat.match8),std(patdat.mismatch8),std(patdat.neutral8);std(patdat.match16),std(patdat.mismatch16),std(patdat.neutral16)];
        stesarray = stdsarray./sqrt(size(patdat.match4,2));   %Assumes equal number of trials in every condition
        figure
        set(gcf,'position',[100,100,1200,800])
        barweb(patdat.meansarray,stesarray,[],{'4 channels';'8 channels';'16 channels'},['Clarity Ratings by Prime Type and Vocoder Channels for Patient ' Patinits(S,:)],[],'Mean Clarity Rating',[],[],{'Match','Neutral','Mismatch'}) ;
        legend('Match','Neutral','Mismatch','location','NorthWest');
        set(gca,'ylim',[1,4])
        
        patdat.stesarray = stesarray;
        %save individual "patdat" structures and figures
        save([fpath 'beh_' PatSes{S}(11:end-4) '.mat'],'patdat');
        saveas(gcf,[DIRS.dataFigs 'Neutral Rating Patient ' Patinits(S,:) '_' Timepoints{this_time} '.jpeg'])
        
        %create concatenated datafile
        if S == 1,
            patdatA{1} = patdat;
            patdatA{1}.filename = [];
            patdatA{1}.filename = patdat.filename;
        else
            patdatA{end+1}.filename = patdat.filename;
            patdatA{end} = patdat;
        end
        
    end
    for S = 1:size(ConSes,2),
        if this_time ~= 1
            ConSes{S} = strrep(ConSes{S},Timepoints{this_time-1},Timepoints{this_time});
        end
        disp(['Working on ' Timepoints{this_time} ' control file: ' fpath ConSes{S}]);
        fname = [fpath ConSes{S}];
        
        fid=fopen(fname,'r','n','Unicode');  % open file
        all_resps = [];
        all_congruencies = {};
        all_bands = [];
        while 1
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            if strfind(tline,'CollectResponse.RESP: ')
                this_resp = strsplit(tline, 'CollectResponse.RESP: ');
                all_resps = [all_resps, str2num(this_resp{end})];
            elseif strfind(tline,'Congruency: ')
                this_congruency = strsplit(tline, 'Congruency: ');
                all_congruencies = [all_congruencies, this_congruency{end}];
            elseif strfind(tline,'Bands: ')
                this_bands = strsplit(tline, 'Bands: ');
                all_bands = [all_bands, str2num(this_bands{end})];
            end
        end
        fclose(fid);
        
        all_bands = all_bands(1:length(all_resps)); %If run trunkated, then response is collected last, get surviving data.
        all_congruencies = all_congruencies(1:length(all_resps));
        
        clear patdat
        condat.filename = fname;
        
        condat.trialarray = zeros(length(all_resps),3);
        condat.trialarray(:,3) = all_resps;
        condat.trialarray(:,2) = all_bands;
        condat.trialarray(strcmp(all_congruencies,'Match'),1) = 1;
        condat.trialarray(strcmp(all_congruencies,'Mismatch'),1) = 2;
        condat.trialarray(strcmp(all_congruencies,'Neutral'),1) = 3;
        
        dgz = load(deblank(fname));
        condat.filename = fname;
        
        if condat.trialarray(1,3) > 4
            condat.trialarray(:,3) = condat.trialarray(:,3)-4;
        end
        
        condat.match4=[];
        condat.match8 = [];
        condat.match16 = [];
        condat.mismatch4 = [];
        condat.mismatch8 = [];
        condat.mismatch16 = [];
        condat.neutral4 = [];
        condat.neutral8 = [];
        condat.neutral16 = [];
        for i = 1:size(condat.trialarray,1)
            if condat.trialarray(i,1)==1 && condat.trialarray(i,2)==4
                condat.match4=[condat.match4, condat.trialarray(i,3)];
            elseif condat.trialarray(i,1)==2 && condat.trialarray(i,2)==4
                condat.mismatch4=[condat.mismatch4, condat.trialarray(i,3)];
            elseif condat.trialarray(i,1)==3 && condat.trialarray(i,2)==4
                condat.neutral4=[condat.neutral4, condat.trialarray(i,3)];
            elseif condat.trialarray(i,1)==1 && condat.trialarray(i,2)==8
                condat.match8=[condat.match8, condat.trialarray(i,3)];
            elseif condat.trialarray(i,1)==2 && condat.trialarray(i,2)==8
                condat.mismatch8=[condat.mismatch8, condat.trialarray(i,3)];
            elseif condat.trialarray(i,1)==3 && condat.trialarray(i,2)==8
                condat.neutral8=[condat.neutral8, condat.trialarray(i,3)];
            elseif condat.trialarray(i,1)==1 && condat.trialarray(i,2)==16
                condat.match16=[condat.match16, condat.trialarray(i,3)];
            elseif condat.trialarray(i,1)==2 && condat.trialarray(i,2)==16
                condat.mismatch16=[condat.mismatch16, condat.trialarray(i,3)];
            elseif condat.trialarray(i,1)==3 && condat.trialarray(i,2)==16
                condat.neutral16=[condat.neutral16, condat.trialarray(i,3)];
            end
        end
        
        condat.meansarray = [mean(condat.match4),mean(condat.neutral4),mean(condat.mismatch4);mean(condat.match8),mean(condat.neutral8),mean(condat.mismatch8);mean(condat.match16),mean(condat.neutral16),mean(condat.mismatch16)];
        stdsarray = [std(condat.match4),std(condat.mismatch4),std(condat.neutral4);std(condat.match8),std(condat.mismatch8),std(condat.neutral8);std(condat.match16),std(condat.mismatch16),std(condat.neutral16)];
        stesarray = stdsarray./sqrt(size(condat.match4,2));   %Assumes equal number of trials in every condition
        figure
        set(gcf,'position',[100,100,1200,800])
        barweb(condat.meansarray,stesarray,[],{'4 channels';'8 channels';'16 channels'},['Clarity Ratings by Prime Type and Vocoder Channels for Control ' Coninits(S,:)],[],'Mean Clarity Rating',[],[],{'Match','Neutral','Mismatch'}) ;
        legend('Match','Neutral','Mismatch','location','NorthWest');
        set(gca,'ylim',[1,4])
        
        condat.stesarray = stesarray;
        %save individual "Condat" structures
        save([fpath 'beh_' ConSes{S}(14:end-4) '.mat'],'condat');
        saveas(gcf,[DIRS.dataFigs 'Neutral Rating Control ' Coninits(S,:) '_' Timepoints{this_time} '.jpeg'])
        
        
        
        if S == 1,
            condatA{1} = condat;
            condatA{1}.filename = [];
            condatA{1}.filename = condat.filename;
        else
            condatA{end+1}.filename = condat.filename;
            condatA{end} = condat;
        end
        
    end
    
    if multi,
        
        patdat.trialarray = [];
        for i = 1:length(patdatA)
            patdat.trialarray = [patdat.trialarray; patdatA{i}.trialarray];
        end
        
        if patdat.trialarray(1,3) > 4
            patdat.trialarray(:,3) = patdat.trialarray(:,3)-4;
        end
        
        patdat.match4=[];
        patdat.match8 = [];
        patdat.match16 = [];
        patdat.mismatch4 = [];
        patdat.mismatch8 = [];
        patdat.mismatch16 = [];
        patdat.neutral4 = [];
        patdat.neutral8 = [];
        patdat.neutral16 = [];
        for i = 1:size(patdat.trialarray,1)
            if patdat.trialarray(i,1)==1 && patdat.trialarray(i,2)==4
                patdat.match4=[patdat.match4, patdat.trialarray(i,3)];
            elseif patdat.trialarray(i,1)==2 && patdat.trialarray(i,2)==4
                patdat.mismatch4=[patdat.mismatch4, patdat.trialarray(i,3)];
            elseif patdat.trialarray(i,1)==3 && patdat.trialarray(i,2)==4
                patdat.neutral4=[patdat.neutral4, patdat.trialarray(i,3)];
            elseif patdat.trialarray(i,1)==1 && patdat.trialarray(i,2)==8
                patdat.match8=[patdat.match8, patdat.trialarray(i,3)];
            elseif patdat.trialarray(i,1)==2 && patdat.trialarray(i,2)==8
                patdat.mismatch8=[patdat.mismatch8, patdat.trialarray(i,3)];
            elseif patdat.trialarray(i,1)==3 && patdat.trialarray(i,2)==8
                patdat.neutral8=[patdat.neutral8, patdat.trialarray(i,3)];
            elseif patdat.trialarray(i,1)==1 && patdat.trialarray(i,2)==16
                patdat.match16=[patdat.match16, patdat.trialarray(i,3)];
            elseif patdat.trialarray(i,1)==2 && patdat.trialarray(i,2)==16
                patdat.mismatch16=[patdat.mismatch16, patdat.trialarray(i,3)];
            elseif patdat.trialarray(i,1)==3 && patdat.trialarray(i,2)==16
                patdat.neutral16=[patdat.neutral16, patdat.trialarray(i,3)];
            end
        end
        
        meansarray = [mean(patdat.match4),mean(patdat.neutral4),mean(patdat.mismatch4);mean(patdat.match8),mean(patdat.neutral8),mean(patdat.mismatch8);mean(patdat.match16),mean(patdat.neutral16),mean(patdat.mismatch16)];
        %stdsarray = [std(patdat.match4),std(patdat.mismatch4),std(patdat.neutral4);std(patdat.match8),std(patdat.mismatch8),std(patdat.neutral8);std(patdat.match16),std(patdat.mismatch16),std(patdat.neutral16)];
        patgroupmeans = zeros(size(patdatA{1}.meansarray,1),size(patdatA{1}.meansarray,2),length(patdatA));
        for i = 1:length(patdatA)
            patgroupmeans(:,:,i) = [patdatA{i}.meansarray];
        end
        stdsarray = zeros(size(patdatA{1}.meansarray));
        for i = 1:size(patdatA{i}.meansarray,1)
            for j = 1:size(patdatA{i}.meansarray,2)
                stdsarray(i,j) = std(patgroupmeans(i,j,:));
            end
        end
        stesarray = stdsarray./sqrt(length(patdatA));
        
        
        figure
        set(gcf,'position',[100,100,1200,800])
        barweb(meansarray,stesarray,[],{'4 channels';'8 channels';'16 channels'},['Clarity Ratings by Prime Type and Vocoder Channels for All Patients '],[],'Mean Clarity Rating',[],[],{'Match','Neutral','Mismatch'}) ;
        legend('Match','Neutral','Mismatch','location','NorthWest');
        set(gca,'ylim',[1,4])
        
        patdat.allmeansarray = meansarray;
        patdat.allstdsarray = stdsarray;
        patdat.allstesarray = stesarray;
        
        %save individual "patdat" structures
        save([fpath 'beh_all'],'patdat');
        saveas(gcf,[DIRS.dataFigs 'Neutral Rating All Patients ' Timepoints{this_time} ' .jpeg'])
        
        condat.trialarray = [];
        for i = 1:length(condatA)
            condat.trialarray = [condat.trialarray; condatA{i}.trialarray];
        end
        
        
        if condat.trialarray(1,3) > 4
            condat.trialarray(:,3) = condat.trialarray(:,3)-4;
        end
        
        condat.match4=[];
        condat.match8 = [];
        condat.match16 = [];
        condat.mismatch4 = [];
        condat.mismatch8 = [];
        condat.mismatch16 = [];
        condat.neutral4 = [];
        condat.neutral8 = [];
        condat.neutral16 = [];
        for i = 1:size(condat.trialarray,1)
            if condat.trialarray(i,1)==1 && condat.trialarray(i,2)==4
                condat.match4=[condat.match4, condat.trialarray(i,3)];
            elseif condat.trialarray(i,1)==2 && condat.trialarray(i,2)==4
                condat.mismatch4=[condat.mismatch4, condat.trialarray(i,3)];
            elseif condat.trialarray(i,1)==3 && condat.trialarray(i,2)==4
                condat.neutral4=[condat.neutral4, condat.trialarray(i,3)];
            elseif condat.trialarray(i,1)==1 && condat.trialarray(i,2)==8
                condat.match8=[condat.match8, condat.trialarray(i,3)];
            elseif condat.trialarray(i,1)==2 && condat.trialarray(i,2)==8
                condat.mismatch8=[condat.mismatch8, condat.trialarray(i,3)];
            elseif condat.trialarray(i,1)==3 && condat.trialarray(i,2)==8
                condat.neutral8=[condat.neutral8, condat.trialarray(i,3)];
            elseif condat.trialarray(i,1)==1 && condat.trialarray(i,2)==16
                condat.match16=[condat.match16, condat.trialarray(i,3)];
            elseif condat.trialarray(i,1)==2 && condat.trialarray(i,2)==16
                condat.mismatch16=[condat.mismatch16, condat.trialarray(i,3)];
            elseif condat.trialarray(i,1)==3 && condat.trialarray(i,2)==16
                condat.neutral16=[condat.neutral16, condat.trialarray(i,3)];
            end
        end
        
        meansarray = [mean(condat.match4),mean(condat.neutral4),mean(condat.mismatch4);mean(condat.match8),mean(condat.neutral8),mean(condat.mismatch8);mean(condat.match16),mean(condat.neutral16),mean(condat.mismatch16)];
        %stdsarray = [std(condat.match4),std(condat.mismatch4),std(condat.neutral4);std(condat.match8),std(condat.mismatch8),std(condat.neutral8);std(condat.match16),std(condat.mismatch16),std(condat.neutral16)];
        congroupmeans = zeros(size(patdatA{1}.meansarray,1),size(patdatA{1}.meansarray,2),length(patdatA));
        for i = 1:length(patdatA)
            congroupmeans(:,:,i) = [patdatA{i}.meansarray];
        end
        stdsarray = zeros(size(patdatA{1}.meansarray));
        for i = 1:size(patdatA{i}.meansarray,1)
            for j = 1:size(patdatA{i}.meansarray,2)
                stdsarray(i,j) = std(congroupmeans(i,j,:));
            end
        end
        stesarray = stdsarray./sqrt(length(condatA));   %Assumes equal number of trials in every condition
        figure
        set(gcf,'position',[100,100,1200,800])
        barweb(meansarray,stesarray,[],{'4 channels';'8 channels';'16 channels'},['Clarity Ratings by Prime Type and Vocoder Channels for All Controls '],[],'Mean Clarity Rating',[],[],{'Match','Neutral','Mismatch'}) ;
        legend('Match','Neutral','Mismatch','location','NorthWest');
        set(gca,'ylim',[1,4])
        
        condat.allmeansarray = meansarray;
        condat.allstdsarray = stdsarray;
        condat.allstesarray = stesarray;
        
        %save individual "condat" structures
        save([fpath 'beh_all'],'condat');
        saveas(gcf,[DIRS.dataFigs 'Neutral Rating All Controls ' Timepoints{this_time} ' .jpeg'])
        
        %compare group means
        allmeansarray = zeros(size(patdat.allmeansarray,1),size(patdat.allmeansarray,2)+size(condat.allmeansarray,2));
        allstesarray = zeros(size(patdat.allstesarray,1),size(patdat.allstesarray,2)+size(condat.allstesarray,2));
        for i = 1:size(allmeansarray,1)
            allmeansarray(i,:) = [condat.allmeansarray(i,1), patdat.allmeansarray(i,1), condat.allmeansarray(i,2), patdat.allmeansarray(i,2), condat.allmeansarray(i,3), patdat.allmeansarray(i,3)];
            allstesarray(i,:) = [condat.allstesarray(i,1), patdat.allstesarray(i,1), condat.allstesarray(i,2), patdat.allstesarray(i,2), condat.allstesarray(i,3), patdat.allstesarray(i,3)];
        end
        
        figure
        set(gcf,'position',[100,100,1200,800])
        barweb(allmeansarray,allstesarray,[],{'4 channels';'8 channels';'16 channels'},['Clarity Ratings by Prime Type and Vocoder Channels for All Subjects '],[],'Mean Clarity Rating',[],[],{'ControlMatch','PatientMatch','ControlNeutral','PatientNeutral','ControlMismatch','PatientMismatch'}) ;
        set(gca,'ylim',[1,4])
        legend('ControlMatch','PatientMatch','ControlNeutral','PatientNeutral','ControlMismatch','PatientMismatch','location','NorthWest');
        saveas(gcf,[DIRS.dataFigs 'Neutral Rating All subjects ' Timepoints{this_time} ' .jpeg'])
        
        figure
        set(gcf,'position',[100,100,1200,800])
        lineplot = tight_subplot(1,2,[0 0],[.1 .1],[.1 .1]);
        axes(lineplot(1));
        errorbar(allmeansarray(:,1),allstesarray(:,1),'r','linewidth',3);
        hold on
        set(gca,'ylim',[1,4],'LineWidth', 2, 'Xtick', [1 2 3], 'XTickLabel',[4,8,16],'Fontsize',[14],'FontName','Tahoma')
        errorbar(allmeansarray(:,3),allstesarray(:,3),'k','linewidth',3);
        errorbar(allmeansarray(:,5),allstesarray(:,5),'b','linewidth',3);
        firstlegend = legend('Match','Neutral','Mismatch','location','NorthWest');
        set(firstlegend,'FontSize',18);
        title('Controls','Color','k','fontsize',20)
        ylabel('Clarity Rating')
        xlabel('Vocode Channels')
        
        axes(lineplot(2));
        set(gca,'ylim',[1,4])
        errorbar(allmeansarray(:,2),allstesarray(:,2),'r','linewidth',3);
        hold on
        errorbar(allmeansarray(:,4),allstesarray(:,4),'k','linewidth',3);
        errorbar(allmeansarray(:,6),allstesarray(:,6),'b','linewidth',3);
        set(gca,'ylim',[1,4],'LineWidth', 2, 'Xtick', [1 2 3], 'XTickLabel',[4,8,16],'Fontsize',[14],'FontName','Tahoma','YAxisLocation','right')
        secondlegend = legend('Match','Neutral','Mismatch','location','NorthWest');
        set(secondlegend,'FontSize',18);
        title('Patients','Color','k','fontsize',20)
        ylabel('Clarity Rating')
        xlabel('Vocode Channels')
        img = getframe(gcf);
        imwrite(img.cdata, [DIRS.dataFigs 'Neutral Rating Linegraphs ' Timepoints{this_time} ' .jpeg']);
    end
    
    %save multiple runs by subject's initials
    disp(['Saving concatenated patient runs as: beh_patient' Mfilename]);
    save([fpath 'beh_patient' Mfilename '_' Timepoints{this_time}],'patdatA');
    disp(['Saving concatenated control runs as: beh_control' Mfilename]);
    save([fpath 'beh_control' Mfilename '_' Timepoints{this_time}],'condatA');
end