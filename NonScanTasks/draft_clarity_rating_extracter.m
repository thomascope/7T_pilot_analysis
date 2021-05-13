fid=fopen('P3E1_modded_2_ONE_RUN-P7C17-Post.txt','r');  % open file 
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