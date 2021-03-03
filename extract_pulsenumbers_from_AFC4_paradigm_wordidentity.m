function [startpulses,stimType,stimNumber,stimName,writtenNumber,rts] = extract_pulsenumbers_from_AFC4_paradigm_wordidentity(fileName,runI)

%A function for outputting the values required for MVPA. 
%Stilltodo: Divide mismatches by vowel transition. 
%Sort out response trials if we want to decode these.

load(fileName)

startpulses = [imputed_pulse_numbers_at_normal_trialstart, imputed_pulse_numbers_at_writtenonly_trialstart, imputed_pulse_numbers_at_response_trialstart];

if length(this_cue_types) ~= length(imputed_pulse_numbers_at_normal_trialstart)
    warning('The length of the cue and imputed pulsenumbers do not match - assuming truncated acquisition')
    this_cue_types = this_cue_types(1:length(imputed_pulse_numbers_at_normal_trialstart));
end

% there were 8 conditions - Match low, Match high, Mismatch low, Mismatch high, Neutral low, Neutral high, Writtenonly, Response
stimType = zeros(1,length(startpulses));

for i = 1:length(this_cue_types);
    if this_cue_types(i) == 1
        if this_vocoder_channels(i)  == 1
            stimType(i) = 1;
        elseif this_vocoder_channels(i)  == 2
            stimType(i) = 2;
        end
    elseif this_cue_types(i) == 2
        if this_vocoder_channels(i)  == 1
            stimType(i) = 3;
        elseif this_vocoder_channels(i)  == 2
            stimType(i) = 4;
        end
    elseif this_cue_types(i) == 3
        if this_vocoder_channels(i)  == 1
            stimType(i) = 5;
        elseif this_vocoder_channels(i)  == 2
            stimType(i) = 6;
        end
    end
end

stimType(i+1:i+length(imputed_pulse_numbers_at_writtenonly_trialstart)) = 7;
stimType(i+1+length(imputed_pulse_numbers_at_writtenonly_trialstart):end) = 8;

% there were 16 words 
stimNumber = [this_word(1:length(imputed_pulse_numbers_at_normal_trialstart)) written_trial_word(1:length(imputed_pulse_numbers_at_writtenonly_trialstart))]; %Easy for the normal and written only trials as recorded by the delivery script

stimNumber = this_word(1:length(imputed_pulse_numbers_at_normal_trialstart)); %Easy for the normal trials as recorded by the delivery script
writtenNumber = this_mismatch_cue; %First mismatch trials
writtenNumber(isnan(this_mismatch_cue))=this_word(isnan(this_mismatch_cue)); %Then match trials

%Now do written trials

writtenNumber = [writtenNumber written_trial_word];
stimNumber = [stimNumber zeros(size(written_trial_word))];

response_indexC = strfind(trialtype,'Response Trial');
response_index = find(not(cellfun('isempty', response_indexC)));
for i = 1:length(response_index)
    j = i+(length(imputed_pulse_numbers_at_response_trialstart)*(runI-1));
    if strcmp('Match',cue_types{response_order.this_cue_types(j)}) 
        writtenNumber(end+1) = response_order.this_word(j);
        stimNumber(end+1) = response_order.this_word(j);
    elseif strcmp('MisMatch',cue_types{response_order.this_cue_types(j)})
        writtenNumber(end+1) = response_order.this_mismatch_cue(j);
        stimNumber(end+1) = response_order.this_word(j);
    elseif strcmp('Neutral',cue_types{response_order.this_cue_types(j)})
        writtenNumber(end+1) = 0;
        stimNumber(end+1) = response_order.this_word(j);
    end
end

stimName = wordlist;
rts = all_rts(all_rts~=0);

end