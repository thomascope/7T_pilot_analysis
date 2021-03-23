% function [msg_length. msg] = display_progress(cfg,cnt,n_decodings,start_time,msg_length,add_to_message)
%
% Function to display progress and the estimated time to go in The Decoding
% Toolbox (i.e. how far is the analysis?).
%
% Call with msg_length = '' on first call, then use returned msg_length as
% input to next call (example below). This removes the previous message.
% Function does not overwrite tdt warnings.
%
% By default, the progress is only shown for selected values. See below
% cfg.display_progress.display_values how to provide your own.
%
% IN
%   cfg: struct with field or empty. See OPTIONAL below.
%   cnt: current loop iteration
%   n_decodings: total number of loop iterations
%   start_time: start time of loop, e.g. start_time = now
%   msg_length: length of last message to delete. Set '' or 0 to not delete
%     a previous message. If a tdt warning (flag global warningv_active~=0) 
%     occured, this will not be overwritten, and flag is reset.
%
% OPTIONAL
%  created text
%   cfg.analysis: string with analysis, e.g. 'ROI', added in message
%   cfg.sn: double, number of current subject, added in message
%   cfg.display_progress.string: string at beginning of message,
%     potentially after subject.
%   add_to_message: extra string to be added to message before progress
%  when to show progress
%   cfg.display_progress.display_values: values for which a message
%     should be shown. If field exist and is empty, progress will always
%     be shown.
%
% OUT
%   msg_length: chars to be deleted in next run
%   message: the generated message, excluding line breaks
%
% EXAMPLE
% msg_length = 0; % init
% n_decodings = 100;
% start_time = now;
% cfg.analysis = 'your_analysis'; % in TDT e.g. searchlight
% for cnt = 1:n_decodings
%    [msg_length] = display_progress(cfg,cnt,n_decodings,start_time,msg_length);
% end

% KG 2020/03/25: introduced option to show fractions of 1 (for large designs)
% changed how msg_length is used. changed when start_time is reported.
% MH 2015/06/18: introduced prev_time to get more accurate estimation of
% time to go if other processes start while analysis is running or if
% analysis is paused on the fly

function [msg_length, message] = display_progress(cfg,cnt,n_decodings,start_time,msg_length,add_to_message)

global warningv_active % was a warning shown in between? (otherwise message will be truncated)
if isempty(warningv_active), warningv_active = 0; end % init if not been done
persistent prev_time   % how much time has elapsed since last call?
persistent prev_start_time % show starting time if this time is new

message = []; % init message
if ~exist('msg_length', 'var'), msg_length = 0; end

if isempty(prev_start_time) ||  prev_start_time ~= start_time
    prev_start_time = start_time;
    fprintf('\nStarting time: %s\n',datestr(start_time));
end

if isfield(cfg, 'display_progress') && isfield(cfg.display_progress, 'display_values')
    if isempty(cfg.display_progress.display_values)
        % add current value to show current step if display_values is empty
        display_values = cnt;
    else
        % use provided display values
        display_values = cfg.display_progress.display_values;
    end
    is1000 = false;
    below1 = false;
else
    if n_decodings > 50
        % display progress of these iterations - afterwards at each 1000 steps
        display_values = [1 2 5 10 15 20 30 40 50 100 200 300 400 500 n_decodings];
    else
        display_values = 1:n_decodings;
    end
    is1000 = mod(cnt,1000) == 0;
    below1 = (cnt > 0 && cnt < 1);
end

if below1 || is1000 || any(display_values == cnt)
    
    if isfield(cfg, 'sn')
        message = [message sprintf('Subject: %02d, ', cfg.sn)];
    end
%         message = sprintf('Subject: %02d %s: %g/%g', cfg.sn, cfg.analysis, cnt, n_decodings);
    if isfield(cfg, 'display_progress') && isfield(cfg.display_progress, 'string')
        message = [message cfg.display_progress.string ' '];
    end
    
    if isfield(cfg, 'analysis')
        message = [message cfg.analysis ', '];
    end
    
    if exist('add_to_message', 'var')
        message = [message add_to_message ', '];
    end
    
    % if text thus far: replace final ', ' by ':' 
    if ~isempty(message)
        
        message(end-1:end) = ': ';
    end
    
    % add standard text cnt/n_decodings
    message = [message sprintf('%g/%g', cnt, n_decodings)];
    
    if cnt == 1
        % dont estimate on first round: message shown before start of
        % computation, so estimate not usefule
        message = [message ', estimating time on next call'];
    else
        % add estimated time to go
        t0 = now;
        el_time = t0 - start_time;
        el_time_str = datestr(el_time, 'dd HH:MM:SS');
        if str2double(el_time_str(1:2)) == 0, el_time_str = el_time_str(4:end); end
        if ~is1000 || cnt == 1000 % if less than 1000 iterations, base estimation on all previous
            est_time =  n_decodings/max(cnt-1, 1) * el_time;
        else % otherwise base on most recent 1000 only
            est_time =  el_time + (n_decodings-cnt)/1000 * (t0 - prev_time);
        end
        prev_time = t0; % update
        est_time_left = est_time - el_time; % how long we think it will still take
        est_time_left_str = datestr(est_time_left, 'dd HH:MM:SS');
        if str2double(est_time_left_str(1:2)) == 0, est_time_left_str = est_time_left_str(4:end); end
        est_finish = start_time + est_time;
        est_finish_str = datestr(est_finish, 'yyyy/mm/dd HH:MM:SS');
        message = [message ', time to go: ' est_time_left_str ', time running: ' el_time_str ', finish: ' est_finish_str];
    end
    
    % print message and delete old message
    
    % get delete characters to delete previous message
    if ~isempty(msg_length) && ~warningv_active
        % msg_length includes potential line breaks etc.
        reverse_str = repmat('\b', 1, msg_length); % to delete old text, for fprintf
    else
        % dont delete previous message (warning or first call)
        reverse_str = [];
    end
    % message_str
    message_str = ['\n' message '\n\n']; % new text, for fprintf
    
    % delete old text and write new text. do in one call to avoid flicker
    fprintf([reverse_str message_str]);

    % get number of printed chars. can differ from length(message)+3 due to interpretation in s/fprintf    
    msg_length = length(sprintf(message_str));
    
    % reset warningv_active flag
    warningv_active = 0; % this not deactivate warningv but resets waringv_active. it is a flag that, if active, prevents that previous (warning) text is deleted here.
end