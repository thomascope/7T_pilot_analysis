function model = distance_train(labels_train,data_train,cfg)

% This function averages samples within the same class. It then passes the
% data for calculations of distance estimate carried out in similarity_test.

% TODO: place calculations in separate function even though they are
% simple.

switch lower(cfg.decoding.method)
    
    case 'classification'
        % Get unique labels
        u_labels = uniqueq(labels_train); % sorts labels!
        n_labels = size(u_labels,1);
        
        for i_label = n_labels:-1:1
            label_ind = labels_train==u_labels(i_label);
            m(:,i_label) = (1/sum(label_ind)) *sum(data_train(label_ind,:),1);
        end
        
        model.vectors_train = m;
        model.labels_train =  u_labels;
        
    case 'classification_kernel'
        error('cfg.decoding.method = ''classification_kernel''. Distance calculations currently don''t work with passed kernels.')
        
    otherwise
        error(...
           ['The "similarity" decoding software used for calculating distances (cfg.decoding.software = ''similarity'') ', ...
           'only takes cfg.decoding.method = ''classification'', to avoid confusions. ', ...
           'The currently set method is ''cfg.decoding.method = %s'' ', ...
           'for cfg.decoding.software = %s. ', ...
           'Please change.'],...
            cfg.decoding.method, cfg.decoding.software)
end