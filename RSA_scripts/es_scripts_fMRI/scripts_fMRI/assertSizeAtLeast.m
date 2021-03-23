function assertSizeAtLeast(expected_size, cellarray, errormsg)
    if ~iscell(cellarray)
        error('assertSizeAtLeast can only be used for cell arrays');
    end
    
    if length(cellarray) < expected_size
        error(errormsg);
    end        
end