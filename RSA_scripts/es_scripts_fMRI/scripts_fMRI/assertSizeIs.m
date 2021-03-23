function assertSizeIs(expected_size, cellarray, errormsg)
    if ~iscell(cellarray)
        error('assertSizeIs can only be used for cell arrays');
    end
    
    if length(cellarray) ~= expected_size
        error(errormsg);
    end        
end