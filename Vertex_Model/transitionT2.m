function [newcells, newlocs] = transitionT2(cells, locs, L, cellix)
%TRANSITIONT2 Remove a cell from the cell set, as well as its nodes.
    cell = cells{cellix};
    newcells = [cells(1:cellix-1); cells(cellix+1:end)];
    cellmid = locs(cell(1), :);
    for ix = cell(2:end)
        dr = mod(locs(ix, :) - locs(cell(1), :) + L/2, L) - L/2;
        cellmid = cellmid + dr / length(cell);
    end
    
    if length(cell) == 3
        celllocix = min(cell);
        remove_ix = cell(cell ~= celllocix);
    elseif length(cell) == 2
        remove_ix = cell;
    else
        error(['Not Implemented: I don''t know how to remove a ' ...
            'cell with %d nodes, only 2 or 3 nodes.'], ...
            length(cell));
    end
    
    ixset = 1:length(locs);
    new_ix_set = ixset - cumsum(ismember(ixset, remove_ix));
    newix = setdiff(ixset, remove_ix);
    newlocs = locs(newix, :);
    
    if length(cell) == 3
        newlocs(celllocix, :) = cellmid;
    end
    
    numcells = zeros(length(newlocs), 1);
    for cix = 1:length(newcells)
        other_cell = newcells{cix};
        if length(cell) == 3
            other_cell(ismember(other_cell, remove_ix)) = celllocix;
            other_cell = unique(other_cell, 'stable');
        elseif length(cell) == 2
            other_cell = other_cell(~ismember(other_cell, remove_ix));
        end
        other_cell = new_ix_set(other_cell);
        newcells{cix} = other_cell;
        numcells(other_cell) = numcells(other_cell) + 1;
    end
    
    assert(all(numcells == 3));
end

