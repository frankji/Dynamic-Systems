function drawpolys(locs, cells, L)
% DRAWPOLYS(LOCS, CELLS, L)
%
% Draw cells with coordinates at LOCS, given cell network CELLS. Uses a
% periodic box of length L.
    Lset = [0,0;0,L(2);L(1),0;L];
    cols = [ ...
        255,255,255; ...
        255,255,255; ...
        0,0,0; ...
        77,175,74; ...
        255,255,51; ...
        180,180,180; ...
        55,126,184; ...
        228,26,28; ...
        152,78,163; ...
        255,127,0; ...
    ] / 255;
    
    for ix = 1:length(cells)
        ixs = cells{ix};
        colix = min(length(ixs), 9);
        color = cols(colix, :);
        cell_locs = locs(ixs, :);
        cell_drs = cell_locs;
        for vix=1:length(cell_drs)
            cell_drs(vix, :) = mod(cell_locs(vix, :) - cell_locs(1, :) + L/2, L) - L/2;
        end
        cell_mid = cell_locs(1, :) + mean(cell_drs, 1);
        cell_mid = mod(cell_mid + L/2, L) - L/2;

        cell_diffs = zeros(length(ixs), 2);
        for pix = 1:length(ixs)
            cell_diffs(pix, :) = mod(cell_locs(pix, :) - cell_mid + L/2, L) - L/2;
        end
        
        for i=1:length(Lset)
            new_mid = cell_mid + Lset(i, :);
            pts = [cell_diffs(:, 1) + new_mid(1), cell_diffs(:, 2) + new_mid(2)];
            
            if (min(pts(:, 1)) < L(1)) && (max(pts(:, 1)) > 0) && (min(pts(:, 2)) < L(2)) && (max(pts(:, 2)) > 0)
                fill(pts(:, 1), pts(:, 2), color);
                hold all;
            end
        end
    end
    
    daspect([1 1 1]);
    xlim([0, L(1)]);
    ylim([0, L(2)]);
end