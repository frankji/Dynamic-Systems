function [newcells, newlocs, T1s, T2s, tries] = transition(cells, locs, L, coeffs, A0, lowdist)
% TRANSITION Attempt T1 and T2 transitions, and return new cells and new
% locs when any transition occurs. It is suggested to attempt energy
% minimization between transition attempts, and keep applying the
% TRANSITION in a WHILE loop until no transitions occur.

% CELLS is a cell array indicating which vertices belong to which cell.
% LOCS is an Nx2 array of the locations of each vertex, with the row in
%     LOCS matching the vertices specified in CELLS.
% L is a 1x2 array of the x- and y- lengths of the periodic boundary
%     conditions.
% COEFFS should be a 1x3 specifying [K, LAMBDA, GAMMA].
% A0 is the preferred area of the cells.
% LOWDIST specifies how small a side should be before a T1 transition is
%     attempted. 
%
% T2 transitions will be attempted for any 2- or 3-sided polygon.
% Transitions are only accepted if they lower the total energy.
%
% Output
% NEWCELLS and NEWLOCS specify the new network and particle positions
%     yielded by the transitions.
% T1S and T2S are the number of T1s and T2s that took place. If both are
%     zero, then you can expect NEWCELLS and NEWLOCS to equal the input
%     LOCS and CELLS.
% TRIES is the number of transitions that were attempted. Only those which
%     lowered the total energy will be accepted.

    %[EA, EL, EP] = cellenergies(A0, cells, locs, L);
    %E0 = [EA, EL, EP] * coeffs';
    [E0,~] = scaleenergy(locs,...
         cells,L,1,coeffs(1),coeffs(2),coeffs(3),A0);
    T1s = 0;
    T2s = 0;
    tries = 0;
    newcells = cells;
    newlocs = locs;
    
    cix = 1;
    while cix <= length(newcells)
        cell = newcells{cix};
        P = 0;
        
        N = length(cell);
        if N <= 3
            tries = tries + 1;
            [newcellsA, newlocsA] = transitionT2(newcells, newlocs, L, cix);
            %[EA, EL, EP] = cellenergies(A0, newcellsA, newlocsA, L);
            %E1 = [EA, EL, EP] * coeffs';
            [E1,~] = scaleenergy(newlocsA,...
                newcellsA,L,1,coeffs(1),coeffs(2),coeffs(3),A0);
            if E1 < E0
                T2s = T2s + 1;
                newcells = newcellsA;
                newlocs = newlocsA;
                return
            end
        end
        
        for i=1:length(cell)
            ix1 = cell(i);
            ix2 = cell(mod(i, length(cell))+1);
            r1 = locs(ix1, :);
            r2 = locs(ix2, :);
            dr = sqrt(sum((mod(r1 - r2 + L/2, L) - L/2).^2));
            P = P + dr;
            if dr < lowdist
                tries = tries + 1;
                [newcellsA, newcellsB] = transitionT1(newcells, ix1, ix2);
                %[EA, EL, EP] = cellenergies(A0, newcellsA, locs, L);
                %E1 = [EA, EL, EP] * coeffs';
                %[EA, EL, EP] = cellenergies(A0, newcellsB, locs, L);
                %E2 = [EA, EL, EP] * coeffs';
                [E1,~] = scaleenergy(locs,...
                    newcellsA,L,1,coeffs(1),coeffs(2),coeffs(3),A0);
                [E2,~] = scaleenergy(locs,...
                    newcellsB,L,1,coeffs(1),coeffs(2),coeffs(3),A0);
                if E1 < E0 && E1 < E2
                    newcells = newcellsA;
                    T1s = T1s + 1;
                    return
                elseif E2 < E0
                    newcells = newcellsB;
                    T1s = T1s + 1;
                    return
                end
            end
        end
        cix = cix + 1;
    end
end

