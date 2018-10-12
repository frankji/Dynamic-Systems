function [newcellsA, newcellsB] = transitionT1(cells, locAix, locBix)
    
    %     \    m1    /
    %      \        /
    %       \A     /
    %     e1 ------ e2
    %       /     B\
    %      /        \
    %     /    m2    \
    ixs1 = -ones(3, 1);
    ixs2 = -ones(3, 1);
    n1=1;
    n2=1;
    for cix=1:length(cells)
        if any(cells{cix} == locAix)
            if length(cells{cix}) <= 2
                % This cell needs to wait for a T2 transition.
                newcellsA = cells;
                newcellsB = cells;
                return
            end
            ixs1(n1)=cix;
            n1 = n1+1;
        end
        if any(cells{cix} == locBix)
            if length(cells{cix}) <= 2
                newcellsA = cells;
                newcellsB = cells;
                return
            end
            ixs2(n2)=cix;
            n2 = n2+1;
        end
    end
    
    assert(length(ixs1) == 3);
    assert(all(ixs1 > 0)); 
    assert(length(ixs2) == 3);
    assert(all(ixs2 > 0));
    
    midcells = intersect(ixs1, ixs2);
    endcells = setxor(ixs1, ixs2);
    assert(length(midcells) == 2);
    assert(length(endcells) == 2);
    % midcells are now the two cells that contain both locAix and locBix.
    % endcells only contain one or the other.
    
    % To keep things simple, let's make sure that endcells(1) has locAix,
    % and endcells(2) has locBix. If that's backward, then flip them.
    if any(cells{endcells(2)} == locAix)
        endcells = endcells(end:-1:1);
    end
    
    celle1 = cells{endcells(1)};
    celle2 = cells{endcells(2)};
    cellm1 = cells{midcells(1)};
    ixAe1 = find(celle1 == locAix);
    ixBe2 = find(celle2 == locBix);
    ixAm1 = find(cellm1 == locAix);
    ixBm1 = find(cellm1 == locBix);
    
    % Also to keep things simple, we need to make sure locAix comes before
    % locBix in midcells(1). If not, switch midcells. This ensures
    % that the order e1-m1-e2-m2 is clockwise.
    
    if mod(ixAm1, length(cellm1)) + 1 ~= ixBm1
        midcells = midcells(end:-1:1);
        cellm1 = cells{midcells(1)};
        ixAm1 = find(cellm1 == locAix);
        ixBm1 = find(cellm1 == locBix);
    end
    
       
    cellm2 = cells{midcells(2)};
    ixAm2 = find(cellm2 == locAix);
    ixBm2 = find(cellm2 == locBix);
    
    % Switch it to
    %     \ m1 /
    %      \  /
    %       \/    B
    %     e1 ------ e2
    %        A    /\
    %            /  \
    %           / m2 \
    
    newcellsA = cells(1:end);
    % insert B into e1 before A
    newcellsA{endcells(1)} = [celle1(1:ixAe1-1), locBix, celle1(ixAe1:end)];
    % insert A into e2 before B
    newcellsA{endcells(2)} = [celle2(1:ixBe2-1), locAix, celle2(ixBe2:end)];
    % remove B from m1
    newcellsA{midcells(1)} = [cellm1(1:ixBm1-1), cellm1(ixBm1+1:end)];
    % remove A from m2
    newcellsA{midcells(2)} = [cellm2(1:ixAm2-1), cellm2(ixAm2+1:end)];
    
%     if length(newcellsA{max(midcells)}) <= 2
%         newcellsA = [newcellsA(1:max(midcells)-1); newcellsA(max(midcells)+1:end)];
%     end
%     if length(newcellsA{min(midcells)}) <= 2
%         newcellsA = [newcellsA(1:min(midcells)-1); newcellsA(min(midcells)+1:end)];
%     end
    
    
    % Switch it to
    %           \ m1 /
    %            \  /
    %        A    \/
    %     e1 ------ e2
    %       /\    B
    %      /  \
    %     / m2 \
    
    newcellsB = cells(1:end);
    % insert B into e1 after A
    newcellsB{endcells(1)} = [celle1(1:ixAe1), locBix, celle1(ixAe1+1:end)];
    % insert A into e2 after B
    newcellsB{endcells(2)} = [celle2(1:ixBe2), locAix, celle2(ixBe2+1:end)];
    % remove A from m1
    newcellsB{midcells(1)} = [cellm1(1:ixAm1-1), cellm1(ixAm1+1:end)];
    % remove B from m2
    newcellsB{midcells(2)} = [cellm2(1:ixBm2-1), cellm2(ixBm2+1:end)];
    
    
%     if length(newcellsB{max(midcells)}) <= 2
%         newcellsB = [newcellsB(1:max(midcells)-1); newcellsB(max(midcells)+1:end)];
%     end
%     if length(newcellsB{min(midcells)}) <= 2
%         newcellsB = [newcellsB(1:min(midcells)-1); newcellsB(min(midcells)+1:end)];
%     end
    
    
%     clf;
%     drawcells(locs, cells([midcells; endcells]), L, 'ko-');
%     drawcells(locs+0.05, newcellsA([midcells; endcells]), L, 'ro-');
%     drawcells(locs-0.05, newcellsB([midcells; endcells]), L, 'bo-');
end