function [area, avg, num]=cellarea(cells,locs,L0)
    L=L0;
    num=unique(cellnum(cells));
    cell_num=cellnum(cells);
    numcount=zeros(1,length(num));
    area=zeros(1,length(num));
    for i = 1:length(cell_num)
        numcount(num==cell_num(i))=numcount(num==cell_num(i))+1;
    end
    for ix = 1:length(cells)
        ixs=cells{ix};
        cell_locs = locs(ixs, :);
        cell_drs = cell_locs;
        for vix=1:length(cell_drs)
            cell_drs(vix, :) = mod(cell_locs(vix, :) - cell_locs(1, :) + L/2, L) - L/2;
        end
        N=length(ixs);
        cell_area=abs(sum(cell_drs(1:N,1).*cell_drs([2:N,1],2))-...
            sum(cell_drs(1:N,2).*cell_drs([2:N,1],1)))/2;
        area(num==length(ixs))=area(num==length(ixs))+cell_area;
    end
    avg=sum(area)./length(cell_num);
    area=area./numcount;
end