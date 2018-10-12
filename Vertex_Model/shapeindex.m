function index=shapeindex(cells,locs,L0)
    L=L0;
    index=zeros(1,length(cells));
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
        cell_edge=cell_drs(1:N,:)-cell_drs([2:N,1],:);
        cell_edge_len=sqrt(sum(cell_edge.^2,2));
        La=sum(cell_edge_len);
        index(ix)=La./sqrt(cell_area);
    end
    index=mean(index);
end