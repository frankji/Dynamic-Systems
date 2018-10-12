function [E, E_elastic, E_tension, E_contractility]=vertexenergy(locs,...
    cells,L0,K,lambda,gamma,A0)
    L=L0;
    E_elastic=zeros(1,length(cells));
    E_contractility=zeros(1,length(cells));
    lij=zeros(length(locs));
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
        E_elastic(ix)=0.5*(cell_area-A0)^2;
        cell_edge=cell_drs(1:N,:)-cell_drs([2:N,1],:);
        cell_edge_len=sqrt(sum(cell_edge.^2,2));
        La=sum(cell_edge_len);
        E_contractility(ix)=0.5*La^2;
        ijs=[ixs(1:N);ixs([2:N,1])]';
        for i = 1:length(ijs)
            lij(ijs(i,1),ijs(i,2))=cell_edge_len(i);
            lij(ijs(i,2),ijs(i,1))=cell_edge_len(i);
        end
    end
    ij_ix = triu(ones(length(locs), length(locs)), 1) > 0;
    E_tension=lambda*sum(lij(ij_ix));
    E_elastic=K.*sum(E_elastic);
    E_contractility=gamma.*sum(E_contractility);
    E=E_elastic+E_tension+E_contractility;
end