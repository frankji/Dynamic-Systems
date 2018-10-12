function [Fs, F_elastic, F_tension, F_contractility]=scaleforces(locs,...
    cells,L0,h, K,lambda,gamma,A0)
    L=L0;
    F_elastic=zeros(length(locs),2);
    F_tension=zeros(length(locs),2);
    F_contractility=zeros(length(locs),2);
    connections=zeros(length(locs));
    for ix = 1:length(cells)
        ixs=cells{ix};
        cell_locs = locs(ixs, :);
        cell_drs = cell_locs;
        for vix=1:length(cell_drs)
            cell_drs(vix, :) = mod(cell_locs(vix, :) - cell_locs(1, :) + L/2, L) - L/2;
        end
        N=length(ixs);
        mat=[0 1;-1 0]';
        cell_area=abs(sum(cell_drs(1:N,1).*cell_drs([2:N,1],2))-...
            sum(cell_drs(1:N,2).*cell_drs([2:N,1],1)))/2;
        cell_edge=cell_drs(1:N,:)-cell_drs([2:N,1],:);
        cell_edge_len=sqrt(sum(cell_edge.^2,2));
        cell_edge(:,1)=cell_edge(:,1)./cell_edge_len;
        cell_edge(:,2)=cell_edge(:,2)./cell_edge_len;
        La=sum(cell_edge_len);
        for i = 1:length(ixs)
            F_elastic(ixs(i),:)=F_elastic(ixs(i),:)+...
                0.5.*((h^2)*cell_area-A0).*...
                (cell_drs((mod(i-1-1,N)+1),:)-...
                cell_drs((mod(i+1-1,N)+1),:))*mat;
            if(connections(ixs(i),ixs(mod(i+1-1,N)+1))==0)
                F_tension(ixs(i),:)=F_tension(ixs(i),:)-...
                    cell_edge(i,:);
                connections(ixs(i),ixs(mod(i+1-1,N)+1))=1;
            end
            F_contractility(ixs(i),:)=F_contractility(ixs(i),:)-...
                    La*cell_edge(i,:);
            if(connections(ixs(i),ixs(mod(i-1-1,N)+1))==0)
                F_tension(ixs(i),:)=F_tension(ixs(i),:)+...
                    cell_edge(mod(i-1-1,N)+1,:);
                connections(ixs(i),ixs(mod(i-1-1,N)+1))=1;
            end
            F_contractility(ixs(i),:)=F_contractility(ixs(i),:)+...
                    La*cell_edge(mod(i-1-1,N)+1,:);
        end
    end
    F_elastic=K.*F_elastic;
    F_tension=h.*lambda.*F_tension;
    F_contractility=(h^2).*gamma.*F_contractility;
    Fs=F_elastic+F_tension+F_contractility;
end