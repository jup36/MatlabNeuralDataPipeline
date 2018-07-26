% working with swc files

%% load the data
[inode,R,X,Y,Z,hD,idpar] = textread ('20100126c1-traced-2.swc', '%d%d%f%f%f%f%d', 'delimiter', ' ');

cellTree.node(1).x = X(1);
cellTree.node(1).y = Y(1);
cellTree.node(1).z = Z(1);
cellTree.node(1).hD= hD(1);
validNodes = 2;

size(inode,1)

for i = 2:size(inode,1)
    % eliminate the orphaned elements
    if idpar(i)>-1
        % this is a valid connection and thus a valid node
        cellTree.cnxns(validNodes,1)    = inode(i);
        cellTree.cnxns(validNodes,2)    = idpar(i);
        cellTree.node(validNodes).x     = X(i);
        cellTree.node(validNodes).y     = Y(i);
        cellTree.node(validNodes).z     = Z(i);
        cellTree.node(validNodes).hD    = hD(i);
        
        validNodes = validNodes+1;
        
    end

end

cellTree.validNodes = validNodes-1;

tmp = [cellTree.node(:).x,cellTree.node(:).y,cellTree.node(:).z];
cellTree.max = max(max(tmp));

%% plot a swc file

zProj       = 0;
thickThresh = 1;
isomorphic  = 0;

figure(1); clf;
hold on;

for i=1:cellTree.validNodes
    
    if thickThresh==1
        
        if cellTree.node(i).hD > 1.0            
            % find the real index of the parent
            parentOriginalID = cellTree.cnxns(i,2);
            parentNewID = find(cellTree.cnxns(:,1)==parentOriginalID,1);

            Xc = [cellTree.node(parentNewID).x , cellTree.node(i).x];
            Yc = [cellTree.node(parentNewID).y , cellTree.node(i).y];
            if zProj==1
                Zc = [1 , 1];
            else
                Zc = [cellTree.node(parentNewID).z , cellTree.node(i).z];
            end

            line(Xc,Yc,Zc,'Color',[1 0.251 0.251],'LineWidth',cellTree.node(i).hD);
        end
        
    else

        % find the real index of the parent
        parentOriginalID = cellTree.cnxns(i,2);
        parentNewID = find(cellTree.cnxns(:,1)==parentOriginalID,1);

        Xc = [cellTree.node(parentNewID).x , cellTree.node(i).x];
        Yc = [cellTree.node(parentNewID).y , cellTree.node(i).y];
        if zProj==1
            Zc = [1 , 1];
        else
            Zc = [cellTree.node(parentNewID).z , cellTree.node(i).z];
        end

        line(Xc,Yc,Zc,'Color',[0 0.53 0.8],'LineWidth',cellTree.node(i).hD);
    end
end

if isomorphic==1
    axis([0 cellTree.max 0 cellTree.max 0 cellTree.max]);
end
grid on;
view([-20 30]);

hold off;