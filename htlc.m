% Calculate HTLC metric
% INPUT: adjacency matrix and coordinates matrix (1st column is nodeID, the rests are the coordinates)
% OUTPUT: RBC = HTLC values of nodes
function RBC = htlc( adjMatrix, coordinatesMatrix )
 
    %% Pre-Processing of Adjacency Matrix
    nodesNumber = size( adjMatrix, 2 );
    pSize = max(sum(adjMatrix));
    testMatrix = zeros(nodesNumber, pSize);
    indexMatrix = zeros(nodesNumber,1);
    for i=1:nodesNumber
        tempArray = find(adjMatrix(i,:));
        for j=1:length( tempArray )
             testMatrix(i,j) = tempArray(j);
        end
        indexMatrix(i) = length( tempArray );
    end
 
    clear adjMatrix tempArray     
 
    tic;
    % RBC initialization
    dimensions = size(coordinatesMatrix,2);
    successorsVi = zeros(1,pSize);
    nodesNumber = ( size( coordinatesMatrix, 1 ) );
    RBC = zeros( nodesNumber,1 );

    for destination=1:nodesNumber
        %destination
        distances = zeros(1,nodesNumber);
        %% STAGE 1 - TOPOLOGICAL SORT
        dst = coordinatesMatrix(destination,:);
        
        % a faster way to calculate nodes distance
        ysum = 1; %Sxi^2
        for j=2:dimensions
            ysum = ysum + dst(j).^2;
        end
        for vertex=1:nodesNumber
            xsum = 1; %Syi^2
            xysum = 0;
            
            for j=2:dimensions
                xsum = xsum + coordinatesMatrix(vertex,j).^2;
                xysum = xysum + coordinatesMatrix(vertex,j)*dst(j);
            end
            
            t = sqrt( ysum*xsum ) - xysum;
            dist = acosh(t);
            distances(vertex) = dist;
            %end
        end
        
        % sort dag by the distance of vertices from destination
        [~,DAG] = sort(distances, 'descend');
        
        %% STAGE 2 - INIT DELTA
        delta = ones(nodesNumber,1);
        %% STAGE 3 - ACCUMULATE d.,.(v)
        %     tempAdjacency = adjMatrix;
        for i=1:nodesNumber
            vi = DAG(i);
            distanceVi = distances(vi);
            
            sizeofSucVi = 0;
            for j=1:indexMatrix(vi)
                vj = testMatrix(vi, j);
                if ( distances(vj) < distanceVi )
                    sizeofSucVi = sizeofSucVi+1;
                    successorsVi( sizeofSucVi) = vj;
                end
            end
            
            
            % set the value of R(vi,vj)
            if ( sizeofSucVi ~=0 )
                R= 1/sizeofSucVi;
            else R=0;
            end
            
            for j=1:sizeofSucVi
                vj = successorsVi(j);
                delta(vj) = delta(vj)+delta(vi)*R;
            end
        end
        RBC = RBC + delta;
    end
    RBC = RBC/nodesNumber;
    toc
end