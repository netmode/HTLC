function HBC = hyperBC( adjMatrix, coordinatesMatrix )
% ALGORITHM FOR HYPERBOLIC BETWEENNESS CENTRALITY
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
nodesNumber = ( size( coordinatesMatrix, 1 ) );
HBC = zeros( nodesNumber,1 );
for destination=1:5
    indexP = zeros( nodesNumber,1 );
    
    sigma = zeros( nodesNumber, 1);
    sigma( destination ) =1;
    distances = zeros(1,nodesNumber);
    %% STAGE 1 - TOPOLOGICAL SORT
    dst = coordinatesMatrix(destination,:);
    
    % a faster way to calculate nodes distance
    ysum = 1; %Sxi^2
    for j=2:9
        ysum = ysum + dst(j).^2;
    end
    for vertex=1:nodesNumber 
        %if ( vertex~=destination )
        xsum = 1; %Syi^2
        xysum = 0;
        
        for j=2:9
            xsum = xsum + coordinatesMatrix(vertex,j).^2;
            xysum = xysum + coordinatesMatrix(vertex,j)*dst(j);
        end

        t = sqrt( ysum*xsum ) - xysum;
        dist = acosh(t);
        distances(vertex) = dist;
        %end
    end
    %DAG = [DAG [destination; 0]];
    % sort dag by the distance of vertices from destination
    %DAG = sortrows(DAG',-2);
    %        DAG = DAG(:,1)';
    [~,DAG] = sort(distances, 'descend');
    %% PART 2 - 
    for i=nodesNumber:-1:1
        v = DAG(i);
        for j=1:indexMatrix(v)
            w = testMatrix(v,j);
            if ( distances(w) > distances(v) + 0.3 )
                sigma(w) = sigma(w)+sigma(v);
				indexP(w) = indexP(w)+1;
                P(w,indexP(w)) = v;
            end
        end
    end
    
    %PART 3
    delta = zeros( nodesNumber, 1 );
    % S returns vertices in order of non-increasing distance from s
    for node=1:nodesNumber
        % pop w<-S
        w = DAG( node );
        if sigma(w)>0
            for j=1:indexP(w)
                v = P(w,j);
                delta(v) = delta(v)+( sigma(v)/sigma(w) )*(1 + delta(w)) ;
            end
        end
        if ( w~=destination )
             HBC(w) = HBC(w)+delta(w);
        end
    end        
    
    
     
end
HBC = HBC/nodesNumber;
toc
end