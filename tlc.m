% Calculate TLC metric
% INPUT: adjacency matrix
% OUTPUT: LC = TLC values of nodes
function LC = tlc(adjMatrix)
% Normal TLC ( Traffic Load Centrality )

nodesNumber = length( adjMatrix );
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
nodesNumber = length( indexMatrix );
tic
P = zeros( nodesNumber,pSize );
indexP = zeros(nodesNumber,1);
Q = zeros( 1, nodesNumber );
distance = ones( nodesNumber, 1).*(-1);
%intialize betw
LC = zeros( nodesNumber, 1 );
% for s c- V do
for destination=1:nodesNumber
    % S<- empty stack;
    Stack = zeros( 1, nodesNumber );
    stackIndex = nodesNumber;
    % P[w]<- empty list, w c- V;   
    distance( destination ) = 0;
    
    % Q is a queue
    qStart = 1;
    qEnd = 1;
    Q( qStart ) =  destination ;
    while ( qStart<=qEnd )
        
        v = Q( qStart );
        qStart = qStart+1;
        
        Stack(stackIndex) = v;
        stackIndex = stackIndex-1;
        
        neighboursV = testMatrix(v,1:indexMatrix(v));
        for n=1:indexMatrix(v) 
            w = neighboursV(n);
            % w found for the first time?
            if distance(w)<0
                qEnd = qEnd + 1;
                Q( qEnd ) = w;
                distance(w) = distance(v)+1;
            end
            % shortest path to w via v?
            if ( distance(w) ==  distance(v) + 1 )
                indexP(w) = indexP(w)+1;
                P(w,indexP(w)) = v;
            end
        end
    end
    delta = ones( nodesNumber, 1 );
    % S returns vertices in order of non-increasing distance from
    % destination
    for node=1:nodesNumber
        % pop w<-S
        w = Stack( node );
        
        for j=1:indexP(w)
            v = P(w,j);
            delta(v) = delta(v)+( delta(w)/indexP(w)) ;
        end
        indexP(w) = 0;
        distance(w) = -1;
       
    end
    LC = LC+delta;
end
%LC = LC / nodesNumber;
clear Stack Q distance w v source distance
toc
end