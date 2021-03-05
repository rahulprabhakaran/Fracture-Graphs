function [chain_length] = compute_chain_length(chain, XY)
% this function computes the cumulative length of a chain that is specified
% by a graph node sequence 

    % reshaping chain into component edges
    chain = chain';
    chain(numel(chain)-1,2) = chain(numel(chain),1);
    chain(end,:) = [];
    chain(1:numel(chain(:,1))-1,2) = chain(2:numel(chain(:,1)),1); 

    for i=1:numel(chain(:,1))
        A(i,1) = distancePoints( [ XY(chain(i,1),1) XY(chain(i,1),2)], ...
         [XY(chain(i,2),1) XY(chain(i,2),2)],2 );
    end

    chain_length = sum(A);


end

