function proposal = block_Markov_chain_proposal(theta,COV,c,block)

% measure parameter dimension
n = length(theta);

% generate block covariance matrix by taking all off-block part zero
V_block = COV;
complementary_index = setdiff(linspace(1,n,n),block);
V_block(:,complementary_index) = 0;        
V_block(complementary_index,:) = 0;

% draw the (blockwise) proposal from normal transition density
proposal = mvnrnd(theta,c^2*V_block)';

end