function [block1 block2] = generate_random_block(n,N_block)

% block size depends on the # of blocks
block_size = floor(n/N_block);

% if block size cannot be equal, add the remainder to the first block
block_remainder = mod(n,block_size);

% randomize block at each iteration
block_full = randsample(linspace(1,n,n),n);

% indexing block in order
block_start = 1;
block_end = block_size+block_remainder;
block1 = sort(block_full(block_start:block_end));

% block baskets
block2=[];

for block_stacking = 2:N_block
    
    % update block index
    block_start = block_end + 1;
    block_end = block_end + block_size;
    
    % pick up block in order
    block_jth = sort(block_full(block_start:block_end));

    % stack up the block and its covariance matrix
    block2 = [block2; block_jth];
end

end