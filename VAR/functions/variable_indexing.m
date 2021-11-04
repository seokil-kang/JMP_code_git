function variable_indexing(data_index)

% full variable list
var_list = ["y" "c" "i" "w" "n" "b" "ps" "g" "z" "Pi" "R" "RB" "x" "ps"];

for j = 1:length(data_index)
    assignin('caller',var_list(data_index(j)),j);
end
end