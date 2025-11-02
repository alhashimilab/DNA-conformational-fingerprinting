function seq_con_T = convSeqtoT(seq_con_G)

seq_cell_split   = cell(length(seq_con_G),3);
seq_cell_split_T = cell(length(seq_con_G),3);
seq_con_T        = cell(length(seq_con_G),1);

for i = 1:length(seq_con_G)
    
    seq_cell_split{i,1} = seq_con_G{i,1}(1);
    seq_cell_split{i,2} = seq_con_G{i,1}(2);
    seq_cell_split{i,3} = seq_con_G{i,1}(3);
    
    seq_cell_split_T{i,2} = 'T';
    
    if seq_cell_split{i,1} == 'A'
        seq_cell_split_T{i,3} = 'T';
    elseif seq_cell_split{i,1} == 'G'
        seq_cell_split_T{i,3} = 'C';
    elseif seq_cell_split{i,1} == 'C'
        seq_cell_split_T{i,3} = 'G';
    elseif seq_cell_split{i,1} == 'T'
        seq_cell_split_T{i,3} = 'A';
    end
    
    if seq_cell_split{i,3} == 'A'
        seq_cell_split_T{i,1} = 'T';
    elseif seq_cell_split{i,3} == 'G'
        seq_cell_split_T{i,1} = 'C';
    elseif seq_cell_split{i,3} == 'C'
        seq_cell_split_T{i,1} = 'G';
    elseif seq_cell_split{i,3} == 'T'
        seq_cell_split_T{i,1} = 'A';
    end
    
    seq_con_T{i,1}(1) = seq_cell_split_T{i,1};
    seq_con_T{i,1}(2) = seq_cell_split_T{i,2};
    seq_con_T{i,1}(3) = seq_cell_split_T{i,3};
    
end