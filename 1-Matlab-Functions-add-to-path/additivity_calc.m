function [additivity_tbl_ranked,additivity_tbl_uniq,parent_seq_cell] = additivity_calc(sorted_seqs_T, cs_sorted, uniq_seqs_am)

seqs_T = sorted_seqs_T;

seqs_T_sort_by_name = {'GTG';'GTA';'GTT';'GTC';'ATG';'ATA';'ATT';'ATC';'TTG';'TTA';'TTT';'TTC';'CTG';'CTA';'CTT';'CTC'};

dw_pred  = zeros(length(seqs_T),9);
dw_exp   = zeros(length(seqs_T),9);
dw_r_mut = zeros(length(seqs_T),9);
dw_l_mut = zeros(length(seqs_T),9);

sequences_for_corr = cell(length(seqs_T),10);
sequences_left     = cell(length(seqs_T),9);
sequences_right    = cell(length(seqs_T),9);

sequences_for_corr(:,1) = seqs_T_sort_by_name;

for i = 1:length(seqs_T)
    
    ref_seq         = seqs_T_sort_by_name{i};
    [~,ind_ref_seq] = ismember(ref_seq,seqs_T);
    cs_ref_seq      = cs_sorted(ind_ref_seq);

    if ref_seq(1) == 'T'
        sub_5p = {'A','C','G'};
    elseif ref_seq(1) == 'A'
        sub_5p = {'T','C','G'};
    elseif ref_seq(1) == 'C'
        sub_5p = {'A','T','G'};
    elseif ref_seq(1) == 'G'
        sub_5p = {'A','T','C'};
    end
    
    if ref_seq(3) == 'T'
        sub_3p = {'A','C','G'};
    elseif ref_seq(3) == 'A'
        sub_3p = {'T','C','G'};
    elseif ref_seq(3) == 'C'
        sub_3p = {'A','T','G'};
    elseif ref_seq(3) == 'G'
        sub_3p = {'A','T','C'};
    end
    
    for j = 1:length(sub_5p)
        
        sub_5p_from_ref    = ref_seq;
        sub_5p_from_ref(1) = sub_5p{j};

        [~,ind_sub_5p] = ismember(sub_5p_from_ref,seqs_T);
        cs_sub_5p    = cs_sorted(ind_sub_5p);
        dw_exp_5p    = cs_sub_5p - cs_ref_seq;

        for k = 1:length(sub_3p)
            
            sub_3p_from_ref      = ref_seq;
            sub_3p_from_ref(3)   = sub_3p{k};

            [~,ind_sub_3p] = ismember(sub_3p_from_ref,seqs_T);
            cs_sub_3p    = cs_sorted(ind_sub_3p);
            dw_exp_3p    = cs_sub_3p - cs_ref_seq;
          
            dw_pred(i,k+3*(j-1))  = dw_exp_5p + dw_exp_3p;

            sub_exp      = ref_seq;
            sub_exp(1)   = sub_5p{j};
            sub_exp(3)   = sub_3p{k};
            
            [~,ind_exp]  = ismember(sub_exp,seqs_T);
            cs_exp     = cs_sorted(ind_exp);

            dw_exp(i,k+3*(j-1))  = cs_exp - cs_ref_seq;

            dw_l_mut(i,k+3*(j-1)) = dw_exp_5p;
            dw_r_mut(i,k+3*(j-1)) = dw_exp_3p;
        
            sequences_for_corr(i,k+3*(j-1)+1) = {sub_exp};
            sequences_left(i,k+3*(j-1))       = {sub_5p_from_ref};
            sequences_right(i,k+3*(j-1))      = {sub_3p_from_ref};

        end
        
    end

end

%% create vectors with 3 letter code

size_seq = size(sequences_for_corr);

parent_seq_cell = cell(size_seq(1),1);
dest_seq_cell   = cell(size_seq(1),(size_seq(2)-1));
right_mut_cell  = cell(size_seq(1),(size_seq(2)-1));
left_mut_cell   = cell(size_seq(1),(size_seq(2)-1));

for i = 1:size_seq(1)
   
    parent_seq_cell{i,1} = sequences_for_corr{i,1};

    for j = 2:size_seq(2)
        
        dest_seq_cell{i,j-1}  = sequences_for_corr{i,j};

        right_mut_cell{i,j-1} = sequences_right{i,j-1};

        left_mut_cell{i,j-1} = sequences_left{i,j-1}; 
        
    end
    
end

parent_seq_vec   = [parent_seq_cell;parent_seq_cell;parent_seq_cell;parent_seq_cell;parent_seq_cell;parent_seq_cell;parent_seq_cell;parent_seq_cell;parent_seq_cell];
dest_seq_vec     = reshape(dest_seq_cell,[],1);
right_mut_vec    = reshape(right_mut_cell,[],1);
left_mut_vec     = reshape(left_mut_cell,[],1);

obs_dw_vec       = reshape(dw_exp,[],1);
pred_dw_vec      = reshape(dw_pred,[],1);
right_mut_dw_vec = reshape(dw_r_mut,[],1);
left_mut_dw_vec  = reshape(dw_l_mut,[],1);

diff_vec         = abs(obs_dw_vec - pred_dw_vec);

additivity_tbl = table(dest_seq_vec,parent_seq_vec,...
    right_mut_vec,left_mut_vec,...
    pred_dw_vec,obs_dw_vec,...
    right_mut_dw_vec,left_mut_dw_vec,diff_vec);

[~,ind_rank]          = sort(additivity_tbl.diff_vec,'descend');
additivity_tbl_ranked = additivity_tbl(ind_rank,:);

%% compare with Akanksha's parent and dest sequences

ind_seq_uniq    = zeros(height(uniq_seqs_am),1);

for i = 1:height(uniq_seqs_am)
    
    for j = 1:length(dest_seq_vec)
        
        if strcmpi(uniq_seqs_am.dest_seq(i),dest_seq_vec(j)) && strcmpi(uniq_seqs_am.parent_seq(i),parent_seq_vec(j))
            ind_seq_uniq(i,1) = j;
        end
        
    end

end

parent_seq           = parent_seq_vec(ind_seq_uniq,1);
dest_seq             = dest_seq_vec(ind_seq_uniq,1);
right_mut            = right_mut_vec(ind_seq_uniq,1);
left_mut             = left_mut_vec(ind_seq_uniq,1);

obs_dw           = obs_dw_vec(ind_seq_uniq,1);
pred_dw          = pred_dw_vec(ind_seq_uniq,1);
right_mut_dw     = right_mut_dw_vec(ind_seq_uniq,1);
left_mut_dw      = left_mut_dw_vec(ind_seq_uniq,1);
diff             = diff_vec(ind_seq_uniq,1);

additivity_tbl_uniq = table(dest_seq,parent_seq,...
    right_mut,left_mut,...
    pred_dw,obs_dw,...
    right_mut_dw,left_mut_dw,diff);

end
