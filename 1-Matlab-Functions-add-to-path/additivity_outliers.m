function [enrch_outlier,num_successes,p_values,x_thresh_per,y_thresh_per] = additivity_outliers(parent_seq_cell,additivity_tbl_ranked,rms)

%% find enriched outliers using all 144 permutations, according to the dest sequence (first column), threshold of outlier = higher than RMSD, my data

seq_for_enrch_calc   = parent_seq_cell;

thresh_outliers      = rms;
outliers_over_thresh = additivity_tbl_ranked(additivity_tbl_ranked.diff_vec > thresh_outliers,:);

sample_size          = height(outliers_over_thresh);
population_size      = height(additivity_tbl_ranked);
num_success_in_pop   = height(additivity_tbl_ranked)./length(seq_for_enrch_calc);

num_successes        = zeros(length(seq_for_enrch_calc),1);

for i = 1:length(seq_for_enrch_calc)

    num_successes(i,1) = sum(strcmp(outliers_over_thresh.dest_seq_vec,seq_for_enrch_calc(i)));

end

enrch_outlier = 100.*num_successes./sample_size;

p_values      = hygepdf(num_successes,population_size,num_success_in_pop,sample_size);
p_val_norm    = p_values.*length(seq_for_enrch_calc);


thresh_per   = 100*(1/length(seq_for_enrch_calc));
x_thresh_per = 0:0.5:(length(seq_for_enrch_calc)+1);
y_thresh_per = thresh_per.*ones(1,length(x_thresh_per));

end
