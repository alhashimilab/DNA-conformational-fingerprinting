% ------------------------------------------------------------
% AUTHORSHIP
% ------------------------------------------------------------
% Author: Hashim Al-Hashimi Lab
% Created: October 2025
% ------------------------------------------------------------

clear;
clc;

%% set colors and axes properties

color_sig_CtoA = [0/255 200/255 255/255];
color_sig_CtoG = [0/255 0/255 0/255];
color_sig_CtoT = [255/255 0/255 0/255];
color_sig_TtoA = [206/255 206/255 206/255];
color_sig_TtoC = [146/255 208/255 80/255];
color_sig_TtoG = [235/255 190/255 190/255];

color_sig_comb = [255/255, 200/255, 80/255];

fontsize  = 28;
ticksize  = 0.015;
axeswidth = 2;
linewidth = 3;
markersize = 18;
err_linewidth = 2;
err_capsize   = 10;

x_line = -3.0:0.05:3.0;
y_line = x_line;

conf_int_per = 0.95;

x_ticks = 1:1:16; % ticks for plots

parent_path = '/Users/orsula1/OrsDocs/OrsDocs_since_Jan2025/2025-Sequence-dependence/2025-Data-Deposition/Mutational-Signatures-JSD/JSD-Matlab';

%% load conformational propensities

path_input     = 'Input-files/conf-sig-input-files-for-Matlab-JSDs';
pathname_input = sprintf('%s/%s',parent_path,path_input);

%% choose data set

dataset = 'GT_Anion'; % choose: 'GT_Anion' / 'AT_HG' / 'GC_HG' / 'AT_BaseOpening'

if strcmp(dataset,'GT_Anion')

    filename   = 'Anion_populations_19F_1C_all_pHs.xlsx';
    filepath   = sprintf('%s/%s',pathname_input,filename);    
    data_table = readtable(filepath);

    pKas_19F              = data_table.pKa_avg(~isnan(data_table.pKa_avg));
    pKas_19F_errs         = data_table.pKa_avg_err(~isnan(data_table.pKa_avg));
    seqs_G_19F            = data_table.Seq_G(~isnan(data_table.pKa_avg));

    [sorted_pKas_19F,ind] = sort(pKas_19F);
    sorted_pKas_19F_errs  = pKas_19F_errs(ind);
    sorted_seqs_G_19F     = seqs_G_19F(ind);

    seqs_T                = convSeqtoT(sorted_seqs_G_19F);

    pH_for_calc           = 7.4;
    conf_prop_T           = (10.^(pH_for_calc-sorted_pKas_19F))./(1+10.^(pH_for_calc-sorted_pKas_19F));
    
    path_output           = 'Output-JSDs/GT-anion-JSDs';

elseif strcmp(dataset,'AT_HG')

    filename   = 'AT_HG_Matlab.csv';
    filepath   = sprintf('%s/%s',pathname_input,filename);
    data_table = readtable(filepath);

    seqs_A      = data_table.Trinucleotide; % sequence context. "Trinucleotide" is the title of the column containing the sequence context names. Change this to the actual title you have in yout data file
    seqs_T      = convSeqAtoT(seqs_A); % sequence context with T/C as the center nucleotide. Comment this line out if your sequences are already centered on the T/C
    conf_prop_T = data_table.pB; % propensities. "pB" is the title of the column containing the conformational propensities. Change this to the actual title you have in yout data file. Also - check if these are normalized, if not - they will be normalized by the script

    path_output = 'Output-JSDs/AT-HG-JSDs';

elseif strcmp(dataset,'GC_HG')

    filename   = 'GC_HG_Matlab.csv';
    filepath   = sprintf('%s/%s',pathname_input,filename);
    data_table = readtable(filepath);

    seqs_G      = data_table.Seq_G; % sequence context. "Seq_G" is the title of the column containing the sequence context names. Change this to the actual title you have in yout data file
    seqs_T      = convSeqtoT(seqs_G); % sequence context with T/C as the center nucleotide. Comment this line out if your sequences are already centered on the T/C
    conf_prop_T = data_table.pB; % propensities. "pB" is the title of the column containing the conformational propensities. Change this to the actual title you have in yout data file. Also - check if these are normalized, if not - they will be normalized by the script

    path_output = 'Output-JSDs/GC-HG-JSDs';

elseif strcmp(dataset,'AT_BaseOpening')

    filename   = 'AT_BaseOpening_Matlab.csv';
    filepath   = sprintf('%s/%s',pathname_input,filename);
    data_table = readtable(filepath);

    seqs_A      = data_table.Trinucleotide; % sequence context. "Trinucleotide" is the title of the column containing the sequence context names. Change this to the actual title you have in yout data file
    seqs_T      = convSeqAtoT(seqs_A); % sequence context with T/C as the center nucleotide. Comment this line out if your sequences are already centered on the T/C
    conf_prop_T = data_table.pB; % propensities. "pB" is the title of the column containing the conformational propensities. Change this to the actual title you have in yout data file. Also - check if these are normalized, if not - they will be normalized by the script

    path_output = 'Output-JSDs/AT-Base-Open-JSDs';

end


%% load mutational signatures - GRCh37

path_input     = 'Input-files/mut-sig-input-files-for-Matlab-JSDs';
pathname_input = sprintf('%s/%s',parent_path,path_input); % file name of SBS mutational signatures (downloaded from the COSMIC data base)
filename_mut_sig = 'GRCh37_SBS_profiles.txt';

filepath_mut_sig = sprintf('%s/%s',pathname_input,filename_mut_sig);

mut_sig_profiles = readtable(filepath_mut_sig);

%% rearrange signatures according to type of substitutions

ind              = zeros(height(mut_sig_profiles),1);
mut_sig_profiles = addvars(mut_sig_profiles,ind,'Before',"Type");
seqs_mut         = cell(height(mut_sig_profiles),1);
num_trans_mut    = 6; % number of single base substitution (transition mutations)

for i = 1:height(mut_sig_profiles)

    temp_char = mut_sig_profiles.Type(i);

    if all((temp_char{1,1}(3:5)) == 'C>A')
        mut_sig_profiles.ind(i) = 1;
    elseif all((temp_char{1,1}(3:5)) == 'C>G')
        mut_sig_profiles.ind(i) = 2;
    elseif all((temp_char{1,1}(3:5)) == 'C>T')
        mut_sig_profiles.ind(i) = 3;
    elseif all((temp_char{1,1}(3:5)) == 'T>A')
        mut_sig_profiles.ind(i) = 4;
    elseif all((temp_char{1,1}(3:5)) == 'T>C')
        mut_sig_profiles.ind(i) = 5;
    elseif all((temp_char{1,1}(3:5)) == 'T>G')
        mut_sig_profiles.ind(i) = 6;
    end

    seqs_mut{i,1} = [temp_char{1,1}(1),'T',temp_char{1,1}(7)];

end

mut_sig_profiles         = addvars(mut_sig_profiles,seqs_mut,'Before',"SBS1");
mut_sig_profiles_ordered = sortrows(mut_sig_profiles,mut_sig_profiles.ind);

CtoA = mut_sig_profiles_ordered((mut_sig_profiles_ordered.ind(:) == 1),3:end);
CtoG = mut_sig_profiles_ordered((mut_sig_profiles_ordered.ind(:) == 2),3:end);
CtoT = mut_sig_profiles_ordered((mut_sig_profiles_ordered.ind(:) == 3),3:end);
TtoA = mut_sig_profiles_ordered((mut_sig_profiles_ordered.ind(:) == 4),3:end);
TtoC = mut_sig_profiles_ordered((mut_sig_profiles_ordered.ind(:) == 5),3:end);
TtoG = mut_sig_profiles_ordered((mut_sig_profiles_ordered.ind(:) == 6),3:end);

seqs_mut_sig_16 = CtoT.seqs_mut;

%% only compare sequence contexts which appear in the conformational propensities list

[~,ind_seqs_T] = ismember(seqs_T,seqs_mut_sig_16);

%% normalize conformational propensities

Q_T = conf_prop_T; % using the pyrimidine as the middle of the triplet
Q_T = Q_T ./sum(Q_T); % normalizing the conformational propensities. Not needed if they are already normalized

%% correcting indices for plot

[~,ind_cell_T]       = ismember(seqs_mut_sig_16,seqs_T);
seqs_T_cell          = [seqs_T;seqs_mut_sig_16(ind_cell_T==0)];
conf_prop_added_T    = [Q_T;zeros(length(CtoT.seqs_mut(ind_cell_T==0)),1)];
[~,ind_for_plot_T]   = ismember(seqs_mut_sig_16,seqs_T_cell);

seqs_T_sort          = seqs_T_cell(ind_for_plot_T);
conf_prop_for_plot_T = conf_prop_added_T(ind_for_plot_T);

seqs_mut_for_plot    = seqs_mut_sig_16;

%% compute JSD values

js_div_T = zeros(num_trans_mut,(width(CtoT)-1));

for j = 1:num_trans_mut

    if j == 1
        mut_sig          = CtoA;
    elseif j == 2
        mut_sig          = CtoG;
    elseif j == 3
        mut_sig          = CtoT;
    elseif j == 4
        mut_sig          = TtoA;
    elseif j == 5
        mut_sig          = TtoC;
    elseif j == 6
        mut_sig          = TtoG;
    end

    for i = 1:width(js_div_T)

        P1 = table2array(mut_sig(:,i+1));
        P  = P1(ind_seqs_T,:);
        P  = P ./sum(P); % normalizing mutational probabilities

        KL_T = kldiv(seqs_T,P,Q_T,'js'); % using kldiv from Matlab Central file exchange to claculate JSD values. Cite as: David, F. KLDIV. MATLAB Central File Exchange. https://www.mathworks.com/matlabcentral/fileexchange/13089-kldiv (accessed 05/17/2024)
        js_div_T(j,i) = KL_T;

    end

end


%% compute JSDs between conformational fingerprint and a random distribution

%% create matrix of all probabilities (normalized), for a random distribution

CtoA_ord_by_conf = table2array(CtoA(ind_seqs_T,2:end));
CtoG_ord_by_conf = table2array(CtoG(ind_seqs_T,2:end));
CtoT_ord_by_conf = table2array(CtoT(ind_seqs_T,2:end));
TtoA_ord_by_conf = table2array(TtoA(ind_seqs_T,2:end));
TtoC_ord_by_conf = table2array(TtoC(ind_seqs_T,2:end));
TtoG_ord_by_conf = table2array(TtoG(ind_seqs_T,2:end));

CtoA_ord_by_conf_norm = CtoA_ord_by_conf./sum(CtoA_ord_by_conf);
CtoG_ord_by_conf_norm = CtoG_ord_by_conf./sum(CtoG_ord_by_conf);
CtoT_ord_by_conf_norm = CtoT_ord_by_conf./sum(CtoT_ord_by_conf);
TtoA_ord_by_conf_norm = TtoA_ord_by_conf./sum(TtoA_ord_by_conf);
TtoC_ord_by_conf_norm = TtoC_ord_by_conf./sum(TtoC_ord_by_conf);
TtoG_ord_by_conf_norm = TtoG_ord_by_conf./sum(TtoG_ord_by_conf);

mut_sig_comb_ord_by_conf_norm = [CtoA_ord_by_conf_norm,CtoG_ord_by_conf_norm,...
    CtoT_ord_by_conf_norm,TtoA_ord_by_conf_norm,TtoC_ord_by_conf_norm,TtoG_ord_by_conf_norm];

%% create the random distribution of probabilities and calcualte random JSDs

num_perms = 1000000; % number of permutations for random distributions

mut_sig_rand = zeros(length(seqs_T),1);
js_rand = zeros(1,num_perms);

tic

for i = 1:num_perms

    for j = 1:length(mut_sig_rand)

        x_seq             = mut_sig_comb_ord_by_conf_norm(j,:);
        mut_sig_rand(j,1) = randsample(x_seq,1);

    end

    mut_sig_rand_norm = mut_sig_rand./sum(mut_sig_rand);

    js_rand(1,i) = kldiv(seqs_T,mut_sig_rand_norm,Q_T,'js');

end

toc

%% check if the randomized JSD vector has a normal distribution

data = js_rand;

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

numBins = floor(sqrt(length(data))); % Number of bins
histogram(data, numBins, 'Normalization', 'pdf'); % Plot normalized histogram
hold on;
mu = mean(data);
sigma = std(data);
x = linspace(min(data), max(data), 1000);
plot(x, normpdf(x, mu, sigma), 'LineWidth', 2); % Plot the theoretical Gaussian PDF
hold off;

propedit


%% unbiased search - compare conformational propensities to individual transition mutations

js_threshold = 0.09; % JS threshold (arbitrary), can be changed
FDR_threshold = 0.05;

%% only C-->A transition mutations (ind = 1)

ind_CtoA                         = 1;
mut_sig_for_plot                 = CtoA(:,2:end);

js_CtoA_T                        = js_div_T(ind_CtoA,:);
[js_CtoA_sorted_T,js_CtoA_ind_T] = sort(js_CtoA_T);

num_sigs = width(mut_sig_for_plot);

JS_pVal_tbl_CtoA_T = cell(num_sigs,9);

for i = 1:num_sigs

    ind_sig_min_T          = js_CtoA_ind_T(i);

    rel_mut_sig_for_plot = mut_sig_for_plot(:,ind_sig_min_T);
    norm_mut_sig         = rel_mut_sig_for_plot{:,:}./sum(rel_mut_sig_for_plot{:,:});

    rel_mut_sig_for_corr_plot = mut_sig_for_plot(ind_seqs_T,ind_sig_min_T);
    norm_mut_sig_for_corr_plot = rel_mut_sig_for_corr_plot{:,:}./sum(rel_mut_sig_for_corr_plot{:,:});

    [r_corr,Rsq,Rsq_det,rms,a_rms_per,bestfit,x2,inBetween,slope,intercept,slope_err,intercept_err,p_corr] = corr_plot(Q_T,norm_mut_sig_for_corr_plot,x_line,conf_int_per);

    p_rand = sum(js_rand<=js_CtoA_sorted_T(i))./num_perms;

    JS_pVal_tbl_CtoA_T{i,1} = mut_sig_for_plot(:,ind_sig_min_T).Properties.VariableNames{1,1};
    JS_pVal_tbl_CtoA_T{i,2} = js_CtoA_sorted_T(i);
    JS_pVal_tbl_CtoA_T{i,3} = r_corr;
    JS_pVal_tbl_CtoA_T{i,4} = p_corr;
    JS_pVal_tbl_CtoA_T{i,5} = rms;
    JS_pVal_tbl_CtoA_T{i,6} = bestfit;
    JS_pVal_tbl_CtoA_T{i,7} = [slope,slope_err];
    JS_pVal_tbl_CtoA_T{i,8} = [intercept,intercept_err];
    JS_pVal_tbl_CtoA_T{i,9} = p_rand;

end

JS_tbl_CtoA_T = cell2table(JS_pVal_tbl_CtoA_T,...
    "VariableNames",["SBS" "JS" "r_pearson" "p_corr" "RMSD" "bestfit" "slope" "intercept" "p_rand"]);
p_values  = JS_tbl_CtoA_T.p_rand;
FDR       = mafdr(p_values,'BHFDR',false);
FDR_BH    = mafdr(p_values,'BHFDR',true);
JS_tbl_CtoA_T = addvars(JS_tbl_CtoA_T,FDR,FDR_BH,'After','p_rand');

best_ind_T = find(js_CtoA_sorted_T<js_threshold&JS_tbl_CtoA_T.FDR_BH'<FDR_threshold); % number of mutationl signatures fulfulling the JS threshold, as well as the FDR threshold
num_mut_sig_T = length(best_ind_T);
best_similarities_CtoA_T = JS_tbl_CtoA_T(best_ind_T,:); % to find which mutational signature has the best JSD scores for this transition mutation

for i = 1:num_mut_sig_T

    ind_sig_min_T          = js_CtoA_ind_T(i);

    rel_mut_sig_for_plot = mut_sig_for_plot(:,ind_sig_min_T);
    norm_mut_sig         = rel_mut_sig_for_plot{:,:}./sum(rel_mut_sig_for_plot{:,:});
    bar_vals             = [norm_mut_sig,conf_prop_for_plot_T];

    figure('units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'Color','w');

    hb1 = bar(bar_vals);

    hb1(1).FaceColor = color_sig_CtoA;
    hb1(1).EdgeColor = 'k';
    hb1(1).LineWidth = 1;

    hb1(2).FaceColor = 'b';
    hb1(2).EdgeColor = 'k';
    hb1(2).LineWidth = 1;

    axis([0.0 17.0 0.0 0.3]);

    set(gca, 'FontSize', fontsize, 'XminorTick', 'off', 'YminorTick', 'off','LineWidth',linewidth,...
        'XTick',x_ticks,'TickLength',[ticksize,ticksize]);

    xticklabels(seqs_mut_for_plot);

    ylabel('Probability', 'FontSize', 36);

    legend({sprintf('C-->A probability from %s',mut_sig_for_plot(:,ind_sig_min_T).Properties.VariableNames{1,1}),...
        'conformational probability from NMR'});

    text(1,0.25,sprintf('JS divergence = %0.3f', js_CtoA_sorted_T(i)),'FontSize',fontsize);

    propedit;



end


%% only C-->G transition mutations (ind = 2)

ind_CtoG                     = 2;
mut_sig_for_plot             = CtoG(:,2:end);

js_CtoG_T                        = js_div_T(ind_CtoG,:);
[js_CtoG_sorted_T,js_CtoG_ind_T] = sort(js_CtoG_T);

num_sigs = width(mut_sig_for_plot);

JS_pVal_tbl_CtoG_T = cell(num_sigs,9);

for i = 1:num_sigs

    ind_sig_min_T          = js_CtoG_ind_T(i);

    rel_mut_sig_for_plot = mut_sig_for_plot(:,ind_sig_min_T);
    norm_mut_sig         = rel_mut_sig_for_plot{:,:}./sum(rel_mut_sig_for_plot{:,:});

    rel_mut_sig_for_corr_plot = mut_sig_for_plot(ind_seqs_T,ind_sig_min_T);
    norm_mut_sig_for_corr_plot = rel_mut_sig_for_corr_plot{:,:}./sum(rel_mut_sig_for_corr_plot{:,:});

    [r_corr,Rsq,Rsq_det,rms,a_rms_per,bestfit,x2,inBetween,slope,intercept,slope_err,intercept_err,p_corr] = corr_plot(Q_T,norm_mut_sig_for_corr_plot,x_line,conf_int_per);

    p_rand = sum(js_rand<=js_CtoG_sorted_T(i))./num_perms;

    JS_pVal_tbl_CtoG_T{i,1} = mut_sig_for_plot(:,ind_sig_min_T).Properties.VariableNames{1,1};
    JS_pVal_tbl_CtoG_T{i,2} = js_CtoG_sorted_T(i);
    JS_pVal_tbl_CtoG_T{i,3} = r_corr;
    JS_pVal_tbl_CtoG_T{i,4} = p_corr;
    JS_pVal_tbl_CtoG_T{i,5} = rms;
    JS_pVal_tbl_CtoG_T{i,6} = bestfit;
    JS_pVal_tbl_CtoG_T{i,7} = [slope,slope_err];
    JS_pVal_tbl_CtoG_T{i,8} = [intercept,intercept_err];
    JS_pVal_tbl_CtoG_T{i,9} = p_rand;

end

JS_tbl_CtoG_T = cell2table(JS_pVal_tbl_CtoG_T,...
    "VariableNames",["SBS" "JS" "r_pearson" "p_corr" "RMSD" "bestfit" "slope" "intercept" "p_rand"]);
p_values  = JS_tbl_CtoG_T.p_rand;
FDR       = mafdr(p_values,'BHFDR',false);
FDR_BH    = mafdr(p_values,'BHFDR',true);
JS_tbl_CtoG_T = addvars(JS_tbl_CtoG_T,FDR,FDR_BH,'After','p_rand');

best_ind_T = find(js_CtoG_sorted_T<js_threshold&JS_tbl_CtoG_T.FDR_BH'<FDR_threshold); % number of mutationl signatures fulfulling the JS threshold, as well as the FDR threshold
num_mut_sig_T = length(best_ind_T);
best_similarities_CtoG_T = JS_tbl_CtoG_T(best_ind_T,:); % to find which mutational signature has the best JSD scores for this transition mutation

for i = 1:num_mut_sig_T

    ind_sig_min_T          = js_CtoG_ind_T(i);

    rel_mut_sig_for_plot = mut_sig_for_plot(:,ind_sig_min_T);
    norm_mut_sig         = rel_mut_sig_for_plot{:,:}./sum(rel_mut_sig_for_plot{:,:});
    bar_vals             = [norm_mut_sig,conf_prop_for_plot_T];

    figure('units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'Color','w');

    hb1 = bar(bar_vals);

    hb1(1).FaceColor = color_sig_CtoG;
    hb1(1).EdgeColor = 'k';
    hb1(1).LineWidth = 1;

    hb1(2).FaceColor = 'b';
    hb1(2).EdgeColor = 'k';
    hb1(2).LineWidth = 1;

    axis([0.0 17.0 0.0 0.3]);

    set(gca, 'FontSize', fontsize, 'XminorTick', 'off', 'YminorTick', 'off','LineWidth',linewidth,...
        'XTick',x_ticks,'TickLength',[ticksize,ticksize]);

    xticklabels(seqs_mut_for_plot);

    ylabel('Probability', 'FontSize', 36);

    legend({sprintf('C-->G probability from %s',mut_sig_for_plot(:,ind_sig_min_T).Properties.VariableNames{1,1}),...
        'conformational probability from NMR'});

    text(1,0.25,sprintf('JS divergence = %0.3f', js_CtoG_sorted_T(i)),'FontSize',fontsize);

    propedit;


end



%% only C-->T transition mutations (ind = 3)

ind_CtoT                     = 3;
mut_sig_for_plot             = CtoT(:,2:end);

js_CtoT_T                        = js_div_T(ind_CtoT,:);
[js_CtoT_sorted_T,js_CtoT_ind_T] = sort(js_CtoT_T);

num_sigs = width(mut_sig_for_plot);

JS_pVal_tbl_CtoT_T = cell(num_sigs,9);

for i = 1:num_sigs

    ind_sig_min_T          = js_CtoT_ind_T(i);

    rel_mut_sig_for_plot = mut_sig_for_plot(:,ind_sig_min_T);
    norm_mut_sig         = rel_mut_sig_for_plot{:,:}./sum(rel_mut_sig_for_plot{:,:});

    rel_mut_sig_for_corr_plot = mut_sig_for_plot(ind_seqs_T,ind_sig_min_T);
    norm_mut_sig_for_corr_plot = rel_mut_sig_for_corr_plot{:,:}./sum(rel_mut_sig_for_corr_plot{:,:});

    [r_corr,Rsq,Rsq_det,rms,a_rms_per,bestfit,x2,inBetween,slope,intercept,slope_err,intercept_err,p_corr] = corr_plot(Q_T,norm_mut_sig_for_corr_plot,x_line,conf_int_per);

    p_rand = sum(js_rand<=js_CtoT_sorted_T(i))./num_perms;

    JS_pVal_tbl_CtoT_T{i,1} = mut_sig_for_plot(:,ind_sig_min_T).Properties.VariableNames{1,1};
    JS_pVal_tbl_CtoT_T{i,2} = js_CtoT_sorted_T(i);
    JS_pVal_tbl_CtoT_T{i,3} = r_corr;
    JS_pVal_tbl_CtoT_T{i,4} = p_corr;
    JS_pVal_tbl_CtoT_T{i,5} = rms;
    JS_pVal_tbl_CtoT_T{i,6} = bestfit;
    JS_pVal_tbl_CtoT_T{i,7} = [slope,slope_err];
    JS_pVal_tbl_CtoT_T{i,8} = [intercept,intercept_err];
    JS_pVal_tbl_CtoT_T{i,9} = p_rand;

end

JS_tbl_CtoT_T = cell2table(JS_pVal_tbl_CtoT_T,...
    "VariableNames",["SBS" "JS" "r_pearson" "p_corr" "RMSD" "bestfit" "slope" "intercept" "p_rand"]);
p_values  = JS_tbl_CtoT_T.p_rand;
FDR       = mafdr(p_values,'BHFDR',false);
FDR_BH    = mafdr(p_values,'BHFDR',true);
JS_tbl_CtoT_T = addvars(JS_tbl_CtoT_T,FDR,FDR_BH,'After','p_rand');

best_ind_T = find(js_CtoT_sorted_T<js_threshold&JS_tbl_CtoT_T.FDR_BH'<FDR_threshold); % number of mutationl signatures fulfulling the JS threshold, as well as the FDR threshold
num_mut_sig_T = length(best_ind_T);
best_similarities_CtoT_T = JS_tbl_CtoT_T(best_ind_T,:); % to find which mutational signature has the best JSD scores for this transition mutation

for i = 1:num_mut_sig_T

    ind_sig_min_T          = js_CtoT_ind_T(i);

    rel_mut_sig_for_plot = mut_sig_for_plot(:,ind_sig_min_T);
    norm_mut_sig         = rel_mut_sig_for_plot{:,:}./sum(rel_mut_sig_for_plot{:,:});
    bar_vals             = [norm_mut_sig,conf_prop_for_plot_T];

    figure('units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'Color','w');

    hb1 = bar(bar_vals);

    hb1(1).FaceColor = color_sig_CtoT;
    hb1(1).EdgeColor = 'k';
    hb1(1).LineWidth = 1;

    hb1(2).FaceColor = 'b';
    hb1(2).EdgeColor = 'k';
    hb1(2).LineWidth = 1;

    axis([0.0 17.0 0.0 0.3]);

    set(gca, 'FontSize', fontsize, 'XminorTick', 'off', 'YminorTick', 'off','LineWidth',linewidth,...
        'XTick',x_ticks,'TickLength',[ticksize,ticksize]);

    xticklabels(seqs_mut_for_plot);

    ylabel('Probability', 'FontSize', 36);

    legend({sprintf('C-->T probability from %s',mut_sig_for_plot(:,ind_sig_min_T).Properties.VariableNames{1,1}),...
        'conformational probability from NMR'});

    text(1,0.25,sprintf('JS divergence = %0.3f', js_CtoT_sorted_T(i)),'FontSize',fontsize);

    propedit;


end



%% only T-->A transition mutations (ind = 4)

ind_TtoA                     = 4;
mut_sig_for_plot             = TtoA(:,2:end);

js_TtoA_T                        = js_div_T(ind_TtoA,:);
[js_TtoA_sorted_T,js_TtoA_ind_T] = sort(js_TtoA_T);

num_sigs = width(mut_sig_for_plot);

JS_pVal_tbl_TtoA_T = cell(num_sigs,9);

for i = 1:num_sigs

    ind_sig_min_T          = js_TtoA_ind_T(i);

    rel_mut_sig_for_plot = mut_sig_for_plot(:,ind_sig_min_T);
    norm_mut_sig         = rel_mut_sig_for_plot{:,:}./sum(rel_mut_sig_for_plot{:,:});

    rel_mut_sig_for_corr_plot = mut_sig_for_plot(ind_seqs_T,ind_sig_min_T);
    norm_mut_sig_for_corr_plot = rel_mut_sig_for_corr_plot{:,:}./sum(rel_mut_sig_for_corr_plot{:,:});

    [r_corr,Rsq,Rsq_det,rms,a_rms_per,bestfit,x2,inBetween,slope,intercept,slope_err,intercept_err,p_corr] = corr_plot(Q_T,norm_mut_sig_for_corr_plot,x_line,conf_int_per);

    p_rand = sum(js_rand<=js_TtoA_sorted_T(i))./num_perms;

    JS_pVal_tbl_TtoA_T{i,1} = mut_sig_for_plot(:,ind_sig_min_T).Properties.VariableNames{1,1};
    JS_pVal_tbl_TtoA_T{i,2} = js_TtoA_sorted_T(i);
    JS_pVal_tbl_TtoA_T{i,3} = r_corr;
    JS_pVal_tbl_TtoA_T{i,4} = p_corr;
    JS_pVal_tbl_TtoA_T{i,5} = rms;
    JS_pVal_tbl_TtoA_T{i,6} = bestfit;
    JS_pVal_tbl_TtoA_T{i,7} = [slope,slope_err];
    JS_pVal_tbl_TtoA_T{i,8} = [intercept,intercept_err];
    JS_pVal_tbl_TtoA_T{i,9} = p_rand;

end

JS_tbl_TtoA_T = cell2table(JS_pVal_tbl_TtoA_T,...
    "VariableNames",["SBS" "JS" "r_pearson" "p_corr" "RMSD" "bestfit" "slope" "intercept" "p_rand"]);
p_values  = JS_tbl_TtoA_T.p_rand;
FDR       = mafdr(p_values,'BHFDR',false);
FDR_BH    = mafdr(p_values,'BHFDR',true);
JS_tbl_TtoA_T = addvars(JS_tbl_TtoA_T,FDR,FDR_BH,'After','p_rand');

best_ind_T = find(js_TtoA_sorted_T<js_threshold&JS_tbl_TtoA_T.FDR_BH'<FDR_threshold); % number of mutationl signatures fulfulling the JS threshold, as well as the FDR threshold
num_mut_sig_T = length(best_ind_T);
best_similarities_TtoA_T = JS_tbl_TtoA_T(best_ind_T,:); % to find which mutational signature has the best JSD scores for this transition mutation

for i = 1:num_mut_sig_T

    ind_sig_min_T          = js_TtoA_ind_T(i);

    rel_mut_sig_for_plot = mut_sig_for_plot(:,ind_sig_min_T);
    norm_mut_sig         = rel_mut_sig_for_plot{:,:}./sum(rel_mut_sig_for_plot{:,:});
    bar_vals             = [norm_mut_sig,conf_prop_for_plot_T];

    figure('units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'Color','w');

    hb1 = bar(bar_vals);

    hb1(1).FaceColor = color_sig_TtoA;
    hb1(1).EdgeColor = 'k';
    hb1(1).LineWidth = 1;

    hb1(2).FaceColor = 'b';
    hb1(2).EdgeColor = 'k';
    hb1(2).LineWidth = 1;

    axis([0.0 17.0 0.0 0.3]);

    set(gca, 'FontSize', fontsize, 'XminorTick', 'off', 'YminorTick', 'off','LineWidth',linewidth,...
        'XTick',x_ticks,'TickLength',[ticksize,ticksize]);

    xticklabels(seqs_mut_for_plot);

    ylabel('Probability', 'FontSize', 36);

    legend({sprintf('T-->A probability from %s',mut_sig_for_plot(:,ind_sig_min_T).Properties.VariableNames{1,1}),...
        'conformational probability from NMR'});

    text(1,0.25,sprintf('JS divergence = %0.3f', js_TtoA_sorted_T(i)),'FontSize',fontsize);

    propedit;


end




%% only T-->C transition mutations (ind = 5)

ind_TtoC                     = 5;
mut_sig_for_plot             = TtoC(:,2:end);

js_TtoC_T                        = js_div_T(ind_TtoC,:);
[js_TtoC_sorted_T,js_TtoC_ind_T] = sort(js_TtoC_T);

num_sigs = width(mut_sig_for_plot);

JS_pVal_tbl_TtoC_T = cell(num_sigs,9);

for i = 1:num_sigs

    ind_sig_min_T          = js_TtoC_ind_T(i);

    rel_mut_sig_for_plot = mut_sig_for_plot(:,ind_sig_min_T);
    norm_mut_sig         = rel_mut_sig_for_plot{:,:}./sum(rel_mut_sig_for_plot{:,:});

    rel_mut_sig_for_corr_plot = mut_sig_for_plot(ind_seqs_T,ind_sig_min_T);
    norm_mut_sig_for_corr_plot = rel_mut_sig_for_corr_plot{:,:}./sum(rel_mut_sig_for_corr_plot{:,:});

    [r_corr,Rsq,Rsq_det,rms,a_rms_per,bestfit,x2,inBetween,slope,intercept,slope_err,intercept_err,p_corr] = corr_plot(Q_T,norm_mut_sig_for_corr_plot,x_line,conf_int_per);

    p_rand = sum(js_rand<=js_TtoC_sorted_T(i))./num_perms;

    JS_pVal_tbl_TtoC_T{i,1} = mut_sig_for_plot(:,ind_sig_min_T).Properties.VariableNames{1,1};
    JS_pVal_tbl_TtoC_T{i,2} = js_TtoC_sorted_T(i);
    JS_pVal_tbl_TtoC_T{i,3} = r_corr;
    JS_pVal_tbl_TtoC_T{i,4} = p_corr;
    JS_pVal_tbl_TtoC_T{i,5} = rms;
    JS_pVal_tbl_TtoC_T{i,6} = bestfit;
    JS_pVal_tbl_TtoC_T{i,7} = [slope,slope_err];
    JS_pVal_tbl_TtoC_T{i,8} = [intercept,intercept_err];
    JS_pVal_tbl_TtoC_T{i,9} = p_rand;

end

JS_tbl_TtoC_T = cell2table(JS_pVal_tbl_TtoC_T,...
    "VariableNames",["SBS" "JS" "r_pearson" "p_corr" "RMSD" "bestfit" "slope" "intercept" "p_rand"]);
p_values  = JS_tbl_TtoC_T.p_rand;
FDR       = mafdr(p_values,'BHFDR',false);
FDR_BH    = mafdr(p_values,'BHFDR',true);
JS_tbl_TtoC_T = addvars(JS_tbl_TtoC_T,FDR,FDR_BH,'After','p_rand');

best_ind_T = find(js_TtoC_sorted_T<js_threshold&JS_tbl_TtoC_T.FDR_BH'<FDR_threshold); % number of mutationl signatures fulfulling the JS threshold, as well as the FDR threshold
num_mut_sig_T = length(best_ind_T);
best_similarities_TtoC_T = JS_tbl_TtoC_T(best_ind_T,:); % to find which mutational signature has the best JSD scores for this transition mutation

for i = 1:num_mut_sig_T

    ind_sig_min_T          = js_TtoC_ind_T(i);

    rel_mut_sig_for_plot = mut_sig_for_plot(:,ind_sig_min_T);
    norm_mut_sig         = rel_mut_sig_for_plot{:,:}./sum(rel_mut_sig_for_plot{:,:});
    bar_vals             = [norm_mut_sig,conf_prop_for_plot_T];

    figure('units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'Color','w');

    hb1 = bar(bar_vals);

    hb1(1).FaceColor = color_sig_TtoC;
    hb1(1).EdgeColor = 'k';
    hb1(1).LineWidth = 1;

    hb1(2).FaceColor = 'b';
    hb1(2).EdgeColor = 'k';
    hb1(2).LineWidth = 1;

    axis([0.0 17.0 0.0 0.3]);

    set(gca, 'FontSize', fontsize, 'XminorTick', 'off', 'YminorTick', 'off','LineWidth',linewidth,...
        'XTick',x_ticks,'TickLength',[ticksize,ticksize]);

    xticklabels(seqs_mut_for_plot);

    ylabel('Probability', 'FontSize', 36);

    legend({sprintf('T-->C probability from %s',mut_sig_for_plot(:,ind_sig_min_T).Properties.VariableNames{1,1}),...
        'conformational probability from NMR'});

    text(1,0.25,sprintf('JS divergence = %0.3f', js_TtoC_sorted_T(i)),'FontSize',fontsize);

    propedit;


end



%% only T-->G transition mutations (ind = 6)

ind_TtoG                     = 6;
mut_sig_for_plot             = TtoG(:,2:end);

js_TtoG_T                        = js_div_T(ind_TtoG,:);
[js_TtoG_sorted_T,js_TtoG_ind_T] = sort(js_TtoG_T);

num_sigs = width(mut_sig_for_plot);

JS_pVal_tbl_TtoG_T = cell(num_sigs,9);

for i = 1:num_sigs

    ind_sig_min_T          = js_TtoG_ind_T(i);

    rel_mut_sig_for_plot = mut_sig_for_plot(:,ind_sig_min_T);
    norm_mut_sig         = rel_mut_sig_for_plot{:,:}./sum(rel_mut_sig_for_plot{:,:});

    rel_mut_sig_for_corr_plot = mut_sig_for_plot(ind_seqs_T,ind_sig_min_T);
    norm_mut_sig_for_corr_plot = rel_mut_sig_for_corr_plot{:,:}./sum(rel_mut_sig_for_corr_plot{:,:});

    [r_corr,Rsq,Rsq_det,rms,a_rms_per,bestfit,x2,inBetween,slope,intercept,slope_err,intercept_err,p_corr] = corr_plot(Q_T,norm_mut_sig_for_corr_plot,x_line,conf_int_per);

    p_rand = sum(js_rand<=js_TtoG_sorted_T(i))./num_perms;

    JS_pVal_tbl_TtoG_T{i,1} = mut_sig_for_plot(:,ind_sig_min_T).Properties.VariableNames{1,1};
    JS_pVal_tbl_TtoG_T{i,2} = js_TtoG_sorted_T(i);
    JS_pVal_tbl_TtoG_T{i,3} = r_corr;
    JS_pVal_tbl_TtoG_T{i,4} = p_corr;
    JS_pVal_tbl_TtoG_T{i,5} = rms;
    JS_pVal_tbl_TtoG_T{i,6} = bestfit;
    JS_pVal_tbl_TtoG_T{i,7} = [slope,slope_err];
    JS_pVal_tbl_TtoG_T{i,8} = [intercept,intercept_err];
    JS_pVal_tbl_TtoG_T{i,9} = p_rand;

end

JS_tbl_TtoG_T = cell2table(JS_pVal_tbl_TtoG_T,...
    "VariableNames",["SBS" "JS" "r_pearson" "p_corr" "RMSD" "bestfit" "slope" "intercept" "p_rand"]);
p_values  = JS_tbl_TtoG_T.p_rand;
FDR       = mafdr(p_values,'BHFDR',false);
FDR_BH    = mafdr(p_values,'BHFDR',true);
JS_tbl_TtoG_T = addvars(JS_tbl_TtoG_T,FDR,FDR_BH,'After','p_rand');

best_ind_T = find(js_TtoG_sorted_T<js_threshold&JS_tbl_TtoG_T.FDR_BH'<FDR_threshold); % number of mutationl signatures fulfulling the JS threshold, as well as the FDR threshold
num_mut_sig_T = length(best_ind_T);
best_similarities_TtoG_T = JS_tbl_TtoG_T(best_ind_T,:); % to find which mutational signature has the best JSD scores for this transition mutation

for i = 1:num_mut_sig_T

    ind_sig_min_T          = js_TtoG_ind_T(i);

    rel_mut_sig_for_plot = mut_sig_for_plot(:,ind_sig_min_T);
    norm_mut_sig         = rel_mut_sig_for_plot{:,:}./sum(rel_mut_sig_for_plot{:,:});
    bar_vals             = [norm_mut_sig,conf_prop_for_plot_T];

    figure('units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'Color','w');

    hb1 = bar(bar_vals);

    hb1(1).FaceColor = color_sig_TtoG;
    hb1(1).EdgeColor = 'k';
    hb1(1).LineWidth = 1;

    hb1(2).FaceColor = 'b';
    hb1(2).EdgeColor = 'k';
    hb1(2).LineWidth = 1;

    axis([0.0 17.0 0.0 0.3]);

    set(gca, 'FontSize', fontsize, 'XminorTick', 'off', 'YminorTick', 'off','LineWidth',linewidth,...
        'XTick',x_ticks,'TickLength',[ticksize,ticksize]);

    xticklabels(seqs_mut_for_plot);

    ylabel('Probability', 'FontSize', 36);

    legend({sprintf('T-->G probability from %s',mut_sig_for_plot(:,ind_sig_min_T).Properties.VariableNames{1,1}),...
        'conformational probability from NMR'});

    text(1,0.25,sprintf('JS divergence = %0.3f', js_TtoG_sorted_T(i)),'FontSize',fontsize);

    propedit;

end


%% saving

pathname_output  = sprintf('%s/%s',parent_path,path_output); % file name of SBS mutational signatures (downloaded from the COSMIC data base)
filename_js_rand = sprintf('%s_js_rand.mat',dataset);
filename_js_rand = sprintf('%s/%s',pathname_output,filename_js_rand);

save(filename_js_rand,'js_rand');

JS_tbl_CtoA_save = JS_tbl_CtoA_T(:,[1:3,5,9:11]);
JS_tbl_CtoG_save = JS_tbl_CtoG_T(:,[1:3,5,9:11]);
JS_tbl_CtoT_save = JS_tbl_CtoT_T(:,[1:3,5,9:11]);
JS_tbl_TtoA_save = JS_tbl_TtoA_T(:,[1:3,5,9:11]);
JS_tbl_TtoC_save = JS_tbl_TtoC_T(:,[1:3,5,9:11]);
JS_tbl_TtoG_save = JS_tbl_TtoG_T(:,[1:3,5,9:11]);

JS_best_sim_CtoA_save = best_similarities_CtoA_T(:,[1:3,5,9:11]);
JS_best_sim_CtoG_save = best_similarities_CtoG_T(:,[1:3,5,9:11]);
JS_best_sim_CtoT_save = best_similarities_CtoT_T(:,[1:3,5,9:11]);
JS_best_sim_TtoA_save = best_similarities_TtoA_T(:,[1:3,5,9:11]);
JS_best_sim_TtoC_save = best_similarities_TtoC_T(:,[1:3,5,9:11]);
JS_best_sim_TtoG_save = best_similarities_TtoG_T(:,[1:3,5,9:11]);


writetable(JS_tbl_CtoA_save,sprintf('%s/%s_JSD_CtoA_T_all.csv',pathname_output,dataset));
writetable(JS_tbl_CtoG_save,sprintf('%s/%s_JSD_CtoG_T_all.csv',pathname_output,dataset));
writetable(JS_tbl_CtoT_save,sprintf('%s/%s_JSD_CtoT_T_all.csv',pathname_output,dataset));
writetable(JS_tbl_TtoA_save,sprintf('%s/%s_JSD_TtoA_T_all.csv',pathname_output,dataset));
writetable(JS_tbl_TtoC_save,sprintf('%s/%s_JSD_TtoC_T_all.csv',pathname_output,dataset));
writetable(JS_tbl_TtoG_save,sprintf('%s/%s_JSD_TtoG_T_all.csv',pathname_output,dataset));

writetable(JS_best_sim_CtoA_save,sprintf('%s/%s_JSD_CtoA_T_top_hits.csv',pathname_output,dataset));
writetable(JS_best_sim_CtoG_save,sprintf('%s/%s_JSD_CtoG_T_top_hits.csv',pathname_output,dataset));
writetable(JS_best_sim_CtoT_save,sprintf('%s/%s_JSD_CtoT_T_top_hits.csv',pathname_output,dataset));
writetable(JS_best_sim_TtoA_save,sprintf('%s/%s_JSD_TtoA_T_top_hits.csv',pathname_output,dataset));
writetable(JS_best_sim_TtoC_save,sprintf('%s/%s_JSD_TtoC_T_top_hits.csv',pathname_output,dataset));
writetable(JS_best_sim_TtoG_save,sprintf('%s/%s_JSD_TtoG_T_top_hits.csv',pathname_output,dataset));




