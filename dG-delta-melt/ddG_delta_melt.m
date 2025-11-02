clear;
clc;

%% comparison to dG delta-melt

parent_path_dG_melt = '/Users/orsula1/OrsDocs/OrsDocs_since_Jan2025/2025-Sequence-dependence/2025-Data-Deposition/dG-delta-melt';
filename_dG_melt    = 'ddGmelt-GC-vs-GT.xlsx';
filepath_dG_melt    = sprintf('%s/%s',parent_path_dG_melt,filename_dG_melt);

%% dG melt, GT

params_dG_melt_GT       = readtable(filepath_dG_melt,'Sheet','GT');
dG_melt_GT_seq_T_all    = params_dG_melt_GT.Seq_T;

dG_melt_GT_25C_all      = - params_dG_melt_GT.dG25; % values in excel file are negative, here we take the positive values, to match with the deltaMelt paper, defined as "energy required for the melting"
dG_melt_GT_25C_all_errs = params_dG_melt_GT.dG25_std;

dG_melt_GT_25C          = dG_melt_GT_25C_all(~isnan(dG_melt_GT_25C_all));
dG_melt_GT_25C_errs     = dG_melt_GT_25C_all_errs(~isnan(dG_melt_GT_25C_all));
dG_melt_GT_seq_T        = dG_melt_GT_seq_T_all(~isnan(dG_melt_GT_25C_all));

%% dG melt, GC

params_dG_melt_GC       = readtable(filepath_dG_melt,'Sheet','GC');
dG_melt_GC_seq_T_all    = params_dG_melt_GC.Seq_T;

dG_melt_GC_25C_all       = - params_dG_melt_GC.dG25; % values in excel file are negative, here we take the positive values, to match with the deltaMelt paper, defined as "energy required for the melting"
dG_melt_GC_25C_all_errs  = params_dG_melt_GC.dG25_std;

dG_melt_GC_25C          = dG_melt_GC_25C_all(~isnan(dG_melt_GC_25C_all));
dG_melt_GC_25C_errs     = dG_melt_GC_25C_all_errs(~isnan(dG_melt_GC_25C_all));
dG_melt_GC_seq_T        = dG_melt_GC_seq_T_all(~isnan(dG_melt_GC_25C_all));


%% ddG melt, GT-GC

[~,ind_dG_melt] = ismember(dG_melt_GC_seq_T,dG_melt_GT_seq_T);

ddG_melt_GT_GC_25C      = dG_melt_GT_25C - dG_melt_GC_25C(ind_dG_melt);
ddG_melt_GT_GC_25C_errs = sqrt(dG_melt_GT_25C_errs.^2 + dG_melt_GC_25C_errs(ind_dG_melt).^2);

%% saving dG melt data

path_to_save_dG = '/Users/orsula1/OrsDocs/OrsDocs_since_Jan2025/2025-Sequence-dependence/2025-Data-Deposition/dG-delta-melt';
filename_dG     = sprintf('%s/dG_GT_GC_delta_melt.mat',path_to_save_dG);

save(filename_dG,'dG_melt_GT_seq_T','dG_melt_GT_25C','dG_melt_GT_25C_errs',...
    'dG_melt_GC_seq_T','dG_melt_GC_25C','dG_melt_GC_25C_errs',...
    'ddG_melt_GT_GC_25C','ddG_melt_GT_GC_25C_errs');