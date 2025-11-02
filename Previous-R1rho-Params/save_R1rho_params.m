clear;
clc;

%% loading parameters, previous data

parent_path = '/Users/orsula1/OrsDocs/OrsDocs_since_Jan2025/2025-Sequence-dependence/2025-Data-Deposition/Previous-R1rho-Params';
filename = '2023-GT-anion-seq-dep-incl-Genol-Tenol.xlsx';
filepath_RD = sprintf('%s/%s',parent_path,filename);

%% get previous data pKas (Issac, Nature 2018), and sort them

params_High_pH_10C_tbl     = readtable(filepath_RD,'Sheet','10C_all_pHs');

construct_name_High_pH_10C = params_High_pH_10C_tbl.Construct;
orig_pH_High_pH_10C        = params_High_pH_10C_tbl.orig_pH;

params_High_pH_25C_tbl     = readtable(filepath_RD,'Sheet','25C_all_pHs');

construct_name_High_pH_25C = params_High_pH_25C_tbl.Construct;
orig_pH_High_pH_25C        = params_High_pH_25C_tbl.orig_pH;

sort_ind_10C = [5,6:7,9,1:4,12,10:11,14:16]; % without TGA (outlier) % without RRE (outlier)
sort_ind_25C = [4:5,1:3,9,8,6:7,10:13]; 

params_High_pH_10C_sort_tbl     = params_High_pH_10C_tbl(sort_ind_10C,:);
construct_name_High_pH_10C_sort = construct_name_High_pH_10C(sort_ind_10C,:);
orig_pH_High_pH_10C_sort        = orig_pH_High_pH_10C(sort_ind_10C,:);

params_25C_High_pH_sort_tbl     = params_High_pH_25C_tbl(sort_ind_25C,:);
construct_name_High_pH_25C_sort = construct_name_High_pH_25C(sort_ind_25C,:);
orig_pH_High_pH_25C_sort        = orig_pH_High_pH_25C(sort_ind_25C,:);

ind_DNA_10C    = find(contains(params_High_pH_10C_sort_tbl.NA_id,'DNA')); 
ind_RNA_10C    = find(contains(params_High_pH_10C_sort_tbl.NA_id,'RNA')); 
ind_Mod_10C    = find(contains(params_High_pH_10C_sort_tbl.NA_id,'mod'));

ind_DNA_25C    = find(contains(params_25C_High_pH_sort_tbl.NA_id,'DNA'));
ind_RNA_25C    = find(contains(params_25C_High_pH_sort_tbl.NA_id,'RNA'));
ind_Hybrid_25C = [10:11,13];

pKas_prev_High_pH_10C        = params_High_pH_10C_sort_tbl.orig_pH - log10(params_High_pH_10C_sort_tbl.pES2_Anion./(1-params_High_pH_10C_sort_tbl.pES2_Anion));
pKas_prev_High_pH_10C_errs   = (1/log(10))*sqrt(((params_High_pH_10C_sort_tbl.pES2_Anion_Error)./(params_High_pH_10C_sort_tbl.pES2_Anion)).^2+(params_High_pH_10C_sort_tbl.pES2_Anion_Error).^2./((1-params_High_pH_10C_sort_tbl.pES2_Anion).^2));

pKas_prev_High_pH_25C        = params_25C_High_pH_sort_tbl.orig_pH - log10(params_25C_High_pH_sort_tbl.pES2_Anion./(1-params_25C_High_pH_sort_tbl.pES2_Anion));
pKas_prev_High_pH_25C_errs   = (1/log(10))*sqrt(((params_25C_High_pH_sort_tbl.pES2_Anion_Error)./(params_25C_High_pH_sort_tbl.pES2_Anion)).^2+(params_25C_High_pH_sort_tbl.pES2_Anion_Error).^2./((1-params_25C_High_pH_sort_tbl.pES2_Anion).^2));

%% Previous data pKas (sorted) - with sequence name

seq_prev_sort_High_pH_10C = cell(length(construct_name_High_pH_10C_sort),1);

for i = 1:length(construct_name_High_pH_10C_sort)

    seq_prev_sort_High_pH_10C{i,1} = construct_name_High_pH_10C_sort{i,1}(end-2:end);

end

seq_prev_sort_High_pH_25C = cell(length(construct_name_High_pH_25C_sort),1);

for i = 1:length(construct_name_High_pH_25C_sort)

    seq_prev_sort_High_pH_25C{i,1} = construct_name_High_pH_25C_sort{i,1}(end-2:end);

end


%% pAnion and pTautomer - high pH Issac

seq_prev_DNA_High_pH_10C         = seq_prev_sort_High_pH_10C(ind_DNA_10C);
seq_prev_DNA_High_pH_10C_T       = convSeqtoT(seq_prev_DNA_High_pH_10C);
pAnion_prev_DNA_High_pH_10C      = 100.*params_High_pH_10C_sort_tbl.pES2_Anion(ind_DNA_10C);
pAnion_prev_DNA_High_pH_10C_errs = 100.*params_High_pH_10C_sort_tbl.pES2_Anion_Error(ind_DNA_10C);
pTaut_prev_DNA_High_pH_10C       = 100.*params_High_pH_10C_sort_tbl.pES1_Tautomer(ind_DNA_10C);
pTaut_prev_DNA_High_pH_10C_errs  = 100.*params_High_pH_10C_sort_tbl.pES1_Tautomer_Error(ind_DNA_10C);
pGenol_prev_DNA_High_pH_10C      = pTaut_prev_DNA_High_pH_10C.*params_High_pH_10C_sort_tbl.Genol(ind_DNA_10C);
pGenol_prev_DNA_High_pH_10C_errs = pGenol_prev_DNA_High_pH_10C.*sqrt((pTaut_prev_DNA_High_pH_10C_errs./pTaut_prev_DNA_High_pH_10C).^2 + (params_High_pH_10C_sort_tbl.Genol_err(ind_DNA_10C)./params_High_pH_10C_sort_tbl.Genol(ind_DNA_10C)).^2);
pTenol_prev_DNA_High_pH_10C      = pTaut_prev_DNA_High_pH_10C.*params_High_pH_10C_sort_tbl.Tenol(ind_DNA_10C);
pTenol_prev_DNA_High_pH_10C_errs = pTenol_prev_DNA_High_pH_10C.*sqrt((pTaut_prev_DNA_High_pH_10C_errs./pTaut_prev_DNA_High_pH_10C).^2 + (params_High_pH_10C_sort_tbl.Tenol_err(ind_DNA_10C)./params_High_pH_10C_sort_tbl.Tenol(ind_DNA_10C)).^2);


seq_prev_RNA_High_pH_10C         = {'CGG';'CGC';'CGC'};
seq_prev_RNA_High_pH_10C_T       = convSeqtoT(seq_prev_RNA_High_pH_10C);
pAnion_prev_RNA_High_pH_10C      = 100.*params_High_pH_10C_sort_tbl.pES2_Anion(ind_RNA_10C);
pAnion_prev_RNA_High_pH_10C_errs = 100.*params_High_pH_10C_sort_tbl.pES2_Anion_Error(ind_RNA_10C);
pTaut_prev_RNA_High_pH_10C       = 100.*params_High_pH_10C_sort_tbl.pES1_Tautomer(ind_RNA_10C);
pTaut_prev_RNA_High_pH_10C_errs  = 100.*params_High_pH_10C_sort_tbl.pES1_Tautomer_Error(ind_RNA_10C);
pGenol_prev_RNA_High_pH_10C      = pTaut_prev_RNA_High_pH_10C.*params_High_pH_10C_sort_tbl.Genol(ind_RNA_10C);
pGenol_prev_RNA_High_pH_10C_errs = pGenol_prev_RNA_High_pH_10C.*sqrt((pTaut_prev_RNA_High_pH_10C_errs./pTaut_prev_RNA_High_pH_10C).^2 + (params_High_pH_10C_sort_tbl.Genol_err(ind_RNA_10C)./params_High_pH_10C_sort_tbl.Genol(ind_RNA_10C)).^2);
pTenol_prev_RNA_High_pH_10C      = pTaut_prev_RNA_High_pH_10C.*params_High_pH_10C_sort_tbl.Tenol(ind_RNA_10C);
pTenol_prev_RNA_High_pH_10C_errs = pTenol_prev_RNA_High_pH_10C.*sqrt((pTaut_prev_RNA_High_pH_10C_errs./pTaut_prev_RNA_High_pH_10C).^2 + (params_High_pH_10C_sort_tbl.Tenol_err(ind_RNA_10C)./params_High_pH_10C_sort_tbl.Tenol(ind_RNA_10C)).^2);

seq_prev_DNA_High_pH_25C         = seq_prev_sort_High_pH_25C(ind_DNA_25C);
seq_prev_DNA_High_pH_25C_T       = convSeqtoT(seq_prev_DNA_High_pH_25C);
pAnion_prev_DNA_High_pH_25C      = 100.*params_25C_High_pH_sort_tbl.pES2_Anion(ind_DNA_25C);
pAnion_prev_DNA_High_pH_25C_errs = 100.*params_25C_High_pH_sort_tbl.pES2_Anion_Error(ind_DNA_25C);
pTaut_prev_DNA_High_pH_25C       = 100.*params_25C_High_pH_sort_tbl.pES1_Tautomer(ind_DNA_25C);
pTaut_prev_DNA_High_pH_25C_errs  = 100.*params_25C_High_pH_sort_tbl.pES1_Tautomer_Error(ind_DNA_25C);
pGenol_prev_DNA_High_pH_25C      = pTaut_prev_DNA_High_pH_25C.*params_25C_High_pH_sort_tbl.Genol(ind_DNA_25C);
pGenol_prev_DNA_High_pH_25C_errs = pGenol_prev_DNA_High_pH_25C.*sqrt((pTaut_prev_DNA_High_pH_25C_errs./pTaut_prev_DNA_High_pH_25C).^2 + (params_25C_High_pH_sort_tbl.Genol_err(ind_DNA_25C)./params_25C_High_pH_sort_tbl.Genol(ind_DNA_25C)).^2);
pTenol_prev_DNA_High_pH_25C      = pTaut_prev_DNA_High_pH_25C.*params_25C_High_pH_sort_tbl.Tenol(ind_DNA_25C);
pTenol_prev_DNA_High_pH_25C_errs = pTenol_prev_DNA_High_pH_25C.*sqrt((pTaut_prev_DNA_High_pH_25C_errs./pTaut_prev_DNA_High_pH_25C).^2 + (params_25C_High_pH_sort_tbl.Tenol_err(ind_DNA_25C)./params_25C_High_pH_sort_tbl.Tenol(ind_DNA_25C)).^2);

seq_prev_RNA_High_pH_25C         = {'CGG';'TGC';'CGC';'CGC'};
seq_prev_RNA_High_pH_25C_T       = convSeqtoT(seq_prev_RNA_High_pH_25C);
pAnion_prev_RNA_High_pH_25C      = 100.*params_25C_High_pH_sort_tbl.pES2_Anion(ind_RNA_25C);
pAnion_prev_RNA_High_pH_25C_errs = 100.*params_25C_High_pH_sort_tbl.pES2_Anion_Error(ind_RNA_25C);
pTaut_prev_RNA_High_pH_25C       = 100.*params_25C_High_pH_sort_tbl.pES1_Tautomer(ind_RNA_25C);
pTaut_prev_RNA_High_pH_25C_errs  = 100.*params_25C_High_pH_sort_tbl.pES1_Tautomer_Error(ind_RNA_25C);
pGenol_prev_RNA_High_pH_25C      = pTaut_prev_RNA_High_pH_25C.*params_25C_High_pH_sort_tbl.Genol(ind_RNA_25C);
pGenol_prev_RNA_High_pH_25C_errs = pGenol_prev_RNA_High_pH_25C.*sqrt((pTaut_prev_RNA_High_pH_25C_errs./pTaut_prev_RNA_High_pH_25C).^2 + (params_25C_High_pH_sort_tbl.Genol_err(ind_RNA_25C)./params_25C_High_pH_sort_tbl.Genol(ind_RNA_25C)).^2);
pTenol_prev_RNA_High_pH_25C      = pTaut_prev_RNA_High_pH_25C.*params_25C_High_pH_sort_tbl.Tenol(ind_RNA_25C);
pTenol_prev_RNA_High_pH_25C_errs = pTenol_prev_RNA_High_pH_25C.*sqrt((pTaut_prev_RNA_High_pH_25C_errs./pTaut_prev_RNA_High_pH_25C).^2 + (params_25C_High_pH_sort_tbl.Tenol_err(ind_RNA_25C)./params_25C_High_pH_sort_tbl.Tenol(ind_RNA_25C)).^2);


seq_prev_Hybrid_High_pH_25C         = {'GGA';'GGA';'TGC'};
seq_prev_Hybrid_High_pH_25C_T       = convSeqtoT(seq_prev_Hybrid_High_pH_25C);
pAnion_prev_Hybrid_High_pH_25C      = 100.*params_25C_High_pH_sort_tbl.pES2_Anion(ind_Hybrid_25C);
pAnion_prev_Hybrid_High_pH_25C_errs = 100.*params_25C_High_pH_sort_tbl.pES2_Anion_Error(ind_Hybrid_25C);
pTaut_prev_Hybrid_High_pH_25C       = 100.*params_25C_High_pH_sort_tbl.pES1_Tautomer(ind_Hybrid_25C);
pTaut_prev_Hybrid_High_pH_25C_errs  = 100.*params_25C_High_pH_sort_tbl.pES1_Tautomer_Error(ind_Hybrid_25C);
pGenol_prev_Hybrid_High_pH_25C      = pTaut_prev_Hybrid_High_pH_25C.*params_25C_High_pH_sort_tbl.Genol(ind_Hybrid_25C);
pGenol_prev_Hybrid_High_pH_25C_errs = pGenol_prev_Hybrid_High_pH_25C.*sqrt((pTaut_prev_Hybrid_High_pH_25C_errs./pTaut_prev_Hybrid_High_pH_25C).^2 + (params_25C_High_pH_sort_tbl.Genol_err(ind_Hybrid_25C)./params_25C_High_pH_sort_tbl.Genol(ind_Hybrid_25C)).^2);
pTenol_prev_Hybrid_High_pH_25C      = pTaut_prev_Hybrid_High_pH_25C.*params_25C_High_pH_sort_tbl.Tenol(ind_Hybrid_25C);
pTenol_prev_Hybrid_High_pH_25C_errs = pTenol_prev_Hybrid_High_pH_25C.*sqrt((pTaut_prev_Hybrid_High_pH_25C_errs./pTaut_prev_Hybrid_High_pH_25C).^2 + (params_25C_High_pH_sort_tbl.Tenol_err(ind_Hybrid_25C)./params_25C_High_pH_sort_tbl.Tenol(ind_Hybrid_25C)).^2);


%% pAnion extrapolated to pH 6.9 vs. pTaut measured at pH 6.9

params_6p9_10C_tbl     = readtable(filepath_RD,'Sheet','10C_6p9');

construct_name_6p9     = params_6p9_10C_tbl.Construct;
orig_pH_6p9_10C        = params_6p9_10C_tbl.orig_pH;

params_6p9_25C_tbl     = readtable(filepath_RD,'Sheet','25C_6p9');

construct_name_6p9_25C = params_6p9_25C_tbl.Construct;
orig_pH_6p9_25C        = params_6p9_25C_tbl.orig_pH;

sort_ind_6p9_10C = [2,3,4,5,1,7,6];
sort_ind_6p9_25C = [2,3,4,5,1,9,12,11,7,6,10,8,13:17]; % RNA sorted, with RRE

params_6p9_10C_sort_tbl     = params_6p9_10C_tbl(sort_ind_6p9_10C,:);
construct_name_6p9_10C_sort = construct_name_6p9(sort_ind_6p9_10C,:);
orig_pH_6p9_10C_sort        = orig_pH_6p9_10C(sort_ind_6p9_10C,:);

params_6p9_25C_sort_tbl     = params_6p9_25C_tbl(sort_ind_6p9_25C,:);
construct_name_6p9_25C_sort = construct_name_6p9_25C(sort_ind_6p9_25C,:);
orig_pH_6p9_25C_sort        = orig_pH_6p9_25C(sort_ind_6p9_25C,:);

ind_DNA_6p9_10C    = [1:3,5]; % without Dickerson-CGA outlier (TTG)
ind_RNA_6p9_10C    = find(contains(params_6p9_10C_sort_tbl.NA_id,'RNA'));
ind_Mod_6p9_10C    = find(contains(params_6p9_10C_sort_tbl.NA_id,'mod'));

ind_DNA_6p9_25C    = find(contains(params_6p9_25C_sort_tbl.NA_id,'DNA'));
ind_RNA_6p9_25C    = find(contains(params_6p9_25C_sort_tbl.NA_id,'RNA'));
ind_Hybrid_6p9_25C = [11,13];

seq_sort_6p9_10C = cell(length(construct_name_6p9_10C_sort),1);

for i = 1:length(construct_name_6p9_10C_sort)

    seq_sort_6p9_10C{i,1} = construct_name_6p9_10C_sort{i,1}(end-2:end);

end

seq_sort_6p9_25C = cell(length(construct_name_6p9_25C_sort),1);

for i = 1:length(construct_name_6p9_25C_sort)

    seq_sort_6p9_25C{i,1} = construct_name_6p9_25C_sort{i,1}(end-2:end);

end

seq_DNA_6p9_10C     = seq_sort_6p9_10C(ind_DNA_6p9_10C);
seq_DNA_6p9_10C_T   = convSeqtoT(seq_DNA_6p9_10C);

pTaut_DNA_6p9_10C       = 100.*params_6p9_10C_sort_tbl.pES1_Tautomer(ind_DNA_6p9_10C);
pTaut_DNA_6p9_10C_errs  = 100.*params_6p9_10C_sort_tbl.pES1_Tautomer_Error(ind_DNA_6p9_10C);
pGenol_DNA_6p9_10C      = pTaut_DNA_6p9_10C.*params_6p9_10C_sort_tbl.Genol(ind_DNA_6p9_10C);
pGenol_DNA_6p9_10C_errs = pGenol_DNA_6p9_10C.*sqrt((pTaut_DNA_6p9_10C_errs./pTaut_DNA_6p9_10C).^2 + (params_6p9_10C_sort_tbl.Genol_err(ind_DNA_6p9_10C)./params_6p9_10C_sort_tbl.Genol(ind_DNA_6p9_10C)).^2);
pTenol_DNA_6p9_10C      = pTaut_DNA_6p9_10C.*params_6p9_10C_sort_tbl.Tenol(ind_DNA_6p9_10C);
pTenol_DNA_6p9_10C_errs = pTenol_DNA_6p9_10C.*sqrt((pTaut_DNA_6p9_10C_errs./pTaut_DNA_6p9_10C).^2 + (params_6p9_10C_sort_tbl.Tenol_err(ind_DNA_6p9_10C)./params_6p9_10C_sort_tbl.Tenol(ind_DNA_6p9_10C)).^2);

seq_RNA_6p9_10C     = seq_sort_6p9_10C(ind_RNA_6p9_10C);
seq_RNA_6p9_10C_T   = convSeqtoT(seq_RNA_6p9_10C);

pTaut_RNA_6p9_10C       = 100.*params_6p9_10C_sort_tbl.pES1_Tautomer(ind_RNA_6p9_10C);
pTaut_RNA_6p9_10C_errs  = 100.*params_6p9_10C_sort_tbl.pES1_Tautomer_Error(ind_RNA_6p9_10C);
pGenol_RNA_6p9_10C      = pTaut_RNA_6p9_10C.*params_6p9_10C_sort_tbl.Genol(ind_RNA_6p9_10C);
pGenol_RNA_6p9_10C_errs = pGenol_RNA_6p9_10C.*sqrt((pTaut_RNA_6p9_10C_errs./pTaut_RNA_6p9_10C).^2 + (params_6p9_10C_sort_tbl.Genol_err(ind_RNA_6p9_10C)./params_6p9_10C_sort_tbl.Genol(ind_RNA_6p9_10C)).^2);
pTenol_RNA_6p9_10C      = pTaut_RNA_6p9_10C.*params_6p9_10C_sort_tbl.Tenol(ind_RNA_6p9_10C);
pTenol_RNA_6p9_10C_errs = pTenol_RNA_6p9_10C.*sqrt((pTaut_RNA_6p9_10C_errs./pTaut_RNA_6p9_10C).^2 + (params_6p9_10C_sort_tbl.Tenol_err(ind_RNA_6p9_10C)./params_6p9_10C_sort_tbl.Tenol(ind_RNA_6p9_10C)).^2);

seq_DNA_6p9_25C     = seq_sort_6p9_25C(ind_DNA_6p9_25C);
seq_DNA_6p9_25C_T   = convSeqtoT(seq_DNA_6p9_25C);

pTaut_DNA_6p9_25C       = 100.*params_6p9_25C_sort_tbl.pES1_Tautomer(ind_DNA_6p9_25C);
pTaut_DNA_6p9_25C_errs  = 100.*params_6p9_25C_sort_tbl.pES1_Tautomer_Error(ind_DNA_6p9_25C);
pGenol_DNA_6p9_25C      = pTaut_DNA_6p9_25C.*params_6p9_25C_sort_tbl.Genol(ind_DNA_6p9_25C);
pGenol_DNA_6p9_25C_errs = pGenol_DNA_6p9_25C.*sqrt((pTaut_DNA_6p9_25C_errs./pTaut_DNA_6p9_25C).^2 + (params_6p9_25C_sort_tbl.Genol_err(ind_DNA_6p9_25C)./params_6p9_25C_sort_tbl.Genol(ind_DNA_6p9_25C)).^2);
pTenol_DNA_6p9_25C      = pTaut_DNA_6p9_25C.*params_6p9_25C_sort_tbl.Tenol(ind_DNA_6p9_25C);
pTenol_DNA_6p9_25C_errs = pTenol_DNA_6p9_25C.*sqrt((pTaut_DNA_6p9_25C_errs./pTaut_DNA_6p9_25C).^2 + (params_6p9_25C_sort_tbl.Tenol_err(ind_DNA_6p9_25C)./params_6p9_25C_sort_tbl.Tenol(ind_DNA_6p9_25C)).^2);

seq_RNA_6p9_25C     = seq_sort_6p9_25C(ind_RNA_6p9_25C);
seq_RNA_6p9_25C_U   = convSeqtoU(seq_RNA_6p9_25C);
seq_RNA_6p9_25C_T   = cell(length(seq_RNA_6p9_25C_U),1);

for i = 1:length(seq_RNA_6p9_25C_U)

    temp = seq_RNA_6p9_25C_U{i,1};
    seq_RNA_6p9_25C_T{i,1} = [temp(1),'T',temp(end)];

end

pTaut_RNA_6p9_25C       = 100.*params_6p9_25C_sort_tbl.pES1_Tautomer(ind_RNA_6p9_25C);
pTaut_RNA_6p9_25C_errs  = 100.*params_6p9_25C_sort_tbl.pES1_Tautomer_Error(ind_RNA_6p9_25C);
pGenol_RNA_6p9_25C      = pTaut_RNA_6p9_25C.*params_6p9_25C_sort_tbl.Genol(ind_RNA_6p9_25C);
pGenol_RNA_6p9_25C_errs = pGenol_RNA_6p9_25C.*sqrt((pTaut_RNA_6p9_25C_errs./pTaut_RNA_6p9_25C).^2 + (params_6p9_25C_sort_tbl.Genol_err(ind_RNA_6p9_25C)./params_6p9_25C_sort_tbl.Genol(ind_RNA_6p9_25C)).^2);
pTenol_RNA_6p9_25C      = pTaut_RNA_6p9_25C.*params_6p9_25C_sort_tbl.Tenol(ind_RNA_6p9_25C);
pTenol_RNA_6p9_25C_errs = pTenol_RNA_6p9_25C.*sqrt((pTaut_RNA_6p9_25C_errs./pTaut_RNA_6p9_25C).^2 + (params_6p9_25C_sort_tbl.Tenol_err(ind_RNA_6p9_25C)./params_6p9_25C_sort_tbl.Tenol(ind_RNA_6p9_25C)).^2);

% pTaut is pH-independent, so the RNA:DNA hybrid tatutomer populations at pH 6.9 are taken from the high pH data
seq_prev_Hybrid_6p9_25C         = {'GGA';'TGC'};
seq_prev_Hybrid_6p9_25C_T       = convSeqtoT(seq_prev_Hybrid_6p9_25C);
pTaut_prev_Hybrid_6p9_25C       = 100.*params_25C_High_pH_sort_tbl.pES1_Tautomer(ind_Hybrid_6p9_25C);
pTaut_prev_Hybrid_6p9_25C_errs  = 100.*params_25C_High_pH_sort_tbl.pES1_Tautomer_Error(ind_Hybrid_6p9_25C);
pGenol_prev_Hybrid_6p9_25C      = pTaut_prev_Hybrid_6p9_25C.*params_25C_High_pH_sort_tbl.Genol(ind_Hybrid_6p9_25C);
pGenol_prev_Hybrid_6p9_25C_errs = pGenol_prev_Hybrid_6p9_25C.*sqrt((pTaut_prev_Hybrid_6p9_25C_errs./pTaut_prev_Hybrid_6p9_25C).^2 + (params_25C_High_pH_sort_tbl.Genol_err(ind_Hybrid_6p9_25C)./params_25C_High_pH_sort_tbl.Genol(ind_Hybrid_6p9_25C)).^2);
pTenol_prev_Hybrid_6p9_25C      = pTaut_prev_Hybrid_6p9_25C.*params_25C_High_pH_sort_tbl.Tenol(ind_Hybrid_6p9_25C);
pTenol_prev_Hybrid_6p9_25C_errs = pTenol_prev_Hybrid_6p9_25C.*sqrt((pTaut_prev_Hybrid_6p9_25C_errs./pTaut_prev_Hybrid_6p9_25C).^2 + (params_25C_High_pH_sort_tbl.Tenol_err(ind_Hybrid_6p9_25C)./params_25C_High_pH_sort_tbl.Tenol(ind_Hybrid_6p9_25C)).^2);


%% saving tautomer populations

Taut_prev_DNA_High_pH_10C = table(seq_prev_DNA_High_pH_10C_T,params_High_pH_10C_sort_tbl.orig_pH(ind_DNA_10C),...
    pTaut_prev_DNA_High_pH_10C,pTaut_prev_DNA_High_pH_10C_errs,...
    params_High_pH_10C_sort_tbl.kGStoES1(ind_DNA_10C),params_High_pH_10C_sort_tbl.kES1toGS(ind_DNA_10C));
Taut_prev_DNA_High_pH_10C = renamevars(Taut_prev_DNA_High_pH_10C,["Var2","Var5","Var6"],["orig_pH","kGStoES1","kES1toGS"]);

Taut_prev_RNA_High_pH_10C = table(seq_prev_RNA_High_pH_10C_T,params_High_pH_10C_sort_tbl.orig_pH(ind_RNA_10C),...
    pTaut_prev_RNA_High_pH_10C,pTaut_prev_RNA_High_pH_10C_errs,...
    params_High_pH_10C_sort_tbl.kGStoES1(ind_RNA_10C),params_High_pH_10C_sort_tbl.kES1toGS(ind_RNA_10C));
Taut_prev_RNA_High_pH_10C = renamevars(Taut_prev_RNA_High_pH_10C,["Var2","Var5","Var6"],["orig_pH","kGStoES1","kES1toGS"]);

Taut_prev_DNA_High_pH_25C = table(seq_prev_DNA_High_pH_25C_T,params_25C_High_pH_sort_tbl.orig_pH(ind_DNA_25C),...
    pTaut_prev_DNA_High_pH_25C,pTaut_prev_DNA_High_pH_25C_errs,...
    params_25C_High_pH_sort_tbl.kGStoES1(ind_DNA_25C),params_25C_High_pH_sort_tbl.kES1toGS(ind_DNA_25C));
Taut_prev_DNA_High_pH_25C = renamevars(Taut_prev_DNA_High_pH_25C,["Var2","Var5","Var6"],["orig_pH","kGStoES1","kES1toGS"]);

Taut_prev_RNA_High_pH_25C = table(seq_prev_RNA_High_pH_25C_T,params_25C_High_pH_sort_tbl.orig_pH(ind_RNA_25C),...
    pTaut_prev_RNA_High_pH_25C,pTaut_prev_RNA_High_pH_25C_errs,...
    params_25C_High_pH_sort_tbl.kGStoES1(ind_RNA_25C),params_25C_High_pH_sort_tbl.kES1toGS(ind_RNA_25C));
Taut_prev_RNA_High_pH_25C = renamevars(Taut_prev_RNA_High_pH_25C,["Var2","Var5","Var6"],["orig_pH","kGStoES1","kES1toGS"]);

Taut_prev_Hybrid_High_pH_25C = table(seq_prev_Hybrid_High_pH_25C_T,params_25C_High_pH_sort_tbl.orig_pH(ind_Hybrid_25C),...
    pTaut_prev_Hybrid_High_pH_25C,pTaut_prev_Hybrid_High_pH_25C_errs,...
    params_25C_High_pH_sort_tbl.kGStoES1(ind_Hybrid_25C),params_25C_High_pH_sort_tbl.kES1toGS(ind_Hybrid_25C));
Taut_prev_Hybrid_High_pH_25C = renamevars(Taut_prev_Hybrid_High_pH_25C,["Var2","Var5","Var6"],["orig_pH","kGStoES1","kES1toGS"]);

Taut_DNA_6p9_10C = table(seq_DNA_6p9_10C_T,params_6p9_10C_sort_tbl.orig_pH(ind_DNA_6p9_10C),...
    pTaut_DNA_6p9_10C,pTaut_DNA_6p9_10C_errs,...
    params_6p9_10C_sort_tbl.kGStoES1(ind_DNA_6p9_10C),params_6p9_10C_sort_tbl.kES1toGS(ind_DNA_6p9_10C));
Taut_DNA_6p9_10C = renamevars(Taut_DNA_6p9_10C,["Var2","Var5","Var6"],["orig_pH","kGStoES1","kES1toGS"]);

Taut_RNA_6p9_10C = table(seq_RNA_6p9_10C_T,params_6p9_10C_sort_tbl.orig_pH(ind_RNA_6p9_10C),...
    pTaut_RNA_6p9_10C,pTaut_RNA_6p9_10C_errs,...
    params_6p9_10C_sort_tbl.kGStoES1(ind_RNA_6p9_10C),params_6p9_10C_sort_tbl.kES1toGS(ind_RNA_6p9_10C));
Taut_RNA_6p9_10C = renamevars(Taut_RNA_6p9_10C,["Var2","Var5","Var6"],["orig_pH","kGStoES1","kES1toGS"]);

Taut_DNA_6p9_25C = table(seq_DNA_6p9_25C_T,params_6p9_25C_sort_tbl.orig_pH(ind_DNA_6p9_25C),...
    pTaut_DNA_6p9_25C,pTaut_DNA_6p9_25C_errs,...
    params_6p9_25C_sort_tbl.kGStoES1(ind_DNA_6p9_25C),params_6p9_25C_sort_tbl.kES1toGS(ind_DNA_6p9_25C));
Taut_DNA_6p9_25C = renamevars(Taut_DNA_6p9_25C,["Var2","Var5","Var6"],["orig_pH","kGStoES1","kES1toGS"]);

Taut_RNA_6p9_25C = table(seq_RNA_6p9_25C_T,params_6p9_25C_sort_tbl.orig_pH(ind_RNA_6p9_25C),...
    pTaut_RNA_6p9_25C,pTaut_RNA_6p9_25C_errs,...
    params_6p9_25C_sort_tbl.kGStoES1(ind_RNA_6p9_25C),params_6p9_25C_sort_tbl.kES1toGS(ind_RNA_6p9_25C));
Taut_RNA_6p9_25C = renamevars(Taut_RNA_6p9_25C,["Var2","Var5","Var6"],["orig_pH","kGStoES1","kES1toGS"]);



Taut_pop_all_10C = [Taut_prev_DNA_High_pH_10C.pTaut_prev_DNA_High_pH_10C;Taut_prev_RNA_High_pH_10C.pTaut_prev_RNA_High_pH_10C;...
    Taut_DNA_6p9_10C.pTaut_DNA_6p9_10C;Taut_RNA_6p9_10C.pTaut_RNA_6p9_10C];
Taut_pop_avg_10C = mean(Taut_pop_all_10C);
Taut_pop_std_10C = std(Taut_pop_all_10C);


Taut_k12_all_10C = [Taut_prev_DNA_High_pH_10C.kGStoES1;Taut_prev_RNA_High_pH_10C.kGStoES1;...
    Taut_DNA_6p9_10C.kGStoES1;Taut_RNA_6p9_10C.kGStoES1];
Taut_k12_avg_10C = mean(Taut_k12_all_10C);
Taut_k12_std_10C = std(Taut_k12_all_10C);


Taut_k21_all_10C = [Taut_prev_DNA_High_pH_10C.kES1toGS;Taut_prev_RNA_High_pH_10C.kES1toGS;...
    Taut_DNA_6p9_10C.kES1toGS;Taut_RNA_6p9_10C.kES1toGS];
Taut_k21_avg_10C = mean(Taut_k21_all_10C);
Taut_k21_std_10C = std(Taut_k21_all_10C);
Taut_k21_se_10C  = std(Taut_k21_all_10C)./sqrt(length(Taut_k21_all_10C));
Taut_k21_med_10C = median(Taut_k21_all_10C);
Taut_k21_med_dev_10C = mad(Taut_k21_all_10C,1);



Taut_pop_all_25C = [Taut_prev_DNA_High_pH_25C.pTaut_prev_DNA_High_pH_25C;Taut_prev_RNA_High_pH_25C.pTaut_prev_RNA_High_pH_25C;...
    Taut_prev_Hybrid_High_pH_25C.pTaut_prev_Hybrid_High_pH_25C;...
    Taut_DNA_6p9_25C.pTaut_DNA_6p9_25C;Taut_RNA_6p9_25C.pTaut_RNA_6p9_25C];
Taut_pop_avg_25C = mean(Taut_pop_all_25C);
Taut_pop_std_25C = std(Taut_pop_all_25C);

Taut_k12_all_25C = [Taut_prev_DNA_High_pH_25C.kGStoES1;Taut_prev_RNA_High_pH_25C.kGStoES1;...
    Taut_prev_Hybrid_High_pH_25C.kGStoES1;...
    Taut_DNA_6p9_25C.kGStoES1;Taut_RNA_6p9_25C.kGStoES1];
Taut_k12_avg_25C = mean(Taut_k12_all_25C);
Taut_k12_std_25C = std(Taut_k12_all_25C);

Taut_k21_all_25C = [Taut_prev_DNA_High_pH_25C.kES1toGS;Taut_prev_RNA_High_pH_25C.kES1toGS;...
    Taut_prev_Hybrid_High_pH_25C.kES1toGS;...
    Taut_DNA_6p9_25C.kES1toGS;Taut_RNA_6p9_25C.kES1toGS];
Taut_k21_avg_25C = mean(Taut_k21_all_25C);
Taut_k21_std_25C = std(Taut_k21_all_25C);

path_to_save_Taut = '/Users/orsula1/OrsDocs/OrsDocs_since_Jan2025/2025-Sequence-dependence/2025-Data-Deposition/Previous-R1rho-Params/Prev-Tautomer-Populations';
filename_Taut_10C = sprintf('%s/Taut_pops_and_kex_10C.mat',path_to_save_Taut);
filename_Taut_25C = sprintf('%s/Taut_pops_and_kex_25C.mat',path_to_save_Taut);

save(filename_Taut_10C,'Taut_prev_DNA_High_pH_10C','Taut_prev_RNA_High_pH_10C',...
    'Taut_DNA_6p9_10C','Taut_RNA_6p9_10C',...
    'Taut_pop_all_10C','Taut_k12_all_10C','Taut_k21_all_10C',...
    'Taut_pop_avg_10C','Taut_pop_std_10C',...
    'Taut_k12_avg_10C','Taut_k12_std_10C',...
    'Taut_k21_avg_10C','Taut_k21_std_10C');

save(filename_Taut_25C,'Taut_prev_DNA_High_pH_25C','Taut_prev_RNA_High_pH_25C',...
    'Taut_DNA_6p9_25C','Taut_RNA_6p9_25C',...
    'Taut_pop_all_25C','Taut_k12_all_25C','Taut_k21_all_25C',...
    'Taut_pop_avg_25C','Taut_pop_std_25C',...
    'Taut_k12_avg_25C','Taut_k12_std_25C',...
    'Taut_k21_avg_25C','Taut_k21_std_25C');


%% Previously meawsured Anion (High pH, interpolation to different pHs)

parent_path_output = '/Users/orsula1/OrsDocs/OrsDocs_since_Jan2025/2025-Sequence-dependence/2025-Data-Deposition/Previous-R1rho-Params/Prev-Anion-Populations';

%% parameters for pH interpolation, DNA, 10C

pKas_prev_High_pH_10C_DNA_all      = pKas_prev_High_pH_10C(ind_DNA_10C);
pKas_prev_High_pH_10C_DNA_all_errs = pKas_prev_High_pH_10C_errs(ind_DNA_10C,:);

% DNA
pKas_prev_High_pH_10C_CGC_mean    = mean(pKas_prev_High_pH_10C_DNA_all(strcmp(seq_prev_DNA_High_pH_10C,'CGC')));
pKas_prev_High_pH_10C_CGC_std     = std(pKas_prev_High_pH_10C_DNA_all(strcmp(seq_prev_DNA_High_pH_10C,'CGC')));
pKas_prev_High_pH_10C_GGC_mean    = mean(pKas_prev_High_pH_10C_DNA_all(strcmp(seq_prev_DNA_High_pH_10C,'GGC')));
pKas_prev_High_pH_10C_GGC_std     = std(pKas_prev_High_pH_10C_DNA_all(strcmp(seq_prev_DNA_High_pH_10C,'GGC')));

pKas_prev_High_pH_10C_DNA         = [pKas_prev_High_pH_10C_DNA_all(1);pKas_prev_High_pH_10C_GGC_mean;...
    pKas_prev_High_pH_10C_DNA_all(4);pKas_prev_High_pH_10C_CGC_mean];
pKas_prev_High_pH_10C_DNA_errs    = [pKas_prev_High_pH_10C_DNA_all_errs(1);pKas_prev_High_pH_10C_GGC_std;...
    pKas_prev_High_pH_10C_DNA_all_errs(4);pKas_prev_High_pH_10C_CGC_std];
seq_prev_High_pH_10C_DNA          = [seq_prev_DNA_High_pH_10C(1);seq_prev_DNA_High_pH_10C(2);...
    seq_prev_DNA_High_pH_10C(4);seq_prev_DNA_High_pH_10C(5)];
seq_prev_High_pH_10C_DNA_T        = convSeqtoT(seq_prev_High_pH_10C_DNA);

params_High_pH_10C_DNA_sort_tbl = params_High_pH_10C_sort_tbl(ind_DNA_10C,:);

% parameters for kex calculation
params_prev_High_pH_10C_DNA       = params_High_pH_10C_DNA_sort_tbl([1,3,4,5],4:end);
params_prev_High_pH_10C_DNA       = addvars(params_prev_High_pH_10C_DNA,seq_prev_High_pH_10C_DNA_T,'Before','orig_pH');
params_prev_High_pH_10C_DNA       = addvars(params_prev_High_pH_10C_DNA,pKas_prev_High_pH_10C_DNA,'Before','orig_pH');
params_prev_High_pH_10C_DNA       = addvars(params_prev_High_pH_10C_DNA,pKas_prev_High_pH_10C_DNA_errs,'Before','orig_pH');

seq_R1rho_T   = params_prev_High_pH_10C_DNA.seq_prev_High_pH_10C_DNA_T;
num_seqs      = height(params_prev_High_pH_10C_DNA);

pH            = params_prev_High_pH_10C_DNA.orig_pH;

pAnion        = 100.*params_prev_High_pH_10C_DNA.pES2_Anion;
pAnion_err    = 100.*params_prev_High_pH_10C_DNA.pES2_Anion_Error;

pTaut         = 100.*params_prev_High_pH_10C_DNA.pES1_Tautomer;
pTaut_err     = 100.*params_prev_High_pH_10C_DNA.pES1_Tautomer_Error;

pGS           = 100.*params_prev_High_pH_10C_DNA.pGS_Wobble;
pGS_err       = 100.*params_prev_High_pH_10C_DNA.pGS_Wobble_Error;

pKas          = params_prev_High_pH_10C_DNA.pKas_prev_High_pH_10C_DNA;
pKas_err      = params_prev_High_pH_10C_DNA.pKas_prev_High_pH_10C_DNA_errs;

kGStoES1      = params_prev_High_pH_10C_DNA.kGStoES1;
kGStoES1_err  = params_prev_High_pH_10C_DNA.kGStoES1_Error;
kES1toGS      = params_prev_High_pH_10C_DNA.kES1toGS;
kES1toGS_err  = params_prev_High_pH_10C_DNA.kES1toGS_Error;

kGStoES2      = params_prev_High_pH_10C_DNA.kGStoES2;
kGStoES2_err  = params_prev_High_pH_10C_DNA.kGStoES2_Error;
kES2toGS      = params_prev_High_pH_10C_DNA.kES2toGS;
kES2toGS_err  = params_prev_High_pH_10C_DNA.kES2toGS_Error;

kES2toES1     = params_prev_High_pH_10C_DNA.kES2toES1;
kES2toES1_err = params_prev_High_pH_10C_DNA.kES2toES1_Error;
kES1toES2     = params_prev_High_pH_10C_DNA.kES1toES2;
kES1toES2_err = params_prev_High_pH_10C_DNA.kES1toES2_Error;

%% load interpolation parameters

parent_path_interp = '/Users/orsula1/OrsDocs/OrsDocs_since_Jan2025/2025-Sequence-dependence/2025-Data-Deposition/Previous-R1rho-Params/Prev-Anion-Populations';
filename_interp    = 'pH_interp_params.mat';
pathname_interp    = sprintf('%s/%s',parent_path_interp,filename_interp);

load(pathname_interp);

%% interpolate original params to pH 6.9

pH_for_interp       = 6.9.*ones(num_seqs,1);
pH_for_interp_errs  = zeros(num_seqs,1);

interp_params_DNA_6p9_10C = pH_interpolation(seq_R1rho_T, pH_for_interp, pH_for_interp_errs, ...
    pH, pKas, pKas_err,...
    pGS, pGS_err, pTaut, pTaut_err, pAnion, pAnion_err, ...
    kGStoES1, kGStoES1_err, kES1toGS, kES1toGS_err, ...
    kGStoES2, kGStoES2_err, kES2toGS, kES2toGS_err, ...
    kES2toES1, kES2toES1_err, kES1toES2, kES1toES2_err, ...
    slopes_CGC_rates_10C, slopes_GGC_kminor_10C);

interp_params_DNA_6p9_10C.kMinor_pH_err(isnan(interp_params_DNA_6p9_10C.kMinor_pH_err)) = 0;
interp_params_DNA_6p9_10C.kES1toES2_pH_err(isnan(interp_params_DNA_6p9_10C.kES1toES2_pH_err)) = 0;
interp_params_DNA_6p9_10C.kES2toES1_pH_err(isnan(interp_params_DNA_6p9_10C.kES2toES1_pH_err)) = 0;

filename_input = 'DNA_interp_to_6p9_10C.csv';
pathname_input = sprintf('%s/%s',parent_path_output,filename_input);

writetable(interp_params_DNA_6p9_10C,pathname_input);

%% interpolate original params to pH 7.0

pH_for_interp       = 7.0.*ones(num_seqs,1);
pH_for_interp_errs  = zeros(num_seqs,1);

interp_params_DNA_7p0_10C = pH_interpolation(seq_R1rho_T, pH_for_interp, pH_for_interp_errs, ...
    pH, pKas, pKas_err,...
    pGS, pGS_err, pTaut, pTaut_err, pAnion, pAnion_err, ...
    kGStoES1, kGStoES1_err, kES1toGS, kES1toGS_err, ...
    kGStoES2, kGStoES2_err, kES2toGS, kES2toGS_err, ...
    kES2toES1, kES2toES1_err, kES1toES2, kES1toES2_err, ...
    slopes_CGC_rates_10C, slopes_GGC_kminor_10C);

interp_params_DNA_7p0_10C.kMinor_pH_err(isnan(interp_params_DNA_7p0_10C.kMinor_pH_err)) = 0;
interp_params_DNA_7p0_10C.kES1toES2_pH_err(isnan(interp_params_DNA_7p0_10C.kES1toES2_pH_err)) = 0;
interp_params_DNA_7p0_10C.kES2toES1_pH_err(isnan(interp_params_DNA_7p0_10C.kES2toES1_pH_err)) = 0;

filename_input = 'DNA_interp_to_7p0_10C.csv';
pathname_input = sprintf('%s/%s',parent_path_output,filename_input);

writetable(interp_params_DNA_7p0_10C,pathname_input);

%% interpolate original params to pH 7.2

pH_for_interp       = 7.2.*ones(num_seqs,1);
pH_for_interp_errs  = zeros(num_seqs,1);

interp_params_DNA_7p2_10C = pH_interpolation(seq_R1rho_T, pH_for_interp, pH_for_interp_errs, ...
    pH, pKas, pKas_err,...
    pGS, pGS_err, pTaut, pTaut_err, pAnion, pAnion_err, ...
    kGStoES1, kGStoES1_err, kES1toGS, kES1toGS_err, ...
    kGStoES2, kGStoES2_err, kES2toGS, kES2toGS_err, ...
    kES2toES1, kES2toES1_err, kES1toES2, kES1toES2_err, ...
    slopes_CGC_rates_10C, slopes_GGC_kminor_10C);

interp_params_DNA_7p2_10C.kMinor_pH_err(isnan(interp_params_DNA_7p2_10C.kMinor_pH_err)) = 0;
interp_params_DNA_7p2_10C.kES1toES2_pH_err(isnan(interp_params_DNA_7p2_10C.kES1toES2_pH_err)) = 0;
interp_params_DNA_7p2_10C.kES2toES1_pH_err(isnan(interp_params_DNA_7p2_10C.kES2toES1_pH_err)) = 0;

filename_input = 'DNA_interp_to_7p2_10C.csv';
pathname_input = sprintf('%s/%s',parent_path_output,filename_input);

writetable(interp_params_DNA_7p2_10C,pathname_input);

%% interpolate original params to pH 7.4

pH_for_interp       = 7.4.*ones(num_seqs,1);
pH_for_interp_errs  = zeros(num_seqs,1);

interp_params_DNA_7p4_10C = pH_interpolation(seq_R1rho_T, pH_for_interp, pH_for_interp_errs, ...
    pH, pKas, pKas_err,...
    pGS, pGS_err, pTaut, pTaut_err, pAnion, pAnion_err, ...
    kGStoES1, kGStoES1_err, kES1toGS, kES1toGS_err, ...
    kGStoES2, kGStoES2_err, kES2toGS, kES2toGS_err, ...
    kES2toES1, kES2toES1_err, kES1toES2, kES1toES2_err, ...
    slopes_CGC_rates_10C, slopes_GGC_kminor_10C);

interp_params_DNA_7p4_10C.kMinor_pH_err(isnan(interp_params_DNA_7p4_10C.kMinor_pH_err)) = 0;
interp_params_DNA_7p4_10C.kES1toES2_pH_err(isnan(interp_params_DNA_7p4_10C.kES1toES2_pH_err)) = 0;
interp_params_DNA_7p4_10C.kES2toES1_pH_err(isnan(interp_params_DNA_7p4_10C.kES2toES1_pH_err)) = 0;

filename_input = 'DNA_interp_to_7p4_10C.csv';
pathname_input = sprintf('%s/%s',parent_path_output,filename_input);

writetable(interp_params_DNA_7p4_10C,pathname_input);


%% interpolate original params to pH 7.5

pH_for_interp       = 7.5.*ones(num_seqs,1);
pH_for_interp_errs  = zeros(num_seqs,1);

interp_params_DNA_7p5_10C = pH_interpolation(seq_R1rho_T, pH_for_interp, pH_for_interp_errs, ...
    pH, pKas, pKas_err,...
    pGS, pGS_err, pTaut, pTaut_err, pAnion, pAnion_err, ...
    kGStoES1, kGStoES1_err, kES1toGS, kES1toGS_err, ...
    kGStoES2, kGStoES2_err, kES2toGS, kES2toGS_err, ...
    kES2toES1, kES2toES1_err, kES1toES2, kES1toES2_err, ...
    slopes_CGC_rates_10C, slopes_GGC_kminor_10C);

interp_params_DNA_7p5_10C.kMinor_pH_err(isnan(interp_params_DNA_7p5_10C.kMinor_pH_err)) = 0;
interp_params_DNA_7p5_10C.kES1toES2_pH_err(isnan(interp_params_DNA_7p5_10C.kES1toES2_pH_err)) = 0;
interp_params_DNA_7p5_10C.kES2toES1_pH_err(isnan(interp_params_DNA_7p5_10C.kES2toES1_pH_err)) = 0;

filename_input = 'DNA_interp_to_7p5_10C.csv';
pathname_input = sprintf('%s/%s',parent_path_output,filename_input);

writetable(interp_params_DNA_7p5_10C,pathname_input);


%% interpolate original params to pH 7.6

pH_for_interp       = 7.6.*ones(num_seqs,1);
pH_for_interp_errs  = zeros(num_seqs,1);

interp_params_DNA_7p6_10C = pH_interpolation(seq_R1rho_T, pH_for_interp, pH_for_interp_errs, ...
    pH, pKas, pKas_err,...
    pGS, pGS_err, pTaut, pTaut_err, pAnion, pAnion_err, ...
    kGStoES1, kGStoES1_err, kES1toGS, kES1toGS_err, ...
    kGStoES2, kGStoES2_err, kES2toGS, kES2toGS_err, ...
    kES2toES1, kES2toES1_err, kES1toES2, kES1toES2_err, ...
    slopes_CGC_rates_10C, slopes_GGC_kminor_10C);

interp_params_DNA_7p6_10C.kMinor_pH_err(isnan(interp_params_DNA_7p6_10C.kMinor_pH_err)) = 0;
interp_params_DNA_7p6_10C.kES1toES2_pH_err(isnan(interp_params_DNA_7p6_10C.kES1toES2_pH_err)) = 0;
interp_params_DNA_7p6_10C.kES2toES1_pH_err(isnan(interp_params_DNA_7p6_10C.kES2toES1_pH_err)) = 0;

filename_input = 'DNA_interp_to_7p6_10C.csv';
pathname_input = sprintf('%s/%s',parent_path_output,filename_input);

writetable(interp_params_DNA_7p6_10C,pathname_input);

%% interpolate original params to pH 8.0

pH_for_interp       = 8.0.*ones(num_seqs,1);
pH_for_interp_errs  = zeros(num_seqs,1);

interp_params_DNA_8p0_10C = pH_interpolation(seq_R1rho_T, pH_for_interp, pH_for_interp_errs, ...
    pH, pKas, pKas_err,...
    pGS, pGS_err, pTaut, pTaut_err, pAnion, pAnion_err, ...
    kGStoES1, kGStoES1_err, kES1toGS, kES1toGS_err, ...
    kGStoES2, kGStoES2_err, kES2toGS, kES2toGS_err, ...
    kES2toES1, kES2toES1_err, kES1toES2, kES1toES2_err, ...
    slopes_CGC_rates_10C, slopes_GGC_kminor_10C);

interp_params_DNA_8p0_10C.kMinor_pH_err(isnan(interp_params_DNA_8p0_10C.kMinor_pH_err)) = 0;
interp_params_DNA_8p0_10C.kES1toES2_pH_err(isnan(interp_params_DNA_8p0_10C.kES1toES2_pH_err)) = 0;
interp_params_DNA_8p0_10C.kES2toES1_pH_err(isnan(interp_params_DNA_8p0_10C.kES2toES1_pH_err)) = 0;

filename_input = 'DNA_interp_to_8p0_10C.csv';
pathname_input = sprintf('%s/%s',parent_path_output,filename_input);

writetable(interp_params_DNA_8p0_10C,pathname_input);

%% interpolate original params to pH 8.4

pH_for_interp       = 8.4.*ones(num_seqs,1);
pH_for_interp_errs  = zeros(num_seqs,1);

interp_params_DNA_8p4_10C = pH_interpolation(seq_R1rho_T, pH_for_interp, pH_for_interp_errs, ...
    pH, pKas, pKas_err,...
    pGS, pGS_err, pTaut, pTaut_err, pAnion, pAnion_err, ...
    kGStoES1, kGStoES1_err, kES1toGS, kES1toGS_err, ...
    kGStoES2, kGStoES2_err, kES2toGS, kES2toGS_err, ...
    kES2toES1, kES2toES1_err, kES1toES2, kES1toES2_err, ...
    slopes_CGC_rates_10C, slopes_GGC_kminor_10C);

interp_params_DNA_8p4_10C.kMinor_pH_err(isnan(interp_params_DNA_8p4_10C.kMinor_pH_err)) = 0;
interp_params_DNA_8p4_10C.kES1toES2_pH_err(isnan(interp_params_DNA_8p4_10C.kES1toES2_pH_err)) = 0;
interp_params_DNA_8p4_10C.kES2toES1_pH_err(isnan(interp_params_DNA_8p4_10C.kES2toES1_pH_err)) = 0;

filename_input = 'DNA_interp_to_8p4_10C.csv';
pathname_input = sprintf('%s/%s',parent_path_output,filename_input);

writetable(interp_params_DNA_8p4_10C,pathname_input);


%% dG calculations

R_constant = 1.9872036E-3; %kCal/molK
Temp_25C   = 298.15; %K
Temp_10C   = 283.15; %K
Temp_1C    = 274.15; %K

%% dG Anion calculation

pH_for_calc = 7.4;

pAnion_prev_pH_7p4_10C      = (10.^(pH_for_calc-pKas_prev_High_pH_10C_DNA))./(1+10.^(pH_for_calc-pKas_prev_High_pH_10C_DNA));
pAnion_prev_pH_7p4_10C_errs = pAnion_prev_pH_7p4_10C.*sqrt((pKas_prev_High_pH_10C_DNA_errs).^2.*(log(10)).^2+((10.^(pH_for_calc-pKas_prev_High_pH_10C_DNA)).^2.*(pKas_prev_High_pH_10C_DNA_errs).^2.*(log(10))^2)./((1+10.^(pH_for_calc-pKas_prev_High_pH_10C_DNA)).^2));

dG_Anion_prev_7p4_10C      = -R_constant*Temp_10C*log(10)*(pH_for_calc - pKas_prev_High_pH_10C_DNA);
dG_Anion_prev_7p4_10C_errs = pKas_prev_High_pH_10C_DNA_errs;

%% dG Tautomer calculation

Taut_DNA_6p9_25C_pops_for_dG = table(seq_DNA_6p9_25C_T,params_6p9_25C_sort_tbl.orig_pH(ind_DNA_6p9_25C),...
    pTaut_DNA_6p9_25C,pTaut_DNA_6p9_25C_errs,...
    pGenol_DNA_6p9_25C,pGenol_DNA_6p9_25C_errs,...
    pTenol_DNA_6p9_25C,pTenol_DNA_6p9_25C_errs);

pTaut_pH_6p9_25C      = Taut_DNA_6p9_25C_pops_for_dG.pTaut_DNA_6p9_25C;
pTaut_pH_6p9_25C_errs = Taut_DNA_6p9_25C_pops_for_dG.pTaut_DNA_6p9_25C_errs;

dG_Taut_6p9_25C       = -R_constant*Temp_25C.*log(pTaut_pH_6p9_25C./(100-pTaut_pH_6p9_25C));
dG_Taut_6p9_25C_errs  = -R_constant*Temp_25C*sqrt(2).*(pTaut_pH_6p9_25C_errs./pTaut_pH_6p9_25C);

path_to_save_dG = '/Users/orsula1/OrsDocs/OrsDocs_since_Jan2025/2025-Sequence-dependence/2025-Data-Deposition/Previous-R1rho-Params/Prev-dG-Anion-and-Tautomer';
filename_dG = sprintf('%s/dG_from_R1rho.mat',path_to_save_dG);

save(filename_dG,'params_prev_High_pH_10C_DNA','Taut_DNA_6p9_25C_pops_for_dG',...
    'pAnion_prev_pH_7p4_10C','pAnion_prev_pH_7p4_10C_errs',...
    'dG_Anion_prev_7p4_10C','dG_Anion_prev_7p4_10C_errs',...
    'pTaut_pH_6p9_25C','pTaut_pH_6p9_25C_errs',...
    'dG_Taut_6p9_25C','dG_Taut_6p9_25C_errs');