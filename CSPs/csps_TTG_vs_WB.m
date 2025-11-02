clear;
clc;

%% colors and axes

sw_HC_HSQC_600_Duke = 28.0025942735608;
sw_HC_HSQC_700_UNC  = 28.0047826032834;

%% loading data file paths

parent_path     = '/Users/orsula1/OrsDocs/OrsDocs_since_Jan2025/2025-Sequence-dependence/2025-Data-Deposition/CSPs';
TTG_pathname    = 'TTG-peak-lists';

filename_GC_WC_TTG_pH6p8   = 'TTG_GC_WC_pH6p8_shifts';
filename_G5FdU_TTG_pH6p8    = 'TTG_G5FdU_pH6p8_shifts';
filename_G5FdU_TTG_pH10p1   = 'TTG_G5FdU_pH10p1_shifts';
filename_A5FdU_TTG_pH6p8    = 'TTG_A5FdU_WC_pH6p8_shifts';
filename_isoG5FdU_TTG_pH6p8 = 'TTG_isoG5FdU_pH6p8_shifts';

%% GC WC TTG, pH 6.8

filepath_GC_WC_TTG_6p8  = sprintf('%s/%s/%s',parent_path,TTG_pathname,filename_GC_WC_TTG_pH6p8);

GC_WC_TTG_6p8_HSQC_list = readtable(filepath_GC_WC_TTG_6p8);

Assign_GC_WC_TTG_6p8    = GC_WC_TTG_6p8_HSQC_list.Atom;

ind_GC_WC_TTG_6p8_C6C8  = zeros(length(Assign_GC_WC_TTG_6p8),1);
ind_GC_WC_TTG_6p8_C2    = zeros(length(Assign_GC_WC_TTG_6p8),1);
ind_GC_WC_TTG_6p8_C1p   = zeros(length(Assign_GC_WC_TTG_6p8),1);
ind_GC_WC_TTG_6p8_C4p   = zeros(length(Assign_GC_WC_TTG_6p8),1);
ind_GC_WC_TTG_6p8_H6H8  = zeros(length(Assign_GC_WC_TTG_6p8),1);
ind_GC_WC_TTG_6p8_H2    = zeros(length(Assign_GC_WC_TTG_6p8),1);
ind_GC_WC_TTG_6p8_H1p   = zeros(length(Assign_GC_WC_TTG_6p8),1);
ind_GC_WC_TTG_6p8_H4p   = zeros(length(Assign_GC_WC_TTG_6p8),1);

for i = 1:length(Assign_GC_WC_TTG_6p8)

    if ismember('C6',Assign_GC_WC_TTG_6p8(i,1)) || ismember('C8',Assign_GC_WC_TTG_6p8(i,1))

        ind_GC_WC_TTG_6p8_C6C8(i,1) = i;

    elseif ismember('C2',Assign_GC_WC_TTG_6p8(i,1))

        ind_GC_WC_TTG_6p8_C2(i,1)   = i;

    elseif ismember("C1'",Assign_GC_WC_TTG_6p8(i,1))

        ind_GC_WC_TTG_6p8_C1p(i,1)   = i;

    elseif ismember("C4'",Assign_GC_WC_TTG_6p8(i,1))

        ind_GC_WC_TTG_6p8_C4p(i,1)   = i;

    elseif ismember('H6',Assign_GC_WC_TTG_6p8(i,1)) || ismember('H8',Assign_GC_WC_TTG_6p8(i,1))

        ind_GC_WC_TTG_6p8_H6H8(i,1)   = i;

    elseif ismember('H2',Assign_GC_WC_TTG_6p8(i,1))

        ind_GC_WC_TTG_6p8_H2(i,1)   = i;

    elseif ismember("H1'",Assign_GC_WC_TTG_6p8(i,1))

        ind_GC_WC_TTG_6p8_H1p(i,1)   = i;

    elseif ismember("H4'",Assign_GC_WC_TTG_6p8(i,1))

        ind_GC_WC_TTG_6p8_H4p(i,1)   = i;
        
    end

end

ind_GC_WC_TTG_6p8_C6C8 = ind_GC_WC_TTG_6p8_C6C8(ind_GC_WC_TTG_6p8_C6C8~=0);
ind_GC_WC_TTG_6p8_C2   = ind_GC_WC_TTG_6p8_C2(ind_GC_WC_TTG_6p8_C2~=0);
ind_GC_WC_TTG_6p8_C1p  = ind_GC_WC_TTG_6p8_C1p(ind_GC_WC_TTG_6p8_C1p~=0);
ind_GC_WC_TTG_6p8_C4p  = ind_GC_WC_TTG_6p8_C4p(ind_GC_WC_TTG_6p8_C4p~=0);
ind_GC_WC_TTG_6p8_H6H8 = ind_GC_WC_TTG_6p8_H6H8(ind_GC_WC_TTG_6p8_H6H8~=0);
ind_GC_WC_TTG_6p8_H2   = ind_GC_WC_TTG_6p8_H2(ind_GC_WC_TTG_6p8_H2~=0);
ind_GC_WC_TTG_6p8_H1p  = ind_GC_WC_TTG_6p8_H1p(ind_GC_WC_TTG_6p8_H1p~=0);
ind_GC_WC_TTG_6p8_H4p  = ind_GC_WC_TTG_6p8_H4p(ind_GC_WC_TTG_6p8_H4p~=0);

GC_WC_TTG_6p8_HSQC_C6C8_list = GC_WC_TTG_6p8_HSQC_list(ind_GC_WC_TTG_6p8_C6C8,1:4);
GC_WC_TTG_6p8_HSQC_C2_list   = GC_WC_TTG_6p8_HSQC_list(ind_GC_WC_TTG_6p8_C2,1:4);
GC_WC_TTG_6p8_HSQC_C1p_list  = GC_WC_TTG_6p8_HSQC_list(ind_GC_WC_TTG_6p8_C1p,1:4);
GC_WC_TTG_6p8_HSQC_C4p_list  = GC_WC_TTG_6p8_HSQC_list(ind_GC_WC_TTG_6p8_C4p,1:4);
GC_WC_TTG_6p8_HSQC_H6H8_list = GC_WC_TTG_6p8_HSQC_list(ind_GC_WC_TTG_6p8_H6H8,1:4);
GC_WC_TTG_6p8_HSQC_H2_list   = GC_WC_TTG_6p8_HSQC_list(ind_GC_WC_TTG_6p8_H2,1:4);
GC_WC_TTG_6p8_HSQC_H1p_list  = GC_WC_TTG_6p8_HSQC_list(ind_GC_WC_TTG_6p8_H1p,1:4);
GC_WC_TTG_6p8_HSQC_H4p_list  = GC_WC_TTG_6p8_HSQC_list(ind_GC_WC_TTG_6p8_H4p,1:4);

GC_WC_TTG_6p8_resname_C6C8   = GC_WC_TTG_6p8_HSQC_C6C8_list.Group;
GC_WC_TTG_6p8_resnum_C6C8    = zeros(height(GC_WC_TTG_6p8_HSQC_C6C8_list),1);

for i = 1:height(GC_WC_TTG_6p8_HSQC_C6C8_list)

    GC_WC_TTG_6p8_resnum_C6C8(i,1) = str2double(regexp(GC_WC_TTG_6p8_resname_C6C8{i,1}, '(+|-)?(\d+)?\.?\d*', 'match'));

end

GC_WC_TTG_6p8_resname_C2     = GC_WC_TTG_6p8_HSQC_C2_list.Group;
GC_WC_TTG_6p8_resnum_C2      = zeros(height(GC_WC_TTG_6p8_HSQC_C2_list),1);

for i = 1:height(GC_WC_TTG_6p8_HSQC_C2_list)

    GC_WC_TTG_6p8_resnum_C2(i,1) = str2double(regexp(GC_WC_TTG_6p8_resname_C2{i,1}, '(+|-)?(\d+)?\.?\d*', 'match'));

end

GC_WC_TTG_6p8_resname_C1p    = GC_WC_TTG_6p8_HSQC_C1p_list.Group;
GC_WC_TTG_6p8_resnum_C1p     = zeros(height(GC_WC_TTG_6p8_HSQC_C1p_list),1);

for i = 1:height(GC_WC_TTG_6p8_HSQC_C1p_list)

    GC_WC_TTG_6p8_resnum_C1p(i,1) = str2double(regexp(GC_WC_TTG_6p8_resname_C1p{i,1}, '(+|-)?(\d+)?\.?\d*', 'match'));

end

GC_WC_TTG_6p8_resname_C4p    = GC_WC_TTG_6p8_HSQC_C4p_list.Group;
GC_WC_TTG_6p8_resnum_C4p     = zeros(height(GC_WC_TTG_6p8_HSQC_C4p_list),1);

for i = 1:height(GC_WC_TTG_6p8_HSQC_C4p_list)

    GC_WC_TTG_6p8_resnum_C4p(i,1) = str2double(regexp(GC_WC_TTG_6p8_resname_C4p{i,1}, '(+|-)?(\d+)?\.?\d*', 'match'));

end


%% G5FdU TTG, pH 6.8

filepath_G5FdU_TTG_6p8  = sprintf('%s/%s/%s',parent_path,TTG_pathname,filename_G5FdU_TTG_pH6p8);

G5FdU_TTG_6p8_HSQC_list = readtable(filepath_G5FdU_TTG_6p8);

Assign_G5FdU_TTG_6p8    = G5FdU_TTG_6p8_HSQC_list.Atom;

ind_G5FdU_TTG_6p8_C6C8  = zeros(length(Assign_G5FdU_TTG_6p8),1);
ind_G5FdU_TTG_6p8_C2    = zeros(length(Assign_G5FdU_TTG_6p8),1);
ind_G5FdU_TTG_6p8_C1p   = zeros(length(Assign_G5FdU_TTG_6p8),1);
ind_G5FdU_TTG_6p8_C4p   = zeros(length(Assign_G5FdU_TTG_6p8),1);
ind_G5FdU_TTG_6p8_H6H8  = zeros(length(Assign_G5FdU_TTG_6p8),1);
ind_G5FdU_TTG_6p8_H2    = zeros(length(Assign_G5FdU_TTG_6p8),1);
ind_G5FdU_TTG_6p8_H1p   = zeros(length(Assign_G5FdU_TTG_6p8),1);
ind_G5FdU_TTG_6p8_H4p   = zeros(length(Assign_G5FdU_TTG_6p8),1);

for i = 1:length(Assign_G5FdU_TTG_6p8)

    if ismember('C6',Assign_G5FdU_TTG_6p8(i,1)) || ismember('C8',Assign_G5FdU_TTG_6p8(i,1))

        ind_G5FdU_TTG_6p8_C6C8(i,1) = i;

    elseif ismember('C2',Assign_G5FdU_TTG_6p8(i,1))

        ind_G5FdU_TTG_6p8_C2(i,1)   = i;

    elseif ismember("C1'",Assign_G5FdU_TTG_6p8(i,1))

        ind_G5FdU_TTG_6p8_C1p(i,1)   = i;

    elseif ismember("C4'",Assign_G5FdU_TTG_6p8(i,1))

        ind_G5FdU_TTG_6p8_C4p(i,1)   = i;
        
    elseif ismember('H6',Assign_G5FdU_TTG_6p8(i,1)) || ismember('H8',Assign_G5FdU_TTG_6p8(i,1))

        ind_G5FdU_TTG_6p8_H6H8(i,1)   = i;

    elseif ismember('H2',Assign_G5FdU_TTG_6p8(i,1))

        ind_G5FdU_TTG_6p8_H2(i,1)   = i;

    elseif ismember("H1'",Assign_G5FdU_TTG_6p8(i,1))

        ind_G5FdU_TTG_6p8_H1p(i,1)   = i;

    elseif ismember("H4'",Assign_G5FdU_TTG_6p8(i,1))

        ind_G5FdU_TTG_6p8_H4p(i,1)   = i;

    end

end

ind_G5FdU_TTG_6p8_C6C8 = ind_G5FdU_TTG_6p8_C6C8(ind_G5FdU_TTG_6p8_C6C8~=0);
ind_G5FdU_TTG_6p8_C2   = ind_G5FdU_TTG_6p8_C2(ind_G5FdU_TTG_6p8_C2~=0);
ind_G5FdU_TTG_6p8_C1p  = ind_G5FdU_TTG_6p8_C1p(ind_G5FdU_TTG_6p8_C1p~=0);
ind_G5FdU_TTG_6p8_C4p  = ind_G5FdU_TTG_6p8_C4p(ind_G5FdU_TTG_6p8_C4p~=0);
ind_G5FdU_TTG_6p8_H6H8 = ind_G5FdU_TTG_6p8_H6H8(ind_G5FdU_TTG_6p8_H6H8~=0);
ind_G5FdU_TTG_6p8_H2   = ind_G5FdU_TTG_6p8_H2(ind_G5FdU_TTG_6p8_H2~=0);
ind_G5FdU_TTG_6p8_H1p  = ind_G5FdU_TTG_6p8_H1p(ind_G5FdU_TTG_6p8_H1p~=0);
ind_G5FdU_TTG_6p8_H4p  = ind_G5FdU_TTG_6p8_H4p(ind_G5FdU_TTG_6p8_H4p~=0);

G5FdU_TTG_6p8_HSQC_C6C8_list = G5FdU_TTG_6p8_HSQC_list(ind_G5FdU_TTG_6p8_C6C8,1:4);
G5FdU_TTG_6p8_HSQC_C2_list   = G5FdU_TTG_6p8_HSQC_list(ind_G5FdU_TTG_6p8_C2,1:4);
G5FdU_TTG_6p8_HSQC_C1p_list  = G5FdU_TTG_6p8_HSQC_list(ind_G5FdU_TTG_6p8_C1p,1:4);
G5FdU_TTG_6p8_HSQC_C4p_list  = G5FdU_TTG_6p8_HSQC_list(ind_G5FdU_TTG_6p8_C4p,1:4);
G5FdU_TTG_6p8_HSQC_H6H8_list = G5FdU_TTG_6p8_HSQC_list(ind_G5FdU_TTG_6p8_H6H8,1:4);
G5FdU_TTG_6p8_HSQC_H2_list   = G5FdU_TTG_6p8_HSQC_list(ind_G5FdU_TTG_6p8_H2,1:4);
G5FdU_TTG_6p8_HSQC_H1p_list  = G5FdU_TTG_6p8_HSQC_list(ind_G5FdU_TTG_6p8_H1p,1:4);
G5FdU_TTG_6p8_HSQC_H4p_list  = G5FdU_TTG_6p8_HSQC_list(ind_G5FdU_TTG_6p8_H4p,1:4);

G5FdU_TTG_6p8_resname_C6C8   = G5FdU_TTG_6p8_HSQC_C6C8_list.Group;
G5FdU_TTG_6p8_resnum_C6C8    = zeros(height(G5FdU_TTG_6p8_HSQC_C6C8_list),1);

for i = 1:height(G5FdU_TTG_6p8_HSQC_C6C8_list)

    G5FdU_TTG_6p8_resnum_C6C8(i,1) = str2double(regexp(G5FdU_TTG_6p8_resname_C6C8{i,1}, '(+|-)?(\d+)?\.?\d*', 'match'));

end

G5FdU_TTG_6p8_U5_C6_folded_shift = G5FdU_TTG_6p8_HSQC_C6C8_list(ismember(G5FdU_TTG_6p8_resname_C6C8,'U5'),:).Shift;
G5FdU_TTG_6p8_U5_C6_new_shift    = G5FdU_TTG_6p8_U5_C6_folded_shift - sw_HC_HSQC_600_Duke;
G5FdU_TTG_6p8_HSQC_C6C8_list(ismember(G5FdU_TTG_6p8_resname_C6C8,'U5'),:).Shift = G5FdU_TTG_6p8_U5_C6_new_shift;

G5FdU_TTG_6p8_resname_C2     = G5FdU_TTG_6p8_HSQC_C2_list.Group;
G5FdU_TTG_6p8_resnum_C2      = zeros(height(G5FdU_TTG_6p8_HSQC_C2_list),1);

for i = 1:height(G5FdU_TTG_6p8_HSQC_C2_list)

    G5FdU_TTG_6p8_resnum_C2(i,1) = str2double(regexp(G5FdU_TTG_6p8_resname_C2{i,1}, '(+|-)?(\d+)?\.?\d*', 'match'));

end

G5FdU_TTG_6p8_resname_C1p    = G5FdU_TTG_6p8_HSQC_C1p_list.Group;
G5FdU_TTG_6p8_resnum_C1p     = zeros(height(G5FdU_TTG_6p8_HSQC_C1p_list),1);

for i = 1:height(G5FdU_TTG_6p8_HSQC_C1p_list)

    G5FdU_TTG_6p8_resnum_C1p(i,1) = str2double(regexp(G5FdU_TTG_6p8_resname_C1p{i,1}, '(+|-)?(\d+)?\.?\d*', 'match'));

end

G5FdU_TTG_6p8_resname_C4p    = G5FdU_TTG_6p8_HSQC_C4p_list.Group;
G5FdU_TTG_6p8_resnum_C4p     = zeros(height(G5FdU_TTG_6p8_HSQC_C4p_list),1);

for i = 1:height(G5FdU_TTG_6p8_HSQC_C4p_list)

    G5FdU_TTG_6p8_resnum_C4p(i,1) = str2double(regexp(G5FdU_TTG_6p8_resname_C4p{i,1}, '(+|-)?(\d+)?\.?\d*', 'match'));

end


%% G5FdU TTG, pH 10.1

filepath_G5FdU_TTG_10p1  = sprintf('%s/%s/%s',parent_path,TTG_pathname,filename_G5FdU_TTG_pH10p1);

G5FdU_TTG_10p1_HSQC_list = readtable(filepath_G5FdU_TTG_10p1);

Assign_G5FdU_TTG_10p1    = G5FdU_TTG_10p1_HSQC_list.Atom;

ind_G5FdU_TTG_10p1_C6C8  = zeros(length(Assign_G5FdU_TTG_10p1),1);
ind_G5FdU_TTG_10p1_C2    = zeros(length(Assign_G5FdU_TTG_10p1),1);
ind_G5FdU_TTG_10p1_C1p   = zeros(length(Assign_G5FdU_TTG_10p1),1);
ind_G5FdU_TTG_10p1_C4p   = zeros(length(Assign_G5FdU_TTG_10p1),1);
ind_G5FdU_TTG_10p1_H6H8  = zeros(length(Assign_G5FdU_TTG_10p1),1);
ind_G5FdU_TTG_10p1_H2    = zeros(length(Assign_G5FdU_TTG_10p1),1);
ind_G5FdU_TTG_10p1_H1p   = zeros(length(Assign_G5FdU_TTG_10p1),1);
ind_G5FdU_TTG_10p1_H4p   = zeros(length(Assign_G5FdU_TTG_10p1),1);

for i = 1:length(Assign_G5FdU_TTG_10p1)

    if ismember('C6',Assign_G5FdU_TTG_10p1(i,1)) || ismember('C8',Assign_G5FdU_TTG_10p1(i,1))

        ind_G5FdU_TTG_10p1_C6C8(i,1) = i;

    elseif ismember('C2',Assign_G5FdU_TTG_10p1(i,1))

        ind_G5FdU_TTG_10p1_C2(i,1)   = i;

    elseif ismember("C1'",Assign_G5FdU_TTG_10p1(i,1))

        ind_G5FdU_TTG_10p1_C1p(i,1)   = i;

    elseif ismember("C4'",Assign_G5FdU_TTG_10p1(i,1))

        ind_G5FdU_TTG_10p1_C4p(i,1)   = i;

    elseif ismember('H6',Assign_G5FdU_TTG_10p1(i,1)) || ismember('H8',Assign_G5FdU_TTG_10p1(i,1))

        ind_G5FdU_TTG_10p1_H6H8(i,1)   = i;

    elseif ismember('H2',Assign_G5FdU_TTG_10p1(i,1))

        ind_G5FdU_TTG_10p1_H2(i,1)   = i;

    elseif ismember("H1'",Assign_G5FdU_TTG_10p1(i,1))

        ind_G5FdU_TTG_10p1_H1p(i,1)   = i;

    elseif ismember("H4'",Assign_G5FdU_TTG_10p1(i,1))

        ind_G5FdU_TTG_10p1_H4p(i,1)   = i;

    end

end

ind_G5FdU_TTG_10p1_C6C8 = ind_G5FdU_TTG_10p1_C6C8(ind_G5FdU_TTG_10p1_C6C8~=0);
ind_G5FdU_TTG_10p1_C2   = ind_G5FdU_TTG_10p1_C2(ind_G5FdU_TTG_10p1_C2~=0);
ind_G5FdU_TTG_10p1_C1p  = ind_G5FdU_TTG_10p1_C1p(ind_G5FdU_TTG_10p1_C1p~=0);
ind_G5FdU_TTG_10p1_C4p  = ind_G5FdU_TTG_10p1_C4p(ind_G5FdU_TTG_10p1_C4p~=0);
ind_G5FdU_TTG_10p1_H6H8 = ind_G5FdU_TTG_10p1_H6H8(ind_G5FdU_TTG_10p1_H6H8~=0);
ind_G5FdU_TTG_10p1_H2   = ind_G5FdU_TTG_10p1_H2(ind_G5FdU_TTG_10p1_H2~=0);
ind_G5FdU_TTG_10p1_H1p  = ind_G5FdU_TTG_10p1_H1p(ind_G5FdU_TTG_10p1_H1p~=0);
ind_G5FdU_TTG_10p1_H4p  = ind_G5FdU_TTG_10p1_H4p(ind_G5FdU_TTG_10p1_H4p~=0);

G5FdU_TTG_10p1_HSQC_C6C8_list = G5FdU_TTG_10p1_HSQC_list(ind_G5FdU_TTG_10p1_C6C8,1:4);
G5FdU_TTG_10p1_HSQC_C2_list   = G5FdU_TTG_10p1_HSQC_list(ind_G5FdU_TTG_10p1_C2,1:4);
G5FdU_TTG_10p1_HSQC_C1p_list  = G5FdU_TTG_10p1_HSQC_list(ind_G5FdU_TTG_10p1_C1p,1:4);
G5FdU_TTG_10p1_HSQC_C4p_list  = G5FdU_TTG_10p1_HSQC_list(ind_G5FdU_TTG_10p1_C4p,1:4);
G5FdU_TTG_10p1_HSQC_H6H8_list = G5FdU_TTG_10p1_HSQC_list(ind_G5FdU_TTG_10p1_H6H8,1:4);
G5FdU_TTG_10p1_HSQC_H2_list   = G5FdU_TTG_10p1_HSQC_list(ind_G5FdU_TTG_10p1_H2,1:4);
G5FdU_TTG_10p1_HSQC_H1p_list  = G5FdU_TTG_10p1_HSQC_list(ind_G5FdU_TTG_10p1_H1p,1:4);
G5FdU_TTG_10p1_HSQC_H4p_list  = G5FdU_TTG_10p1_HSQC_list(ind_G5FdU_TTG_10p1_H4p,1:4);

G5FdU_TTG_10p1_resname_C6C8   = G5FdU_TTG_10p1_HSQC_C6C8_list.Group;
G5FdU_TTG_10p1_resnum_C6C8    = zeros(height(G5FdU_TTG_10p1_HSQC_C6C8_list),1);

for i = 1:height(G5FdU_TTG_10p1_HSQC_C6C8_list)

    G5FdU_TTG_10p1_resnum_C6C8(i,1) = str2double(regexp(G5FdU_TTG_10p1_resname_C6C8{i,1}, '(+|-)?(\d+)?\.?\d*', 'match'));

end

G5FdU_TTG_10p1_U5_C6_folded_shift = G5FdU_TTG_10p1_HSQC_C6C8_list(ismember(G5FdU_TTG_10p1_resname_C6C8,'U5'),:).Shift;
G5FdU_TTG_10p1_U5_C6_new_shift    = G5FdU_TTG_10p1_U5_C6_folded_shift - sw_HC_HSQC_600_Duke;
G5FdU_TTG_10p1_HSQC_C6C8_list(ismember(G5FdU_TTG_10p1_resname_C6C8,'U5'),:).Shift = G5FdU_TTG_10p1_U5_C6_new_shift;

G5FdU_TTG_10p1_resname_C2     = G5FdU_TTG_10p1_HSQC_C2_list.Group;
G5FdU_TTG_10p1_resnum_C2      = zeros(height(G5FdU_TTG_10p1_HSQC_C2_list),1);

for i = 1:height(G5FdU_TTG_10p1_HSQC_C2_list)

    G5FdU_TTG_10p1_resnum_C2(i,1) = str2double(regexp(G5FdU_TTG_10p1_resname_C2{i,1}, '(+|-)?(\d+)?\.?\d*', 'match'));

end

G5FdU_TTG_10p1_resname_C1p    = G5FdU_TTG_10p1_HSQC_C1p_list.Group;
G5FdU_TTG_10p1_resnum_C1p     = zeros(height(G5FdU_TTG_10p1_HSQC_C1p_list),1);

for i = 1:height(G5FdU_TTG_10p1_HSQC_C1p_list)

    G5FdU_TTG_10p1_resnum_C1p(i,1) = str2double(regexp(G5FdU_TTG_10p1_resname_C1p{i,1}, '(+|-)?(\d+)?\.?\d*', 'match'));

end


G5FdU_TTG_10p1_resname_C4p    = G5FdU_TTG_10p1_HSQC_C4p_list.Group;
G5FdU_TTG_10p1_resnum_C4p     = zeros(height(G5FdU_TTG_10p1_HSQC_C4p_list),1);

for i = 1:height(G5FdU_TTG_10p1_HSQC_C4p_list)

    G5FdU_TTG_10p1_resnum_C4p(i,1) = str2double(regexp(G5FdU_TTG_10p1_resname_C4p{i,1}, '(+|-)?(\d+)?\.?\d*', 'match'));

end

%% A5FdU TTG WC, pH 6.8

filepath_A5FdU_TTG_6p8  = sprintf('%s/%s/%s',parent_path,TTG_pathname,filename_A5FdU_TTG_pH6p8);

A5FdU_TTG_6p8_HSQC_list = readtable(filepath_A5FdU_TTG_6p8);

Assign_A5FdU_TTG_6p8    = A5FdU_TTG_6p8_HSQC_list.Atom;

ind_A5FdU_TTG_6p8_C6C8  = zeros(length(Assign_A5FdU_TTG_6p8),1);
ind_A5FdU_TTG_6p8_C2    = zeros(length(Assign_A5FdU_TTG_6p8),1);
ind_A5FdU_TTG_6p8_C1p   = zeros(length(Assign_A5FdU_TTG_6p8),1);
ind_A5FdU_TTG_6p8_C4p   = zeros(length(Assign_A5FdU_TTG_6p8),1);
ind_A5FdU_TTG_6p8_H6H8  = zeros(length(Assign_A5FdU_TTG_6p8),1);
ind_A5FdU_TTG_6p8_H2    = zeros(length(Assign_A5FdU_TTG_6p8),1);
ind_A5FdU_TTG_6p8_H1p   = zeros(length(Assign_A5FdU_TTG_6p8),1);
ind_A5FdU_TTG_6p8_H4p   = zeros(length(Assign_A5FdU_TTG_6p8),1);

for i = 1:length(Assign_A5FdU_TTG_6p8)

    if ismember('C6',Assign_A5FdU_TTG_6p8(i,1)) || ismember('C8',Assign_A5FdU_TTG_6p8(i,1))

        ind_A5FdU_TTG_6p8_C6C8(i,1) = i;

    elseif ismember('C2',Assign_A5FdU_TTG_6p8(i,1))

        ind_A5FdU_TTG_6p8_C2(i,1)   = i;

    elseif ismember("C1'",Assign_A5FdU_TTG_6p8(i,1))

        ind_A5FdU_TTG_6p8_C1p(i,1)   = i;

    elseif ismember("C4'",Assign_A5FdU_TTG_6p8(i,1))

        ind_A5FdU_TTG_6p8_C4p(i,1)   = i;

    elseif ismember('H6',Assign_A5FdU_TTG_6p8(i,1)) || ismember('H8',Assign_A5FdU_TTG_6p8(i,1))

        ind_A5FdU_TTG_6p8_H6H8(i,1)   = i;

    elseif ismember('H2',Assign_A5FdU_TTG_6p8(i,1))

        ind_A5FdU_TTG_6p8_H2(i,1)   = i;

    elseif ismember("H1'",Assign_A5FdU_TTG_6p8(i,1))

        ind_A5FdU_TTG_6p8_H1p(i,1)   = i;

    elseif ismember("H4'",Assign_A5FdU_TTG_6p8(i,1))

        ind_A5FdU_TTG_6p8_H4p(i,1)   = i;

    end

end

ind_A5FdU_TTG_6p8_C6C8 = ind_A5FdU_TTG_6p8_C6C8(ind_A5FdU_TTG_6p8_C6C8~=0);
ind_A5FdU_TTG_6p8_C2   = ind_A5FdU_TTG_6p8_C2(ind_A5FdU_TTG_6p8_C2~=0);
ind_A5FdU_TTG_6p8_C1p  = ind_A5FdU_TTG_6p8_C1p(ind_A5FdU_TTG_6p8_C1p~=0);
ind_A5FdU_TTG_6p8_C4p  = ind_A5FdU_TTG_6p8_C4p(ind_A5FdU_TTG_6p8_C4p~=0);
ind_A5FdU_TTG_6p8_H6H8 = ind_A5FdU_TTG_6p8_H6H8(ind_A5FdU_TTG_6p8_H6H8~=0);
ind_A5FdU_TTG_6p8_H2   = ind_A5FdU_TTG_6p8_H2(ind_A5FdU_TTG_6p8_H2~=0);
ind_A5FdU_TTG_6p8_H1p  = ind_A5FdU_TTG_6p8_H1p(ind_A5FdU_TTG_6p8_H1p~=0);
ind_A5FdU_TTG_6p8_H4p  = ind_A5FdU_TTG_6p8_H4p(ind_A5FdU_TTG_6p8_H4p~=0);

A5FdU_TTG_6p8_HSQC_C6C8_list = A5FdU_TTG_6p8_HSQC_list(ind_A5FdU_TTG_6p8_C6C8,1:4);
A5FdU_TTG_6p8_HSQC_C2_list   = A5FdU_TTG_6p8_HSQC_list(ind_A5FdU_TTG_6p8_C2,1:4);
A5FdU_TTG_6p8_HSQC_C1p_list  = A5FdU_TTG_6p8_HSQC_list(ind_A5FdU_TTG_6p8_C1p,1:4);
A5FdU_TTG_6p8_HSQC_C4p_list  = A5FdU_TTG_6p8_HSQC_list(ind_A5FdU_TTG_6p8_C4p,1:4);
A5FdU_TTG_6p8_HSQC_H6H8_list = A5FdU_TTG_6p8_HSQC_list(ind_A5FdU_TTG_6p8_H6H8,1:4);
A5FdU_TTG_6p8_HSQC_H2_list   = A5FdU_TTG_6p8_HSQC_list(ind_A5FdU_TTG_6p8_H2,1:4);
A5FdU_TTG_6p8_HSQC_H1p_list  = A5FdU_TTG_6p8_HSQC_list(ind_A5FdU_TTG_6p8_H1p,1:4);
A5FdU_TTG_6p8_HSQC_H4p_list  = A5FdU_TTG_6p8_HSQC_list(ind_A5FdU_TTG_6p8_H4p,1:4);

A5FdU_TTG_6p8_resname_C6C8   = A5FdU_TTG_6p8_HSQC_C6C8_list.Group;
A5FdU_TTG_6p8_resnum_C6C8    = zeros(height(A5FdU_TTG_6p8_HSQC_C6C8_list),1);

for i = 1:height(A5FdU_TTG_6p8_HSQC_C6C8_list)

    A5FdU_TTG_6p8_resnum_C6C8(i,1) = str2double(regexp(A5FdU_TTG_6p8_resname_C6C8{i,1}, '(+|-)?(\d+)?\.?\d*', 'match'));

end

A5FdU_TTG_6p8_U5_C6_folded_shift = A5FdU_TTG_6p8_HSQC_C6C8_list(ismember(A5FdU_TTG_6p8_resname_C6C8,'U5'),:).Shift;
A5FdU_TTG_6p8_U5_C6_new_shift    = A5FdU_TTG_6p8_U5_C6_folded_shift - sw_HC_HSQC_700_UNC;
A5FdU_TTG_6p8_HSQC_C6C8_list(ismember(A5FdU_TTG_6p8_resname_C6C8,'U5'),:).Shift = A5FdU_TTG_6p8_U5_C6_new_shift;


A5FdU_TTG_6p8_resname_C2     = A5FdU_TTG_6p8_HSQC_C2_list.Group;
A5FdU_TTG_6p8_resnum_C2      = zeros(height(A5FdU_TTG_6p8_HSQC_C2_list),1);

for i = 1:height(A5FdU_TTG_6p8_HSQC_C2_list)

    A5FdU_TTG_6p8_resnum_C2(i,1) = str2double(regexp(A5FdU_TTG_6p8_resname_C2{i,1}, '(+|-)?(\d+)?\.?\d*', 'match'));

end

A5FdU_TTG_6p8_resname_C1p    = A5FdU_TTG_6p8_HSQC_C1p_list.Group;
A5FdU_TTG_6p8_resnum_C1p     = zeros(height(A5FdU_TTG_6p8_HSQC_C1p_list),1);

for i = 1:height(A5FdU_TTG_6p8_HSQC_C1p_list)

    A5FdU_TTG_6p8_resnum_C1p(i,1) = str2double(regexp(A5FdU_TTG_6p8_resname_C1p{i,1}, '(+|-)?(\d+)?\.?\d*', 'match'));

end

A5FdU_TTG_6p8_resname_C4p    = A5FdU_TTG_6p8_HSQC_C4p_list.Group;
A5FdU_TTG_6p8_resnum_C4p     = zeros(height(A5FdU_TTG_6p8_HSQC_C4p_list),1);

for i = 1:height(A5FdU_TTG_6p8_HSQC_C4p_list)

    A5FdU_TTG_6p8_resnum_C4p(i,1) = str2double(regexp(A5FdU_TTG_6p8_resname_C4p{i,1}, '(+|-)?(\d+)?\.?\d*', 'match'));

end


%% isoG5FdU TTG, pH 6.8

filepath_isoG5FdU_TTG_6p8  = sprintf('%s/%s/%s',parent_path,TTG_pathname,filename_isoG5FdU_TTG_pH6p8);

isoG5FdU_TTG_6p8_HSQC_list = readtable(filepath_isoG5FdU_TTG_6p8);

Assign_isoG5FdU_TTG_6p8    = isoG5FdU_TTG_6p8_HSQC_list.Atom;

ind_isoG5FdU_TTG_6p8_C6C8  = zeros(length(Assign_isoG5FdU_TTG_6p8),1);
ind_isoG5FdU_TTG_6p8_C2    = zeros(length(Assign_isoG5FdU_TTG_6p8),1);
ind_isoG5FdU_TTG_6p8_C1p   = zeros(length(Assign_isoG5FdU_TTG_6p8),1);
ind_isoG5FdU_TTG_6p8_C4p   = zeros(length(Assign_isoG5FdU_TTG_6p8),1);
ind_isoG5FdU_TTG_6p8_H6H8  = zeros(length(Assign_isoG5FdU_TTG_6p8),1);
ind_isoG5FdU_TTG_6p8_H2    = zeros(length(Assign_isoG5FdU_TTG_6p8),1);
ind_isoG5FdU_TTG_6p8_H1p   = zeros(length(Assign_isoG5FdU_TTG_6p8),1);
ind_isoG5FdU_TTG_6p8_H4p   = zeros(length(Assign_isoG5FdU_TTG_6p8),1);

for i = 1:length(Assign_isoG5FdU_TTG_6p8)

    if ismember('C6',Assign_isoG5FdU_TTG_6p8(i,1)) || ismember('C8',Assign_isoG5FdU_TTG_6p8(i,1))

        ind_isoG5FdU_TTG_6p8_C6C8(i,1) = i;

    elseif ismember('C2',Assign_isoG5FdU_TTG_6p8(i,1))

        ind_isoG5FdU_TTG_6p8_C2(i,1)   = i;

    elseif ismember("C1'",Assign_isoG5FdU_TTG_6p8(i,1))

        ind_isoG5FdU_TTG_6p8_C1p(i,1)   = i;

    elseif ismember("C4'",Assign_isoG5FdU_TTG_6p8(i,1))

        ind_isoG5FdU_TTG_6p8_C4p(i,1)   = i;

    elseif ismember('H6',Assign_isoG5FdU_TTG_6p8(i,1)) || ismember('H8',Assign_isoG5FdU_TTG_6p8(i,1))

        ind_isoG5FdU_TTG_6p8_H6H8(i,1)   = i;

    elseif ismember('H2',Assign_isoG5FdU_TTG_6p8(i,1))

        ind_isoG5FdU_TTG_6p8_H2(i,1)   = i;

    elseif ismember("H1'",Assign_isoG5FdU_TTG_6p8(i,1))

        ind_isoG5FdU_TTG_6p8_H1p(i,1)   = i;

    elseif ismember("H4'",Assign_isoG5FdU_TTG_6p8(i,1))

        ind_isoG5FdU_TTG_6p8_H4p(i,1)   = i;

    end

end

ind_isoG5FdU_TTG_6p8_C6C8 = ind_isoG5FdU_TTG_6p8_C6C8(ind_isoG5FdU_TTG_6p8_C6C8~=0);
ind_isoG5FdU_TTG_6p8_C2   = ind_isoG5FdU_TTG_6p8_C2(ind_isoG5FdU_TTG_6p8_C2~=0);
ind_isoG5FdU_TTG_6p8_C1p  = ind_isoG5FdU_TTG_6p8_C1p(ind_isoG5FdU_TTG_6p8_C1p~=0);
ind_isoG5FdU_TTG_6p8_C4p  = ind_isoG5FdU_TTG_6p8_C4p(ind_isoG5FdU_TTG_6p8_C4p~=0);
ind_isoG5FdU_TTG_6p8_H6H8 = ind_isoG5FdU_TTG_6p8_H6H8(ind_isoG5FdU_TTG_6p8_H6H8~=0);
ind_isoG5FdU_TTG_6p8_H2   = ind_isoG5FdU_TTG_6p8_H2(ind_isoG5FdU_TTG_6p8_H2~=0);
ind_isoG5FdU_TTG_6p8_H1p  = ind_isoG5FdU_TTG_6p8_H1p(ind_isoG5FdU_TTG_6p8_H1p~=0);
ind_isoG5FdU_TTG_6p8_H4p  = ind_isoG5FdU_TTG_6p8_H4p(ind_isoG5FdU_TTG_6p8_H4p~=0);

isoG5FdU_TTG_6p8_HSQC_C6C8_list = isoG5FdU_TTG_6p8_HSQC_list(ind_isoG5FdU_TTG_6p8_C6C8,1:4);
isoG5FdU_TTG_6p8_HSQC_C2_list   = isoG5FdU_TTG_6p8_HSQC_list(ind_isoG5FdU_TTG_6p8_C2,1:4);
isoG5FdU_TTG_6p8_HSQC_C1p_list  = isoG5FdU_TTG_6p8_HSQC_list(ind_isoG5FdU_TTG_6p8_C1p,1:4);
isoG5FdU_TTG_6p8_HSQC_C4p_list  = isoG5FdU_TTG_6p8_HSQC_list(ind_isoG5FdU_TTG_6p8_C4p,1:4);
isoG5FdU_TTG_6p8_HSQC_H6H8_list = isoG5FdU_TTG_6p8_HSQC_list(ind_isoG5FdU_TTG_6p8_H6H8,1:4);
isoG5FdU_TTG_6p8_HSQC_H2_list   = isoG5FdU_TTG_6p8_HSQC_list(ind_isoG5FdU_TTG_6p8_H2,1:4);
isoG5FdU_TTG_6p8_HSQC_H1p_list  = isoG5FdU_TTG_6p8_HSQC_list(ind_isoG5FdU_TTG_6p8_H1p,1:4);
isoG5FdU_TTG_6p8_HSQC_H4p_list  = isoG5FdU_TTG_6p8_HSQC_list(ind_isoG5FdU_TTG_6p8_H4p,1:4);

isoG5FdU_TTG_6p8_resname_C6C8   = isoG5FdU_TTG_6p8_HSQC_C6C8_list.Group;
isoG5FdU_TTG_6p8_resnum_C6C8    = zeros(height(isoG5FdU_TTG_6p8_HSQC_C6C8_list),1);

for i = 1:height(isoG5FdU_TTG_6p8_HSQC_C6C8_list)

    isoG5FdU_TTG_6p8_resnum_C6C8(i,1) = str2double(regexp(isoG5FdU_TTG_6p8_resname_C6C8{i,1}, '(+|-)?(\d+)?\.?\d*', 'match'));

end

isoG5FdU_TTG_6p8_U5_C6_folded_shift = isoG5FdU_TTG_6p8_HSQC_C6C8_list(ismember(isoG5FdU_TTG_6p8_resname_C6C8,'U5'),:).Shift;
isoG5FdU_TTG_6p8_U5_C6_new_shift    = isoG5FdU_TTG_6p8_U5_C6_folded_shift - sw_HC_HSQC_600_Duke;
isoG5FdU_TTG_6p8_HSQC_C6C8_list(ismember(isoG5FdU_TTG_6p8_resname_C6C8,'U5'),:).Shift = isoG5FdU_TTG_6p8_U5_C6_new_shift;

isoG5FdU_TTG_6p8_resname_C2     = isoG5FdU_TTG_6p8_HSQC_C2_list.Group;
isoG5FdU_TTG_6p8_resnum_C2      = zeros(height(isoG5FdU_TTG_6p8_HSQC_C2_list),1);

for i = 1:height(isoG5FdU_TTG_6p8_HSQC_C2_list)

    isoG5FdU_TTG_6p8_resnum_C2(i,1) = str2double(regexp(isoG5FdU_TTG_6p8_resname_C2{i,1}, '(+|-)?(\d+)?\.?\d*', 'match'));

end

isoG5FdU_TTG_6p8_resname_C1p    = isoG5FdU_TTG_6p8_HSQC_C1p_list.Group;
isoG5FdU_TTG_6p8_resnum_C1p     = zeros(height(isoG5FdU_TTG_6p8_HSQC_C1p_list),1);

for i = 1:height(isoG5FdU_TTG_6p8_HSQC_C1p_list)

    isoG5FdU_TTG_6p8_resnum_C1p(i,1) = str2double(regexp(isoG5FdU_TTG_6p8_resname_C1p{i,1}, '(+|-)?(\d+)?\.?\d*', 'match'));

end


isoG5FdU_TTG_6p8_resname_C4p    = isoG5FdU_TTG_6p8_HSQC_C4p_list.Group;
isoG5FdU_TTG_6p8_resnum_C4p     = zeros(height(isoG5FdU_TTG_6p8_HSQC_C4p_list),1);

for i = 1:height(isoG5FdU_TTG_6p8_HSQC_C4p_list)

    isoG5FdU_TTG_6p8_resnum_C4p(i,1) = str2double(regexp(isoG5FdU_TTG_6p8_resname_C4p{i,1}, '(+|-)?(\d+)?\.?\d*', 'match'));

end

%% compare chemical shifts, C6C8

delta_C6C8_GC_WC_6p8_vs_G5FdU_6p8 = GC_WC_TTG_6p8_HSQC_C6C8_list.Shift - G5FdU_TTG_6p8_HSQC_C6C8_list.Shift;
delta_C6C8_A5FdU_vs_G5FdU_6p8      = A5FdU_TTG_6p8_HSQC_C6C8_list.Shift - G5FdU_TTG_6p8_HSQC_C6C8_list.Shift;
delta_C6C8_G5FdU_10p1_vs_G5FdU_6p8 = G5FdU_TTG_10p1_HSQC_C6C8_list.Shift - G5FdU_TTG_6p8_HSQC_C6C8_list.Shift;
delta_C6C8_isoG5FdU_vs_G5FdU_6p8   = isoG5FdU_TTG_6p8_HSQC_C6C8_list.Shift - G5FdU_TTG_6p8_HSQC_C6C8_list.Shift;


%% compare chemical shifts, C2

delta_C2_GC_WC_6p8_vs_G5FdU_6p8 = GC_WC_TTG_6p8_HSQC_C2_list.Shift - G5FdU_TTG_6p8_HSQC_C2_list.Shift;
delta_C2_A5FdU_vs_G5FdU_6p8      = A5FdU_TTG_6p8_HSQC_C2_list.Shift(1:2) - G5FdU_TTG_6p8_HSQC_C2_list.Shift(1:2);
delta_C2_G5FdU_10p1_vs_G5FdU_6p8 = G5FdU_TTG_10p1_HSQC_C2_list.Shift - G5FdU_TTG_6p8_HSQC_C2_list.Shift;
delta_C2_isoG5FdU_vs_G5FdU_6p8   = isoG5FdU_TTG_6p8_HSQC_C2_list.Shift - G5FdU_TTG_6p8_HSQC_C2_list.Shift;



%% compare chemical shifts, C1'

delta_C1p_GC_WC_6p8_vs_G5FdU_6p8 = GC_WC_TTG_6p8_HSQC_C1p_list.Shift - G5FdU_TTG_6p8_HSQC_C1p_list.Shift;
delta_C1p_A5FdU_vs_G5FdU_6p8      = A5FdU_TTG_6p8_HSQC_C1p_list.Shift - G5FdU_TTG_6p8_HSQC_C1p_list.Shift;
delta_C1p_G5FdU_10p1_vs_G5FdU_6p8 = G5FdU_TTG_10p1_HSQC_C1p_list.Shift - G5FdU_TTG_6p8_HSQC_C1p_list.Shift;
delta_C1p_isoG5FdU_vs_G5FdU_6p8   = isoG5FdU_TTG_6p8_HSQC_C1p_list.Shift - G5FdU_TTG_6p8_HSQC_C1p_list.Shift;

%% compare chemical shifts, C4'

delta_C4p_GC_WC_6p8_vs_G5FdU_6p8 = GC_WC_TTG_6p8_HSQC_C4p_list.Shift - G5FdU_TTG_6p8_HSQC_C4p_list.Shift;
delta_C4p_A5FdU_vs_G5FdU_6p8      = A5FdU_TTG_6p8_HSQC_C4p_list.Shift - G5FdU_TTG_6p8_HSQC_C4p_list.Shift;
delta_C4p_G5FdU_10p1_vs_G5FdU_6p8 = G5FdU_TTG_10p1_HSQC_C4p_list.Shift - G5FdU_TTG_6p8_HSQC_C4p_list.Shift;
delta_C4p_isoG5FdU_vs_G5FdU_6p8   = isoG5FdU_TTG_6p8_HSQC_C4p_list.Shift - G5FdU_TTG_6p8_HSQC_C4p_list.Shift;

%% compare chemical shifts, H6H8

delta_H6H8_GC_WC_6p8_vs_G5FdU_6p8 = GC_WC_TTG_6p8_HSQC_H6H8_list.Shift - G5FdU_TTG_6p8_HSQC_H6H8_list.Shift;
delta_H6H8_A5FdU_vs_G5FdU_6p8      = A5FdU_TTG_6p8_HSQC_H6H8_list.Shift - G5FdU_TTG_6p8_HSQC_H6H8_list.Shift;
delta_H6H8_G5FdU_10p1_vs_G5FdU_6p8 = G5FdU_TTG_10p1_HSQC_H6H8_list.Shift - G5FdU_TTG_6p8_HSQC_H6H8_list.Shift;
delta_H6H8_isoG5FdU_vs_G5FdU_6p8   = isoG5FdU_TTG_6p8_HSQC_H6H8_list.Shift - G5FdU_TTG_6p8_HSQC_H6H8_list.Shift;

%% compare chemical shifts, H2

delta_H2_GC_WC_6p8_vs_G5FdU_6p8 = GC_WC_TTG_6p8_HSQC_H2_list.Shift - G5FdU_TTG_6p8_HSQC_H2_list.Shift;
delta_H2_A5FdU_vs_G5FdU_6p8      = A5FdU_TTG_6p8_HSQC_H2_list.Shift(1:2) - G5FdU_TTG_6p8_HSQC_H2_list.Shift(1:2);
delta_H2_G5FdU_10p1_vs_G5FdU_6p8 = G5FdU_TTG_10p1_HSQC_H2_list.Shift - G5FdU_TTG_6p8_HSQC_H2_list.Shift;
delta_H2_isoG5FdU_vs_G5FdU_6p8   = isoG5FdU_TTG_6p8_HSQC_H2_list.Shift - G5FdU_TTG_6p8_HSQC_H2_list.Shift;


%% compare chemical shifts, H1'

delta_H1p_GC_WC_6p8_vs_G5FdU_6p8 = GC_WC_TTG_6p8_HSQC_H1p_list.Shift - G5FdU_TTG_6p8_HSQC_H1p_list.Shift;
delta_H1p_A5FdU_vs_G5FdU_6p8      = A5FdU_TTG_6p8_HSQC_H1p_list.Shift - G5FdU_TTG_6p8_HSQC_H1p_list.Shift;
delta_H1p_G5FdU_10p1_vs_G5FdU_6p8 = G5FdU_TTG_10p1_HSQC_H1p_list.Shift - G5FdU_TTG_6p8_HSQC_H1p_list.Shift;
delta_H1p_isoG5FdU_vs_G5FdU_6p8   = isoG5FdU_TTG_6p8_HSQC_H1p_list.Shift - G5FdU_TTG_6p8_HSQC_H1p_list.Shift;

%% compare chemical shifts, H4'

delta_H4p_GC_WC_6p8_vs_G5FdU_6p8 = GC_WC_TTG_6p8_HSQC_H4p_list.Shift - G5FdU_TTG_6p8_HSQC_H4p_list.Shift;
delta_H4p_A5FdU_vs_G5FdU_6p8      = A5FdU_TTG_6p8_HSQC_H4p_list.Shift - G5FdU_TTG_6p8_HSQC_H4p_list.Shift;
delta_H4p_G5FdU_10p1_vs_G5FdU_6p8 = G5FdU_TTG_10p1_HSQC_H4p_list.Shift - G5FdU_TTG_6p8_HSQC_H4p_list.Shift;
delta_H4p_isoG5FdU_vs_G5FdU_6p8   = isoG5FdU_TTG_6p8_HSQC_H4p_list.Shift - G5FdU_TTG_6p8_HSQC_H4p_list.Shift;


%% G, 5FdU/T/C, and immidiate neighbors

neigh_resnum = [(4:6)';(14:16)'];

TTG_neigh_G5FdU_cs(:,1) = neigh_resnum;
TTG_neigh_G5FdU_cs(:,2) = G5FdU_TTG_6p8_HSQC_C6C8_list.Shift([4:6,13:15],1);
TTG_neigh_G5FdU_cs(:,3) = G5FdU_TTG_6p8_HSQC_C1p_list.Shift([4:6,11:13],1);
TTG_neigh_G5FdU_cs(:,4) = G5FdU_TTG_6p8_HSQC_C4p_list.Shift([4:6,12:14],1);
TTG_neigh_G5FdU_cs(:,5) = G5FdU_TTG_6p8_HSQC_H6H8_list.Shift([4:6,13:15],1);
TTG_neigh_G5FdU_cs(:,6) = G5FdU_TTG_6p8_HSQC_H1p_list.Shift([4:6,11:13],1);
TTG_neigh_G5FdU_cs(:,7) = G5FdU_TTG_6p8_HSQC_H4p_list.Shift([4:6,12:14],1);

TTG_neigh_delta_cs        = cell(5,7);
TTG_neigh_delta_cs(1,:)   = {'delta','C6C8','C1p','C4p','H6H8','H1','H4p'};
TTG_neigh_delta_cs(2:5,1) = {'high_pH';'isoG';'GC_WC';'UA_WC'};

TTG_neigh_delta_cs{2,2}(:,1) = neigh_resnum;
TTG_neigh_delta_cs{2,2}(:,2) = delta_C6C8_G5FdU_10p1_vs_G5FdU_6p8([4:6,13:15],1);
TTG_neigh_delta_cs{3,2}(:,1) = neigh_resnum;
TTG_neigh_delta_cs{3,2}(:,2) = delta_C6C8_isoG5FdU_vs_G5FdU_6p8([4:6,13:15],1);
TTG_neigh_delta_cs{4,2}(:,1) = neigh_resnum;
TTG_neigh_delta_cs{4,2}(:,2) = delta_C6C8_GC_WC_6p8_vs_G5FdU_6p8([4:6,13:15],1);
TTG_neigh_delta_cs{5,2}(:,1) = neigh_resnum;
TTG_neigh_delta_cs{5,2}(:,2) = delta_C6C8_A5FdU_vs_G5FdU_6p8([4:6,13:15],1);

TTG_neigh_delta_cs{2,3}(:,1) = neigh_resnum;
TTG_neigh_delta_cs{2,3}(:,2) = delta_C1p_G5FdU_10p1_vs_G5FdU_6p8([4:6,11:13],1);
TTG_neigh_delta_cs{3,3}(:,1) = neigh_resnum;
TTG_neigh_delta_cs{3,3}(:,2) = delta_C1p_isoG5FdU_vs_G5FdU_6p8([4:6,11:13],1);
TTG_neigh_delta_cs{4,3}(:,1) = neigh_resnum;
TTG_neigh_delta_cs{4,3}(:,2) = delta_C1p_GC_WC_6p8_vs_G5FdU_6p8([4:6,11:13],1);
TTG_neigh_delta_cs{5,3}(:,1) = neigh_resnum;
TTG_neigh_delta_cs{5,3}(:,2) = delta_C1p_A5FdU_vs_G5FdU_6p8([4:6,11:13],1);

TTG_neigh_delta_cs{2,4}(:,1) = neigh_resnum;
TTG_neigh_delta_cs{2,4}(:,2) = delta_C4p_G5FdU_10p1_vs_G5FdU_6p8([4:6,12:14],1);
TTG_neigh_delta_cs{3,4}(:,1) = neigh_resnum;
TTG_neigh_delta_cs{3,4}(:,2) = delta_C4p_isoG5FdU_vs_G5FdU_6p8([4:6,12:14],1);
TTG_neigh_delta_cs{4,4}(:,1) = neigh_resnum;
TTG_neigh_delta_cs{4,4}(:,2) = delta_C4p_GC_WC_6p8_vs_G5FdU_6p8([4:6,12:14],1);
TTG_neigh_delta_cs{5,4}(:,1) = neigh_resnum;
TTG_neigh_delta_cs{5,4}(:,2) = delta_C4p_A5FdU_vs_G5FdU_6p8([4:6,12:14],1);

TTG_neigh_delta_cs{2,5}(:,1) = neigh_resnum;
TTG_neigh_delta_cs{2,5}(:,2) = delta_H6H8_G5FdU_10p1_vs_G5FdU_6p8([4:6,13:15],1);
TTG_neigh_delta_cs{3,5}(:,1) = neigh_resnum;
TTG_neigh_delta_cs{3,5}(:,2) = delta_H6H8_isoG5FdU_vs_G5FdU_6p8([4:6,13:15],1);
TTG_neigh_delta_cs{4,5}(:,1) = neigh_resnum;
TTG_neigh_delta_cs{4,5}(:,2) = delta_H6H8_GC_WC_6p8_vs_G5FdU_6p8([4:6,13:15],1);
TTG_neigh_delta_cs{5,5}(:,1) = neigh_resnum;
TTG_neigh_delta_cs{5,5}(:,2) = delta_H6H8_A5FdU_vs_G5FdU_6p8([4:6,13:15],1);

TTG_neigh_delta_cs{2,6}(:,1) = neigh_resnum;
TTG_neigh_delta_cs{2,6}(:,2) = delta_H1p_G5FdU_10p1_vs_G5FdU_6p8([4:6,11:13],1);
TTG_neigh_delta_cs{3,6}(:,1) = neigh_resnum;
TTG_neigh_delta_cs{3,6}(:,2) = delta_H1p_isoG5FdU_vs_G5FdU_6p8([4:6,11:13],1);
TTG_neigh_delta_cs{4,6}(:,1) = neigh_resnum;
TTG_neigh_delta_cs{4,6}(:,2) = delta_H1p_GC_WC_6p8_vs_G5FdU_6p8([4:6,11:13],1);
TTG_neigh_delta_cs{5,6}(:,1) = neigh_resnum;
TTG_neigh_delta_cs{5,6}(:,2) = delta_H1p_A5FdU_vs_G5FdU_6p8([4:6,11:13],1);

TTG_neigh_delta_cs{2,7}(:,1) = neigh_resnum;
TTG_neigh_delta_cs{2,7}(:,2) = delta_H4p_G5FdU_10p1_vs_G5FdU_6p8([4:6,12:14],1);
TTG_neigh_delta_cs{3,7}(:,1) = neigh_resnum;
TTG_neigh_delta_cs{3,7}(:,2) = delta_H4p_isoG5FdU_vs_G5FdU_6p8([4:6,12:14],1);
TTG_neigh_delta_cs{4,7}(:,1) = neigh_resnum;
TTG_neigh_delta_cs{4,7}(:,2) = delta_H4p_GC_WC_6p8_vs_G5FdU_6p8([4:6,12:14],1);
TTG_neigh_delta_cs{5,7}(:,1) = neigh_resnum;
TTG_neigh_delta_cs{5,7}(:,2) = delta_H4p_A5FdU_vs_G5FdU_6p8([4:6,12:14],1);

save('TTG_CSPs.mat','TTG_neigh_delta_cs','TTG_neigh_G5FdU_cs');

