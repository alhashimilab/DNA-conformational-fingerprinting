clear;
clc;

%% colors and plot properties

x_min_1H = 9.0;
x_max_1H = 14.5;
y_min_1H = 0.5E8;
y_max_1H = 1.5E8;

x_min_19F = -166.0;
x_max_19F = -160.0; %-161.5
y_min_19F = -1E6;
y_max_19F = 1E7;


y_axis_min_pKa = 7.5;
y_axis_max_pKa = 10.0;

fontsize      = 28;
ticksize      = 0.015;
axeswidth     = 2;
linewidth     = 3;
markersize    = 18;
err_linewidth = 2;
err_capsize   = 10;
font_title    = 36;
markersize_cs = 19;
markersize_pKa = 22;
markersize_ddG = 22.5;

color_GT                = 'k';
color_G5FdU_iminos      = '#cc6677';
color_G5FdU_plots       = 'k';
color_G5FdU_Anion_plots = 'k';

color_G5FdU_lowpH = 'k';
color_G5FdU_medpH = '#b59a87';
color_G5FdU_Anion = '#ff9600';
color_G5FdU_ES    = '#44aa99';

color_A5FdU_WC    = '#5b8ea7';
color_GC_WC       = '#1f7d1f';
color_isoG5FdU    = '#5a2a23';

color_box_plot = 'k';
color_box_mean = '#808080';

color_1C  = '#003f5c';
color_5C  = '#444e86';
color_10C = '#955196';
color_15C = '#dd5182';
color_20C = '#ff6e54';
color_25C = '#ffa600';

colors_CTC = {color_1C;color_5C;color_10C;color_15C;color_20C};
colors_ATC = {color_1C;color_10C};


color_shade_CI = [220/255,220/255,220/255];


x_line = -6.0:0.05:6.0;
y_line = x_line;

conf_int_per = 0.95;

parent_path = '/Users/orsula1/OrsDocs/OrsDocs_since_Jan2025/2025-Sequence-dependence/2025-Data-Deposition';

%% Figure 2B + Figure S1 - 1D 1H imino spectra of G-T and G-5FdU hairpins

%% load 1D 1H imino spectra

path_1D_imino = 'NMR-Spectra/1D-1H-imino-spectra/1H-1D-G5FdU-vs-GT';
parent_path_1D_imino = sprintf('%s/%s',parent_path,path_1D_imino);

sequence_tri_1D_imino = {'GTC';'GTT';'GTA';'GTG';'ATC';'ATT';'ATA';'ATG';...
    'TTC';'TTT';'TTA';'TTG';'CTC';'CTT';'CTA';'CTG'};

mimics = {'GT';'G5FdU'};

F2_left_1D_imino  = 14.5;
F2_right_1D_imino = 9.0;

specs_1D_imino = cell(length(sequence_tri_1D_imino),3);

for i = 1:length(sequence_tri_1D_imino)

    specs_1D_imino{i,1} = sequence_tri_1D_imino{i};

    for j = 1:length(mimics)

        filename = sprintf('%s_%s_1D_imino_pH_Low.txt',sequence_tri_1D_imino{i},mimics{j});

        [axis_1H,spec_1H] = load_1D_spec(parent_path_1D_imino,filename,F2_left_1D_imino,F2_right_1D_imino);

        if ~isempty(axis_1H)

            specs_1D_imino{i,j+1} = [axis_1H,spec_1H];

        end

    end

end

specs_1D_imino_tbl = cell2table(specs_1D_imino,"VariableNames",["seq_T" "GT" "G5FdU"]);


%% Figure S1 - Plotting the data - 1D 1H imino, offset vertically

%% 1D 1H imino, unmodified GT

vert_offset = 2E8;
plot_ind = 0;

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S1A 1",NumberTitle="off");
set(gcf,'Color','w');

for i = 1:4

    plot(specs_1D_imino{i,2}(:,1),(specs_1D_imino{i,2}(:,2))-vert_offset*i,'Color',color_GT,'LineWidth',linewidth);
    hold on;
    text(9.8,1.0E8-vert_offset*i,sprintf('%s',sequence_tri_1D_imino{i}),'FontSize',fontsize,'FontWeight','bold');

end

hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([x_min_1H x_max_1H (-vert_offset*(i+1)+y_max_1H) (-vert_offset*plot_ind+y_min_1H)]);

propedit;

plot_ind = i;

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S1A 2",NumberTitle="off");
set(gcf,'Color','w');

for i = 5:8

    plot(specs_1D_imino{i,2}(:,1),(specs_1D_imino{i,2}(:,2))-vert_offset*i,'Color',color_GT,'LineWidth',linewidth);
    hold on;
    text(9.8,1.0E8-vert_offset*i,sprintf('%s',sequence_tri_1D_imino{i}),'FontSize',fontsize,'FontWeight','bold');

end

hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([x_min_1H x_max_1H (-vert_offset*(i+1)+y_max_1H) (-vert_offset*plot_ind+y_min_1H)]);

propedit;

plot_ind = i;

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S1A 3",NumberTitle="off");
set(gcf,'Color','w');

for i = 9:12

    plot(specs_1D_imino{i,2}(:,1),(specs_1D_imino{i,2}(:,2))-vert_offset*i,'Color',color_GT,'LineWidth',linewidth);
    hold on;
    text(9.8,1.0E8-vert_offset*i,sprintf('%s',sequence_tri_1D_imino{i}),'FontSize',fontsize,'FontWeight','bold');

end

hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([x_min_1H x_max_1H (-vert_offset*(i+1)+y_max_1H) (-vert_offset*plot_ind+y_min_1H)]);

propedit;

plot_ind = i;

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S1A 4",NumberTitle="off");
set(gcf,'Color','w');

for i = 13:16

    plot(specs_1D_imino{i,2}(:,1),(specs_1D_imino{i,2}(:,2))-vert_offset*i,'Color',color_GT,'LineWidth',linewidth);
    hold on;
    text(9.8,1.0E8-vert_offset*i,sprintf('%s',sequence_tri_1D_imino{i}),'FontSize',fontsize,'FontWeight','bold');

end

hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([x_min_1H x_max_1H (-vert_offset*(i+1)+y_max_1H) (-vert_offset*plot_ind+y_min_1H)]);

propedit;


%% 1D 1H imino, modified G5FdU

scale_G5FdU = 4;
plot_ind = 0;

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S1B 1",NumberTitle="off");
set(gcf,'Color','w');

for i = 1:4

    if i == 2

        plot(specs_1D_imino{i,3}(:,1),(specs_1D_imino{i,3}(:,2)).*scale_G5FdU*2.3-vert_offset*i,...
            'Color',color_G5FdU_iminos,'LineWidth',linewidth);
        hold on;

    elseif i == 4

        plot(specs_1D_imino{i,3}(:,1),(specs_1D_imino{i,3}(:,2)).*scale_G5FdU*1.2-vert_offset*i,...
            'Color',color_G5FdU_iminos,'LineWidth',linewidth);
        hold on;

    else

        plot(specs_1D_imino{i,3}(:,1),(specs_1D_imino{i,3}(:,2)).*scale_G5FdU*1.0-vert_offset*i,...
            'Color',color_G5FdU_iminos,'LineWidth',linewidth);
        hold on;

    end

    text(9.8,1.0E8-vert_offset*i,sprintf('%s',sequence_tri_1D_imino{i}),'FontSize',fontsize,'FontWeight','bold');

end

hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([x_min_1H x_max_1H (-vert_offset*(i+1)+y_max_1H) (-vert_offset*plot_ind+y_min_1H)]);

propedit;

plot_ind = i;

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S1B 2",NumberTitle="off");
set(gcf,'Color','w');

for i = 5:8

    if i == 6

        plot(specs_1D_imino{i,3}(:,1),(specs_1D_imino{i,3}(:,2)).*scale_G5FdU*2.0-vert_offset*i,...
            'Color',color_G5FdU_iminos,'LineWidth',linewidth);
        hold on;

    else

        plot(specs_1D_imino{i,3}(:,1),(specs_1D_imino{i,3}(:,2)).*scale_G5FdU*1.2-vert_offset*i,...
            'Color',color_G5FdU_iminos,'LineWidth',linewidth);
        hold on;

    end

    text(9.8,1.0E8-vert_offset*i,sprintf('%s',sequence_tri_1D_imino{i}),'FontSize',fontsize,'FontWeight','bold');

end

hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([x_min_1H x_max_1H (-vert_offset*(i+1)+y_max_1H) (-vert_offset*plot_ind+y_min_1H)]);

propedit;

plot_ind = i;

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S1B 3",NumberTitle="off");
set(gcf,'Color','w');

for i = 9:12

    if i == 9

        plot(specs_1D_imino{i,3}(:,1),(specs_1D_imino{i,3}(:,2)).*scale_G5FdU*2.2-vert_offset*i,...
            'Color',color_G5FdU_iminos,'LineWidth',linewidth);
        hold on;

    elseif i == 11

        plot(specs_1D_imino{i,3}(:,1),(specs_1D_imino{i,3}(:,2)).*scale_G5FdU*1.5-vert_offset*i,...
            'Color',color_G5FdU_iminos,'LineWidth',linewidth);
        hold on;

    else

        plot(specs_1D_imino{i,3}(:,1),(specs_1D_imino{i,3}(:,2)).*scale_G5FdU*2.0-vert_offset*i,...
            'Color',color_G5FdU_iminos,'LineWidth',linewidth);
        hold on;

    end

    text(9.8,1.0E8-vert_offset*i,sprintf('%s',sequence_tri_1D_imino{i}),'FontSize',fontsize,'FontWeight','bold');

end

hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([x_min_1H x_max_1H (-vert_offset*(i+1)+y_max_1H) (-vert_offset*plot_ind+y_min_1H)]);

propedit;

plot_ind = i;

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S1B 4",NumberTitle="off");
set(gcf,'Color','w');

for i = 13:16

    if i == 13

        plot(specs_1D_imino{i,3}(:,1),(specs_1D_imino{i,3}(:,2)).*scale_G5FdU*0.9-vert_offset*i,...
            'Color',color_G5FdU_iminos,'LineWidth',linewidth);
        hold on;

    elseif i == 14

        plot(specs_1D_imino{i,3}(:,1),(specs_1D_imino{i,3}(:,2)).*scale_G5FdU*1.7-vert_offset*i,...
            'Color',color_G5FdU_iminos,'LineWidth',linewidth);
        hold on;

    else

        plot(specs_1D_imino{i,3}(:,1),(specs_1D_imino{i,3}(:,2)).*scale_G5FdU*1.2-vert_offset*i,...
            'Color',color_G5FdU_iminos,'LineWidth',linewidth);
        hold on;
    end

    text(9.8,1.0E8-vert_offset*i,sprintf('%s',sequence_tri_1D_imino{i}),'FontSize',fontsize,'FontWeight','bold');

end

hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([x_min_1H x_max_1H (-vert_offset*(i+1)+y_max_1H) (-vert_offset*plot_ind+y_min_1H)]);

propedit;

%% Figure 2B - 1D 1H imino for two sequence contexts, offset vertically

sel_seqs = {'ATG';'TTC'};

%% 1D 1H imino, unmodified GT

specs_1D_imino_top = specs_1D_imino(strcmp(sel_seqs{1,1},specs_1D_imino(:,1)),:);
specs_1D_imino_bot = specs_1D_imino(strcmp(sel_seqs{2,1},specs_1D_imino(:,1)),:);


figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S2B GT",NumberTitle="off");
set(gcf,'Color','w');

plot(specs_1D_imino_top{1,2}(:,1),(specs_1D_imino_top{1,2}(:,2))-vert_offset*1.0,'Color',color_GT,'LineWidth',linewidth);
hold on;
text(9.8,1.0E8-vert_offset*1.0,sprintf('%s',sel_seqs{1,1}),'FontSize',fontsize,'FontWeight','bold');

plot(specs_1D_imino_bot{1,2}(:,1),(specs_1D_imino_bot{1,2}(:,2))-vert_offset*2.0,'Color',color_GT,'LineWidth',linewidth);
text(9.8,1.0E8-vert_offset*2.0,sprintf('%s',sel_seqs{2,1}),'FontSize',fontsize,'FontWeight','bold');

hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([x_min_1H x_max_1H (-vert_offset*3.0+y_max_1H) y_min_1H]);

propedit;


figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S2B G5FdU",NumberTitle="off");
set(gcf,'Color','w');

plot(specs_1D_imino_top{1,3}(:,1),(specs_1D_imino_top{1,3}(:,2)).*scale_G5FdU*1.2-vert_offset*1.0,'Color',color_G5FdU_iminos,'LineWidth',linewidth);
hold on;
text(9.8,1.0E8-vert_offset*1.0,sprintf('%s',sel_seqs{1,1}),'FontSize',fontsize,'FontWeight','bold');

plot(specs_1D_imino_bot{1,3}(:,1),(specs_1D_imino_bot{1,3}(:,2)).*scale_G5FdU*2.2-vert_offset*2.0,'Color',color_G5FdU_iminos,'LineWidth',linewidth);
text(9.8,1.0E8-vert_offset*2.0,sprintf('%s',sel_seqs{2,1}),'FontSize',fontsize,'FontWeight','bold');

hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([x_min_1H x_max_1H (-vert_offset*3.0+y_max_1H) y_min_1H]);

propedit;

%% Figure 2C + Figure S3A + Figure S4A - Chemical shift heatmaps (1H, 19F)

%% loading parameters from 19F measurements

path_anion_pops     = 'Populations-and-Chemical-Shifts';
parent_path_pops_cs = sprintf('%s/%s',parent_path,path_anion_pops);
filename_pops_19F   = 'Anion_populations_19F_1C_all_pHs.xlsx';
filepath            = sprintf('%s/%s',parent_path_pops_cs,filename_pops_19F);

%% 19F-NMR-derived pKas

params_pH_pKas        = readtable(filepath);

pKas_19F              = params_pH_pKas.pKa_avg(~isnan(params_pH_pKas.pKa_avg));
pKas_19F_errs         = params_pH_pKas.pKa_avg_err(~isnan(params_pH_pKas.pKa_avg));
seqs_G_19F            = params_pH_pKas.Seq_G(~isnan(params_pH_pKas.pKa_avg));

[sorted_pKas_19F,ind] = sort(pKas_19F);
sorted_pKas_19F_errs  = pKas_19F_errs(ind);
sorted_seqs_G_19F     = seqs_G_19F(ind);

sorted_seqs_T_19F     = convSeqtoT(sorted_seqs_G_19F);

x_axis_seqs           = 1:length(sorted_seqs_T_19F);

ref_seq = 'CTA';

dpKa_T_N3_G5FdU_Anion_os_sorted      = sorted_pKas_19F - sorted_pKas_19F(ismember(sorted_seqs_T_19F,ref_seq));
dpKa_T_N3_G5FdU_Anion_os_sorted_errs = sqrt(sorted_pKas_19F_errs.^2 + sorted_pKas_19F_errs(ismember(sorted_seqs_T_19F,ref_seq)).^2);

%% 1H dw T(H3) in G-T wobble

filename_cs = '19F-and-1H-Chemical-Shifts.xlsx';
filepath    = sprintf('%s/%s',parent_path_pops_cs,filename_cs);

params_GT_1H_cs = readtable(filepath,'Sheet','GT 1H chemical shifts');

cs_T_H3_GT      = params_GT_1H_cs.T5_H3_ppm(~isnan(params_GT_1H_cs.T5_H3_ppm));
cs_T_H3_GT_err  = params_GT_1H_cs.T5_H3_err(~isnan(params_GT_1H_cs.T5_H3_ppm));
cs_G_H1_GT      = params_GT_1H_cs.G15_H1_ppm(~isnan(params_GT_1H_cs.G15_H1_ppm));
cs_G_H1_GT_err  = params_GT_1H_cs.G15_H1_err(~isnan(params_GT_1H_cs.G15_H1_err));

[~,ind_cs_GT] = ismember(sorted_seqs_G_19F,params_GT_1H_cs.Seq_G);

cs_T_H3_GT_sort_by_pKa      = cs_T_H3_GT(ind_cs_GT);
cs_T_H3_GT_sort_by_pKa_err  = cs_T_H3_GT_err(ind_cs_GT);
cs_G_H1_GT_sort_by_pKa      = cs_G_H1_GT(ind_cs_GT);
cs_G_H1_GT_sort_by_pKa_err  = cs_G_H1_GT_err(ind_cs_GT);

dw_T_H3_GT_WB_os_sort_by_pKa      = cs_T_H3_GT_sort_by_pKa - cs_T_H3_GT_sort_by_pKa(ismember(sorted_seqs_T_19F,ref_seq));
dw_T_H3_GT_WB_os_sort_by_pKa_errs = sqrt(cs_T_H3_GT_sort_by_pKa_err.^2 + cs_T_H3_GT_sort_by_pKa_err(ismember(sorted_seqs_T_19F,ref_seq)).^2);

dw_G_H1_GT_WB_os_sort_by_pKa      = cs_G_H1_GT_sort_by_pKa - cs_G_H1_GT_sort_by_pKa(ismember(sorted_seqs_T_19F,ref_seq));
dw_G_H1_GT_WB_os_sort_by_pKa_errs = sqrt(cs_G_H1_GT_sort_by_pKa_err.^2 + cs_G_H1_GT_sort_by_pKa_err(ismember(sorted_seqs_T_19F,ref_seq)).^2);

%% 1H dw 5FdU(H3) in G-5FdU wobble

params_G5FdU_1H_cs = readtable(filepath,'Sheet','G5FdU 1H chemical shifts');

cs_U_H3_G5FdU      = params_G5FdU_1H_cs.U5_H3_ppm(~isnan(params_G5FdU_1H_cs.U5_H3_ppm));
cs_U_H3_G5FdU_err  = params_G5FdU_1H_cs.U5_H3_err(~isnan(params_G5FdU_1H_cs.U5_H3_ppm));
cs_G_H1_G5FdU      = params_G5FdU_1H_cs.G15_H1_ppm(~isnan(params_G5FdU_1H_cs.G15_H1_ppm));
cs_G_H1_G5FdU_err  = params_G5FdU_1H_cs.G15_H1_err(~isnan(params_G5FdU_1H_cs.G15_H1_err));

[~,ind_cs_G5FdU_1H]   = ismember(sorted_seqs_G_19F,params_G5FdU_1H_cs.Seq_G);

cs_U_H3_G5FdU_sort_by_pKa      = cs_U_H3_G5FdU(ind_cs_G5FdU_1H);
cs_U_H3_G5FdU_sort_by_pKa_err  = cs_U_H3_G5FdU_err(ind_cs_G5FdU_1H);
cs_G_H1_G5FdU_sort_by_pKa      = cs_G_H1_G5FdU(ind_cs_G5FdU_1H);
cs_G_H1_G5FdU_sort_by_pKa_err  = cs_G_H1_G5FdU_err(ind_cs_G5FdU_1H);

dw_T_H3_G5FdU_WB_os_sort_by_pKa      = cs_U_H3_G5FdU_sort_by_pKa - cs_U_H3_G5FdU_sort_by_pKa(ismember(sorted_seqs_T_19F,ref_seq));
dw_T_H3_G5FdU_WB_os_sort_by_pKa_errs = sqrt(cs_U_H3_G5FdU_sort_by_pKa_err.^2 + cs_U_H3_G5FdU_sort_by_pKa_err(ismember(sorted_seqs_T_19F,ref_seq)).^2);

dw_G_H1_G5FdU_WB_os_sort_by_pKa      = cs_G_H1_G5FdU_sort_by_pKa - cs_G_H1_G5FdU_sort_by_pKa(ismember(sorted_seqs_T_19F,ref_seq));
dw_G_H1_G5FdU_WB_os_sort_by_pKa_errs = sqrt(cs_G_H1_G5FdU_sort_by_pKa_err.^2 + cs_G_H1_G5FdU_sort_by_pKa_err(ismember(sorted_seqs_T_19F,ref_seq)).^2);


%% 19F dw 5FdU(F5) in G-5FdU wobble

params_G5FdU_19F_cs = readtable(filepath,'Sheet','G5FdU 19F chemical shifts');

cs_F5_G5FdU_WB         = params_G5FdU_19F_cs.cs_WB_ppm(~isnan(params_G5FdU_19F_cs.cs_WB_ppm));
cs_F5_G5FdU_WB_err     = params_G5FdU_19F_cs.cs_WB_err_ppm(~isnan(params_G5FdU_19F_cs.cs_WB_ppm));

[~,ind_cs_G5FdU_19F]   = ismember(sorted_seqs_G_19F,params_G5FdU_19F_cs.Seq_G);

cs_F5_G5FdU_WB_sort_by_pKa         = cs_F5_G5FdU_WB(ind_cs_G5FdU_19F);
cs_F5_G5FdU_WB_sort_by_pKa_err     = cs_F5_G5FdU_WB_err(ind_cs_G5FdU_19F);

dw_T_F5_G5FdU_WB_os_sort_by_pKa = cs_F5_G5FdU_WB_sort_by_pKa - cs_F5_G5FdU_WB_sort_by_pKa(ismember(sorted_seqs_T_19F,ref_seq));
dw_T_F5_G5FdU_WB_os_sort_by_pKa_errs = sqrt(cs_F5_G5FdU_WB_sort_by_pKa_err.^2 + cs_F5_G5FdU_WB_sort_by_pKa_err(ismember(sorted_seqs_T_19F,ref_seq)).^2);

%% 19F dw 5FdU(F5) in G-5FdU anion

cs_F5_G5FdU_Anion      = params_G5FdU_19F_cs.cs_Anion_ppm(~isnan(params_G5FdU_19F_cs.cs_WB_ppm));
cs_F5_G5FdU_Anion_err  = params_G5FdU_19F_cs.cs_Anion_err_ppm(~isnan(params_G5FdU_19F_cs.cs_WB_ppm));

cs_F5_G5FdU_Anion_sort_by_pKa      = cs_F5_G5FdU_Anion(ind_cs_G5FdU_19F);
cs_F5_G5FdU_Anion_sort_by_pKa_err  = cs_F5_G5FdU_Anion_err(ind_cs_G5FdU_19F);

dw_T_F5_G5FdU_Anion_os_sort_by_pKa = cs_F5_G5FdU_Anion_sort_by_pKa - cs_F5_G5FdU_Anion_sort_by_pKa(ismember(sorted_seqs_T_19F,ref_seq));
dw_T_F5_G5FdU_Anion_os_sort_by_pKa_errs = sqrt(cs_F5_G5FdU_Anion_sort_by_pKa_err.^2 + cs_F5_G5FdU_Anion_sort_by_pKa_err(ismember(sorted_seqs_T_19F,ref_seq)).^2);


%% Heatmaps

seq_mat                 = zeros(length(sorted_seqs_G_19F),2);

dw_T_H3_GT_WB_os_sort_by_pKa_heatmap        = NaN(4);
dw_T_H3_G5FdU_WB_os_sort_by_pKa_heatmap     = NaN(4);

dw_G_H1_GT_WB_os_sort_by_pKa_heatmap        = NaN(4);
dw_G_H1_G5FdU_WB_os_sort_by_pKa_heatmap     = NaN(4);

dw_T_F5_G5FdU_WB_os_sort_by_pKa_heatmap     = NaN(4);
dw_T_F5_G5FdU_Anion_os_sort_by_pKa_heatmap  = NaN(4);

dpKa_T_N3_G5FdU_Anion_os_sorted_heatmap     = NaN(4);

for i = 1:length(seq_mat)

    if sorted_seqs_G_19F{i,1}(1) == 'C' % 3'-neighbor of dU
        seq_mat(i,1) = 1;
    elseif sorted_seqs_G_19F{i,1}(1) == 'T'
        seq_mat(i,1) = 2;
    elseif sorted_seqs_G_19F{i,1}(1) == 'A'
        seq_mat(i,1) = 3;
    elseif sorted_seqs_G_19F{i,1}(1) == 'G'
        seq_mat(i,1) = 4;
    end

    if sorted_seqs_G_19F{i,1}(3) == 'C' % 5'-neighbor of dU
        seq_mat(i,2) = 1;
    elseif sorted_seqs_G_19F{i,1}(3) == 'T'
        seq_mat(i,2) = 2;
    elseif sorted_seqs_G_19F{i,1}(3) == 'A'
        seq_mat(i,2) = 3;
    elseif sorted_seqs_G_19F{i,1}(3) == 'G'
        seq_mat(i,2) = 4;
    end

    dw_T_H3_GT_WB_os_sort_by_pKa_heatmap(seq_mat(i,1),seq_mat(i,2))       = dw_T_H3_GT_WB_os_sort_by_pKa(i,1);
    dw_T_H3_G5FdU_WB_os_sort_by_pKa_heatmap(seq_mat(i,1),seq_mat(i,2))    = dw_T_H3_G5FdU_WB_os_sort_by_pKa(i,1);

    dw_G_H1_GT_WB_os_sort_by_pKa_heatmap(seq_mat(i,1),seq_mat(i,2))       = dw_G_H1_GT_WB_os_sort_by_pKa(i,1);
    dw_G_H1_G5FdU_WB_os_sort_by_pKa_heatmap(seq_mat(i,1),seq_mat(i,2))    = dw_G_H1_G5FdU_WB_os_sort_by_pKa(i,1);

    dw_T_F5_G5FdU_WB_os_sort_by_pKa_heatmap(seq_mat(i,1),seq_mat(i,2))    = dw_T_F5_G5FdU_WB_os_sort_by_pKa(i,1);
    dw_T_F5_G5FdU_Anion_os_sort_by_pKa_heatmap(seq_mat(i,1),seq_mat(i,2)) = dw_T_F5_G5FdU_Anion_os_sort_by_pKa(i,1);

    dpKa_T_N3_G5FdU_Anion_os_sorted_heatmap(seq_mat(i,1),seq_mat(i,2))    = dpKa_T_N3_G5FdU_Anion_os_sorted(i,1);

end

dw_G_H1_GT_WB_os_sort_by_pKa_heatmap_by_G    = transpose(dw_G_H1_GT_WB_os_sort_by_pKa_heatmap);
dw_G_H1_G5FdU_WB_os_sort_by_pKa_heatmap_by_G = transpose(dw_G_H1_G5FdU_WB_os_sort_by_pKa_heatmap);

%% neighbors for heatmaps

seq_3p_dU = {'G';'A';'T';'C'}; % 1;2;3;4
seq_5p_dU = {'G';'A';'T';'C'}; % 4;3;2;1

seq_5p_G = {'G';'A';'T';'C'}; % 1;2;3;4
seq_3p_G = {'G';'A';'T';'C'}; % 4;3;2;1

%% find 3' and 5' neighbors of T (of 5FdU)

neighbor_3p_to_T = zeros(length(sorted_seqs_T_19F),1);
neighbor_5p_to_T = zeros(length(sorted_seqs_T_19F),1);

neighbor_3p_to_G = zeros(length(sorted_seqs_T_19F),1);
neighbor_5p_to_G = zeros(length(sorted_seqs_T_19F),1);

for i = 1:length(neighbor_3p_to_T)

    if sorted_seqs_T_19F{i,1}(3) == 'T' | sorted_seqs_T_19F{i,1}(3) == 'C'
        neighbor_3p_to_T(i,1) = 1; % 1 = pyrimidine, 2 = purine
        neighbor_5p_to_G(i,1) = 2;
    elseif sorted_seqs_T_19F{i,1}(3) == 'A' | sorted_seqs_T_19F{i,1}(3) == 'G'
        neighbor_3p_to_T(i,1) = 2;
        neighbor_5p_to_G(i,1) = 1;
    end

    if sorted_seqs_T_19F{i,1}(1) == 'T' | sorted_seqs_T_19F{i,1}(1) == 'C'
        neighbor_5p_to_T(i,1) = 1;
        neighbor_3p_to_G(i,1) = 2;
    elseif sorted_seqs_T_19F{i,1}(1) == 'A' | sorted_seqs_T_19F{i,1}(1) == 'G'
        neighbor_5p_to_T(i,1) = 2;
        neighbor_3p_to_G(i,1) = 1;
    end

end

[neighbor_5p_to_T_sorted,ind_5p]   = sort(neighbor_5p_to_T);
[neighbor_3p_to_G_sorted,ind_3p_G] = sort(neighbor_3p_to_G,'descend');

dpKa_T_N3_G5FdU_Anion_os_sort_by_5p = dpKa_T_N3_G5FdU_Anion_os_sorted(ind_5p);

dw_T_H3_GT_WB_os_sort_by_5p         = dw_T_H3_GT_WB_os_sort_by_pKa(ind_5p);
dw_G_H1_GT_WB_os_sort_by_3p_G       = dw_G_H1_GT_WB_os_sort_by_pKa(ind_3p_G);
dw_T_H3_G5FdU_WB_os_sort_by_5p      = dw_T_H3_G5FdU_WB_os_sort_by_pKa(ind_5p);
dw_G_H1_G5FdU_WB_os_sort_by_3p_G    = dw_G_H1_G5FdU_WB_os_sort_by_pKa(ind_3p_G);
dw_T_F5_G5FdU_WB_os_sort_by_5p      = dw_T_F5_G5FdU_WB_os_sort_by_pKa(ind_5p);
dw_T_F5_G5FdU_Anion_os_sort_by_5p   = dw_T_F5_G5FdU_Anion_os_sort_by_pKa(ind_5p);

%% Wilcoxon rank sum test (non-parametric test, equivalent to the Mann-Whitney U-test

[pWR_T_H3_GT_WB_3p,hWR_T_H3_GT_WB_3p] = ranksum(dw_T_H3_GT_WB_os_sort_by_pKa(1:8),dw_T_H3_GT_WB_os_sort_by_pKa(9:16));
[pWR_T_H3_GT_WB_5p,hWR_T_H3_GT_WB_5p] = ranksum(dw_T_H3_GT_WB_os_sort_by_5p(1:8),dw_T_H3_GT_WB_os_sort_by_5p(9:16));

[pWR_G_H1_GT_WB_5p_G,hWR_G_H1_GT_WB_5p_G] = ranksum(dw_G_H1_GT_WB_os_sort_by_pKa(1:8),dw_G_H1_GT_WB_os_sort_by_pKa(9:16));
[pWR_G_H1_GT_WB_3p_G,hWR_G_H1_GT_WB_3p_G] = ranksum(dw_G_H1_GT_WB_os_sort_by_3p_G(1:8),dw_G_H1_GT_WB_os_sort_by_3p_G(9:16));

[pWR_T_H3_G5FdU_WB_3p,hWR_T_H3_G5FdU_WB_3p] = ranksum(dw_T_H3_G5FdU_WB_os_sort_by_pKa(1:8),dw_T_H3_G5FdU_WB_os_sort_by_pKa(9:16));
[pWR_T_H3_G5FdU_WB_5p,hWR_T_H3_G5FdU_WB_5p] = ranksum(dw_T_H3_G5FdU_WB_os_sort_by_5p(1:8),dw_T_H3_G5FdU_WB_os_sort_by_5p(9:16));

[pWR_G_H1_G5FdU_WB_5p_G,hWR_G_H1_G5FdU_WB_5p_G] = ranksum(dw_G_H1_G5FdU_WB_os_sort_by_pKa(1:8),dw_G_H1_G5FdU_WB_os_sort_by_pKa(9:16));
[pWR_G_H1_G5FdU_WB_3p_G,hWR_G_H1_G5FdU_WB_3p_G] = ranksum(dw_G_H1_G5FdU_WB_os_sort_by_3p_G(1:8),dw_G_H1_G5FdU_WB_os_sort_by_3p_G(9:16));

[pWR_T_F5_G5FdU_WB_3p,hWR_T_F5_G5FdU_WB_3p] = ranksum(dw_T_F5_G5FdU_WB_os_sort_by_pKa(1:8),dw_T_F5_G5FdU_WB_os_sort_by_pKa(9:16));
[pWR_T_F5_G5FdU_WB_5p,hWR_T_F5_G5FdU_WB_5p] = ranksum(dw_T_F5_G5FdU_WB_os_sort_by_5p(1:8),dw_T_F5_G5FdU_WB_os_sort_by_5p(9:16));

[pWR_T_F5_G5FdU_Anion_3p,hWR_T_F5_G5FdU_Anion_3p] = ranksum(dw_T_F5_G5FdU_Anion_os_sort_by_pKa(1:8),dw_T_F5_G5FdU_Anion_os_sort_by_pKa(9:16));
[pWR_T_F5_G5FdU_Anion_5p,hWR_T_F5_G5FdU_Anion_5p] = ranksum(dw_T_F5_G5FdU_Anion_os_sort_by_5p(1:8),dw_T_F5_G5FdU_Anion_os_sort_by_5p(9:16));

[pWR_dpKa_T_N3_G5FdU_Anion_3p,hWR_dpKa_T_N3_G5FdU_Anion_3p] = ranksum(dpKa_T_N3_G5FdU_Anion_os_sorted(1:8),dpKa_T_N3_G5FdU_Anion_os_sorted(9:16));
[pWR_dpKa_T_N3_G5FdU_Anion_5p,hWR_dpKa_T_N3_G5FdU_Anion_5p] = ranksum(dpKa_T_N3_G5FdU_Anion_os_sort_by_5p(1:8),dpKa_T_N3_G5FdU_Anion_os_sort_by_5p(9:16));


%% Figure 2C, top - heatmap for dw T(H3) GT WB 

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig 2C top",NumberTitle="off");
set(gcf,'Color','w');

h3 = heatmap(seq_5p_dU,flip(seq_3p_dU),flip(dw_T_H3_GT_WB_os_sort_by_pKa_heatmap),'CellLabelColor','none');
h3.MissingDataColor = [0.8 0.8 0.8];
h3.Colormap = earthyjet;

set(gca, 'FontSize', fontsize);

propedit;

%% Figure 2C, bottom - heatmap for dw G(H1) GT WB

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig 2C bottom",NumberTitle="off");
set(gcf,'Color','w');

h3 = heatmap(seq_5p_G,flip(seq_3p_G),flip(dw_G_H1_GT_WB_os_sort_by_pKa_heatmap_by_G,2),'CellLabelColor','none');
h3.MissingDataColor = [0.8 0.8 0.8];
h3.Colormap = earthyjet;

set(gca, 'FontSize', fontsize);

propedit;

%% Figure S3A, top - heatmap for dw 5FdU(H3) G5FdU WB

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S3A top",NumberTitle="off");
set(gcf,'Color','w');

h3 = heatmap(seq_5p_dU,flip(seq_3p_dU),flip(dw_T_H3_G5FdU_WB_os_sort_by_pKa_heatmap),'CellLabelColor','none');
h3.MissingDataColor = [0.8 0.8 0.8];
h3.Colormap = earthyjet;

set(gca, 'FontSize', fontsize);

propedit;


%% Figure S3A, bottom - heatmap for dw G(H1) G5FdU WB

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S3A bottom",NumberTitle="off");
set(gcf,'Color','w');

h3 = heatmap(seq_5p_G,flip(seq_3p_G),flip(dw_G_H1_G5FdU_WB_os_sort_by_pKa_heatmap_by_G,2),'CellLabelColor','none');
h3.MissingDataColor = [0.8 0.8 0.8];
h3.Colormap = earthyjet;

set(gca, 'FontSize', fontsize);

propedit;


%% Figure S4A, left - heatmap for dw T(F5) G5FdU WB

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S4A left",NumberTitle="off");
set(gcf,'Color','w');

h3 = heatmap(seq_5p_dU,flip(seq_3p_dU),flip(dw_T_F5_G5FdU_WB_os_sort_by_pKa_heatmap),'CellLabelColor','none');
h3.MissingDataColor = [0.8 0.8 0.8];
h3.Colormap = earthyjet;

set(gca, 'FontSize', fontsize);

propedit;

%% Figure S4B, left - heatmap for dw T(F5) G5FdU Anion

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S4B left",NumberTitle="off");
set(gcf,'Color','w');

h3 = heatmap(seq_5p_dU,flip(seq_3p_dU),flip(dw_T_F5_G5FdU_Anion_os_sort_by_pKa_heatmap),'CellLabelColor','none');
h3.MissingDataColor = [0.8 0.8 0.8];
h3.Colormap = earthyjet;

set(gca, 'FontSize', fontsize);

propedit;

%% Figure 4B - pKa vs. sequence, 19F NMR

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig 4B",NumberTitle="off");
set(gcf,'Color','w');

errorbar(x_axis_seqs,sorted_pKas_19F,sorted_pKas_19F_errs,...
    'o','MarkerEdgeColor',color_G5FdU_Anion_plots,'MarkerFaceColor',color_G5FdU_Anion_plots,...
    'MarkerSize',markersize_pKa,'Color',color_G5FdU_Anion_plots,...
    'CapSize',err_capsize,'Linewidth',axeswidth);

xticklabels(sorted_seqs_T_19F);

axis([0 (length(sorted_seqs_T_19F)+1) y_axis_min_pKa y_axis_max_pKa]);

set(gca, 'FontSize', fontsize, 'XminorTick', 'off', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize],'XTick',x_axis_seqs);

ylabel('pKa', 'FontSize', font_title);

propedit;


%% Figure 4C - heatmap for dpKa T(N3) G5FdU Anion

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig 4C",NumberTitle="off");
set(gcf,'Color','w');

h3 = heatmap(seq_5p_dU,flip(seq_3p_dU),flip(dpKa_T_N3_G5FdU_Anion_os_sorted_heatmap),'CellLabelColor','none');
h3.MissingDataColor = [0.8 0.8 0.8];
h3.Colormap = earthyjet;
h3.ColorLimits = [-1.08 0.9];

set(gca, 'FontSize', fontsize);

propedit;


%% Figure 4D - box plots for dpKa T(N3) G5FdU Anion

vector_for_boxplot = [dpKa_T_N3_G5FdU_Anion_os_sorted(1:8),dpKa_T_N3_G5FdU_Anion_os_sorted(9:16)];
neighbor_3p_to_T_for_plot = [{'pyrimidine'},{'purine'}];

% figure('units','normalized','outerposition',[0 0 1 1],Name="Fig 4D",NumberTitle="off");
% set(gcf,'Color','w');
% 
% h = boxplot(vector_for_boxplot);
% set(h,'LineWidth',2);
% err1_upper_line = h(1);
% err1_lower_line = h(2);
% err1_top        = h(3);
% err1_bottom     = h(4);
% box1            = h(5);
% box1_median     = h(6);
% 
% err2_upper_line = h(8);
% err2_lower_line = h(9);
% err2_top        = h(10);
% err2_bottom     = h(11);
% box2            = h(12);
% box2_median     = h(13);
% 
% set(err1_upper_line, 'LineStyle','-','Color', color_box_plot,'LineWidth',linewidth);
% set(err1_lower_line, 'LineStyle','-','Color', color_box_plot,'LineWidth',linewidth);
% set(err1_top, 'LineStyle','-','Color', color_box_plot,'LineWidth',linewidth);
% set(err1_bottom, 'LineStyle','-','Color', color_box_plot,'LineWidth',linewidth);
% set(box1, 'LineStyle','-','Color', color_box_plot,'LineWidth',linewidth);
% set(box1_median, 'LineStyle','-','Color', color_box_mean,'LineWidth',linewidth);
% 
% set(err2_upper_line, 'LineStyle','-','Color', color_box_plot,'LineWidth',linewidth);
% set(err2_lower_line, 'LineStyle','-','Color', color_box_plot,'LineWidth',linewidth);
% set(err2_top, 'LineStyle','-','Color', color_box_plot,'LineWidth',linewidth);
% set(err2_bottom, 'LineStyle','-','Color', color_box_plot,'LineWidth',linewidth);
% set(box2, 'LineStyle','-','Color', color_box_plot,'LineWidth',linewidth);
% set(box2_median, 'LineStyle','-','Color', color_box_mean,'LineWidth',linewidth);
% 
% set(gca, 'FontSize', fontsize,'XminorTick', 'off', 'YminorTick', 'on',...
%     'LineWidth',axeswidth,'TickLength',[ticksize,ticksize]);
% 
% axis([0.5 (length(neighbor_3p_to_T_for_plot)+0.5) -1.3 1.3]);
% 
% xticklabels(neighbor_3p_to_T_for_plot);
% 
% xlabel('3p neighbor', 'FontSize', font_title);
% ylabel('dpKa', 'FontSize', font_title);
% 
% propedit;


%% Figure 2D - ddw of modified G5FdU vs. unmodified GT

%% correlation plot 5FdU(H3) in G5FdU WB vs T(H3) in GT WB

x_data      = dw_T_H3_GT_WB_os_sort_by_pKa;
x_data_errs = dw_T_H3_GT_WB_os_sort_by_pKa_errs;
y_data      = dw_T_H3_G5FdU_WB_os_sort_by_pKa;
y_data_errs = dw_T_H3_G5FdU_WB_os_sort_by_pKa_errs;

[H3_G5FdU_WB_vs_H3_GT_WB_r,H3_G5FdU_WB_vs_H3_GT_WB_Rsq,H3_G5FdU_WB_vs_H3_GT_WB_Rsq_det,...
    H3_G5FdU_WB_vs_H3_GT_WB_rms,H3_G5FdU_WB_vs_H3_GT_WB_rms_per,...
    bestfit_H3,x2_H3,inBetween_H3,...
    H3_G5FdU_WB_vs_H3_GT_WB_slope,H3_G5FdU_WB_vs_H3_GT_WB_inter,...
    H3_G5FdU_WB_vs_H3_GT_WB_slope_err,H3_G5FdU_WB_vs_H3_GT_WB_inter_err] = corr_plot(x_data,y_data,x_line,conf_int_per);

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig 2D top",NumberTitle="off");
set(gcf,'Color','w');

fill(x2_H3, inBetween_H3, 'k', 'FaceColor', color_shade_CI, 'FaceAlpha',.5,'LineStyle','none');

hold on;

errorbar(x_data,y_data,y_data_errs,y_data_errs,...
    x_data_errs,x_data_errs,...
    'o','MarkerEdgeColor',color_G5FdU_plots,'MarkerFaceColor',color_G5FdU_plots,'MarkerSize',22,...
    'Color',color_G5FdU_plots,'CapSize',err_capsize,'Linewidth',axeswidth);

plot(x_line,y_line,'k-','LineWidth',linewidth);

plot(x_line,bestfit_H3,'k--','LineWidth',err_linewidth);

hold off;

text(-0.45,0.48,sprintf('r = %0.2f',H3_G5FdU_WB_vs_H3_GT_WB_r),'FontSize',fontsize,'FontWeight','bold');
text(-0.45,0.30,sprintf('RMSD = %0.2f ppm',H3_G5FdU_WB_vs_H3_GT_WB_rms),'FontSize',fontsize,'FontWeight','bold');

axis([-0.5 0.5 -0.7 0.6]);

set(gca, 'FontSize', fontsize, 'XminorTick', 'on', 'YminorTick', 'on',...
    'XTick',-1.0:0.5:1.0,'YTick',-1.0:0.5:1.0,'LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

xlabel('dw T(H3) GT WB [ppm]', 'FontSize', font_title);
ylabel('dw U(H3) G5FdU WB [ppm]', 'FontSize', font_title);

propedit;


%% correlation plot G(H1) in G5FdU WB vs G(H1) in GT WB

x_data      = dw_G_H1_GT_WB_os_sort_by_pKa;
x_data_errs = dw_G_H1_GT_WB_os_sort_by_pKa_errs;
y_data      = dw_G_H1_G5FdU_WB_os_sort_by_pKa;
y_data_errs = dw_G_H1_G5FdU_WB_os_sort_by_pKa_errs;

[H1_G5FdU_WB_vs_H1_GT_WB_r,H1_G5FdU_WB_vs_H1_GT_WB_Rsq,H1_G5FdU_WB_vs_H1_GT_WB_Rsq_det,...
    H1_G5FdU_WB_vs_H1_GT_WB_rms,H1_G5FdU_WB_vs_H1_GT_WB_rms_per,...
    bestfit_H1,x2_H1,inBetween_H1,...
    H1_G5FdU_WB_vs_H1_GT_WB_slope,H1_G5FdU_WB_vs_H1_GT_WB_inter,...
    H1_G5FdU_WB_vs_H1_GT_WB_slope_err,H1_G5FdU_WB_vs_H1_GT_WB_inter_err] = corr_plot(x_data,y_data,x_line,conf_int_per);

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig 2D bottom",NumberTitle="off");
set(gcf,'Color','w');

fill(x2_H1, inBetween_H1, 'k', 'FaceColor', color_shade_CI, 'FaceAlpha',.5,'LineStyle','none');

hold on;

errorbar(x_data,y_data,y_data_errs,y_data_errs,...
    x_data_errs,x_data_errs,...
    'o','MarkerEdgeColor',color_G5FdU_plots,'MarkerFaceColor',color_G5FdU_plots,'MarkerSize',22,...
    'Color',color_G5FdU_plots,'CapSize',err_capsize,'Linewidth',axeswidth);

plot(x_line,y_line,'k-','LineWidth',linewidth);

plot(x_line,bestfit_H1,'k--','LineWidth',err_linewidth);

hold off;

text(-0.45,0.48,sprintf('r = %0.2f',H1_G5FdU_WB_vs_H1_GT_WB_r),'FontSize',fontsize,'FontWeight','bold');
text(-0.45,0.30,sprintf('RMSD = %0.2f ppm',H1_G5FdU_WB_vs_H1_GT_WB_rms),'FontSize',fontsize,'FontWeight','bold');

axis([-0.55 0.5 -0.8 0.65]);

set(gca, 'FontSize', fontsize, 'XminorTick', 'on', 'YminorTick', 'on',...
    'XTick',-1.0:0.5:1.0,'YTick',-1.0:0.5:1.0,'LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

xlabel('dw G(H1) GT WB [ppm]', 'FontSize', font_title);
ylabel('dw G(H1) G5FdU WB [ppm]', 'FontSize', font_title);

propedit;


diff_TH3 = abs(dw_T_H3_GT_WB_os_sort_by_pKa - dw_T_H3_G5FdU_WB_os_sort_by_pKa);
a = find(diff_TH3 > H3_G5FdU_WB_vs_H3_GT_WB_rms);
diff_TH3_a = diff_TH3(a);
seq_TH3_a = sorted_seqs_T_19F(a);

diff_GH1 = abs(dw_G_H1_GT_WB_os_sort_by_pKa - dw_G_H1_G5FdU_WB_os_sort_by_pKa);
a = find(diff_GH1 > H1_G5FdU_WB_vs_H1_GT_WB_rms);
diff_GH1_a = diff_GH1(a);
seq_GH1_a = sorted_seqs_T_19F(a);

%% Addivity calculations

%% loading list of unique sequences

filename_addit_uniq = 'parent-dest-seqs-unique.xlsx';
filepath_addit_uniq = sprintf('%s/%s',parent_path_pops_cs,filename_addit_uniq);

uniq_seqs_am = readtable(filepath_addit_uniq);

%% Figure S3B, left - Additivity T(H3) in G-T WB

[T_H3_GT_WB_additivity_tbl_ranked,T_H3_GT_WB_additivity_tbl_uniq,...
    parent_seq_cell] = additivity_calc(sorted_seqs_T_19F, cs_T_H3_GT_sort_by_pKa, uniq_seqs_am);

x_data = T_H3_GT_WB_additivity_tbl_uniq.pred_dw;
y_data = T_H3_GT_WB_additivity_tbl_uniq.obs_dw;

[T_H3_GT_WB_additivity_r,T_H3_GT_WB_additivity_Rsq,T_H3_GT_WB_additivity_Rsq_det,...
    T_H3_GT_WB_additivity_rms,T_H3_GT_WB_additivity_rms_per,...
    T_H3_GT_WB_additivity_bestfit,T_H3_GT_WB_additivity_x2,T_H3_GT_WB_additivity_inBetween,...
    T_H3_GT_WB_additivity_slope,T_H3_GT_WB_additivity_inter,...
    T_H3_GT_WB_additivity_slope_err,T_H3_GT_WB_additivity_inter_err] = corr_plot(x_data,y_data,x_line,conf_int_per);

[T_H3_GT_WB_enrch_outlier,T_H3_GT_WB_num_successes,T_H3_GT_WB_p_values,...
    T_H3_GT_WB_x_thresh_per,T_H3_GT_WB_y_thresh_per] = additivity_outliers(parent_seq_cell,...
    T_H3_GT_WB_additivity_tbl_ranked,T_H3_GT_WB_additivity_rms); % Outliers additivity T(H3) in G-T WB

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S3B left",NumberTitle="off");
set(gcf,'Color','w');

fill(T_H3_GT_WB_additivity_x2, T_H3_GT_WB_additivity_inBetween, 'k', 'FaceColor', color_shade_CI, 'FaceAlpha',.5,'LineStyle','none');

hold on;

plot(x_data,y_data,...
    'o','MarkerEdgeColor',color_GT,'MarkerFaceColor',color_GT,'MarkerSize',markersize_cs);

plot(x_line,y_line,'k-','LineWidth',err_linewidth);

plot(x_line,T_H3_GT_WB_additivity_bestfit,'k--','LineWidth',err_linewidth);

hold off;

axis([-1.15 1.15 -1.15 1.15]);

text(-1.0,0.98,sprintf('r = %0.2f',T_H3_GT_WB_additivity_r),'FontSize',fontsize,'FontWeight','bold');
text(-1.0,0.80,sprintf('RMSD = %0.2f ppm',T_H3_GT_WB_additivity_rms),'FontSize',fontsize,'FontWeight','bold');

set(gca, 'FontSize', fontsize, 'FontName', 'Arial',...
    'XminorTick', 'on', 'YminorTick', 'on',...
    'LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

xlabel('dw pred T(H3) GT [ppm]', 'FontSize', font_title);
ylabel('dw exp T(H3) GT [ppm]', 'FontSize', font_title);

propedit;

%% Figure S3B, right - Additivity G(H1) in G-T WB

[G_H1_GT_WB_additivity_tbl_ranked,G_H1_GT_WB_additivity_tbl_uniq,...
    parent_seq_cell] = additivity_calc(sorted_seqs_T_19F, cs_G_H1_GT_sort_by_pKa, uniq_seqs_am);

x_data = G_H1_GT_WB_additivity_tbl_uniq.pred_dw;
y_data = G_H1_GT_WB_additivity_tbl_uniq.obs_dw;

[G_H1_GT_WB_additivity_r,G_H1_GT_WB_additivity_Rsq,G_H1_GT_WB_additivity_Rsq_det,...
    G_H1_GT_WB_additivity_rms,G_H1_GT_WB_additivity_rms_per,...
    G_H1_GT_WB_additivity_bestfit,G_H1_GT_WB_additivity_x2,G_H1_GT_WB_additivity_inBetween,...
    G_H1_GT_WB_additivity_slope,G_H1_GT_WB_additivity_inter,...
    G_H1_GT_WB_additivity_slope_err,G_H1_GT_WB_additivity_inter_err] = corr_plot(x_data,y_data,x_line,conf_int_per);

[G_H1_GT_WB_enrch_outlier,G_H1_GT_WB_num_successes,G_H1_GT_WB_p_values,...
    G_H1_GT_WB_x_thresh_per,G_H1_GT_WB_y_thresh_per] = additivity_outliers(parent_seq_cell,...
    G_H1_GT_WB_additivity_tbl_ranked,G_H1_GT_WB_additivity_rms); % Outliers additivity G(H1) in G-T WB

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S3B right",NumberTitle="off");
set(gcf,'Color','w');

fill(G_H1_GT_WB_additivity_x2, G_H1_GT_WB_additivity_inBetween, 'k', 'FaceColor', color_shade_CI, 'FaceAlpha',.5,'LineStyle','none');

hold on;

plot(x_data,y_data,...
    'o','MarkerEdgeColor',color_GT,'MarkerFaceColor',color_GT,'MarkerSize',markersize_cs);

plot(x_line,y_line,'k-','LineWidth',err_linewidth);

plot(x_line,G_H1_GT_WB_additivity_bestfit,'k--','LineWidth',err_linewidth);

hold off;

axis([-1.15 1.15 -1.15 1.15]);

text(-1.0,0.98,sprintf('r = %0.2f',G_H1_GT_WB_additivity_r),'FontSize',fontsize,'FontWeight','bold');
text(-1.0,0.80,sprintf('RMSD = %0.2f ppm',G_H1_GT_WB_additivity_rms),'FontSize',fontsize,'FontWeight','bold');

set(gca, 'FontSize', fontsize, 'FontName', 'Arial',...
    'XminorTick', 'on', 'YminorTick', 'on',...
    'LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

xlabel('dw pred G(H1) GT [ppm]', 'FontSize', font_title);
ylabel('dw exp G(H1) GT [ppm]', 'FontSize', font_title);

propedit;

%% Figure S3C, left - Additivity 5FdU(H3) in G-5FdU WB

[U_H3_G5FdU_WB_additivity_tbl_ranked,U_H3_G5FdU_WB_additivity_tbl_uniq,...
    parent_seq_cell] = additivity_calc(sorted_seqs_T_19F, cs_U_H3_G5FdU_sort_by_pKa, uniq_seqs_am);

x_data = U_H3_G5FdU_WB_additivity_tbl_uniq.pred_dw;
y_data = U_H3_G5FdU_WB_additivity_tbl_uniq.obs_dw;

[U_H3_G5FdU_WB_additivity_r,U_H3_G5FdU_WB_additivity_Rsq,U_H3_G5FdU_WB_additivity_Rsq_det,...
    U_H3_G5FdU_WB_additivity_rms,U_H3_G5FdU_WB_additivity_rms_per,...
    U_H3_G5FdU_WB_additivity_bestfit,U_H3_G5FdU_WB_additivity_x2,U_H3_G5FdU_WB_additivity_inBetween,...
    U_H3_G5FdU_WB_additivity_slope,U_H3_G5FdU_WB_additivity_inter,...
    U_H3_G5FdU_WB_additivity_slope_err,U_H3_G5FdU_WB_additivity_inter_err] = corr_plot(x_data,y_data,x_line,conf_int_per);

[U_H3_G5FdU_WB_enrch_outlier,U_H3_G5FdU_WB_num_successes,U_H3_G5FdU_WB_p_values,...
    U_H3_G5FdU_WB_x_thresh_per,U_H3_G5FdU_WB_y_thresh_per] = additivity_outliers(parent_seq_cell,...
    U_H3_G5FdU_WB_additivity_tbl_ranked,U_H3_G5FdU_WB_additivity_rms); % Outliers additivity 5FdU(H3) in G-5FdU WB

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S3C left",NumberTitle="off");
set(gcf,'Color','w');

fill(U_H3_G5FdU_WB_additivity_x2, U_H3_G5FdU_WB_additivity_inBetween, 'k', 'FaceColor', color_shade_CI, 'FaceAlpha',.5,'LineStyle','none');

hold on;

plot(x_data,y_data,...
    'o','MarkerEdgeColor',color_G5FdU_plots,'MarkerFaceColor',color_G5FdU_plots,'MarkerSize',markersize_cs);

plot(x_line,y_line,'k-','LineWidth',err_linewidth);

plot(x_line,U_H3_G5FdU_WB_additivity_bestfit,'k--','LineWidth',err_linewidth);

hold off;

axis([-1.15 1.15 -1.15 1.15]);

text(-1.0,0.98,sprintf('r = %0.2f',U_H3_G5FdU_WB_additivity_r),'FontSize',fontsize,'FontWeight','bold');
text(-1.0,0.80,sprintf('RMSD = %0.2f ppm',U_H3_G5FdU_WB_additivity_rms),'FontSize',fontsize,'FontWeight','bold');

set(gca, 'FontSize', fontsize, 'FontName', 'Arial',...
    'XminorTick', 'on', 'YminorTick', 'on',...
    'LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

xlabel('dw pred U(H3) G5FdU [ppm]', 'FontSize', font_title);
ylabel('dw exp U(H3) G5FdU [ppm]', 'FontSize', font_title);

propedit;

%% Figure S3C, right - Additivity G(H1) in G-5FdU WB

[G_H1_G5FdU_WB_additivity_tbl_ranked,G_H1_G5FdU_WB_additivity_tbl_uniq,...
    parent_seq_cell] = additivity_calc(sorted_seqs_T_19F, cs_G_H1_G5FdU_sort_by_pKa, uniq_seqs_am);

x_data = G_H1_G5FdU_WB_additivity_tbl_uniq.pred_dw;
y_data = G_H1_G5FdU_WB_additivity_tbl_uniq.obs_dw;

[G_H1_G5FdU_WB_additivity_r,G_H1_G5FdU_WB_additivity_Rsq,G_H1_G5FdU_WB_additivity_Rsq_det,...
    G_H1_G5FdU_WB_additivity_rms,G_H1_G5FdU_WB_additivity_rms_per,...
    G_H1_G5FdU_WB_additivity_bestfit,G_H1_G5FdU_WB_additivity_x2,G_H1_G5FdU_WB_additivity_inBetween,...
    G_H1_G5FdU_WB_additivity_slope,G_H1_G5FdU_WB_additivity_inter,...
    G_H1_G5FdU_WB_additivity_slope_err,G_H1_G5FdU_WB_additivity_inter_err] = corr_plot(x_data,y_data,x_line,conf_int_per);

[G_H1_G5FdU_WB_enrch_outlier,G_H1_G5FdU_WB_num_successes,G_H1_G5FdU_WB_p_values,...
    G_H1_G5FdU_WB_x_thresh_per,G_H1_G5FdU_WB_y_thresh_per] = additivity_outliers(parent_seq_cell,...
    G_H1_G5FdU_WB_additivity_tbl_ranked,G_H1_G5FdU_WB_additivity_rms); % Outliers additivity G(H1) in G-5FdU WB

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S3C right",NumberTitle="off");
set(gcf,'Color','w');

fill(G_H1_G5FdU_WB_additivity_x2, G_H1_G5FdU_WB_additivity_inBetween, 'k', 'FaceColor', color_shade_CI, 'FaceAlpha',.5,'LineStyle','none');

hold on;

plot(x_data,y_data,...
    'o','MarkerEdgeColor',color_G5FdU_plots,'MarkerFaceColor',color_G5FdU_plots,'MarkerSize',markersize_cs);

plot(x_line,y_line,'k-','LineWidth',err_linewidth);

plot(x_line,G_H1_G5FdU_WB_additivity_bestfit,'k--','LineWidth',err_linewidth);

hold off;

axis([-1.15 1.15 -1.15 1.15]);

text(-1.0,0.98,sprintf('r = %0.2f',G_H1_G5FdU_WB_additivity_r),'FontSize',fontsize,'FontWeight','bold');
text(-1.0,0.80,sprintf('RMSD = %0.2f ppm',G_H1_G5FdU_WB_additivity_rms),'FontSize',fontsize,'FontWeight','bold');

set(gca, 'FontSize', fontsize, 'FontName', 'Arial',...
    'XminorTick', 'on', 'YminorTick', 'on',...
    'LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

xlabel('dw pred G(H1) G5FdU [ppm]', 'FontSize', font_title);
ylabel('dw exp G(H1) G5FdU [ppm]', 'FontSize', font_title);

propedit;

%% Figure S4A, right - Additivity 5FdU(F5) in G-5FdU WB

[U_F5_G5FdU_WB_additivity_tbl_ranked,U_F5_G5FdU_WB_additivity_tbl_uniq,...
    parent_seq_cell] = additivity_calc(sorted_seqs_T_19F, cs_F5_G5FdU_WB_sort_by_pKa, uniq_seqs_am);

x_data = U_F5_G5FdU_WB_additivity_tbl_uniq.pred_dw;
y_data = U_F5_G5FdU_WB_additivity_tbl_uniq.obs_dw;

[U_F5_G5FdU_WB_additivity_r,U_F5_G5FdU_WB_additivity_Rsq,U_F5_G5FdU_WB_additivity_Rsq_det,...
    U_F5_G5FdU_WB_additivity_rms,U_F5_G5FdU_WB_additivity_rms_per,...
    U_F5_G5FdU_WB_additivity_bestfit,U_F5_G5FdU_WB_additivity_x2,U_F5_G5FdU_WB_additivity_inBetween,...
    U_F5_G5FdU_WB_additivity_slope,U_F5_G5FdU_WB_additivity_inter,...
    U_F5_G5FdU_WB_additivity_slope_err,U_F5_G5FdU_WB_additivity_inter_err] = corr_plot(x_data,y_data,x_line,conf_int_per);

[U_F5_G5FdU_WB_enrch_outlier,U_F5_G5FdU_WB_num_successes,U_F5_G5FdU_WB_p_values,...
    U_F5_G5FdU_WB_x_thresh_per,U_F5_G5FdU_WB_y_thresh_per] = additivity_outliers(parent_seq_cell,...
    U_F5_G5FdU_WB_additivity_tbl_ranked,U_F5_G5FdU_WB_additivity_rms); % Outliers additivity 5FdU(F5) in G-5FdU WB

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S4A right",NumberTitle="off");
set(gcf,'Color','w');

fill(U_F5_G5FdU_WB_additivity_x2, U_F5_G5FdU_WB_additivity_inBetween, 'k', 'FaceColor', color_shade_CI, 'FaceAlpha',.5,'LineStyle','none');

hold on;

plot(x_data,y_data,...
    'o','MarkerEdgeColor',color_G5FdU_plots,'MarkerFaceColor',color_G5FdU_plots,'MarkerSize',markersize_cs);

plot(x_line,y_line,'k-','LineWidth',err_linewidth);

plot(x_line,U_F5_G5FdU_WB_additivity_bestfit,'k--','LineWidth',err_linewidth);

hold off;

axis([-2.0 2.0 -2.0 2.0]);

text(-1.8,1.5,sprintf('r = %0.2f',U_F5_G5FdU_WB_additivity_r),'FontSize',fontsize,'FontWeight','bold');
text(-1.8,1.2,sprintf('RMSD = %0.2f ppm',U_F5_G5FdU_WB_additivity_rms),'FontSize',fontsize,'FontWeight','bold');

set(gca, 'FontSize', fontsize, 'FontName', 'Arial',...
    'XminorTick', 'on', 'YminorTick', 'on',...
    'LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

xlabel('dw pred U(F5) G5FdU WB [ppm]', 'FontSize', font_title);
ylabel('dw exp U(F5) G5FdU WB [ppm]', 'FontSize', font_title);

propedit;


%% Figure S4B, right - Additivity 5FdU(F5) in G-5FdU Anion

[U_F5_G5FdU_Anion_additivity_tbl_ranked,U_F5_G5FdU_Anion_additivity_tbl_uniq,...
    parent_seq_cell] = additivity_calc(sorted_seqs_T_19F, cs_F5_G5FdU_Anion_sort_by_pKa, uniq_seqs_am);

x_data = U_F5_G5FdU_Anion_additivity_tbl_uniq.pred_dw;
y_data = U_F5_G5FdU_Anion_additivity_tbl_uniq.obs_dw;

[U_F5_G5FdU_Anion_additivity_r,U_F5_G5FdU_Anion_additivity_Rsq,U_F5_G5FdU_Anion_additivity_Rsq_det,...
    U_F5_G5FdU_Anion_additivity_rms,U_F5_G5FdU_Anion_additivity_rms_per,...
    U_F5_G5FdU_Anion_additivity_bestfit,U_F5_G5FdU_Anion_additivity_x2,U_F5_G5FdU_Anion_additivity_inBetween,...
    U_F5_G5FdU_Anion_additivity_slope,U_F5_G5FdU_Anion_additivity_inter,...
    U_F5_G5FdU_Anion_additivity_slope_err,U_F5_G5FdU_Anion_additivity_inter_err] = corr_plot(x_data,y_data,x_line,conf_int_per);

[U_F5_G5FdU_Anion_enrch_outlier,U_F5_G5FdU_Anion_num_successes,U_F5_G5FdU_Anion_p_values,...
    U_F5_G5FdU_Anion_x_thresh_per,U_F5_G5FdU_Anion_y_thresh_per] = additivity_outliers(parent_seq_cell,...
    U_F5_G5FdU_Anion_additivity_tbl_ranked,U_F5_G5FdU_Anion_additivity_rms); % Outliers additivity 5FdU(F5) in G-5FdU Anion

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S4B right",NumberTitle="off");
set(gcf,'Color','w');

fill(U_F5_G5FdU_Anion_additivity_x2, U_F5_G5FdU_Anion_additivity_inBetween, 'k', 'FaceColor', color_shade_CI, 'FaceAlpha',.5,'LineStyle','none');

hold on;

plot(x_data,y_data,...
    'o','MarkerEdgeColor',color_G5FdU_Anion_plots,'MarkerFaceColor',color_G5FdU_Anion_plots,'MarkerSize',markersize_cs);

plot(x_line,y_line,'k-','LineWidth',err_linewidth);

plot(x_line,U_F5_G5FdU_Anion_additivity_bestfit,'k--','LineWidth',err_linewidth);

hold off;

axis([-2.0 2.0 -2.0 2.0]);

text(-1.8,1.5,sprintf('r = %0.2f',U_F5_G5FdU_Anion_additivity_r),'FontSize',fontsize,'FontWeight','bold');
text(-1.8,1.2,sprintf('RMSD = %0.2f ppm',U_F5_G5FdU_Anion_additivity_rms),'FontSize',fontsize,'FontWeight','bold');

set(gca, 'FontSize', fontsize, 'FontName', 'Arial',...
    'XminorTick', 'on', 'YminorTick', 'on',...
    'LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

xlabel('dw pred U(F5) G5FdU Anion [ppm]', 'FontSize', font_title);
ylabel('dw exp U(F5) G5FdU Anion [ppm]', 'FontSize', font_title);

propedit;


%% Figure 4E - Additivity pKa U(N3) in G-5FdU Anion

[pKa_U_N3_G5FdU_Anion_additivity_tbl_ranked,pKa_U_N3_G5FdU_Anion_additivity_tbl_uniq,...
    parent_seq_cell] = additivity_calc(sorted_seqs_T_19F, sorted_pKas_19F, uniq_seqs_am);

x_data = pKa_U_N3_G5FdU_Anion_additivity_tbl_uniq.pred_dw;
y_data = pKa_U_N3_G5FdU_Anion_additivity_tbl_uniq.obs_dw;

[pKa_U_N3_G5FdU_Anion_additivity_r,pKa_U_N3_G5FdU_Anion_additivity_Rsq,pKa_U_N3_G5FdU_Anion_additivity_Rsq_det,...
    pKa_U_N3_G5FdU_Anion_additivity_rms,pKa_U_N3_G5FdU_Anion_additivity_rms_per,...
    pKa_U_N3_G5FdU_Anion_additivity_bestfit,pKa_U_N3_G5FdU_Anion_additivity_x2,pKa_U_N3_G5FdU_Anion_additivity_inBetween,...
    pKa_U_N3_G5FdU_Anion_additivity_slope,pKa_U_N3_G5FdU_Anion_additivity_inter,...
    pKa_U_N3_G5FdU_Anion_additivity_slope_err,pKa_U_N3_G5FdU_Anion_additivity_inter_err] = corr_plot(x_data,y_data,x_line,conf_int_per);

[pKa_U_N3_G5FdU_Anion_enrch_outlier,pKa_U_N3_G5FdU_Anion_num_successes,pKa_U_N3_G5FdU_Anion_p_values,...
    pKa_U_N3_G5FdU_Anion_x_thresh_per,pKa_U_N3_G5FdU_Anion_y_thresh_per] = additivity_outliers(parent_seq_cell,...
    pKa_U_N3_G5FdU_Anion_additivity_tbl_ranked,pKa_U_N3_G5FdU_Anion_additivity_rms); % Outliers additivity 5FdU(F5) in G-5FdU Anion

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig 4E",NumberTitle="off");
set(gcf,'Color','w');

fill(pKa_U_N3_G5FdU_Anion_additivity_x2, pKa_U_N3_G5FdU_Anion_additivity_inBetween, 'k', 'FaceColor', color_shade_CI, 'FaceAlpha',.5,'LineStyle','none');

hold on;

plot(x_data,y_data,...
    'o','MarkerEdgeColor',color_G5FdU_Anion_plots,'MarkerFaceColor',color_G5FdU_Anion_plots,'MarkerSize',markersize_cs);

plot(x_line,y_line,'k-','LineWidth',err_linewidth);

plot(x_line,pKa_U_N3_G5FdU_Anion_additivity_bestfit,'k--','LineWidth',err_linewidth);

hold off;

axis([-2.0 2.0 -2.0 2.0]);

text(-1.8,1.5,sprintf('r = %0.2f',pKa_U_N3_G5FdU_Anion_additivity_r),'FontSize',fontsize,'FontWeight','bold');
text(-1.8,1.2,sprintf('RMSD = %0.2f',pKa_U_N3_G5FdU_Anion_additivity_rms),'FontSize',fontsize,'FontWeight','bold');

set(gca, 'FontSize', fontsize, 'FontName', 'Arial',...
    'XminorTick', 'on', 'YminorTick', 'on',...
    'LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

xlabel('dpKa pred [ppm]', 'FontSize', font_title);
ylabel('dpKa exp [ppm]', 'FontSize', font_title);

propedit;


%% Plot 1D 19F spectra

%% load 1D 19F spectra

path_1D_19F = 'NMR-Spectra/1D-19F-spectra';
parent_path_1D_19F = sprintf('%s/%s',parent_path,path_1D_19F);

sequence_tri_1D_19F = {'GTC';'GTT';'GTA';'GTG';'ATC';'ATT';'ATA';'ATG';...
    'TTC';'TTT';'TTA';'TTG';'CTC';'CTT';'CTA';'CTG'};

pHs_1D_19F = {'Low';'Med1';'Med2';'Med3';'High'};

pHs_1D_19F_Low  = [6.88;6.84;6.87;7.01;6.88;6.91;6.92;6.82;...
    NaN;6.84;6.86;6.85;6.80;6.86;6.87;6.85];
pHs_1D_19F_Med1 = [8.10;8.07;8.26;8.20;8.11;8.08;8.38;8.14;...
    8.12;8.05;8.07;8.08;7.97;8.02;8.37;8.13];
pHs_1D_19F_Med2 = [NaN;NaN;8.96;8.91;9.22;NaN;8.95;9.055;...
    NaN;8.75;8.92;9.10;NaN;9.14;8.92;8.99];
pHs_1D_19F_Med3 = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;9.36;...
    NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];
pHs_1D_19F_High = [10.4;10.9;10.7;10.2;10.3;10.1;11.0;10.8;...
    NaN;10.4;10.7;10.9;NaN;10.8;10.9;10.6];

pHs_1D_19F_all_tbl = table(sequence_tri_1D_19F,pHs_1D_19F_Low,...
    pHs_1D_19F_Med1,pHs_1D_19F_Med2,pHs_1D_19F_Med3,pHs_1D_19F_High,...
    'VariableNames',{'seq_T', 'Low', 'Med1', 'Med2', 'Med3', 'High'});

F2_left_1D_19F  = -158.0;
F2_right_1D_19F = -170.0;

specs_1D_19F = cell(length(sequence_tri_1D_19F),6);

for i = 1:length(sequence_tri_1D_19F)

    specs_1D_19F{i,1} = sequence_tri_1D_19F{i};

    for j = 1:length(pHs_1D_19F)

        filename = sprintf('%s_G5FdU_1D_19F_pH_%s.txt',sequence_tri_1D_19F{i},pHs_1D_19F{j});

        [axis_19F,spec_19F] = load_1D_spec(parent_path_1D_19F,filename,F2_left_1D_19F,F2_right_1D_19F);

        if ~isempty(axis_19F)

            specs_1D_19F{i,j+1} = [axis_19F,spec_19F];

        end

    end

end

specs_1D_19F_tbl = cell2table(specs_1D_19F,"VariableNames",["seq_T" "Low" "Med1" "Med2" "Med3" "High"]);


%% Plotting 1D 19F spectra at different pHs - offset vertically

%% scaling for plotting Med pH spectra

scale_GTC = 2.0;
scale_GTT = 1.8;
scale_GTA = 1.2;
scale_GTG = 2.0;
scale_ATC = 1.6;
scale_ATT = 2.2;
scale_ATA = 2.0;
scale_ATG = 2.5;
scale_TTC = 4.8;
scale_TTT = 3.0;
scale_TTA = 1.7;
scale_TTG = 2.1;
scale_CTC = 1.0;
scale_CTT = 1.5;
scale_CTA = 2.5;
scale_CTG = 2.1;

scale_GTA_for_all_pHs = scale_GTA./scale_GTA;
scale_CTG_for_all_pHs = scale_CTG./scale_GTA;

scales_1D_19F_for_fits = [scale_GTC; scale_GTT; scale_GTA; scale_GTG;...
    scale_ATC; scale_ATT; scale_ATA; scale_ATG;...
    scale_TTC; scale_TTT; scale_TTA; scale_TTG;...
    scale_CTC; scale_CTT; scale_CTA; scale_CTG];

vert_offset_19F = 8E6;

%% Figure 3A, left - 1D 19F spectra GTA different pHs

[~,ind_GTA] = ismember('GTA',specs_1D_19F_tbl.seq_T);

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig 3A left",NumberTitle="off");
set(gcf,'Color','w');

plot(specs_1D_19F_tbl.Low{ind_GTA,1}(:,1),(specs_1D_19F_tbl.Low{ind_GTA,1}(:,2))'+vert_offset_19F*2.0,...
    'Color',color_G5FdU_lowpH,'LineWidth',linewidth);

hold on;

plot(specs_1D_19F_tbl.Med1{ind_GTA,1}(:,1),(specs_1D_19F_tbl.Med1{ind_GTA,1}(:,2)*scale_GTA_for_all_pHs)'+vert_offset_19F,...
    'Color',color_G5FdU_medpH,'LineWidth',linewidth);

plot(specs_1D_19F_tbl.High{ind_GTA,1}(:,1),(specs_1D_19F_tbl.High{ind_GTA,1}(:,2)*0.6)',...
    'Color',color_G5FdU_Anion,'LineWidth',linewidth);

hold off;

text(-160.3,2.6E7,'GTA',...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);

text(-161.8,2.3E7,'G-5FdU',...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);

text(-163.9,0.4E7,'G-5FdU^{-}',...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_Anion);

text(-160.3,1.7E7,sprintf('pH %0.1f',pHs_1D_19F_all_tbl.Low(ind_GTA)),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);

text(-160.3,(1.7E7-vert_offset_19F),sprintf('pH %0.1f',pHs_1D_19F_all_tbl.Med1(ind_GTA)),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_medpH);

text(-160.3,(1.7E7-vert_offset_19F*2.0),sprintf('pH %0.1f',pHs_1D_19F_all_tbl.High(ind_GTA)),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_Anion);

set(gca, 'FontSize', fontsize, 'FontName', 'Arial', 'Xdir', 'rev',...
    'XminorTick', 'on', 'YTick',[],'YminorTick', 'off','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([x_min_19F x_max_19F y_min_19F y_max_19F+18E6]);

propedit;


%% Figure 3A, right - 1D 19F spectra CTG different pHs

[~,ind_CTG] = ismember('CTG',specs_1D_19F_tbl.seq_T);

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig 3A right",NumberTitle="off");
set(gcf,'Color','w');

plot(specs_1D_19F_tbl.Low{ind_CTG,1}(:,1),(specs_1D_19F_tbl.Low{ind_CTG,1}(:,2))'+vert_offset_19F*2.0,...
    'Color',color_G5FdU_lowpH,'LineWidth',linewidth);

hold on;

plot(specs_1D_19F_tbl.Med1{ind_CTG,1}(:,1),(specs_1D_19F_tbl.Med2{ind_CTG,1}(:,2)*scale_CTG_for_all_pHs)'+vert_offset_19F,...
    'Color',color_G5FdU_medpH,'LineWidth',linewidth);

plot(specs_1D_19F_tbl.High{ind_CTG,1}(:,1),(specs_1D_19F_tbl.High{ind_CTG,1}(:,2)*0.8)',...
    'Color',color_G5FdU_Anion,'LineWidth',linewidth);

hold off;

text(-160.3,2.6E7,'CTG',...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);

text(-164.2,2.3E7,'G-5FdU',...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);

text(-161.6,0.4E7,'G-5FdU^{-}',...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_Anion);

text(-160.3,1.7E7,sprintf('pH %0.1f',pHs_1D_19F_all_tbl.Low(ind_CTG)),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);

text(-160.3,(1.7E7-vert_offset_19F),sprintf('pH %0.1f',pHs_1D_19F_all_tbl.Med1(ind_CTG)),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_medpH);

text(-160.3,(1.7E7-vert_offset_19F*2.0),sprintf('pH %0.1f',pHs_1D_19F_all_tbl.High(ind_CTG)),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_Anion);

set(gca, 'FontSize', fontsize, 'FontName', 'Arial', 'Xdir', 'rev',...
    'XminorTick', 'on', 'YTick',[],'YminorTick', 'off','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([x_min_19F x_max_19F y_min_19F y_max_19F+18E6]);

propedit;

%% Lorenzian fits

filename_init_fit_params  = 'fitting-params-19F-nmr.xlsx';
filepath                  = sprintf('%s/%s',parent_path_1D_19F,filename_init_fit_params);

init_fit_params_Med1 = readtable(filepath,'Sheet','Med1');
init_fit_params_Med2 = readtable(filepath,'Sheet','Med2');
init_fit_params_Med3 = readtable(filepath,'Sheet','Med3');

noise_lims   = [-160.0, -161.5];
x_min_fits = -166.0;
x_max_fits = -161.5;
x_fit = -160:-0.01:-170;

specs_1D_19F_Med1_fits_params      = cell(height(specs_1D_19F_tbl),11);
specs_1D_19F_Med1_fits_params(:,1) = specs_1D_19F_tbl.seq_T;

specs_1D_19F_Med2_fits_params      = cell(height(specs_1D_19F_tbl),11);
specs_1D_19F_Med2_fits_params(:,1) = specs_1D_19F_tbl.seq_T;

specs_1D_19F_Med3_fits_params      = cell(height(specs_1D_19F_tbl),11);
specs_1D_19F_Med3_fits_params(:,1) = specs_1D_19F_tbl.seq_T;

for i = 1:height(specs_1D_19F_tbl)

    if ~isnan(init_fit_params_Med1.pH(i))

        spec_Med1         = cell2mat(specs_1D_19F_tbl.Med1(i));
        pH_Med1           = pHs_1D_19F_Med1(i);
        x_data_Med1_all   = spec_Med1(:,1);
        noise_region_Med1 = spec_Med1(x_data_Med1_all>noise_lims(2)&x_data_Med1_all<noise_lims(1),2);
        noise_mean_Med1   = mean(noise_region_Med1);

        y_data_Med1  = spec_Med1(x_data_Med1_all>(x_min_fits)&x_data_Med1_all<x_max_fits,2) - noise_mean_Med1;
        x_data_Med1  = x_data_Med1_all(x_data_Med1_all>(x_min_fits)&x_data_Med1_all<x_max_fits);

        if ismember('fit1Lorentzian',init_fit_params_Med1.fitFunc(i))

            x0_WB_amp   = init_fit_params_Med1.x0_WB_amp(i);
            x0_WB_width = init_fit_params_Med1.x0_WB_width(i);
            x0_WB_cs    = init_fit_params_Med1.x0_WB_cs(i);

            x0 = [x0_WB_amp x0_WB_width x0_WB_cs];

            [FitLorentzians_Med1,resnorm_Med1,exitflag_Med1,output_Med1,...
                FitLorentz1_Med1,...
                cs_WB_Med1] = fit1Lorentzian(x_data_Med1, y_data_Med1, x0, x_fit, pH_Med1);

            FitLorentz2_Med1 = [];
            FitLorentz3_Med1 = [];
            pWB_Med1 = [];
            pAnion_Med1 = [];
            pES_Med1 = [];
            cs_Anion_Med1 = [];
            cs_ES_Med1 = [];
            pKa_Med1 = [];

        elseif ismember('fit2Lorentzians',init_fit_params_Med1.fitFunc(i))

            x0_anion_amp   = init_fit_params_Med1.x0_anion_amp(i);
            x0_anion_width = init_fit_params_Med1.x0_anion_width(i);
            x0_anion_cs    = init_fit_params_Med1.x0_anion_cs(i);
            x0_WB_amp      = init_fit_params_Med1.x0_WB_amp(i);
            x0_WB_width    = init_fit_params_Med1.x0_WB_width(i);
            x0_WB_cs       = init_fit_params_Med1.x0_WB_cs(i);

            x0 = [x0_anion_amp x0_anion_width x0_anion_cs x0_WB_amp x0_WB_width x0_WB_cs];

            [FitLorentzians_Med1,resnorm_Med1,exitflag_Med1,output_Med1,...
                FitLorentz1_Med1,FitLorentz2_Med1,...
                pWB_Med1,pAnion_Med1,pKa_Med1,...
                cs_WB_Med1,cs_Anion_Med1,dcs_Anion_Med1] = fit2Lorentzians(x_data_Med1, y_data_Med1, x0, x_fit, pH_Med1);

            FitLorentz3_Med1 = [];
            pES_Med1 = [];
            cs_ES_Med1 = [];

        elseif ismember('fit3Lorentzians',init_fit_params_Med1.fitFunc(i))

            x0_anion_amp   = init_fit_params_Med1.x0_anion_amp(i);
            x0_anion_width = init_fit_params_Med1.x0_anion_width(i);
            x0_anion_cs    = init_fit_params_Med1.x0_anion_cs(i);
            x0_WB_amp      = init_fit_params_Med1.x0_WB_amp(i);
            x0_WB_width    = init_fit_params_Med1.x0_WB_width(i);
            x0_WB_cs       = init_fit_params_Med1.x0_WB_cs(i);
            x0_ES_amp      = init_fit_params_Med1.x0_ES_amp(i);
            x0_ES_width    = init_fit_params_Med1.x0_ES_width(i);
            x0_ES_cs       = init_fit_params_Med1.x0_ES_cs(i);

            x0 = [x0_anion_amp x0_anion_width x0_anion_cs x0_WB_amp x0_WB_width x0_WB_cs x0_ES_amp x0_ES_width x0_ES_cs];

            [FitLorentzians_Med1,resnorm_Med1,exitflag_Med1,output_Med1,...
                FitLorentz1_Med1,FitLorentz2_Med1,FitLorentz3_Med1,...
                pWB_Med1,pAnion_Med1,pES_Med1,pKa_Med1,...
                cs_WB_Med1,cs_Anion_Med1,dcs_Anion_Med1,cs_ES_Med1] = fit3Lorentzians(x_data_Med1, y_data_Med1, x0, x_fit, pH_Med1);


        elseif ismember('fit2LorentziansConst',init_fit_params_Med1.fitFunc(i))

            x0_anion_amp   = init_fit_params_Med1.x0_anion_amp(i);
            x0_anion_width = init_fit_params_Med1.x0_anion_width(i);
            x0_anion_cs    = init_fit_params_Med1.x0_anion_cs(i);
            x0_WB_amp      = init_fit_params_Med1.x0_WB_amp(i);
            x0_WB_width    = init_fit_params_Med1.x0_WB_width(i);
            x0_WB_cs       = init_fit_params_Med1.x0_WB_cs(i);

            x0 = [x0_anion_amp x0_anion_width x0_anion_cs x0_WB_amp x0_WB_width x0_WB_cs];

            if i == 7
                lb = [-inf,-inf,x0_anion_cs,-inf,-inf,-inf];
                ub = [inf,inf,x0_anion_cs,inf,inf,inf];
            elseif i == 12 || i == 15
                lb = [-inf,x0_anion_width,x0_anion_cs,-inf,-inf,-inf];
                ub = [inf,x0_anion_width,x0_anion_cs,inf,inf,inf];
            end

            [FitLorentzians_Med1,resnorm_Med1,exitflag_Med1,output_Med1,...
                FitLorentz1_Med1,FitLorentz2_Med1,...
                pWB_Med1,pAnion_Med1,pKa_Med1,...
                cs_WB_Med1,cs_Anion_Med1,dcs_Anion_Med1] = fit2LorentziansConst(x_data_Med1, y_data_Med1, x0, lb, ub, x_fit, pH_Med1);

            FitLorentz3_Med1 = [];
            pES_Med1 = [];
            cs_ES_Med1 = [];

        elseif ismember('fit3LorentziansConst',init_fit_params_Med1.fitFunc(i))

            x0_anion_amp   = init_fit_params_Med1.x0_anion_amp(i);
            x0_anion_width = init_fit_params_Med1.x0_anion_width(i);
            x0_anion_cs    = init_fit_params_Med1.x0_anion_cs(i);
            x0_WB_amp      = init_fit_params_Med1.x0_WB_amp(i);
            x0_WB_width    = init_fit_params_Med1.x0_WB_width(i);
            x0_WB_cs       = init_fit_params_Med1.x0_WB_cs(i);
            x0_ES_amp      = init_fit_params_Med1.x0_ES_amp(i);
            x0_ES_width    = init_fit_params_Med1.x0_ES_width(i);
            x0_ES_cs       = init_fit_params_Med1.x0_ES_cs(i);

            x0 = [x0_anion_amp x0_anion_width x0_anion_cs x0_WB_amp x0_WB_width x0_WB_cs x0_ES_amp x0_ES_width x0_ES_cs];

            if i == 2
                lb = [-inf,-inf,-inf,-inf,-inf,-inf,x0_ES_amp,x0_ES_width,-inf];
                ub = [inf,inf,inf,inf,inf,inf,x0_ES_amp,x0_ES_width,inf];
            end

            [FitLorentzians_Med1,resnorm_Med1,exitflag_Med1,output_Med1,...
                FitLorentz1_Med1,FitLorentz2_Med1,FitLorentz3_Med1,...
                pWB_Med1,pAnion_Med1,pES_Med1,pKa_Med1,...
                cs_WB_Med1,cs_Anion_Med1,dcs_Anion_Med1,cs_ES_Med1] = fit3LorentziansConst(x_data_Med1, y_data_Med1, x0, lb, ub, x_fit, pH_Med1);

        end

        specs_1D_19F_Med1_fits_params{i,2} = FitLorentz1_Med1;
        specs_1D_19F_Med1_fits_params{i,3} = FitLorentz2_Med1;
        specs_1D_19F_Med1_fits_params{i,4} = FitLorentz3_Med1;

        specs_1D_19F_Med1_fits_params{i,5} = cs_WB_Med1;
        specs_1D_19F_Med1_fits_params{i,6} = cs_Anion_Med1;
        specs_1D_19F_Med1_fits_params{i,7} = cs_ES_Med1;
        specs_1D_19F_Med1_fits_params{i,8} = pWB_Med1;
        specs_1D_19F_Med1_fits_params{i,9} = pAnion_Med1;
        specs_1D_19F_Med1_fits_params{i,10} = pES_Med1;
        specs_1D_19F_Med1_fits_params{i,11} = pKa_Med1;

    end

end


for i = 1:height(specs_1D_19F_tbl)

    if ~isnan(init_fit_params_Med2.pH(i))

        spec_Med2         = cell2mat(specs_1D_19F_tbl.Med2(i));
        pH_Med2           = pHs_1D_19F_Med2(i);
        x_data_Med2_all   = spec_Med2(:,1);
        noise_region_Med2 = spec_Med2(x_data_Med2_all>noise_lims(2)&x_data_Med2_all<noise_lims(1),2);
        noise_mean_Med2   = mean(noise_region_Med2);

        y_data_Med2  = spec_Med2(x_data_Med2_all>(x_min_fits)&x_data_Med2_all<x_max_fits,2) - noise_mean_Med2;
        x_data_Med2  = x_data_Med2_all(x_data_Med2_all>(x_min_fits)&x_data_Med2_all<x_max_fits);

        if ismember('fit1Lorentzian',init_fit_params_Med2.fitFunc(i))

            x0_WB_amp   = init_fit_params_Med2.x0_WB_amp(i);
            x0_WB_width = init_fit_params_Med2.x0_WB_width(i);
            x0_WB_cs    = init_fit_params_Med2.x0_WB_cs(i);

            x0 = [x0_WB_amp x0_WB_width x0_WB_cs];

            [FitLorentzians_Med2,resnorm_Med2,exitflag_Med2,output_Med2,...
                FitLorentz1_Med2,...
                cs_WB_Med2] = fit1Lorentzian(x_data_Med2, y_data_Med2, x0, x_fit, pH_Med2);

            FitLorentz2_Med2 = [];
            FitLorentz3_Med2 = [];
            pWB_Med2 = [];
            pAnion_Med2 = [];
            pES_Med2 = [];
            cs_Anion_Med2 = [];
            cs_ES_Med2 = [];
            pKa_Med2 = [];

        elseif ismember('fit2Lorentzians',init_fit_params_Med2.fitFunc(i))

            x0_anion_amp   = init_fit_params_Med2.x0_anion_amp(i);
            x0_anion_width = init_fit_params_Med2.x0_anion_width(i);
            x0_anion_cs    = init_fit_params_Med2.x0_anion_cs(i);
            x0_WB_amp      = init_fit_params_Med2.x0_WB_amp(i);
            x0_WB_width    = init_fit_params_Med2.x0_WB_width(i);
            x0_WB_cs       = init_fit_params_Med2.x0_WB_cs(i);

            x0 = [x0_anion_amp x0_anion_width x0_anion_cs x0_WB_amp x0_WB_width x0_WB_cs];

            [FitLorentzians_Med2,resnorm_Med2,exitflag_Med2,output_Med2,...
                FitLorentz1_Med2,FitLorentz2_Med2,...
                pWB_Med2,pAnion_Med2,pKa_Med2,...
                cs_WB_Med2,cs_Anion_Med2,dcs_Anion_Med2] = fit2Lorentzians(x_data_Med2, y_data_Med2, x0, x_fit, pH_Med2);

            FitLorentz3_Med2 = [];
            pES_Med2 = [];
            cs_ES_Med2 = [];

        elseif ismember('fit3Lorentzians',init_fit_params_Med2.fitFunc(i))

            x0_anion_amp   = init_fit_params_Med2.x0_anion_amp(i);
            x0_anion_width = init_fit_params_Med2.x0_anion_width(i);
            x0_anion_cs    = init_fit_params_Med2.x0_anion_cs(i);
            x0_WB_amp      = init_fit_params_Med2.x0_WB_amp(i);
            x0_WB_width    = init_fit_params_Med2.x0_WB_width(i);
            x0_WB_cs       = init_fit_params_Med2.x0_WB_cs(i);
            x0_ES_amp      = init_fit_params_Med2.x0_ES_amp(i);
            x0_ES_width    = init_fit_params_Med2.x0_ES_width(i);
            x0_ES_cs       = init_fit_params_Med2.x0_ES_cs(i);

            x0 = [x0_anion_amp x0_anion_width x0_anion_cs x0_WB_amp x0_WB_width x0_WB_cs x0_ES_amp x0_ES_width x0_ES_cs];

            [FitLorentzians_Med2,resnorm_Med2,exitflag_Med2,output_Med2,...
                FitLorentz1_Med2,FitLorentz2_Med2,FitLorentz3_Med2,...
                pWB_Med2,pAnion_Med2,pES_Med2,pKa_Med2,...
                cs_WB_Med2,cs_Anion_Med2,dcs_Anion_Med2,cs_ES_Med2] = fit3Lorentzians(x_data_Med2, y_data_Med2, x0, x_fit, pH_Med2);


        elseif ismember('fit2LorentziansConst',init_fit_params_Med2.fitFunc(i))

            x0_anion_amp   = init_fit_params_Med2.x0_anion_amp(i);
            x0_anion_width = init_fit_params_Med2.x0_anion_width(i);
            x0_anion_cs    = init_fit_params_Med2.x0_anion_cs(i);
            x0_WB_amp      = init_fit_params_Med2.x0_WB_amp(i);
            x0_WB_width    = init_fit_params_Med2.x0_WB_width(i);
            x0_WB_cs       = init_fit_params_Med2.x0_WB_cs(i);

            x0 = [x0_anion_amp x0_anion_width x0_anion_cs x0_WB_amp x0_WB_width x0_WB_cs];

            if i == 5
                lb = [-inf,-inf,-inf,-inf,x0_WB_width,x0_WB_cs];
                ub = [inf,inf,inf,inf,x0_WB_width,x0_WB_cs];
            end

            [FitLorentzians_Med2,resnorm_Med2,exitflag_Med2,output_Med2,...
                FitLorentz1_Med2,FitLorentz2_Med2,...
                pWB_Med2,pAnion_Med2,pKa_Med2,...
                cs_WB_Med2,cs_Anion_Med2,dcs_Anion_Med2] = fit2LorentziansConst(x_data_Med2, y_data_Med2, x0, lb, ub, x_fit, pH_Med2);

            FitLorentz3_Med2 = [];
            pES_Med2 = [];
            cs_ES_Med2 = [];

        elseif ismember('fit3LorentziansConst',init_fit_params_Med2.fitFunc(i))

            x0_anion_amp   = init_fit_params_Med2.x0_anion_amp(i);
            x0_anion_width = init_fit_params_Med2.x0_anion_width(i);
            x0_anion_cs    = init_fit_params_Med2.x0_anion_cs(i);
            x0_WB_amp      = init_fit_params_Med2.x0_WB_amp(i);
            x0_WB_width    = init_fit_params_Med2.x0_WB_width(i);
            x0_WB_cs       = init_fit_params_Med2.x0_WB_cs(i);
            x0_ES_amp      = init_fit_params_Med2.x0_ES_amp(i);
            x0_ES_width    = init_fit_params_Med2.x0_ES_width(i);
            x0_ES_cs       = init_fit_params_Med2.x0_ES_cs(i);

            x0 = [x0_anion_amp x0_anion_width x0_anion_cs x0_WB_amp x0_WB_width x0_WB_cs x0_ES_amp x0_ES_width x0_ES_cs];

            if i == 4
                lb = [-inf,-inf,-inf,-inf,-inf,-inf,-inf,x0_ES_width,-inf];
                ub = [inf,inf,inf,inf,inf,inf,inf,x0_ES_width,inf];
            elseif i == 8
                lb = [-inf,-inf,x0_anion_cs,-inf,-inf,-inf,-inf,x0_ES_width,x0_ES_cs];
                ub = [inf,inf,x0_anion_cs,inf,inf,inf,inf,x0_ES_width,x0_ES_cs];
            elseif i == 14
                lb = [-inf, -inf, -inf, x0_WB_amp, x0_WB_width, x0_WB_cs, -inf, x0_ES_width, x0_ES_cs];
                ub = [inf, inf, inf, x0_WB_amp, x0_WB_width, x0_WB_cs, inf, x0_ES_width, x0_ES_cs];
            end

            [FitLorentzians_Med2,resnorm_Med2,exitflag_Med2,output_Med2,...
                FitLorentz1_Med2,FitLorentz2_Med2,FitLorentz3_Med2,...
                pWB_Med2,pAnion_Med2,pES_Med2,pKa_Med2,...
                cs_WB_Med2,cs_Anion_Med2,dcs_Anion_Med2,cs_ES_Med2] = fit3LorentziansConst(x_data_Med2, y_data_Med2, x0, lb, ub, x_fit, pH_Med2);

        end

        specs_1D_19F_Med2_fits_params{i,2} = FitLorentz1_Med2;
        specs_1D_19F_Med2_fits_params{i,3} = FitLorentz2_Med2;
        specs_1D_19F_Med2_fits_params{i,4} = FitLorentz3_Med2;

        specs_1D_19F_Med2_fits_params{i,5} = cs_WB_Med2;
        specs_1D_19F_Med2_fits_params{i,6} = cs_Anion_Med2;
        specs_1D_19F_Med2_fits_params{i,7} = cs_ES_Med2;
        specs_1D_19F_Med2_fits_params{i,8} = pWB_Med2;
        specs_1D_19F_Med2_fits_params{i,9} = pAnion_Med2;
        specs_1D_19F_Med2_fits_params{i,10} = pES_Med2;
        specs_1D_19F_Med2_fits_params{i,11} = pKa_Med2;

    end

end



for i = 1:height(specs_1D_19F_tbl)

    if ~isnan(init_fit_params_Med3.pH(i))

        spec_Med3         = cell2mat(specs_1D_19F_tbl.Med3(i));
        pH_Med3           = pHs_1D_19F_Med3(i);
        x_data_Med3_all   = spec_Med3(:,1);
        noise_region_Med3 = spec_Med3(x_data_Med3_all>noise_lims(2)&x_data_Med3_all<noise_lims(1),2);
        noise_mean_Med3   = mean(noise_region_Med3);

        y_data_Med3  = spec_Med3(x_data_Med3_all>(x_min_fits)&x_data_Med3_all<x_max_fits,2) - noise_mean_Med3;
        x_data_Med3  = x_data_Med3_all(x_data_Med3_all>(x_min_fits)&x_data_Med3_all<x_max_fits);

        if ismember('fit1Lorentzian',init_fit_params_Med3.fitFunc(i))

            x0_WB_amp   = init_fit_params_Med3.x0_WB_amp(i);
            x0_WB_width = init_fit_params_Med3.x0_WB_width(i);
            x0_WB_cs    = init_fit_params_Med3.x0_WB_cs(i);

            x0 = [x0_WB_amp x0_WB_width x0_WB_cs];

            [FitLorentzians_Med3,resnorm_Med3,exitflag_Med3,output_Med3,...
                FitLorentz1_Med3,...
                cs_WB_Med3] = fit1Lorentzian(x_data_Med3, y_data_Med3, x0, x_fit, pH_Med3);

            FitLorentz2_Med3 = [];
            FitLorentz3_Med3 = [];
            pWB_Med3 = [];
            pAnion_Med3 = [];
            pES_Med3 = [];
            cs_Anion_Med3 = [];
            cs_ES_Med3 = [];
            pKa_Med3 = [];

        elseif ismember('fit2Lorentzians',init_fit_params_Med3.fitFunc(i))

            x0_anion_amp   = init_fit_params_Med3.x0_anion_amp(i);
            x0_anion_width = init_fit_params_Med3.x0_anion_width(i);
            x0_anion_cs    = init_fit_params_Med3.x0_anion_cs(i);
            x0_WB_amp      = init_fit_params_Med3.x0_WB_amp(i);
            x0_WB_width    = init_fit_params_Med3.x0_WB_width(i);
            x0_WB_cs       = init_fit_params_Med3.x0_WB_cs(i);

            x0 = [x0_anion_amp x0_anion_width x0_anion_cs x0_WB_amp x0_WB_width x0_WB_cs];

            [FitLorentzians_Med3,resnorm_Med3,exitflag_Med3,output_Med3,...
                FitLorentz1_Med3,FitLorentz2_Med3,...
                pWB_Med3,pAnion_Med3,pKa_Med3,...
                cs_WB_Med3,cs_Anion_Med3,dcs_Anion_Med3] = fit2Lorentzians(x_data_Med3, y_data_Med3, x0, x_fit, pH_Med3);

            FitLorentz3_Med3 = [];
            pES_Med3 = [];
            cs_ES_Med3 = [];

        elseif ismember('fit3Lorentzians',init_fit_params_Med3.fitFunc(i))

            x0_anion_amp   = init_fit_params_Med3.x0_anion_amp(i);
            x0_anion_width = init_fit_params_Med3.x0_anion_width(i);
            x0_anion_cs    = init_fit_params_Med3.x0_anion_cs(i);
            x0_WB_amp      = init_fit_params_Med3.x0_WB_amp(i);
            x0_WB_width    = init_fit_params_Med3.x0_WB_width(i);
            x0_WB_cs       = init_fit_params_Med3.x0_WB_cs(i);
            x0_ES_amp      = init_fit_params_Med3.x0_ES_amp(i);
            x0_ES_width    = init_fit_params_Med3.x0_ES_width(i);
            x0_ES_cs       = init_fit_params_Med3.x0_ES_cs(i);

            x0 = [x0_anion_amp x0_anion_width x0_anion_cs x0_WB_amp x0_WB_width x0_WB_cs x0_ES_amp x0_ES_width x0_ES_cs];

            [FitLorentzians_Med3,resnorm_Med3,exitflag_Med3,output_Med3,...
                FitLorentz1_Med3,FitLorentz2_Med3,FitLorentz3_Med3,...
                pWB_Med3,pAnion_Med3,pES_Med3,pKa_Med3,...
                cs_WB_Med3,cs_Anion_Med3,dcs_Anion_Med3,cs_ES_Med3] = fit3Lorentzians(x_data_Med3, y_data_Med3, x0, x_fit, pH_Med3);


        elseif ismember('fit2LorentziansConst',init_fit_params_Med3.fitFunc(i))

            x0_anion_amp   = init_fit_params_Med3.x0_anion_amp(i);
            x0_anion_width = init_fit_params_Med3.x0_anion_width(i);
            x0_anion_cs    = init_fit_params_Med3.x0_anion_cs(i);
            x0_WB_amp      = init_fit_params_Med3.x0_WB_amp(i);
            x0_WB_width    = init_fit_params_Med3.x0_WB_width(i);
            x0_WB_cs       = init_fit_params_Med3.x0_WB_cs(i);

            x0 = [x0_anion_amp x0_anion_width x0_anion_cs x0_WB_amp x0_WB_width x0_WB_cs];

            if i == 5
                lb = [-inf,-inf,-inf,-inf,x0_WB_width,x0_WB_cs];
                ub = [inf,inf,inf,inf,x0_WB_width,x0_WB_cs];
            end

            [FitLorentzians_Med3,resnorm_Med3,exitflag_Med3,output_Med3,...
                FitLorentz1_Med3,FitLorentz2_Med3,...
                pWB_Med3,pAnion_Med3,pKa_Med3,...
                cs_WB_Med3,cs_Anion_Med3,dcs_Anion_Med3] = fit2LorentziansConst(x_data_Med3, y_data_Med3, x0, lb, ub, x_fit, pH_Med3);

            FitLorentz3_Med3 = [];
            pES_Med3 = [];
            cs_ES_Med3 = [];

        elseif ismember('fit3LorentziansConst',init_fit_params_Med3.fitFunc(i))

            x0_anion_amp   = init_fit_params_Med3.x0_anion_amp(i);
            x0_anion_width = init_fit_params_Med3.x0_anion_width(i);
            x0_anion_cs    = init_fit_params_Med3.x0_anion_cs(i);
            x0_WB_amp      = init_fit_params_Med3.x0_WB_amp(i);
            x0_WB_width    = init_fit_params_Med3.x0_WB_width(i);
            x0_WB_cs       = init_fit_params_Med3.x0_WB_cs(i);
            x0_ES_amp      = init_fit_params_Med3.x0_ES_amp(i);
            x0_ES_width    = init_fit_params_Med3.x0_ES_width(i);
            x0_ES_cs       = init_fit_params_Med3.x0_ES_cs(i);

            x0 = [x0_anion_amp x0_anion_width x0_anion_cs x0_WB_amp x0_WB_width x0_WB_cs x0_ES_amp x0_ES_width x0_ES_cs];
            lb = [-inf,-inf,x0_anion_cs,-inf,-inf,-inf,-inf,x0_ES_width,x0_ES_cs];
            ub = [inf,inf,x0_anion_cs,inf,inf,inf,inf,x0_ES_width,x0_ES_cs];

            [FitLorentzians_Med3,resnorm_Med3,exitflag_Med3,output_Med3,...
                FitLorentz1_Med3,FitLorentz2_Med3,FitLorentz3_Med3,...
                pWB_Med3,pAnion_Med3,pES_Med3,pKa_Med3,...
                cs_WB_Med3,cs_Anion_Med3,dcs_Anion_Med3,cs_ES_Med3] = fit3LorentziansConst(x_data_Med3, y_data_Med3, x0, lb, ub, x_fit, pH_Med3);

        end

        specs_1D_19F_Med3_fits_params{i,2} = FitLorentz1_Med3;
        specs_1D_19F_Med3_fits_params{i,3} = FitLorentz2_Med3;
        specs_1D_19F_Med3_fits_params{i,4} = FitLorentz3_Med3;

        specs_1D_19F_Med3_fits_params{i,5} = cs_WB_Med3;
        specs_1D_19F_Med3_fits_params{i,6} = cs_Anion_Med3;
        specs_1D_19F_Med3_fits_params{i,7} = cs_ES_Med3;
        specs_1D_19F_Med3_fits_params{i,8} = pWB_Med3;
        specs_1D_19F_Med3_fits_params{i,9} = pAnion_Med3;
        specs_1D_19F_Med3_fits_params{i,10} = pES_Med3;
        specs_1D_19F_Med3_fits_params{i,11} = pKa_Med3;

    end
end



specs_1D_19F_fits_plots = cell(height(specs_1D_19F_tbl),9);
specs_1D_19F_fits_plots(:,1) = specs_1D_19F_tbl.seq_T;

for i = 1:height(specs_1D_19F_tbl)

    if init_fit_params_Med1.used_in_plot(i) == 1

        specs_1D_19F_fits_plots{i,2}   = pHs_1D_19F_Med1(i);
        specs_1D_19F_fits_plots(i,3:5) = specs_1D_19F_Med1_fits_params(i,2:4);
        specs_1D_19F_fits_plots(i,6:8) = specs_1D_19F_Med1_fits_params(i,8:10);
        specs_1D_19F_fits_plots(i,9)   = specs_1D_19F_tbl.Med1(i);

    elseif init_fit_params_Med1.used_in_plot(i) == 0 && init_fit_params_Med2.used_in_plot(i) == 1

        specs_1D_19F_fits_plots{i,2}   = pHs_1D_19F_Med2(i);
        specs_1D_19F_fits_plots(i,3:5) = specs_1D_19F_Med2_fits_params(i,2:4);
        specs_1D_19F_fits_plots(i,6:8) = specs_1D_19F_Med2_fits_params(i,8:10);
        specs_1D_19F_fits_plots(i,9)   = specs_1D_19F_tbl.Med2(i);

    elseif init_fit_params_Med1.used_in_plot(i) == 0 && init_fit_params_Med2.used_in_plot(i) == 0 && init_fit_params_Med3.used_in_plot(i) == 1

        specs_1D_19F_fits_plots{i,2}   = pHs_1D_19F_Med3(i);
        specs_1D_19F_fits_plots(i,3:5) = specs_1D_19F_Med3_fits_params(i,2:4);
        specs_1D_19F_fits_plots(i,6:8) = specs_1D_19F_Med3_fits_params(i,8:10);
        specs_1D_19F_fits_plots(i,9)   = specs_1D_19F_tbl.Med3(i);

    end

end



%% Figure S6A - 1D 19F spectra at medium pH with fits

vert_offset_19F_1 = 7.8E6;
y_min_19F_all = 1E6;
y_max_19F_all = 6E6;

plot_ind = 0;

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S6A 1",NumberTitle="off");
set(gcf,'Color','w');

for i = 1:4

    plot(specs_1D_19F_fits_plots{i,9}(:,1),...
        (specs_1D_19F_fits_plots{i,9}(:,2).*scales_1D_19F_for_fits(i)) - vert_offset_19F_1*i,...
        'Color',color_G5FdU_medpH,'LineWidth',linewidth);

    hold on;

    if ~isempty(specs_1D_19F_fits_plots{i,5})

        plot(x_fit,specs_1D_19F_fits_plots{i,5}.*scales_1D_19F_for_fits(i) - vert_offset_19F_1*i,...
            '-','Color',color_G5FdU_ES,'LineWidth',linewidth);

    end

    plot(x_fit,specs_1D_19F_fits_plots{i,4}.*scales_1D_19F_for_fits(i) - vert_offset_19F_1*i,...
        '-','Color',color_G5FdU_lowpH,'LineWidth',linewidth);
    plot(x_fit,specs_1D_19F_fits_plots{i,3}.*scales_1D_19F_for_fits(i) - vert_offset_19F_1*i,...
        '-','Color',color_G5FdU_Anion,'LineWidth',linewidth);

    text(-160.3,(1E6-vert_offset_19F_1*i),sprintf('pH %0.1f',specs_1D_19F_fits_plots{i,2}),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);
    text(-165.0,(1E6-vert_offset_19F_1*i),sprintf('%s',specs_1D_19F_fits_plots{i,1}),...
        'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);

    text(-161.8,(5E6-vert_offset_19F_1*i),sprintf('%0.0f%%',specs_1D_19F_fits_plots{i,6}(1)),...
        'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);

    text(-163.9,(5E6-vert_offset_19F_1*i),sprintf('%0.0f%%',specs_1D_19F_fits_plots{i,7}(1)),...
        'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_Anion);

end

hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([x_min_19F x_max_19F (-vert_offset_19F_1*(i+1)+y_max_19F_all) (-vert_offset_19F_1*plot_ind+y_min_19F_all)]);

propedit;

plot_ind = i;

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S6A 2",NumberTitle="off");
set(gcf,'Color','w');

for i = 5:8

    plot(specs_1D_19F_fits_plots{i,9}(:,1),...
        (specs_1D_19F_fits_plots{i,9}(:,2).*scales_1D_19F_for_fits(i)) - vert_offset_19F_1*i,...
        'Color',color_G5FdU_medpH,'LineWidth',linewidth);

    hold on;

    if ~isempty(specs_1D_19F_fits_plots{i,5})

        plot(x_fit,specs_1D_19F_fits_plots{i,5}.*scales_1D_19F_for_fits(i) - vert_offset_19F_1*i,...
            '-','Color',color_G5FdU_ES,'LineWidth',linewidth);

    end

    plot(x_fit,specs_1D_19F_fits_plots{i,4}.*scales_1D_19F_for_fits(i) - vert_offset_19F_1*i,...
        '-','Color',color_G5FdU_lowpH,'LineWidth',linewidth);
    plot(x_fit,specs_1D_19F_fits_plots{i,3}.*scales_1D_19F_for_fits(i) - vert_offset_19F_1*i,...
        '-','Color',color_G5FdU_Anion,'LineWidth',linewidth);

    text(-160.3,(1E6-vert_offset_19F_1*i),sprintf('pH %0.1f',specs_1D_19F_fits_plots{i,2}),...
        'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);
    text(-165.0,(1E6-vert_offset_19F_1*i),sprintf('%s',specs_1D_19F_fits_plots{i,1}),...
        'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);

end

hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([x_min_19F x_max_19F (-vert_offset_19F_1*(i+1)+y_max_19F_all) (-vert_offset_19F_1*plot_ind+y_min_19F_all)]);

propedit;

plot_ind = i;

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S6A 3",NumberTitle="off");
set(gcf,'Color','w');

for i = 9:12

    plot(specs_1D_19F_fits_plots{i,9}(:,1),...
        (specs_1D_19F_fits_plots{i,9}(:,2).*scales_1D_19F_for_fits(i)) - vert_offset_19F_1*i,...
        'Color',color_G5FdU_medpH,'LineWidth',linewidth);

    hold on;

    if ~isempty(specs_1D_19F_fits_plots{i,5})

        plot(x_fit,specs_1D_19F_fits_plots{i,5}.*scales_1D_19F_for_fits(i) - vert_offset_19F_1*i,...
            '-','Color',color_G5FdU_ES,'LineWidth',linewidth);

    end

    plot(x_fit,specs_1D_19F_fits_plots{i,4}.*scales_1D_19F_for_fits(i) - vert_offset_19F_1*i,...
        '-','Color',color_G5FdU_lowpH,'LineWidth',linewidth);
    plot(x_fit,specs_1D_19F_fits_plots{i,3}.*scales_1D_19F_for_fits(i) - vert_offset_19F_1*i,...
        '-','Color',color_G5FdU_Anion,'LineWidth',linewidth);

    text(-160.3,(1E6-vert_offset_19F_1*i),sprintf('pH %0.1f',specs_1D_19F_fits_plots{i,2}),...
        'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);
    text(-165.0,(1E6-vert_offset_19F_1*i),sprintf('%s',specs_1D_19F_fits_plots{i,1}),...
        'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);

end

hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([x_min_19F x_max_19F (-vert_offset_19F_1*(i+1)+y_max_19F_all) (-vert_offset_19F_1*plot_ind+y_min_19F_all)]);

propedit;

plot_ind = i;

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S6A 4",NumberTitle="off");
set(gcf,'Color','w');

for i = 13:16

    plot(specs_1D_19F_fits_plots{i,9}(:,1),...
        (specs_1D_19F_fits_plots{i,9}(:,2).*scales_1D_19F_for_fits(i)) - vert_offset_19F_1*i,...
        'Color',color_G5FdU_medpH,'LineWidth',linewidth);

    hold on;

    if ~isempty(specs_1D_19F_fits_plots{i,5})

        plot(x_fit,specs_1D_19F_fits_plots{i,5}.*scales_1D_19F_for_fits(i) - vert_offset_19F_1*i,...
            '-','Color',color_G5FdU_ES,'LineWidth',linewidth);

    end

    plot(x_fit,specs_1D_19F_fits_plots{i,4}.*scales_1D_19F_for_fits(i) - vert_offset_19F_1*i,...
        '-','Color',color_G5FdU_lowpH,'LineWidth',linewidth);
    plot(x_fit,specs_1D_19F_fits_plots{i,3}.*scales_1D_19F_for_fits(i) - vert_offset_19F_1*i,...
        '-','Color',color_G5FdU_Anion,'LineWidth',linewidth);

    text(-160.3,(1E6-vert_offset_19F_1*i),sprintf('pH %0.1f',specs_1D_19F_fits_plots{i,2}),...
        'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);
    text(-165.0,(1E6-vert_offset_19F_1*i),sprintf('%s',specs_1D_19F_fits_plots{i,1}),...
        'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);

end

hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([x_min_19F x_max_19F (-vert_offset_19F_1*(i+1)+y_max_19F_all) (-vert_offset_19F_1*plot_ind+y_min_19F_all)]);

propedit;

%% variable temperature 

path_1D_19F_VT = 'NMR-Spectra/1D-19F-spectra/VT';
parent_path_1D_19F_VT = sprintf('%s/%s',parent_path,path_1D_19F_VT);

sequence_tri_1D_19F_VT = {'CTC';'ATC'};

pHs_19F_VT = {'Med1'};

temperature_19F_VT = {'1C'; '5C'; '10C'; '15C'; '20C'};

specs_1D_19F_VT = cell(length(sequence_tri_1D_19F_VT),6);

for i = 1:length(sequence_tri_1D_19F_VT)

    specs_1D_19F_VT{i,1} = sequence_tri_1D_19F_VT{i};

    for j = 1:length(temperature_19F_VT)

        filename = sprintf('%s_G5FdU_1D_19F_pH_Med1_%s.txt',sequence_tri_1D_19F_VT{i},temperature_19F_VT{j});

        [axis_19F,spec_19F] = load_1D_spec(parent_path_1D_19F_VT,filename,F2_left_1D_19F,F2_right_1D_19F);

        if ~isempty(axis_19F)

            specs_1D_19F_VT{i,j+1} = [axis_19F,spec_19F];

        end

    end

end

specs_1D_19F_VT_tbl = cell2table(specs_1D_19F_VT,"VariableNames",["seq_T" "1C" "5C" "10C" "15C" "20C"]);

%% fitting 1D 19F VT - CTC

pH_CTC       = 7.97;
specs_CTC    = specs_1D_19F_VT(1,2:end)';

fitparams_CTC = cell(length(temperature_19F_VT),17);
fitparams_CTC(:,1) = temperature_19F_VT;
fitparams_CTC(:,2) = specs_CTC;

for i = 1:length(temperature_19F_VT)

    x_data_all = specs_CTC{i,1}(:,1);
    noise_region = specs_CTC{i,1}(x_data_all>noise_lims(2)&x_data_all<noise_lims(1),2);
    noise_mean   = mean(noise_region);

    y_data_CTC = specs_CTC{i,1}(x_data_all>(x_min_fits)&x_data_all<x_max_fits,2) - noise_mean;
    x_data_CTC = x_data_all(x_data_all>(x_min_fits)&x_data_all<x_max_fits);

    if i == 5

        x0_anion_cs    = -162.6;
        x0_anion_width = 0.3;
        x0_anion_amp   = 1000000;
        x0_WB_cs = -162.9;
        x0_WB_width = 0.3;
        x0_WB_amp = 100000;
        x0_ES_cs = -163.7;
        x0_ES_width = 0.15;
        x0_ES_amp = 170000;

    else
        
        x0_anion_cs    = -162.639;
        x0_anion_width = 0.150;
        x0_anion_amp   = 1770000;
        x0_WB_cs = -163.222;
        x0_WB_width = 0.183;
        x0_WB_amp = 1360000;
        x0_ES_cs = -164.126;
        x0_ES_width = 0.5;
        x0_ES_amp = 800000;

    end

    x0 = [x0_anion_amp x0_anion_width x0_anion_cs x0_WB_amp x0_WB_width x0_WB_cs x0_ES_amp x0_ES_width x0_ES_cs];

    [fitparams_CTC{i,3},fitparams_CTC{i,4},fitparams_CTC{i,5},fitparams_CTC{i,6},...
        fitparams_CTC{i,7},fitparams_CTC{i,8},fitparams_CTC{i,9},...
        fitparams_CTC{i,10},fitparams_CTC{i,11},fitparams_CTC{i,12},fitparams_CTC{i,13},...
        fitparams_CTC{i,14},fitparams_CTC{i,15},fitparams_CTC{i,16},fitparams_CTC{i,17}] = fit3Lorentzians(x_data_CTC, y_data_CTC, x0, x_fit, pH_CTC);

end

pAnion_CTC = [mean([fitparams_CTC{1,11}(1),fitparams_CTC{2,11}(1),fitparams_CTC{3,11}(1)]),mean([fitparams_CTC{1,11}(2),fitparams_CTC{2,11}(2),fitparams_CTC{3,11}(2)])];
pWB_CTC = [mean([fitparams_CTC{1,10}(1),fitparams_CTC{2,10}(1),fitparams_CTC{3,10}(1)]),mean([fitparams_CTC{1,10}(2),fitparams_CTC{2,10}(2),fitparams_CTC{3,10}(2)])];
pES_CTC = [mean([fitparams_CTC{1,12}(1),fitparams_CTC{2,12}(1),fitparams_CTC{3,12}(1)]),mean([fitparams_CTC{1,12}(2),fitparams_CTC{2,12}(2),fitparams_CTC{3,12}(2)])];

pKa_CTC(1,1) = pH_CTC - log10(pAnion_CTC(1,1)./(100-pAnion_CTC(1,1)));
pKa_CTC(1,2) = 1/(log(10))*sqrt((pAnion_CTC(1,2)./pAnion_CTC(1,1)).^2+((pAnion_CTC(1,2))./(100-pAnion_CTC(1,1))).^2);

%% fitting 1D 19F VT - ATC

pH_ATC       = 8.11;
specs_ATC    = specs_1D_19F_VT(2,[2,4])';
temperatures_ATC = temperature_19F_VT([1,3],1);

fitparams_ATC = cell(length(temperatures_ATC),17);
fitparams_ATC(:,1) = temperatures_ATC;
fitparams_ATC(:,2) = specs_ATC;


for i = 1:length(temperatures_ATC)

    x_data_all = specs_ATC{i,1}(:,1);
    noise_region = specs_ATC{i,1}(x_data_all>noise_lims(2)&x_data_all<noise_lims(1),2);
    noise_mean   = mean(noise_region);

    y_data_ATC = specs_ATC{i,1}(x_data_all>(x_min_fits)&x_data_all<x_max_fits,2) - noise_mean;
    x_data_ATC = x_data_all(x_data_all>(x_min_fits)&x_data_all<x_max_fits);

    if i == 1

        x0_anion_cs    = -163.27;
        x0_anion_width = 0.2;
        x0_anion_amp   = 1000000;
        x0_WB_cs = -162.78;
        x0_WB_width = 0.1;
        x0_WB_amp = 1000000;
        x0_ES_cs = -164.6;
        x0_ES_width = 0.1;
        x0_ES_amp = 10000;

    else

        x0_anion_cs    = -163.34;
        x0_anion_width = 0.2;
        x0_anion_amp   = 300000;
        x0_WB_cs = -162.78;
        x0_WB_width = 0.1;
        x0_WB_amp = 1000000;
        x0_ES_cs = -164.5;
        x0_ES_width = 0.5;
        x0_ES_amp = 80000;

    end

    x0 = [x0_anion_amp x0_anion_width x0_anion_cs x0_WB_amp x0_WB_width x0_WB_cs x0_ES_amp x0_ES_width x0_ES_cs];

    [fitparams_ATC{i,3},fitparams_ATC{i,4},fitparams_ATC{i,5},fitparams_ATC{i,6},...
        fitparams_ATC{i,7},fitparams_ATC{i,8},fitparams_ATC{i,9},...
        fitparams_ATC{i,10},fitparams_ATC{i,11},fitparams_ATC{i,12},fitparams_ATC{i,13},...
        fitparams_ATC{i,14},fitparams_ATC{i,15},fitparams_ATC{i,16},fitparams_ATC{i,17}] = fit3Lorentzians(x_data_ATC, y_data_ATC, x0, x_fit, pH_ATC);

end


%% Figure S6C - 1D 19F VT, offset vertically, CTC and ATC

scale_5FdU = 7;

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S6C left",NumberTitle="off");
set(gcf,'Color','w');

for i = 1:length(temperature_19F_VT)

    plot(specs_CTC{i,1}(:,1),(specs_CTC{i,1}(:,2))'.*scale_5FdU-4E7*i,...
        'Color',colors_CTC{i},'LineWidth',linewidth);

    hold on;

    plot(x_fit,fitparams_CTC{i,9}.*scale_5FdU-4E7*i,'-','Color',color_G5FdU_ES,'LineWidth',3);
    plot(x_fit,fitparams_CTC{i,8}.*scale_5FdU-4E7*i,'-','Color',color_G5FdU_lowpH,'LineWidth',3);
    plot(x_fit,fitparams_CTC{i,7}.*scale_5FdU-4E7*i,'-','Color',color_G5FdU_Anion,'LineWidth',3);


end

hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([x_min_19F x_max_19F -2.1E8 9E6]);

propedit;

scale_5FdU = 15;

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S6C right",NumberTitle="off");
set(gcf,'Color','w');

for i = 1:length(temperatures_ATC)

    plot(specs_ATC{i,1}(:,1),(specs_ATC{i,1}(:,2))'.*scale_5FdU-8E7*i,...
        'Color',colors_ATC{i},'LineWidth',linewidth);
    hold on;

    plot(x_fit,fitparams_ATC{i,9}.*scale_5FdU-8E7*i,'-','Color',color_G5FdU_ES,'LineWidth',3);
    plot(x_fit,fitparams_ATC{i,8}.*scale_5FdU-8E7*i,'-','Color',color_G5FdU_lowpH,'LineWidth',3);
    plot(x_fit,fitparams_ATC{i,7}.*scale_5FdU-8E7*i,'-','Color',color_G5FdU_Anion,'LineWidth',3);


end

hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([x_min_19F x_max_19F -2.1E8 9E6]);

propedit;

%% Figure 4A - 1D 19F CTT

[~,ind_CTT] = ismember('CTT',specs_1D_19F_tbl.seq_T);

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig 4A left",NumberTitle="off");
set(gcf,'Color','w');

plot(specs_1D_19F_fits_plots{ind_CTT,9}(:,1),...
    (specs_1D_19F_fits_plots{ind_CTT,9}(:,2).*scales_1D_19F_for_fits(ind_CTT)) - vert_offset_19F_1*ind_CTT,...
    'Color',color_G5FdU_medpH,'LineWidth',linewidth);

hold on;

plot(x_fit,specs_1D_19F_fits_plots{ind_CTT,4}.*scales_1D_19F_for_fits(ind_CTT) - vert_offset_19F_1*ind_CTT,...
    '-','Color',color_G5FdU_lowpH,'LineWidth',linewidth);
plot(x_fit,specs_1D_19F_fits_plots{ind_CTT,3}.*scales_1D_19F_for_fits(ind_CTT) - vert_offset_19F_1*ind_CTT,...
    '-','Color',color_G5FdU_Anion,'LineWidth',linewidth);

text(-161.2,(1E6-vert_offset_19F_1*ind_CTT),sprintf('pH %0.1f',specs_1D_19F_fits_plots{ind_CTT,2}),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);
text(-163.5,(1E6-vert_offset_19F_1*ind_CTT),sprintf('%s',specs_1D_19F_fits_plots{ind_CTT,1}),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);

hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([-164.0 -161.0 -11.3E7 -9.7E7]);

propedit;

%% Figure 4A - 1D 19F TTG


[~,ind_TTG] = ismember('TTG',specs_1D_19F_tbl.seq_T);


figure('units','normalized','outerposition',[0 0 1 1],Name="Fig 4A right 1",NumberTitle="off");
set(gcf,'Color','w');

plot(specs_1D_19F_tbl.Med1{ind_TTG,1}(:,1),...
    (specs_1D_19F_tbl.Med1{ind_TTG,1}(:,2).*scales_1D_19F_for_fits(ind_TTG)) - vert_offset_19F_1*ind_TTG,...
    'Color',color_G5FdU_medpH,'LineWidth',linewidth);

hold on;

plot(x_fit,specs_1D_19F_Med1_fits_params{ind_TTG,2}.*scales_1D_19F_for_fits(ind_TTG) - vert_offset_19F_1*ind_TTG,...
    '-','Color',color_G5FdU_Anion,'LineWidth',linewidth);
plot(x_fit,specs_1D_19F_Med1_fits_params{ind_TTG,3}.*scales_1D_19F_for_fits(ind_TTG) - vert_offset_19F_1*ind_TTG,...
    '-','Color',color_G5FdU_lowpH,'LineWidth',linewidth);

text(-162.25,(1E6-vert_offset_19F_1*ind_TTG),sprintf('pH %0.1f',pHs_1D_19F_all_tbl.Med1(ind_TTG,1)),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);
text(-164.8,(1E6-vert_offset_19F_1*ind_TTG),sprintf('%s',specs_1D_19F_Med1_fits_params{ind_TTG,1}),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);

hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([-165.2 -162.2 -9.95E7 -7.5E7]);

propedit;



figure('units','normalized','outerposition',[0 0 1 1],Name="Fig 4B right 2",NumberTitle="off");
set(gcf,'Color','w');

plot(specs_1D_19F_fits_plots{ind_TTG,9}(:,1),...
    (specs_1D_19F_fits_plots{ind_TTG,9}(:,2).*scales_1D_19F_for_fits(ind_TTG)) - vert_offset_19F_1*ind_TTG,...
    'Color',color_G5FdU_medpH,'LineWidth',linewidth);

hold on;

plot(x_fit,specs_1D_19F_fits_plots{ind_TTG,4}.*scales_1D_19F_for_fits(ind_TTG) - vert_offset_19F_1*ind_TTG,...
    '-','Color',color_G5FdU_lowpH,'LineWidth',linewidth);
plot(x_fit,specs_1D_19F_fits_plots{ind_TTG,3}.*scales_1D_19F_for_fits(ind_TTG) - vert_offset_19F_1*ind_TTG,...
    '-','Color',color_G5FdU_Anion,'LineWidth',linewidth);

text(-162.25,(1E6-vert_offset_19F_1*ind_TTG),sprintf('pH %0.1f',specs_1D_19F_fits_plots{ind_TTG,2}),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);
text(-164.8,(1E6-vert_offset_19F_1*ind_TTG),sprintf('%s',specs_1D_19F_fits_plots{ind_TTG,1}),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);

hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([-165.2 -162.2 -9.8E7 -8E7]);

propedit;





%% Figure S6B - GTA and TTT - pKa at different pHs

[~,ind_GTA] = ismember('GTA',specs_1D_19F_tbl.seq_T);

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S6B left 1",NumberTitle="off");
set(gcf,'Color','w');

plot(specs_1D_19F_fits_plots{ind_GTA,9}(:,1),...
    (specs_1D_19F_fits_plots{ind_GTA,9}(:,2).*scales_1D_19F_for_fits(ind_GTA)) - vert_offset_19F_1*ind_GTA,...
    'Color',color_G5FdU_medpH,'LineWidth',linewidth);

hold on;

plot(x_fit,specs_1D_19F_fits_plots{ind_GTA,4}.*scales_1D_19F_for_fits(ind_GTA) - vert_offset_19F_1*ind_GTA,...
    '-','Color',color_G5FdU_lowpH,'LineWidth',linewidth);
plot(x_fit,specs_1D_19F_fits_plots{ind_GTA,3}.*scales_1D_19F_for_fits(ind_GTA) - vert_offset_19F_1*ind_GTA,...
    '-','Color',color_G5FdU_Anion,'LineWidth',linewidth);

text(-162.1,(1E6-vert_offset_19F_1*ind_GTA),sprintf('pH %0.1f',specs_1D_19F_fits_plots{ind_GTA,2}),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);
text(-164.0,(1E6-vert_offset_19F_1*ind_GTA),sprintf('%s',specs_1D_19F_fits_plots{ind_GTA,1}),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);

text(-162.5,(5E6-vert_offset_19F_1*ind_GTA),sprintf('%0.0f%%',specs_1D_19F_fits_plots{ind_GTA,6}(1)),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);

text(-163.9,(5E6-vert_offset_19F_1*ind_GTA),sprintf('%0.0f%%',specs_1D_19F_fits_plots{ind_GTA,7}(1)),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_Anion);


hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([-164.6 -162.0 -2.7E7 -1.1E7]);

propedit;


figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S6B left 2",NumberTitle="off");
set(gcf,'Color','w');

plot(specs_1D_19F_tbl.Med2{ind_GTA,1}(:,1),...
    (specs_1D_19F_tbl.Med2{ind_GTA,1}(:,2).*scales_1D_19F_for_fits(ind_GTA)) - vert_offset_19F_1*ind_GTA,...
    'Color',color_G5FdU_medpH,'LineWidth',linewidth);

hold on;

plot(x_fit,specs_1D_19F_Med2_fits_params{ind_GTA,3}.*scales_1D_19F_for_fits(ind_GTA) - vert_offset_19F_1*ind_GTA,...
    '-','Color',color_G5FdU_lowpH,'LineWidth',linewidth);
plot(x_fit,specs_1D_19F_Med2_fits_params{ind_GTA,2}.*scales_1D_19F_for_fits(ind_GTA) - vert_offset_19F_1*ind_GTA,...
    '-','Color',color_G5FdU_Anion,'LineWidth',linewidth);

text(-162.1,(1E6-vert_offset_19F_1*ind_GTA),sprintf('pH %0.1f',pHs_1D_19F_all_tbl.Med2(ind_GTA,1)),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);
text(-164.0,(1E6-vert_offset_19F_1*ind_GTA),sprintf('%s',specs_1D_19F_Med2_fits_params{ind_GTA,1}),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);

text(-162.5,(5E6-vert_offset_19F_1*ind_GTA),sprintf('%0.0f%%',specs_1D_19F_Med2_fits_params{ind_GTA,8}(1)),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);

text(-163.9,(5E6-vert_offset_19F_1*ind_GTA),sprintf('%0.0f%%',specs_1D_19F_Med2_fits_params{ind_GTA,9}(1)),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_Anion);


hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([-164.6 -162.0 -2.7E7 -1.1E7]);

propedit;



[~,ind_TTT] = ismember('TTT',specs_1D_19F_tbl.seq_T);

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S6B right 1",NumberTitle="off");
set(gcf,'Color','w');

plot(specs_1D_19F_fits_plots{ind_TTT,9}(:,1),...
    (specs_1D_19F_fits_plots{ind_TTT,9}(:,2).*scales_1D_19F_for_fits(ind_TTT)) - vert_offset_19F_1*ind_TTT,...
    'Color',color_G5FdU_medpH,'LineWidth',linewidth);

hold on;

plot(x_fit,specs_1D_19F_fits_plots{ind_TTT,4}.*scales_1D_19F_for_fits(ind_TTT) - vert_offset_19F_1*ind_TTT,...
    '-','Color',color_G5FdU_lowpH,'LineWidth',linewidth);
plot(x_fit,specs_1D_19F_fits_plots{ind_TTT,3}.*scales_1D_19F_for_fits(ind_TTT) - vert_offset_19F_1*ind_TTT,...
    '-','Color',color_G5FdU_Anion,'LineWidth',linewidth);

text(-161.9,(1E6-vert_offset_19F_1*ind_TTT),sprintf('pH %0.1f',specs_1D_19F_fits_plots{ind_TTT,2}),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);
text(-163.9,(1E6-vert_offset_19F_1*ind_TTT),sprintf('%s',specs_1D_19F_fits_plots{ind_TTT,1}),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);

text(-163.4,(5E6-vert_offset_19F_1*ind_TTT),sprintf('%0.0f%%',specs_1D_19F_fits_plots{ind_TTT,6}(1)),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);

text(-162.5,(5E6-vert_offset_19F_1*ind_TTT),sprintf('%0.0f%%',specs_1D_19F_fits_plots{ind_TTT,7}(1)),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_Anion);


hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([-164.2 -161.8 -8.15E7 -6.6E7]);

propedit;



figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S6B right 2",NumberTitle="off");
set(gcf,'Color','w');

plot(specs_1D_19F_tbl.Med2{ind_TTT,1}(:,1),...
    (specs_1D_19F_tbl.Med2{ind_TTT,1}(:,2).*scales_1D_19F_for_fits(ind_TTT)) - vert_offset_19F_1*ind_TTT,...
    'Color',color_G5FdU_medpH,'LineWidth',linewidth);

hold on;

plot(x_fit,specs_1D_19F_Med2_fits_params{ind_TTT,3}.*scales_1D_19F_for_fits(ind_TTT) - vert_offset_19F_1*ind_TTT,...
    '-','Color',color_G5FdU_lowpH,'LineWidth',linewidth);
plot(x_fit,specs_1D_19F_Med2_fits_params{ind_TTT,2}.*scales_1D_19F_for_fits(ind_TTT) - vert_offset_19F_1*ind_TTT,...
    '-','Color',color_G5FdU_Anion,'LineWidth',linewidth);

text(-161.9,(1E6-vert_offset_19F_1*ind_TTT),sprintf('pH %0.1f',pHs_1D_19F_all_tbl.Med2(ind_TTT,1)),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);
text(-163.9,(1E6-vert_offset_19F_1*ind_TTT),sprintf('%s',specs_1D_19F_Med2_fits_params{ind_TTT,1}),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);

text(-163.4,(5E6-vert_offset_19F_1*ind_TTT),sprintf('%0.0f%%',specs_1D_19F_Med2_fits_params{ind_TTT,8}(1)),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_lowpH);

text(-162.5,(5E6-vert_offset_19F_1*ind_TTT),sprintf('%0.0f%%',specs_1D_19F_Med2_fits_params{ind_TTT,9}(1)),...
    'FontSize',fontsize,'FontWeight','bold','Color',color_G5FdU_Anion);


hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([-164.2 -161.8 -8.35E7 -6E7]);

propedit;



%% CSPs (Figure S2B + Figure 3D)

%% loading data file paths CSPs

path_CSPs        = 'CSPs';
parent_path_CSPs = sprintf('%s/%s',parent_path,path_CSPs);

%parent_path     = '/Users/orsula1/OrsDocs/OrsDocs_since_Jan2025/2025-Sequence-dependence/2025-Data-Deposition/CSPs';

filename_GTA_19F = 'GTA_CSPs.mat';
filename_ATT_19F = 'ATT_CSPs.mat';
filename_CTG_19F = 'CTG_CSPs.mat';
filename_TTG_19F = 'TTG_CSPs.mat';

filepath_GTA_19F  = sprintf('%s/%s',parent_path_CSPs,filename_GTA_19F);
load(filepath_GTA_19F);
filepath_ATT_19F  = sprintf('%s/%s',parent_path_CSPs,filename_ATT_19F);
load(filepath_ATT_19F);
filepath_CTG_19F  = sprintf('%s/%s',parent_path_CSPs,filename_CTG_19F);
load(filepath_CTG_19F);
filepath_TGG_19F  = sprintf('%s/%s',parent_path_CSPs,filename_TTG_19F);
load(filepath_TGG_19F);

seq_con_CSPs = {'GTA';'ATT';'CTG';'TTG'};
atom_name = {'C4p';'H4p'};

neighbor_name = {'5FdU';'5'' of 5FdU';'3'' of 5FdU';...
    'G';'5'' of G';'3'' of G'};

color_C6C8 = '#ff961e';
color_C1p = '#ff463c';
color_C4p = '#ff008c';


%% Figure S2B - CSPs G5FdU vs GT WB 

%% bar plots - GTA (GTA) only

%% GTA CSPs, 5FdU

bar_vals_GTA = [GTA_neigh_G5FdU_vs_GT([2,1,3,5,4,6],:)];

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S4B GTA U 13C",NumberTitle="off");
set(gcf,'Color','w');

hb = bar(bar_vals_GTA(1:3,2:4));

hb(1).FaceColor = color_C6C8;
hb(1).EdgeColor = 'k';
hb(1).LineWidth = 1;

hb(2).FaceColor = color_C1p;
hb(2).EdgeColor = 'k';
hb(2).LineWidth = 1;

hb(3).FaceColor = color_C4p;
hb(3).EdgeColor = 'k';
hb(3).LineWidth = 1;

axis([0.5 3.5 -0.05 0.7]);

set(gca, 'FontSize', fontsize, 'XminorTick', 'off', 'YminorTick', 'on',...
    'LineWidth', axeswidth, 'TickLength', [ticksize,ticksize]);

xticklabels(neighbor_name(1:3));

xlabel('probe name', 'FontSize', 36);
ylabel('GTA 5FdU 13C CSPs [ppm]', 'FontSize', 36);

propedit;


figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S4B GTA U 1H",NumberTitle="off");
set(gcf,'Color','w');

hb = bar(bar_vals_GTA(1:3,5:7));

hb(1).FaceColor = color_C6C8;
hb(1).EdgeColor = 'k';
hb(1).LineWidth = 1;

hb(2).FaceColor = color_C1p;
hb(2).EdgeColor = 'k';
hb(2).LineWidth = 1;

hb(3).FaceColor = color_C4p;
hb(3).EdgeColor = 'k';
hb(3).LineWidth = 1;

axis([0.5 3.5 -0.03 0.5]);

set(gca, 'FontSize', fontsize, 'XminorTick', 'off', 'YminorTick', 'on',...
    'LineWidth', axeswidth, 'TickLength', [ticksize,ticksize]);

xticklabels(neighbor_name(1:3));

xlabel('probe name', 'FontSize', 36);
ylabel('GTA 5FdU 1H CSPs [ppm]', 'FontSize', 36);

propedit;


%% GTA CSPs, G

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S4B GTA G 13C",NumberTitle="off");
set(gcf,'Color','w');

hb = bar(bar_vals_GTA(4:6,2:4));

hb(1).FaceColor = color_C6C8;
hb(1).EdgeColor = 'k';
hb(1).LineWidth = 1;

hb(2).FaceColor = color_C1p;
hb(2).EdgeColor = 'k';
hb(2).LineWidth = 1;

hb(3).FaceColor = color_C4p;
hb(3).EdgeColor = 'k';
hb(3).LineWidth = 1;

axis([0.5 3.5 -0.05 0.12]);

set(gca, 'FontSize', fontsize, 'XminorTick', 'off', 'YminorTick', 'on',...
    'LineWidth', axeswidth, 'TickLength', [ticksize,ticksize]);

xticklabels(neighbor_name(4:6));

xlabel('probe name', 'FontSize', 36);
ylabel('GTA G 13C CSPs [ppm]', 'FontSize', 36);

propedit;


figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S4B GTA G 1H",NumberTitle="off");
set(gcf,'Color','w');

hb = bar(bar_vals_GTA(4:6,5:7));

hb(1).FaceColor = color_C6C8;
hb(1).EdgeColor = 'k';
hb(1).LineWidth = 1;

hb(2).FaceColor = color_C1p;
hb(2).EdgeColor = 'k';
hb(2).LineWidth = 1;

hb(3).FaceColor = color_C4p;
hb(3).EdgeColor = 'k';
hb(3).LineWidth = 1;

axis([0.5 3.5 -0.015 0.015]);

set(gca, 'FontSize', fontsize, 'XminorTick', 'off', 'YminorTick', 'on',...
    'LineWidth', axeswidth, 'TickLength', [ticksize,ticksize]);

xticklabels(neighbor_name(4:6));

xlabel('probe name', 'FontSize', 36);
ylabel('GTA G 1H CSPs [ppm]', 'FontSize', 36);

propedit;



%% bar plots - ATT (ATT) only

%% ATT CSPs, 5FdU

bar_vals_ATT = [ATT_neigh_G5FdU_vs_GT([2,1,3,5,4,6],:)];

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S4B ATT U 13C",NumberTitle="off");
set(gcf,'Color','w');

hb = bar(bar_vals_ATT(1:3,2:4));

hb(1).FaceColor = color_C6C8;
hb(1).EdgeColor = 'k';
hb(1).LineWidth = 1;

hb(2).FaceColor = color_C1p;
hb(2).EdgeColor = 'k';
hb(2).LineWidth = 1;

hb(3).FaceColor = color_C4p;
hb(3).EdgeColor = 'k';
hb(3).LineWidth = 1;

axis([0.5 3.5 -0.05 0.8]);

set(gca, 'FontSize', fontsize, 'XminorTick', 'off', 'YminorTick', 'on',...
    'LineWidth', axeswidth, 'TickLength', [ticksize,ticksize]);

xticklabels(neighbor_name(1:3));

xlabel('probe name', 'FontSize', 36);
ylabel('ATT 5FdU 13C CSPs [ppm]', 'FontSize', 36);

propedit;


figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S4B ATT U 1H",NumberTitle="off");
set(gcf,'Color','w');

hb = bar(bar_vals_ATT(1:3,5:7));

hb(1).FaceColor = color_C6C8;
hb(1).EdgeColor = 'k';
hb(1).LineWidth = 1;

hb(2).FaceColor = color_C1p;
hb(2).EdgeColor = 'k';
hb(2).LineWidth = 1;

hb(3).FaceColor = color_C4p;
hb(3).EdgeColor = 'k';
hb(3).LineWidth = 1;

axis([0.5 3.5 -0.03 0.5]);

set(gca, 'FontSize', fontsize, 'XminorTick', 'off', 'YminorTick', 'on',...
    'LineWidth', axeswidth, 'TickLength', [ticksize,ticksize]);

xticklabels(neighbor_name(1:3));

xlabel('probe name', 'FontSize', 36);
ylabel('ATT 5FdU 1H CSPs [ppm]', 'FontSize', 36);

propedit;


%% ATT CSPs, G

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S4B ATT G 13C",NumberTitle="off");
set(gcf,'Color','w');

hb = bar(bar_vals_ATT(4:6,2:4));

hb(1).FaceColor = color_C6C8;
hb(1).EdgeColor = 'k';
hb(1).LineWidth = 1;

hb(2).FaceColor = color_C1p;
hb(2).EdgeColor = 'k';
hb(2).LineWidth = 1;

hb(3).FaceColor = color_C4p;
hb(3).EdgeColor = 'k';
hb(3).LineWidth = 1;

axis([0.5 3.5 -0.16 0.22]);

set(gca, 'FontSize', fontsize, 'XminorTick', 'off', 'YminorTick', 'on',...
    'LineWidth', axeswidth, 'TickLength', [ticksize,ticksize]);

xticklabels(neighbor_name(4:6));

xlabel('probe name', 'FontSize', 36);
ylabel('ATT G 13C CSPs [ppm]', 'FontSize', 36);

propedit;


figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S4B ATT G 1H",NumberTitle="off");
set(gcf,'Color','w');

hb = bar(bar_vals_ATT(4:6,5:7));

hb(1).FaceColor = color_C6C8;
hb(1).EdgeColor = 'k';
hb(1).LineWidth = 1;

hb(2).FaceColor = color_C1p;
hb(2).EdgeColor = 'k';
hb(2).LineWidth = 1;

hb(3).FaceColor = color_C4p;
hb(3).EdgeColor = 'k';
hb(3).LineWidth = 1;

axis([0.5 3.5 -0.024 0.05]);

set(gca, 'FontSize', fontsize, 'XminorTick', 'off', 'YminorTick', 'on',...
    'LineWidth', axeswidth, 'TickLength', [ticksize,ticksize]);

xticklabels(neighbor_name(4:6));

xlabel('probe name', 'FontSize', 36);
ylabel('ATT G 1H CSPs [ppm]', 'FontSize', 36);

propedit;





%% bar plots - CTG (CTG) only

%% CTG CSPs, 5FdU

bar_vals_CTG = [CTG_neigh_G5FdU_vs_GT([2,1,3,5,4,6],:)];

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S4B CTG U 13C",NumberTitle="off");
set(gcf,'Color','w');

hb = bar(bar_vals_CTG(1:3,2:4));

hb(1).FaceColor = color_C6C8;
hb(1).EdgeColor = 'k';
hb(1).LineWidth = 1;

hb(2).FaceColor = color_C1p;
hb(2).EdgeColor = 'k';
hb(2).LineWidth = 1;

hb(3).FaceColor = color_C4p;
hb(3).EdgeColor = 'k';
hb(3).LineWidth = 1;

axis([0.5 3.5 -0.05 0.8]);

set(gca, 'FontSize', fontsize, 'XminorTick', 'off', 'YminorTick', 'on',...
    'LineWidth', axeswidth, 'TickLength', [ticksize,ticksize]);

xticklabels(neighbor_name(1:3));

xlabel('probe name', 'FontSize', 36);
ylabel('CTG 5FdU 13C CSPs [ppm]', 'FontSize', 36);

propedit;


figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S4B CTG U 1H",NumberTitle="off");
set(gcf,'Color','w');

hb = bar(bar_vals_CTG(1:3,5:7));

hb(1).FaceColor = color_C6C8;
hb(1).EdgeColor = 'k';
hb(1).LineWidth = 1;

hb(2).FaceColor = color_C1p;
hb(2).EdgeColor = 'k';
hb(2).LineWidth = 1;

hb(3).FaceColor = color_C4p;
hb(3).EdgeColor = 'k';
hb(3).LineWidth = 1;

axis([0.5 3.5 -0.03 0.53]);

set(gca, 'FontSize', fontsize, 'XminorTick', 'off', 'YminorTick', 'on',...
    'LineWidth', axeswidth, 'TickLength', [ticksize,ticksize]);

xticklabels(neighbor_name(1:3));

xlabel('probe name', 'FontSize', 36);
ylabel('CTG 5FdU 1H CSPs [ppm]', 'FontSize', 36);

propedit;


%% CTG CSPs, G

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S4B CTG G 13C",NumberTitle="off");
set(gcf,'Color','w');

hb = bar(bar_vals_CTG(4:6,2:4));

hb(1).FaceColor = color_C6C8;
hb(1).EdgeColor = 'k';
hb(1).LineWidth = 1;

hb(2).FaceColor = color_C1p;
hb(2).EdgeColor = 'k';
hb(2).LineWidth = 1;

hb(3).FaceColor = color_C4p;
hb(3).EdgeColor = 'k';
hb(3).LineWidth = 1;

axis([0.5 3.5 -0.04 0.18]);

set(gca, 'FontSize', fontsize, 'XminorTick', 'off', 'YminorTick', 'on',...
    'LineWidth', axeswidth, 'TickLength', [ticksize,ticksize]);

xticklabels(neighbor_name(4:6));

xlabel('probe name', 'FontSize', 36);
ylabel('CTG G 13C CSPs [ppm]', 'FontSize', 36);

propedit;


figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S4B CTG G 1H",NumberTitle="off");
set(gcf,'Color','w');

hb = bar(bar_vals_CTG(4:6,5:7));

hb(1).FaceColor = color_C6C8;
hb(1).EdgeColor = 'k';
hb(1).LineWidth = 1;

hb(2).FaceColor = color_C1p;
hb(2).EdgeColor = 'k';
hb(2).LineWidth = 1;

hb(3).FaceColor = color_C4p;
hb(3).EdgeColor = 'k';
hb(3).LineWidth = 1;

axis([0.5 3.5 -0.012 0.018]);

set(gca, 'FontSize', fontsize, 'XminorTick', 'off', 'YminorTick', 'on',...
    'LineWidth', axeswidth, 'TickLength', [ticksize,ticksize]);

xticklabels(neighbor_name(4:6));

xlabel('probe name', 'FontSize', 36);
ylabel('CTG G 1H CSPs [ppm]', 'FontSize', 36);

propedit;

%% CSPs with mimics

%% Figure 3D, left - C4p CSPs bar plot

bar_vals_C4p_5FdU = [[GTA_neigh_delta_cs{2,4}(2,2),GTA_neigh_delta_cs{5,4}(2,2),GTA_neigh_delta_cs{4,4}(2,2),GTA_neigh_delta_cs{3,4}(2,2)];...
    [ATT_neigh_delta_cs{2,4}(2,2),ATT_neigh_delta_cs{5,4}(2,2),ATT_neigh_delta_cs{4,4}(2,2),ATT_neigh_delta_cs{3,4}(2,2)];...
    [CTG_neigh_delta_cs{2,4}(2,2),CTG_neigh_delta_cs{5,4}(2,2),CTG_neigh_delta_cs{4,4}(2,2),CTG_neigh_delta_cs{3,4}(2,2)];...
    [TTG_neigh_delta_cs{2,4}(2,2),TTG_neigh_delta_cs{5,4}(2,2),TTG_neigh_delta_cs{4,4}(2,2),TTG_neigh_delta_cs{3,4}(2,2)]];

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S3D left",NumberTitle="off");
set(gcf,'Color','w');

hb = bar(bar_vals_C4p_5FdU);

hb(1).FaceColor = color_G5FdU_Anion;
hb(1).EdgeColor = 'k';
hb(1).LineWidth = 1;

hb(2).FaceColor = color_A5FdU_WC;
hb(2).EdgeColor = 'k';
hb(2).LineWidth = 1;

hb(3).FaceColor = color_GC_WC;
hb(3).EdgeColor = 'k';
hb(3).LineWidth = 1;

hb(4).FaceColor = color_isoG5FdU;
hb(4).EdgeColor = 'k';
hb(4).LineWidth = 1;

axis([0.5 4.5 -1.0 2.1]);

set(gca, 'FontSize', fontsize, 'XminorTick', 'off', 'YminorTick', 'on',...
    'LineWidth', axeswidth, 'TickLength', [ticksize,ticksize]);
xticklabels(seq_con_CSPs);

xlabel('Sequence Context', 'FontSize', 36);
ylabel('5FdU delta C4p [ppm]', 'FontSize', 36);

propedit;


%% Figure 3D, right - H4p CSPs bar plot

bar_vals_H4p_5FdU = [[GTA_neigh_delta_cs{2,7}(2,2),GTA_neigh_delta_cs{5,7}(2,2),GTA_neigh_delta_cs{4,7}(2,2),GTA_neigh_delta_cs{3,7}(2,2)];...
    [ATT_neigh_delta_cs{2,7}(2,2),ATT_neigh_delta_cs{5,7}(2,2),ATT_neigh_delta_cs{4,7}(2,2),ATT_neigh_delta_cs{3,7}(2,2)];...
    [CTG_neigh_delta_cs{2,7}(2,2),CTG_neigh_delta_cs{5,7}(2,2),CTG_neigh_delta_cs{4,7}(2,2),CTG_neigh_delta_cs{3,7}(2,2)];...
    [TTG_neigh_delta_cs{2,7}(2,2),TTG_neigh_delta_cs{5,7}(2,2),TTG_neigh_delta_cs{4,7}(2,2),TTG_neigh_delta_cs{3,7}(2,2)]];

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S3D right",NumberTitle="off");
set(gcf,'Color','w');

hb = bar(bar_vals_H4p_5FdU);

hb(1).FaceColor = color_G5FdU_Anion;
hb(1).EdgeColor = 'k';
hb(1).LineWidth = 1;

hb(2).FaceColor = color_A5FdU_WC;
hb(2).EdgeColor = 'k';
hb(2).LineWidth = 1;

hb(3).FaceColor = color_GC_WC;
hb(3).EdgeColor = 'k';
hb(3).LineWidth = 1;

hb(4).FaceColor = color_isoG5FdU;
hb(4).EdgeColor = 'k';
hb(4).LineWidth = 1;

axis([0.5 4.5 0.0 0.35]);

set(gca, 'FontSize', fontsize, 'XminorTick', 'off', 'YminorTick', 'on',...
    'LineWidth', axeswidth, 'TickLength', [ticksize,ticksize]);
xticklabels(seq_con_CSPs);

xlabel('Sequence Context', 'FontSize', 36);
ylabel('5FdU delta H4p [ppm]', 'FontSize', 36);

propedit;


%% 1D imino for mimics (A-5FdU, G-C, isoG-5FdU, Figure S5A)

path_1D_imino_mimics = 'NMR-Spectra/1D-1H-imino-spectra/1H-1D-invWB-and-WC-mimics';
filepath_mimics_1D_imino = sprintf('%s/%s',parent_path,path_1D_imino_mimics);

mimics = {'isoG5FdU';'A5FdU_WC';'GC_WC'};


%% loading all 1H spectra of mimics

specs_mimics_1D_imino = cell(length(seq_con_CSPs),length(mimics)+1);

for i = 1:length(seq_con_CSPs)

    specs_mimics_1D_imino{i,1} = seq_con_CSPs{i};

    for j = 1:length(mimics)

        filename = sprintf('%s_%s_1D_imino_pH_Low.txt',seq_con_CSPs{i},mimics{j});

        [axis_1H,spec_1H] = load_1D_spec(filepath_mimics_1D_imino,filename,F2_left_1D_imino,F2_right_1D_imino);

        if ~isempty(axis_1H)

            specs_mimics_1D_imino{i,j+1} = [axis_1H,spec_1H];

        end

    end

end

specs_mimics_1D_imino_tbl = cell2table(specs_mimics_1D_imino,"VariableNames",["seq_T" "isoG5FdU" "A5FdU_WC" "GC_WC"]);



scale_5FdU = 6;

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S5A top left",NumberTitle="off");
set(gcf,'Color','w');

for i = 1:4
    
    if i == 1

        plot(specs_mimics_1D_imino{i,2}(:,1),(specs_mimics_1D_imino{i,2}(:,2))'.*scale_5FdU*0.8-2.5E8*i,...
            'Color',color_isoG5FdU,'LineWidth',linewidth);
        hold on;

    else

        plot(specs_mimics_1D_imino{i,2}(:,1),(specs_mimics_1D_imino{i,2}(:,2))'.*scale_5FdU-2.5E8*i,...
            'Color',color_isoG5FdU,'LineWidth',linewidth);
        hold on;

    end
    
end

hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([9.0 14.5 -2.5E8*(i+1)+2E8 0.5E8]);

propedit;


scale_5FdU = 5.5;

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S5A top right",NumberTitle="off");
set(gcf,'Color','w');

for i = 1:4
    
    if i == 1

        plot(specs_mimics_1D_imino{i,4}(:,1),(specs_mimics_1D_imino{i,4}(:,2))'.*scale_5FdU*0.3-2.5E8*i,...
            'Color',color_GC_WC,'LineWidth',linewidth);
        hold on;

    elseif i == 3

        plot(specs_mimics_1D_imino{i,4}(:,1),(specs_mimics_1D_imino{i,4}(:,2))'.*scale_5FdU*0.3-2.5E8*i,...
            'Color',color_GC_WC,'LineWidth',linewidth);
        hold on;

    else

        plot(specs_mimics_1D_imino{i,4}(:,1),(specs_mimics_1D_imino{i,4}(:,2))'.*scale_5FdU*2.1-2.5E8*i,...
            'Color',color_GC_WC,'LineWidth',linewidth);
        hold on;

    end

    
end

hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([9.0 14.5 -2.5E8*(i+1)+2E8 0.5E8]);

propedit;


scale_5FdU = 6;

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S5A bottom left",NumberTitle="off");
set(gcf,'Color','w');

for i = 1:4
    
    if i == 1

        plot(specs_mimics_1D_imino{i,3}(:,1),(specs_mimics_1D_imino{i,3}(:,2))'.*scale_5FdU*0.4-2.5E8*i,...
            'Color',color_A5FdU_WC,'LineWidth',linewidth);
        hold on;

    else

        plot(specs_mimics_1D_imino{i,3}(:,1),(specs_mimics_1D_imino{i,3}(:,2))'.*scale_5FdU*1.2-2.5E8*i,...
            'Color',color_A5FdU_WC,'LineWidth',linewidth);
        hold on;

    end

    
end

hold off;

set(gca, 'FontSize', fontsize, 'Xdir', 'rev', 'YTick',[],...
    'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

axis([9.0 14.5 -2.5E8*(i+1)+2E8 0.5E8]);

propedit;


%% Figure S7 - pAnion and pTautomer vs. sequence

path_Taut_pops = 'Previous-R1rho-Params/Prev-Tautomer-Populations';
path_Taut_10C = sprintf('%s/%s',parent_path,path_Taut_pops);

filename_Taut_10C = sprintf('%s/Taut_pops_and_kex_10C.mat',path_Taut_10C);
load(filename_Taut_10C);

pKas_diff_average = 2.705399674455320;
pKas_diff_stdv    = 0.209225409074581;

pKas_adj_from_19F      = sorted_pKas_19F + pKas_diff_average;
pKas_adj_from_19F_errs = sqrt(sorted_pKas_19F_errs.^2 + pKas_diff_stdv.^2);


pH_for_calc = 7.6;

pAnion_adj_from_19F_pH_7p6          = (10.^(pH_for_calc-pKas_adj_from_19F))./(1+10.^(pH_for_calc-pKas_adj_from_19F));
pAnion_adj_from_19F_pH_7p6_errs     = pAnion_adj_from_19F_pH_7p6.*sqrt((pKas_adj_from_19F_errs).^2.*(log(10)).^2+((10.^(pH_for_calc-pKas_adj_from_19F)).^2.*(pKas_adj_from_19F_errs).^2.*(log(10))^2)./((1+10.^(pH_for_calc-pKas_adj_from_19F)).^2));

pTaut_min = min(Taut_pop_all_10C);
pTaut_max = max(Taut_pop_all_10C);

x_rec = 0:0.2:18;
y_Taut_min = pTaut_min.*(ones(length(x_rec)));
y_Taut_max = pTaut_max.*(ones(length(x_rec)));


ind_CTC_Taut_DNA_High_pH = find(strcmp('CTC',Taut_prev_DNA_High_pH_10C.seq_prev_DNA_High_pH_10C_T));
ind_GTC_Taut_DNA_High_pH = find(strcmp('GTC',Taut_prev_DNA_High_pH_10C.seq_prev_DNA_High_pH_10C_T));
ind_TTG_Taut_DNA_High_pH = find(strcmp('TTG',Taut_prev_DNA_High_pH_10C.seq_prev_DNA_High_pH_10C_T));
ind_GTG_Taut_DNA_High_pH = find(strcmp('GTG',Taut_prev_DNA_High_pH_10C.seq_prev_DNA_High_pH_10C_T));

ind_CTC_Taut_DNA_6p9 = find(strcmp('CTC',Taut_DNA_6p9_10C.seq_DNA_6p9_10C_T));
ind_GTC_Taut_DNA_6p9 = find(strcmp('GTC',Taut_DNA_6p9_10C.seq_DNA_6p9_10C_T));
ind_TTA_Taut_DNA_6p9 = find(strcmp('TTA',Taut_DNA_6p9_10C.seq_DNA_6p9_10C_T));
ind_GTG_Taut_DNA_6p9 = find(strcmp('GTG',Taut_DNA_6p9_10C.seq_DNA_6p9_10C_T));


ind_CTG_Taut_RNA_High_pH = find(strcmp('CTG',Taut_prev_RNA_High_pH_10C.seq_prev_RNA_High_pH_10C_T));
ind_GTG_Taut_RNA_High_pH = find(strcmp('GTG',Taut_prev_RNA_High_pH_10C.seq_prev_RNA_High_pH_10C_T));

ind_CTT_Taut_RNA_6p9 = find(strcmp('CTT',Taut_RNA_6p9_10C.seq_RNA_6p9_10C_T));
ind_GTG_Taut_RNA_6p9 = find(strcmp('GTG',Taut_RNA_6p9_10C.seq_RNA_6p9_10C_T));


pB_Taut_per_CTC_10C     = mean([Taut_DNA_6p9_10C.pTaut_DNA_6p9_10C(ind_CTC_Taut_DNA_6p9);Taut_prev_DNA_High_pH_10C.pTaut_prev_DNA_High_pH_10C(ind_CTC_Taut_DNA_High_pH)]);
pB_Taut_per_CTC_10C_std = std([Taut_DNA_6p9_10C.pTaut_DNA_6p9_10C(ind_CTC_Taut_DNA_6p9);Taut_prev_DNA_High_pH_10C.pTaut_prev_DNA_High_pH_10C(ind_CTC_Taut_DNA_High_pH)]);

pB_Taut_per_GTC_10C     = mean([Taut_DNA_6p9_10C.pTaut_DNA_6p9_10C(ind_GTC_Taut_DNA_6p9);Taut_prev_DNA_High_pH_10C.pTaut_prev_DNA_High_pH_10C(ind_GTC_Taut_DNA_High_pH)]);
pB_Taut_per_GTC_10C_std = std([Taut_DNA_6p9_10C.pTaut_DNA_6p9_10C(ind_GTC_Taut_DNA_6p9);Taut_prev_DNA_High_pH_10C.pTaut_prev_DNA_High_pH_10C(ind_GTC_Taut_DNA_High_pH)]);

pB_Taut_per_TTA_10C     = Taut_DNA_6p9_10C.pTaut_DNA_6p9_10C(ind_TTA_Taut_DNA_6p9);
pB_Taut_per_TTA_10C_std = Taut_DNA_6p9_10C.pTaut_DNA_6p9_10C_errs(ind_TTA_Taut_DNA_6p9);

pB_Taut_per_TTG_10C     = Taut_prev_DNA_High_pH_10C.pTaut_prev_DNA_High_pH_10C(ind_TTG_Taut_DNA_High_pH);
pB_Taut_per_TTG_10C_std = Taut_prev_DNA_High_pH_10C.pTaut_prev_DNA_High_pH_10C_errs(ind_TTG_Taut_DNA_High_pH);

pB_Taut_per_GTG_10C     = mean([Taut_DNA_6p9_10C.pTaut_DNA_6p9_10C(ind_GTG_Taut_DNA_6p9);Taut_prev_DNA_High_pH_10C.pTaut_prev_DNA_High_pH_10C(ind_GTG_Taut_DNA_High_pH);...
    Taut_RNA_6p9_10C.pTaut_RNA_6p9_10C(ind_GTG_Taut_RNA_6p9);Taut_prev_RNA_High_pH_10C.pTaut_prev_RNA_High_pH_10C(ind_GTG_Taut_RNA_High_pH)]);
pB_Taut_per_GTG_10C_std = std([Taut_DNA_6p9_10C.pTaut_DNA_6p9_10C(ind_GTG_Taut_DNA_6p9);Taut_prev_DNA_High_pH_10C.pTaut_prev_DNA_High_pH_10C(ind_GTG_Taut_DNA_High_pH);...
    Taut_RNA_6p9_10C.pTaut_RNA_6p9_10C(ind_GTG_Taut_RNA_6p9);Taut_prev_RNA_High_pH_10C.pTaut_prev_RNA_High_pH_10C(ind_GTG_Taut_RNA_High_pH)]);

pB_Taut_per_CTG_10C     = Taut_prev_RNA_High_pH_10C.pTaut_prev_RNA_High_pH_10C(ind_CTG_Taut_RNA_High_pH);
pB_Taut_per_CTG_10C_std = Taut_prev_RNA_High_pH_10C.pTaut_prev_RNA_High_pH_10C_errs(ind_CTG_Taut_RNA_High_pH);

pB_Taut_per_CTT_10C     = Taut_RNA_6p9_10C.pTaut_RNA_6p9_10C(ind_CTT_Taut_RNA_6p9);
pB_Taut_per_CTT_10C_std = Taut_RNA_6p9_10C.pTaut_RNA_6p9_10C(ind_CTT_Taut_RNA_6p9);

seqs_for_taut_plot_10C = {'CTC';'GTC';'TTA';'TTG';'GTG';'CTG';'CTT'};

[~,ind_10C_for_plot] = ismember(seqs_for_taut_plot_10C,sorted_seqs_T_19F);

pB_Taut_per_for_plot_10C      = [pB_Taut_per_CTC_10C;pB_Taut_per_GTC_10C;...
    pB_Taut_per_TTA_10C;pB_Taut_per_TTG_10C;pB_Taut_per_GTG_10C;...
    pB_Taut_per_CTG_10C;pB_Taut_per_CTT_10C];
pB_Taut_per_for_plot_10C_errs = [pB_Taut_per_CTC_10C_std;pB_Taut_per_GTC_10C_std;...
    pB_Taut_per_TTA_10C_std;pB_Taut_per_TTG_10C_std;pB_Taut_per_GTG_10C_std;...
    pB_Taut_per_CTG_10C_std;pB_Taut_per_CTT_10C_std];


path_Anion_pops = 'Previous-R1rho-Params/Prev-Anion-Populations';
path_interp_anion = sprintf('%s/%s',parent_path,path_Anion_pops);

filename_interp_anion   = 'DNA_interp_to_7p6_10C.csv';
interp_orig_DNA_7p6_10C = readtable(sprintf('%s/%s',path_interp_anion,filename_interp_anion));

pAnion_prev_DNA_7p6_10C       = 100.*interp_orig_DNA_7p6_10C.pA_pH;
pAnion_prev_DNA_7p6_10C_errs  = 100.*interp_orig_DNA_7p6_10C.pA_pH_err;
pAnion_prev_DNA_7p6_10C_seq_T = interp_orig_DNA_7p6_10C.seq;

[~,ind_DNA_7p6_10C_for_plot] = ismember(pAnion_prev_DNA_7p6_10C_seq_T,sorted_seqs_T_19F);

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig S7",NumberTitle="off");
set(gcf,'Color','w');

errorbar(1:16,100.*pAnion_adj_from_19F_pH_7p6,100.*pAnion_adj_from_19F_pH_7p6_errs,...
    'o','MarkerEdgeColor',color_G5FdU_Anion,'MarkerFaceColor',color_G5FdU_Anion,...
    'MarkerSize',markersize_pKa,...
    'Color',color_G5FdU_Anion,'CapSize',err_capsize,'Linewidth',axeswidth);

hold on;

rectangle('Position',[0 pTaut_min 17 (pTaut_max - pTaut_min)],...
    'FaceColor',color_G5FdU_ES,'FaceAlpha',0.2,'EdgeColor','k','LineStyle','none');

errorbar(ind_10C_for_plot,pB_Taut_per_for_plot_10C,pB_Taut_per_for_plot_10C_errs,...
    'o','MarkerEdgeColor',color_G5FdU_ES,...
    'MarkerSize',markersize_pKa,...
    'Color',color_G5FdU_ES,'CapSize',err_capsize,'Linewidth',axeswidth);

errorbar(ind_DNA_7p6_10C_for_plot,pAnion_prev_DNA_7p6_10C,pAnion_prev_DNA_7p6_10C_errs,...
    'o','MarkerEdgeColor',color_G5FdU_medpH,...
    'MarkerSize',markersize_pKa,...
    'Color',color_G5FdU_medpH,'CapSize',err_capsize,'Linewidth',axeswidth);

hold off;

axis([0.0 17.0 0.0 0.3]);

set(gca, 'FontSize', fontsize, 'XTick', 1:17, 'XminorTick', 'off', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

set(gca, 'YScale', 'linear');

xticklabels(sorted_seqs_T_19F);

xlabel('Sequence', 'FontSize', 36);
ylabel('pB [%]', 'FontSize', 36);

propedit;


%% Figure 5 - ddG 19F vs. R1rho

R_constant = 1.9872036E-3; %kCal/molK
Temp_25C   = 298.15; %K
Temp_10C   = 283.15; %K
Temp_4C    = 277.15; %K
Temp_1C    = 274.15; %K

path_dG = 'Previous-R1rho-Params/Prev-dG-Anion-and-Tautomer';
path_dG_R1rho = sprintf('%s/%s',parent_path,path_dG);
filename_dG_R1rho = sprintf('%s/dG_from_R1rho.mat',path_dG_R1rho);

load(filename_dG_R1rho);

seq_Anion_prev_T               = params_prev_High_pH_10C_DNA.seq_prev_High_pH_10C_DNA_T;
pKas_prev_High_pH_10C_DNA      = params_prev_High_pH_10C_DNA.pKas_prev_High_pH_10C_DNA;
pKas_prev_High_pH_10C_DNA_errs = params_prev_High_pH_10C_DNA.pKas_prev_High_pH_10C_DNA_errs;

pH_for_calc = 7.4;

pAnion_19F_pH_7p4        = (10.^(pH_for_calc-sorted_pKas_19F))./(1+10.^(pH_for_calc-sorted_pKas_19F));
pAnion_19F_pH_7p4_errs   = pAnion_19F_pH_7p4.*sqrt((sorted_pKas_19F_errs).^2.*(log(10)).^2+((10.^(pH_for_calc-sorted_pKas_19F)).^2.*(sorted_pKas_19F_errs).^2.*(log(10))^2)./((1+10.^(pH_for_calc-sorted_pKas_19F)).^2));

dG_Anion_19F_7p4_1C      = -R_constant*Temp_1C*log(10)*(pH_for_calc - sorted_pKas_19F);
dG_Anion_19F_7p4_1C_errs = sorted_pKas_19F_errs;

ref_seq_dG = 'CTC';

ddG_Anion_19F_7p4_1C      = dG_Anion_19F_7p4_1C-dG_Anion_19F_7p4_1C(strcmp(sorted_seqs_T_19F,ref_seq_dG));
ddG_Anion_19F_7p4_1C_errs = sqrt(dG_Anion_19F_7p4_1C_errs.^2+(dG_Anion_19F_7p4_1C_errs(strcmp(sorted_seqs_T_19F,ref_seq_dG))).^2);

ddG_Anion_prev_7p4_10C      = dG_Anion_prev_7p4_10C-dG_Anion_prev_7p4_10C(strcmp(seq_Anion_prev_T,ref_seq_dG));
ddG_Anion_prev_7p4_10C_errs = sqrt(dG_Anion_prev_7p4_10C_errs.^2+(dG_Anion_prev_7p4_10C_errs(strcmp(seq_Anion_prev_T,ref_seq_dG))).^2);

[~,ind_dG_R1rho] = ismember(seq_Anion_prev_T,sorted_seqs_T_19F);

%% GTA from R1rho measurement

seq_for_R1rho = 'GTA';
pH_for_calc   = 8.52;

pAnion_GTA_R1rho_pH8p52_per = 0.062;
pAnion_GTA_R1rho_pH8p52_per_err = 0.009;

pKa_GTA_R1rho_4C     = pH_for_calc - log10(pAnion_GTA_R1rho_pH8p52_per/(100-pAnion_GTA_R1rho_pH8p52_per));
pKa_GTA_R1rho_4C_err = (1/log(10))*sqrt(((pAnion_GTA_R1rho_pH8p52_per_err)./(pAnion_GTA_R1rho_pH8p52_per)).^2+(pAnion_GTA_R1rho_pH8p52_per_err).^2./((100-pAnion_GTA_R1rho_pH8p52_per).^2));

dG_Anion_GTA_R1rho_pH8p52_4C      = -R_constant*Temp_4C*log(10)*(pH_for_calc - pKa_GTA_R1rho_4C);
dG_Anion_GTA_R1rho_pH8p52_4C_err  = pKa_GTA_R1rho_4C_err;

dG_Anion_ref_seq_pH8p52_10C       = -R_constant*Temp_10C*log(10)*(pH_for_calc - pKas_prev_High_pH_10C_DNA(strcmp(seq_Anion_prev_T,ref_seq_dG)));
dG_Anion_ref_seq_pH8p52_10C_err   = pKas_prev_High_pH_10C_DNA_errs(strcmp(seq_Anion_prev_T,ref_seq_dG));

ddG_Anion_GTA_R1rhp_pH8p52_4C     = dG_Anion_GTA_R1rho_pH8p52_4C - dG_Anion_ref_seq_pH8p52_10C;
ddG_Anion_GTA_R1rhp_pH8p52_4C_err = sqrt((dG_Anion_GTA_R1rho_pH8p52_4C_err)^2+(dG_Anion_ref_seq_pH8p52_10C_err)^2);

ddG_Anion_GTA_19F_7p4_1C          = ddG_Anion_19F_7p4_1C(strcmp(sorted_seqs_T_19F,seq_for_R1rho));
ddG_Anion_GTA_19F_7p4_1C_err      = ddG_Anion_19F_7p4_1C_errs(strcmp(sorted_seqs_T_19F,seq_for_R1rho));

%% Figure 5A - ddG 19F vs. R1rho

x_data      = ddG_Anion_prev_7p4_10C;
x_data_errs = ddG_Anion_prev_7p4_10C_errs;
y_data      = ddG_Anion_19F_7p4_1C(ind_dG_R1rho);
y_data_errs = ddG_Anion_19F_7p4_1C_errs(ind_dG_R1rho);

[ddG_Anion_19F_vs_R1rho_r,ddG_Anion_19F_vs_R1rho_Rsq,ddG_Anion_19F_vs_R1rho_Rsq_det,...
    ddG_Anion_19F_vs_R1rho_rms,ddG_Anion_19F_vs_R1rho_rms_per,...
    ddG_Anion_19F_vs_R1rho_bestfit,x2_ddG_R1rho,inBetween_ddG_R1rho,...
    ddG_Anion_19F_vs_R1rho_slope,ddG_Anion_19F_vs_R1rho_inter,...
    ddG_Anion_19F_vs_R1rho_slope_err,ddG_Anion_19F_vs_R1rho_inter_err] = corr_plot(x_data,y_data,x_line,conf_int_per);

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig 5A",NumberTitle="off");
set(gcf,'Color','w');

errorbar(x_data,y_data,y_data_errs,y_data_errs,...
    x_data_errs,x_data_errs,...
    'o','MarkerEdgeColor',color_GT,'MarkerFaceColor',color_GT,'MarkerSize',markersize_ddG,...
    'Color',color_GT,'CapSize',err_capsize,'Linewidth',axeswidth);

hold on;

plot(x_line,y_line,'k-','LineWidth',2);

plot(x_line,ddG_Anion_19F_vs_R1rho_bestfit,'k--','LineWidth',2);

errorbar(ddG_Anion_GTA_R1rhp_pH8p52_4C,ddG_Anion_GTA_19F_7p4_1C,...
    ddG_Anion_GTA_19F_7p4_1C_err,ddG_Anion_GTA_19F_7p4_1C_err,...
    ddG_Anion_GTA_R1rhp_pH8p52_4C_err,ddG_Anion_GTA_R1rhp_pH8p52_4C_err,...
    'o','MarkerEdgeColor','r','MarkerFaceColor','r',...
    'MarkerSize',markersize_ddG,'Color','r','CapSize',err_capsize,'Linewidth',2);

hold off;

axis([-0.5 2.5 -0.5 2.5]);

set(gca, 'FontSize', fontsize, 'XminorTick', 'on', 'YminorTick', 'on',...
    'XTick',-1.0:0.5:2.0,'YTick',-1.0:0.5:2.0,'LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

text(-0.4,2.2,sprintf('RMSD = %0.1f kcal/mol',ddG_Anion_19F_vs_R1rho_rms),'FontSize',fontsize,'FontWeight','bold');
text(-0.4,2.0,sprintf('r = %0.2f',ddG_Anion_19F_vs_R1rho_r),'FontSize',fontsize,'FontWeight','bold');

xlabel('DDG R1rho', 'FontSize', 36);
ylabel('DDG 19F', 'FontSize', 36);

propedit;



%% loading R1rho data

path_GTA_R1rho = 'R1rho-GTA/R1rho-Figure';
parent_path_GTA_R1rho = sprintf('%s/%s',parent_path,path_GTA_R1rho);

filename     = 'Data_GN1_1.csv';
filepath_GN1 = sprintf('%s/%s',parent_path_GTA_R1rho,filename);
R1rho_GN1    = readtable(filepath_GN1);

filename     = 'Data_TN3_1.csv';
filepath_TN3 = sprintf('%s/%s',parent_path_GTA_R1rho,filename);
R1rho_TN3    = readtable(filepath_TN3);

SLPs_GN1 = [500,800,1200,1500,2000];
SLPs_TN3 = [500,800,1200,1500,2000,2500];

colors_for_R1rho_plot = {'#2d55aa','#9ebd5d','#ff9600','#a52a2a','#800080','#cc6677'};

%% loading simulations

filename     = 'GN1-3state-triangle-sim/sim-r1p.csv';
filepath_GN1 = sprintf('%s/%s',parent_path_GTA_R1rho,filename);
sim_GN1      = readtable(filepath_GN1);

filename     = 'TN3-3state-triangle-sim/sim-r1p.csv';
filepath_TN3 = sprintf('%s/%s',parent_path_GTA_R1rho,filename);
sim_TN3      = readtable(filepath_TN3);


%% R1rho GTA pH 8.52 G(N1) - Figure 5B left

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig 5B left",NumberTitle="off");
set(gcf,'Color','w');

errorbar(R1rho_GN1.Offset(R1rho_GN1.SLP==SLPs_GN1(1))./1000,...
    R1rho_GN1.R2eff(R1rho_GN1.SLP==SLPs_GN1(1)),...
    R1rho_GN1.R2effErr(R1rho_GN1.SLP==SLPs_GN1(1)),...
    'o','MarkerEdgeColor',colors_for_R1rho_plot{1},'MarkerFaceColor',colors_for_R1rho_plot{1},...
    'MarkerSize',markersize,'Color',colors_for_R1rho_plot{1},'CapSize',err_capsize,'Linewidth',err_linewidth);

hold on;

plot(sim_GN1.offset(sim_GN1.slp==SLPs_GN1(1))./1000,...
    sim_GN1.r2eff(sim_GN1.slp==SLPs_GN1(1)),...
    '-','Color',colors_for_R1rho_plot{1},'LineWidth',linewidth);

errorbar(R1rho_GN1.Offset(R1rho_GN1.SLP==SLPs_GN1(2))./1000,...
    R1rho_GN1.R2eff(R1rho_GN1.SLP==SLPs_GN1(2)),...
    R1rho_GN1.R2effErr(R1rho_GN1.SLP==SLPs_GN1(2)),...
    'o','MarkerEdgeColor',colors_for_R1rho_plot{2},'MarkerFaceColor',colors_for_R1rho_plot{2},...
    'MarkerSize',markersize,'Color',colors_for_R1rho_plot{2},'CapSize',err_capsize,'Linewidth',err_linewidth);

plot(sim_GN1.offset(sim_GN1.slp==SLPs_GN1(2))./1000,...
    sim_GN1.r2eff(sim_GN1.slp==SLPs_GN1(2)),...
    '-','Color',colors_for_R1rho_plot{2},'LineWidth',linewidth);

errorbar(R1rho_GN1.Offset(R1rho_GN1.SLP==SLPs_GN1(3))./1000,...
    R1rho_GN1.R2eff(R1rho_GN1.SLP==SLPs_GN1(3)),...
    R1rho_GN1.R2effErr(R1rho_GN1.SLP==SLPs_GN1(3)),...
    'o','MarkerEdgeColor',colors_for_R1rho_plot{3},'MarkerFaceColor',colors_for_R1rho_plot{3},...
    'MarkerSize',markersize,'Color',colors_for_R1rho_plot{3},'CapSize',err_capsize,'Linewidth',err_linewidth);

plot(sim_GN1.offset(sim_GN1.slp==SLPs_GN1(3))./1000,...
    sim_GN1.r2eff(sim_GN1.slp==SLPs_GN1(3)),...
    '-','Color',colors_for_R1rho_plot{3},'LineWidth',linewidth);

errorbar(R1rho_GN1.Offset(R1rho_GN1.SLP==SLPs_GN1(4))./1000,...
    R1rho_GN1.R2eff(R1rho_GN1.SLP==SLPs_GN1(4)),...
    R1rho_GN1.R2effErr(R1rho_GN1.SLP==SLPs_GN1(4)),...
    'o','MarkerEdgeColor',colors_for_R1rho_plot{4},'MarkerFaceColor',colors_for_R1rho_plot{4},...
    'MarkerSize',markersize,'Color',colors_for_R1rho_plot{4},'CapSize',err_capsize,'Linewidth',err_linewidth);

plot(sim_GN1.offset(sim_GN1.slp==SLPs_GN1(4))./1000,...
    sim_GN1.r2eff(sim_GN1.slp==SLPs_GN1(4)),...
    '-','Color',colors_for_R1rho_plot{4},'LineWidth',linewidth);

errorbar(R1rho_GN1.Offset(R1rho_GN1.SLP==SLPs_GN1(5))./1000,...
    R1rho_GN1.R2eff(R1rho_GN1.SLP==SLPs_GN1(5)),...
    R1rho_GN1.R2effErr(R1rho_GN1.SLP==SLPs_GN1(5)),...
    'o','MarkerEdgeColor',colors_for_R1rho_plot{5},'MarkerFaceColor',colors_for_R1rho_plot{5},...
    'MarkerSize',markersize,'Color',colors_for_R1rho_plot{5},'CapSize',err_capsize,'Linewidth',err_linewidth);

plot(sim_GN1.offset(sim_GN1.slp==SLPs_GN1(5))./1000,...
    sim_GN1.r2eff(sim_GN1.slp==SLPs_GN1(5)),...
    '-','Color',colors_for_R1rho_plot{5},'LineWidth',linewidth);

hold off;

set(gca, 'FontSize', fontsize,...
    'XminorTick', 'off', 'YminorTick', 'off',...
    'LineWidth',axeswidth,'TickLength',[ticksize,ticksize]);

axis([-6.200 6.200 6 41]);

xlabel('Offset (kHz)', 'FontSize', 36);
ylabel('R2+Rex (s-1)', 'FontSize', 36);

propedit;


%% R1rho GTA pH 8.52 T(N3) - Figure 5B right

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig 5B right",NumberTitle="off");
set(gcf,'Color','w');

errorbar(R1rho_TN3.Offset(R1rho_TN3.SLP==SLPs_TN3(1))./1000,...
    R1rho_TN3.R2eff(R1rho_TN3.SLP==SLPs_TN3(1)),...
    R1rho_TN3.R2effErr(R1rho_TN3.SLP==SLPs_TN3(1)),...
    'o','MarkerEdgeColor',colors_for_R1rho_plot{1},'MarkerFaceColor',colors_for_R1rho_plot{1},...
    'MarkerSize',markersize,'Color',colors_for_R1rho_plot{1},'CapSize',err_capsize,'Linewidth',err_linewidth);

hold on;

plot(sim_TN3.offset(sim_TN3.slp==SLPs_TN3(1))./1000,...
    sim_TN3.r2eff(sim_TN3.slp==SLPs_TN3(1)),...
    '-','Color',colors_for_R1rho_plot{1},'LineWidth',linewidth);

errorbar(R1rho_TN3.Offset(R1rho_TN3.SLP==SLPs_TN3(2))./1000,...
    R1rho_TN3.R2eff(R1rho_TN3.SLP==SLPs_TN3(2)),...
    R1rho_TN3.R2effErr(R1rho_TN3.SLP==SLPs_TN3(2)),...
    'o','MarkerEdgeColor',colors_for_R1rho_plot{2},'MarkerFaceColor',colors_for_R1rho_plot{2},...
    'MarkerSize',markersize,'Color',colors_for_R1rho_plot{2},'CapSize',err_capsize,'Linewidth',err_linewidth);

plot(sim_TN3.offset(sim_TN3.slp==SLPs_TN3(2))./1000,...
    sim_TN3.r2eff(sim_TN3.slp==SLPs_TN3(2)),...
    '-','Color',colors_for_R1rho_plot{2},'LineWidth',linewidth);

errorbar(R1rho_TN3.Offset(R1rho_TN3.SLP==SLPs_TN3(3))./1000,...
    R1rho_TN3.R2eff(R1rho_TN3.SLP==SLPs_TN3(3)),...
    R1rho_TN3.R2effErr(R1rho_TN3.SLP==SLPs_TN3(3)),...
    'o','MarkerEdgeColor',colors_for_R1rho_plot{3},'MarkerFaceColor',colors_for_R1rho_plot{3},...
    'MarkerSize',markersize,'Color',colors_for_R1rho_plot{3},'CapSize',err_capsize,'Linewidth',err_linewidth);

plot(sim_TN3.offset(sim_TN3.slp==SLPs_TN3(3))./1000,...
    sim_TN3.r2eff(sim_TN3.slp==SLPs_TN3(3)),...
    '-','Color',colors_for_R1rho_plot{3},'LineWidth',linewidth);

errorbar(R1rho_TN3.Offset(R1rho_TN3.SLP==SLPs_TN3(4))./1000,...
    R1rho_TN3.R2eff(R1rho_TN3.SLP==SLPs_TN3(4)),...
    R1rho_TN3.R2effErr(R1rho_TN3.SLP==SLPs_TN3(4)),...
    'o','MarkerEdgeColor',colors_for_R1rho_plot{4},'MarkerFaceColor',colors_for_R1rho_plot{4},...
    'MarkerSize',markersize,'Color',colors_for_R1rho_plot{4},'CapSize',err_capsize,'Linewidth',err_linewidth);

plot(sim_TN3.offset(sim_TN3.slp==SLPs_TN3(4))./1000,...
    sim_TN3.r2eff(sim_TN3.slp==SLPs_TN3(4)),...
    '-','Color',colors_for_R1rho_plot{4},'LineWidth',linewidth);
 
errorbar(R1rho_TN3.Offset(R1rho_TN3.SLP==SLPs_TN3(6))./1000,...
    R1rho_TN3.R2eff(R1rho_TN3.SLP==SLPs_TN3(6)),...
    R1rho_TN3.R2effErr(R1rho_TN3.SLP==SLPs_TN3(6)),...
    'o','MarkerEdgeColor',colors_for_R1rho_plot{6},'MarkerFaceColor',colors_for_R1rho_plot{6},...
    'MarkerSize',markersize,'Color',colors_for_R1rho_plot{6},'CapSize',err_capsize,'Linewidth',err_linewidth);

plot(sim_TN3.offset(sim_TN3.slp==SLPs_TN3(6))./1000,...
    sim_TN3.r2eff(sim_TN3.slp==SLPs_TN3(6)),...
    '-','Color',colors_for_R1rho_plot{6},'LineWidth',linewidth);

hold off;

set(gca, 'FontSize', fontsize,'XTick',-10:2:10,...
    'XminorTick', 'off', 'YminorTick', 'off',...
    'LineWidth',axeswidth,'TickLength',[ticksize,ticksize]);

axis([-9.300 9.300 7 34]);

xlabel('Offset (kHz)', 'FontSize', 36);
ylabel('R2+Rex (s-1)', 'FontSize', 36);

propedit;


%% Figure 6B - ddG plots

%% ddG delta melt vs. ddG Anion

pathname_dG_melt = 'dG-delta-melt';
path_dG_melt = sprintf('%s/%s',parent_path,pathname_dG_melt);

filename_dG_melt = sprintf('%s/dG_GT_GC_delta_melt.mat',path_dG_melt);

load(filename_dG_melt);

[~,ind_dG_melt]  = ismember(dG_melt_GT_seq_T,sorted_seqs_T_19F);

ddG_melt_GT_25C      = dG_melt_GT_25C - dG_melt_GT_25C(strcmp(dG_melt_GT_seq_T,ref_seq_dG));
ddG_melt_GT_25C_errs = sqrt(dG_melt_GT_25C_errs.^2 + (dG_melt_GT_25C_errs(strcmp(dG_melt_GT_seq_T,ref_seq_dG))).^2);

ddG_melt_GC_25C      = dG_melt_GC_25C - dG_melt_GC_25C(strcmp(dG_melt_GC_seq_T,ref_seq_dG));
ddG_melt_GC_25C_errs = sqrt(dG_melt_GC_25C_errs.^2 + (dG_melt_GC_25C_errs(strcmp(dG_melt_GC_seq_T,ref_seq_dG))).^2);

dddG_melt_GT_GC_25C      = ddG_melt_GT_GC_25C - ddG_melt_GT_GC_25C(strcmp(dG_melt_GT_seq_T,ref_seq_dG));
dddG_melt_GT_GC_25C_errs = sqrt(ddG_melt_GT_GC_25C_errs.^2 + (ddG_melt_GT_GC_25C_errs(strcmp(dG_melt_GT_seq_T,ref_seq_dG))).^2);


%% Figure 6B left - ddG Anion vs. dddG melt (GT-GC)

x_data      = dddG_melt_GT_GC_25C;
x_data_errs = dddG_melt_GT_GC_25C_errs;
y_data      = ddG_Anion_19F_7p4_1C(ind_dG_melt);
y_data_errs = ddG_Anion_19F_7p4_1C_errs(ind_dG_melt);

[ddG_Anion_19F_vs_dddG_melt_r,ddG_Anion_19F_vs_dddG_melt_Rsq,ddG_Anion_19F_vs_dddG_melt_Rsq_det,...
    ddG_Anion_19F_vs_dddG_melt_rms,ddG_Anion_19F_vs_dddG_melt_rms_per,...
    ddG_Anion_19F_vs_dddG_melt_bestfit,x2_dddG,inBetween_dddG,...
    ddG_Anion_19F_vs_dddG_melt_slope,ddG_Anion_19F_vs_dddG_melt_inter,...
    ddG_Anion_19F_vs_dddG_melt_slope_err,ddG_Anion_19F_vs_dddG_melt_inter_err] = corr_plot(x_data,y_data,x_line,conf_int_per);

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig 6B left",NumberTitle="off");
set(gcf,'Color','w');

errorbar(x_data,y_data,y_data_errs,y_data_errs,...
    x_data_errs,x_data_errs,...
    'o','MarkerEdgeColor',color_GT,'MarkerFaceColor',color_GT,'MarkerSize',markersize_ddG,...
    'Color',color_GT,'CapSize',err_capsize,'Linewidth',axeswidth);

hold on;

plot(x_line,ddG_Anion_19F_vs_dddG_melt_bestfit,'k--','LineWidth',2);

hold off;

axis([-2.3 1.3 -0.5 2.9]); % for ddG, ref = CTC

set(gca, 'FontSize', fontsize, 'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

text(-2.2,2.6,sprintf('RMSD = %0.1f kcal/mol',ddG_Anion_19F_vs_dddG_melt_rms),'FontSize',fontsize,'FontWeight','bold');
text(-2.2,2.4,sprintf('r = %0.2f',ddG_Anion_19F_vs_dddG_melt_r),'FontSize',fontsize,'FontWeight','bold');

xlabel('DDDG melt (GT-GC)', 'FontSize', 36);
ylabel('DDG Anion', 'FontSize', 36);

propedit;


%% Figure 6B middle - ddG Anion vs. ddG melt (GT)

x_data      = ddG_melt_GT_25C;
x_data_errs = ddG_melt_GT_25C_errs;
y_data      = ddG_Anion_19F_7p4_1C(ind_dG_melt);
y_data_errs = ddG_Anion_19F_7p4_1C_errs(ind_dG_melt);

[ddG_Anion_19F_vs_ddG_melt_GT_r,ddG_Anion_19F_vs_ddG_melt_GT_Rsq,ddG_Anion_19F_vs_ddG_melt_GT_Rsq_det,...
    ddG_Anion_19F_vs_ddG_melt_GT_rms,ddG_Anion_19F_vs_ddG_melt_GT_rms_per,...
    ddG_Anion_19F_vs_ddG_melt_GT_bestfit,x2_ddG,inBetween_ddG,...
    ddG_Anion_19F_vs_ddG_melt_GT_slope,ddG_Anion_19F_vs_ddG_melt_GT_inter,...
    ddG_Anion_19F_vs_ddG_melt_GT_slope_err,ddG_Anion_19F_vs_ddG_melt_GT_inter_err] = corr_plot(x_data,y_data,x_line,conf_int_per);

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig 6B middle",NumberTitle="off");
set(gcf,'Color','w');

errorbar(x_data,y_data,y_data_errs,y_data_errs,...
    x_data_errs,x_data_errs,...
    'o','MarkerEdgeColor',color_GT,'MarkerFaceColor',color_GT,'MarkerSize',markersize_ddG,...
    'Color',color_GT,'CapSize',err_capsize,'Linewidth',axeswidth);

hold on;

plot(x_line,ddG_Anion_19F_vs_ddG_melt_GT_bestfit,'k--','LineWidth',2);

hold off;

axis([-5.4 1.1 -0.5 2.9]); % for ddG, ref = CTC

set(gca, 'FontSize', fontsize, 'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

text(-5.2,2.6,sprintf('RMSD = %0.1f kcal/mol',ddG_Anion_19F_vs_ddG_melt_GT_rms),'FontSize',fontsize,'FontWeight','bold');
text(-5.2,2.4,sprintf('r = %0.2f',ddG_Anion_19F_vs_ddG_melt_GT_r),'FontSize',fontsize,'FontWeight','bold');

xlabel('DDDG melt (GT)', 'FontSize', 36);
ylabel('DDG Anion', 'FontSize', 36);

propedit;


%% Figure 6B right - ddG_Anion (19F) vs. ddG_Tautomer (R1rho)

seq_DNA_6p9_25C_T = Taut_DNA_6p9_25C_pops_for_dG.seq_DNA_6p9_25C_T;
[~,ind_dG_Taut] = ismember(seq_DNA_6p9_25C_T,sorted_seqs_T_19F);

ddG_Taut_25C_6p9           = dG_Taut_6p9_25C - dG_Taut_6p9_25C(ismember(seq_DNA_6p9_25C_T,ref_seq_dG));
ddG_Taut_25C_6p9_errs      = sqrt(dG_Taut_6p9_25C_errs.^2 + (0.1).^2);

x_data      = ddG_Taut_25C_6p9;
x_data_errs = ddG_Taut_25C_6p9_errs;
y_data      = ddG_Anion_19F_7p4_1C(ind_dG_Taut);
y_data_errs = ddG_Anion_19F_7p4_1C_errs(ind_dG_Taut);

[ddG_Anion_19F_vs_ddG_Taut_r,ddG_Anion_19F_vs_ddG_Taut_Rsq,ddG_Anion_19F_vs_ddG_Taut_Rsq_det,...
    ddG_Anion_19F_vs_ddG_Taut_rms,ddG_Anion_19F_vs_ddG_Taut_rms_per,...
    ddG_Anion_19F_vs_ddG_Taut_bestfit,x2_ddG_Taut,inBetween_ddG_Taut,...
    ddG_Anion_19F_vs_ddG_Taut_slope,ddG_Anion_19F_vs_ddG_Taut_inter,...
    ddG_Anion_19F_vs_ddG_Taut_slope_err,ddG_Anion_19F_vs_ddG_Taut_inter_err] = corr_plot(x_data,y_data,x_line,conf_int_per);

figure('units','normalized','outerposition',[0 0 1 1],Name="Fig 6B right",NumberTitle="off");
set(gcf,'Color','w');

errorbar(x_data,y_data,y_data_errs,y_data_errs,...
    x_data_errs,x_data_errs,...
    'o','MarkerEdgeColor',color_GT,'MarkerFaceColor',color_GT,'MarkerSize',markersize,...
    'Color',color_GT,'CapSize',err_capsize,'Linewidth',axeswidth);

hold on;

plot(x_line,ddG_Anion_19F_vs_ddG_Taut_bestfit,'k--','LineWidth',err_linewidth);

hold off;

axis([-0.2 1.1 -0.5 2.9]);

set(gca, 'FontSize', fontsize, 'XminorTick', 'on', 'YminorTick', 'on','LineWidth',axeswidth,...
    'TickLength',[ticksize,ticksize]);

text(-0.15,2.6,sprintf('RMSD = %0.1f kcal/mol',ddG_Anion_19F_vs_ddG_Taut_rms),'FontSize',fontsize,'FontWeight','bold');
text(-0.15,2.4,sprintf('r = %0.2f',ddG_Anion_19F_vs_ddG_Taut_r),'FontSize',fontsize,'FontWeight','bold');

xlabel('DDG Tautomer', 'FontSize', 36);
ylabel('DDG Anion', 'FontSize', 36);

propedit;

