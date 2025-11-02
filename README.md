#########################################################################
# DNA-conformational-fingerprinting
#  Processed Data and Scripts to analyze and plot data for manuscript
#  entitled:
#  Assessing the contribution of rare DNA states to cancer
#  mutational signatures using sequence-specific conformational fingerprinting
#  Or Szekely
#  2025-11-02
#########################################################################

This repository contains data folders and scripts for analysis and plotting.

1. 1-Matlab-Functions-add-to-path
A folder containing Matlab in-house functions necessary for running scripts and analyzing data. Add this folder to the Matlab path using the command: addpath('pathname'); 

2. NMR Spectra
2.1. Process 1D NMR spectra using NMRPipe (Version 10.6 Revision 2020.007.13.34) and convert spectra to txt format (.txt).
2.2. Spectra .txt files are saved in NMR-Spectra/1D-1H-imino-spectra and NMR-Spectra/1D-19F-spectra
2.3. Use figures_all.m to load the .txt files and plot 1D spectra
2.4. Process 2D NMR spectra (NOESY, HSQC) using NMRPipe (Version 10.6 Revision 2020.007.13.34) and convert spectrum to SPARKY format (.ucsf). 
2.5. 2D spectra are available upon request. Open ucsf files using SPARKY software. Use SPARKY to create figure of overlayed spectra.
2.6. Peak the peaks and save the peak lists in the folder CSPs (CSPs/GTA-peak-lists; CSPs/ATT-peak-lists; CSPs/CTG-peak-lists; CSPs/TTG-peak-lists).
2.7. Use CSPs/csps_GTA_vs_WB; CSPs/csps_ATT_vs_WB; CSPs/csps_CTG_vs_WB; CSPs/csps_TTG_vs_WB to read the peak lists and save .mat files with the chemical shift differences (CSPs).

3. Populations-and-Chemical-Shifts
3.1. Extract chemical shifts from the 1D spectra using NMRPipe (Version 10.6 Revision 2020.007.13.34). Chemical shifts are saved in Populations-and-Chemical-Shifts.19F-and-1H-Chemical-Shifts.xlsx
3.2. Use figures_all.m to load the 1D 19F spectra .txt files, fit the peak shapes to Lorentzian functions, and extract the GT anion populations. Populations are saved in Populations-and-Chemical-Shifts/Anion_populations_19F_1C_all_pHs.xlsx
3.3. Also included is a list of unique permutations of parent and destination sequences for additivity calculations taken from Manghrani et al (DOI: 10.1021/acs.biochem.4c00820)
3.4. Chemical shift additivity (as well as pKa additivity) can be calculated using figures_all.m.

4. Free energy difference (delta G) calculations
4.1. Previous exchange parameters measured using R1rho Relaxation Dispersion (RD) by Kimsey et al (DOI: 10.1038/nature25487), are saved in Previous-R1rho-Params/2023-GT-anion-seq-dep-incl-Genol-Tenol.xlsx.
4.2. Use Previous-R1rho-Params/save_R1rho_params.m to read the RD parameters and save both anion and tautomer data, including delta G, and save in: Previous-R1rho-Params/Prev-Anion-Populations; Previous-R1rho-Params/Prev-Tautomer-Populations; Previous-R1rho-Params/Prev-dG-Anion-and-Tautomer.
4.3. Melting energies are saved in dG-delta-melt/ddGmelt-GC-vs-GT.xlsx
4.4. Use dG-delta-melt/ddG_delta_melt.m to read the melting deltaG values and save them as a Matlab .mat file.
4.5. Use figures_all.m to load the dG values from above, compute the difference in free energy with sequence (ddG), and plot against ddG(anion).

5. NMR R1rho data
5.1. Raw R1rho data was processed using Disprun (https://github.com/alhashimilab/Disprun.git) based on NMRPipe (Version 10.6 Revision 2020.007.13.34) to extract peak intensities and fit the relaxation curves using mono-exponential decays.
5.2. R1rho values are then fit to Bloch-McConnell equations using BMNS (https://github.com/alhashimilab/BMNS.git). The data from two probes are globally fit, to a 3-state exchange model with triangular topology. Resulting fits are saved as csv files in R1rho-GTA/3state-triangle-bestfit.
5.3. The fitting parameters are used to create an RD curve using BMNS simulation (https://github.com/alhashimilab/BMNS.git). The simulations are saved in R1rho-GTA/R1rho-Figure.
5.4. Use the script figures_all.m to plot the R1rho profiles (R2eff vs. omega_obs).

6. Calculating Jenssen-Shannon-Divergence between conformational fingerprints and cancer mutational signatures - can be done either using an in-house Matlab script or an in-house python script
6.1. Matlab:
6.1.1. Located at Mutational-Signatures-JSD/JSD-Matlab/mut_sig_vs_conf_fingerprint.m
6.1.2. In the script choose which dataset to compare to mutational signatures (GT_Anion / AT_HG / GC_HG / AT_BaseOpening). The script will read the input file containing the conformational fingerprint (sequence contexts and rare state populations), and compute the Jenssen-Shannon-Divergence (JSD) with the cancer mutational signatures file (.txt file downloaded from the COSMIC database).
6.1.3. Output .csv files contain the full table of JSDs per substitution (transition mutation C>A, C>G, C>T, T>A, T>C, or T>G); as well as the best similarities between the SBS signature and the conformational fingerprint.
6.1.4. Another output file is a JSD table computed against a random distribution, for calculation of p-values and FDRs.

6.2. Python:
6.2.1. Located at Mutational-Signatures-JSD/JSD-python/mut_sig_vs_conf_fingerprint.py
Requirements: Python 3.8 or higher. See inside for package version requirements. 
6.2.2. In the script, under main() choose which dataset to compare to mutational signatures (conf_filename = "GT_Anion.csv" / "AT_HG.csv" / "GC_HG.csv" / "AT_BaseOpening.csv"). See inside for file requirements. The script will read the input file containing the conformational fingerprint (sequence contexts and rare state populations), and compute the Jenssen-Shannon-Divergence (JSD) with the cancer mutational signatures file (.txt file downloaded from the COSMIC database).
6.2.3. Output .csv files contain the full table of JSDs (JSD_summary_unfiltered.csv); as well as the best similarities between the SBS signature and the conformational fingerprint (JSD_summary_filtered.csv).
6.2.4. p-values and FDRs are calculated using against a random distribution.
