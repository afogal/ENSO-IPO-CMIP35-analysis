%% Licensing (BSD 3 Clause)
% Copyright 2022 Alexander Fogal
%
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
%
% 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% Code not written by me may carry different licensing
% CMIP data does not fall under this license

%% init
close all;
reprocess_data=0;
redo_eofs=0;
redo_tpi=1;

%% main sequence

% define the box we want
latlim = [-50 50];
lonlim = [120 300];

% either freshly process the data or load it
if reprocess_data
    % load CMIP3 data
    CSIRO35_sst = load_CSIRO35(latlim,lonlim);
    HadGem1_sst = load_HadGem1(latlim,lonlim);
    CCSM3_sst = load_CCSM3(latlim,lonlim);

    % load CMIP5 data
    CSIRO36_sst = load_CSIRO36(latlim,lonlim);
    HadGem2_sst = load_HadGem2(latlim,lonlim);
    CCSM4_sst = load_CCSM4(latlim,lonlim);
    
    % load observational and reanalysis datasets
    hadsst_sst = load_hadsst(latlim,lonlim);
    noaa_sst = load_noaa(latlim,lonlim);
    
    % save amalgam
    save("CMIPdata/allmodels_sst.mat", "CSIRO35_sst","CSIRO36_sst", "HadGem1_sst", "HadGem2_sst","CCSM3_sst","CCSM4_sst","hadsst_sst","noaa_sst");
   
else
    load("CMIPdata/allmodels_sst.mat", "CSIRO35_sst","CSIRO36_sst", "HadGem1_sst", "HadGem2_sst","CCSM3_sst","CCSM4_sst","hadsst_sst","noaa_sst");
end    

% either save or load eof analysis
if redo_eofs
    % comptue all EOFs and save
    [CCSM3_eofs,CCSM3_eigvs,CCSM3_pcs] = get_eofs(CCSM3_sst);
    [CCSM4_eofs,CCSM4_eigvs,CCSM4_pcs] = get_eofs(CCSM4_sst);

    [HadGem1_eofs,HadGem1_eigvs,HadGem1_pcs] = get_eofs(HadGem1_sst);
    [HadGem2_eofs,HadGem2_eigvs,HadGem2_pcs] = get_eofs(HadGem2_sst);

    [CSIRO35_eofs,CSIRO35_eigvs,CSIRO35_pcs] = get_eofs(CSIRO35_sst);
    [CSIRO36_eofs,CSIRO36_eigvs,CSIRO36_pcs] = get_eofs(CSIRO36_sst);

    [hadsst_eofs,hadsst_eigvs,hadsst_pcs] = get_eofs(hadsst_sst);
    [noaa_eofs,noaa_eigvs,noaa_pcs] = get_eofs(noaa_sst);
    
    save("CMIPdata/allmodels_eofs.mat","CCSM3_eofs","CCSM3_eigvs","CCSM3_pcs","CCSM4_eofs","CCSM4_eigvs","CCSM4_pcs","HadGem1_eofs","HadGem1_eigvs","HadGem1_pcs","HadGem2_eofs","HadGem2_eigvs","HadGem2_pcs","CSIRO35_eofs","CSIRO35_eigvs","CSIRO35_pcs","CSIRO36_eofs","CSIRO36_eigvs","CSIRO36_pcs","hadsst_eofs","hadsst_eigvs","hadsst_pcs","noaa_eofs","noaa_eigvs","noaa_pcs");
else
    load("CMIPdata/allmodels_eofs.mat","CCSM3_eofs","CCSM3_eigvs","CCSM3_pcs","CCSM4_eofs","CCSM4_eigvs","CCSM4_pcs","HadGem1_eofs","HadGem1_eigvs","HadGem1_pcs","HadGem2_eofs","HadGem2_eigvs","HadGem2_pcs","CSIRO35_eofs","CSIRO35_eigvs","CSIRO35_pcs","CSIRO36_eofs","CSIRO36_eigvs","CSIRO36_pcs","hadsst_eofs","hadsst_eigvs","hadsst_pcs","noaa_eofs","noaa_eigvs","noaa_pcs");
end

if redo_tpi
   
    [CCSM3_TPI, CCSM3_TPI_filt] = get_tpi(CCSM3_sst,latlim,lonlim);
    [CCSM4_TPI, CCSM4_TPI_filt] = get_tpi(CCSM4_sst,latlim,lonlim);
    
    [CSIRO35_TPI, CSIRO35_TPI_filt] = get_tpi(CSIRO35_sst,latlim,lonlim);
    [CSIRO36_TPI, CSIRO36_TPI_filt] = get_tpi(CSIRO36_sst,latlim,lonlim);
    
    [HadGem1_TPI, HadGem1_TPI_filt] = get_tpi(HadGem1_sst,latlim,lonlim);
    [HadGem2_TPI, HadGem2_TPI_filt] = get_tpi(HadGem2_sst,latlim,lonlim);
    
    [hadsst_TPI, hadsst_TPI_filt] = get_tpi(hadsst_sst,latlim,lonlim);
    [noaa_TPI, noaa_TPI_filt] = get_tpi(noaa_sst,latlim,lonlim);
    
    save("CMIPdata/allmodels_tpi.mat","CCSM3_TPI", "CCSM3_TPI_filt","CCSM4_TPI", "CCSM4_TPI_filt","CSIRO35_TPI", "CSIRO35_TPI_filt","CSIRO36_TPI", "CSIRO36_TPI_filt","HadGem1_TPI", "HadGem1_TPI_filt","HadGem2_TPI", "HadGem2_TPI_filt","hadsst_TPI", "hadsst_TPI_filt","noaa_TPI", "noaa_TPI_filt");
else
    load("CMIPdata/allmodels_tpi.mat","CCSM3_TPI", "CCSM3_TPI_filt","CCSM4_TPI", "CCSM4_TPI_filt","CSIRO35_TPI", "CSIRO35_TPI_filt","CSIRO36_TPI", "CSIRO36_TPI_filt","HadGem1_TPI", "HadGem1_TPI_filt","HadGem2_TPI", "HadGem2_TPI_filt","hadsst_TPI", "hadsst_TPI_filt","noaa_TPI", "noaa_TPI_filt");
end

% Figure 1
%plot_eigspec(CCSM3_sst,"CCSM3",50);

% EOF Figures
% plot_eofs(CCSM3_eofs,CCSM3_eigvs,4,latlim,lonlim,1.5,"CCSM3");
% plot_eofs(CCSM4_eofs,CCSM4_eigvs,4,latlim,lonlim,1.5,"CCSM4");
% 
% plot_eofs(CSIRO35_eofs,CSIRO35_eigvs,4,latlim,lonlim,1.5,"CSIRO mk3.5");
% plot_eofs(CSIRO36_eofs,CSIRO36_eigvs,4,latlim,lonlim,1.5,"CSIRO mk3.6");
% 
% plot_eofs(HadGem1_eofs,HadGem1_eigvs,4,latlim,lonlim,1.5,"HadGEM1");
% plot_eofs(HadGem2_eofs,HadGem2_eigvs,4,latlim,lonlim,1.5,"HadGEM2-ES");
% 
% plot_eofs(hadsst_eofs,hadsst_eigvs,4,latlim,lonlim,1.5,"HadSST3");
% plot_eofs(noaa_eofs,noaa_eigvs,4,latlim,lonlim,1.5,"NOAA NCEP OI v2");

% FFT figures
% plot_fftpcs(CCSM3_pcs,4,"CCSM3");
% plot_fftpcs(CCSM4_pcs,4,"CCSM4");
% 
% plot_fftpcs(CSIRO35_pcs,4,"CSIRO mk3.5");
% plot_fftpcs(CSIRO36_pcs,4,"CSIRO mk3.6");
% 
% plot_fftpcs(HadGem1_pcs,4,"HadGEM1");
% plot_fftpcs(HadGem2_pcs,4,"HadGEM2-ES");
% 
%  plot_fftpcs(hadsst_pcs,4,"HadSST3");
% plot_fftpcs(noaa_pcs,4,"NOAA NCEP OI v2");

% TPI figures
plot_tpi(CCSM3_TPI, CCSM3_TPI_filt,"CCSM3",1870*12);
plot_tpi(CCSM4_TPI, CCSM4_TPI_filt,"CCSM4",1850*12);

plot_tpi(CSIRO35_TPI, CSIRO35_TPI_filt,"CSIRO mk3.5",1871*12);
plot_tpi(CSIRO36_TPI, CSIRO36_TPI_filt,"CSIRO mk3.6",1850*12);

plot_tpi(HadGem1_TPI, HadGem1_TPI_filt,"HadGEM1",1860*12);
plot_tpi(HadGem2_TPI, HadGem2_TPI_filt,"HadGEM2-ES",1850*12-1);

plot_tpi(hadsst_TPI,hadsst_TPI_filt,"HadSST3",1850*12);
plot_tpi(noaa_TPI,noaa_TPI_filt,"NOAA NCEP OI v2",1981*12+8);


% end main sequence

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HERE BE FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% functions for loading CMIP3 data
% takes ensemble means if multiple members available
% interpolates to a common 1.5x1.5 grid
% masks out land points

function [CSIRO35_sst] = load_CSIRO35(latlim,lonlim)
    %% CSIRO mk3.5 raw data & ensemble mean
    ssts1 = ncread("CMIPdata/tas_A1_csiro_mk3_5_00_120-300E_-50-50N_1871_2020_anom.nc","tas");
    ssts2 = ncread("CMIPdata/tas_A1_csiro_mk3_5_01_120-300E_-50-50N_1871_2020_anom.nc","tas");
    ssts3 = ncread("CMIPdata/tas_A1_csiro_mk3_5_02_120-300E_-50-50N_1871_2020_anom.nc","tas");
    
    % get size, make new array, take ensemble mean
    original_size = size(ssts1);
    CSIRO_sst = zeros(original_size(1),original_size(2),original_size(3),3);
    CSIRO_sst(:,:,:,1) = ssts1;
    CSIRO_sst(:,:,:,2) = ssts2;
    CSIRO_sst(:,:,:,3) = ssts3;
    CSIRO_sst = squeeze(mean(CSIRO_sst,4));
    
    %% interp to 1.5x1.5
    % first need a geo reference object
    R = georefcells(latlim,lonlim,[original_size(2),original_size(1)],'ColumnsStartFrom','south', 'RowsStartFrom','west');

    % define new lat/lon arrays and grids
    interp_lats = latlim(1):1.5:latlim(2);
    interp_lons = lonlim(1):1.5:lonlim(2);
    [interp_lats_grid,interp_lons_grid] = meshgrid(interp_lats,interp_lons);

    % define a place to put the interp results
    interp_size = size(interp_lats_grid);
    CSIRO35_sst = zeros(interp_size(1),interp_size(2),original_size(3));

    % do for all time (can do this in a better way?)
    for ti=1:original_size(3)
        curr = squeeze(CSIRO_sst(:,:,ti));
        % runs fasters with these transposes
        CSIRO35_sst(:,:,ti) = geointerp(curr',R,interp_lats_grid',interp_lons_grid')';
    end
    
    %% landmask
    mask = landmask(interp_lats,interp_lons,8,0);
    
    mask = repmat(mask, [1,1,original_size(3)]);
    CSIRO35_sst(~mask) = nan;
    
    %save("CMIPdata/CSIRO35_sst.mat", "CSIRO35_sst");

end

function [HadGem1_sst] = load_HadGem1(latlim,lonlim)
    %% HadGem1 raw data & ensemble mean
    ssts1 = ncread("CMIPdata/tas_A1_ukmo_hadgem1_00_120-300E_-50-50N_1860_2020_anom.nc","tas");
    ssts2 = ncread("CMIPdata/tas_A1_ukmo_hadgem1_01_120-300E_-50-50N_1860_2020_anom.nc","tas");
    
    % get size, make new array, take ensemble mean
    original_size = size(ssts1);
    HadGem_sst = zeros(original_size(1),original_size(2),original_size(3),2);
    HadGem_sst(:,:,:,1) = ssts1;
    HadGem_sst(:,:,:,2) = ssts2;
    HadGem_sst = squeeze(mean(HadGem_sst,4));
    
    %% interp to 1.5x1.5
    % first need a geo reference object
    R = georefcells(latlim,lonlim,[original_size(2),original_size(1)],'ColumnsStartFrom','south', 'RowsStartFrom','west');

    % define new lat/lon arrays and grids
    interp_lats = latlim(1):1.5:latlim(2);
    interp_lons = lonlim(1):1.5:lonlim(2);
    [interp_lats_grid,interp_lons_grid] = meshgrid(interp_lats,interp_lons);

    % define a place to put the interp results
    interp_size = size(interp_lats_grid);
    HadGem1_sst = zeros(interp_size(1),interp_size(2),original_size(3));

    % do for all time (can do this in a better way?)
    for ti=1:original_size(3)
        curr = squeeze(HadGem_sst(:,:,ti));
        % runs fasters with these transposes
        HadGem1_sst(:,:,ti) = geointerp(curr',R,interp_lats_grid',interp_lons_grid')';
    end
    
    %% landmask
    mask = landmask(interp_lats,interp_lons,8,0);
    
    mask = repmat(mask, [1,1,original_size(3)]);
    HadGem1_sst(~mask) = nan;
    
    %save("CMIPdata/HadGem1_sst.mat", "HadGem1_sst");

end

function [CCSM3_sst] = load_CCSM3(latlim,lonlim)
    %% CSIRO mk3.5 raw data & ensemble mean
    ssts1 = ncread("CMIPdata/tas_A1.20C3M_00.CCSM.atmm.1870-01_cat_1999-12_120-300E_-50-50N_1870_2020_anom.nc","tas");
    ssts2 = ncread("CMIPdata/tas_A1.20C3M_01.CCSM.atmm.1870-01_cat_1999-12_120-300E_-50-50N_1870_2020_anom.nc","tas");
    ssts3 = ncread("CMIPdata/tas_A1.20C3M_02.CCSM.atmm.1870-01_cat_1999-12_120-300E_-50-50N_1870_2020_anom.nc","tas");
    ssts4 = ncread("CMIPdata/tas_A1.20C3M_03.CCSM.atmm.1870-01_cat_1999-12_120-300E_-50-50N_1870_2020_anom.nc","tas");
    ssts5 = ncread("CMIPdata/tas_A1.20C3M_04.CCSM.atmm.1870-01_cat_1999-12_120-300E_-50-50N_1870_2020_anom.nc","tas");
    ssts6 = ncread("CMIPdata/tas_A1.20C3M_05.CCSM.atmm.1870-01_cat_1999-12_120-300E_-50-50N_1870_2020_anom.nc","tas");
    
    % get size, make new array, take ensemble mean
    original_size = size(ssts1);
    CCSM_sst = zeros(original_size(1),original_size(2),original_size(3),6);
    CCSM_sst(:,:,:,1) = ssts1;
    CCSM_sst(:,:,:,2) = ssts2;
    CCSM_sst(:,:,:,3) = ssts3;
    CCSM_sst(:,:,:,4) = ssts4;
    CCSM_sst(:,:,:,5) = ssts5;
    CCSM_sst(:,:,:,6) = ssts6;
    CCSM_sst = squeeze(mean(CCSM_sst,4));
    
    %% interp to 1.5x1.5
    % first need a geo reference object
    R = georefcells(latlim,lonlim,[original_size(2),original_size(1)],'ColumnsStartFrom','south', 'RowsStartFrom','west');

    % define new lat/lon arrays and grids
    interp_lats = latlim(1):1.5:latlim(2);
    interp_lons = lonlim(1):1.5:lonlim(2);
    [interp_lats_grid,interp_lons_grid] = meshgrid(interp_lats,interp_lons);

    % define a place to put the interp results
    interp_size = size(interp_lats_grid);
    CCSM3_sst = zeros(interp_size(1),interp_size(2),original_size(3));

    % do for all time (can do this in a better way?)
    for ti=1:original_size(3)
        curr = squeeze(CCSM_sst(:,:,ti));
        % runs fasters with these transposes
        CCSM3_sst(:,:,ti) = geointerp(curr',R,interp_lats_grid',interp_lons_grid')';
    end
    
    %% landmask
    mask = landmask(interp_lats,interp_lons,8,0);
    
    mask = repmat(mask, [1,1,original_size(3)]);
    CCSM3_sst(~mask) = nan;
    
    %save("CMIPdata/CCSM3_sst.mat", "CCSM3_sst");

end

%% functions for loading CMIP5 data
% same as above really

function [CSIRO36_sst] = load_CSIRO36(latlim,lonlim)
    %% CSIRO mk3.6 raw data & ensemble mean
    ssts1 = ncread("CMIPdata/tas_Amon_CSIRO-Mk3-6-0_historical_000_120-300E_-50-50N_1850_2020_anom.nc","tas");
    ssts2 = ncread("CMIPdata/tas_Amon_CSIRO-Mk3-6-0_historical_001_120-300E_-50-50N_1850_2020_anom.nc","tas");
    ssts3 = ncread("CMIPdata/tas_Amon_CSIRO-Mk3-6-0_historical_002_120-300E_-50-50N_1850_2020_anom.nc","tas");
    ssts4 = ncread("CMIPdata/tas_Amon_CSIRO-Mk3-6-0_historical_003_120-300E_-50-50N_1850_2020_anom.nc","tas");
    ssts5 = ncread("CMIPdata/tas_Amon_CSIRO-Mk3-6-0_historical_004_120-300E_-50-50N_1850_2020_anom.nc","tas");
    ssts6 = ncread("CMIPdata/tas_Amon_CSIRO-Mk3-6-0_historical_005_120-300E_-50-50N_1850_2020_anom.nc","tas");
    ssts7 = ncread("CMIPdata/tas_Amon_CSIRO-Mk3-6-0_historical_006_120-300E_-50-50N_1850_2020_anom.nc","tas");
    ssts8 = ncread("CMIPdata/tas_Amon_CSIRO-Mk3-6-0_historical_007_120-300E_-50-50N_1850_2020_anom.nc","tas");
    ssts9 = ncread("CMIPdata/tas_Amon_CSIRO-Mk3-6-0_historical_008_120-300E_-50-50N_1850_2020_anom.nc","tas");
    ssts10 = ncread("CMIPdata/tas_Amon_CSIRO-Mk3-6-0_historical_009_120-300E_-50-50N_1850_2020_anom.nc","tas");
    
    % get size, make new array, take ensemble mean
    original_size = size(ssts1);
    CSIRO_sst = zeros(original_size(1),original_size(2),original_size(3),10);
    CSIRO_sst(:,:,:,1) = ssts1;
    CSIRO_sst(:,:,:,2) = ssts2;
    CSIRO_sst(:,:,:,3) = ssts3;
    CSIRO_sst(:,:,:,4) = ssts4;
    CSIRO_sst(:,:,:,5) = ssts5;
    CSIRO_sst(:,:,:,6) = ssts6;
    CSIRO_sst(:,:,:,7) = ssts7;
    CSIRO_sst(:,:,:,8) = ssts8;
    CSIRO_sst(:,:,:,9) = ssts9;
    CSIRO_sst(:,:,:,10) = ssts10;
    CSIRO_sst = squeeze(mean(CSIRO_sst,4));
    
    %% interp to 1.5x1.5
    % first need a geo reference object
    R = georefcells(latlim,lonlim,[original_size(2),original_size(1)],'ColumnsStartFrom','south', 'RowsStartFrom','west');

    % define new lat/lon arrays and grids
    interp_lats = latlim(1):1.5:latlim(2);
    interp_lons = lonlim(1):1.5:lonlim(2);
    [interp_lats_grid,interp_lons_grid] = meshgrid(interp_lats,interp_lons);

    % define a place to put the interp results
    interp_size = size(interp_lats_grid);
    CSIRO36_sst = zeros(interp_size(1),interp_size(2),original_size(3));

    % do for all time (can do this in a better way?)
    for ti=1:original_size(3)
        curr = squeeze(CSIRO_sst(:,:,ti));
        % runs fasters with these transposes
        CSIRO36_sst(:,:,ti) = geointerp(curr',R,interp_lats_grid',interp_lons_grid')';
    end
    
    %% landmask
    mask = landmask(interp_lats,interp_lons,8,0);
    
    mask = repmat(mask, [1,1,original_size(3)]);
    CSIRO36_sst(~mask) = nan;
    
    %save("CMIPdata/CSIRO36_sst.mat", "CSIRO36_sst");

end

function [HadGem2_sst] = load_HadGem2(latlim,lonlim)
    %% HadGem2 raw data & ensemble mean
    ssts1 = ncread("CMIPdata/tas_Amon_HadGEM2-ES_historical_000_120-300E_-50-50N_1859_2020_anom.nc","tas");
    ssts2 = ncread("CMIPdata/tas_Amon_HadGEM2-ES_historical_001_120-300E_-50-50N_1859_2020_anom.nc","tas");
    ssts3 = ncread("CMIPdata/tas_Amon_HadGEM2-ES_historical_002_120-300E_-50-50N_1859_2020_anom.nc","tas");
    ssts4 = ncread("CMIPdata/tas_Amon_HadGEM2-ES_historical_003_120-300E_-50-50N_1859_2020_anom.nc","tas");
    
    % get size, make new array, take ensemble mean
    original_size = size(ssts1);
    HadGem_sst = zeros(original_size(1),original_size(2),original_size(3),4);
    HadGem_sst(:,:,:,1) = ssts1;
    HadGem_sst(:,:,:,2) = ssts2(:,:,1:(end-1)); % for some reason these two ensemble members have an extra month at the end
    HadGem_sst(:,:,:,3) = ssts3(:,:,1:(end-1));
    HadGem_sst(:,:,:,4) = ssts4;
    HadGem_sst = squeeze(mean(HadGem_sst,4));
    
    %% interp to 1.5x1.5
    % first need a geo reference object
    R = georefcells(latlim,lonlim,[original_size(2),original_size(1)],'ColumnsStartFrom','south', 'RowsStartFrom','west');

    % define new lat/lon arrays and grids
    interp_lats = latlim(1):1.5:latlim(2);
    interp_lons = lonlim(1):1.5:lonlim(2);
    [interp_lats_grid,interp_lons_grid] = meshgrid(interp_lats,interp_lons);

    % define a place to put the interp results
    interp_size = size(interp_lats_grid);
    HadGem2_sst = zeros(interp_size(1),interp_size(2),original_size(3));

    % do for all time (can do this in a better way?)
    for ti=1:original_size(3)
        curr = squeeze(HadGem_sst(:,:,ti));
        % runs fasters with these transposes
        HadGem2_sst(:,:,ti) = geointerp(curr',R,interp_lats_grid',interp_lons_grid')';
    end
    
    %% landmask
    mask = landmask(interp_lats,interp_lons,8,0);
    
    mask = repmat(mask, [1,1,original_size(3)]);
    HadGem2_sst(~mask) = nan;
    
    %save("CMIPdata/HadGem2_sst.mat", "HadGem2_sst");

end

function [CCSM4_sst] = load_CCSM4(latlim,lonlim)
    %% CCSM4 raw data & ensemble mean
    ssts1 = ncread("CMIPdata/tas_Amon_CCSM4_historical_000_120-300E_-50-50N_1850_2020_anom.nc","tas");
    ssts2 = ncread("CMIPdata/tas_Amon_CCSM4_historical_001_120-300E_-50-50N_1850_2020_anom.nc","tas");
    ssts3 = ncread("CMIPdata/tas_Amon_CCSM4_historical_002_120-300E_-50-50N_1850_2020_anom.nc","tas");
    ssts4 = ncread("CMIPdata/tas_Amon_CCSM4_historical_003_120-300E_-50-50N_1850_2020_anom.nc","tas");
    ssts5 = ncread("CMIPdata/tas_Amon_CCSM4_historical_004_120-300E_-50-50N_1850_2020_anom.nc","tas");
    ssts6 = ncread("CMIPdata/tas_Amon_CCSM4_historical_005_120-300E_-50-50N_1850_2020_anom.nc","tas");
    
    
    % get size, make new array, take ensemble mean
    original_size = size(ssts1);
    CCSM_sst = zeros(original_size(1),original_size(2),original_size(3),6);
    CCSM_sst(:,:,:,1) = ssts1;
    CCSM_sst(:,:,:,2) = ssts2;
    CCSM_sst(:,:,:,3) = ssts3;
    CCSM_sst(:,:,:,4) = ssts4;
    CCSM_sst(:,:,:,5) = ssts5;
    CCSM_sst(:,:,:,6) = ssts6;
    CCSM_sst = squeeze(mean(CCSM_sst,4));
    
    %% interp to 1.5x1.5
    % first need a geo reference object
    R = georefcells(latlim,lonlim,[original_size(2),original_size(1)],'ColumnsStartFrom','south', 'RowsStartFrom','west');

    % define new lat/lon arrays and grids
    interp_lats = latlim(1):1.5:latlim(2);
    interp_lons = lonlim(1):1.5:lonlim(2);
    [interp_lats_grid,interp_lons_grid] = meshgrid(interp_lats,interp_lons);

    % define a place to put the interp results
    interp_size = size(interp_lats_grid);
    CCSM4_sst = zeros(interp_size(1),interp_size(2),original_size(3));

    % do for all time (can do this in a better way?)
    for ti=1:original_size(3)
        curr = squeeze(CCSM_sst(:,:,ti));
        % runs fasters with these transposes
        CCSM4_sst(:,:,ti) = geointerp(curr',R,interp_lats_grid',interp_lons_grid')';
    end
    
    %% landmask
    mask = landmask(interp_lats,interp_lons,8,0);
    
    mask = repmat(mask, [1,1,original_size(3)]);
    CCSM4_sst(~mask) = nan;
    
    %save("CMIPdata/CCSM4_sst.mat", "CCSM4_sst");

end

%% functions for loading obs/reanalysis data
function [hadsst_sst] = load_hadsst(latlim,lonlim)
    %% HadSST raw data & ensemble mean
    ssts1 = ncread("CMIPdata/HadSST.3.1.1.0.median_120-300E_-50-50N_120-300E_-50-50N.nc","sst");
    
    original_size = size(ssts1);
    
    %% interp to 1.5x1.5
    % first need a geo reference object
    R = georefcells(latlim,lonlim,[original_size(2),original_size(1)],'ColumnsStartFrom','south', 'RowsStartFrom','west');

    % define new lat/lon arrays and grids
    interp_lats = latlim(1):1.5:latlim(2);
    interp_lons = lonlim(1):1.5:lonlim(2);
    [interp_lats_grid,interp_lons_grid] = meshgrid(interp_lats,interp_lons);

    % define a place to put the interp results
    interp_size = size(interp_lats_grid);
    hadsst_sst = zeros(interp_size(1),interp_size(2),original_size(3));

    % do for all time (can do this in a better way?)
    for ti=1:original_size(3)
        curr = squeeze(ssts1(:,:,ti));
        % runs fasters with these transposes
        hadsst_sst(:,:,ti) = geointerp(curr',R,interp_lats_grid',interp_lons_grid')';
    end
    
    %% landmask
    mask = landmask(interp_lats,interp_lons,8,0);
    
    mask = repmat(mask, [1,1,original_size(3)]);
    hadsst_sst(~mask) = nan;
    
    %save("CMIPdata/hadsst_sst.mat", "hadsst_sst");

end

function [noaa_sst] = load_noaa(latlim,lonlim)
    %% NOAA raw data & ensemble mean
    ssts1 = ncread("CMIPdata/sstoi_v2_120-300E_-50-50N_1981_2020_anom.nc","sst");
    
    original_size = size(ssts1);
    
    %% interp to 1.5x1.5
    % first need a geo reference object
    R = georefcells(latlim,lonlim,[original_size(2),original_size(1)],'ColumnsStartFrom','south', 'RowsStartFrom','west');

    % define new lat/lon arrays and grids
    interp_lats = latlim(1):1.5:latlim(2);
    interp_lons = lonlim(1):1.5:lonlim(2);
    [interp_lats_grid,interp_lons_grid] = meshgrid(interp_lats,interp_lons);

    % define a place to put the interp results
    interp_size = size(interp_lats_grid);
    noaa_sst = zeros(interp_size(1),interp_size(2),original_size(3));

    % do for all time (can do this in a better way?)
    for ti=1:original_size(3)
        curr = squeeze(ssts1(:,:,ti));
        % runs fasters with these transposes
        noaa_sst(:,:,ti) = geointerp(curr',R,interp_lats_grid',interp_lons_grid')';
    end
    
    %% landmask
    mask = landmask(interp_lats,interp_lons,8,0);
    
    mask = repmat(mask, [1,1,original_size(3)]);
    noaa_sst(~mask) = nan;

end

%% Visualization functions

% quick n dirty visualization of the data
function plotDatasetSlice(dataset,latlim,lonlim,gridsize,time)
    % pass in time or else use default value
    switch nargin
        case 4
            time=10;
        case 5
            % got everything
    end    

    % make grid of lats
    lats = latlim(1):gridsize:latlim(2);
    lons = lonlim(1):gridsize:lonlim(2);
    [lats_grid,lons_grid] = meshgrid(lats,lons);

    % figure
    figure();
    worldmap(latlim,lonlim);
    hold on;
    
    load coastlines
    sst = squeeze(dataset(:,:,time));
    geoshow(lats_grid,lons_grid,sst,'DisplayType','surface');
    geoshow(coastlat,coastlon,'Color','k');
    colorbar;
    setm(gca,MLabelParallel="south",MLabelLocation=30);
end

% plot the top N EOFs
function plot_eofs(eofs,eigvs,N,latlim,lonlim,gdsize,model)
    for i=1:N
       plotDatasetSlice(eofs,latlim,lonlim,gdsize,i);
       title(sprintf("%s: EOF %d, (Explains %.2f%% of Variance)",model,i,eigvs(i)));
       mkdir(sprintf("figs/%s",model));
       caxis([-1 1])
       saveas(gcf,sprintf("figs/%s/eof_%d.png",model,i),"png");
    end
end

% plot the top N FFTs of PCs
function plot_fftpcs(pc,N,model)

    % init figure and set log scaled x axis
    figure();
    set(gca, 'XScale', 'log');
    hold on;
    
    % line styles
    ls = ["k-","r-","b-","g-","c-","m-"];
    for i=1:N
        ff = abs(fft(squeeze(pc(:,i)))); % fourier transform each PC
        s = size(ff);
        ff = ff./s(1); % normalize
        ff = ff(1:floor(s(1)/2)+1); % take the first half (FFT is symmetric, back half is "negative" frequencies)
        ff(2:end-1) = 2*ff(2:end-1); % correct normalization
        fs = (0:floor(s(1)/2))/s(1); % define frequency vector
        fs = 1./fs; % transform into periods
        semilogx(fs./12, ff, ls(i),'DisplayName',sprintf("PC %d",i)); % note that its plotted against period in years
        grid on;
    end
    xlabel("Period (Years)");
    ylabel("Magnitude of Fourier Coefficient");
    title(sprintf("%s: Fourier Transform of the first %d PCs",model,N));
    xlim([1,inf]); % sets lower bound at 1yr, upper bound auto determined by data
    legend('Location','northwest');
    saveas(gcf,sprintf("figs/%s/fft_pc.png",model),"png");
end

% plot eigenvalue spectrum
% this is figure 1 in the report
function plot_eigspec(dataset,model,trunc)
    figure();
    hold on;
    
    [~,eigvs] = get_eofs(dataset);
    
    plot(1:trunc, eigvs(1:trunc), "ko");
    xlabel("Eigenvalue");
    ylabel("Variance Explained (%)");
    title(sprintf("Top %d Eigenvalues from %s",trunc,model));
    
    saveas(gcf,"figs/eigenvalue_spectrum_ccsm3.png",'png');

end

% plot TPI
function plot_tpi(TPI,TPI_filt,model,start)
    f=figure();
    hold on;
    s = size(TPI,1);
    plot(((0:s-1)+start)/12, TPI,"-",'DisplayName',"TPI");
    plot(((0:s-1)+start)/12, TPI_filt, "-",'DisplayName',"Filtered TPI",'LineWidth',1.2);
    f.Position = [f.Position(1),f.Position(2),f.Position(3)*1.3,f.Position(4)*(2/3)];
    xlabel("Year")
    ylabel("TPI (°C)");
    title(sprintf("%s: IPO Tripole Index",model));
    legend('Location','northwest');
    xlim([start/12,inf]);
    saveas(gcf,sprintf("figs/%s/tpi.png",model),"png");
end

%% EOF function
function [eofs,eigvs,U] = get_eofs(dataset)
    % Based on a presentation by Shane Elipot, UMiami
    
    % we're taking in data as a lat x lon x time  array
    % but we want a N x time array
    % using reshape will let us put it into this format
    % and reshape is one to one invertible so itll put everything back exactly
    original_size = size(dataset);
    dataset = reshape(dataset, [original_size(1)*original_size(2),original_size(3)]);
    
    dataset(isnan(dataset)) = 0; % SVD doesn't like NaNs or INFs

    dataset = detrend(dataset,'constant');
    [P,L,U] = svd(dataset,'econ'); % computes the SVD
    eofs = P*L; % EOF matrix
    G = ctranspose(L)*L/(original_size(3)-1); % eigenvalue matrix
    eigvs = diag(G);
    [eigvs,inds] = sort((eigvs./sum(eigvs)).*100,"descend");
    eofs = eofs(:,inds); % rearranges EOFs to match sorted eigenvalues
    eofs = reshape(eofs, original_size);
    
    
    % normalize EOFs
    eofs = eofs./max(eofs,[],'all');

end

%% TPI function
function [TPI,TPI_filt] = get_tpi(dataset,latlim,lonlim)
    % going to "interpolate" to the desired boxes to make life easy
    
    % first need a geo reference object
    original_size = size(dataset);
    R = georefcells(latlim,lonlim,[original_size(2),original_size(1)],'ColumnsStartFrom','south', 'RowsStartFrom','west');

    % limits of each box
    latlim1 = [25 45];
    lonlim1 = [140 215];
    
    latlim2 = [-10 10];
    lonlim2 = [170 270];
    
    latlim3 = [-50 -15];
    lonlim3 = [150 200];
    
    dataset(isnan(dataset)) = 0; % mean doesn't like NaNs or INFs
    
    %% box one
    % define new lat/lon arrays and grids
    interp_lats = latlim1(1):1.5:latlim1(2);
    interp_lons = lonlim1(1):1.5:lonlim1(2);
    [interp_lats_grid,interp_lons_grid] = meshgrid(interp_lats,interp_lons);

    % define a place to put the interp results
    T1 = zeros(original_size(3),1);

    % do for all time (can do this in a better way?)
    for ti=1:original_size(3)
        curr = squeeze(dataset(:,:,ti));
        % runs fasters with these transposes
        T1(ti) = squeeze(mean(mean(geointerp(curr',R,interp_lats_grid',interp_lons_grid'),1),2));
    end
    
%     temp = geointerp(curr',R,interp_lats_grid',interp_lons_grid')';
%     temp = repmat(temp,[1,1,3]);
%     plotDatasetSlice(temp,latlim1,lonlim1,1.5,1);
    
    %% box two
    % define new lat/lon arrays and grids
    interp_lats = latlim2(1):1.5:latlim2(2);
    interp_lons = lonlim2(1):1.5:lonlim2(2);
    [interp_lats_grid,interp_lons_grid] = meshgrid(interp_lats,interp_lons);

    % define a place to put the interp results
    T2 = zeros(original_size(3),1);

    % do for all time (can do this in a better way?)
    for ti=1:original_size(3)
        curr = squeeze(dataset(:,:,ti));
        % runs fasters with these transposes
        T2(ti) = squeeze(mean(mean(geointerp(curr',R,interp_lats_grid',interp_lons_grid'),1),2));
    end
    
%     temp = geointerp(curr',R,interp_lats_grid',interp_lons_grid')';
%     temp = repmat(temp,[1,1,3]);
%     plotDatasetSlice(temp,latlim2,lonlim2,1.5,1);
    
    %% box three
    % define new lat/lon arrays and grids
    interp_lats = latlim3(1):1.5:latlim3(2);
    interp_lons = lonlim3(1):1.5:lonlim3(2);
    [interp_lats_grid,interp_lons_grid] = meshgrid(interp_lats,interp_lons);

    % define a place to put the interp results
    T3 = zeros(original_size(3),1);

    % do for all time (can do this in a better way?)
    for ti=1:original_size(3)
        curr = squeeze(dataset(:,:,ti));
        % runs fasters with these transposes
        T3(ti) = squeeze(mean(mean(geointerp(curr',R,interp_lats_grid',interp_lons_grid'),1),2));
    end
    
%     temp = geointerp(curr',R,interp_lats_grid',interp_lons_grid')';
%     temp = repmat(temp,[1,1,3]);
%     plotDatasetSlice(temp,latlim3,lonlim3,1.5,1);
    
    %% combine
    TPI = T2 - 0.5 *(T1+T3);
    
    [b,a] = cheby1(4,3,2*pi/(13*12),'low'); % define a 4th order low pass chebyshev filter with frequency 13 years
    TPI_filt = filter(b,a,TPI);

end
