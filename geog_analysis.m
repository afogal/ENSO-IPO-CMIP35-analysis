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

%% main sequence

% define the box we want
latlim = [-50 50];
lonlim = [120 300];

% load CSIRO CMIP3 data
CSIRO35_sst = load_CSIRO35(latlim,lonlim);

% plot it
plotDatasetSlice(CSIRO35_sst,latlim,lonlim,1.5,10);


% end main sequence

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HERE BE FUNCTIONS
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
    
    save("CMIPdata/CSIRO35_sst.mat", "CSIRO35_sst");

end

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
    figure(1);
    worldmap(latlim,lonlim);
    hold on;
    
    load coastlines
    sst = squeeze(dataset(:,:,time));
    geoshow(lats_grid,lons_grid,sst,'DisplayType','surface');
    geoshow(coastlat,coastlon,'Color','k');
    colorbar;
    setm(gca,MLabelParallel="south",MLabelLocation=30);
end
