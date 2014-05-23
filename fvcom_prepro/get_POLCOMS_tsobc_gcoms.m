function Mobj = get_POLCOMS_tsobc_gcoms(Mobj, ts)
% Read temperature and salinity from the PML POLCOMS-ERSEM NetCDF model
% output files and interpolate onto the open boundaries in Mobj.
%
% function Mobj = get_POLCOMS_tsobc(Mobj, ts, polcoms_bathy, varlist)
%
% DESCRIPTION:
%    Interpolate temperature and salinity values onto the FVCOM open
%    boundaries at all sigma levels.
%
% INPUT:
%   Mobj    = MATLAB mesh structure which must contain:
%               - Mobj.siglayz - sigma layer depths for all model nodes.
%               - Mobj.lon, Mobj.lat - node coordinates (lat/long).
%               - Mobj.obc_nodes - list of open boundary node inidices.
%               - Mobj.nObcNodes - number of nodes in each open boundary.
%   ts      = Cell array of PML POLCOMS-ERSEM NetCDF file(s) in which 4D
%             variables of temperature and salinity (called 'ETWD' and
%             'x1XD') exist. Their shape should be (y, x, sigma, time).
%
% NOTES:
%
%   - If you supply multiple files in ts, there are a few assumptions:
%
%       - Variables are only appended if there are 4 dimensions; fewer than
%       that, and the values are assumed to be static across all the given
%       files (e.g. longitude, latitude etc.). The fourth dimension is
%       time.
%       - The order of the files given should be chronological.
% 
%   - The NetCDF files used here are those from the PML POLCOMS-ERSEM model
%   output.
%
% OUTPUT:
%    Mobj = MATLAB structure in which three new fields (called temperature,
%           salinity and ts_time). temperature and salinity have sizes
%           (sum(Mobj.nObcNodes), sigma, time). The time dimension is
%           determined based on the input NetCDF file. The ts_time variable
%           is just the input file times in Modified Julian Day.
%
% EXAMPLE USAGE
%    Mobj = get_POLCOMS_tsobc(Mobj, ts, depth)
%
% Author(s):
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2013-02-07 First version based on get_POLCOMS_tsobc.m.
%
%==========================================================================

subname = 'get_POLCOMS_tsobc_gcoms';

global ftbverbose;
if ftbverbose
    fprintf('\n')
    fprintf(['begin : ' subname '\n'])
    fprintf(['Expecting two files  : ' ts{1} '\n' ts{2} '\n'])
end

varlist = {'lon', 'lat', 'ETW', 'x1X', 'time'};

varlistdmeaneco3d = {'depth', 'pdepth'};

% Data format:
% 
%   pc.ETWD.data and pc.x1XD.data are y, x, sigma, time
% 
dd = 1;
catstr = 'pc=catstruct(';
for ff = 1:2:length(ts)-1
    pc = get_POLCOMS_netCDF(ts{ff}, varlist);
    pcdmean = get_POLCOMS_netCDF(ts{ff+1}, varlistdmeaneco3d);
    eval(['dump',num2str(dd),'=catstruct(pc,pcdmean);'])
    dd = dd+1;
    catstr = [catstr,'dump',num2str(dd),','];
end
eval([catstr(1:end-1),');'])
clear pcdmean dump*
[~, ~, nz, nt] = size(pc.ETW.data);

% Make rectangular arrays for the nearest point lookup.
[lon, lat] = meshgrid(pc.lon.data, pc.lat.data);

fvlon = Mobj.lon(Mobj.obc_nodes(Mobj.obc_nodes ~= 0));
fvlat = Mobj.lat(Mobj.obc_nodes(Mobj.obc_nodes ~= 0));

% Number of boundary nodes
nf = sum(Mobj.nObcNodes);
% Number of sigma layers.
fz = size(Mobj.siglayz, 2);

fvtemp = nan(nf, fz, nt); % FVCOM interpolated temperatures
fvsal = nan(nf, fz, nt); % FVCOM interpolated salinities

if ftbverbose
    tic
end
for t = 1:nt
    if ftbverbose
        fprintf('%s : %i of %i timesteps... ', subname, t, nt)
    end
    % Get the current 3D array of PML POLCOMS-ERSEM results.
    pctemp3 = pc.ETW.data(:, :, :, t);
    pcsalt3 = pc.x1X.data(:, :, :, t);
    
    % Preallocate the intermediate results arrays.
    itempz = nan(nf, nz);
    isalz = nan(nf, nz);
    idepthz = nan(nf, nz);
    
    for j = 1:nz
        % Now extract the relevant layer from the 3D subsets. Transpose the
        % data to be (x, y) rather than (y, x).
        pctemp2 = pctemp3(:, :, j)';
        pcsalt2 = pcsalt3(:, :, j)';
        pcdepth2 = squeeze(pc.depth.data(:, :, j, t))';
       
        % Create new arrays which will be flattened when masking (below).
        tpctemp2 = pctemp2;
        tpcsalt2 = pcsalt2;
        tpcdepth2 = pcdepth2;
        tlon = lon;
        tlat = lat;
        
        % Create and apply a mask to remove values outside the domain. This
        % inevitably flattens the arrays, but it shouldn't be a problem
        % since we'll be searching for the closest values in such a manner
        % as is appropriate for an unstructured grid (i.e. we're assuming
        % the PML POLCOMS-ERSEM data is irregularly spaced).
        mask = tpcdepth2 < -10e+10;
        tpctemp2(mask) = [];
        tpcsalt2(mask) = [];
        tpcdepth2(mask) = [];
        % Also apply the masks to the position arrays so we can't even find
        % positions outside the domain, effectively meaning if a value is
        % outside the domain, the nearest value to the boundary node will
        % be used.
        tlon(mask) = [];
        tlat(mask) = [];
        
        % Preallocate the intermediate results arrays.
        itempobc = nan(nf, 1);
        isalobc = nan(nf, 1);
        idepthobc = nan(nf, 1);
        
        % Speed up the tightest loop with a parallelized loop.
        parfor i = 1:nf
%        for i = 1:nf
            % Now we can do each position within the 2D layer.

            fx = fvlon(i);
            fy = fvlat(i);

            [~, ii] = sort(sqrt((tlon - fx).^2 + (tlat - fy).^2));
            % Get the n nearest nodes from PML POLCOMS-ERSEM data (more?
            % fewer?).
            ixy = ii(1:16);

            % Get the variables into static variables for the
            % parallelisation.
            plon = tlon(ixy);
            plat = tlat(ixy);
            ptemp = tpctemp2(ixy);
            psal = tpcsalt2(ixy);
            pdepth = tpcdepth2(ixy);
            
            % Use a triangulation to do the horizontal interpolation.
            tritemp = TriScatteredInterp(plon', plat', ptemp', 'natural');
            trisal = TriScatteredInterp(plon', plat', psal', 'natural');
            triz = TriScatteredInterp(plon', plat', pdepth', 'natural');
            itempobc(i) = tritemp(fx, fy);
            isalobc(i) = trisal(fx, fy);
            idepthobc(i) = triz(fx, fy);
            
            % Check all three, though if one is NaN, they all will be.
            if isnan(itempobc(i)) || isnan(isalobc(i)) || isnan(idepthobc(i))
                warning('FVCOM boundary node at %f, %f is outside the PML POLCOMS-ERSEM domain. Setting to the closest PML POLCOMS-ERSEM value.', fx, fy)
                itempobc(i) = tpctemp2(ii(1));
                isalobc(i) = tpcsalt2(ii(1));
                idepthobc(i) = tpcdepth2(ii(1));
            end
        end
        
        % Put the results in this intermediate array.
        itempz(:, j) = itempobc;
        isalz(:, j) = isalobc;
        idepthz(:, j) = idepthobc;
    end

    % Now we've interpolated in space, we can interpolate the z-values
    % to the sigma depths.
    oNodes = Mobj.obc_nodes(Mobj.obc_nodes ~= 0);

        % Preallocate the output arrays
        fvtempz = nan(nf, fz);
        fvsalz = nan(nf, fz);

        for pp = 1:nf
            % Get the FVCOM depths at this node
            tfz = Mobj.siglayz(oNodes(pp), :);
            % Now get the interpolated PML POLCOMS-ERSEM depth at this node
            tpz = idepthz(pp, :);

        % To ensure we get the full vertical expression of the vertical
        % profiles, we need to normalise the POLCOMS-ERSEM and FVCOM
        % depths to the same range. This is because in instances where
        % FVCOM depths are shallower (e.g. in coastal regions), if we
        % don't normalise the depths, we end up truncating the vertical
        % profile. This approach ensures we always use the full
        % vertical profile, but we're potentially squeezing it into a
        % smaller depth.
        A = max(tpz);
        B = min(tpz);
        C = max(tfz);
        D = min(tfz);
        norm_tpz = (((D - C) * (tpz - A)) / (B - A)) + C;

        % Get the temperature and salinity values for this node and
        % interpolate down the water column (from PML POLCOMS-ERSEM to
        % FVCOM). I had originally planned to use csaps for the
        % vertical interplation/subsampling at each location. However,
        % the demo of csaps in the MATLAB documentation makes the
        % interpolation look horrible (shaving off extremes). interp1
        % provides a range of interpolation schemes, of which pchip
        % seems to do a decent job of the interpolation (at least
        % qualitatively).
        if ~isnan(tpz)
            fvtempz(pp, :) = interp1(norm_tpz, itempz(pp, :), tfz, 'linear', 'extrap');
            fvsalz(pp, :) = interp1(norm_tpz, isalz(pp, :), tfz, 'linear', 'extrap');
        else
            warning('Should never see this... ') % because we test for NaNs when fetching the values.
            warning('FVCOM boundary node at %f, %f is outside the PML POLCOMS-ERSEM domain. Skipping.', fvlon(pp), fvlat(pp))
            continue
        end
    end
    
    % The horizontally- and vertically-interpolated values in the final
    % FVCOM results array.
    fvtemp(:, :, t) = fvtempz;
    fvsal(:, :, t) = fvsalz;
    
    if ftbverbose
        fprintf('done.\n')
    end
end
if ftbverbose
    toc
end

Mobj.temperature = fvtemp;
Mobj.salt = fvsal;

% Convert the current times to Modified Julian Day (this is a bit ugly).
pc.time.all = strtrim(regexp(pc.time.units, 'since', 'split'));
pc.time.datetime = strtrim(regexp(pc.time.all{end}, ' ', 'split'));
pc.time.ymd = str2double(strtrim(regexp(pc.time.datetime{1}, '-', 'split')));
pc.time.hms = str2double(strtrim(regexp(datestr(datenum(pc.time.ymd),'HH:MM:SS'), ':', 'split')));

Mobj.ts_times = greg2mjulian(...
    pc.time.ymd(1), ...
    pc.time.ymd(2), ...
    pc.time.ymd(3), ...
    pc.time.hms(1), ...
    pc.time.hms(2), ...
    pc.time.hms(3)) + (pc.time.data / 3600 / 24);

if ftbverbose
    fprintf(['end   : ' subname '\n'])
end
