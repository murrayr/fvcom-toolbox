function Mobj = get_POLCOMS_tsobc(Mobj, ts)
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
%       files (e.g. longitude, latitude). The fourth dimension is time.
%       - The order of the files given should be chronological.
% 
%   - The NetCDF files used here are those from the PML POLCOMS-ERSEM model
%   output.
%
% OUTPUT:
%    Mobj = MATLAB structure in which three new fields (called temperature,
%           salinity and ts_times). temperature and salinity have sizes
%           (sum(Mobj.nObcNodes), sigma, time). The time dimension is
%           determined based on the input NetCDF file. The ts_time variable
%           is just the input file times in Modified Julian Day.
%
% EXAMPLE USAGE
%    Mobj = get_POLCOMS_tsobc(Mobj, ts)
%
% Author(s):
%    Pierre Cazenave (Plymouth Marine Laboratory)
%
% Revision history
%    2013-02-07 First version.
%    2013-02-27 Change the vertical interpolation to be scaled within the
%    POLCOMS-ERSEM depth range for the current node. The net result is that
%    the vertical profiles are squashed or stretched to fit within the
%    FVCOM depths. This means the full profile structure is maintained in
%    the resulting FVCOM boundary input despite the differing depths at the
%    FVCOM boundary node.
%    2013-06-03 Fix some bugs in the way the open boundary node values were
%    stored (the order in which they were stored did not match the order of
%    the nodes in Casename_obc.dat). Also fix the order of the vertically
%    interpolated values so that FVCOM starts at the surface instead of
%    mirroring POLCOMS' approach (where the first value is the seabed). The
%    effect of these two fixes (nodes and vertical) should match what FVCOM
%    expects. Also add a set of figures (commented out) at the end for
%    diagnostic purposes.
%
%==========================================================================

subname = 'get_POLCOMS_tsobc';

global ftbverbose;
if ftbverbose
    fprintf('\n')
    fprintf(['begin : ' subname '\n'])
end

wasOpened = false;
if license('test', 'Distrib_Computing_Toolbox')
    % We have the Parallel Computing Toolbox, so launch a bunch of workers.
    if matlabpool('size') == 0
        % Force pool to be local in case we have remote pools available.
        matlabpool open local
        wasOpened = true;
    end
end

varlist = {'lon', 'lat', 'ETWD', 'x1XD', 'time', 'depth', 'pdepthD'};

% Data format:
% 
%   pc.ETWD.data and pc.x1XD.data are y, x, sigma, time
% 
pc = get_POLCOMS_netCDF(ts, varlist);

[~, ~, nz, nt] = size(pc.ETWD.data);

% Make rectangular arrays for the nearest point lookup.
[lon, lat] = meshgrid(pc.lon.data, pc.lat.data);

%oNodes = Mobj.obc_nodes(Mobj.obc_nodes ~= 0);
% Change the way the nodes are listed to match the order in the
% Casename_obc.dat file.
tmpObcNodes = Mobj.obc_nodes';
oNodes = tmpObcNodes(tmpObcNodes ~= 0)';

fvlon = Mobj.lon(oNodes);
fvlat = Mobj.lat(oNodes);

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
    pctemp3 = pc.ETWD.data(:, :, :, t);
    pcsalt3 = pc.x1XD.data(:, :, :, t);

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
        mask = tpcdepth2 > 20000;
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

        % Put the results in the intermediate array.
        itempz(:, j) = itempobc;
        isalz(:, j) = isalobc;
        idepthz(:, j) = idepthobc;

    end

    % Now we've interpolated in space, we can interpolate the z-values
    % to the sigma depths.

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
            fvtempz(pp, :) = interp1(norm_tpz, itempz(pp, :), tfz, 'pchip', 'extrap');
            fvsalz(pp, :) = interp1(norm_tpz, isalz(pp, :), tfz, 'pchip', 'extrap');
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
pc.time.hms = str2double(strtrim(regexp(pc.time.datetime{2}, ':', 'split')));

Mobj.ts_times = greg2mjulian(...
    pc.time.ymd(1), ...
    pc.time.ymd(2), ...
    pc.time.ymd(3), ...
    pc.time.hms(1), ...
    pc.time.hms(2), ...
    pc.time.hms(3)) + (pc.time.data / 3600 / 24);

% Close the MATLAB pool if we opened it.
if wasOpened
    matlabpool close
end

if ftbverbose
    fprintf(['end   : ' subname '\n'])
end


%%
% Plot a vertical profile for a boundary node (for my Irish Sea case, this
% is one of the ones along the Celtic Sea boundary). Also plot the
% distribution of interpolated values over the POLCOMS data. Add the
% location of the vertical profile (both FVCOM and POLCOMS) to the plot.
% nn = 55;   % open boundary index
% tt = 1;    % time index
%
% % Get the corresponding indices for the POLCOMS data
% [~, xidx] = min(abs(lon(1, :) - fvlon(nn)));
% [~, yidx] = min(abs(lat(:, 1) - fvlat(nn)));
%
% % close all
%
% figure
% clf
% subplot(2,2,1)
% plot(Mobj.temperature(nn, :, 1), Mobj.siglayz(oNodes(nn), :), 'x-')
% xlabel('Temperature (^{\circ}C)')
% ylabel('Depth (m)')
% title('FVCOM')
%
% subplot(2,2,2)
% % Although POLCOMS stores its temperature values from seabed to surface,
% % the depths are stored surface to seabed. Nice.
% plot(squeeze(pc.ETWD.data(xidx, yidx, :, 1)), squeeze(pc.depth.data(xidx, yidx, :, 1)), 'rx-')
% xlabel('Temperature (^{\circ}C)')
% ylabel('Depth (m)')
% title('POLCOMS')
%
% subplot(2,2,3)
% plot(Mobj.temperature(nn, :, tt), 1:fz, 'x-')
% xlabel('Temperature (^{\circ}C)')
% ylabel('Array index')
% title('FVCOM')
%
% subplot(2,2,4)
% plot(squeeze(pc.ETWD.data(xidx, yidx, :, tt)), 1:nz, 'rx-')
% xlabel('Temperature (^{\circ}C)')
% ylabel('Array index')
% title('POLCOMS')
%
% % Figure to check everything's as we'd expect. Plot first time step with
% % the POLCOMS surface temperature as a background with the interpolated
% % boundary node surface values on top.
%
% figure
% clf
% % Plot POLCOMS surface data (last sigma layer)
% dx = mean(diff(pc.lon.data));
% dy = mean(diff(pc.lat.data));
% pcolor(pc.lon.data - (dx / 2), pc.lat.data - (dy / 2), ...
%     squeeze(pc.ETWD.data(:, :, 1, tt))')
% shading flat
% axis('equal', 'tight')
% daspect([1.5, 1, 1])
% hold on
% % Add the interpolated surface data (first sigma layer)
% scatter(Mobj.lon(oNodes), Mobj.lat(oNodes), 40, Mobj.temperature(:, 1, tt), 'filled', 'MarkerEdgeColor', 'k')
% axis([min(Mobj.lon(oNodes)), max(Mobj.lon(oNodes)), min(Mobj.lat(oNodes)), max(Mobj.lat(oNodes))])
% caxis([6, 20])
% plot(lon(yidx, xidx), lat(yidx, xidx), 'rs') % polcoms is all backwards
% plot(Mobj.lon(oNodes(nn)), Mobj.lat(oNodes(nn)), 'wo')
% colorbar
