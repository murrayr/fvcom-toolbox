% Find the closest elements/nodes/time to a location/time
%
% [nID eID tID] = find_FVCOM_location(varargin)
%
% DESCRIPTION:
%
% INPUT
%   Pass the variable names that you want to extract
%   filename, the netCDF filename
%   location, the location
%   time, the time
%
% EXAMPLE USAGE
%   
%   [nID eID tID] = read_netcdf_vars('filename', 'the_filename.nc', 'time', 'x', 'y', 'temp', 'location', [xloc yloc]);
%
%   
%   [nID eID tID] = read_netcdf_vars('filename', 'the_filename.nc', 'time', 'x', 'y', 'temp', 'time_val', the_time);
%
% Author(s)
%   Rory O'Hara Murray, Marine Scotland Science
%
% Revision history
%   v0 May 2014
%
function [nID eID tID] = find_FVCOM_location(varargin)

get_time = false; tID = 0;
get_space = false; nID = 0; eID = 0;

for ii=1:1:length(varargin)
    keyword  = lower(varargin{ii});
    if length(keyword)<8, continue; end
    switch(keyword(1:8))
        case 'filename'
            netcdf_filename = varargin{ii+1};
        case 'location'
            xy = varargin{ii+1};
            get_space = true;
        case 'time_val'
            model_time = varargin{ii+1};
            get_time = true;
    end
end

% Load model spatial data, and time series data
if get_space & get_time
    M = read_netcdf_vars('filename', netcdf_filename, 'x', 'y', 'xc', 'yc', 'time');
elseif get_space
    M = read_netcdf_vars('filename', netcdf_filename, 'x', 'y', 'xc', 'yc');
elseif get_time
    M = read_netcdf_vars('filename', netcdf_filename, 'time');
end

if get_space
    nID = fun_nearest2D_internal(xy(1), xy(2), M.x, M.y);
    eID = fun_nearest2D_internal(xy(1), xy(2), M.xc, M.yc);
end

if get_time
    tID = fun_nearest_internal(M.time', model_time);
end

% define function fun_nearest2D_internal
function PointID = fun_nearest2D_internal(xloc, yloc, x, y)
    for ii=1:length(xloc)
        radvec = sqrt( (xloc(ii)-x).^2 + (yloc(ii)-y).^2);
        [Distance(ii),PointID(ii)] = min(radvec);
    end
end

% define function fun_nearest_internal
function ib = fun_nearest_internal(b, a)
    m = size(a,2); n = size(b,2);
    [c,p] = sort([a,b]);
    q = 1:m+n; q(p) = q;
    t = cumsum(p>m);
    r = 1:n; r(t(q(m+1:m+n))) = r;
    s = t(q(1:m));
    id = r(max(s,1));
    iu = r(min(s+1,n));
    [d,it] = min([abs(a-b(id));abs(b(iu)-a)]);
    ib = id+(it-1).*(iu-id);
end

end