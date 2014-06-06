function [out_east,out_north] = my_project(in_east,in_north,direction) 

% Sample user-defined projection and inverse projection of (lon,lat) to (x,y) 
% Copy to my_project (not a member of the toolbox) and modify to suite you
%
% function [out_east,out_north] = my_project(in_east,in_north,direction) 
%
% DESCRIPTION:
%    Define projections between geographical and Euclidean coordinates 
%
% INPUT: 
%   in_east   = 1D vector containing longitude (forward) x (reverse)
%   in_north  = 1D vector containing latitude  (forward) y (reverse)
%   direction = ['forward' ;  'inverse']
%           
% OUTPUT:
%   (lon,lat) or (x,y) depending on choice of forward or reverse projection
%
% EXAMPLE USAGE
%    [lon,lat] = my_project(x,y,'reverse') 
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

%subname = 'my_project';
%fprintf('\n')
%fprintf(['begin : ' subname '\n'])

%------------------------------------------------------------------------------
% Parse input arguments
%------------------------------------------------------------------------------

ProjectDirection = 'forward';

if(direction == 'forward')
	ProjectDirection = 'forward';
        lon = in_east;
        lat = in_north;
else
	ProjectDirection = 'inverse';
        x = in_east;
        y = in_north;
end;


if(ProjectDirection == 'forward')
	fprintf('Projecting from (lon,lat) to (x,y) UTM 30 N\n');
    [x, y] = ll2utm(lon, lat, 30);
else
	fprintf('Inverse Projecting from (x,y) UTM 30 N to (lon,lat)\n')
	[lon, lat]=utm2ll(x,y,30); 
end;


% set the output
if(ProjectDirection == 'forward')
  out_east = x;
  out_north = y;
else
  out_east = lon;
  out_north = lat;
end;

