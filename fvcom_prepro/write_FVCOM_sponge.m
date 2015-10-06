function write_FVCOM_sponge(Mobj,filename) 

% Write FVCOM format sponge layer file  
%
% function write_FVCOM_sponge(Mobj,filename)
%
% DESCRIPTION:
%    Generate an ascii FVCOM 3.x format sponge file from Mesh object
%
% INPUT 
%   Mobj     = Mesh object
%   filename = FVCOM sponge file name
%
% OUTPUT:
%    FVCOM sponge file: filename
%
% EXAMPLE USAGE
%    write_FVCOM_sponge(Mobj,'tst_spg.dat')   
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%    Karen Thurston (National Oceanography Centre, Liverpool)
%
% Revision history
%   2013-01-18  Added support for variable sponge radius
%   2014-10-28 Added support for variable sponge damping coefficient, by
%   assuming the size of the sponge_fac and sponge_rad arrays are equal the
%   number of sponge nodes by default.
%   
%==============================================================================
subname = 'write_FVCOM_sponge';
global ftbverbose 
if(ftbverbose)
  fprintf('\n'); fprintf(['begin : ' subname '\n']);
end;

%------------------------------------------------------------------------------
% Parse input arguments
%------------------------------------------------------------------------------
if(exist('Mobj')*exist('filename')==0)
	error('arguments to write_FVCOM_sponge are incorrect')
end;

%------------------------------------------------------------------------------
% Correct possible errors arrising from previouse Mobj storage methods
%------------------------------------------------------------------------------

% make sure sponge_fac and sponge_rad are the right size
if size(Mobj.sponge_fac,2)==1
    Mobj.sponge_fac = Mobj.sponge_fac(:,1)*ones(1,max(Mobj.nSpongeNodes));
elseif size(Mobj.sponge_fac,2)<max(Mobj.nSpongeNodes)
    error('sponge_fac is an incompatible size, check it''s been written correctly.')
end
if size(Mobj.sponge_rad,2)==1
    Mobj.sponge_rad = Mobj.sponge_rad(:,1)*ones(1,max(Mobj.nSpongeNodes));
elseif size(Mobj.sponge_rad,2)<max(Mobj.nSpongeNodes)
    error('sponge_rad is an incompatible size, check it''s been written correctly.')
end

% If there are zeros across the sponge_fac and sponge_rad arrays then
% assume they are constant values and fill. This may have happened if one
% sponge layer has constant values, and the second has variable values.
for n=1:length(Mobj.nSpongeNodes)
    if sum(Mobj.sponge_fac(n,2:Mobj.nSpongeNodes(n)))==0
        Mobj.sponge_fac(n,2:Mobj.nSpongeNodes(n)) = Mobj.sponge_fac(n,1);
    end
    if sum(Mobj.sponge_rad(n,2:Mobj.nSpongeNodes(n)))==0
        Mobj.sponge_rad(n,2:Mobj.nSpongeNodes(n)) = Mobj.sponge_rad(n,1);
    end
end

%------------------------------------------------------------------------------
% Dump the file
%------------------------------------------------------------------------------
if(ftbverbose); fprintf('writing FVCOM spongefile %s\n',filename); end;
fid = fopen(filename,'w');

if(Mobj.nSponge==0)
	fprintf(fid,'Sponge Node Number = %d\n',0);
else
	Total_Sponge = sum(Mobj.nSpongeNodes(1:Mobj.nSponge));
	fprintf(fid,'Sponge Node Number = %d\n',Total_Sponge);
	for i=1:Mobj.nSponge
        for j=1:Mobj.nSpongeNodes(i)
            fprintf(fid,'%d %f %f \n',Mobj.sponge_nodes(i,j),Mobj.sponge_rad(i,j),Mobj.sponge_fac(i,j));
        end;
	end;
end;
fclose(fid);
		
if(ftbverbose)
  fprintf(['end   : ' subname '\n'])
end;

