% Function to define the metadata for the physics variables included in the
% nesting forcing files. This is can be used to define phsyics metadata
% before running write_FVCOM_nested_forcing.m or it is called by this
% function if no metadata are supplied.
%
% Rory O'Hara Murray, 08/12/2023
%
function var_names = define_physics_nesting_var_names
    name = {'zeta'; 'ua'; 'va'; 'u'; 'v'; 'temp'; 'salinity'; 'hyw'};
    long_name = {'Water Surface Elevation'; 'Vertically Averaged x-velocity'; 'Vertically Averaged y-velocity'; 'Eastward Water Velocity'; 'Northward Water Velocity'; 'Temperature'; 'Salinity'; 'hydro static vertical velocity'};
    standard_name = {'sea_surface_height_above_geoid'; 'vertically-Averaged-x-velocity'; 'vertically-Averaged-y-velocity'; 'eastward_sea_water_velocity'; 'northward_sea_water_velocity'; 'sea_water_temperature'; 'sea_water_salinity'; 'hydro_static_vertical_velocity'};
    units = {'meters'; 'meters s-1'; 'meters s-1'; 'meters s-1'; 'meters s-1'; 'degrees Celcius'; '1e-3'; 'meters s-1'};
    grid = {'Bathymetry_Mesh'; 'fvcom_grid'; 'fvcom_grid'; 'fvcom_grid'; 'fvcom_grid'; 'fvcom_grid'; 'fvcom_grid'; 'fvcom_grid'};
    type = {'data'; 'data'; 'data'; 'data'; 'data'; 'data'; 'data'; 'data'};
    location = {'node'; 'face'; 'face'; 'face'; 'face'; 'node'; 'node'; 'node'};
    coordinates = {'time lat lon'; 'time latc lonc'; 'time latc lonc'; 'time siglay latc lonc'; 'time siglay latc lonc'; 'time siglay lat lon'; 'time siglay lat lon'; 'time siglev lat lon'};
%    dimension = {'2D'; '2D'; '2D'; '3D'; '3D'; '3D'; '3D'; '3D'};

    var_names_table = table(name, long_name, standard_name, units, grid, type, location, coordinates);
    var_names = table2struct(var_names_table);
end
