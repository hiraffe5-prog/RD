function [center_lon, center_lat, center_h, row_c, col_c, lon1, lat1] = get_dem_center(DEM, lon_lat_range)

[Nrow, Ncol] = size(DEM);

lon = linspace(lon_lat_range(1), lon_lat_range(2), Ncol);
lat = linspace(lon_lat_range(3), lon_lat_range(4), Nrow);

row_c = round((Nrow+1)/2);
col_c = round((Ncol+1)/2);

center_lon = lon(col_c);
center_lat = lat(row_c);
center_h   = DEM(row_c, col_c);

[lon1, lat1] = meshgrid(lon, lat);

end