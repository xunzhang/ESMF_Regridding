#include "grid.h"
#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include <stdlib.h>
#include <math.h>


#define NC_ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(2);}


#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define ZERO  0.0
#define TWO 2
#define HALF 0.5
#define PI2   TWO*PI
#define PIH   HALF*PI
#define MAX_VALUE    9999
#define MIN_VALUE   -9999
#define LATITUDE 1
#define LATLON    2

grid::grid(int out_npnts, int out_nlats, int out_nlons, double *out_center_lat_coords, double *out_center_lon_coords, bool *out_mask)
{
	double last_lat;
	double *tmp_lon_array1;
	double *tmp_lon_array2;
	double *tmp_double;
	int i, j, k, m, n, count, indx;
	int *tmp_nlat_each_lon1;
	int *tmp_nlat_each_lon2;
	int *tmp_int;
	int *tmp_current_indx_each_lon;
	double max_lat;
	double min_lat;

	//ndim = out_ndim;
	npnts = out_npnts;
	nlats = out_nlats;
	nlons = out_nlons;
	center_lat_coords = out_center_lat_coords;
	center_lon_coords = out_center_lon_coords;
	mask = out_mask;


	/* calculate the maximum and minimum latitudes of each row */
	max_lat_each_row = new double [nlats];
	min_lat_each_row = new double [nlats];
	nlon_each_lat = new int [nlats];
	nlat_each_lon = new int [nlons];
	different_lons = new double [nlons];
	different_lats = new double [nlats];
	reverse_indx = new int [npnts];


	//printf("nlons is %d\n", nlons);

	if (npnts == nlats * nlons) {		// logical rectangular grid
		for (i = 0; i < nlats; i ++) {
			nlon_each_lat[i] = nlons;
			max_lat = -100000;
			min_lat = 100000;
			for (j = 0; j < nlons; j ++) {
				if (center_lat_coords[i * nlons + j] > max_lat)
					max_lat = center_lat_coords[i * nlons + j];
				if (center_lat_coords[i * nlons + j] < min_lat)
					min_lat = center_lat_coords[i * nlons + j];
			}
			max_lat_each_row[i] = max_lat;
			min_lat_each_row[i] = min_lat;
		}
		for (i = 0; i < nlats; i ++)
			if (max_lat_each_row[i] != min_lat_each_row[i])
				break;
		if (i != nlats) {
			type_grid = 1;
			//printf("irregular grid\n");
		}
		else {
			type_grid = 0;
			//printf("regular grid\n");
		}

		//printf("grid init done!~!~\n"); 
	}
	else {
		/*calculate the number of longitudes in each latitude*/
		for (i = 0; i < nlats; i ++)
			nlon_each_lat[i] = 0;
		last_lat = center_lat_coords[0];
		indx = 0;
		for (i = 0, count = 0; i < npnts; i ++) {
			if (last_lat == center_lat_coords[i])
				count ++;
			else {
				nlon_each_lat[indx ++] = count; 
				last_lat = center_lat_coords[i];
				count = 1;
			}
		}
		nlon_each_lat[indx ++] = count; 

		/*calculate the different latitudes in the grid*/
		for (i = 0, count = 0; i < nlats; i ++) {
			different_lats[i] = center_lat_coords[count];
			count += nlon_each_lat[i];
		}

		/*calculate the number of latitudes in each longitude and the different longitudes in the grid*/
		tmp_lon_array1 = new double [nlons];
		tmp_lon_array2 = new double [nlons];
		tmp_nlat_each_lon1 = new int [nlons];
		tmp_nlat_each_lon2 = new int [nlons];
		for (i = 0; i < nlons; i ++)
			tmp_nlat_each_lon1[i] = 0;
		for (i = 0; i < nlon_each_lat[0]; i ++) {
			tmp_nlat_each_lon1[i] = 1;
			tmp_lon_array1[i] = center_lon_coords[i];
		}
		count = nlon_each_lat[0];
		for (i = 1, j = nlon_each_lat[0]; i < nlats; i ++) {
			k = 0;
			m = 0;
			n = 0;
			while (k < nlon_each_lat[i] && m < count) {
				if (tmp_lon_array1[m] == center_lon_coords[j + k]) {
					tmp_nlat_each_lon2[n] = tmp_nlat_each_lon1[m] + 1;
					tmp_lon_array2[n ++] = tmp_lon_array1[m];
					m ++;
					k ++;
				}
				else if (tmp_lon_array1[m] < center_lon_coords[j + k]) {
					tmp_nlat_each_lon2[n] = tmp_nlat_each_lon1[m];
					tmp_lon_array2[n ++] = tmp_lon_array1[m];
					m ++;
				}
				else {
					tmp_nlat_each_lon2[n] = 1;
					tmp_lon_array2[n ++] = center_lon_coords[j + k];
					k ++;
				}
			}
			for (; k < nlon_each_lat[i]; k ++) {
				tmp_nlat_each_lon2[n] = 1;
				tmp_lon_array2[n ++] = center_lon_coords[j + k];
			}
			for (; m < count; m ++) {
				tmp_nlat_each_lon2[n] = tmp_nlat_each_lon1[m];
				tmp_lon_array2[n ++] = tmp_lon_array1[m];
			}
			count = n;
			tmp_double = tmp_lon_array1;
			tmp_lon_array1 = tmp_lon_array2;
			tmp_lon_array2 = tmp_double;
			tmp_int = tmp_nlat_each_lon1;
			tmp_nlat_each_lon1 = tmp_nlat_each_lon2;
			tmp_nlat_each_lon2 = tmp_int;
			j += nlon_each_lat[i];
		}
		for (i = 0; i < nlons; i ++) {
			nlat_each_lon[i] = tmp_nlat_each_lon1[i];
			different_lons[i] = tmp_lon_array1[i];
		}

		/*check whether the grid is retangular*/
		is_rect = true;
		for (i = 0; i < nlats; i ++)
			if (nlons != nlon_each_lat[i])
				is_rect = false;
		for (i = 0; i < nlons; i ++)
			if (nlats != nlat_each_lon[i])
				is_rect = false;

		/*calculate the reverse index of each point when transforming the lat-major matrix into lon-major matrix*/
		tmp_current_indx_each_lon = new int [nlons];
		tmp_current_indx_each_lon[0] = 0;
		for (i = 1; i < nlons; i ++)
			tmp_current_indx_each_lon[i] = tmp_current_indx_each_lon[i - 1] + nlat_each_lon[i];
		for (i = 0, j = 0; i < nlats; i ++) {
			k = 0;
			m = 0;
			n = 0;
			while (k < nlon_each_lat[i] && m < nlons) {
				if (tmp_lon_array1[m] == center_lon_coords[j + k]) {
					reverse_indx[j + k] = tmp_current_indx_each_lon[m];
					tmp_current_indx_each_lon[m] ++;
					m ++;
					k ++;
				}
				else if (tmp_lon_array1[m] < center_lon_coords[j + k]) 
					m ++;
				else {
					printf("never happen case when building reverse indexes for grid\n");
				}
			}
			j += nlon_each_lat[i];
		}	
		delete [] tmp_current_indx_each_lon;
		delete [] tmp_lon_array1;
		delete [] tmp_lon_array2;
		delete [] tmp_nlat_each_lon1;
		delete [] tmp_nlat_each_lon2;
	}
}

grid::~grid()
{
	//delete [] center_lon_coords;
	//delete [] center_lat_coords;
	//	if (center_heights != NULL)
	//		delete [] center_heights;
	//	delete [] vertex_lat_coords;
	//	delete [] vertex_lon_coords;
	//	if (vertex_heights != NULL)
	//		delete [] vertex_heights;
	//	delete [] area_in;
	//	delete [] areas;
	//	for (int i = 0; i < 1; i ++) 
	//		if (mask[i] != NULL) {
	//			delete [] mask[i];
	//			delete [] model_frac[i];
	//		}
	//	delete [] mask;
	//	delete [] model_frac;
	delete [] nlat_each_lon;
	delete [] nlon_each_lat;
	delete [] reverse_indx;
	delete [] different_lats;
	delete [] max_lat_each_row;
	delete [] min_lat_each_row;
	delete [] different_lons;

	//adding
	// delete [] bound_box; 						
	//	delete [] dims;	                                  
	//delete [] bin_addr;                             
	//delete [] bin_lats;           
	//delete [] bin_lons;	
}

#ifdef DEBUG_GRID
void grid::check_grid_info(int check_type)
{
	int i, j;
	double *tmp_lat_coords;
	double *tmp_lon_coords;

	/*check whether the grid is lat-major and sorted after loading it*/
	if (check_type == 0) {
		for (i = 1; i < npnts; i ++) {
			if (center_lat_coords[i] < center_lat_coords[i - 1]) {
				printf("grid is not lat-major and sorted: case 1\n");
				break;
			}
			if (center_lat_coords[i] == center_lat_coords[i - 1] &&
					center_lon_coords[i] < center_lon_coords[i - 1]) {
				printf("grid is not lat-major and sorted: case 2 %d\n", i);
				break;
			}
		}
	}

	/*check the correctness of nlats, nlons and reverse matrix*/
	if (check_type == 1) {
		if (!is_rect)
			printf("%s is not a rect grid\n", my_name);
		tmp_lat_coords = new double [npnts];
		tmp_lon_coords = new double [npnts];
		for (i = 0; i < nlats; i ++)
			if (nlon_each_lat[i] == 0) {
				printf("nlats of grid is wrong\n");
				break;
			}
		for (i = 0; i < nlons; i ++)
			if (nlat_each_lon[i] == 0) {
				printf("nlons of grid is wrong\n");
				break;
			}
		for (i = 0; i < npnts; i ++) {
			tmp_lat_coords[reverse_indx[i]] = tmp_lat_coords[i];
			tmp_lon_coords[reverse_indx[i]] = tmp_lon_coords[i];
		}
		/*check whether the reversed grid is lon-major and sorted*/
		for (i = 1; i < npnts; i ++) {
			if (tmp_lon_coords[i] < tmp_lon_coords[i - 1]) {
				printf("the reversed grid is not lon-major and sorted: case 1\n");
				break;
			}
			if (tmp_lon_coords[i] == tmp_lon_coords[i - 1] &&
					tmp_lat_coords[i] < tmp_lat_coords[i - 1]) {
				printf("the reversed grid is not lon-major and sorted: case 2\n");
				break;
			}
		}
		//delete [] tmp_lat_coords;
		//delete [] tmp_lon_coords;
	}
}
#endif


