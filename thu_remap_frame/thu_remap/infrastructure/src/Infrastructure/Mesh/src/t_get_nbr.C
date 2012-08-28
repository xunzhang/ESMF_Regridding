#include "t_get_nbr.h"
#include "remap_base.h"
#include "t_sort.h"
#include "grid.h"
#include <stdio.h>
#include <math.h>

#define EARTH_R ((double) 6378.137)

double cal_pnts_dist(double lon_point1, 
                             double lat_point1, 
                             double lon_point2, 
                             double lat_point2)
{
	double rad_lat1 = RADIAN(lat_point1);
	double rad_lat2 = RADIAN(lat_point2);
	double lon_diff = RADIAN(lon_point1) - RADIAN(lon_point2);
	double s = acos(cos(rad_lat1) * cos(rad_lat2) * cos(lon_diff) + sin(rad_lat1) * sin(rad_lat2));
	
	return s;
}


void cal_rect_border_of_dist(double lon_coords_src, 
                                          double lat_coords_src, 
                                          double &lon_dist_border_low, 
                                          double &lon_dist_border_high, 
                                          double &lat_dist_border_low, 
                                          double &lat_dist_border_high,
                                          double dist_threshold)
{
	double r;		// radius of current lat
	
	lat_dist_border_low = lat_coords_src - (180 * dist_threshold) / (PI * EARTH_R);
	lat_dist_border_high = lat_coords_src + (180 * dist_threshold) / (PI * EARTH_R);

	if (lat_dist_border_high > 90) 	{
		lat_dist_border_high = 90;
		lon_dist_border_low = 0;
		lon_dist_border_high = 360;
	}
	else if (lat_dist_border_low < -90) 	{
		lat_dist_border_low = -90;
		lon_dist_border_low = 0;
		lon_dist_border_high = 360;
	}
	else {
		r = EARTH_R * cos(lat_coords_src*PI / 180);
		if (r < 0)
			r = -r;
		lon_dist_border_low = lon_coords_src - (360 * dist_threshold) / (2 * PI * r);
		if (lon_dist_border_low < 0) 
			lon_dist_border_low += 360;
		lon_dist_border_high = lon_coords_src + (360 * dist_threshold) / (2 * PI * r);
		if (lon_dist_border_high > 360) 
			lon_dist_border_high -= 360; 
	} 
}

bool search_nearest_pnts(double lat_dst,
                                   double lon_dst,
                                   grid *grid_src,
                                   double dist_threshold,
                                   int num_nearest_pnts,
                                   int *point_indx_src,
                                   double *weights,
                                   double power_p)
{
	double lon_dist_border_low, lon_dist_border_high;
	double lat_dist_border_low, lat_dist_border_high;
	int i, j, m, old_j;
	double dist;
	bool potential_neighbor;
	double sum;
	int npts_src = grid_src->get_num_points();
	double *lat_coords_src = grid_src->get_lat_coords();
	double *lon_coords_src = grid_src->get_lon_coords();
	int *nlon_each_lat_src = grid_src->get_nlon_each_lat();
	static int *local_point_indx_src = NULL;
	static double *local_weights = NULL;
	static int buf_size = 4096;
	int count = 0;
	bool *mask_src = grid_src->get_mymask();
	double *max_lat_each_row_src = grid_src->get_max_lat_each_row();
	double *min_lat_each_row_src = grid_src->get_min_lat_each_row();
	int nlats_src = grid_src->get_num_lats();
	int row_src;

	if (local_point_indx_src == NULL) {
		local_point_indx_src = new int [buf_size];
		local_weights = new double [buf_size];	
	}
	
	cal_rect_border_of_dist(lon_dst, lat_dst, 
	                                    lon_dist_border_low, 
	                                    lon_dist_border_high, 
	                                    lat_dist_border_low, 
	                                    lat_dist_border_high,
	                                    dist_threshold);

	for (j = 0, row_src = 0; row_src < nlats_src; row_src ++) {
		old_j = j;
		if (max_lat_each_row_src[row_src] >= lat_dist_border_low && min_lat_each_row_src[row_src] <= lat_dist_border_high) {
			for (; j < old_j + nlon_each_lat_src[row_src]; j ++) {
				if (!mask_src[j]) 
					continue;
				potential_neighbor = false;
				if (lon_dist_border_low > lon_dist_border_high &&
				    (lon_coords_src[j] >= lon_dist_border_low 
				    || lon_coords_src[j] <= lon_dist_border_high))
					potential_neighbor = true;
				if (lon_dist_border_low < lon_dist_border_high &&
				    lon_coords_src[j] >= lon_dist_border_low 
				    && lon_coords_src[j] <= lon_dist_border_high)
					potential_neighbor = true;
				if (potential_neighbor)	{
					dist = cal_pnts_dist(lon_coords_src[j], 
					                           lat_coords_src[j], 
					                           lon_dst, 
					                           lat_dst);

					if (dist <= dist_threshold / EARTH_R) {
						if (dist == 0)	{
							for (m = 0; m < num_nearest_pnts; m ++) {
								weights[m] = 0.0;
								point_indx_src[m] = j;
							}
							weights[0] = 1.0;
							return true;
						}
						local_point_indx_src[count] = j;
						local_weights[count] = 1.0 / pow(dist, power_p);
						count ++;
						if (count >= buf_size) {
							double *new_weight_buf = new double [buf_size * 2];
							int *new_indx_buf = new int [buf_size * 2];
							for (i = 0; i < buf_size; i ++) {
								new_weight_buf[i] = local_weights[i];
								new_indx_buf[i] = local_point_indx_src[i];
							}
							delete [] local_weights;
							delete [] local_point_indx_src;
							local_weights = new_weight_buf;
							local_point_indx_src = new_indx_buf;
							buf_size *= 2;
						}
					}
				}
			}
		}
		j = old_j + nlon_each_lat_src[row_src];
	}
	if (count < num_nearest_pnts) {
		return false;
	}
	d_quick_sort(local_weights, local_point_indx_src, 0, count - 1);
	for (sum = 0.0, j = 0; j < num_nearest_pnts; j ++) 
		sum += local_weights[j];
	for (j = 0; j < num_nearest_pnts; j ++) {
		point_indx_src[j] = local_point_indx_src[j];
		weights[j] = local_weights[j] / sum;
	}
	
	return true;
}

double cal_dist_threshold(double lat_dst,
                                   grid *grid_src,
                                   int num_nearest_pnts)
{
	int nlats_src, nlons_src;
	double *max_lat_each_row_src, *min_lat_each_row_src;
	int i;
	double min_lat_diff, second_lat_diff, cur_lat_diff;
	double candidate_dist;

	nlats_src = grid_src->get_num_lats();
	nlons_src = grid_src->get_num_lons();
	max_lat_each_row_src = grid_src->get_max_lat_each_row();
	min_lat_each_row_src = grid_src->get_min_lat_each_row();

	min_lat_diff = second_lat_diff = 10000;

	for (i = 0; i < nlats_src; i ++) {
		cur_lat_diff = fabs(max_lat_each_row_src[i] - lat_dst);
		if (cur_lat_diff < min_lat_diff) {
			second_lat_diff = min_lat_diff;
			min_lat_diff = cur_lat_diff;
		}
		else if (cur_lat_diff < second_lat_diff) 
			second_lat_diff = cur_lat_diff;
		cur_lat_diff = fabs(min_lat_each_row_src[i] - lat_dst);
		if (cur_lat_diff < min_lat_diff) {
			second_lat_diff = min_lat_diff;
			min_lat_diff = cur_lat_diff;
		}
		else if (cur_lat_diff < second_lat_diff) 
			second_lat_diff = cur_lat_diff;
	}

	candidate_dist = cos(lat_dst * PI / 180) * num_nearest_pnts * PI / nlons_src;
	if (candidate_dist < second_lat_diff * PI / 180)
		candidate_dist = second_lat_diff * PI / 180;
	return candidate_dist * 6378.137;
}
                                   

void get_nearest_nbrs(double lat_dst,
                               double lon_dst,
                               grid *grid_src,
                               double dist_threshold,
                               int num_nearest_pnts,
                               int *point_indx_src,
                               double *weights,
                               double power_p)
{
	bool flag;
	double caled_thresh;
	double increment;
	static double last_thresh;
	static double last_lat = 10000;

	if (last_lat == lat_dst)
		caled_thresh = last_thresh;
	else caled_thresh = cal_dist_threshold(lat_dst, grid_src, num_nearest_pnts);

	increment = 1.2;

	do {
		flag = search_nearest_pnts(lat_dst,
		                                    lon_dst,
		                                    grid_src,
		                                    caled_thresh,
		                                    num_nearest_pnts, 
		                                    point_indx_src,
		                                    weights, 
		                                    power_p);
		caled_thresh *= increment;
	}
	while (!flag);
	last_thresh = caled_thresh/increment;
	last_lat = lat_dst;
}


