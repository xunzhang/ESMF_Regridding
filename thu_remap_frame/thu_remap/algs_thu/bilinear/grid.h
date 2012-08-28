#ifndef GRID_H
#define GRID_H

#include "name_cfg.h"

class grid
{
	private:
		//public:
		char my_name[NAME_STR_SIZE];
		unsigned int ndim;							// The number of dimensions of the grid
		int npnts;								// The number of points (cells) in the grid
		int nlats;								// The number of different latitudes
		int nlons;								// The number of different longitudes
		int nhghts;								// The number of different heights, 3D grid only
		int nvertexes;							// The maximum number of vertexes of each cell
		double *center_lat_coords;				// The latitude coordination of the center point of each cell
		double *center_lon_coords;				// The longitude coordination of the center point of each cell
		double *center_heights;					// The height coordination of the center point of each cell
		double *vertex_lat_coords;				// The latitude coordination of the vertexes of each cell
		double *vertex_lon_coords;				// The longitude coordination of the vertexes of each cell
		double *vertex_heights;					// The longitude coordination of the vertexes of each cell, 3D grid only
		double *area_in;							// The area of each cell read from file
		double *areas;							// The area of each cell after caculation
		bool *mask;								// The mask of each cell for each model
		double **model_frac;						// The area fraction of each cell for each model
		int *nlon_each_lat;						// The number of different longitudes in each latitude
		int *nlat_each_lon;						// The number of different latitudes in each longitude
		int *reverse_indx;						// The indexes for reversing grid from lat-major to lon-major
		double *different_lats;					// The different latitudes of the center points
		double *max_lat_each_row;				// The maximum latitude in each row of center latitudes
		double *min_lat_each_row;				// The minimum latitude in each row of center latitudes
		double *different_lons;					// The different longitudes of the center points
		bool is_rect;								// Whether the grid is regular 
		bool type_grid;							// The type of grid
		int my_model_set_id;						// The model set id of the grid corresponding to

		// adding
		double *bound_box;                 // lat/lon bounding box for use in restricting grid searches
		int    *dims;                      // 
		bool   is_gcenter;                 // Whether use centers for bounding boxes
		int    type_bin;                   // type of bins to use
		int    nsbin;                      // The number of bins for restricted search
		int    *bin_addr;                  // min,max adds for grid cells in this lat bin
		double *bin_lats;                  // min,max latitude for each search bin
		double *bin_lons;                  // min,max longitude for each search bin

	public:
		grid(unsigned int, int, int, int, double *, double *, bool *);
		~grid();
		char *get_grid_name() {return my_name; }
		int get_num_points() {return npnts; }
		int get_num_lons() {return nlons; }
		int get_num_lats() {return nlats; }
		double *get_lat_coords() {return center_lat_coords; }
		double *get_lon_coords() {return center_lon_coords; }
		double *get_cell_area() {return areas; }
		bool is_rect_grid() {return is_rect; }
		int *get_nlon_each_lat() {return nlon_each_lat; }
		int *get_nlat_each_lon() {return nlat_each_lon; }
		int *get_reverse_indx() {return reverse_indx; }
		double *get_diff_lats() { return different_lats; }
		double *get_diff_lons() {return different_lons; }
		bool *get_mymask() {return mask;}
		void cal_cell_area();
		int get_type_grid() {return 1; }
		double* get_max_lat_each_row() {return max_lat_each_row; }
		double* get_min_lat_each_row() {return min_lat_each_row; }


		//adding
		int get_num_vertexes() {return nvertexes; }
		int get_num_srch_bins() {return nsbin; }

		double *get_vertex_lat_coords() {return vertex_lat_coords; }
		double *get_vertex_lon_coords() {return vertex_lon_coords; }
		double *get_bound_box() {return bound_box; }
		int    *get_bin_addr() {return bin_addr; }
		double *get_single_frac(int n) {return model_frac[n]; }
		//bool   *get_single_mask(int n) {return mask[n]; }

		void test_bug();

#ifdef DEBUG_GRID
		void check_grid_info(int);
#endif
};

#endif
