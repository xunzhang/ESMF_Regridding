#ifndef REMAP_BASE_H
#define REMAP_BASE_H

#include "grid.h"

#define PI ((double) 3.141592653589793)
#define RADIAN(d) (d*(PI/180.0))

class Common_remap
{
	protected:
		/* Variables for search a remap algorithm */
		grid *grid_src;
		grid *grid_dst;
		char grid1_name[NAME_STR_SIZE];
		char grid2_name[NAME_STR_SIZE];
		char algorithm_name[NAME_STR_SIZE];
		char file_to_cal_derivative[NAME_STR_SIZE];
		//Common_remap *next;

		/* Localized fields of the source grid */
		int npts_src;
		int num_vertexes_src;
		int num_lats_src;	
		int num_lons_src;
		double *center_lons_src;	
		double *center_lats_src;		
		double *vertex_lons_src;
		double *vertex_lats_src;
		double *frac_src;
		double *area_src;
		bool *mask_src;

		/* Localized fields of the destination grid */		
		int npts_dst;
		int num_vertexes_dst;
		int num_lats_dst;
		int num_lons_dst;
		double *center_lons_dst;
		double *center_lats_dst;
		double *vertex_lons_dst;
		double *vertex_lats_dst;
		double *frac_dst;
		double *area_dst;
		bool *mask_dst;

		/* Variables for doing remap functions */
		int num_wgts;
		int *test;
		int *pine;
	public:
		int *wgt_indx_src;
		int *wgt_indx_dst;
		double *wgt_values;
		
		
	public:
		int *get_npts_src() { return &npts_src; }
		int *get_npts_dst() { return &npts_dst; }
		int *get_num_wgts() { return &num_wgts; }
		//int *get_wgt_indx_src() { return wgt_indx_src; }
		//int *get_wgt_indx_dst() { return wgt_indx_dst; }
		//double *get_wgt_values() { return wgt_values; }
		void localize_grid_fields();
		Common_remap(grid*, grid*);
		virtual ~Common_remap(){;}
		virtual void init_remap(){}
		virtual void cal_remap(double *, double *){}
		//Common_remap *search_remap_operator(char *, char *, char *);
};


#endif
