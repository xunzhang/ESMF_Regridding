#include "bilinear_remap.h"
//#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "grid.h"
#include <netcdf.h>
#include "t_get_nbr.h"

#define ERRCODE 2

#define ERR(e) {printf("Error: %s/n", nc_strerror(e)); exit(ERRCODE);}

bilinear_remap_2D_operator::bilinear_remap_2D_operator(grid *out_grid_src, grid *out_grid_dst):Common_remap(out_grid_src, out_grid_dst)
{
	int i,j;

	grid_src_local = out_grid_src; 
	grid_dst_local = out_grid_dst;

	nlon_each_lat_src=grid_src_local->get_nlon_each_lat();
	nlon_each_lat_dst=grid_dst_local->get_nlon_each_lat();

	pine = new int [1];

	lat_begindx_src = new int [num_lats_src];
	lat_begindx_dst = new int [num_lats_dst];
	num_wgts = npts_dst * 4;
	wgt_indx_src = NULL;
	//wgt_indx_src = new int [num_wgts];
	//printf(" is is %x\n", wgt_indx_src);
	//delete wgt_indx_src;
	wgt_indx_src = new int [num_wgts];
	wgt_indx_dst = NULL;
	wgt_values = new double [num_wgts];

	lat_begindx_src[0] = 0;
	for(i = 1; i < num_lats_src; i ++) 
		lat_begindx_src[i] = lat_begindx_src[i - 1] + nlon_each_lat_src[i - 1];

	lat_begindx_dst[0] = 0;
	for(i = 1; i < num_lats_dst; i ++)
		lat_begindx_dst[i] = lat_begindx_dst[i - 1] + nlon_each_lat_dst[i - 1];
}


bilinear_remap_2D_operator::~bilinear_remap_2D_operator()
{
	delete [] lat_begindx_src;
	delete [] lat_begindx_dst;
	delete [] wgt_indx_src;
	//delete [] wgt_indx_dst;
	delete [] wgt_values;
}


void bilinear_remap_2D_operator::wgt_bilinear_regular( double *lons_src,
		double *lats_src,
		double lon_dst,
		double lat_dst,
		double *wgts) 
{
	double t1,t2;
	double u;

	if (!((lon_dst >= lons_src[0]) &&
				(lons_src[1] >= lons_src[0]))) {
		if (lon_dst < lons_src[0])
			lon_dst += 360;
		if (lons_src[1] < lons_src[0])
			lons_src[1] += 360;
	}
	if (!((lon_dst >= lons_src[2]) &&
				(lons_src[3] >= lons_src[2]))) {
		if (lon_dst < lons_src[2])
			lon_dst += 360;
		if (lons_src[3] < lons_src[2])
			lons_src[3] += 360;
	}
	t1=(lon_dst-lons_src[0])/(lons_src[1]-lons_src[0]);
	t2=(lon_dst-lons_src[2])/(lons_src[3]-lons_src[2]);
	u=(lat_dst-lats_src[0])/(lats_src[1]-lats_src[0]);
	wgts[0] = (1 - u) * (1 - t1);
	wgts[1] = (1 - u) * t1;
	wgts[2] = u * (1 - t2);
	wgts[3] = u * t2;
}			 



bool bilinear_remap_2D_operator::wgt_bilinear_irregular( double *lons_src,
		double *lats_src,
		double lon_dst,
		double lat_dst,
		double *wgts) 
{
	double dth1, dth2, dth3, dthp;
	double dph1, dph2, dph3, dphp;
	double iguess, jguess;
	double mat1, mat2, mat3, mat4;
	int iter, max_iter = 100;
	double determinant, deli, delj;
	double converge = 0.0000000001;

	dth1 = lats_src[1] - lats_src[0];
	dth2 = lats_src[3] - lats_src[0];
	dth3 = lats_src[2] - lats_src[1] - dth2;
	dph1 = lons_src[1] - lons_src[0];
	dph2 = lons_src[3] - lons_src[0];
	dph3 = lons_src[2] - lons_src[1];
	if (dph1 >  540) 
		dph1 = dph1 - 360;
	if (dph2 >  540) 
		dph2 = dph2 - 360;
	if (dph3 >  540) 
		dph3 = dph3 - 360;
	if (dph1 < -540) 
		dph1 = dph1 + 360;
	if (dph2 < -540) 
		dph2 = dph2 + 360;
	if (dph3 < -540) 
		dph3 = dph3 + 360;
	dph3 = dph3 - dph2;
	iguess = 0.5;
	jguess = 0.5;

	for (iter = 0; iter < max_iter; iter ++) {
		dthp = lat_dst - lats_src[0] - dth1*iguess - dth2*jguess - dth3*iguess*jguess;
		dphp = lon_dst - lons_src[0];
		if (dphp >  540) dphp = dphp - 360;
		if (dphp < -540) dphp = dphp + 360;
		dphp = dphp - dph1*iguess - dph2*jguess - dph3*iguess*jguess;
		mat1 = dth1 + dth3*jguess;
		mat2 = dth2 + dth3*iguess;
		mat3 = dph1 + dph3*jguess;
		mat4 = dph2 + dph3*iguess;

		determinant = mat1*mat4 - mat2*mat3;
		deli = (dthp*mat4 - mat2*dphp) / determinant;
		delj = (mat1*dphp - dthp*mat3) / determinant;
		if (fabs(deli) < converge && fabs(delj) < converge)
			break;

		iguess = iguess + deli;
		jguess = jguess + delj;
	}

	if (iter < max_iter) {
		wgts[0] = (1.0-iguess) * (1.0-jguess);
		wgts[1] = iguess * (1.0-jguess); 
		wgts[2] = iguess * jguess; 
		wgts[3] = (1.0-iguess) * jguess;
		return true;
	}
	else return false;
}			 


void bilinear_remap_2D_operator::init_remap()
{

	int i, j, k, point_dst = 0;
	int bound_lat_indx1, bound_lat_indx2;
	int bound_lon_indx1, bound_lon_indx2;
	int left_nbr, right_nbr;
	double left_lon, right_lon;
	double left_lat, right_lat;
	double bilinear_lons[4];
	double bilinear_lats[2];
	bool *mask_src, *mask_dst;
	double *max_lat_each_row_src, *min_lat_each_row_src;
	double last_lat = -10000;
	bool find_bilinear_point;
	int nlons_src;
	bool lat_in, lon_in;
	double lons_src_box[4], lats_src_box[4];
	int corner_indx[4];
	double max_lat, min_lat;
	double max_lon, min_lon, tmp_lon;
	double *bound_boxes_src;
	double *cornner_boxes_lat_src, *cornner_boxes_lon_src;

	mask_src = grid_src->get_mymask();
	mask_dst = grid_dst->get_mymask();
	max_lat_each_row_src = grid_src_local->get_max_lat_each_row();
	min_lat_each_row_src = grid_src_local->get_min_lat_each_row();
	nlons_src = grid_src_local->get_num_lons();

	//printf("test mask is %d\n", mask_dst[0]);
	for (point_dst = 0; point_dst < npts_dst; point_dst ++) {
		//printf("he yu 1 is %d\n", mask_src[point_dst]);
		if (last_lat != center_lats_dst[point_dst]) {
			last_lat = center_lats_dst[point_dst];
			bound_lat_indx1 = -1;
			for(i = 0; i < num_lats_src; i ++)
				if(last_lat - max_lat_each_row_src[i] < 0) {
					bound_lat_indx1 = i - 1;
					break;
				}
			bound_lat_indx2 = bound_lat_indx1 + 1;
		}

		if (!mask_dst[point_dst]) {
			wgt_indx_src[point_dst * 4 + 0] = 0;
			wgt_indx_src[point_dst * 4 + 1] = 0;
			wgt_indx_src[point_dst * 4 + 2] = 0;
			wgt_indx_src[point_dst * 4 + 3] = 0;
			wgt_values[point_dst * 4 + 0] = 0;
			wgt_values[point_dst * 4 + 1] = 0;
			wgt_values[point_dst * 4 + 2] = 0;
			wgt_values[point_dst * 4 + 3] = 0;
			continue;
		}

		if (bound_lat_indx1 == -1) {
			get_nearest_nbrs(center_lats_dst[point_dst], 
					center_lons_dst[point_dst],
					grid_src, 
					120.0, 
					4, 
					wgt_indx_src + point_dst * 4,
					wgt_values + point_dst * 4, 
					1.0);

			continue;
		}

		bound_lon_indx2 = 0;
		bound_lon_indx1 = 0;
		for(i = 0; i < nlons_src; i ++) {
			left_nbr = (i - 1 + nlons_src) % nlons_src;
			right_nbr = i + lat_begindx_src[bound_lat_indx1];
			left_nbr += lat_begindx_src[bound_lat_indx1];
			left_lon = center_lons_src[left_nbr];
			right_lon = center_lons_src[right_nbr];
			if (right_lon >= left_lon) {
				if (center_lons_dst[point_dst] - right_lon < 0 &&
						center_lons_dst[point_dst] - left_lon >= 0) {
					bound_lon_indx1 = i - 1;
					break;
				}
			}
			else if (center_lons_dst[point_dst] - right_lon < 0) {
				bound_lon_indx1 = i - 1;
				break;
			}
		}
		bound_lon_indx2 = i;
		if (i == 0 || i == nlons_src) {
			bound_lon_indx1 = nlons_src - 1;
			bound_lon_indx2 = 0;
		} 

		if (!mask_src[bound_lon_indx1 + lat_begindx_src[bound_lat_indx1]] ||
				!mask_src[bound_lon_indx2 + lat_begindx_src[bound_lat_indx1]] ||
				!mask_src[bound_lon_indx1 + lat_begindx_src[bound_lat_indx2]] ||
				!mask_src[bound_lon_indx2 + lat_begindx_src[bound_lat_indx2]]) {
			get_nearest_nbrs(center_lats_dst[point_dst], 
					center_lons_dst[point_dst],
					grid_src, 
					120.0, 
					4, 
					wgt_indx_src + point_dst * 4,
					wgt_values + point_dst * 4, 
					1.0);
		}
		else {

			wgt_indx_src[point_dst * 4 + 0] = bound_lon_indx1 + lat_begindx_src[bound_lat_indx1];
			wgt_indx_src[point_dst * 4 + 1] = bound_lon_indx2 + lat_begindx_src[bound_lat_indx1];
			wgt_indx_src[point_dst * 4 + 2] = bound_lon_indx1 + lat_begindx_src[bound_lat_indx2];
			wgt_indx_src[point_dst * 4 + 3] = bound_lon_indx2 + lat_begindx_src[bound_lat_indx2];
			bilinear_lons[0]= center_lons_src[bound_lon_indx1 + lat_begindx_src[bound_lat_indx1]];
			bilinear_lons[1] = center_lons_src[bound_lon_indx2 + lat_begindx_src[bound_lat_indx1]];
			bilinear_lons[2] = center_lons_src[bound_lon_indx1 + lat_begindx_src[bound_lat_indx2]];
			bilinear_lons[3] = center_lons_src[bound_lon_indx2 + lat_begindx_src[bound_lat_indx2]];
			bilinear_lats[0] = center_lats_src[lat_begindx_src[bound_lat_indx1]]; 
			bilinear_lats[1] = center_lats_src[lat_begindx_src[bound_lat_indx2]];

			//for debug
			/*
			printf("bililiner_lons[0] is %lf \n", bilinear_lons[0]);
			printf("bililiner_lons[1] is %lf \n", bilinear_lons[1]);
			printf("bililiner_lons[2] is %lf \n", bilinear_lons[2]);
			printf("bililiner_lons[3] is %lf \n", bilinear_lons[3]);
			printf("bililiner_lats[0] is %lf \n", bilinear_lats[0]);
			printf("bililiner_lats[1] is %lf \n", bilinear_lats[1]);
			printf("center_lons_dst[%d] is %lf \n", point_dst, center_lons_dst[point_dst]);
			printf("center_lats_dst[%d] is %lf \n", point_dst, center_lats_dst[point_dst]);
			*/

			// end of debug
			wgt_bilinear_regular(bilinear_lons, bilinear_lats, center_lons_dst[point_dst], center_lats_dst[point_dst], wgt_values + 4*point_dst);
			for(int kk = 0; kk < 4; ++kk)
				printf("%f ", wgt_values[4*point_dst+kk]);
		}
	}
/*
	int retval, ncid;
	int nc_src_grid_size_id, nc_dst_grid_size_id;
	int nc_src_grid_corner_id, nc_dst_grid_corner_id;
	int nc_src_grid_rank_id, nc_dst_grid_rank_id;
	int nc_src_grid_dims_id, nc_dst_grid_dims_id;
	int nc_src_grid_corner_lat_id, nc_dst_grid_corner_lat_id;
	int nc_src_grid_corner_lon_id, nc_dst_grid_corner_lon_id;
	int nc_num_links_id, nc_num_wgts_id;
	int nc_remap_matrix_id;
	int nc_dst_mask_id;

	if ((retval = nc_create("weights.nc", NC_CLOBBER, &ncid)))
		ERR(retval);
	
printf("nininimamamam\n");	
	if((retval =  nc_def_dim(ncid, "src_grid_size", npts_src, &nc_src_grid_size_id)))
		ERR(retval);
	if((retval =  nc_def_dim(ncid, "dst_grid_size", npts_dst, &nc_dst_grid_size_id)))
		ERR(retval);
		
	if((retval =  nc_def_dim(ncid, "src_grid_corners", 4, &nc_src_grid_corner_id)))
		ERR(retval);
	if((retval =  nc_def_dim(ncid, "dst_grid_corners", 4, &nc_dst_grid_corner_id)))
		ERR(retval);
		
	if((retval =  nc_def_dim(ncid, "src_grid_rank", 2, &nc_src_grid_rank_id)))
		ERR(retval);
	if((retval =  nc_def_dim(ncid, "dst_grid_rank", 2, &nc_dst_grid_rank_id)))
		ERR(retval);	
	
	if((retval =  nc_def_dim(ncid, "num_links", 4 * npts_dst, &nc_num_links_id)))
		ERR(retval);
	
	if((retval =  nc_def_dim(ncid, "num_wgts", 1, &nc_num_wgts_id)))
		ERR(retval);
		
	if((retval =  nc_def_var(ncid, "src_grid_dims", NC_INT, 1, &nc_src_grid_rank_id, &nc_src_grid_dims_id)))
		ERR(retval);
	if((retval =  nc_def_var(ncid, "dst_grid_dims", NC_INT, 1, &nc_dst_grid_rank_id, &nc_dst_grid_dims_id)))
		ERR(retval);	
		
	if((retval =  nc_def_var(ncid, "src_grid_center_lat", NC_DOUBLE, 1, &nc_src_grid_size_id, &nc_src_grid_corner_lat_id)))
		ERR(retval);
	if((retval =  nc_def_var(ncid, "dst_grid_center_lat", NC_DOUBLE, 1, &nc_dst_grid_size_id, &nc_dst_grid_corner_lat_id)))
		ERR(retval);
	if((retval =  nc_def_var(ncid, "src_grid_center_lon", NC_DOUBLE, 1, &nc_src_grid_size_id, &nc_src_grid_corner_lon_id)))
		ERR(retval);
	if((retval =  nc_def_var(ncid, "dst_grid_center_lon", NC_DOUBLE, 1, &nc_dst_grid_size_id, &nc_dst_grid_corner_lon_id)))
		ERR(retval);
	if((retval = nc_def_var(ncid, "dst_mask", NC_INT, 1, &nc_dst_grid_size_id, &nc_dst_mask_id)))
	        ERR(retval);
	//if((retval =  nc_def_var(ncid, 'src_address', NC_INT, 1, nc_dst_grid_rank_id, nc_dst_grid_conrner_lon_id)))
	//	ERR(retval);
	int nc_dims2 = nc_num_links_id;
	int nc_src_add_id;

	if((retval = nc_def_var(ncid, "src_address", NC_INT, 1, &nc_num_links_id, &nc_src_add_id)))
		ERR(retval);
	if((retval =  nc_def_var(ncid, "remap_matrix", NC_DOUBLE, 1, &nc_dims2, &nc_remap_matrix_id)))
		ERR(retval);

	if((retval = nc_enddef(ncid)))
		ERR(retval);
		
	int itmp1, itmp2, itmp3, itmp4;
	itmp1 = nc_src_grid_corner_lat_id;
	itmp2 = nc_dst_grid_corner_lat_id;
	itmp3 = nc_src_grid_corner_lon_id;
	itmp4 = nc_dst_grid_corner_lon_id;
	if((retval = nc_put_var_double(ncid, itmp1, center_lats_src)))
		ERR(retval);
	if((retval = nc_put_var_double(ncid, itmp2, center_lats_dst)))
		ERR(retval);
	if((retval = nc_put_var_double(ncid, itmp3, center_lons_src)))
		ERR(retval);
	if((retval = nc_put_var_double(ncid, itmp4, center_lons_dst)))
		ERR(retval);
	if((retval = nc_put_var_int(ncid, nc_src_add_id, wgt_indx_src)))
		ERR(retval);
	if((retval = nc_put_var_double(ncid, nc_remap_matrix_id, wgt_values)))
		ERR(retval);
	if((retval = nc_put_var_int(ncid, nc_dst_mask_id, (int *)mask_dst)))
	        ERR(retval);
	if((retval = nc_close(ncid)))
		ERR(retval);

	//for(i = 0; i < npts_dst;i++)
	//	printf("%dth is %lf\n", i, wgt_values[i]);
*/
}

