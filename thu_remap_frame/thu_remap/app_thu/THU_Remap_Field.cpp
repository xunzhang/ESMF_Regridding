//------------------------------------------------------------------------------------
// $Id: THU_Remap.cxx,v 1.00 2011/12 -- 2012/01
//-------------------------------------------------------------------------------------
#include "netcdf.h"
#include <string.h>
#include <iostream>
#include <math.h>
using namespace std;

void nc_err(int rcode)
{
	if(rcode != 0)
	{
		fprintf(stderr, "Netcdf IO error\n");
		exit(1);
	}
}

int main(int argc, char *argv[])
{
	// arguments
	char    srcfldfile[256], dstfldfile[256], funcname[256], flag[40];
	int     i,j,k;

	// nc return code and ncid
	int     rcode,fid,did,vid;

	// nc var id
	int     nc_grdsize_id,nc_grdcntrlat_id,nc_grdcntrlon_id,nc_physical_id;
	int     nc_nlons_src_id, nc_nlats_src_id, nc_nlons_dst_id, nc_nlats_dst_id;
	int     nc_nlonlat_size_id, nc_nlons_id, nc_nlats_id;
	int     nlons_dst_field, nlats_dst_field;
	// src&dst cell num, weight num, nlons, nlats
	size_t    npts_dst,npts_src,num_wgts;
	size_t    nlons_src, nlons_dst, nlats_src, nlats_dst;

	// weight info
	int	   *wgt_indx_src;
	int    *dst_mask;
	int    *src_mask;
	double *wgt_values;

	// field info
	double *data_dst;
	double *data_src;
	double *true_dst;
	double *diff_dst;
	double *src_lat;
	double *src_lon;
	double *dst_lat;
	double *dst_lon;

	// rmsd
	double diff_square_sum = 0;

	// transfer arguments
	if(argc != 3 && argc != 4 ){
		fprintf(stderr, "Arguments error! [ THU_Remap_Field -function functionname ]\n");
		fprintf(stderr, "                 [ THU_Remap_Field -field srcfieldfile dstfieldfile]\n");
		exit(1);
	}
	strcpy(flag, argv[1]);
	if(strcmp(flag,"-function")==0)
	{
		if(argc != 3) 
		{
			fprintf(stderr, "Arguments error! [ THU_Remap_Field -function functionname ]\n");
			exit(1);
		}
		strcpy(funcname, argv[2]);
		if(strcmp(funcname,"Y22")!=0 && strcmp(funcname,"Y16")!=0)
		{
			fprintf(stderr, "Function name error! [ Y22 or Y16 ]\n");
			exit(1);	
		}
	}
	else if(strcmp(flag,"-field")==0)
	{
		if(argc != 4)
		{
			fprintf(stderr, "Arguments error! [ THU_Remap_Field -field srcfieldfile dstfieldfile ]\n");
			exit(1);
		}
		strcpy(srcfldfile, argv[2]);
		strcpy(dstfldfile, argv[3]);
	}
	else
	{
		fprintf(stderr, "Argument error! [ -function ARG or -field ARG1 ARG2]\n");
		exit(1);
	}


	//  read weights.nc
	//---------------------------------------------
	fprintf(stdout, "read file \"weights.nc\" .........\n");
	// open ncfile
	fprintf(stdout, "open file...\n");
	rcode = nc_open("weights.nc", NC_NOWRITE, &fid);
	nc_err(rcode);

	// inquire dimlen
	fprintf(stdout, "inquire dimension...\n");
	rcode = nc_inq_dimid(fid, "n_a", &did);
	nc_err(rcode);
	rcode = nc_inq_dimlen(fid, did, &npts_src);
	nc_err(rcode);
	rcode = nc_inq_dimid(fid, "n_b", &did);
	nc_err(rcode);
	rcode = nc_inq_dimlen(fid, did, &npts_dst);
	nc_err(rcode);

	// allocate 
	data_dst = new double[npts_dst];
	data_src = new double[npts_src];
	true_dst = new double[npts_dst];
	diff_dst = new double[npts_dst];
	src_lat  = new double[npts_src];
	src_lon  = new double[npts_src];
	dst_lat  = new double[npts_dst];
	dst_lon  = new double[npts_dst];
	dst_mask = new int[npts_dst];
	src_mask = new int[npts_src];
	num_wgts = npts_dst * 4;
	wgt_indx_src = new int[num_wgts];
	wgt_values   = new double[num_wgts];

	//  read wgt_indx_src
	fprintf(stdout, "get variables...\n");

	rcode = nc_inq_varid(fid, "col", &vid);
	nc_err(rcode);
	rcode = nc_get_var_int(fid, vid, wgt_indx_src);
	nc_err(rcode);

	// read wgt_values
	rcode = nc_inq_varid(fid, "S", &vid);
	nc_err(rcode);
	rcode = nc_get_var_double(fid, vid, wgt_values);
	nc_err(rcode);

	// read grid_center_lat
	rcode = nc_inq_varid(fid, "yc_a", &vid);
	nc_err(rcode);
	rcode = nc_get_var_double(fid, vid, src_lat);
	nc_err(rcode);
	rcode = nc_inq_varid(fid, "yc_b", &vid);
	nc_err(rcode);
	rcode = nc_get_var_double(fid, vid, dst_lat);
	nc_err(rcode);

	// read grid_center_lon
	rcode = nc_inq_varid(fid, "xc_a", &vid);
	nc_err(rcode);
	rcode = nc_get_var_double(fid, vid, src_lon);
	nc_err(rcode);
	rcode = nc_inq_varid(fid, "xc_b", &vid);
	nc_err(rcode);
	rcode = nc_get_var_double(fid, vid, dst_lon);
	nc_err(rcode);

	// read mask
	rcode = nc_inq_varid(fid, "mask_b", &vid);
	nc_err(rcode);
	rcode = nc_get_var_int(fid, vid, dst_mask);
	nc_err(rcode);
	rcode = nc_inq_varid(fid, "mask_a", &vid);
	nc_err(rcode);
	rcode = nc_get_var_int(fid, vid, src_mask);
	nc_err(rcode);

	// close ncfile
	fprintf(stdout, "close file...\n");
	rcode = nc_close(fid);
	nc_err(rcode);

	fprintf(stdout, "read weights.nc ok..........\n");

	//compute remap
	//------------------------------------------------------

	// field file method
	if (strcmp(flag,"-field")==0)
	{	
		//read srcfield 
		fprintf(stdout, "read source field file.........\n");
		rcode = nc_open(srcfldfile, NC_NOWRITE, &fid);
		nc_err(rcode);

		rcode = nc_inq_varid(fid, "physical_variable", &nc_physical_id);
		nc_err(rcode);

		rcode = nc_get_var_double(fid, nc_physical_id, data_src);
		nc_err(rcode);

		rcode = nc_close(fid);
		nc_err(rcode);

		//read dstfield
		fprintf(stdout, "read destination field file.........\n");
		rcode = nc_open(dstfldfile, NC_NOWRITE, &fid);
		nc_err(rcode);

		rcode = nc_inq_varid(fid, "physical_variable", &nc_physical_id);
		nc_err(rcode);

		rcode = nc_get_var_double(fid, nc_physical_id, true_dst);
		nc_err(rcode);

		rcode = nc_close(fid);
		nc_err(rcode);

		// SMM & get difference
		fprintf(stdout, "calculate remapping.........\n");
		j = 0;
		diff_square_sum = 0;
		for(i=0; i<npts_dst; i++)
		{
			if(dst_mask[i] == 0) 
			{
				data_dst[i] = 0;
				true_dst[i] = 0;
			}
			else
			{
				data_dst[i] =  data_src[wgt_indx_src[j]]*wgt_values[j]+
					data_src[wgt_indx_src[j+1]]*wgt_values[j+1]+ 
					data_src[wgt_indx_src[j+2]]*wgt_values[j+2]+ 
					data_src[wgt_indx_src[j+3]]*wgt_values[j+3];
			}
			// get difference
			diff_dst[i] = data_dst[i]-true_dst[i];
//			if(fabs(diff_dst[i]) > 1.0e-16 && fabs(true_dst[i]) < 1.0e-16) diff_dst[i] = 1;
//			else if(fabs(diff_dst[i]) > 1.0e-16) diff_dst[i] /= true_dst[i];		

    		if (fabs(true_dst[i]) < 1.0e-16)
			{
				if (fabs(diff_dst[i]) < 1.0e-16)
					diff_dst[i] = 0;
				else
					diff_dst[i] = 1;
			}
			else
				diff_dst[i] /= true_dst[i];	

			diff_square_sum += diff_dst[i] * diff_dst[i];

			j = j + 4;
		}

		//output remapped field
		rcode = nc_create ("dstfield_file_remapped.nc", NC_WRITE, &fid);
		nc_err(rcode);
		//define dimension
		rcode = nc_def_dim (fid, "grid_size", npts_dst, &nc_grdsize_id);
		nc_err(rcode);
		rcode = nc_def_dim (fid, "nlonlat_size", 1, &nc_nlonlat_size_id);
		nc_err(rcode);
		//define variable
		rcode = nc_def_var (fid, "grid_center_lat", NC_DOUBLE, 1, &nc_grdsize_id, &nc_grdcntrlat_id);
		nc_err(rcode);
		rcode = nc_def_var (fid, "grid_center_lon", NC_DOUBLE, 1, &nc_grdsize_id, &nc_grdcntrlon_id);
		nc_err(rcode);
		rcode = nc_def_var (fid, "physical_variable", NC_DOUBLE, 1, &nc_grdsize_id, &nc_physical_id);
		nc_err(rcode);
		rcode = nc_def_var (fid, "nlons", NC_INT, 1, &nc_nlonlat_size_id, &nc_nlons_id);
		nc_err(rcode);
		rcode = nc_def_var (fid, "nlats", NC_INT, 1, &nc_nlonlat_size_id, &nc_nlats_id);

		rcode = nc_enddef(fid);
		nc_err(rcode);
		//put variable
		rcode = nc_put_var_double (fid, nc_grdcntrlat_id, dst_lat);
		nc_err(rcode);
		rcode = nc_put_var_double (fid, nc_grdcntrlon_id, dst_lon);
		nc_err(rcode);
		rcode = nc_put_var_double (fid, nc_physical_id, data_dst);
		nc_err(rcode);
		rcode = nc_put_var_int (fid, nc_nlons_id, &nlons_dst_field);
		nc_err(rcode);
		rcode = nc_put_var_int (fid, nc_nlats_id, &nlats_dst_field);		
		nc_err(rcode);

		rcode = nc_close (fid);
		nc_err(rcode);

		fprintf(stdout, "output dstfield_func_remapped.nc\n");

		//output "true" dst field
		rcode = nc_create ("dstfield_file_input.nc", NC_WRITE, &fid);
		nc_err(rcode);
		//define dimension
		rcode = nc_def_dim (fid, "grid_size", npts_dst, &nc_grdsize_id);
		nc_err(rcode);
		rcode = nc_def_dim (fid, "nlonlat_size", 1, &nc_nlonlat_size_id);
		nc_err(rcode);
		//define variable
		rcode = nc_def_var (fid, "grid_center_lat", NC_DOUBLE, 1, &nc_grdsize_id, &nc_grdcntrlat_id);
		nc_err(rcode);
		rcode = nc_def_var (fid, "grid_center_lon", NC_DOUBLE, 1, &nc_grdsize_id, &nc_grdcntrlon_id);
		nc_err(rcode);
		rcode = nc_def_var (fid, "physical_variable", NC_DOUBLE, 1, &nc_grdsize_id, &nc_physical_id);
		nc_err(rcode);
		rcode = nc_def_var (fid, "nlons", NC_INT, 1, &nc_nlonlat_size_id, &nc_nlons_id);
		nc_err(rcode);
		rcode = nc_def_var (fid, "nlats", NC_INT, 1, &nc_nlonlat_size_id, &nc_nlats_id);

		rcode = nc_enddef(fid);
		nc_err(rcode);
		//put variable
		rcode = nc_put_var_double (fid, nc_grdcntrlat_id, dst_lat);
		nc_err(rcode);
		rcode = nc_put_var_double (fid, nc_grdcntrlon_id, dst_lon);
		nc_err(rcode);
		rcode = nc_put_var_double (fid, nc_physical_id, true_dst);
		nc_err(rcode);
		rcode = nc_put_var_int (fid, nc_nlons_id, &nlons_dst_field);
		nc_err(rcode);
		rcode = nc_put_var_int (fid, nc_nlats_id, &nlats_dst_field);		
		nc_err(rcode);

		rcode = nc_close (fid);
		nc_err(rcode);

		fprintf(stdout, "output dstfield_func_input.nc\n");



		// output difference field
		rcode = nc_create ("dstfield_file_difference.nc", NC_WRITE, &fid);
		nc_err(rcode);
		//define dimension
		rcode = nc_def_dim (fid, "grid_size", npts_dst, &nc_grdsize_id);
		nc_err(rcode);
		rcode = nc_def_dim (fid, "nlonlat_size", 1, &nc_nlonlat_size_id);
		nc_err(rcode);
		//define variable
		rcode = nc_def_var (fid, "grid_center_lat", NC_DOUBLE, 1, &nc_grdsize_id, &nc_grdcntrlat_id);
		nc_err(rcode);
		rcode = nc_def_var (fid, "grid_center_lon", NC_DOUBLE, 1, &nc_grdsize_id, &nc_grdcntrlon_id);
		nc_err(rcode);
		rcode = nc_def_var (fid, "physical_variable", NC_DOUBLE, 1, &nc_grdsize_id, &nc_physical_id);
		nc_err(rcode);
		rcode = nc_def_var (fid, "nlons", NC_INT, 1, &nc_nlonlat_size_id, &nc_nlons_id);
		nc_err(rcode);
		rcode = nc_def_var (fid, "nlats", NC_INT, 1, &nc_nlonlat_size_id, &nc_nlats_id);

		rcode = nc_enddef(fid);
		nc_err(rcode);
		//put variable
		rcode = nc_put_var_double (fid, nc_grdcntrlat_id, dst_lat);
		nc_err(rcode);
		rcode = nc_put_var_double (fid, nc_grdcntrlon_id, dst_lon);
		nc_err(rcode);
		rcode = nc_put_var_double (fid, nc_physical_id, diff_dst);
		nc_err(rcode);
		rcode = nc_put_var_int (fid, nc_nlons_id, &nlons_dst_field);
		nc_err(rcode);
		rcode = nc_put_var_int (fid, nc_nlats_id, &nlats_dst_field);		
		nc_err(rcode);

		rcode = nc_close (fid);
		nc_err(rcode);

		fprintf(stdout, "output dstfield_file_difference.nc\n");
		fprintf(stdout, "\nRoot mean square deviation: %lf\n", sqrt(diff_square_sum/npts_dst));
	}	
	// function method
	else 
	{
		//generate source field and true dst field

		// deg to rad 
		double deg2rad = 3.14159265359 / 180;
		if (strcmp(funcname, "Y22") == 0)
		{
			fprintf(stdout, "generate source field data with funtion 'Y22'\n");
			k = 0;
			for(i=0; i<npts_src; i++)
			{
				if(src_mask[k] == 0)
					data_src[k] = 0;
				else
					data_src[k] = 2+cos(src_lat[i]*deg2rad)*cos(src_lat[i]*deg2rad)*cos(2*src_lon[i]*deg2rad);
				k++;
			}
			k = 0;
			for(i=0; i<npts_dst; i++)
			{	
				if(dst_mask[k] == 0)
					true_dst[k] = 0;
				else
					true_dst[k] = 2+cos(dst_lat[i]*deg2rad)*cos(dst_lat[i]*deg2rad)*cos(2*dst_lon[i]*deg2rad);
				k++;
			}
		}
		else if (strcmp(funcname, "Y16") == 0)
		{
			fprintf(stdout, "generate source field data with funtion 'Y16'\n");
			k = 0;
			for(i=0; i<npts_src; i++)
			{
				if(src_mask[i] == 0)
					data_src[k] = 0;
				else
					data_src[k] = 2+pow(cos(src_lat[i]*deg2rad),16)*cos(2*src_lon[i]*deg2rad);
				k++;
			}
			k = 0;
			for(i=0; i<npts_dst; i++)
			{
				if(dst_mask[k] == 0)
					true_dst[k] = 0;
				else
					true_dst[k] = 2+pow(cos(dst_lat[i]*deg2rad),16)*cos(2*dst_lon[i]*deg2rad);
				k++;
			}
		}

		// SMM & get diff
		fprintf(stdout, "calculate remapping.........\n");
		j = 0;
		diff_square_sum = 0;
		for(i=0; i<npts_dst; i++)
		{
			if(dst_mask[i] == 0) 
			{
				data_dst[i] = 0;
			}
			else
			{
				data_dst[i] =  data_src[wgt_indx_src[j]]*wgt_values[j]+
					data_src[wgt_indx_src[j+1]]*wgt_values[j+1]+ 
					data_src[wgt_indx_src[j+2]]*wgt_values[j+2]+ 
					data_src[wgt_indx_src[j+3]]*wgt_values[j+3];
			}
			diff_dst[i] = data_dst[i]-true_dst[i];
			//if(fabs(diff_dst[i]) > 1.0e-16 && fabs(true_dst[i]) < 1.0e-16) diff_dst[i] = 1;
			//else if(fabs(diff_dst[i]) > 1.0e-16 && fabs(truedst[i]) > 1.0e-16) diff_dst[i] /= true_dst[i];

    		if (fabs(true_dst[i]) < 1.0e-16)
			{
				if (fabs(diff_dst[i]) < 1.0e-16)
					diff_dst[i] = 0;
				else
					diff_dst[i] = 1;
			}
			else
				diff_dst[i] /= true_dst[i];	

			diff_square_sum += diff_dst[i] * diff_dst[i];

			j = j + 4;
		}

		// output dst field, true field, difference field
		//-----------------------------------------------

		//output source field

		rcode = nc_create ("dstfield_func_remapped.nc", NC_WRITE, &fid);
		nc_err(rcode);
		// define dimension
		rcode = nc_def_dim (fid, "grid_size", npts_dst, &nc_grdsize_id);
		nc_err(rcode);
		rcode = nc_def_dim (fid, "nlonlat_size", 1, &nc_nlonlat_size_id);
		nc_err(rcode);
		// define variable
		rcode = nc_def_var (fid, "grid_center_lat", NC_DOUBLE, 1, &nc_grdsize_id, &nc_grdcntrlat_id);
		nc_err(rcode);
		rcode = nc_def_var (fid, "grid_center_lon", NC_DOUBLE, 1, &nc_grdsize_id, &nc_grdcntrlon_id);
		nc_err(rcode);
		rcode = nc_def_var (fid, "physical_variable", NC_DOUBLE, 1, &nc_grdsize_id, &nc_physical_id);
		nc_err(rcode);
		rcode = nc_def_var (fid, "nlons", NC_INT, 1, &nc_nlonlat_size_id, &nc_nlons_id);
		nc_err(rcode);
		rcode = nc_def_var (fid, "nlats", NC_INT, 1, &nc_nlonlat_size_id, &nc_nlats_id);

		rcode = nc_enddef(fid);
		nc_err(rcode);
		//put variables
		rcode = nc_put_var_double (fid, nc_grdcntrlat_id, dst_lat);
		nc_err(rcode);
		rcode = nc_put_var_double (fid, nc_grdcntrlon_id, dst_lon);
		nc_err(rcode);
		rcode = nc_put_var_double (fid, nc_physical_id, data_dst);
		nc_err(rcode);
		rcode = nc_put_var_int (fid, nc_nlons_id, &nlons_dst_field);
		nc_err(rcode);
		rcode = nc_put_var_int (fid, nc_nlats_id, &nlats_dst_field);		
		nc_err(rcode);

		rcode = nc_close (fid);
		nc_err(rcode);

		fprintf(stdout, "output dstfield_func_remapped.nc\n");

		// output true field
		rcode = nc_create ("dstfield_func_generated.nc", NC_WRITE, &fid);
		nc_err(rcode);
		//define dimension
		rcode = nc_def_dim (fid, "grid_size", npts_dst, &nc_grdsize_id);
		nc_err(rcode);
		rcode = nc_def_dim (fid, "nlonlat_size", 1, &nc_nlonlat_size_id);
		nc_err(rcode);
		//define variable
		rcode = nc_def_var (fid, "grid_center_lat", NC_DOUBLE, 1, &nc_grdsize_id, &nc_grdcntrlat_id);
		nc_err(rcode);
		rcode = nc_def_var (fid, "grid_center_lon", NC_DOUBLE, 1, &nc_grdsize_id, &nc_grdcntrlon_id);
		nc_err(rcode);
		rcode = nc_def_var (fid, "physical_variable", NC_DOUBLE, 1, &nc_grdsize_id, &nc_physical_id);
		nc_err(rcode);
		rcode = nc_def_var (fid, "nlons", NC_INT, 1, &nc_nlonlat_size_id, &nc_nlons_id);
		nc_err(rcode);
		rcode = nc_def_var (fid, "nlats", NC_INT, 1, &nc_nlonlat_size_id, &nc_nlats_id);

		rcode = nc_enddef(fid);
		nc_err(rcode);
		// put variable
		rcode = nc_put_var_double (fid, nc_grdcntrlat_id, dst_lat);
		nc_err(rcode);
		rcode = nc_put_var_double (fid, nc_grdcntrlon_id, dst_lon);
		nc_err(rcode);
		rcode = nc_put_var_double (fid, nc_physical_id, true_dst);
		nc_err(rcode);
		rcode = nc_put_var_int (fid, nc_nlons_id, &nlons_dst_field);
		nc_err(rcode);
		rcode = nc_put_var_int (fid, nc_nlats_id, &nlats_dst_field);		
		nc_err(rcode);

		rcode = nc_close (fid);
		nc_err(rcode);

		fprintf(stdout, "output dstfield_func_generated.nc\n");

		// output difference field
		rcode = nc_create ("dstfield_func_difference.nc", NC_WRITE, &fid);
		nc_err(rcode);
		// define dimension
		rcode = nc_def_dim (fid, "grid_size", npts_dst, &nc_grdsize_id);
		nc_err(rcode);
		rcode = nc_def_dim (fid, "nlonlat_size", 1, &nc_nlonlat_size_id);
		nc_err(rcode);
		//define variable
		rcode = nc_def_var (fid, "grid_center_lat", NC_DOUBLE, 1, &nc_grdsize_id, &nc_grdcntrlat_id);
		nc_err(rcode);
		rcode = nc_def_var (fid, "grid_center_lon", NC_DOUBLE, 1, &nc_grdsize_id, &nc_grdcntrlon_id);
		nc_err(rcode);
		rcode = nc_def_var (fid, "physical_variable", NC_DOUBLE, 1, &nc_grdsize_id, &nc_physical_id);
		nc_err(rcode);
		rcode = nc_def_var (fid, "nlons", NC_INT, 1, &nc_nlonlat_size_id, &nc_nlons_id);
		nc_err(rcode);
		rcode = nc_def_var (fid, "nlats", NC_INT, 1, &nc_nlonlat_size_id, &nc_nlats_id);

		rcode = nc_enddef(fid);
		nc_err(rcode);
		//put variable
		rcode = nc_put_var_double (fid, nc_grdcntrlat_id, dst_lat);
		nc_err(rcode);
		rcode = nc_put_var_double (fid, nc_grdcntrlon_id, dst_lon);
		nc_err(rcode);
		rcode = nc_put_var_double (fid, nc_physical_id, diff_dst);
		nc_err(rcode);
		rcode = nc_put_var_int (fid, nc_nlons_id, &nlons_dst_field);
		nc_err(rcode);
		rcode = nc_put_var_int (fid, nc_nlats_id, &nlats_dst_field);		
		nc_err(rcode);

		rcode = nc_close (fid);
		nc_err(rcode);

		fprintf(stdout, "output dstfield_func_difference.nc\n");
		fprintf(stdout, "\nRoot mean square deviation: %lf\n", sqrt(diff_square_sum/npts_dst));
	}

	//  delete
	delete []data_dst;
	delete []data_src;
	delete []true_dst;
	delete []diff_dst;
	delete []src_lat;
	delete []src_lon;
	delete []dst_lat;
	delete []dst_lon;
	delete []wgt_indx_src;
	delete []wgt_values;
	delete []dst_mask;
	delete []src_mask;

	fprintf(stdout, "Completed remap successfully\n");

	return 0;
}
