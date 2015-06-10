#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mfhdf.h>
#include "modisresam.h"

// define band names for all bands
const char *bandNames[38] = {
	"1",		/*  0. "EV_250_Aggr1km_RefSB" */
	"2",		/*  1. "EV_250_Aggr1km_RefSB" */
	"3",		/*  2. "EV_500_Aggr1km_RefSB" */
	"4",		/*  3. "EV_500_Aggr1km_RefSB" */
	"5",		/*  4. "EV_500_Aggr1km_RefSB" */
	"6",		/*  5. "EV_500_Aggr1km_RefSB" */
	"7",		/*  6. "EV_500_Aggr1km_RefSB" */
	"8",		/*  7. "EV_1KM_RefSB" */
	"9",		/*  8. "EV_1KM_RefSB" */
	"10",		/*  9. "EV_1KM_RefSB" */
	"11",		/* 10. "EV_1KM_RefSB" */
	"12",		/* 11. "EV_1KM_RefSB" */
	"13lo",		/* 12. "EV_1KM_RefSB" */
	"13hi",		/* 13. "EV_1KM_RefSB" */
	"14lo",		/* 14. "EV_1KM_RefSB" */
	"14hi",		/* 15. "EV_1KM_RefSB" */
	"15",		/* 16. "EV_1KM_RefSB" */
	"16",		/* 17. "EV_1KM_RefSB" */
	"17",		/* 18. "EV_1KM_RefSB" */
	"18",		/* 19. "EV_1KM_RefSB" */
	"19",		/* 20. "EV_1KM_RefSB" */
	"26",		/* 21. "EV_1KM_RefSB" (sic) */
	"20",		/* 22. "EV_1KM_Emissive" */
	"21",		/* 23. "EV_1KM_Emissive" */
	"22",		/* 24. "EV_1KM_Emissive" */
	"23",		/* 25. "EV_1KM_Emissive" */
	"24",		/* 26. "EV_1KM_Emissive" */
	"25",		/* 27. "EV_1KM_Emissive" */
	"27",		/* 28. "EV_1KM_Emissive" (sic) */
	"28",		/* 29. "EV_1KM_Emissive" */
	"29",		/* 30. "EV_1KM_Emissive" */
	"30",		/* 31. "EV_1KM_Emissive" */
	"31",		/* 32. "EV_1KM_Emissive" */
	"32",		/* 33. "EV_1KM_Emissive" */
	"33",		/* 34. "EV_1KM_Emissive" */
	"34",		/* 35. "EV_1KM_Emissive" */
	"35",		/* 36. "EV_1KM_Emissive" */
	"36",		/* 37. "EV_1KM_Emissive" */
};

// define array of wavelengths and assign correct values for emissive bands
double lambda[38] = {
	0,	// 0
	0,	// 1
	0,	// 2
	0,	// 3
	0,	// 4
	0,	// 5
	0,	// 6
	0,	// 7
	0,	// 8
	0,	// 9
	0,	// 10
	0,	// 11
	0,	// 12
	0,	// 13
	0,	// 14
	0,	// 15
	0,	// 16
	0,	// 17
	0,	// 18
	0,	// 19
	0,	// 20
	0,	// 21
	0.5*( 3.660 +  3.840)*1.0E-6, // 22. band 20
	0.5*( 3.929 +  3.989)*1.0E-6, // 23. band 21
	0.5*( 3.929 +  3.989)*1.0E-6, // 24. band 22
	0.5*( 4.020 +  4.080)*1.0E-6, // 25. band 23
	0.5*( 4.433 +  4.498)*1.0E-6, // 26. band 24
	0.5*( 4.482 +  4.549)*1.0E-6, // 27. band 25
	0.5*( 6.535 +  6.895)*1.0E-6, // 28. band 27 (sic)
	0.5*( 7.175 +  7.475)*1.0E-6, // 29. band 28
	0.5*( 8.400 +  8.700)*1.0E-6, // 30. band 29
	0.5*( 9.580 +  9.880)*1.0E-6, // 31. band 30
	0.5*(10.780 + 11.280)*1.0E-6, // 32. band 31
	0.5*(11.770 + 12.270)*1.0E-6, // 33. band 32
	0.5*(13.185 + 13.485)*1.0E-6, // 34. band 33
	0.5*(13.485 + 13.785)*1.0E-6, // 35. band 34
	0.5*(13.785 + 14.085)*1.0E-6, // 36. band 35
	0.5*(14.085 + 14.385)*1.0E-6, // 37. band 36
};

// k = Boltzmann gas constant (joules/Kelvin)
const float k_Boltz = 1.3806488e-23;
// h = Planck’s constant (joule * second)
const float h_Planck = 6.62606957e-34;
// c = speed of light in vacuum (m/s)
const float c_light = 299792458.0;


// Transfrom integers to radience and then to brightness temperature for emissive bands.
//
// is -- band number
// nx -- number of columns
// ny -- number of rows
// buff1 -- input image (1d)
// offset -- offset value for this band
// scale -- scale factor for this band
// maskNaN -- mask of pixels with negative radiance (1d) (preallocated output)
// inp_img -- brightness temperature output (2d) (preallocated output)
//
void
int2bt(int is, int nx, int ny, unsigned short *buff1, float offset, float scale,
	int *maskNaN, float **inp_img)
{
	float r1, r2;
	int nmask, ix, j;
	
	r1 = h_Planck*c_light/(k_Boltz*lambda[is]);
	r2 = lambda[is];
	r2 = 1.0e-6*(2.0*h_Planck*c_light*c_light)/(r2*r2*r2*r2*r2);

	// find the minimum valid radiance > 0, mask all pixels with radiance <= 0
	unsigned short jmin = 65535;
	nmask = 0;
	for(ix=0; ix<nx*ny; ix++) {
		maskNaN[ix] = 0;      // originally assume data are physically valid

		if(buff1[ix] <= offset) {
			// found a pixel with negative radiance
			maskNaN[ix] = 1;  // set the mask for unphysical data
			nmask++;          // count such pixels
			continue;
		}

		// update minimum physical radiance - used later as fill in value for unphysical values
		if(buff1[ix]<jmin) {
			jmin = buff1[ix];
		}
	}
	// printf("nmask = %i nx*ny = %i jmin = %i\n", nmask, nx*ny, jmin);

	// convert radiance to brightness temperature
	nmask = 0;
	double avebt = 0.0;

	for(ix=0; ix<nx*ny; ix++) {
		j = buff1[ix];

		// if radiance less than smallest physical value, set it to fill in value
		if(j<jmin) {
			j = jmin;
			nmask++;
		}

		// scale integers to physical radiance values
		inp_img[0][ix] = scale*(j - offset);

		// calculate Brightness Temperature from Radiance
		inp_img[0][ix] = r1/log(1.0 + r2/inp_img[0][ix]);

		// counter for average BT - just to check sanity of data
		avebt +=  inp_img[0][ix];
	}

	printf("Number of pixels with negative radiances on input = %i\n", nmask);
	// printf("Average Brightness Temperature = %e \n", avebt/(nx*ny));
}


// convert output brightness temperature back to radiance and then back to integer
//
// is -- band number
// nx -- number of columns
// ny -- number of rows
// outp_img -- brightness temperature image (2d)
// offset -- offset value for this band
// scale -- scale factor for this band
// maskNaN -- mask of pixels with negative radiance (1d)
// buff1 -- output integers (preallocated output)
//
void
bt2int(int is, int nx, int ny, float **outp_img, float offset, float scale, int *maskNaN, unsigned short *buff1)
{
	float r1, r2;
	int ix, j;
	float z;

	r1 = h_Planck*c_light/(k_Boltz*lambda[is]);
	r2 = lambda[is];
	r2 = 1.0e-6*(2.0*h_Planck*c_light*c_light)/(r2*r2*r2*r2*r2);
	
	for(ix=0; ix<nx*ny; ix++) {
		// preserve the original data
		// in case of unphysical negative radiance in original data (produces NaN in Brightness Temperature)
		if(maskNaN[ix] == 1) continue;

		// get the radiance from brightness temperature
		z = r2/(exp(r1/outp_img[0][ix]) - 1.0);

		// scale the radiance back to integer
		j = (int) round(z/scale + offset);

		// check that integer is within valid bounds
		if((j<0) || (j>65535)) {
			printf("Value outside range %i\n", j);
			j = 65535;
		}

		// copy value to output buffer
		buff1[ix] = (unsigned short) j;
	}
}


// convert scaled integers to physical reflectance
//
// nx -- number of columns
// ny -- number of rows
// buff1 -- input image (1d)
// offset -- offset value for this band
// scale -- scale factor for this band
// inp_img -- reflectance output (2d) (preallocated output)
//
void
int2ref(int nx, int ny, unsigned short *buff1, float offset, float scale, float **inp_img)
{
	int ix;
	
	for(ix=0; ix<nx*ny; ix++) {
		inp_img[0][ix] = scale*( ((float) (buff1[ix])) - offset);
	}
}

// convert reflectance back to scaled integers
//
// nx -- number of columns
// ny -- number of rows
// outp_img -- reflectance image (2d)
// offset -- offset value for this band
// scale -- scale factor for this band
// buff1 -- output integers (preallocated output)
//
void
ref2int(int nx, int ny, float **outp_img, float offset, float scale, unsigned short *buff1)
{
	int ix, j;
	
	for(ix=0; ix<nx*ny; ix++) {
		// scale the reflectance back to integer
		j = (int) round(outp_img[0][ix]/scale + offset);

		// check that integer is within valid bounds
		if((j<0) || (j>65535)) {
			printf("Value outside range %i\n", j);
			j = 65535;
		}

		// copy value to output buffer
		buff1[ix] = (unsigned short) j;
	}
}

char *progname;

static void
usage()
{
	printf("usage: %s MODIS_hdf_file MOD03_hdf_file destriping_param_file.txt\n", progname);
	printf("	-m	mask out overlapping region before resampling\n");
	printf("	-s	write sorted output image\n");
	exit(2);
}

#define GETARG(x)	do{\
		(x) = *argv++;\
		argc--;\
	}while(0);

int main(int argc, char** argv)
{
	char *flag;

	// there are four data fields in MODIS HDF files that contain all bands this code can process
	// here are the data field names for the four data fields
	char EV_250_Ref[] = "EV_250_Aggr1km_RefSB";
	char EV_500_Ref[] = "EV_500_Aggr1km_RefSB";
	char EV_1KM_Ref[] = "EV_1KM_RefSB";
	char EV_1KM_Emi[] = "EV_1KM_Emissive";
	char * dataFieldNames[4] = { EV_250_Ref, EV_500_Ref, EV_1KM_Ref, EV_1KM_Emi };

	// attribute names for the four data fields
	char AttrRef[] = "reflectance";
	char AttrRad[] = "radiance";
	char * attrBaseNames[4] = { AttrRef, AttrRef, AttrRef, AttrRad };

	// indices in bandNames array of first band in a data field
	int bandIndex[5] = { 0, 2, 7, 22, 38};

	unsigned short *buffer1 = NULL, *buff1;
	float **inp_img  = NULL;
	int    *maskNaN  = NULL;

	int is, status, i, j;

	int ib, nb, nx, ny, iband,  iDataField;
	int Ndet_arr[40], Niter_arr[40], isBand[40];
	float Qmin_arr[40], Qmax_arr[40], Tx_arr[40], Ty_arr[40], NEdQ_arr[40], Scale_arr[40], Offset_arr[40];


	// parse arguments
	GETARG(progname);
	bool maskoverlap = false;
	bool sortoutput = false;
	while(argc > 0 && strlen(argv[0]) == 2 && argv[0][0] == '-') {
		GETARG(flag);

		switch(flag[1]) {
		default:
			usage();
			break;
		case '-':
			goto argdone;
		case 'm':
			maskoverlap = true;
			break;
		case 's':
			sortoutput = true;
			break;
		}
	}
argdone:
	if(argc != 3)
		usage();
	char *hdfpath = argv[0];
	char *geopath = argv[1];
	char *parampath = argv[2];
	printf("modisresam %s%s%s %s %s\n",
		maskoverlap ? "-m " : "",
		sortoutput ? "-s " : "",
		hdfpath, geopath, parampath);

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// read input parameters from parameter file
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	printf("Input parameters\n");

	// open parameter file
	FILE *fp = fopen(parampath,"r");
	if(fp==NULL) {
		printf("ERROR: Cannot open parameter file\n");
		return -8;
	}

	// read parameters
	char tmpbandname[128];
	for(i=0; i<40; i++) {
		isBand[i] = 0;    // initially, set all band parameters as "absent"
	}
	for(i=0; i<40; i++) {                    // read at most 40 lines in parameter file

		j = fscanf(fp,"%s ", tmpbandname);   // read the first item in a line = band name
		if(j!=1) break;                      // if not read, then break

		// check if the band name is sensible
		is = -1;
		for(j=0; j<38; j++) {

			// compare the band name to all possible band names in an array
			if(strcmp(tmpbandname, bandNames[j]) == 0) {
				is = j;  // if found a match, assign band index
				break;   // and stop comparing
			}
		}
		if(is==-1) break; // if band name not valid, break

		// read all destriping parameters for this band
		j = fscanf(fp,"%i %i %f %f %f %f %f\n", &(Ndet_arr[is]), &(Niter_arr[is]), &(NEdQ_arr[is]),
		           &(Tx_arr[is]), &(Ty_arr[is]), &(Qmin_arr[is]), &(Qmax_arr[is]));

		if(j==7) {
			isBand[is]=1;
		} else break;  // if destriping parameters not read correctly, break

		// echo read destriping parameters
		printf(" %s \t %i %i %f %f %f %f %f\n", bandNames[is], (Ndet_arr[is]), (Niter_arr[is]), (NEdQ_arr[is]),
		       (Tx_arr[is]), (Ty_arr[is]), (Qmin_arr[is]), (Qmax_arr[is]));
	}

	fclose(fp); // close parameter file
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// done reading parameters
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// read latitude
	int latrows, latcols;
	float *lat;
	status = readlatitude(&lat, &latcols, &latrows, geopath);
	if(status<0) {
		printf("ERROR: Cannot read data Latitude data\n");
		return 10*status;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// loop over all 4 data fields
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for(iDataField=0; iDataField<4; iDataField++) {
		printf("========================================================================\n");
		printf("Data_field number = %i   name = %s\n", iDataField, dataFieldNames[iDataField]);
		printf("------------------------------------------------------------------------\n");
		ib = bandIndex[iDataField];                              // index of first band of this data field in bandNames[] array
		nb = bandIndex[iDataField+1] - bandIndex[iDataField];    // number of bands in this data field

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// make sure we have at least one band to resample in this data field
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		int nreadwrite = 0;
		for(iband=0; iband<nb; iband++) {
			if(isBand[ib+iband]>0) nreadwrite++;
		}
		if(nreadwrite==0) continue;   // if no bands to resample in this data field, continue

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// read data to resample in this data field
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		status = readwrite_modis( &buffer1, &nx, &ny, nb, &(Scale_arr[ib]), &(Offset_arr[ib]), &(isBand[ib]),
		                          dataFieldNames[iDataField], attrBaseNames[iDataField], hdfpath, 0);
		if(status<0) {
			printf("ERROR: Cannot read data field %s\n", dataFieldNames[iDataField]);
			return 10*status;
		}
		if(latrows != ny || latcols != nx){
			printf("ERROR: latitude image dimensions agree with band image\n");
			return 2;
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// allocate temporary arrays
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		inp_img  = allocate_2d_f(ny,nx);
		maskNaN  = (int *) malloc(ny*nx*sizeof(int));
		if( (maskNaN==NULL) || (inp_img==NULL)) {
			printf("ERROR: Cannot allocate memory\n");
			return -1;
		}


		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// go through all the bands in the current data field
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		int iBandIndx = 0;
		for(iband=0; iband<nb; iband++) {

			// index of the current band in bandNames[] array is a sum of
			//    index of first band of this data field in bandNames[] array and
			//    index of band within the current data field
			is = ib + iband;

			// printf("iband = %i isband = %i\n", iband, isBand[is]);
			if(isBand[is]==0) continue; // if no parameters for this band, then pass

			buff1 = &(buffer1[iBandIndx*nx*ny]); // location of the current band data to resample
			iBandIndx++;                         // increment the index of next band data to resample

			printf("Band = %i  MODIS_band_number = %s   scale = %e  offset = %e\n", iband, bandNames[is], Scale_arr[is], Offset_arr[is]);

			if(iDataField==3) {
				int2bt(is, nx, ny, buff1, Offset_arr[is], Scale_arr[is], maskNaN, inp_img);

				printf("Parameters:\n");
				printf("%i %i %f   %f %f   %f %f\n",  Ndet_arr[is], Niter_arr[is], NEdQ_arr[is], Tx_arr[is], Ty_arr[is], Qmin_arr[is], Qmax_arr[is]);

				resample_modis(inp_img, lat, nx, ny, Qmin_arr[is], Qmax_arr[is], maskoverlap, sortoutput);
				printf("Resampling done\n");

				bt2int(is, nx, ny, inp_img, Offset_arr[is], Scale_arr[is], maskNaN, buff1);


			}  //  if(iDataField==3)
			else {
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				// reflective band
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

				int2ref(nx, ny, buff1, Offset_arr[is], Scale_arr[is], inp_img);

				printf("Parameters:\n");
				printf("%i %i %f   %f %f   %f %f\n",  Ndet_arr[is], Niter_arr[is], NEdQ_arr[is], Tx_arr[is], Ty_arr[is], Qmin_arr[is], Qmax_arr[is]);
				resample_modis(inp_img, lat, nx, ny, Qmin_arr[is], Qmax_arr[is], maskoverlap, sortoutput);
				printf("Resampling done\n");

				ref2int(nx, ny, inp_img, Offset_arr[is], Scale_arr[is], buff1);

			}  //  if(iDataField!=3) meaning reflective band

			printf("------------------------------------------------------------------------\n");

		} // for iband

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// write resampled data back to hdf file, as well as set resampling attribute
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		status = readwrite_modis( &buffer1, &nx, &ny, nb, &(Scale_arr[ib]), &(Offset_arr[ib]), &(isBand[ib]),
		                          dataFieldNames[iDataField], attrBaseNames[iDataField], hdfpath, 1);

		if(status<0) {
			printf("ERROR: Failed to write data\n");
			return 10*status;
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// clean up
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(buffer1!=NULL) {
			free(buffer1);
			buffer1 = NULL;
		}
		free(inp_img[0]);
		free(inp_img);
		free(maskNaN);

	} // for(iDataField = ...

	return 0;
}
