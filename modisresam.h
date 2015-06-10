#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <opencv2/opencv.hpp>

using namespace cv;

#define	nelem(x)	(sizeof(x)/sizeof((x)[0]))

// allocate_2d.cc
float	** allocate_2d_f(int n1, int n2);
int	**allocate_2d_i(int n1, int n2);

// readwrite_modis.cc
int	readwrite_modis(unsigned short ** buffer, int * nx, int * ny, int nband, float *scales, float *offsets,
                    int *isband, char * sds_name, char * attr_name, char * filename, int readwrite);
int	readlatitude(float ** buffer, int *nx, int *ny, const char *filename);

// utils.cc
const char	*type2str(int type);
void	eprintf(const char *fmt, ...);
void	dumpmat(const char *filename, Mat &m);
void	dumpfloat(const char *filename, float *buf, int nbuf);

// resample_modis.cc
void	resample_modis(float **_img, float *_lat, int nx, int ny, float min, float max, bool maskoverlap);
