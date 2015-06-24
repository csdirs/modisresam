//
// Resampling of image based on latitude
//

#include "modisresam.h"
#include "sort.h"

#define SIGN(A)   ((A) > 0 ? 1 : ((A) < 0 ? -1 : 0 ))

// b in between a and c
#define INBETWEEN(a, b, c) (((a) <= (b) && (b) <= (c)) || ((c) <= (b) && (b) <= (a)))

enum {
	WIDTH_1KM = 1354,
	DEBUG = false,
};

// Generate a image of latitude sorting indices.
//
// sind -- sorting indices (output)
// swaths -- number of swaths in the output
//
void
getsortingind(Mat &sind, int swaths)
{
	int height = swaths*SWATH_SIZE;
	sind = Mat::zeros(height, WIDTH_1KM, CV_32SC1);
	
	int x = 0;
	for(int i = 0; i < (int)nelem(SORT_WIDTHS); i++){
		int xe = x + SORT_WIDTHS[i];
		for(; x < xe; x++){
			for(int y = 0; y < SWATH_SIZE; y++){
				sind.at<int>(y, x) = SORT_FIRST[i][y];
			}
			for(int y = SWATH_SIZE; y < height-SWATH_SIZE; y++){
				sind.at<int>(y, x) =
					(y/SWATH_SIZE)*SWATH_SIZE + SORT_MID[i][y%SWATH_SIZE];
			}
			for(int y = height-SWATH_SIZE; y < height; y++){
				sind.at<int>(y, x) =
					(y/SWATH_SIZE)*SWATH_SIZE + SORT_LAST[i][y%SWATH_SIZE];
			}
		}
	}
}

template <class T>
static Mat
resample_unsort_(const Mat &sind, const Mat &img)
{
	Mat newimg;
	int i, j, k;
	int32_t *sp;
	T *ip;

	CHECKMAT(sind, CV_32SC1);
	CV_Assert(img.channels() == 1);

	newimg = Mat::zeros(img.rows, img.cols, img.type());
	sp = (int32_t*)sind.data;
	ip = (T*)img.data;
	k = 0;
	for(i = 0; i < newimg.rows; i++) {
		for(j = 0; j < newimg.cols; j++) {
			newimg.at<T>(sp[k], j) = ip[k];
			k++;
		}
	}
	return newimg;
}

// Returns the unsorted image of the sorted image img.
// Sind is the image of sort indices.
static Mat
resample_unsort(const Mat &sind, const Mat &img)
{
	switch(img.type()) {
	default:
		eprintf("unsupported type %s\n", type2str(img.type()));
		break;
	case CV_8UC1:
		return resample_unsort_<uchar>(sind, img);
		break;
	case CV_32FC1:
		return resample_unsort_<float>(sind, img);
		break;
	case CV_64FC1:
		return resample_unsort_<double>(sind, img);
		break;
	}
	// not reached
	return Mat();
}

template <class T>
static Mat
resample_sort_(const Mat &sind, const Mat &img)
{
	Mat newimg;
	int i, j, k;
	int32_t *sp;
	T *np;

	CHECKMAT(sind, CV_32SC1);
	CV_Assert(img.channels() == 1);

	newimg = Mat::zeros(img.rows, img.cols, img.type());
	sp = (int*)sind.data;
	np = (T*)newimg.data;
	k = 0;
	for(i = 0; i < newimg.rows; i++){
		for(j = 0; j < newimg.cols; j++){
			np[k] = img.at<T>(sp[k], j);
			k++;
		}
	}
	return newimg;
}

// Returns the sorted image of the unsorted image img.
// Sind is the image of sort indices.
Mat
resample_sort(const Mat &sind, const Mat &img)
{
	switch(img.type()){
	default:
		eprintf("unsupported type %s\n", type2str(img.type()));
		break;
	case CV_8UC1:
		return resample_sort_<uchar>(sind, img);
		break;
	case CV_32FC1:
		return resample_sort_<float>(sind, img);
		break;
	case CV_64FC1:
		return resample_sort_<double>(sind, img);
		break;
	}
	// not reached
	return Mat();
}

// Returns the average of 2 values which can be NAN.
static double
avg2(double a, double b)
{
        if(isnan(a))
                return b;
        if(isnan(b))
                return a;
        return (a+b) / 2.0;
}

// Resample 1D data.
//
// sind -- sorting indices
// slat -- sorted latitude
// sval -- sorted values
// dir -- diff direction (-1 or 1)
// rval -- resampled values (intput & output)
//
static void
resample1d(const int *sind, const float *slat, const float *sval, int n, int stride, float *rval)
{
	int i;
	
	// copy first non-nan value for first row
	for(i = 0; i < n-stride; i += stride){
		if(!isnan(sval[i])){
			rval[0] = sval[i];
			break;
		}
	}
	
	// interpolate the middle values
	for(i = stride; i < n-stride; i += stride){
		if(SIGN(sind[i+stride] - sind[i]) == 1 && !isnan(sval[i])){
			rval[i] = sval[i];
			continue;
		}
		if(isnan(sval[i]) && isnan(sval[i-stride]) && isnan(sval[i+stride])){
			printf("unable to resample at row %d\n", i/stride);
			rval[i] = NAN;
			continue;
		}
		double x1 = (slat[i] + slat[i-stride]) / 2;
		double y1 = avg2(sval[i], sval[i-stride]);
		double x2 = (slat[i] + slat[i+stride]) / 2;
		double y2 = avg2(sval[i], sval[i+stride]);
		
		if(isnan(y1)){
			rval[i] = y2;
		}else if(isnan(y2)){
			rval[i] = y1;
		}else if(x2 == x1 || !INBETWEEN(x1, slat[i], x2)){
			// slat[i] might not be in between x1 and x2 because we're
			// using universal sorting indices
			rval[i] = (y1+y2) / 2;
		}else{
			double lam = (slat[i] - x1) / (x2 - x1);
			rval[i] = (1-lam)*y1 + lam*y2;
		}
	}
	
	// copy last non-nan value to last row
	for(int k = i; k >= 0; k -= stride){
		if(!isnan(sval[k])){
			rval[i] = sval[k];
			break;
		}
	}
}


// Resample a 2D image.
//
// ssrc -- image to resample already sorted
// slat -- sorted latitude
// sortidx -- lat sorting indices
// dst -- resampled image (output)
// 
static void
resample2d(const Mat &ssrc, const Mat &slat, const Mat &sortidx, Mat &dst)
{
	int j, width, height;
	Mat col, idx, botidx;
	Range colrg, toprg, botrg;

	CHECKMAT(ssrc, CV_32FC1);
	CHECKMAT(slat, CV_32FC1);
	CHECKMAT(sortidx, CV_32SC1);
	CV_Assert(ssrc.data != dst.data);

	width = ssrc.cols;
	height = ssrc.rows;
	int total = ssrc.total();
	dst = Mat::zeros(height, width, CV_32FC1);	// resampled values

	for(j = 0; j < width; j++) {
		// resample this column
		resample1d(&sortidx.ptr<int>(0)[j],
			&slat.ptr<float>(0)[j],
			&ssrc.ptr<float>(0)[j],
			total, width,
			&dst.ptr<float>(0)[j]);
	}
}


// Set overlapping regions to NAN.
//
static void
setoverlaps1km(Mat &dst, float value)
{
	enum {
		WIDTH_1KM = 1354,
		C0 = 0,
		C1 = 70,
		C2 = 70+130,
		C3 = WIDTH_1KM - C2,
		C4 = WIDTH_1KM - C1,
		C5 = WIDTH_1KM,
	};
	
	CHECKMAT(dst, CV_32FC1);
	
	if(dst.cols != WIDTH_1KM){
		eprintf("width of image is not %d\n", WIDTH_1KM);
	}

	for(int y = 0; y < dst.rows; y += SWATH_SIZE) {
		for(int x = C0; x < C2; x++) {
			dst.at<float>(y+0, x) = value;
		}
		for(int x = C0; x < C1; x++) {
			dst.at<float>(y+1, x) = value;
		}
		for(int x = C0; x < C1; x++) {
			dst.at<float>(y+8, x) = value;
		}
		for(int x = C0; x < C2; x++) {
			dst.at<float>(y+9, x) = value;
		}

		for(int x = C3; x < C5; x++) {
			dst.at<float>(y+0, x) = value;
		}
		for(int x = C4; x < C5; x++) {
			dst.at<float>(y+1, x) = value;
		}
		for(int x = C4; x < C5; x++) {
			dst.at<float>(y+8, x) = value;
		}
		for(int x = C3; x < C5; x++) {
			dst.at<float>(y+9, x) = value;
		}
	}
}

// Resample VIIRS swatch image _img with corresponding
// latitude image _lat.
// _img[0..ny][0..nx]  - original image (brightness temperature)
// _lat[0..ny][0..nx]  - original latitude
// nx = width of image (should be 3200 for VIIRS)
// ny = height of image ( 5408 or 5392 for ~10 min VIIRS granule)
//
// TODO: use min, max arguments
void
resample_modis(float **_img, float *_lat, int nx, int ny, float minvalid, float maxvalid,
	bool maskoverlap, bool sortoutput)
{
	Mat sind, dst;
	
	if(DEBUG) printf("resampling debugging is turned on!\n");
	
	// Mat wrapper around external buffer.
	// Caller of this function still reponsible for freeing the buffers.
	Mat img(ny, nx, CV_32FC1, &_img[0][0]);
	Mat lat(ny, nx, CV_32FC1, _lat);
	if(DEBUG)dumpmat("before.bin", img);
	if(DEBUG)dumpmat("lat.bin", lat);
	
	if(maskoverlap){
		// Set overlapping regions to NAN.
		// Those pixels are interpolated when resampling.
		setoverlaps1km(img, NAN);
		if(DEBUG)dumpmat("masked.bin", img);
	}
	
	getsortingind(sind, lat.rows/SWATH_SIZE);
	Mat slat = resample_sort(sind, lat);
	Mat simg = resample_sort(sind, img);
	if(DEBUG)dumpmat("sind.bin", sind);
	if(DEBUG)dumpmat("simg.bin", simg);
	if(DEBUG)dumpmat("slat.bin", slat);
	
	resample2d(simg, slat, sind, dst);
	if(DEBUG)dumpmat("after.bin", dst);

	if(!sortoutput){
		dst = resample_unsort(sind, dst);
	}
		
	CV_Assert(dst.size() == img.size() && dst.type() == img.type());
	dst.copyTo(img);
	if(DEBUG)dumpfloat("final.bin", &_img[0][0], nx*ny);
	if(DEBUG)exit(3);
}
