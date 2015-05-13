//
// Resampling of image based on latitude
//

#include "modisresam.h"

#define SGN(A)   ((A) > 0 ? 1 : ((A) < 0 ? -1 : 0 ))
#define CHECKMAT(M, T)	CV_Assert((M).type() == (T) && (M).isContinuous())

enum {
	VIIRS_SWATH_SIZE = 16,
	MODIS_SWATH_SIZE = 10,
};

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
	for(i = 0; i < newimg.rows; i++) {
		for(j = 0; j < newimg.cols; j++) {
			np[k] = img.at<T>(sp[k], j);
			k++;
		}
	}
	return newimg;
}

// Returns the sorted image of the unsorted image img.
// Sind is the image of sort indices.
static Mat
resample_sort(const Mat &sind, const Mat &img)
{
	switch(img.type()) {
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

// Returns the average of 3 pixels.
static double
avg3(double a, double b, double c)
{
	if(isnan(b))
		return NAN;
	if(isnan(a) || isnan(c))
		return b;
	return (a+b+c)/3.0;
}

// Returns the average filter of image 'in' with a window of 3x1
// where sorted order is not the same as the original order.
// Sind is the sort indices giving the sort order.
static Mat
avgfilter3(const Mat &in, const Mat &sind)
{
	const int *sindp;
	const float *ip;
	Mat out;
	int i, j, rows, cols;
	float *op;

	CHECKMAT(in, CV_32FC1);
	CHECKMAT(sind, CV_32SC1);
	rows = in.rows;
	cols = in.cols;

	out.create(rows, cols, CV_32FC1);
	in.row(0).copyTo(out.row(0));
	in.row(rows-1).copyTo(out.row(rows-1));

	for(i = 1; i < rows-1; i++) {
		ip = in.ptr<float>(i);
		op = out.ptr<float>(i);
		sindp = sind.ptr<int>(i);
		for(j = 0; j < cols; j++) {
			if(sindp[j] != i)
				op[j] = avg3(ip[j-cols], ip[j], ip[j+cols]);
			else
				op[j] = ip[j];
		}
	}
	return out;
}

// Interpolate the missing values in image simg and returns the result.
// Slat is the latitude image, and slandmask is the land mask image.
// All input arguments must already be sorted.
static Mat
resample_interp(const Mat &simg, const Mat &slat)
{
	int i, j, k, nbuf, *buf;
	Mat newimg, bufmat;
	double x, llat, rlat, lval, rval;

	CHECKMAT(simg, CV_32FC1);
	CHECKMAT(slat, CV_32FC1);

	newimg = simg.clone();
	bufmat = Mat::zeros(simg.rows, 1, CV_32SC1);
	buf = (int*)bufmat.data;

	for(j = 0; j < simg.cols; j++) {
		nbuf = 0;
		llat = -999;
		lval = NAN;
		for(i = 0; i < simg.rows; i++) {
			// valid pixel
			if(!isnan(simg.at<float>(i, j))) {
				// first pixel is not valid, so extrapolate
				if(llat == -999) {
					for(k = 0; k < nbuf; k++) {
						newimg.at<float>(buf[k],j) = simg.at<float>(i, j);
					}
					nbuf = 0;
				}

				// interpolate pixels in buffer
				for(k = 0; k < nbuf; k++) {
					rlat = slat.at<float>(i, j);
					rval = simg.at<float>(i, j);
					x = slat.at<float>(buf[k], j);
					newimg.at<float>(buf[k],j) =
					    lval + (rval - lval)*(x - llat)/(rlat - llat);
				}

				llat = slat.at<float>(i, j);
				lval = simg.at<float>(i, j);
				nbuf = 0;
				continue;
			}

			// no valid pixel
			buf[nbuf++] = i;
		}
		// extrapolate the last pixels
		if(llat != -999) {
			for(k = 0; k < nbuf; k++) {
				newimg.at<float>(buf[k],j) = lval;
			}
		}
	}
	return newimg;
}

enum Pole {
	NORTHPOLE,
	SOUTHPOLE,
	NOPOLE,
};
typedef enum Pole Pole;

// Argsort latitude image 'lat' with given swath size.
// Image of sort indices are return in 'sortidx'.
static void
argsortlat(const Mat &lat, int swathsize, Mat &sortidx)
{
	int i, j, off, width, height, dir, d, split;
	Pole pole;
	Mat col, idx, botidx;
	Range colrg, toprg, botrg;

	CHECKMAT(lat, CV_32FC1);
	CV_Assert(swathsize >= 2);
	CV_Assert(lat.data != sortidx.data);

	width = lat.cols;
	height = lat.rows;
	sortidx.create(height, width, CV_32SC1);

	// For a column in latitude image, look at every 'swathsize' pixels
	// starting from 'off'. If they increases and then decreases, or
	// decreases and then increases, we're at the polar region.
	off = swathsize/2;

	pole = NOPOLE;

	for(j = 0; j < width; j++) {
		col = lat.col(j);

		// find initial direction -- increase, decrease or no change
		dir = 0;
		for(i = off+swathsize; i < height; i += swathsize) {
			dir = SGN(col.at<float>(i) - col.at<float>(i-swathsize));
			if(dir != 0)
				break;
		}

		// find change in direction if there is one
		for(; i < height; i += swathsize) {
			d = SGN(col.at<float>(i) - col.at<float>(i-swathsize));
			if(dir == 1 && d == -1) {
				CV_Assert(pole == NOPOLE || pole == NORTHPOLE);
				pole = NORTHPOLE;
				break;
			}
			if(dir == -1 && d == 1) {
				CV_Assert(pole == NOPOLE || pole == SOUTHPOLE);
				pole = SOUTHPOLE;
				break;
			}
		}

		if(i >= height) {
			pole = NOPOLE;
			if(dir >= 0)
				sortIdx(col, sortidx.col(j), CV_SORT_EVERY_COLUMN + CV_SORT_ASCENDING);
			else
				sortIdx(col, sortidx.col(j), CV_SORT_EVERY_COLUMN + CV_SORT_DESCENDING);
			continue;
		}

		split = i-swathsize;	// split before change in direction
		colrg = Range(j, j+1);
		toprg = Range(0, split);
		botrg = Range(split, height);

		if(pole == NORTHPOLE) {
			botidx = sortidx(botrg, colrg);
			sortIdx(col.rowRange(toprg), sortidx(toprg, colrg),
			        CV_SORT_EVERY_COLUMN + CV_SORT_ASCENDING);
			sortIdx(col.rowRange(botrg), botidx,
			        CV_SORT_EVERY_COLUMN + CV_SORT_DESCENDING);
			botidx += split;
		} else {	// pole == SOUTHPOLE
			botidx = sortidx(botrg, colrg);
			sortIdx(col.rowRange(toprg), sortidx(toprg, colrg),
			        CV_SORT_EVERY_COLUMN + CV_SORT_DESCENDING);
			sortIdx(col.rowRange(botrg), botidx,
			        CV_SORT_EVERY_COLUMN + CV_SORT_ASCENDING);
			botidx += split;
		}
	}
}

void
setoverlaps1km(Mat &dst, float value)
{
	enum {
		Width1KM = 1354,
		C0 = 0,
		C1 = 70,
		C2 = 70+130,
		C3 = Width1KM - C2,
		C4 = Width1KM - C1,
		C5 = Width1KM,
	};
	
	CHECKMAT(dst, CV_32FC1);
	
	if(dst.cols != Width1KM){
		eprintf("width of image is not %d\n", Width1KM);
	}

	for(int y = 0; y < dst.rows; y += MODIS_SWATH_SIZE) {
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
resample_modis(float **_img, float *_lat, int nx, int ny, float minvalid, float maxvalid)
{
	Mat sind, simg;

	// Mat wrapper around external buffer.
	// Caller of this function still reponsible for freeing the buffers.
	Mat img(ny, nx, CV_32FC1, &_img[0][0]);
	Mat lat(ny, nx, CV_32FC1, _lat);
	if(false)dumpmat("before.bin", img);
	
	// set overlapping regions to NAN,
	// so those pixels are interpolated when resampling
	setoverlaps1km(img, NAN);
	if(false)dumpmat("masked.bin", img);

	argsortlat(lat, MODIS_SWATH_SIZE, sind);
	simg = resample_sort(sind, img);
	if(false)dumpmat("simg1.bin", simg);
	simg = avgfilter3(simg, sind);
	if(false)dumpmat("simg2.bin", simg);
	lat = resample_sort(sind, lat);
	simg = resample_interp(simg, lat);
	if(false)dumpmat("simg3.bin", simg);
	simg = resample_unsort(sind, simg);
	if(false)dumpmat("simg4.bin", simg);

	CV_Assert(simg.size() == img.size() && simg.type() == img.type());
	simg.copyTo(img);
	if(false)dumpfloat("final.bin", &_img[0][0], nx*ny);
}
