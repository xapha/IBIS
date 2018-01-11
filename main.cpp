#include <iostream>
#include "ibis.h"
#include <opencv2/opencv.hpp>
#include <unistd.h>
#include <cmath>
#include <fstream>
#include <highgui.h>

using namespace std;
//=================================================================================
/// DrawContoursAroundSegments
///
/// Internal contour drawing option exists. One only needs to comment the if
/// statement inside the loop that looks at neighbourhood.
//=================================================================================
void DrawContoursAroundSegments(
    unsigned int*&			ubuff,
    int*&					labels,
    const int&				width,
    const int&				height,
    const unsigned int&				color )
{
    const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
    const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

    int sz = width*height;
    vector<bool> istaken(sz, false);
    vector<int> contourx(sz);vector<int> contoury(sz);
    int mainindex(0);int cind(0);

    for( int j = 0; j < height; j++ )
    {
        for( int k = 0; k < width; k++ )
        {
            int np(0);
            for( int i = 0; i < 8; i++ )
            {
                int x = k + dx8[i];
                int y = j + dy8[i];

                if( (x >= 0 && x < width) && (y >= 0 && y < height) )
                {
                    int index = y*width + x;

                    //if( false == istaken[index] )//comment this to obtain internal contours
                    {
                        if( labels[mainindex] != labels[index] ) np++;
                    }
                }
            }
            if( np > 1 )
            {
                contourx[cind] = k;
                contoury[cind] = j;
                istaken[mainindex] = true;
                //img[mainindex] = color;
                cind++;
            }
            mainindex++;
        }
    }

    int numboundpix = cind;//int(contourx.size());
    for( int j = 0; j < numboundpix; j++ )
    {
        int ii = contoury[j]*width + contourx[j];
        ubuff[ii] = 0xffffff;

        for( int n = 0; n < 8; n++ )
        {
            int x = contourx[j] + dx8[n];
            int y = contoury[j] + dy8[n];
            if( (x >= 0 && x < width) && (y >= 0 && y < height) )
            {
                int ind = y*width + x;
                if(!istaken[ind])
                    ubuff[ind] = 0;
            }
        }
    }
}


void write_labels(IplImage* input, const std::string& output_labels)
{
    std::ofstream file;
    file.open(output_labels.c_str());

    unsigned char* data = (unsigned char*)input->imageData;
    for (int y=0 ; y<input->height ; y++)
    {
        for (int x=0 ; x<input->width-1 ; x++)
        {
            file << (int) data[y*input->widthStep + x*input->nChannels] << " ";
        }
        file << (int) data[y*input->widthStep + (input->width -1)*input->nChannels] << std::endl;
    }

    file.close();
}

std::string get_name(const std::string& path_with_ext)
{
    int deb = path_with_ext.find_last_of("/");
    int fin = path_with_ext.find_last_of(".");
    return path_with_ext.substr(deb+1, fin-deb-1);
}

int main( int argc, char* argv[] )
{
    printf(" - Iterative Boundary Identification Segmentation - \n\n");

    int K;
    int compa;

    if( argc != 4 ) {
        printf(" --> usage ./ibis SP_number Compacity File_path \n");
        return -1;

    }
    else {
        K = atoi( argv[ 1 ] );
        compa = atoi( argv[ 2 ] );

        if( K < 0 || compa < 0 ) {
            printf(" --> usage ./ibis SP_number Compacity File_path \n");
            return -1;

        }

    }

    // IBIS
    IBIS Super_Pixel( K, compa );

    // get picture
    cv::Mat img = cv::imread( argv[ 3 ] );

    // process IBIS
    Super_Pixel.process( &img );

    // convert int* labels to Mat* labels in gray scale
    int* labels = Super_Pixel.getLabels();

    const int width = img.cols;
    const int height = img.rows;
    std::string output_basename = argv[3];
    output_basename = get_name(output_basename);
    const int color = 0xFFFFFFFF;
    std::string output_labels = output_basename;
    output_labels = std::string("results/") + output_labels + std::string(".seg");

    ofstream outfile;
    outfile.open(output_labels.c_str());
    for (int y=0 ; y<height ; y++)
    {
        for (int x=0 ; x<width-1 ; x++)
        {
            outfile << labels[y*width + x] << " ";
        }
        outfile << labels[y*width + width-1] << std::endl;
    }
    outfile.close();

    IplImage* output_bounds_alpha = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 4);
    cvSet(output_bounds_alpha, cvScalar(0,0,0,0));
    IplImage* output_bounds = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 3);
    unsigned int* ubuff = (unsigned int*)output_bounds_alpha->imageData;
    DrawContoursAroundSegments(ubuff, labels, width, height, color);

    cvCvtColor(output_bounds_alpha, output_bounds, CV_RGBA2RGB);

    std::string output_boundaries_name = output_basename;
    output_boundaries_name = std::string("results/") + output_boundaries_name + std::string("_boundary.png");
    cvSaveImage(output_boundaries_name.c_str(), output_bounds);

    cvReleaseImage(&output_bounds_alpha);
    cvReleaseImage(&output_bounds);

    // print computation times
    float computation_time = Super_Pixel.getComputationTime();
    float post_processing_time = Super_Pixel.getPostProcessingTime();

    printf( " --> Computation time : IBIS             : \t %lf \t ms \n", computation_time );
    printf( " --> Computation time : Post processing  : \t %lf \t ms \n", post_processing_time );
    printf( " --> Computation time : total            : \t %lf \t ms \n", ( post_processing_time + computation_time ) );
    printf( " --> Pixels processed : total            : \t %lf \t \% \n", Super_Pixel.get_complexity()*100 );

}
