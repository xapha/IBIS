#include <iostream>
#include "ibis.h"
#include <opencv2/opencv.hpp>
#include <unistd.h>

int main( int argc, char* argv[] )
{
    printf(" - Iterative Boundary Identification Segmentation - \n\n");

    //printf( " omp threads : %i \n", omp_get_max_threads() );

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

    //K = 200;
    //compa = 50;
    //char* file_path = "16004.jpg";

    // IBIS
    IBIS Super_Pixel( K, compa );

    // get picture
    cv::Mat img = cv::imread( argv[ 3 ] );
    //cv::Mat img = cv::imread( file_path );

    // process IBIS
    Super_Pixel.process( &img );

    // convert int* labels to Mat* labels in gray scale
    cv::Mat labels_mat = cv::Mat( img.rows, img.cols, CV_8UC1 );
    int* labels = Super_Pixel.getLabels();

    for( int i=0; i < img.rows*img.cols; i++ ) {
        labels_mat.ptr()[ i ] = labels[ i ];

    }

    // save labels decomposition
    cv::imwrite( "labels.ppm", labels_mat );

    // print computation times
    float computation_time = Super_Pixel.getComputationTime();
    float post_processing_time = Super_Pixel.getPostProcessingTime();

    printf( " --> Computation time : IBIS             : \t %lf \t ms \n", computation_time );
    printf( " --> Computation time : Post processing  : \t %lf \t ms \n", post_processing_time );
    printf( " --> Computation time : total            : \t %lf \t ms \n", ( post_processing_time + computation_time ) );

    char command[ 1024 ];
    long pid = ( long ) getpid(); // use long to ensure correct format specifier
    sprintf( command, "pmap %ld > pmap_log.txt", pid );
    system( command );

}
