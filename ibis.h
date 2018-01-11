#ifndef IBIS_H
#define IBIS_H

#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <omp.h>
#include <atomic>
#include <thread>

#define RENDER_SP_CONTOURS  0
#define THREAD_count        4
#define size_roi            9 //9 25 49

class IBIS
{

public:

    IBIS(int maxSPNum, int compacity);
    virtual ~IBIS();

    void process( cv::Mat* img );

    int getMaxSPNumber() { return maxSPNumber;}
    int getActualSPNumber() { return SPNumber; }

    float getComputationTime() { return st3; }
    float getPostProcessingTime() { return st4; }
    int* getLabels() { return labels; }

    float get_complexity() { return count_px_processed / size; }

protected:

    void initSeeds();
    void mask_propagate_SP();
    void mean_seeds();
    void generate_mask();
    bool angular_assign( int mask_index, int y, int x, int* angular, int* unique_angular, int index_unique );
    void fill_mask(int x_min, int x_max, int y_min, int y_max, int value);
    void assign_last( int y, int x, int x_min, int x_max, int y_min, int y_max, int* unique_angular, int index_unique );
    void apply_mask( int y, int x, int mask_index );
    int assign_px( int y, int x, int index_xy, int* unique_angular, int index_unique  );
    double now_ms(void);
    void enforceConnectivity();

private:

    struct Mask {
        int* x_var;
        int* y_var;
        int* angular;
        int size;
        int index_angular;

        void Mask_init( int size_in ) {
            x_var = new int[ size_in ];
            y_var = new int[ size_in ];
            angular = new int[ size_in ];
            size = size_in;

        }

        ~Mask() {
            delete[] x_var;
            delete[] y_var;
            delete[] angular;

        }

    };

    std::thread* thread_buff;
    int* vertical_index;
    int* adj_label;
    Mask* mask_buffer;
    int TPSP_type;
    bool iterate_FLAG;
    bool splitting_FLAG;
    int* mask_size;
    int* start_xy;
    int index_mask;
    int* unique_parent;
    bool* updated_px;
    int* adjacent_sp;
    int* initial_repartition;
    //int index_unique;

    float count_px_processed;

    int* x_vec;
    int* y_vec;

    int y_limit;
    int x_limit;

    int nb_iteration;
    float max_xy_dist;

    // V.1
    int count_reset;
    int size;
    int width;
    int height;

    int SPNumber;			// number of Sp actually created
    int SPTypicalLength;		// typical size of the width or height for a SP

    int* labels = nullptr;
    int* looking_area;
    int* count_looking_area;

    float* Xseeds;
    float* Yseeds;
    float* lseeds;
    float* aseeds;
    float* bseeds;

    float* lvec;
    float* avec;
    float* bvec;

    float* inv;

    bool bis_buffer;

    float* Xseeds_Sum;
    float* Yseeds_Sum;
    float* lseeds_Sum;
    float* aseeds_Sum;
    float* bseeds_Sum;

    float* countPx;
    float invwt;
    int* index_propagation;
    bool* elligible;
    float* count_diff;

    bool creation_deletion_FLAG;
    cv::Mat1b contourMask;
    float sensitivity;            // % of SP that needs to appear changed to be updated
    int inertia;                  // min || [ L, a, b ] || value to consider changes

    int minSPSizeThreshold;
    int compacity;                      // compacity factor
    int maxSPNumber;			// number of Sp passed by user

// RENDER_SNR_VALUES

public:
    double slicTime;
    double st1, st2, st3, st4, st5, st6;

    int slicNum;
    int selectedSp = -1;		// Rq : -1 is used in labels for pixels that are not included in any superpixel.
};

#endif // IBIS_H
