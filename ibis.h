#ifndef IBIS_H
#define IBIS_H

#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <omp.h>
#include <atomic>
#include <thread>

#define RENDER_SP_CONTOURS  0
#define THREAD_count        1

class IBIS
{
public:

    IBIS(int maxSPNum, int compacity);
    virtual ~IBIS();

    void process( cv::Mat* img );
    cv::Mat &getLabelContourMask();

    int getMaxSPNumber() { return maxSPNumber;}
    int getActualSPNumber() { return SPNumber; }
    float *getRGBValues(int &numSP) { numSP = SPNumber; return labMeanSeeds; }
    int *getPropagation() { return index_propagation; }
    void renderMeanRGB(cv::Mat* pImg);

    float getComputationTime() { return st3; }
    float getPostProcessingTime() { return st4; }
    int* getLabels() { return labels; }

protected:

    void initSeeds();
    void boundaries();

    // STEP 1: propagate SP through the new frame
    void mask_propagate_SP();

    // STEP 2: check for creation deletion
    void assure_contiguity();

    // STEP 2: check for creation deletion
    void creation_deletion();

    // STEP 4: update seeds values
    void mean_seeds();

    // sub functions
    void generate_mask();
    bool angular_assign( int mask_index, int y, int x, int* angular, int* unique_angular, int index_unique );
    void fill_mask(int x_min, int x_max, int y_min, int y_max, int value);
    void assign_last( int y, int x, int x_min, int x_max, int y_min, int y_max, int* unique_angular, int index_unique );
    void apply_mask( int y, int x, int mask_index );
    int assign_px( int y, int x, int index_xy, int* unique_angular, int index_unique  );
    void update_seeds( int y, int x, int current_sp, int previous_sp, int index_xy );
    static int compare (const void * a, const void * b);
    void find_unique_parent( int mask_index, int* angular, int* unique_angular, int* index_unique );
    void propagate_SP();
    static void thread_mask_propagate( IBIS* , int thread_index, int* thread_limit, int step, int mask_index );
    double now_ms(void);
    void enforceConnectivity();

private:

    // V.2
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
    //int index_unique;
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
    int* previous_labels;

    float* Xseeds;
    float* Yseeds;
    float* lseeds;
    float* aseeds;
    float* bseeds;

    float* lvec;
    float* avec;
    float* bvec;

    float* bis_lvec;
    float* bis_avec;
    float* bis_bvec;

    float* labMeanSeeds;

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
