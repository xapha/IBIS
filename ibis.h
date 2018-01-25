#ifndef IBIS_H
#define IBIS_H

#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <omp.h>
#include <atomic>
#include <thread>
#include "utils.h"

// algo debug parameters
#define THREAD_count        4       // deactivate if <= 1
#define size_roi            9       // 9 25 49 : consider adjacent seeds for pixels assignation
#define MATLAB_lab          0       // 0=>RGB2LAB : l,a,b -> 0-255. 1=>RGB2LAB : l -> 0-100; a,b -> -127:127
#define MASK_chrono         0       // 0:1 for evaluation purpose : slow down the process !
#define VISU                0       // show processed pixels for each iterations
#define VISU_all            0       // for mask by mask visu of the processing : THREAD_count = 0, very slow
#define OUTPUT_log          1       // 0:1 print output log

class IBIS
{

public:
    friend class MASK;

    IBIS(int maxSPNum, int compacity);
    virtual ~IBIS();

    void process( cv::Mat* img );
    void init();
    void reset();

    int getMaxSPNumber() { return maxSPNumber;}
    int getActualSPNumber() { return SPNumber; }

    float getComputationTime() { return st3; }
    float getPostProcessingTime() { return st4; }
    int* getLabels() { return labels; }

    float get_complexity();

protected:
    void getLAB(cv::Mat *img);
    void initSeeds();
    void mask_propagate_SP();
    void mean_seeds();

    double now_ms(void);
    void enforceConnectivity();

private:

    class MASK {

    private:
        int* x_var;
        int* y_var;
        int* xy_var;
        int* last_px_x;
        int* last_px_y;
        int* last_px_xy;
        int* last_parent;

        int count_var;
        int count_last;
        int count_last_parent;

        int* limit_value_fill;

        int x;
        int y;
        int size;
        int mask_index;
        int angular;
        bool filled;
        MASK* sub_mask;
        IBIS* IBIS_data;

        void generate_mask();

    public:
        MASK() {

        }
        ~MASK();

        void init(int size_in, int y_in, int x_in , int mask_level, IBIS *ptr_IBIS);
        void reset();
        void process( int mask_level );
        int assign_px(int y, int x, int index_xy);
        int assign_last_px(int y, int x, int index_xy);
        void assign_last();
        void fill_mask();
        bool angular_assign();
        void assign_labels( int x, int y, int index_xy, int value );

    };

    // input parameters
    int size;
    int width;
    int height;

    // internals buffer
    int* vertical_index;
    MASK* mask_buffer;
    int* mask_size;
    int* adjacent_sp;
    int* count_adjacent;
    int* initial_repartition;
    int* processed;
    int* x_vec;
    int* y_vec;
    int* labels;
    float* inv;
    float* Xseeds_Sum;
    float* Yseeds_Sum;
    float* lseeds_Sum;
    float* aseeds_Sum;
    float* bseeds_Sum;
    float* countPx;

    // inner parameter
    int count_mask;
    int start_xy;
    int index_mask;
    int y_limit;
    int x_limit;
    float invwt;
    int minSPSizeThreshold;

    int SPNumber;               // number of Sp actually created
    int SPTypicalLength;		// typical size of the width or height for a SP
    int compacity;              // compacity factor
    int maxSPNumber;            // number of Sp passed by user

    // seeds value
    float* Xseeds;
    float* Yseeds;
    float* lseeds;
    float* aseeds;
    float* bseeds;

    // initial seeds repartition
    float* Xseeds_init;
    float* Yseeds_init;

    // image data
    float* lvec;
    float* avec;
    float* bvec;

public:
    double slicTime;
    double st1, st2, st3, st4, st5, st6;
    double Mt1, Mt2, Mt3, Mt4;

    int slicNum;
    int selectedSp = -1;		// Rq : -1 is used in labels for pixels that are not included in any superpixel.
};

#endif // IBIS_H
