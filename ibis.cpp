#include "ibis.h"

IBIS::IBIS(int _maxSPNum, int _compacity ) : slicTime(0.), slicNum(0) {
    maxSPNumber = _maxSPNum;
    compacity = _compacity;

    // memory allocation
    Xseeds = new float[maxSPNumber];
    Yseeds = new float[maxSPNumber];
    lseeds = new float[maxSPNumber];
    aseeds = new float[maxSPNumber];
    bseeds = new float[maxSPNumber];

    Xseeds_Sum = new float[maxSPNumber];
    Yseeds_Sum = new float[maxSPNumber];
    lseeds_Sum = new float[maxSPNumber];
    aseeds_Sum = new float[maxSPNumber];
    bseeds_Sum = new float[maxSPNumber];

    countPx = new float[maxSPNumber];
    count_diff = new float[maxSPNumber];

    elligible = new bool[maxSPNumber];
    labMeanSeeds = new float[maxSPNumber * 3];

    // hard coded inner params
    inertia = 10;         // absolute distance in [ 0 255 ]
    sensitivity = 1.f;    // in percent

}

IBIS::~IBIS() {
    // allocated in setSize()
    //if ( labels )
    //{
        delete[] labels;
        delete[] previous_labels;
        delete[] avec;
        delete[] lvec;
        delete[] bvec;

        delete[] bis_avec;
        delete[] bis_lvec;
        delete[] bis_bvec;
    //}

    // allocated in constructor
    delete[] Xseeds;
    delete[] Yseeds;
    delete[] lseeds;
    delete[] aseeds;
    delete[] bseeds;

    delete[] Xseeds_Sum;
    delete[] Yseeds_Sum;
    delete[] lseeds_Sum;
    delete[] aseeds_Sum;
    delete[] bseeds_Sum;

    delete[] countPx;
    delete[] labMeanSeeds;

    delete[] elligible;
    delete[] count_diff;
    delete[] index_propagation;

    if( TPSP_type == 2 ) {
        delete[] mask_buffer;
        delete[] mask_size;
        delete[] start_xy;
        delete[] unique_parent;
        delete[] updated_px;
        delete[] x_vec;
        delete[] y_vec;
        delete[] adj_label;
        delete[] vertical_index;

    }

}

void IBIS::initSeeds() {
    const bool hexgrid = true;
    int n;
    int xstrips, ystrips;
    int xerr, yerr;
    double xerrperstrip, yerrperstrip;
    int xoff, yoff;
    int x, y;
    int xe, ye;

    xstrips = width / SPTypicalLength;
    ystrips = height / SPTypicalLength;

    xerr = width - SPTypicalLength * xstrips;
    yerr = height - SPTypicalLength * ystrips;

    xerrperstrip = (double)xerr / xstrips;
    yerrperstrip = (double)yerr / ystrips;

    xoff = SPTypicalLength / 2;
    yoff = SPTypicalLength / 2;

    n = 0;
    for (y = 0; y < ystrips; y++)
    {
        ye = (int)(y*yerrperstrip);
        for (x = 0; x < xstrips; x++)
        {
            int seedx;
            xe = (int)(x*xerrperstrip);

            if (hexgrid)
            {
                seedx = x * SPTypicalLength + (xoff << (y & 0x1)) + xe;
                if (seedx >= width)
                    seedx = width - 1;
            }
            else
            {
                seedx = (x * SPTypicalLength + xoff + xe);
            }

            int seedy = (y * SPTypicalLength + yoff + ye);
            int i = seedy * width + seedx;

            Xseeds[n] = (float) seedx;
            Yseeds[n] = (float) seedy;
            lseeds[n] = lvec[i];
            aseeds[n] = avec[i];
            bseeds[n] = bvec[i];

            Xseeds_Sum[n] = Xseeds[n];
            Yseeds_Sum[n] = Yseeds[n];
            lseeds_Sum[n] = lseeds[n];
            aseeds_Sum[n] = aseeds[n];
            bseeds_Sum[n] = bseeds[n];

            n++;
        }
    }
    SPNumber = n;

}

void IBIS::propagate_SP( ) {

    //MathLib::TimeDiffAccum tda(slicTime, !PERFORMANCE_MEASURES);

    int current_index;
    int previous_sp;
    float dist_xy, closest_dist;
    float dist_x, dist_y, dist_l, dist_a, dist_b;
    float dist_lab;
    float total_dist;
    float D;
    int best_sp, closest_sp;
    bool reassign;

    for(int y=0; y<height; y++) {
        for(int x=0; x<width; x++) {
            // current index
            reassign = false;
            current_index = y*width + x;

            previous_sp = labels[ current_index ];

            if( previous_sp < 0 )
                reassign = true;
            else if( elligible[ previous_sp ] == 1 )
                reassign = true;

            if( reassign ) {
                D = -1;
                best_sp = -1;
                closest_dist = -1;

                // get the closest SP in both dist_xy & dist_lab
                for( int index_sp=0; index_sp<SPNumber; index_sp++ ) {

                    dist_x = Xseeds[ index_sp ] - x;
                    dist_y = Yseeds[ index_sp ] - y;

                    dist_xy = dist_x * dist_x + dist_y * dist_y;

                    // save closest SP index
                    if( dist_xy < closest_dist || closest_dist < 0 ) {
                        closest_dist = dist_xy;
                        closest_sp = index_sp;
                    }

                    if( dist_xy < max_xy_dist ) {

                        if( bis_buffer ) {
                            dist_l = lseeds[ index_sp ] - bis_lvec[ current_index ];
                            dist_a = aseeds[ index_sp ] - bis_avec[ current_index ];
                            dist_b = bseeds[ index_sp ] - bis_bvec[ current_index ];

                        }
                        else {
                            dist_l = lseeds[ index_sp ] - lvec[ current_index ];
                            dist_a = aseeds[ index_sp ] - avec[ current_index ];
                            dist_b = bseeds[ index_sp ] - bvec[ current_index ];

                        }

                        dist_lab = dist_l * dist_l + dist_a * dist_a + dist_b * dist_b;
                        total_dist = dist_lab + ( dist_xy * invwt );
                        if( total_dist < D || D < 0) {
                            best_sp = index_sp;
                            D = total_dist;
                        }

                    }

                }

                if( best_sp < 0 )
                    best_sp = closest_sp;

                // update labels and seeds
                labels[ current_index ] = best_sp;

                //update seeds
                if( bis_buffer ) {
                    lseeds_Sum[ best_sp ] += bis_lvec[ current_index ];
                    aseeds_Sum[ best_sp ] += bis_avec[ current_index ];
                    bseeds_Sum[ best_sp ] += bis_bvec[ current_index ];

                }
                else {
                    lseeds_Sum[ best_sp ] += lvec[ current_index ];
                    aseeds_Sum[ best_sp ] += avec[ current_index ];
                    bseeds_Sum[ best_sp ] += bvec[ current_index ];

                }

                Xseeds_Sum[ best_sp ] += x;
                Yseeds_Sum[ best_sp ] += y;
                countPx[ best_sp ]++;

                // update seeds
                for( int i=0; i<SPNumber; i++ ) {
                    Xseeds[ i ] = Xseeds_Sum[ i ] / countPx[ i ];
                    Yseeds[ i ] = Yseeds_Sum[ i ] / countPx[ i ];
                    lseeds[ i ] = lseeds_Sum[ i ] / countPx[ i ];
                    aseeds[ i ] = aseeds_Sum[ i ] / countPx[ i ];
                    bseeds[ i ] = bseeds_Sum[ i ] / countPx[ i ];

                }

                /*mathLib->div(lseeds_Sum, countPx, lseeds, SPNumber);
                mathLib->div(aseeds_Sum, countPx, aseeds, SPNumber);
                mathLib->div(bseeds_Sum, countPx, bseeds, SPNumber);
                mathLib->div(Xseeds_Sum, countPx, Xseeds, SPNumber);
                mathLib->div(Yseeds_Sum, countPx, Yseeds, SPNumber);*/

            }

        }

    }

}

void IBIS::creation_deletion() {
    std::fill( elligible, elligible + maxSPNumber, 0 );
    creation_deletion_FLAG = false;
    int biggest_sp = 0;
    int x, y;
    int best_value;
    int current_index;

    for( int i=0; i<SPNumber; i++ ) {
        if( countPx[ i ] < minSPSizeThreshold && !elligible[ i ] ) {
            //get the index of the biggest sp
            best_value = 0;
            for( int j=0; j<SPNumber; j++ ) {
                if( countPx[ j ] > best_value ) {
                    best_value = countPx[ j ];
                    biggest_sp = j;
                }
            }

            // copy seeds
            Xseeds[ i ] = Xseeds[ biggest_sp ];
            Yseeds[ i ] = Yseeds[ biggest_sp ];
            lseeds[ i ] = lseeds[ biggest_sp ];
            aseeds[ i ] = aseeds[ biggest_sp ];
            bseeds[ i ] = bseeds[ biggest_sp ];

            elligible[ i ] = true;
            elligible[ biggest_sp ] = true;

            // historical inheritance
            index_propagation[ i ] = biggest_sp;

            // reinit SP
            countPx[ i ] = countPx[ biggest_sp ] = 1;

            Xseeds_Sum[ i ] = x = round( Xseeds[ biggest_sp ] );
            Yseeds_Sum[ i ] = y = round( Yseeds[ biggest_sp ] );

            current_index = vertical_index[ y ] + x;

            if( bis_buffer ) {
                lseeds_Sum[ i ] = round( bis_lvec[ current_index ] );
                aseeds_Sum[ i ] = round( bis_avec[ current_index ] );
                bseeds_Sum[ i ] = round( bis_bvec[ current_index ] );

            }
            else {
                lseeds_Sum[ i ] = round( lvec[ current_index ] );
                aseeds_Sum[ i ] = round( avec[ current_index ] );
                bseeds_Sum[ i ] = round( bvec[ current_index ] );

            }

            Xseeds_Sum[ biggest_sp ] = Xseeds_Sum[ i ];
            Yseeds_Sum[ biggest_sp ] = Yseeds_Sum[ i ];
            lseeds_Sum[ biggest_sp ] = lseeds_Sum[ i ];
            aseeds_Sum[ biggest_sp ] = aseeds_Sum[ i ];
            bseeds_Sum[ biggest_sp ] = bseeds_Sum[ i ];

            creation_deletion_FLAG = true;
            splitting_FLAG = true;
        }
        else
            index_propagation[ i ] = i;

    }

}

void IBIS::mean_seeds() {
    memset( lseeds_Sum, 0, sizeof( float ) * maxSPNumber );
    memset( aseeds_Sum, 0, sizeof( float ) * maxSPNumber );
    memset( bseeds_Sum, 0, sizeof( float ) * maxSPNumber );
    memset( Xseeds_Sum, 0, sizeof( float ) * maxSPNumber );
    memset( Yseeds_Sum, 0, sizeof( float ) * maxSPNumber );
    memset( countPx, 0, sizeof( float ) * maxSPNumber );

    int current_sp;
    int i;

    for(int y=0; y<height; y++) {
        for(int x=0; x<width; x++) {
            // current index
            i = vertical_index[ y ] + x;
            current_sp = labels[ i ];

            if( bis_buffer ) {
                lseeds_Sum[ current_sp ] += bis_lvec[ i ] ;
                aseeds_Sum[ current_sp ] += bis_avec[ i ];
                bseeds_Sum[ current_sp ] += bis_bvec[ i ];

            }
            else {
                lseeds_Sum[ current_sp ] += lvec[ i ] ;
                aseeds_Sum[ current_sp ] += avec[ i ] ;
                bseeds_Sum[ current_sp ] += bvec[ i ] ;

            }

            Xseeds_Sum[ current_sp ] += x;
            Yseeds_Sum[ current_sp ] += y;
            countPx[ current_sp ]++;

        }

    }

    for( int i=0; i<SPNumber; i++ ) {
        Xseeds[ i ] = Xseeds_Sum[ i ] / countPx[ i ];
        Yseeds[ i ] = Yseeds_Sum[ i ] / countPx[ i ];
        lseeds[ i ] = lseeds_Sum[ i ] / countPx[ i ];
        aseeds[ i ] = aseeds_Sum[ i ] / countPx[ i ];
        bseeds[ i ] = bseeds_Sum[ i ] / countPx[ i ];

    }

    /*mathLib->div(lseeds_Sum, countPx, lseeds, SPNumber);
    mathLib->div(aseeds_Sum, countPx, aseeds, SPNumber);
    mathLib->div(bseeds_Sum, countPx, bseeds, SPNumber);
    mathLib->div(Xseeds_Sum, countPx, Xseeds, SPNumber);
    mathLib->div(Yseeds_Sum, countPx, Yseeds, SPNumber);*/

}

void IBIS::renderMeanRGB(cv::Mat* pImg)
{
    int ii = 0, i;
    for (i = 0; i < 3 * size; i += 3, ii++)
    {
        int sp = labels[ii];

        if (sp >= 0)
        {
            pImg->ptr()[i] = (unsigned char)(labMeanSeeds[3 * sp]);
            pImg->ptr()[i + 1] = (unsigned char)(labMeanSeeds[3 * sp + 1]);
            pImg->ptr()[i + 2] = (unsigned char)(labMeanSeeds[3 * sp + 2]);
        }
    }
}

cv::Mat &IBIS::getLabelContourMask()
{
    return contourMask;
}

void IBIS::boundaries()
{
    int j, k, n, x, y;
    const int dx4[4] = { -1,  0,  1,  0 };
    const int dy4[4] = { 0, -1,  0,  1 };

    int buff_mask[4] = { 0 };
    bool inerPixel;

    //get the boundaries
    for (j = 0; j < height; j++) {
        for (k = 0; k < width; k++) {
            int ref_index = j*width + k;

            //get the 2 x 2 values
            inerPixel = true;

            for (n = 0; n < 4; n++)
            {
                x = k + dx4[n];
                y = j + dy4[n];

                if ((x >= 0 && x < width) && (y >= 0 && y < height))
                {
                    int nindex = y*width + x;
                    buff_mask[n] = labels[nindex];

                }
                else
                {
                    buff_mask[n] = -1;
                    inerPixel = false;
                    break;

                }

            }

            inerPixel = inerPixel && buff_mask[0] == buff_mask[2] && buff_mask[1] == buff_mask[3];

            //mark the boundaries and get the relative pixels for correspondent sp
            //if (flag0 && flag1) //if surrounding pixels are in the pics
            //	boundaries[ref_index] = 1;
            //else
            //	boundaries[ref_index] = 0;

            contourMask.ptr()[ref_index] = (inerPixel) ? 0 : 255;

        }

    }

}

int IBIS::compare (const void * a, const void * b) {
    return ( *(int*)a - *(int*)b );
}

void IBIS::find_unique_parent(int mask_index, int* angular, int* unique_angular , int *index_unique) {
    int local_index_unique = mask_buffer[ mask_index ].size;
    memcpy( unique_angular, angular, sizeof( int )*local_index_unique );
    qsort( unique_angular, local_index_unique, sizeof(int), IBIS::compare);

    // suppress duplicate
    int count_multi = 0;
    for( int i=1; i < local_index_unique; i++ ) {
        if( unique_angular[ i ] == unique_angular[ i - 1 ] ) {
            unique_angular[ i ] = SPNumber + 1; // invalid one
            count_multi++;
        }

    }

    // re order SP list
    qsort(unique_angular, local_index_unique, sizeof(int), IBIS::compare);

    // remove useless SP check
    *index_unique = local_index_unique - count_multi;

}

int IBIS::assign_px( int y, int x, int index_xy, int* unique_angular, int index_unique ) {
    int best_sp = -1;
    int closest_sp = -1;
    int i=0;
    float dist_xy, dist_x, dist_y, dist_l, dist_a, dist_b;
    float dist_lab;
    int index_sp;
    float closest_dist=-1.f;
    float D=-1.f;
    //int index_xy;
    float total_dist;

    for( i=0; i<index_unique; i++ ) {
        index_sp = unique_angular[ i ];

        dist_x = Xseeds[ index_sp ] - x;
        dist_y = Yseeds[ index_sp ] - y;

        dist_xy = dist_x * dist_x + dist_y * dist_y;

        // save closest SP index
        if( dist_xy < closest_dist || closest_dist < 0 ) {
            closest_dist = dist_xy;
            closest_sp = index_sp;
        }

        if( bis_buffer ) {
            dist_l = lseeds[ index_sp ] - bis_lvec[ index_xy ];
            dist_a = aseeds[ index_sp ] - bis_avec[ index_xy ];
            dist_b = bseeds[ index_sp ] - bis_bvec[ index_xy ];

        }
        else {
            dist_l = lseeds[ index_sp ] - lvec[ index_xy ];
            dist_a = aseeds[ index_sp ] - avec[ index_xy ];
            dist_b = bseeds[ index_sp ] - bvec[ index_xy ];

        }

        dist_lab = dist_l * dist_l + dist_a * dist_a + dist_b * dist_b;
        total_dist = dist_lab + ( dist_xy * invwt );

        if( total_dist < D || D < 0) {
            best_sp = index_sp;
            D = total_dist;
        }

    }

    if( best_sp < 0 )
        best_sp = closest_sp;

    return best_sp;
}

void IBIS::update_seeds( int y, int x, int current_sp, int previous_sp, int index_xy ) {

    if( bis_buffer ) {
        lseeds_Sum[ current_sp ] += bis_lvec[ index_xy ];
        aseeds_Sum[ current_sp ] += bis_avec[ index_xy ];
        bseeds_Sum[ current_sp ] += bis_bvec[ index_xy ];

    }
    else {
        lseeds_Sum[ current_sp ] += lvec[ index_xy ];
        aseeds_Sum[ current_sp ] += avec[ index_xy ];
        bseeds_Sum[ current_sp ] += bvec[ index_xy ];

    }

    Xseeds_Sum[ current_sp ] += x;
    Yseeds_Sum[ current_sp ] += y;

    countPx[ current_sp ]++;

}

void IBIS::assign_last( int y, int x, int x_min, int x_max, int y_min, int y_max, int* unique_angular, int index_unique ) {

    int j, k;
    int index_xy;
    int value;
    int previous_sp;
    int index_y;

    for( int index_var_y = y_min; index_var_y <= y_max; index_var_y++ ) {
        j = y + index_var_y;
        index_y = vertical_index[ j ];

        for( int index_var_x = x_min; index_var_x <= x_max; index_var_x++ ) {

            k = x + index_var_x;
            index_xy = index_y + k;

            if( j < height && k < width ) {

                if( labels[ index_xy ] < 0 ) {

                    previous_sp = previous_labels[ index_xy ];

                    if( previous_sp < 0 ) {
                        // assign px
                        value = assign_px( j, k, index_xy, unique_angular, index_unique );

                        if( !iterate_FLAG ) {
                            // update seeds
                            update_seeds( j, k, value, 0, index_xy );
                        }

                        labels[ index_xy ] = value;

                    }
                    else if( elligible[ previous_sp ] == 1 ) {
                        // assign px
                        value = assign_px( j, k, index_xy, unique_angular, index_unique );

                        if( !iterate_FLAG ) {
                            // update seeds
                            update_seeds( j, k, value, labels[ index_xy ], index_xy );
                        }

                        labels[ index_xy ] = value;

                    }
                    else {
                        if( !iterate_FLAG ) {
                            // update seeds
                            update_seeds( j, k, previous_sp, labels[ index_xy ], index_xy );
                        }

                        labels[ index_xy ] = previous_sp;

                    }

                }

            }

        }

    }

}

void IBIS::fill_mask( int x_min, int x_max, int y_min, int y_max, int value ) {

    int j, k;
    int index_xy;
    int previous_sp;
    int index_y;

    for( int index_var_y = y_min; index_var_y <= y_max; index_var_y++ ) {
        j = index_var_y;
        index_y = vertical_index[ index_var_y ]; //index_var_y*width;

        for( int index_var_x = x_min; index_var_x <= x_max; index_var_x++ ) {

            k = index_var_x;
            index_xy = index_y + k;

            if( labels[ index_xy ] < 0 ) {
                previous_sp = previous_labels[ index_xy ];

                if( previous_sp < 0 ) {
                    if( !iterate_FLAG ) {
                        // update seeds
                        update_seeds( j, k, value, 0, index_xy );
                    }

                    labels[ index_xy ] = value;

                }
                else if( elligible[ previous_sp ] == 1 ) {
                    if( !iterate_FLAG ) {
                        // update seeds
                        update_seeds( j, k, value, 0, index_xy );
                    }

                    labels[ index_xy ] = value;

                }
                else {
                    if( !iterate_FLAG ) {
                        // update seeds
                        update_seeds( j, k, previous_sp, 0, index_xy );
                    }

                    labels[ index_xy ] = previous_sp;

                }

            }

        }

    }

}

bool IBIS::angular_assign( int mask_index, int y, int x, int* angular, int* unique_angular, int index_unique ) {

    int index_angular = 0;
    int j, k;
    int previous_sp;
    bool unique_borders = true;
    int best_sp;
    int index_xy;

    for( int index_var=0; index_var < mask_buffer[ mask_index ].size; index_var++ ) {
        j = y + mask_buffer[ mask_index ].y_var[ index_var ];
        k = x + mask_buffer[ mask_index ].x_var[ index_var ];

        if( j < height && k < width ) {
            index_xy = vertical_index[ j ] + k;
            previous_sp = previous_labels[ index_xy ];

            if( labels[ index_xy ] < 0 ) {

                if( previous_sp < 0 ) { // first iteration
                    // assign px
                    best_sp = assign_px( j, k, index_xy, unique_angular, index_unique );

                    angular[ index_angular ] = best_sp;

                    if( iterate_FLAG ) {
                        labels[ index_xy ] = best_sp;

                    }

                }
                else if( elligible[ previous_sp ] == 1 ) { // update from motion
                    // assign px
                    best_sp = assign_px( j, k, index_xy, unique_angular, index_unique );

                    angular[ index_angular ] = best_sp;

                    if( iterate_FLAG ) {
                        labels[ index_xy ] = best_sp;

                    }

                }
                else { // conservation
                    angular[ index_angular ] = previous_sp;

                    if( iterate_FLAG ) {
                        labels[ index_xy ] = previous_sp;

                    }

                }

            }
            else {
                angular[ index_angular ] = labels[ index_xy ];

            }

            // check angular values
            if( angular[ 0 ] != angular[ index_angular ] ) {
                unique_borders = false;

                if( !iterate_FLAG && mask_index > 0 )
                    return unique_borders;

            }

            index_angular++;
        }

    }

    return unique_borders;

}

void IBIS::generate_mask() {
    mask_buffer = new Mask[ index_mask ];
    int limit_val, vertical_val;
    int value_assign;

    for( int k=0; k<index_mask; k++ ) {
        int table_index = 0;
        mask_buffer[ k ].Mask_init( pow(2.0, k+2) );

        if( k > 0 ) {
            limit_val = mask_size[ k - 1 ] - 1;
            vertical_val = -limit_val;

            for( int index_var=1; index_var <= mask_size[ k ]; index_var++ ) {
                if( index_var == 1 ) { //top border

                    value_assign = -limit_val;
                    for( int i=table_index; i<=table_index+limit_val; i++ ) {
                        mask_buffer[ k ].x_var[ i ] = value_assign;
                        mask_buffer[ k ].y_var[ i ] = vertical_val;

                        value_assign += 2;
                    }

                    table_index += limit_val + 1;
                    vertical_val += 2;

                }
                else if( index_var > 1 && index_var < mask_size[ k - 1 ] ) { // vertical border

                    value_assign = -limit_val;
                    for( int i=table_index; i<table_index+2; i++ ) {
                        mask_buffer[ k ].x_var[ i ] = value_assign;
                        mask_buffer[ k ].y_var[ i ] = vertical_val;

                        value_assign = limit_val;
                    }

                    table_index += 2;
                    vertical_val += 2;

                }
                else { // bot border

                    value_assign = -limit_val;
                    for( int i=table_index; i<=table_index+limit_val; i++ ) {
                        mask_buffer[ k ].x_var[ i ] = value_assign;
                        mask_buffer[ k ].y_var[ i ] = limit_val;

                        value_assign += 2;
                    }

                }

            }

        }
        else {

            mask_buffer[ k ].x_var[ 0 ] = -1;
            mask_buffer[ k ].y_var[ 0 ] = -1;

            mask_buffer[ k ].x_var[ 1 ] = 1;
            mask_buffer[ k ].y_var[ 1 ] = -1;

            mask_buffer[ k ].x_var[ 2 ] = -1;
            mask_buffer[ k ].y_var[ 2 ] = 1;

            mask_buffer[ k ].x_var[ 3 ] = 1;
            mask_buffer[ k ].y_var[ 3 ] = 1;

        }

    }

}

void IBIS::apply_mask( int y, int x, int mask_index ) {

    int limit_value;
    int x_iter, y_iter;
    int offset;

    int* angular = new int[ mask_buffer[ mask_index ].size ]; //(int*)malloc( sizeof( int ) * mask_buffer[ mask_index ].size );
    int* unique_angular = new int[ 20 ];//(int*)malloc( sizeof( int ) * 20 );
    int index_unique = 0;

    // identify possible seeds assignement
    float dist_x, dist_y;
    index_unique = 0;
    for( int index_sp = 0; index_sp < SPNumber; index_sp++ ) {

        dist_x = Xseeds[ index_sp ] - x;
        dist_y = Yseeds[ index_sp ] - y;

        if( ( dist_x * dist_x + dist_y * dist_y ) < max_xy_dist ) {
            unique_angular[ index_unique ] = index_sp;
            index_unique++;

        }

    }

    if( index_unique == 0 )
        std::cout << "alerte aux gogoles !!!!" << std::endl;

    if( angular_assign( mask_index, y, x, angular, unique_angular, index_unique ) ) {
        limit_value = ( mask_size[ mask_index ] - 1 )/2;

        if( x - limit_value < width && y - limit_value < height ) {
            fill_mask( ( x - limit_value < 0 ) ? 0 : x - limit_value ,
                       ( x + limit_value >= width ) ? width - 1 : x + limit_value ,
                       ( y - limit_value < 0 ) ? 0 : y - limit_value,
                       ( y + limit_value >= height ) ? height - 1 : y + limit_value,
                        angular[ 0 ] );


        }

    }
    else {
        if( mask_index == 0 ) {
            //find_unique_parent( mask_index, angular, unique_angular, &index_unique );
            assign_last( y, x, -1, 1, -1, 1, unique_angular, index_unique );

        }
        else if( iterate_FLAG ) {
            mask_index--;
            offset = 1 << mask_index;
            for( int iter=0; iter < 4; iter ++ ) {
                if( iter == 0 ) {
                    x_iter = x - offset;
                    y_iter = y - offset;

                }
                else if( iter == 1 ) {
                    x_iter = x + offset;
                    y_iter = y - offset;

                }
                else if( iter == 2 ) {
                    x_iter = x - offset;
                    y_iter = y + offset;

                }
                else {
                    x_iter = x + offset;
                    y_iter = y + offset;

                }

                apply_mask( y_iter, x_iter, mask_index );

            }

        }

    }

    delete[] angular;
    delete[] unique_angular;

    //free( angular );
    //free( unique_angular );

}

double IBIS::now_ms(void)
{
    double milliseconds_since_epoch = (double) (std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now().time_since_epoch()).count());
    //double milliseconds_since_epoch = (double)(std::chrono::high_resolution_clock::now().time_since_epoch() / std::chrono::milliseconds(1));
    //or double milliseconds_since_epoch = (double) (std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    return milliseconds_since_epoch;


    // (clock_t) clock() * 1000 / CLOCKS_PER_SEC
}

// return the current actual SP number (the method may reduce the actual Sp number).
void IBIS::enforceConnectivity()
{
    //local var
    int label;
    int i, j, k;
    int n, c, count;
    int x, y;
    int ind;
    int oindex, adjlabel;
    int nindex;
    const int dx4[4] = { -1,  0,  1,  0 };
    const int dy4[4] = { 0, -1,  0,  1 };
    int* nlabels = new int[ size ];
    std::fill( nlabels, nlabels + size, -1 );

    oindex = 0;
    adjlabel = 0;

    for (j = 0; j < height; j++)
    {
        for (k = 0; k < width; k++)
        {
            if (nlabels[oindex] < 0)
            {
                nlabels[oindex] = label;// !! labels[oindex] --> label

                x_vec[0] = k;
                y_vec[0] = j;

                for (n = 0; n < 4; n++)
                {
                    x = x_vec[0] + dx4[n];
                    y = y_vec[0] + dy4[n];

                    if ((x >= 0 && x < width) && (y >= 0 && y < height))
                    {
                        nindex = y*width + x;

                        if (nlabels[nindex] >= 0)
                            adjlabel = nlabels[nindex];
                    }
                }

                count = 1;
                for (c = 0; c < count; c++)
                {
                    for (n = 0; n < 4; n++)
                    {
                        x = x_vec[c] + dx4[n];
                        y = y_vec[c] + dy4[n];
                        if ((x >= 0 && x < width) && (y >= 0 && y < height))
                        {
                            nindex = y*width + x;

                            if (nlabels[nindex] < 0 && labels[oindex] == labels[nindex])
                            {
                                x_vec[count] = x;
                                y_vec[count] = y;
                                nlabels[nindex] = label;
                                count++;
                            }
                        }
                    }
                }

                if (count <= minSPSizeThreshold)
                {
                    for (c = 0; c < count; c++)
                    {
                        ind = y_vec[c] * width + x_vec[c];
                        nlabels[ind] = adjlabel;
                    }
                    label--;
                }
                label++;
            }
            oindex++;
        }
    }

    for (i = 0; i < size; i++)
        labels[i] = nlabels[i];

    delete[] nlabels;

}

void IBIS::mask_propagate_SP() {

    //MathLib::TimeDiffAccum tda(slicTime, !PERFORMANCE_MEASURES);

    if( iterate_FLAG ) {

        y_limit = height + mask_size[ index_mask-1 ];
        x_limit = width + mask_size[ index_mask-1 ];

//#pragma omp parallel for num_threads( 4 )
        for( int y=start_xy[ index_mask-1 ]; y<y_limit; y+=mask_size[ index_mask-1 ] ) {
            for( int x=start_xy[ index_mask-1 ]; x<x_limit; x+=mask_size[ index_mask-1 ] ) {

                apply_mask( y, x, index_mask-1 );

            }

        }

    }
    else {
        int step;

        for( int mask_index=index_mask-1; mask_index>=0; mask_index--) {

            y_limit = height + mask_size[ mask_index ];
            x_limit = width + mask_size[ mask_index ];

            if( mask_index == index_mask-1 )
                step = mask_size[ mask_index ];
            else
                step = mask_size[ mask_index ] - 1;

            //double lap, diff_lap;
            //lap = now_ms();


#if THREAD_count > 1
#pragma omp parallel for num_threads( THREAD_count )
#endif
            for( int y=start_xy[ mask_index ]; y<y_limit; y+=step ) { // splitted loop
                for( int x=start_xy[ mask_index ]; x<x_limit; x+=step ) {

                    if( y < height && x < width ) {
                        if( labels[ vertical_index[ y ] + x ] < 0 )
                            apply_mask( y, x, mask_index );

                    }
                    else
                        apply_mask( y, x, mask_index );

                }

            }
            //diff_lap = now_ms() - lap;
            //printf( " iteration time : %lf ms \n", diff_lap );

            /*if( THREAD_count >= 1 ) {
                // thread division
                int y_limit_thread;
                int thread_limit[ THREAD_count + 1 ] = { 0 };

                //#pragma omp parallel for num_threads( 4 )

                y_limit_thread = y_limit / THREAD_count;

                for( int index_thread=0; index_thread<THREAD_count+1; index_thread++ ) {
                    if( index_thread == 0 )
                        thread_limit[ index_thread ] = start_xy[ mask_index ];
                    else if( index_thread == THREAD_count )
                        thread_limit[ index_thread ] = y_limit;
                    else {
                        while( thread_limit[ index_thread ] <= y_limit_thread * ( index_thread ) )
                            thread_limit[ index_thread ] += step;

                    }

                    if( index_thread > 0 ) {
                        thread_buff[ index_thread-1 ] = std::thread( IBIS::thread_mask_propagate,
                                                                   this, index_thread-1, thread_limit, step, mask_index );

                    }

                }

                //join threads
                for( int index_thread=0; index_thread<THREAD_count; index_thread++ )
                    thread_buff[ index_thread ].join();

            }*/

            // update seeds
            //#pragma omp simd
            for( int i=0; i<SPNumber; i++ ) {
                Xseeds[ i ] = Xseeds_Sum[ i ] / countPx[ i ];
                Yseeds[ i ] = Yseeds_Sum[ i ] / countPx[ i ];
                lseeds[ i ] = lseeds_Sum[ i ] / countPx[ i ];
                aseeds[ i ] = aseeds_Sum[ i ] / countPx[ i ];
                bseeds[ i ] = bseeds_Sum[ i ] / countPx[ i ];

            }

            /*mathLib->div(lseeds_Sum, countPx, lseeds, SPNumber);
            mathLib->div(aseeds_Sum, countPx, aseeds, SPNumber);
            mathLib->div(bseeds_Sum, countPx, bseeds, SPNumber);
            mathLib->div(Xseeds_Sum, countPx, Xseeds, SPNumber);
            mathLib->div(Yseeds_Sum, countPx, Yseeds, SPNumber);*/

        }

    }

}

void IBIS::thread_mask_propagate( IBIS* SP, int thread_index, int* thread_limit, int step, int mask_index ) {

    double lap, diff_lap;
    lap = SP->now_ms();
    printf(" THREAD NUM : %i => start \n", thread_index);

    for( int y=thread_limit[ thread_index ]; y<thread_limit[ thread_index + 1 ]; y+=step ) { // splitted loop

        for( int x=SP->start_xy[ mask_index ]; x<SP->x_limit; x+=step ) {

            if( y < SP->height && x < SP->width ) {
                if( SP->labels[ SP->vertical_index[ y ] + x ] < 0 )
                    SP->apply_mask( y, x, mask_index );

            }
            else
                SP->apply_mask( y, x, mask_index );

        }

    }

    diff_lap = SP->now_ms() - lap;
    printf( " THREAD NUM : %i => end, time : %lf ms \n", thread_index, diff_lap );

}

void IBIS::assure_contiguity() {
    std::fill( updated_px, updated_px + size, false );

    int dx4[] = { -1, 0, 1, 0 };
    int dy4[] = { 0, -1, 0, 1 };

    int current_sp, prospect_sp;
    int j, k;
    int current_index, prospect_index;

    // assure that a SP is a single 4-connected cluster
    int nb_adj;
    int count_vec;

    bool inscription;

    for(int y=0; y<height; y++) {
        for(int x=0; x<width; x++) {
            // current index
            current_index = vertical_index[ y ] + x;

            if( !updated_px[ current_index ] ) {
                nb_adj = 1;
                current_sp = labels[  current_index  ];
                adj_label[ 0 ] = current_sp;

                // get adjacent px for direct neighbour
                x_vec[ 0 ] = x;
                y_vec[ 0 ] = y;
                count_vec = 1;

                for( int index_vec = 0; index_vec <= count_vec; index_vec++  ) {
                    for( int index_var=0; index_var<4; index_var++ ) {
                        j = y_vec[ index_vec ] + dy4[ index_var ];
                        k = x_vec[ index_vec ] + dx4[ index_var ];
                        prospect_index = vertical_index[ j ] + k;

                        if( j >= 0 && j < height && k >= 0 && k < width ) {
                            prospect_sp = labels[ prospect_index ];

                            if( !updated_px[ prospect_index ] && prospect_sp == current_sp ) {

                                count_vec++;
                                x_vec[ count_vec ] = k;
                                y_vec[ count_vec ] = j;

                                updated_px[ prospect_index ] = true;

                            }
                            else {
                                if( prospect_sp != current_sp ) {
                                    inscription = true;

                                    // check if we know already about this one
                                    for( int index_adj=0; index_adj<nb_adj; index_adj++ ) {
                                        if( adj_label[ index_adj ] == prospect_sp ) {
                                            inscription = false;
                                            break;
                                        }

                                    }

                                    if( inscription ) {
                                        adj_label[ nb_adj ] = prospect_sp;
                                        nb_adj++;

                                    }

                                }

                            }

                        }

                    }

                }

                // single pixel
                if( count_vec == 1 )
                    updated_px[ current_index ] = true;

                // check if only one adjacent SP, if so, reallocate current SP pixels to the global one
                if( nb_adj <= 3 ) {

                    for( int index_vec = 0; index_vec < count_vec; index_vec++ ) {

                        prospect_index = vertical_index[ y_vec[ index_vec ] ] + x_vec[ index_vec ];

                        if( nb_adj == 1 )
                            labels[ prospect_index ] = adj_label[ 0 ];
                        else
                            labels[ prospect_index ] = adj_label[ 1 ];

                    }

                }

            }

        }

    }

}

void IBIS::process( cv::Mat* img ) {

    double t0, t1, t2, t3, t4, t5, t6;
    //t0 = mathLib->now_ms();

    if (labels == nullptr) { // first frame

        count_reset = 0;
        bis_buffer = false;
        size = img->cols * img->rows;
        width = img->cols;
        height = img->rows;

        previous_labels = new int[size];
        avec = new float[size];
        bvec = new float[size];
        lvec = new float[size];

        bis_avec = new float[size];
        bis_bvec = new float[size];
        bis_lvec = new float[size];
        //boundaries = new int[size];
        index_propagation = new int[maxSPNumber];

        contourMask = cv::Mat1b(height, width);

        //store mean distance between 2 seeds info
        SPTypicalLength = (int)(std::sqrt((float)(size) / (float)(maxSPNumber))) + 1;
        invwt = (float)SPTypicalLength / compacity;
        invwt = 1.0f / (invwt * invwt);
        minSPSizeThreshold = (size / maxSPNumber) / 4;
        max_xy_dist = SPTypicalLength * SPTypicalLength * 4;

        // select TPSP version = 1: old version, temporal propag opti / 2: new version, iterative mask
        TPSP_type = 2;

        // prepare mask
        if( TPSP_type == 2 ) {
            iterate_FLAG = false;
            index_mask = 1;

            while( SPTypicalLength > ( pow( 2.0, index_mask ) + 1 ) )
                index_mask++;
            index_mask--;

            mask_size = new int[ index_mask ];
            start_xy = new int[ index_mask ];

            for( int k=0; k<index_mask; k++ ) {
                mask_size[ k ] = pow( 2.0, k+1 ) + 1;

                if( k == 0 )
                    start_xy[ k ] = 1;
                else
                    start_xy[ k ] = mask_size[ k - 1 ] - 1;
            }

            unique_parent = new int[ mask_size[ index_mask - 1 ] ];
            updated_px = new bool[ size ];

            x_vec = new int[ size ];
            y_vec = new int[ size ];

            vertical_index = new int[ height + mask_size[ index_mask - 1 ] ];
            for( int k=0; k < height + mask_size[ index_mask - 1 ]; k++ ) {
                vertical_index[ k ] = k*width;
            }

            adj_label = new int[ 20 ];

            // generates mask coordinates
            generate_mask();

        }

    }

    // STEP 0 : convert to Lab
    cv::Mat lab_image;
    cv::cvtColor(*img, lab_image, CV_BGR2Lab, 0);

    int ii = 0;
    if( bis_buffer ) {
        for (int i = 0; i < size * 3; i += 3) {
            bis_lvec[ii] = lab_image.ptr()[i];
            bis_avec[ii] = lab_image.ptr()[i + 1];
            bis_bvec[ii] = lab_image.ptr()[i + 2];
            ii++;
        }
    }
    else {
        for (int i = 0; i < size * 3; i += 3) {
            lvec[ii] = lab_image.ptr()[i];
            avec[ii] = lab_image.ptr()[i + 1];
            bvec[ii] = lab_image.ptr()[i + 2];
            ii++;
        }
    }

    // reset variable
    creation_deletion_FLAG = false;
    splitting_FLAG = false;

    //SP_PERF_MEASURE(1, 0);
    // STEP 1: compute ||diff|| from consecutive frames
    if( labels != nullptr ) {

        //set_elligible();
        std::fill( elligible, elligible + maxSPNumber, true );
        std::fill( countPx, countPx + maxSPNumber, 1.f );

        // copy labels info into previous labels
        memcpy( previous_labels, labels, sizeof( int ) * size );

        std::fill( labels, labels + size, -1 );
        std::fill( updated_px, updated_px + size, false );

        for( int n=0; n<SPNumber; n++) {
            Xseeds_Sum[n] = Xseeds[n];
            Yseeds_Sum[n] = Yseeds[n];
            lseeds_Sum[n] = lseeds[n];
            aseeds_Sum[n] = aseeds[n];
            bseeds_Sum[n] = bseeds[n];
        }

    }
    else {
        labels = new int[size];

        std::fill( elligible, elligible + maxSPNumber, false );
        std::fill( updated_px, updated_px + size, false );
        std::fill( countPx, countPx + maxSPNumber, 1.f );
        std::fill( labels, labels + size, -1 );
        std::fill( previous_labels, previous_labels + size, -1 );
        memset( count_diff, 0, sizeof( int ) * maxSPNumber );

        initSeeds();

    }

    //SP_PERF_MEASURE(2, 1);
    // STEP 2: propagate SP through the new frame

    double lap, diff_lap;
    lap = now_ms();

    if( TPSP_type == 2 )
        mask_propagate_SP();
    else
        propagate_SP();

    diff_lap = now_ms() - lap;
    st3 = diff_lap;

    //SP_PERF_MEASURE(3, 2);
    // STEP 3: check for creation deletion, ==> STEP 4
    // assure contiguity and compute mean_seeds

    lap = now_ms();
    if( iterate_FLAG ) {
        //assure_contiguity();
        mean_seeds();

    }

    enforceConnectivity();

    //creation_deletion();

    //if( creation_deletion_FLAG )
    //    propagate_SP();

    //assure_contiguity();

    diff_lap = now_ms() - lap;
    st4 = diff_lap;
}
