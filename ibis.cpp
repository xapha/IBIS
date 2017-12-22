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
    inv = new float[maxSPNumber];
    adjacent_sp = new int[size_roi*maxSPNumber];

}

IBIS::~IBIS() {
    delete[] labels;
    delete[] avec;
    delete[] lvec;
    delete[] bvec;

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

    delete[] elligible;
    delete[] count_diff;
    delete[] index_propagation;
    delete[] mask_buffer;
    delete[] mask_size;
    delete[] start_xy;
    delete[] unique_parent;
    delete[] updated_px;
    delete[] x_vec;
    delete[] y_vec;
    delete[] adj_label;
    delete[] vertical_index;
    delete[] looking_area;
    delete[] count_looking_area;

    delete[] adjacent_sp;
    delete[] initial_repartition;

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

void IBIS::get_looking_area() {
    double dist_xy;

    for( int y=0; y < height; y++ ) {
        for( int x=0; x< width; x++ ) {
            count_looking_area[ y*width + x ]=0;

            for( int index_sp = 0; index_sp < SPNumber; index_sp++ ) {
                dist_xy = (x - Xseeds[index_sp])*(x - Xseeds[index_sp]) +
                          (y - Yseeds[index_sp])*(y - Yseeds[index_sp]);

                if(dist_xy < max_xy_dist) {
                    looking_area[ (count_looking_area[ y*width + x ] * size) + (y*width + x ) ] = index_sp;
                    count_looking_area[ y*width + x ]++;

                }

            }

        }

    }

}

void IBIS::initSeeds() {
    const bool hexgrid = false;
    int n;
    int xstrips, ystrips;
    int xerr, yerr;
    double xerrperstrip, yerrperstrip;
    int xoff, yoff;
    int x, y;
    int xe, ye, xe_1, ye_1;
    int start, final, start_y, final_y, start_x, final_x;

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
        ye_1 = (int)((y+1)*yerrperstrip);

        int seedy = y * SPTypicalLength + yoff + ye;

        if( y == 0 ) {
            start_y = 0;
            final_y = SPTypicalLength + ye_1;

        }
        else {
            start_y = y * SPTypicalLength + ye;
            final_y = ( (y + 1) * SPTypicalLength + ye_1 >= height ) ? height-1 : (y + 1) * SPTypicalLength + ye_1;

        }

        for (x = 0; x < xstrips; x++)
        {
            int seedx;
            xe = (int)(x*xerrperstrip);
            xe_1 = (int)((x+1)*xerrperstrip);
            seedx = x * SPTypicalLength + xoff + xe;

            if( x == 0 ) {
                start_x = 0;
                final_x = SPTypicalLength + xe_1;

            }
            else {
                start_x = x * SPTypicalLength + xe;
                final_x = ( (x + 1) * SPTypicalLength + xe_1 >= width ) ? width-1 : (x + 1) * SPTypicalLength + xe_1;

            }

            int i = seedy * width + seedx;

            Xseeds[n] = (float) seedx;
            Yseeds[n] = (float) seedy;
            lseeds[n] = lvec[i];
            aseeds[n] = avec[i];
            bseeds[n] = bvec[i];

            // fill line by line
            for( int index_y=start_y; index_y<=final_y; index_y++ ) {
                std::fill( initial_repartition + index_y*width + start_x, initial_repartition + index_y*width + final_x, n );

            }

            // list adjacents seeds
            int ii=0;
            for( int roi_y=-(sqrt(size_roi) - 1)/2; roi_y <= (sqrt(size_roi) - 1)/2; roi_y++ ) {
                for( int roi_x=-(sqrt(size_roi) - 1)/2; roi_x <= (sqrt(size_roi) - 1)/2; roi_x++ ) {
                    if( y + roi_y < 0 || y + roi_y >= ystrips || x + roi_x < 0 || x + roi_x >= xstrips )
                        adjacent_sp[size_roi*n+ii] = -1;
                    else
                        adjacent_sp[size_roi*n+ii] = n + roi_y*xstrips + roi_x;

                    ii++;

                }
            }

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

void IBIS::mean_seeds() {
    memset( lseeds_Sum, 0, sizeof( float ) * maxSPNumber );
    memset( aseeds_Sum, 0, sizeof( float ) * maxSPNumber );
    memset( bseeds_Sum, 0, sizeof( float ) * maxSPNumber );
    memset( Xseeds_Sum, 0, sizeof( float ) * maxSPNumber );
    memset( Yseeds_Sum, 0, sizeof( float ) * maxSPNumber );
    std::fill( countPx, countPx + maxSPNumber , 0 );

    int current_sp;
    int i;

    for(int y=0; y<height; y++) {
        for(int x=0; x<width; x++) {
            // current index
            i = vertical_index[ y ] + x;
            current_sp = labels[ i ];

            if( current_sp >= 0 ) {
                lseeds_Sum[ current_sp ] += lvec[ i ] ;
                aseeds_Sum[ current_sp ] += avec[ i ] ;
                bseeds_Sum[ current_sp ] += bvec[ i ] ;

                Xseeds_Sum[ current_sp ] += x;
                Yseeds_Sum[ current_sp ] += y;
                countPx[ current_sp ]++;

            }

        }

    }

    for( int i=0; i<SPNumber; i++ ) {
        inv[ i ] = 1.f / countPx[ i ];
    }

    for( int i=0; i<SPNumber; i++ ) {
        if( countPx[ i ] > 0 ) {
            Xseeds[ i ] = Xseeds_Sum[ i ] * inv[ i ];
            Yseeds[ i ] = Yseeds_Sum[ i ] * inv[ i ];
            lseeds[ i ] = lseeds_Sum[ i ] * inv[ i ];
            aseeds[ i ] = aseeds_Sum[ i ] * inv[ i ];
            bseeds[ i ] = bseeds_Sum[ i ] * inv[ i ];

        }

    }

}

int IBIS::assign_px( int y, int x, int index_xy, int* unique_angular, int index_unique ) {
    int best_sp = -1;
    float dist_xy, dist_x, dist_y;
    float l, a, b;
    float dist_lab;
    int index_sp;
    float D=-1.f;
    float total_dist;

    count_px_processed++;

    for( int i=0; i<size_roi; i++ ) {
//    for( i=0; i<index_unique; i++ ) {
        index_sp = adjacent_sp[ size_roi*initial_repartition[index_xy] + i ];

//        index_sp = unique_angular[ i ];
        if( index_sp >= 0 && index_sp < SPNumber) {

            dist_x = Xseeds[ index_sp ] - x;
            dist_y = Yseeds[ index_sp ] - y;

            dist_xy = dist_x * dist_x + dist_y * dist_y;

            //if(dist_xy < max_xy_dist) {
                l = lvec[ index_xy ];
                a = avec[ index_xy ];
                b = bvec[ index_xy ];

                dist_lab = ( l - lseeds[ index_sp ]) * ( l - lseeds[ index_sp ]) +
                           ( a - aseeds[ index_sp ]) * ( a - aseeds[ index_sp ]) +
                           ( b - bseeds[ index_sp ]) * ( b - bseeds[ index_sp ]);

                dist_xy = ( x - Xseeds[ index_sp ] ) * ( x - Xseeds[ index_sp ] ) +
                          ( y - Yseeds[ index_sp ] ) * ( y - Yseeds[ index_sp ] );


                total_dist = dist_lab + dist_xy * invwt;

                if( total_dist < D || D < 0) {
                    best_sp = index_sp;
                    D = total_dist;

                }

            //}

        }

    }

    if( best_sp < 0 )
        printf("impossible value\n");

    return best_sp;
}

void IBIS::assign_last( int y, int x, int x_min, int x_max, int y_min, int y_max, int* unique_angular, int index_unique ) {

    int j, k;
    int index_xy;
    int value;
    int index_y;

    for( int index_var_y = y_min; index_var_y <= y_max; index_var_y++ ) {
        j = y + index_var_y;
        index_y = vertical_index[ j ];

        for( int index_var_x = x_min; index_var_x <= x_max; index_var_x++ ) {

            k = x + index_var_x;
            index_xy = index_y + k;

            if( j < height && k < width ) {

                if( labels[ index_xy ] < 0 ) {

                    // assign px
                    value = assign_px( j, k, index_xy, unique_angular, index_unique );

                    labels[ index_xy ] = value;

                }

            }

        }

    }

}

void IBIS::fill_mask( int x_min, int x_max, int y_min, int y_max, int value ) {

    int j;
    int index_y;

    for( int index_var_y = y_min; index_var_y <= y_max; index_var_y++ ) {
        j = index_var_y;
        index_y = vertical_index[ index_var_y ]; //index_var_y*width;

        std::fill( labels + index_y + x_min, labels + index_y + x_max, value );
        std::fill( initial_repartition + index_y + x_min, initial_repartition + index_y + x_max, value );


    }

}

bool IBIS::angular_assign( int mask_index, int y, int x, int* angular, int* unique_angular, int index_unique ) {

    int index_angular = 0;
    int j, k;
    bool unique_borders = true;
    int best_sp;
    int index_xy;

    for( int index_var=0; index_var < mask_buffer[ mask_index ].size; index_var++ ) {
        j = y + mask_buffer[ mask_index ].y_var[ index_var ];
        k = x + mask_buffer[ mask_index ].x_var[ index_var ];

        if( j < height && k < width ) {
            index_xy = vertical_index[ j ] + k;

            if( labels[ index_xy ] < 0 ) {
                // assign px
                best_sp = assign_px( j, k, index_xy, unique_angular, index_unique );
                angular[ index_angular ] = best_sp;

            }
            else {
                angular[ index_angular ] = labels[ index_xy ];

            }


            // check angular values
            if( angular[ 0 ] != angular[ index_angular ] ) {
                unique_borders = false;

                return unique_borders;

            }

            index_angular++;
        }

    }

    return unique_borders;

}



void IBIS::apply_mask( int y, int x, int mask_index ) {

    int limit_value;

    int* angular = new int[ mask_buffer[ mask_index ].size ];
    int* unique_angular = new int[ 20 ];
    int index_unique = 0;

    // identify possible seeds assignement
    float dist_x, dist_y;
    index_unique = 0;

    /*double lap = now_ms();
    for( int index_sp = 0; index_sp < SPNumber; index_sp++ ) {

        dist_x = Xseeds[ index_sp ] - x;
        dist_y = Yseeds[ index_sp ] - y;

        if( ( dist_x * dist_x + dist_y * dist_y ) < max_xy_dist ) {
            unique_angular[ index_unique ] = index_sp;
            index_unique++;

        }

    }
    st1 += now_ms() - lap;*/

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
            assign_last( y, x, -1, 1, -1, 1, unique_angular, index_unique );

        }

    }

    delete[] angular;
    delete[] unique_angular;

}

double IBIS::now_ms(void)
{
    double milliseconds_since_epoch = (double) (std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now().time_since_epoch()).count());
    //double milliseconds_since_epoch = (double)(std::chrono::high_resolution_clock::now().time_since_epoch() / std::chrono::milliseconds(1));
    //or double milliseconds_since_epoch = (double) (std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    return milliseconds_since_epoch;

}

// return the current actual SP number (the method may reduce the actual Sp number).
void IBIS::enforceConnectivity()
{
    //local var
    int label = 0;
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

    int step;
    for( int mask_index=index_mask-1; mask_index>=0; mask_index--) {

        y_limit = height + mask_size[ mask_index ];
        x_limit = width + mask_size[ mask_index ];

        if( mask_index == index_mask-1 )
            step = mask_size[ mask_index ];
        else
            step = mask_size[ mask_index ] - 1;


#if THREAD_count > 1
#pragma omp parallel for
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

        mean_seeds();

    }

}

void IBIS::process( cv::Mat* img ) {
    count_reset = 0;
    size = img->cols * img->rows;
    width = img->cols;
    height = img->rows;

    avec = new float[size];
    bvec = new float[size];
    lvec = new float[size];

    //store mean distance between 2 seeds info
    SPTypicalLength = (int)(std::sqrt((float)(size) / (float)(maxSPNumber))) + 1;
    invwt = (float)SPTypicalLength / compacity;
    invwt = 1.0f / (invwt * invwt);
    minSPSizeThreshold = (size / maxSPNumber) / 4;
    max_xy_dist = SPTypicalLength * SPTypicalLength * 4;

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
    count_px_processed = 0;

    labels = new int[size];
    looking_area = new int[10 * size];
    count_looking_area = new int[size];
    initial_repartition = new int[size];
    st1 = 0;

    std::fill( updated_px, updated_px + size, false );
    std::fill( countPx, countPx + maxSPNumber, 1.f );
    std::fill( labels, labels + size, -1 );

    memset( count_diff, 0, sizeof( int ) * maxSPNumber );

    // generates mask coordinates
    generate_mask();

    // STEP 0 : convert to Lab
    cv::Mat lab_image;
    cv::cvtColor(*img, lab_image, CV_BGR2Lab, 0);

    int ii = 0;
    for (int i = 0; i < size * 3; i += 3) {
        lvec[ii] = lab_image.ptr()[i];
        avec[ii] = lab_image.ptr()[i + 1];
        bvec[ii] = lab_image.ptr()[i + 2];
        ii++;
    }

    double lap;

    //lap = now_ms();
    initSeeds();
    //st1 = now_ms() - lap;

    lap = now_ms();
    mask_propagate_SP();
    st3 = now_ms() - lap;

    lap = now_ms();
    enforceConnectivity();
    st4 = now_ms() - lap;

    printf("PERF_T %lf\n", st3+st4);
    //printf("Seeds id \t %lf\n", st1);
}
