/********************************************************
 * Kernels to be optimized for the CS:APP Performance Lab
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
/*
 * Please fill in the following team_t struct
 */
team_t team = {
    "TEAM", /* Team Name */

    "e244842",            /* First student ID */
    "Taha Emir Gökçegöz", /* First student name */

    "e230974",               /* Second student ID */
    "Muhammed Batuhan Berk", /* Second student name */

    "e244846",        /* Third student ID */
    "Kerem Recep Gür" /* Third student Name */
};

/********************
 * CONVOLUTION KERNEL
 ********************/

/***************************************************************
 * Your different versions of the convolution functions  go here
 ***************************************************************/

/*
 * naive_conv - The naive baseline version of convolution
 */
char naive_conv_descr[] = "naive_conv: Naive baseline implementation";
void naive_conv(int dim, pixel *src2, pixel *ker, unsigned *dst)
{
    int i, j, k, l;

    for (i = 0; i < dim - 8 + 1; i++)
        for (j = 0; j < dim - 8 + 1; j++)
        {
            dst[RIDX(i, j, dim)] = 0;
            for (k = 0; k < 8; k++)
                for (l = 0; l < 8; l++)
                {
                    dst[RIDX(i, j, dim)] += src2[RIDX((i + k), (j + l), dim)].red * ker[RIDX(k, l, 8)].red;
                    dst[RIDX(i, j, dim)] += src2[RIDX((i + k), (j + l), dim)].green * ker[RIDX(k, l, 8)].green;
                    dst[RIDX(i, j, dim)] += src2[RIDX((i + k), (j + l), dim)].blue * ker[RIDX(k, l, 8)].blue;
                }
        }
}

/*
 * convolution - Your current working version of convolution
 * IMPORTANT: This is the version you will be graded on
 */
char convolution_descr[] = "Convolution: Current working version";
void convolution(int dim, pixel *src, pixel *ker, unsigned *dst)
{
    int i, j, acc, dim2, num, dim3, dim4;
    pixel *src2, *src3;
    unsigned *dst2;

    int k0, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15;
    int k16, k17, k18, k19, k20, k21, k22, k23, k24, k25, k26, k27, k28, k29, k30, k31;
    int k32, k33, k34, k35, k36, k37, k38, k39, k40, k41, k42, k43, k44, k45, k46, k47;
    int k48, k49, k50, k51, k52, k53, k54, k55, k56, k57, k58, k59, k60, k61, k62, k63;
    int l0, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15;
    int l16, l17, l18, l19, l20, l21, l22, l23, l24, l25, l26, l27, l28, l29, l30, l31;
    int l32, l33, l34, l35, l36, l37, l38, l39, l40, l41, l42, l43, l44, l45, l46, l47;
    int l48, l49, l50, l51, l52, l53, l54, l55, l56, l57, l58, l59, l60, l61, l62, l63;
    int m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15;
    int m16, m17, m18, m19, m20, m21, m22, m23, m24, m25, m26, m27, m28, m29, m30, m31;
    int m32, m33, m34, m35, m36, m37, m38, m39, m40, m41, m42, m43, m44, m45, m46, m47;
    int m48, m49, m50, m51, m52, m53, m54, m55, m56, m57, m58, m59, m60, m61, m62, m63;
    num = dim - 7;
    dim2 = 8;
    dim3 = 8;
    dim4 = 8;
    k0 = ker[0].red;
    k1 = ker[1].red;
    k2 = ker[2].red;
    k3 = ker[3].red;
    k4 = ker[4].red;
    k5 = ker[5].red;
    k6 = ker[6].red;
    k7 = ker[7].red;
    k8 = ker[0 + dim2].red;
    k9 = ker[1 + dim2].red;
    k10 = ker[2 + dim2].red;
    k11 = ker[3 + dim2].red;
    k12 = ker[4 + dim2].red;
    k13 = ker[5 + dim2].red;
    k14 = ker[6 + dim2].red;
    k15 = ker[7 + dim2].red;
    dim2 += 8;
    k16 = ker[0 + dim2].red;
    k17 = ker[1 + dim2].red;
    k18 = ker[2 + dim2].red;
    k19 = ker[3 + dim2].red;
    k20 = ker[4 + dim2].red;
    k21 = ker[5 + dim2].red;
    k22 = ker[6 + dim2].red;
    k23 = ker[7 + dim2].red;
    dim2 += 8;
    k24 = ker[0 + dim2].red;
    k25 = ker[1 + dim2].red;
    k26 = ker[2 + dim2].red;
    k27 = ker[3 + dim2].red;
    k28 = ker[4 + dim2].red;
    k29 = ker[5 + dim2].red;
    k30 = ker[6 + dim2].red;
    k31 = ker[7 + dim2].red;
    dim2 += 8;
    k32 = ker[0 + dim2].red;
    k33 = ker[1 + dim2].red;
    k34 = ker[2 + dim2].red;
    k35 = ker[3 + dim2].red;
    k36 = ker[4 + dim2].red;
    k37 = ker[5 + dim2].red;
    k38 = ker[6 + dim2].red;
    k39 = ker[7 + dim2].red;
    dim2 += 8;
    k40 = ker[0 + dim2].red;
    k41 = ker[1 + dim2].red;
    k42 = ker[2 + dim2].red;
    k43 = ker[3 + dim2].red;
    k44 = ker[4 + dim2].red;
    k45 = ker[5 + dim2].red;
    k46 = ker[6 + dim2].red;
    k47 = ker[7 + dim2].red;
    dim2 += 8;
    k48 = ker[0 + dim2].red;
    k49 = ker[1 + dim2].red;
    k50 = ker[2 + dim2].red;
    k51 = ker[3 + dim2].red;
    k52 = ker[4 + dim2].red;
    k53 = ker[5 + dim2].red;
    k54 = ker[6 + dim2].red;
    k55 = ker[7 + dim2].red;
    dim2 += 8;
    k56 = ker[0 + dim2].red;
    k57 = ker[1 + dim2].red;
    k58 = ker[2 + dim2].red;
    k59 = ker[3 + dim2].red;
    k60 = ker[4 + dim2].red;
    k61 = ker[5 + dim2].red;
    k62 = ker[6 + dim2].red;
    k63 = ker[7 + dim2].red;

    l0 = ker[0].green;
    l1 = ker[1].green;
    l2 = ker[2].green;
    l3 = ker[3].green;
    l4 = ker[4].green;
    l5 = ker[5].green;
    l6 = ker[6].green;
    l7 = ker[7].green;
    l8 = ker[0 + dim3].green;
    l9 = ker[1 + dim3].green;
    l10 = ker[2 + dim3].green;
    l11 = ker[3 + dim3].green;
    l12 = ker[4 + dim3].green;
    l13 = ker[5 + dim3].green;
    l14 = ker[6 + dim3].green;
    l15 = ker[7 + dim3].green;
    dim3 += 8;
    l16 = ker[0 + dim3].green;
    l17 = ker[1 + dim3].green;
    l18 = ker[2 + dim3].green;
    l19 = ker[3 + dim3].green;
    l20 = ker[4 + dim3].green;
    l21 = ker[5 + dim3].green;
    l22 = ker[6 + dim3].green;
    l23 = ker[7 + dim3].green;
    dim3 += 8;
    l24 = ker[0 + dim3].green;
    l25 = ker[1 + dim3].green;
    l26 = ker[2 + dim3].green;
    l27 = ker[3 + dim3].green;
    l28 = ker[4 + dim3].green;
    l29 = ker[5 + dim3].green;
    l30 = ker[6 + dim3].green;
    l31 = ker[7 + dim3].green;
    dim3 += 8;
    l32 = ker[0 + dim3].green;
    l33 = ker[1 + dim3].green;
    l34 = ker[2 + dim3].green;
    l35 = ker[3 + dim3].green;
    l36 = ker[4 + dim3].green;
    l37 = ker[5 + dim3].green;
    l38 = ker[6 + dim3].green;
    l39 = ker[7 + dim3].green;
    dim3 += 8;
    l40 = ker[0 + dim3].green;
    l41 = ker[1 + dim3].green;
    l42 = ker[2 + dim3].green;
    l43 = ker[3 + dim3].green;
    l44 = ker[4 + dim3].green;
    l45 = ker[5 + dim3].green;
    l46 = ker[6 + dim3].green;
    l47 = ker[7 + dim3].green;
    dim3 += 8;
    l48 = ker[0 + dim3].green;
    l49 = ker[1 + dim3].green;
    l50 = ker[2 + dim3].green;
    l51 = ker[3 + dim3].green;
    l52 = ker[4 + dim3].green;
    l53 = ker[5 + dim3].green;
    l54 = ker[6 + dim3].green;
    l55 = ker[7 + dim3].green;
    dim3 += 8;
    l56 = ker[0 + dim3].green;
    l57 = ker[1 + dim3].green;
    l58 = ker[2 + dim3].green;
    l59 = ker[3 + dim3].green;
    l60 = ker[4 + dim3].green;
    l61 = ker[5 + dim3].green;
    l62 = ker[6 + dim3].green;
    l63 = ker[7 + dim3].green;

    m0 = ker[0].blue;
    m1 = ker[1].blue;
    m2 = ker[2].blue;
    m3 = ker[3].blue;
    m4 = ker[4].blue;
    m5 = ker[5].blue;
    m6 = ker[6].blue;
    m7 = ker[7].blue;
    m8 = ker[0 + dim4].blue;
    m9 = ker[1 + dim4].blue;
    m10 = ker[2 + dim4].blue;
    m11 = ker[3 + dim4].blue;
    m12 = ker[4 + dim4].blue;
    m13 = ker[5 + dim4].blue;
    m14 = ker[6 + dim4].blue;
    m15 = ker[7 + dim4].blue;
    dim4 += 8;
    m16 = ker[0 + dim4].blue;
    m17 = ker[1 + dim4].blue;
    m18 = ker[2 + dim4].blue;
    m19 = ker[3 + dim4].blue;
    m20 = ker[4 + dim4].blue;
    m21 = ker[5 + dim4].blue;
    m22 = ker[6 + dim4].blue;
    m23 = ker[7 + dim4].blue;
    dim4 += 8;
    m24 = ker[0 + dim4].blue;
    m25 = ker[1 + dim4].blue;
    m26 = ker[2 + dim4].blue;
    m27 = ker[3 + dim4].blue;
    m28 = ker[4 + dim4].blue;
    m29 = ker[5 + dim4].blue;
    m30 = ker[6 + dim4].blue;
    m31 = ker[7 + dim4].blue;
    dim4 += 8;
    m32 = ker[0 + dim4].blue;
    m33 = ker[1 + dim4].blue;
    m34 = ker[2 + dim4].blue;
    m35 = ker[3 + dim4].blue;
    m36 = ker[4 + dim4].blue;
    m37 = ker[5 + dim4].blue;
    m38 = ker[6 + dim4].blue;
    m39 = ker[7 + dim4].blue;
    dim4 += 8;
    m40 = ker[0 + dim4].blue;
    m41 = ker[1 + dim4].blue;
    m42 = ker[2 + dim4].blue;
    m43 = ker[3 + dim4].blue;
    m44 = ker[4 + dim4].blue;
    m45 = ker[5 + dim4].blue;
    m46 = ker[6 + dim4].blue;
    m47 = ker[7 + dim4].blue;
    dim4 += 8;
    m48 = ker[0 + dim4].blue;
    m49 = ker[1 + dim4].blue;
    m50 = ker[2 + dim4].blue;
    m51 = ker[3 + dim4].blue;
    m52 = ker[4 + dim4].blue;
    m53 = ker[5 + dim4].blue;
    m54 = ker[6 + dim4].blue;
    m55 = ker[7 + dim4].blue;
    dim4 += 8;
    m56 = ker[0 + dim4].blue;
    m57 = ker[1 + dim4].blue;
    m58 = ker[2 + dim4].blue;
    m59 = ker[3 + dim4].blue;
    m60 = ker[4 + dim4].blue;
    m61 = ker[5 + dim4].blue;
    m62 = ker[6 + dim4].blue;
    m63 = ker[7 + dim4].blue;

    for (i = 0; i < num; i++)
    {
        src3 = src;
        dst2 = dst;
        
        for (j = 0; j < num; j++)
        {
            src2 = src3;
            acc = 0;

            acc += src2[0].red * k0;
            acc += src2[0].green * l0;
            acc += src2[0].blue * m0;

            acc += src2[1].red * k1;
            acc += src2[1].green * l1;
            acc += src2[1].blue * m1;

            acc += src2[2].red * k2;
            acc += src2[2].green * l2;
            acc += src2[2].blue * m2;

            acc += src2[3].red * k3;
            acc += src2[3].green * l3;
            acc += src2[3].blue * m3;

            acc += src2[4].red * k4;
            acc += src2[4].green * l4;
            acc += src2[4].blue * m4;

            acc += src2[5].red * k5;
            acc += src2[5].green * l5;
            acc += src2[5].blue * m5;

            acc += src2[6].red * k6;
            acc += src2[6].green * l6;
            acc += src2[6].blue * m6;

            acc += src2[7].red * k7;
            acc += src2[7].green * l7;
            acc += src2[7].blue * m7;

            src2 += dim;

            acc += src2[0].red * k8;
            acc += src2[0].green * l8;
            acc += src2[0].blue * m8;

            acc += src2[1].red * k9;
            acc += src2[1].green * l9;
            acc += src2[1].blue * m9;

            acc += src2[2].red * k10;
            acc += src2[2].green * l10;
            acc += src2[2].blue * m10;

            acc += src2[3].red * k11;
            acc += src2[3].green * l11;
            acc += src2[3].blue * m11;

            acc += src2[4].red * k12;
            acc += src2[4].green * l12;
            acc += src2[4].blue * m12;

            acc += src2[5].red * k13;
            acc += src2[5].green * l13;
            acc += src2[5].blue * m13;

            acc += src2[6].red * k14;
            acc += src2[6].green * l14;
            acc += src2[6].blue * m14;

            acc += src2[7].red * k15;
            acc += src2[7].green * l15;
            acc += src2[7].blue * m15;

            src2 += dim;

            acc += src2[0].red * k16;
            acc += src2[0].green * l16;
            acc += src2[0].blue * m16;

            acc += src2[1].red * k17;
            acc += src2[1].green * l17;
            acc += src2[1].blue * m17;

            acc += src2[2].red * k18;
            acc += src2[2].green * l18;
            acc += src2[2].blue * m18;

            acc += src2[3].red * k19;
            acc += src2[3].green * l19;
            acc += src2[3].blue * m19;

            acc += src2[4].red * k20;
            acc += src2[4].green * l20;
            acc += src2[4].blue * m20;

            acc += src2[5].red * k21;
            acc += src2[5].green * l21;
            acc += src2[5].blue * m21;

            acc += src2[6].red * k22;
            acc += src2[6].green * l22;
            acc += src2[6].blue * m22;

            acc += src2[7].red * k23;
            acc += src2[7].green * l23;
            acc += src2[7].blue * m23;

            src2 += dim;

            acc += src2[0].red * k24;
            acc += src2[0].green * l24;
            acc += src2[0].blue * m24;

            acc += src2[1].red * k25;
            acc += src2[1].green * l25;
            acc += src2[1].blue * m25;

            acc += src2[2].red * k26;
            acc += src2[2].green * l26;
            acc += src2[2].blue * m26;

            acc += src2[3].red * k27;
            acc += src2[3].green * l27;
            acc += src2[3].blue * m27;

            acc += src2[4].red * k28;
            acc += src2[4].green * l28;
            acc += src2[4].blue * m28;

            acc += src2[5].red * k29;
            acc += src2[5].green * l29;
            acc += src2[5].blue * m29;

            acc += src2[6].red * k30;
            acc += src2[6].green * l30;
            acc += src2[6].blue * m30;

            acc += src2[7].red * k31;
            acc += src2[7].green * l31;
            acc += src2[7].blue * m31;

            src2 += dim;

            acc += src2[0].red * k32;
            acc += src2[0].green * l32;
            acc += src2[0].blue * m32;

            acc += src2[1].red * k33;
            acc += src2[1].green * l33;
            acc += src2[1].blue * m33;

            acc += src2[2].red * k34;
            acc += src2[2].green * l34;
            acc += src2[2].blue * m34;

            acc += src2[3].red * k35;
            acc += src2[3].green * l35;
            acc += src2[3].blue * m35;

            acc += src2[4].red * k36;
            acc += src2[4].green * l36;
            acc += src2[4].blue * m36;

            acc += src2[5].red * k37;
            acc += src2[5].green * l37;
            acc += src2[5].blue * m37;

            acc += src2[6].red * k38;
            acc += src2[6].green * l38;
            acc += src2[6].blue * m38;

            acc += src2[7].red * k39;
            acc += src2[7].green * l39;
            acc += src2[7].blue * m39;

            src2 += dim;

            acc += src2[0].red * k40;
            acc += src2[0].green * l40;
            acc += src2[0].blue * m40;

            acc += src2[1].red * k41;
            acc += src2[1].green * l41;
            acc += src2[1].blue * m41;

            acc += src2[2].red * k42;
            acc += src2[2].green * l42;
            acc += src2[2].blue * m42;

            acc += src2[3].red * k43;
            acc += src2[3].green * l43;
            acc += src2[3].blue * m43;

            acc += src2[4].red * k44;
            acc += src2[4].green * l44;
            acc += src2[4].blue * m44;

            acc += src2[5].red * k45;
            acc += src2[5].green * l45;
            acc += src2[5].blue * m45;

            acc += src2[6].red * k46;
            acc += src2[6].green * l46;
            acc += src2[6].blue * m46;

            acc += src2[7].red * k47;
            acc += src2[7].green * l47;
            acc += src2[7].blue * m47;

            src2 += dim;

            acc += src2[0].red * k48;
            acc += src2[0].green * l48;
            acc += src2[0].blue * m48;

            acc += src2[1].red * k49;
            acc += src2[1].green * l49;
            acc += src2[1].blue * m49;

            acc += src2[2].red * k50;
            acc += src2[2].green * l50;
            acc += src2[2].blue * m50;

            acc += src2[3].red * k51;
            acc += src2[3].green * l51;
            acc += src2[3].blue * m51;

            acc += src2[4].red * k52;
            acc += src2[4].green * l52;
            acc += src2[4].blue * m52;

            acc += src2[5].red * k53;
            acc += src2[5].green * l53;
            acc += src2[5].blue * m53;

            acc += src2[6].red * k54;
            acc += src2[6].green * l54;
            acc += src2[6].blue * m54;

            acc += src2[7].red * k55;
            acc += src2[7].green * l55;
            acc += src2[7].blue * m55;

            src2 += dim;

            acc += src2[0].red * k56;
            acc += src2[0].green * l56;
            acc += src2[0].blue * m56;

            acc += src2[1].red * k57;
            acc += src2[1].green * l57;
            acc += src2[1].blue * m57;

            acc += src2[2].red * k58;
            acc += src2[2].green * l58;
            acc += src2[2].blue * m58;

            acc += src2[3].red * k59;
            acc += src2[3].green * l59;
            acc += src2[3].blue * m59;

            acc += src2[4].red * k60;
            acc += src2[4].green * l60;
            acc += src2[4].blue * m60;

            acc += src2[5].red * k61;
            acc += src2[5].green * l61;
            acc += src2[5].blue * m61;

            acc += src2[6].red * k62;
            acc += src2[6].green * l62;
            acc += src2[6].blue * m62;

            acc += src2[7].red * k63;
            acc += src2[7].green * l63;
            acc += src2[7].blue * m63;

            ++src3;
            *dst2++ = acc;
        }
        src += dim;
        dst += dim;
    }
}

/*********************************************************************
 * register_conv_functions - Register all of your different versions
 *     of the convolution functions  with the driver by calling the
 *     add_conv_function() for each test function. When you run the
 *     driver program, it will test and report the performance of each
 *     registered test function.
 *********************************************************************/

void register_conv_functions()
{
    add_conv_function(&naive_conv, naive_conv_descr);
    add_conv_function(&convolution, convolution_descr);
    /* ... Register additional test functions here */
}

/************************
 * AVERAGE POOLING KERNEL
 ************************/

/*********************************************************
 * Your different versions of the average pooling  go here
 *********************************************************/

/*
 * naive_average_pooling - The naive baseline version of average pooling
 */
char naive_average_pooling_descr[] = "Naive Average Pooling: Naive baseline implementation";
void naive_average_pooling(int dim, pixel *src2, pixel *dst)
{
    int i, j, k, l;

    for (i = 0; i < dim / 2; i++)
        for (j = 0; j < dim / 2; j++)
        {
            dst[RIDX(i, j, dim / 2)].red = 0;
            dst[RIDX(i, j, dim / 2)].green = 0;
            dst[RIDX(i, j, dim / 2)].blue = 0;
            for (k = 0; k < 2; k++)
            {
                for (l = 0; l < 2; l++)
                {
                    dst[RIDX(i, j, dim / 2)].red += src2[RIDX(i*2 + k, j*2 + l, dim)].red;
                    dst[RIDX(i, j, dim / 2)].green += src2[RIDX(i*2 + k, j*2 + l, dim)].green;
                    dst[RIDX(i, j, dim / 2)].blue += src2[RIDX(i*2 + k, j*2 + l, dim)].blue;
                }
            }
            dst[RIDX(i, j, dim / 2)].red /= 4;
            dst[RIDX(i, j, dim / 2)].green /= 4;
            dst[RIDX(i, j, dim / 2)].blue /= 4;
        }
}

/*
 * average_pooling - Your current working version of average_pooling
 * IMPORTANT: This is the version you will be graded on
 */
char average_pooling_descr[] = "Average Pooling: Current working version";
void average_pooling(int dim, pixel *src2, pixel *dst)
{

     int i, j;

int src=src2;

     for (i = 0; i < dim>>1; i+=2)
        for (j = 0; j < dim>>1; j+=8)
        {
     
     
                    dst[RIDX(i,j,dim>>1)].red = (src2[RIDX(i<<1,j<<1,dim)].red +src2[RIDX(i<<1,j<<1,dim)+1].red+src2[RIDX(i<<1,j<<1,dim)+dim].red+src2[RIDX(i<<1,j<<1,dim)+dim+1].red)>>2;
                    dst[RIDX(i,j,dim>>1)].green += (src2[RIDX(i<<1,j<<1,dim)].green+ src2[RIDX(i<<1,j<<1,dim)+1].green+src2[RIDX(i<<1,j<<1,dim)+dim].green+src2[RIDX(i<<1,j<<1,dim)+dim+1].green)>>2;
                    dst[RIDX(i,j,dim>>1)].blue += (src2[RIDX(i<<1,j<<1,dim)].blue+ 
                    src2[RIDX(i<<1,j<<1,dim)+1].blue+
                     src2[RIDX(i<<1,j<<1,dim)+dim].blue+
                    src2[RIDX(i<<1,j<<1,dim)+dim+1].blue)>>2;
                                
                                  
  
                    
                    src2+=2;
                         dst[RIDX(i,j+1,dim>>1)].red = (src2[RIDX(i<<1,j<<1,dim)].red +src2[RIDX(i<<1,j<<1,dim)+1].red+src2[RIDX(i<<1,j<<1,dim)+dim].red+src2[RIDX(i<<1,j<<1,dim)+dim+1].red)>>2;
                    dst[RIDX(i,j+1,dim>>1)].green += (src2[RIDX(i<<1,j<<1,dim)].green+ src2[RIDX(i<<1,j<<1,dim)+1].green+src2[RIDX(i<<1,j<<1,dim)+dim].green+src2[RIDX(i<<1,j<<1,dim)+dim+1].green)>>2;
                    dst[RIDX(i,j+1,dim>>1)].blue += (src2[RIDX(i<<1,j<<1,dim)].blue+ 
                    src2[RIDX(i<<1,j<<1,dim)+1].blue+
                     src2[RIDX(i<<1,j<<1,dim)+dim].blue+
                    src2[RIDX(i<<1,j<<1,dim)+dim+1].blue)>>2;
                    
                    
                    src2+=2;
                         dst[RIDX(i,j+2,dim>>1)].red = (src2[RIDX(i<<1,j<<1,dim)].red +src2[RIDX(i<<1,j<<1,dim)+1].red+src2[RIDX(i<<1,j<<1,dim)+dim].red+src2[RIDX(i<<1,j<<1,dim)+dim+1].red)>>2;
                    dst[RIDX(i,j+2,dim>>1)].green += (src2[RIDX(i<<1,j<<1,dim)].green+ src2[RIDX(i<<1,j<<1,dim)+1].green+src2[RIDX(i<<1,j<<1,dim)+dim].green+src2[RIDX(i<<1,j<<1,dim)+dim+1].green)>>2;
                    dst[RIDX(i,j+2,dim>>1)].blue += (src2[RIDX(i<<1,j<<1,dim)].blue+ 
                    src2[RIDX(i<<1,j<<1,dim)+1].blue+
                     src2[RIDX(i<<1,j<<1,dim)+dim].blue+
                    src2[RIDX(i<<1,j<<1,dim)+dim+1].blue)>>2;
                    
                    src2+=2;
                         dst[RIDX(i,j+3,dim>>1)].red = (src2[RIDX(i<<1,j<<1,dim)].red +src2[RIDX(i<<1,j<<1,dim)+1].red+src2[RIDX(i<<1,j<<1,dim)+dim].red+src2[RIDX(i<<1,j<<1,dim)+dim+1].red)>>2;
                    dst[RIDX(i,j+3,dim>>1)].green += (src2[RIDX(i<<1,j<<1,dim)].green+ src2[RIDX(i<<1,j<<1,dim)+1].green+src2[RIDX(i<<1,j<<1,dim)+dim].green+src2[RIDX(i<<1,j<<1,dim)+dim+1].green)>>2;
                    dst[RIDX(i,j+3,dim>>1)].blue += (src2[RIDX(i<<1,j<<1,dim)].blue+ 
                    src2[RIDX(i<<1,j<<1,dim)+1].blue+
                     src2[RIDX(i<<1,j<<1,dim)+dim].blue+
                    src2[RIDX(i<<1,j<<1,dim)+dim+1].blue)>>2;
       
                    src2+=2;
                         dst[RIDX(i,j+4,dim>>1)].red = (src2[RIDX(i<<1,j<<1,dim)].red +src2[RIDX(i<<1,j<<1,dim)+1].red+src2[RIDX(i<<1,j<<1,dim)+dim].red+src2[RIDX(i<<1,j<<1,dim)+dim+1].red)>>2;
                    dst[RIDX(i,j+4,dim>>1)].green += (src2[RIDX(i<<1,j<<1,dim)].green+ src2[RIDX(i<<1,j<<1,dim)+1].green+src2[RIDX(i<<1,j<<1,dim)+dim].green+src2[RIDX(i<<1,j<<1,dim)+dim+1].green)>>2;
                    dst[RIDX(i,j+4,dim>>1)].blue += (src2[RIDX(i<<1,j<<1,dim)].blue+ 
                    src2[RIDX(i<<1,j<<1,dim)+1].blue+
                     src2[RIDX(i<<1,j<<1,dim)+dim].blue+
                    src2[RIDX(i<<1,j<<1,dim)+dim+1].blue)>>2;
                    
                    src2+=2;
                         dst[RIDX(i,j+5,dim>>1)].red = (src2[RIDX(i<<1,j<<1,dim)].red +src2[RIDX(i<<1,j<<1,dim)+1].red+src2[RIDX(i<<1,j<<1,dim)+dim].red+src2[RIDX(i<<1,j<<1,dim)+dim+1].red)>>2;
                    dst[RIDX(i,j+5,dim>>1)].green += (src2[RIDX(i<<1,j<<1,dim)].green+ src2[RIDX(i<<1,j<<1,dim)+1].green+src2[RIDX(i<<1,j<<1,dim)+dim].green+src2[RIDX(i<<1,j<<1,dim)+dim+1].green)>>2;
                    dst[RIDX(i,j+5,dim>>1)].blue += (src2[RIDX(i<<1,j<<1,dim)].blue+ 
                    src2[RIDX(i<<1,j<<1,dim)+1].blue+
                     src2[RIDX(i<<1,j<<1,dim)+dim].blue+
                    src2[RIDX(i<<1,j<<1,dim)+dim+1].blue)>>2;
                    
                    src2+=2;
                         dst[RIDX(i,j+6,dim>>1)].red = (src2[RIDX(i<<1,j<<1,dim)].red +src2[RIDX(i<<1,j<<1,dim)+1].red+src2[RIDX(i<<1,j<<1,dim)+dim].red+src2[RIDX(i<<1,j<<1,dim)+dim+1].red)>>2;
                    dst[RIDX(i,j+6,dim>>1)].green += (src2[RIDX(i<<1,j<<1,dim)].green+ src2[RIDX(i<<1,j<<1,dim)+1].green+src2[RIDX(i<<1,j<<1,dim)+dim].green+src2[RIDX(i<<1,j<<1,dim)+dim+1].green)>>2;
                    dst[RIDX(i,j+6,dim>>1)].blue += (src2[RIDX(i<<1,j<<1,dim)].blue+ 
                    src2[RIDX(i<<1,j<<1,dim)+1].blue+
                     src2[RIDX(i<<1,j<<1,dim)+dim].blue+
                    src2[RIDX(i<<1,j<<1,dim)+dim+1].blue)>>2;
                    
                    src2+=2;
                         dst[RIDX(i,j+7,dim>>1)].red = (src2[RIDX(i<<1,j<<1,dim)].red +src2[RIDX(i<<1,j<<1,dim)+1].red+src2[RIDX(i<<1,j<<1,dim)+dim].red+src2[RIDX(i<<1,j<<1,dim)+dim+1].red)>>2;
                    dst[RIDX(i,j+7,dim>>1)].green += (src2[RIDX(i<<1,j<<1,dim)].green+ src2[RIDX(i<<1,j<<1,dim)+1].green+src2[RIDX(i<<1,j<<1,dim)+dim].green+src2[RIDX(i<<1,j<<1,dim)+dim+1].green)>>2;
                    dst[RIDX(i,j+7,dim>>1)].blue += (src2[RIDX(i<<1,j<<1,dim)].blue+ 
                    src2[RIDX(i<<1,j<<1,dim)+1].blue+
                     src2[RIDX(i<<1,j<<1,dim)+dim].blue+
                    src2[RIDX(i<<1,j<<1,dim)+dim+1].blue)>>2;
         
         
                    
                    src2+=dim*2-14;
                                 dst[RIDX(i+1,j,dim>>1)].red = (src2[RIDX(i<<1,j<<1,dim)].red +src2[RIDX(i<<1,j<<1,dim)+1].red+src2[RIDX(i<<1,j<<1,dim)+dim].red+src2[RIDX(i<<1,j<<1,dim)+dim+1].red)>>2;
                    dst[RIDX(i+1,j,dim>>1)].green += (src2[RIDX(i<<1,j<<1,dim)].green+ src2[RIDX(i<<1,j<<1,dim)+1].green+src2[RIDX(i<<1,j<<1,dim)+dim].green+src2[RIDX(i<<1,j<<1,dim)+dim+1].green)>>2;
                    dst[RIDX(i+1,j,dim>>1)].blue += (src2[RIDX(i<<1,j<<1,dim)].blue+ 
                    src2[RIDX(i<<1,j<<1,dim)+1].blue+
                     src2[RIDX(i<<1,j<<1,dim)+dim].blue+
                    src2[RIDX(i<<1,j<<1,dim)+dim+1].blue)>>2;
                    
                    src2+=2;
                                     dst[RIDX(i+1,j+1,dim>>1)].red = (src2[RIDX(i<<1,j<<1,dim)].red +src2[RIDX(i<<1,j<<1,dim)+1].red+src2[RIDX(i<<1,j<<1,dim)+dim].red+src2[RIDX(i<<1,j<<1,dim)+dim+1].red)>>2;
                    dst[RIDX(i+1,j+1,dim>>1)].green += (src2[RIDX(i<<1,j<<1,dim)].green+ src2[RIDX(i<<1,j<<1,dim)+1].green+src2[RIDX(i<<1,j<<1,dim)+dim].green+src2[RIDX(i<<1,j<<1,dim)+dim+1].green)>>2;
                    dst[RIDX(i+1,j+1,dim>>1)].blue += (src2[RIDX(i<<1,j<<1,dim)].blue+ 
                    src2[RIDX(i<<1,j<<1,dim)+1].blue+
                     src2[RIDX(i<<1,j<<1,dim)+dim].blue+
                    src2[RIDX(i<<1,j<<1,dim)+dim+1].blue)>>2;
          
               src2+=2;
                                     dst[RIDX(i+1,j+2,dim>>1)].red = (src2[RIDX(i<<1,j<<1,dim)].red +src2[RIDX(i<<1,j<<1,dim)+1].red+src2[RIDX(i<<1,j<<1,dim)+dim].red+src2[RIDX(i<<1,j<<1,dim)+dim+1].red)>>2;
                    dst[RIDX(i+1,j+2,dim>>1)].green += (src2[RIDX(i<<1,j<<1,dim)].green+ src2[RIDX(i<<1,j<<1,dim)+1].green+src2[RIDX(i<<1,j<<1,dim)+dim].green+src2[RIDX(i<<1,j<<1,dim)+dim+1].green)>>2;
                    dst[RIDX(i+1,j+2,dim>>1)].blue += (src2[RIDX(i<<1,j<<1,dim)].blue+ 
                    src2[RIDX(i<<1,j<<1,dim)+1].blue+
                     src2[RIDX(i<<1,j<<1,dim)+dim].blue+
                    src2[RIDX(i<<1,j<<1,dim)+dim+1].blue)>>2;
                    
                         src2+=2;
                                     dst[RIDX(i+1,j+3,dim>>1)].red = (src2[RIDX(i<<1,j<<1,dim)].red +src2[RIDX(i<<1,j<<1,dim)+1].red+src2[RIDX(i<<1,j<<1,dim)+dim].red+src2[RIDX(i<<1,j<<1,dim)+dim+1].red)>>2;
                    dst[RIDX(i+1,j+3,dim>>1)].green += (src2[RIDX(i<<1,j<<1,dim)].green+ src2[RIDX(i<<1,j<<1,dim)+1].green+src2[RIDX(i<<1,j<<1,dim)+dim].green+src2[RIDX(i<<1,j<<1,dim)+dim+1].green)>>2;
                    dst[RIDX(i+1,j+3,dim>>1)].blue += (src2[RIDX(i<<1,j<<1,dim)].blue+ 
                    src2[RIDX(i<<1,j<<1,dim)+1].blue+
                     src2[RIDX(i<<1,j<<1,dim)+dim].blue+
                    src2[RIDX(i<<1,j<<1,dim)+dim+1].blue)>>2;
                    
                          src2+=2;
                                     dst[RIDX(i+1,j+4,dim>>1)].red = (src2[RIDX(i<<1,j<<1,dim)].red +src2[RIDX(i<<1,j<<1,dim)+1].red+src2[RIDX(i<<1,j<<1,dim)+dim].red+src2[RIDX(i<<1,j<<1,dim)+dim+1].red)>>2;
                    dst[RIDX(i+1,j+4,dim>>1)].green += (src2[RIDX(i<<1,j<<1,dim)].green+ src2[RIDX(i<<1,j<<1,dim)+1].green+src2[RIDX(i<<1,j<<1,dim)+dim].green+src2[RIDX(i<<1,j<<1,dim)+dim+1].green)>>2;
                    dst[RIDX(i+1,j+4,dim>>1)].blue += (src2[RIDX(i<<1,j<<1,dim)].blue+ 
                    src2[RIDX(i<<1,j<<1,dim)+1].blue+
                     src2[RIDX(i<<1,j<<1,dim)+dim].blue+
                    src2[RIDX(i<<1,j<<1,dim)+dim+1].blue)>>2;
                    
                          src2+=2;
                                     dst[RIDX(i+1,j+5,dim>>1)].red = (src2[RIDX(i<<1,j<<1,dim)].red +src2[RIDX(i<<1,j<<1,dim)+1].red+src2[RIDX(i<<1,j<<1,dim)+dim].red+src2[RIDX(i<<1,j<<1,dim)+dim+1].red)>>2;
                    dst[RIDX(i+1,j+5,dim>>1)].green += (src2[RIDX(i<<1,j<<1,dim)].green+ src2[RIDX(i<<1,j<<1,dim)+1].green+src2[RIDX(i<<1,j<<1,dim)+dim].green+src2[RIDX(i<<1,j<<1,dim)+dim+1].green)>>2;
                    dst[RIDX(i+1,j+5,dim>>1)].blue += (src2[RIDX(i<<1,j<<1,dim)].blue+ 
                    src2[RIDX(i<<1,j<<1,dim)+1].blue+
                     src2[RIDX(i<<1,j<<1,dim)+dim].blue+
                    src2[RIDX(i<<1,j<<1,dim)+dim+1].blue)>>2;
                    
                          src2+=2;
                                     dst[RIDX(i+1,j+6,dim>>1)].red = (src2[RIDX(i<<1,j<<1,dim)].red +src2[RIDX(i<<1,j<<1,dim)+1].red+src2[RIDX(i<<1,j<<1,dim)+dim].red+src2[RIDX(i<<1,j<<1,dim)+dim+1].red)>>2;
                    dst[RIDX(i+1,j+6,dim>>1)].green += (src2[RIDX(i<<1,j<<1,dim)].green+ src2[RIDX(i<<1,j<<1,dim)+1].green+src2[RIDX(i<<1,j<<1,dim)+dim].green+src2[RIDX(i<<1,j<<1,dim)+dim+1].green)>>2;
                    dst[RIDX(i+1,j+6,dim>>1)].blue += (src2[RIDX(i<<1,j<<1,dim)].blue+ 
                    src2[RIDX(i<<1,j<<1,dim)+1].blue+
                     src2[RIDX(i<<1,j<<1,dim)+dim].blue+
                    src2[RIDX(i<<1,j<<1,dim)+dim+1].blue)>>2;
                    
                          src2+=2;
                                     dst[RIDX(i+1,j+7,dim>>1)].red = (src2[RIDX(i<<1,j<<1,dim)].red +src2[RIDX(i<<1,j<<1,dim)+1].red+src2[RIDX(i<<1,j<<1,dim)+dim].red+src2[RIDX(i<<1,j<<1,dim)+dim+1].red)>>2;
                    dst[RIDX(i+1,j+7,dim>>1)].green += (src2[RIDX(i<<1,j<<1,dim)].green+ src2[RIDX(i<<1,j<<1,dim)+1].green+src2[RIDX(i<<1,j<<1,dim)+dim].green+src2[RIDX(i<<1,j<<1,dim)+dim+1].green)>>2;
                    dst[RIDX(i+1,j+7,dim>>1)].blue += (src2[RIDX(i<<1,j<<1,dim)].blue+ 
                    src2[RIDX(i<<1,j<<1,dim)+1].blue+
                     src2[RIDX(i<<1,j<<1,dim)+dim].blue+
                    src2[RIDX(i<<1,j<<1,dim)+dim+1].blue)>>2;
                    
                    
                    
                    
                    
          
src2=src;
        }
}

/******************************************************************************
 * register_average_pooling_functions - Register all of your different versions
 *     of the average pooling  with the driver by calling the
 *     add_average_pooling_function() for each test function. When you run the
 *     driver program, it will test and report the performance of each
 *     registered test function.
 ******************************************************************************/

void register_average_pooling_functions()
{
    add_average_pooling_function(&naive_average_pooling, naive_average_pooling_descr);
    add_average_pooling_function(&average_pooling, average_pooling_descr);
    /* ... Register additional test functions here */
}
