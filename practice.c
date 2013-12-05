/***
 *  DETAILED EXPLANATION OF THE ALGORITHM
 *  For simplicity, this example will use very small numbers and systematic float values.
 *
 *  fMRI takes 3D images of the brain, recording a float value at each point in the brain.
 *  This can be repeated many times over a time series, say 1 image of the brain per second.
 *
 *  Suppose that a simple fMRI image records 8 places in the brain (think of a 2*2*2 3D matrix).
 *
 *  Suppose for the first fMRI image (time = 0 ), we have the following float values: 
 *  (note that it is hard to represent 3D in this 2dimensional document, so this 
 *   is done by representing 2D matrices at z = 0 and z = 1)
 *  Time = 0:
 *  
 *     z = 0:                  z = 1:                (i.e. image[0][1][1] = 6.00, for example)
 *     _   y0       y1 _        _   y0       y1 _
 *  x0|   1.00 |   2.00 |    x0|   5.00 |   6.00 |
 *  x1|_  3.00 |   4.00_|    x1|_  7.00 |   8.00_|
 *
 *    We can also represent this image data as one long vector at time t = 0:
 *    image(0) = {0.00, 1.00, 2.00. 3.00, 4.00, 5.00, 6.00, 7.00}
 *
 *  Now, suppose that we are taking images at 5 different times:  time = 0, 1, 2, 3, 4 and we have:
 *    image(0) = {  1.00,   2.00,   3.00.   4.00,   5.00,   6.00,   7.00,   8.00}
 *    image(1) = {101.00, 102.00, 103.00. 104.00, 105.00, 106.00, 107.00, 108.00}
 *    image(2) = {201.00, 202.00, 203.00. 204.00, 205.00, 206.00, 207.00, 208.00}
 *    image(3) = {301.00, 302.00, 303.00. 304.00, 305.00, 306.00, 307.00, 308.00}
 *    image(4) = {401.00, 402.00, 403.00. 404.00, 405.00, 406.00, 407.00, 408.00}
 *
 * The challenge is, that this 'raw' fMRI data contains 'noise' created by the
 * patient's heartbeat, head motion, breathing, etc, and we need to clean
 * out that noise to get a more accurate fMRI scan.  
 *
 * To do this, a 'covariate file' is produced containing calculate estimates for
 * how much each element such as heartbeat, etc, has influenced the fMRI values  
 * at a given time stamp.  
 *
 *  Suppose we have a covariate file as follows: 
 *                   heartRate     head Motion    Respiration
 *  for time(0) = {    1.00,          2.00,          3.00   }
 *  for time(1) = {   11.00,         12.00,         13.00   }
 *  for time(2) = {   21.00,         22.00,         23.00   }
 *  for time(3) = {   31.00,         32.00,         33.00   }
 *  for time(4) = {   41.00,         42.00,         43.00   }
 * 
 *
 *  Now consider a position of the brain that is measured, say image[0][1][1]. 
 *  Recall that this position is imaged 5 times in our case.  We want to make
 *  a vector for position image[0][1][1] over time, so we have a vector of 
 *  all the values for image[0][1][1] at each of the time stamps, and we need
 *  to do this for each position:  
 *                                                              time0   time1   time2   time3   time4
 *  position [0][0][0]  (or [0] in the long vecotor) we have:    1.00  101.00  201.00  301.00  401.00
 *  position [0][1][0]  (or [1] in the long vecotor) we have:    2.00  102.00  202.00  302.00  402.00
 *  position [1][0][0]  (or [2] in the long vecotor) we have:    3.00  103.00  203.00  303.00  403.00
 *  position [1][1][0]  (or [3] in the long vecotor) we have:    4.00  104.00  204.00  304.00  404.00
 *  position [0][0][1]  (or [4] in the long vecotor) we have:    5.00  105.00  205.00  305.00  405.00
 *  position [0][1][1]  (or [5] in the long vecotor) we have:    6.00  106.00  206.00  306.00  406.00
 *  position [1][0][1]  (or [6] in the long vecotor) we have:    7.00  107.00  207.00  307.00  407.00
 *  position [1][1][1]  (or [7] in the long vecotor) we have:    8.00  108.00  208.00  308.00  408.00
 *
 *
 *  Once we have these vectors, we can use them with the covariate data to find an estimate for
 *  how much the extraneous elements such as heart rate, etc, influenced the value at each point
 *  in the brain over time.  We want to solve:
 *
 *         Y = b1*X1 + b2*X2 + b3*X3 + U   // where b1, b2, b3 are scalars
 *
 *  Where Yt is a vector over time for a given point (such as position [0][0][0] above:  {1.00 101.00... 401.00})
 *  X1 is the vector of covariates over time for heart rate:   {1.00, 11.00, 21.00, 31.00, 41.00}
 *  X2 is the vector of covariates over time for head motion:  {2.00, 12.00, 22.00, 32.00, 42.00}  
 *  X3 is the vector of covariates over time for respiration:  {3.00, 13.00, 23.00, 33.00, 43.00}
 *
 *           consider example for position [0][1][1] (5 in long vector)
 * 
 *                   Y         b1   x1        b2     X2        b3     X3     U
 *                  6.00  =  (b1*  1.00)  +  (b2 *  2.00)  +  (b3 *  3.00) + u0
 *                106.00  =  (b1* 11.00)  +  (b2 * 12.00)  +  (b3 * 13.00) + u1
 *                206.00  =  (b1* 21.00)  +  (b2 * 22.00)  +  (b3 * 23.00) + u2
 *                306.00  =  (b1* 31.00)  +  (b2 * 32.00)  +  (b3 * 33.00) + u3
 *                406.00  =  (b1* 41.00)  +  (b2 * 42.00)  +  (b3 * 43.00) + u4
 *                
 *
 *   We need a way to estimate the vectors B1, B2, and B3, so that we can figure out what the remaining
 *  'U' vector data is, and this new vector will be our 'cleaned' data. 
 *
 *  To solve/estimate the Beta b1, b2, b3 values, we will use the Ordinary Least Squares (OLS) method:
 *
 *              B = (X^T * X)^-1  *  X^T * y
 *
 *       B     _          X^T                            X        _^-1   _           X^T                    Y  _ 
 *      b1    |1.00 11.00 21.00 31.00 41.00     1.00   2.00   3.00 |    |1.00 11.00 21.00 31.00 41.00      6.00 |
 *      b2    |2.00 12.00 22.00 32.00 42.00    11.00  12.00  13.00 |    |2.00 12.00 22.00 32.00 42.00  * 106.00 |
 *      b3  = |3.00 13.00 23.00 33.00 43.00 *  21.00  22.00  23.00 |  * |3.00 13.00 23.00 33.00 43.00    206.00 |
 *            |                                31.00  32.00  33.00 |    |                                306.00 |
 *            |_                               41.00  42.00  43.00_|    |_                               406.00_|
 *
 *
 *   Once we have B = {b1, b2, b3}, we can turn around and solve
 *
 *                   U = Y - b1X1 + b2X2 + b3X3
 *
 *     U         Y         b1   x1         b2     X2        b3     X3
 *     u0       6.00      (b1*  1.00)  +  (b2 *  2.00)  +  (b3 *  3.00)
 *     u1     106.00      (b1* 11.00)  +  (b2 * 12.00)  +  (b3 * 13.00) 
 *     u2  =  206.00   -  (b1* 21.00)  +  (b2 * 22.00)  +  (b3 * 23.00)
 *     u3     306.00      (b1* 31.00)  +  (b2 * 32.00)  +  (b3 * 33.00)
 *     u4     406.00      (b1* 41.00)  +  (b2 * 42.00)  +  (b3 * 43.00)
 * 
 *  These new cleaned U values will need to be written back into new 3D image matrices.  
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// for simplicity, we will use a very small 3D brain image matrix
#define x 2
#define y 2
#define z 2
// time Count is the number of times we repeatedly took an image of the brain
// in this case we will have brain images over 5 points in time
#define timeCount 5
// number of Covariates is the number of extra elements we are using that
// influence the value of the fMRI image:  In this case, we will 
// use 3 covariates (think of them as "hearRate, headMotion, respiration")
#define numCovariates 3

//pre-declare methods
void fillDataStructures(float** dirtyBrainFiles, float** cleanBrainFiles, float** covariates, float** brainLocations, int print);
void createTranspose(float** covariates, float** covariatesTranspose, int print);
void multiplyCovariatesAndTranspose(float** covariatesTranspose, float** covariates, float** covariatesXtranspose, int print);
// methods for finding matrix inverse
float determinant(float** inputMatrix, int dimension);
void cofactor(float** inputMatrix, int dimension);
void transpose(float** inputMatrix, int dimension);

/****************************************************************************
 *
 *             MAIN
 *
 ***************************************************************************/
void main(int argc, char* argv[])
{
    int i, j;
   
    //"dirty" brain files are the raw data coming out of the fMRI
    // For example, at time = 0, we will have 
    // dirtyBrainFiles[0] = {float, float, float, float, float, float, float, float}
    // with a 'float' value for each of the 8 points (2*2*2) of the brain
    float* dirtyBrainFiles[timeCount];
    float* cleanBrainFiles[timeCount]; // holds the cleaned values

    //for each time period an image was taken, we need an array of covariates
    //for example, at time = 0, we need a heartRate, a headMotion, and a respiratory value
    float* covariates[timeCount];

    // this will be the transpose of the covariate matrix
    float* covariatesTranspose[numCovariates];
    //this will hold covariatesTranxpose * covariates - a numCovariates x numCovariates matrix
    float* covariatesXtranspose[numCovariates];
    //This will hold the inverse of the covariatesXtranspose matrix
    float* cXtInverse[numCovariates];

    //Instead of vectors across space (i.e. at time zero we have a value for each point in the 
    // brain), we need vectors across time (i.e. for a specific location in the brain, 
    // we want a vector of values for that space at each time shot.
    // We will have a vector for each brain location (so x*y*z vectors).
    // Each of these vectors wil have "timeCount" entries, one for each time shot.
    float* brainLocations[x*y*z];

    // malloc the data structures
    for (i = 0; i < timeCount; i++)
    {
	dirtyBrainFiles[i] = (float*)malloc(sizeof(float) * (x*y*z));
	cleanBrainFiles[i] = (float*)malloc(sizeof(float) * (x*y*z));
	covariates[i] = (float*)malloc(sizeof(float) * numCovariates);
    }
    for (i = 0; i < (x*y*z); i++)
    	brainLocations[i] = (float*)malloc(sizeof(float) * timeCount);
    for (i = 0; i < numCovariates; i++)
    {
	covariatesTranspose[i]  = (float*)malloc(sizeof(float) * timeCount);
	covariatesXtranspose[i] = (float*)malloc(sizeof(float) * numCovariates);
	cXtInverse[i]           = (float*)malloc(sizeof(float) * numCovariates);
    }
    
    //fill the data sturctures with dummy data
    // brainLocations will be filled from the dirtyBrainFiles data        // 1 for print, 0 for no print
    fillDataStructures(dirtyBrainFiles, cleanBrainFiles, covariates, brainLocations, 1); 
    
    

    /***************meat of the problem***********************************************/

    //create the transpose matrix
    createTranspose(covariates, covariatesTranspose, 1); // 1 for print
    //create the ( numCovariates x numCovariates )matrix X^T * X  (covariatesTranspose * covariates)
    multiplyCovariatesAndTranspose(covariatesTranspose, covariates, covariatesXtranspose, 1);
    // create the inverse of the covariatesXtranspose matrix
    createInverse(covariatesXtranspose, cXtInverse, 1);


    /*******************end meat of the problem*********************************/

    // clean up
    for (i = 0; i < timeCount; i++)
    {
	free(dirtyBrainFiles[i]);
	free(cleanBrainFiles[i]);
	free(covariates[i]);
    }//end for
    for (i = 0; i < x*y*z; i++)
	free (brainLocations[i]);
    for (i = 0; i < numCovariates; i++)
    {
	free (covariatesTranspose[i]);
	free (covariatesXtranspose[i]);
	free (cXtInverse[i]);
    }
	    
	   

}// end main

/****************************************************************************
 *
 *             Helper methods
 *
 ***************************************************************************/
/***
 *  This method fills the data strutures with systematic floats for 
 *  testing purposes.  
 *
 */
void fillDataStructures(float** dirtyBrainFiles, float** cleanBrainFiles, float** covariates, float** brainLocations, int print)
{
    int i, j;

    // fill the dirtyBrainFiles, cleanBrainFiles
    if (print)
	printf("PRINTING DIRTY BRAIN FILES\n");
    for (i = 0; i < timeCount; i++)
    {
	if (print)
	    printf("     printing file %d:  ", i);
	for (j = 0; j < x*y*z; j++)
	{
	    cleanBrainFiles[i][j] = 0.0; //initialize
	    dirtyBrainFiles[i][j] = (i*100) + j + 1;
	    if (print)
		printf("%3.2f, ", dirtyBrainFiles[i][j]);
	}//end for j
	if (print)
	    printf("\n");
    }//end for i

    if (print)
	printf("\n\nPRINTING COVARIATES\n");
    //fill covariates
    for (i = 0; i < timeCount; i++)
    {
	if (print)
	    printf("    printing file %d:   ", i);  
	for (j = 0; j < numCovariates; j++)
        {
	    covariates[i][j] = (i*10) + j + 1;
	    if (print)
		printf("%3.2f, ", covariates[i][j]);	    
	}//end for j
	if (print)
	    printf("\n");
    }//end for i



    //fill the brainLocations vectors
    if (print)
	printf("\n\nPRINTING BRAINLOCATIONS\n");
    for (i = 0; i< (x*y*z); i++)
    {
	if (print)
	    printf("     printing file %d:   ", i);
	for (j = 0; j < timeCount; j++)
	{
	    brainLocations[i][j] = dirtyBrainFiles[j][i];
	    if (print)
		printf("%3.2f, ", brainLocations[i][j]);
	}//end 
	if (print)
	    printf("\n");	    
    }//end for
    if (print)
	printf("\n\n");

}//end fillDataStructures

/***
 *  This method takes in a matrix and loads a transpose matrix from
 *  that data. 
 */
void createTranspose(float** covariates, float** covariatesTranspose, int print)
{
    // recall that covariates has  'timeCount' rows, and 'numCovariates' columns
    // the transpose will have 'numCovariates' rows and 'timeCount' columns
    int i, j;

    if (print)
	printf("PRINTING COVARIATES TRANSPOSE\n");
    for (i = 0; i < numCovariates; i++)
    {
	if (print)
	    printf("     printing row %d:    ", i);
	for (j = 0; j < timeCount; j++)
	{
	    covariatesTranspose[i][j] = covariates[j][i];
	    if (print)
		printf("%3.2f, ", covariatesTranspose[i][j]);
	}//end for j
	if(print)
	    printf("\n");
    }//end for i
    if(print)
	printf("\n\n");

}//end createTranspose

/***
 *  This method creates the numCovariates x numCovarites matrix
 *  by multiplying covariatesTranspose x covariates
 */
void multiplyCovariatesAndTranspose(float** covariatesTranspose, float** covariates, float** covariatesXtranspose, int print)
{
    //dense matrix multiply
    // replace with cannon's algorithm???
    int i, j, k;

    if (print)
	printf("PRINTING COVARIATES X TRANSPOSE\n");
    for (i = 0; i < numCovariates; i++)
    {
	if (print)
	    printf("     printing row %d:    ", i);
	for (j = 0; j < numCovariates; j++)
	{
	    covariatesXtranspose[i][j] = 0.0; //initialize
	    for (k = 0; k < timeCount; k++)
	    {
		covariatesXtranspose[i][j] += covariatesTranspose[i][k] * covariates[k][j];
	    }//end for k
	    if (print)
		printf("%5.2f, ", covariatesXtranspose[i][j]);
	}//end for j
	if (print)
	    printf("\n");
    }//end for i
    if (print)
	printf("\n\n");

    

}//end multiplyCovariatesAndTranspose

/**
 * Create the inverse of the covaratiates * transpose matrix
 */
void createInverse(float** covariatesXtranspose, float** cXtInverse, int dim, int print)
{
    int i, j, k;
    float augmentedMatrix[dim][2*dim]; //stores matrix augmented with identity matrix
    

    //first, augment with an identity matrix of similar dimensions
    for (i = 0; i < dim; i++)
    {
	for (j = dim, j < 2*dim; j++)
	{
	    if (i
	}
    }//end for i

    //use gauss-jordan elimination

}//end createInverse


