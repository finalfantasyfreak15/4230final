/***
*  Keri Anderson
*  CS4230  Assignment 6
*  Due Friday, Nov 22, 2013
*
*
*  Parallelize Sparse Matrix Vector Mutliplication
*  on a GPU using CUDA
*
*   - Use lab1-x.eng.utah.edu, where x can be 1, 2, 3, 4 ...
*   - The makefile provided should be used to compile 
*     the program.  Use "make" command
*   - To run:  ./spmv appu.mtx
*
*   You might have to set the LD_LIBRARY_PATH
*
*    -export
*    LD_LIBRARY_PATH = /usr/local/apps/cuda/3.2/cuda/lib64   (if using bash)
*
*    or 
*
*    setenv  LD_LIBRARY_PATH /usr/local/apps/cuda/3.2/cuda/lib64
*
*    Expect speedup to be around 2
*
*
*/

#include <stdio.h>
#include <cutil.h>
extern int cudaMalloc();


/***********************************************************
 *
 *         GLOBAL FUNCTIONS
 *
 **********************************************************/
__global__ void sparse_GPU(float *data, float *t, float *b, int *ptr, 
                           int *indices) {

// ToDo: implement the cuda SpMV code here

//  Individual threads run this code.  
//  We need to figure out where to tell
//  *this* thread to start processing. 
//
//  Suppose A was an array of 16 rows, and we 
//  have 4 blocks with 4 threads each.  
//  Then if *this* thread is in block 2 and 
//  thread number 3 within that block, 
//  Then the row we want to process is 
//  (2*4) + 3, or row 11.  
//  Block 0:                           Block 1:
//        thread 0 processes row 0;          thread 0 processes row 4; (1*4 + 0 = 4)
//        thread 1 processes row 1;          thread 1 processes row 5;
//        thread 2 processes row 2;          thread 2 processes row 6;
//        thread 3 processes row 3;          thread 3 processes row 7;
//
//  Block 2:                           Block 3:
//        thread 0 processes row 8;          thread 0 processes row 12; (3*4 + 0 = 12)
//        thread 1 processes row 9;          thread 1 processes row 13;
//        thread 2 processes row 10;         thread 2 processes row 14;
//        thread 3 processes row 11;         thread 3 processes row 15;

 
    // i tells the program to access a specific row, 
    // j advances the column

  //parallelize the outer loop
  int i = blockIdx.x *blockDim.x + threadIdx.x;
  int j;

  for (j = ptr[i]; j<ptr[i+1]; j++) {
      t[i] = t[i] + data[j] * b[indices[j]];

  }//end for j


}//end global sparse_GPU


/***********************************************************
 *
 *         MAIN
 *
 **********************************************************/
main (int argc, char **argv) {
  FILE *fp;
  char line[1024]; 
  int *ptr, *indices;
  float *data, *b, *t;
  float *d_a, *d_b, *d_c, *h_a;
  int *d_ptr, *d_indices; 
  int i,j;
  int n; // number of nonzero elements in data
  int nr; // number of rows in matrix
  int nc; // number of columns in matrix

  // Open input file and read to end of comments
  if (argc !=2) abort(); 

  if ((fp = fopen(argv[1], "r")) == NULL) {
    abort();
  }

  fgets(line, 128, fp);
  while (line[0] == '%') {
    fgets(line, 128, fp); 
  }

  // Read number of rows (nr), number of columns (nc) and
  // number of elements and allocate memory for ptr, indices, data, b and t.
  sscanf(line,"%d %d %d\n", &nr, &nc, &n);
  ptr = (int *) malloc ((nr+1)*sizeof(int));  // nr+1 because the last entry holds
  indices = (int *) malloc(n*sizeof(int));    //    info about where last row ends
  data = (float *) malloc(n*sizeof(float));
  b = (float *) malloc(nc*sizeof(float));
  t = (float *) malloc(nr*sizeof(float));

  //  The matrix will be loaded in CSR Format
  //  ptr array:  holds the 'row' information
  //  indices:    holds the 'col' information
  //  data:       holds the actual values.

  //   example:  suppose we had sparse matrix 

  //          2 0 0 5
  //          0 3 0 0
  //          1 0 1 0
  //          0 0 0 7

  //  Then our CSR representation would be

  //      ptr      indices     data       note:  for the ptr array, the first entry
  //      [0]        [0]        [2]              represents row '0' and where in the 
  //      [2]        [3]        [5]              indices array the row 0 elements start.
  //      [3]        [1]        [3]              The sencond entry represents row 1
  //      [5]        [0]        [1]              and the value there is where the elements
  //      [6]        [2]        [1]              of row 1 start:  in this case, 2, since
  //                 [3]        [7]              we had 2 elements in row 0.

  // Read data in coordinate format and initialize sparse matrix
  int lastr=0;
  for (i=0; i<n; i++) {
    int r;
    fscanf(fp,"%d %d %f\n", &r, &(indices[i]), &(data[i]));  
    indices[i]--;  // start numbering at 0
    if (r!=lastr) { 
      ptr[r-1] = i; 
      lastr = r; 
    }
  }
  ptr[nr] = n;

  // initialize t to 0 and b with random data  
  for (i=0; i<nr; i++) {
    t[i] = 0.0;
  }

  for (i=0; i<nc; i++) {
    b[i] = (float) rand()/1111111111;
  }

  // create CUDA event handles for timing purposes
  cudaEvent_t start_event, stop_event;
  float elapsed_time_seq, elapsed_time_gpu;

  CUDA_SAFE_CALL( cudaEventCreate(&start_event) );
  CUDA_SAFE_CALL( cudaEventCreate(&stop_event) );
  cudaEventRecord(start_event, 0);   

  // MAIN COMPUTATION, SEQUENTIAL VERSION
  for (i=0; i<nr; i++) {                                                      
    for (j = ptr[i]; j<ptr[i+1]; j++) {
      t[i] = t[i] + data[j] * b[indices[j]];
    }
  }
  cudaEventRecord(stop_event, 0);
  cudaEventSynchronize(stop_event);
  CUDA_SAFE_CALL( cudaEventElapsedTime(&elapsed_time_seq,start_event, stop_event) )

  // The result from the GPU should be copied to h_a
  // By convention, 'h_' is used to indiate memory
  // allocated on the Host CPU
  h_a = (float *) malloc(nr*sizeof(float));

  CUDA_SAFE_CALL( cudaEventCreate(&start_event) );
  CUDA_SAFE_CALL( cudaEventCreate(&stop_event) );
  cudaEventRecord(start_event, 0);   

  // ToDo: cuda memory allocation and copy
  //  cudaMalloc creates space on the GPU, with
  //  a memory-mapped pointer on the CPU host
  cudaMalloc((void**) &d_a, n*sizeof(float));
  cudaMalloc((void**) &d_b, nc*sizeof(float)); //number of cols
  cudaMalloc((void**) &d_c, nr*sizeof(float)); //number of rows
  cudaMalloc((void**) &d_ptr, (nr+1)*sizeof(float));
  cudaMalloc((void**) &d_indices, n*sizeof(float));

  //now copy
  cudaMemcpy(d_a, data, n*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_b, b, nc*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_c, h_a, nr*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_ptr, ptr, (nr+1)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_indices, indices, n*sizeof(float), cudaMemcpyHostToDevice);


  // ToDo: declare grid and block dimensions
  // We want each thread to perform a row's calculation.
  // Essentially, each thread will be performing a 
  // vector x vector computation (or dot product).
  // Each block has a certain amount of threads.

  // We want:  number of blocks * number of threads (in each block) = total number of rows
  // Example:  suppose we had 16 rows.  If we create 4 blocks with 4 threads in 
  //           each block, then 4*4 = 16 and we will have enough threads.
  //           Each block would be assigned 4 rows of data to compute. 
  
  // dim3 is a CUDA type
  // we are parallelizing 1 loop, so we only need one dimension
  // in our case nr/100 will divide nicely
  // grids are usually 2 dimensional
  dim3 gridDim((nr/100), 1);  //pass in 1, 1, for the 1 and z dimension so they are not being used
  dim3 blockDim(100, 1, 1);   //dim3 gridDim.x = 100;  another way to express the assignment

  // ToDo: call the sparse_GPU cuda kernel  (this is the thread code)
  sparse_GPU<<<gridDim, blockDim>>>(d_a, d_c, d_b, d_ptr, d_indices);

  // ToDo: copy back the result
  // We only need to copy back the results
  cudaMemcpy(h_a, d_c, nr*sizeof(float), cudaMemcpyDeviceToHost);  //these are the results
 

  cudaThreadSynchronize();
  cudaEventRecord(stop_event, 0);
  cudaEventSynchronize(stop_event);
  CUDA_SAFE_CALL( cudaEventElapsedTime(&elapsed_time_gpu,start_event, stop_event) )

  CUTBoolean res = cutComparefe( h_a, t, nr, 1.0);
  if (res == 1) {
    printf("VALID!\n  Sequential Time: %.2f msec\n  Parallel Time: %.2f msec\n Speedup = %.2f\n", elapsed_time_seq, elapsed_time_gpu, elapsed_time_seq/elapsed_time_gpu);
  }// end if res == 1
  else {
   printf("INVALID...\n");
   for (i=0; i<nr; i++) {
     if (abs(h_a[i]-t[i]) > 1.0) {
       printf("i=%d, h_a[i]=%f, t[i]=%f\n", i, h_a[i], t[i]);  
       break;
     }//end if
   }//end for i 
  }//end else

  //CPU free
  free(h_a);
  free(ptr);
  free(indices);
  free(data);
  free(b);
  free(t);


  //GPU free
  CUDA_SAFE_CALL(cudaFree(d_a));
  CUDA_SAFE_CALL(cudaFree(d_b));
  CUDA_SAFE_CALL(cudaFree(d_c));
  CUDA_SAFE_CALL(cudaFree(d_ptr));
  CUDA_SAFE_CALL(cudaFree(d_indices));



}//end main


