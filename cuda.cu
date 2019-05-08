// Include required header files
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#include <cuda_profiler_api.h>

// Constants used for determining what print statements are executed.
#define PRINT_ERR 0
#define PRINT_TIME 0
#define PRINT_RESULT 1

// Function prototypes
char *create_input(int);
char *sobel(char *);

__global__ void VectorComputeSobel(char *, char *);

// Code for calculating elapsed time between two points in the code adpated from PCA
//    (EEL 6763) Lecture - Spring 2019 semester.
typedef double ttype;

ttype tdiff(struct timespec a, struct timespec b)
{
  ttype dt = (( b.tv_sec - a.tv_sec) + ( b.tv_nsec - a.tv_nsec ) / 1E9);
  return dt;
}

struct timespec now()
{
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  return t;
}

// Create a vector of size image_size filled wit random integers between 0 and 255, inclusive.
char *create_input(int image_size) {
  
  int i;

  char *image;

  image = (char *) malloc(image_size*sizeof(char));

  for (i = 0; i < image_size; ++i) {
    image[i] = (int) (rand() % 256);
  }

  return image;

}

// Create a vector of size image_size filled with zeros.
char *initialize_output(int image_size) {

  int i;

  char *image;

  image = (char *) malloc(image_size*sizeof(char));

  for (i = 0; i < image_size; ++i) {
    image[i] = 0;
  }

  return image;

}

// Serial implementation of the 3x3 Sobel filter.
char *sobel_serial(int width, char *image) {

  int i;
  int image_size = width*width;

  char *sobel_image;

  int Gx[9] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
  int Gy[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};

  int S_1, S_2;
  double G_b;

  sobel_image = (char *) malloc(image_size*sizeof(char));

  for (i = 0; i < image_size; ++i) {

    if (i/width == 0 || i/width == width-1 || i%width == 0 || i%width == width-1) {
      sobel_image[(i/width)*width + (i%width)] = 0;
      continue;
    }

    S_1 = (image[((i/width)-1)*width + ((i%width)-1)] & 0x000000FF)*Gx[0] +
          //(image[i-1][j  ] & 0x000000FF)*Gx[1] +
          (image[((i/width)-1)*width + ((i%width)+1)] & 0x000000FF)*Gx[2] +
          (image[((i/width)  )*width + ((i%width)-1)] & 0x000000FF)*Gx[3] +
          //(image[i  ][j  ] & 0x000000FF)*Gx[4] +
          (image[((i/width)  )*width + ((i%width)+1)] & 0x000000FF)*Gx[5] +
          (image[((i/width)+1)*width + ((i%width)-1)] & 0x000000FF)*Gx[6] +
          //(image[i+1][j  ] & 0x000000FF)*Gx[7] +
          (image[((i/width)+1)*width + ((i%width)+1)] & 0x000000FF)*Gx[8];

    S_2 = (image[((i/width)-1)*width + ((i%width)-1)] & 0x000000FF)*Gy[0] +
          (image[((i/width)-1)*width + ((i%width)  )] & 0x000000FF)*Gy[1] +
          (image[((i/width)-1)*width + ((i%width)+1)] & 0x000000FF)*Gy[2] +
          //(image[i  ][j-1] & 0x000000FF)*Gy[3] +
          //(image[i  ][j  ] & 0x000000FF)*Gy[4] +
          //(image[i  ][j+1] & 0x000000FF)*Gy[5] +
          (image[((i/width)+1)*width + ((i%width)-1)] & 0x000000FF)*Gy[6] +
          (image[((i/width)+1)*width + ((i%width)  )] & 0x000000FF)*Gy[7] +
          (image[((i/width)+1)*width + ((i%width)+1)] & 0x000000FF)*Gy[8];
      
    G_b = sqrt((double) (S_1*S_1 + S_2*S_2));
      
    if (G_b < 0.0)
      G_b = 0.0;

    sobel_image[(i/width)*width + (i%width)] = (int) round(G_b);
  }

  return sobel_image;

}

// GPU implementation of the 3x3 Sobel filter.
__global__ void VectorComputeSobel(char *image, char *sobel_image, int width) {

  int Gx[9] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
  int Gy[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};

  int i = threadIdx.x + (blockIdx.x * blockDim.x);

  int S_1, S_2;
  double G_b;

  if (i/width == 0 || i/width == width-1 || i%width == 0 || i%width == width-1) {
      sobel_image[(i/width)*width + (i%width)] = 0;
      return;
  }

  S_1 = (image[((i/width)-1)*width + ((i%width)-1)] & 0x000000FF)*Gx[0] +
        //(image[i-1][j  ] & 0x000000FF)*Gx[1] +
        (image[((i/width)-1)*width + ((i%width)+1)] & 0x000000FF)*Gx[2] +
        (image[((i/width)  )*width + ((i%width)-1)] & 0x000000FF)*Gx[3] +
        //(image[i  ][j  ] & 0x000000FF)*Gx[4] +
        (image[((i/width)  )*width + ((i%width)+1)] & 0x000000FF)*Gx[5] +
        (image[((i/width)+1)*width + ((i%width)-1)] & 0x000000FF)*Gx[6] +
        //(image[i+1][j  ] & 0x000000FF)*Gx[7] +
        (image[((i/width)+1)*width + ((i%width)+1)] & 0x000000FF)*Gx[8];

  S_2 = (image[((i/width)-1)*width + ((i%width)-1)] & 0x000000FF)*Gy[0] +
        (image[((i/width)-1)*width + ((i%width)  )] & 0x000000FF)*Gy[1] +
        (image[((i/width)-1)*width + ((i%width)+1)] & 0x000000FF)*Gy[2] +
        //(image[i  ][j-1] & 0x000000FF)*Gy[3] +
        //(image[i  ][j  ] & 0x000000FF)*Gy[4] +
        //(image[i  ][j+1] & 0x000000FF)*Gy[5] +
        (image[((i/width)+1)*width + ((i%width)-1)] & 0x000000FF)*Gy[6] +
        (image[((i/width)+1)*width + ((i%width)  )] & 0x000000FF)*Gy[7] +
        (image[((i/width)+1)*width + ((i%width)+1)] & 0x000000FF)*Gy[8];

  G_b = sqrt((double) (S_1*S_1 + S_2*S_2));

  if (G_b < 0.0)
    G_b = 0.0;

  sobel_image[(i/width)*width + (i%width)] = (int) round(G_b);

  return;

}

int main(int argc, char *argv[])
{
  // Receive command line arguments.
  if (argc != 2) {
    printf("usage: ./hw4_a1 <image size>\n\n");
    return 1;
  }

  // Calculate image size.
  int width = atoi(argv[1]);
  int image_size = width*width;

  if (PRINT_TIME)
    printf("image size = %d\n", image_size);

  // Declare useful variables.
  struct timespec b_init, e_init, b_cpyHD, e_cpyHD, b_kernel, e_kernel, b_cpyDH, e_cpyDH, b_serial, e_serial;
  char *image, *sobel_image, *image_serial, *sobel_image_serial;
  FILE *output_fp, *output_serial_fp;

  // Start CUDA profiler.
  cudaProfilerStart();

  // Initialize input data.
  b_init = now();
  srand(time(NULL));

  image = create_input(image_size);
  sobel_image = initialize_output(image_size);

  //image_serial = create_input(image_size);
  sobel_image_serial = initialize_output(image_size);

  int byte_size = image_size*sizeof(char);
  e_init = now();

  // Allocate device memory.
  char *d_image, *d_sobel_image;

  cudaError_t err_image = cudaMalloc((void**) &d_image, byte_size);
  if (PRINT_ERR)
    printf("CUDA malloc d_image: %s\n",cudaGetErrorString(err_image));

  cudaError_t err_sobel_image = cudaMalloc((void**) &d_sobel_image, byte_size);
  if (PRINT_ERR)
    printf("CUDA malloc d_sobel_image: %s\n",cudaGetErrorString(err_sobel_image));  

  // Copy data from host memory to device memory.
  b_cpyHD = now();
  err_image = cudaMemcpy(d_image, image, byte_size, cudaMemcpyHostToDevice);
  if (PRINT_ERR)
    printf("CUDA Memcpy image->d_image: %s\n",cudaGetErrorString(err_image));

  err_sobel_image = cudaMemcpy(d_sobel_image, sobel_image, byte_size, cudaMemcpyHostToDevice);
  if (PRINT_ERR)
    printf("CUDA Memcpy sobel_image->d_sobel_image: %s\n",cudaGetErrorString(err_sobel_image));
  e_cpyHD = now();

  // Launch the device kernel.
  b_kernel = now();
  VectorComputeSobel<<<width,width>>>(d_image, d_sobel_image, width);
  cudaError_t err_VCS = cudaDeviceSynchronize();
  e_kernel = now();

  // Copy data from device memory to host memory.
  b_cpyDH = now();
  err_sobel_image = cudaMemcpy(sobel_image, d_sobel_image, byte_size, cudaMemcpyDeviceToHost);
  if (PRINT_ERR)
    printf("CUDA Memcpy d_sobel_image->sobel_image: %s\n",cudaGetErrorString(err_sobel_image));
  e_cpyDH = now();

  // Free the device memory.
  cudaFree(d_image);
  cudaFree(d_sobel_image);

  // Stop the CUDA profiler.
  cudaProfilerStop();

  // Run the serial implementation of the 3x3 Sobel filter.
  b_serial = now();
  sobel_image_serial = sobel_serial(width, image);
  e_serial = now();
  
  // Store the output from the serial and GPU implementations of the 3x3 Sobel filter.
  //    Used for comparison to check for correctness of the code using the linux diff command.
  output_fp = fopen("output.txt", "w");
  output_serial_fp = fopen("output_serial.txt", "w");

  int i;
  for (i = 0; i < image_size; ++i) {
    fprintf(output_fp, "%d\n", (sobel_image[i] & 0x000000FF));
    fprintf(output_serial_fp, "%d\n", (sobel_image_serial[i] & 0x000000FF));
  }

  // Close the files and free the pointers.
  fclose(output_fp);
  fclose(output_serial_fp);

  free(image);
  free(sobel_image);

  //free(image_serial);
  free(sobel_image_serial);
  
  // Print elapsed times for the different portions of the code.
  if (PRINT_ERR) {
    printf("Time elapsed for initialization: %.8f sec\n", tdiff(b_init, e_init));
    printf("Time elapsed for memory copy from host to device: %.8f sec\n", tdiff(b_cpyHD, e_cpyHD));
  }
  if (PRINT_RESULT) {
    printf("Time elapsed for kernel: %.8f sec\n", tdiff(b_kernel, e_kernel));
//    printf("%.8f ", tdiff(b_kernel, e_kernel));
  }
  if (PRINT_ERR)
    printf("Time elapsed for memory copy from device to host: %.8f sec\n", tdiff(b_cpyDH, e_cpyDH));
  if (PRINT_RESULT) {
    printf("Time elapsed for serial execution: %.8f sec\n", tdiff(b_serial, e_serial));
//    printf("%.8f\n", tdiff(b_serial, e_serial));
  }

  return 0;

}
