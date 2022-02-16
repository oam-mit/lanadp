#include <stdio.h>

__global__ void addKernel(int* c, const int* a, const int* b, int size) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size) {
        c[i] = a[i] + b[i];
    }
}

// Helper function for using CUDA to add vectors in parallel.
void addWithCuda(int* c, const int* a, const int* b, int size) {
    int* dev_a = NULL;
    int* dev_b = NULL;
    int* dev_c = NULL;

    // Allocate GPU buffers for three vectors (two input, one output)
    cudaMalloc((void**)&dev_c, size * sizeof(int));
    cudaMalloc((void**)&dev_a, size * sizeof(int));
    cudaMalloc((void**)&dev_b, size * sizeof(int));

    // Copy input vectors from host memory to GPU buffers.
    cudaMemcpy(dev_a, a, size * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_b, b, size * sizeof(int), cudaMemcpyHostToDevice);

    // Launch a kernel on the GPU with one thread for each element.
    // 2 is number of computational blocks and (size + 1) / 2 is a number of threads in a block


    addKernel<<<2, 64>>>(dev_c, dev_a, dev_b, size);
    
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaDeviceSynchronize();

    // Copy output vector from GPU buffer to host memory.
    cudaMemcpy(c, dev_c, size * sizeof(int), cudaMemcpyDeviceToHost);

    cudaFree(dev_c);
    cudaFree(dev_a);
    cudaFree(dev_b);
}

int main(int argc, char** argv) {
    int a[20];
    int b[20];
    int n;
    printf("Enter the size of array\n");
    scanf("%d",&n);
    printf("Enter the Array A\n");
    for(int i=0;i<n;i++){   
    scanf("%d",&a[i]);
}
printf("Enter the Array B\n");
    for(int i=0;i<n;i++){   
    scanf("%d",&b[i]);
}
    int c[20] = { 0 };

    addWithCuda(c, a, b, n);
    for(int i=0;i<n;i++){
    printf("A[%d]= %d + B[%d]=%d = %d \n", i,a[i],i,b[i],c[i]);
}
    cudaDeviceReset();

    return 0;
}
