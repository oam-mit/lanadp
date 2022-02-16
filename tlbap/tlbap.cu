#include <stdio.h>
#include <cuda.h>
#include "dcmtk/dcmimgle/dcmimage.h"
#include <fstream>

__constant__ int powers[3][3]={
    {8,4,2},
    {16,0,1},
    {32,64,128}
};



//dichom, dcmtk
//arr to unsigned short int
__global__ void TLBAP (unsigned int *arr, unsigned char *out,unsigned int rows,unsigned int cols,float threshold_1) {
    
    
    int global_row = threadIdx.y + blockDim.y * blockIdx.y;
    int global_col = threadIdx.x + blockDim.x * blockIdx.x;

    if(global_row <rows && global_col<cols){
         if(global_row == 0 ||global_col ==0 || global_col==cols-1 || global_row==rows-1 ){
        //do nothing
        }

          else {
              
              short int max = arr[(global_row-1)*cols+global_col];
              
              for(int i=global_row-1;i<=global_row+1;i++) {
                  for(int i1=global_col-1;i1<=global_col+1;i1++){
                      if(i!=global_row || i1!=global_col){
                          if(arr[i*cols+i1]>max){
                              max=arr[i*cols+i1];
                          }
                      }
                  }
              }
            
              float threshold_2 = threshold_1 * (float)max;

              short int ele = arr[global_row*cols + global_col];
              //printf("row:%d\tcol:%d\tThres:%f\nMax:%d\n",global_row,global_col,threshold_2,max);

              int ans=0;
              int power_i=0;

              for(int i=global_row-1;i<=global_row+1;i++){
                  int power_i1=0;
                  for(int i1=global_col-1; i1<=global_col+1 ;i1++) {
                      if(i!=global_row || i1!=global_col) {
                          if(arr[i*cols+i1]>ele && (float)arr[i*cols+i1]>threshold_2) {
                              ans+=powers[power_i][power_i1];
                          }
                      }
                      power_i1++;
                  }
                  power_i++;
              }

                out[(global_row-1)*(cols-2)+(global_col-1)] = ans;
                
        }
    }

   
  
  
}


int main(int argc, char *argv[]) {

    int row,column;
    
    DicomImage *img=new DicomImage(argv[1]);
    

    if(img != NULL && img->getStatus()==EIS_Normal) {

        if(img->isMonochrome()) {

            img->setMinMaxWindow();
			
			row=img->getHeight();
			column=img->getWidth();

           

            unsigned int *arr = new unsigned int [row*column];

            printf("Rows: %d\tColumns:%d\n",row,column);

            img->getOutputData(arr,row*column*sizeof(unsigned int));

            unsigned char output[row-2][column-2];

            int output_rows=row-2,output_colums=column-2;

            int size_arr = row*column*sizeof(unsigned int);
            int size_out = output_rows*output_colums*sizeof(unsigned char);

            //memory variables
            unsigned int *d_arr; 
            unsigned char *d_output;

            //allocating memory
            cudaMalloc((void **)&d_arr,size_arr);
            cudaMalloc((void **)&d_output,size_out);
        
            //copying data
            cudaMemcpy(d_arr,arr,size_arr,cudaMemcpyHostToDevice); 

	//     float no_threads;
	// 	printf("Enter threads\n");
	// scanf("%f",&no_threads);


            dim3 numberOfBlocks(ceil((float)row/atof(argv[2])),ceil((float)column/atof(argv[2])),1);
            dim3 numberOfThreads(atoi(argv[2]),atoi(argv[2]),1);

        
            float elapsed=0;
            cudaEvent_t start, stop;

            cudaEventCreate(&start);
            cudaEventCreate(&stop);

            cudaEventRecord(start, 0);

            for(int i=1;i<=100;i++)
                TLBAP<<<numberOfBlocks,numberOfThreads>>>(d_arr,d_output,row,column,0.9);
        
            cudaEventRecord(stop, 0);
            cudaEventSynchronize (stop);

            cudaEventElapsedTime(&elapsed, start, stop);
        
            printf("Time Elapsed: %f\n",elapsed);

            cudaMemcpy(output,d_output,size_out,cudaMemcpyDeviceToHost);
              

            // for(int i=0;i<row-2;i++) {
                
            //     for(int i1=0;i1<column-2;i1++) {
                    
            //         printf("%d  ",output[i][i1]);
            //     }
            //     printf("\n");
            // }
            std::ofstream outfile;

            outfile.open("readings/tlbap/tlbap.csv", std::ios_base::app); // append instead of overwrite
            outfile <<argv[1]<<","<< argv[2]<<","<<elapsed<<"\n"; 


        }
    }

    
 
   

   

  



}


