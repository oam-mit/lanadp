#include <stdio.h>
#include <cuda.h>
#include "dcmtk/dcmimgle/dcmimage.h"
#include <fstream>

int powers[3][3]={
    {8,4,2},
    {16,0,1},
    {32,64,128}
};

//dichom, dcmtk
void TLBAP (int *arr, unsigned char *out,unsigned int rows,unsigned int cols,float threshold_1) {
    
    for(int global_row=0;global_row<rows;global_row++){
        
        for(int global_col=0;global_col<cols;global_col++){

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

    }
   



   
  
  
}


int main(int argc, char *argv[]) {

    DicomImage *img=new DicomImage(argv[1]);
    int row,column;
    int *arr;

    if(img != NULL && img->getStatus()==EIS_Normal) {

        if(img->isMonochrome()) {

            img->setMinMaxWindow();
			
			row=img->getHeight();
			column=img->getWidth();

            arr = (int *)new int[row*column];

            printf("Rows: %d\tColumns:%d\n",row,column);

            img->getOutputData(arr,row*column*sizeof(int));


            int output_rows=row-2,output_colums=column-2;
            unsigned char output[output_rows][output_colums];
            printf("Output Rows: %d\tOutput Columns:%d\n",output_rows,output_colums);

            int size_arr = row*column*sizeof(short int);
            int size_out = output_rows*output_colums*sizeof(unsigned char);



        
            float elapsed=0;
            cudaEvent_t start, stop;

            cudaEventCreate(&start);
            cudaEventCreate(&stop);

            cudaEventRecord(start, 0);


            TLBAP(arr,(unsigned char *)output,row,column,0.9);
        
            cudaEventRecord(stop, 0);
            cudaEventSynchronize (stop);

            cudaEventElapsedTime(&elapsed, start, stop);
        
            printf("Time Elapsed: %f\n",elapsed);

            // cudaMemcpy(output,d_output,size_out,cudaMemcpyDeviceToHost);

            // for(int i=0;i<row-2;i++) {
                
            //     for(int i1=0;i1<column-2;i1++) {
                    
            //         printf("%d  ",output[i][i1]);
            //     }
            //     printf("\n");
            // }

            std::ofstream outfile;

            outfile.open("readings/tlbap/sequential_tlbap.csv", std::ios_base::app); // append instead of overwrite
            outfile <<argv[1]<<","<< elapsed<<"\n"; 



        }
    }

   
   //short int arr[9]= {3,5,4,10,5,3,12,2,6};

   
 



}