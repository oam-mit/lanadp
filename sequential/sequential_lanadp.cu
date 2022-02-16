#include <stdio.h>
#include <cuda.h>
#include "dcmtk/dcmimgle/dcmimage.h"

int powers[9]={0,1,2,4,8,16,32,64,128};

int getGlobalID(int row,int col,int cols){
    return row*cols+col;
}

int getCircularIndex(int i){
    
    if (i<1) {
        return 8 + i;
    }
    else if (i>8 )
      return (i%9 )+1;
    else 
      return i;
  
}

//dichom, dcmtk
void LANADP (int *arr, unsigned char *out,unsigned int rows,unsigned int cols) {
    
    for(int global_row=0;global_row<rows;global_row++){
        
        for(int global_col=0;global_col<cols;global_col++){

            if(global_row <rows && global_col<cols){
                if(global_row == 0 ||global_col ==0 || global_col==cols-1 || global_row==rows-1 ){
               //do nothing
               }
       
                 else {
                     
                     int mapping[9];
                     mapping[1]=arr[getGlobalID(global_row,global_col+1,cols)];
                     mapping[2]=arr[getGlobalID(global_row-1,global_col+1,cols)];
                     mapping[3]=arr[getGlobalID(global_row-1,global_col,cols)];
                     mapping[4]=arr[getGlobalID(global_row-1,global_col-1,cols)];
                     mapping[5]=arr[getGlobalID(global_row,global_col-1,cols)];
                     mapping[6]=arr[getGlobalID(global_row+1,global_col-1,cols)];
                     mapping[7]=arr[getGlobalID(global_row+1,global_col,cols)];
                     mapping[8]=arr[getGlobalID(global_row+1,global_col+1,cols)];
       
       
                     int ans=0;
                     for (int i=1;i<=8;i++) {
                         
                         float avg1 = ((float)(mapping[getCircularIndex(i+1)] + mapping[getCircularIndex(i+2)]))/2.0;
                         float avg2 = ((float)(mapping[getCircularIndex(i-1)] + mapping[getCircularIndex(i-2)]))/2.0;
       
                         //printf("Number: %d, Avg1: %f, Avg2: %f\n",mapping[i],avg1,avg2);
       
                         if (avg1>=(float)mapping[i] && avg2>=(float)mapping[i] || (avg1<=(float)mapping[i] && avg2<=(float)mapping[i]) ){
                             
                             //printf("Adding:%d\n",powers[i]);
                             ans+=powers[i];
                         }
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


            LANADP(arr,(unsigned char *)output,row,column);
        
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

            outfile.open("readings/sequential_lanadp.csv", std::ios_base::app); // append instead of overwrite
            outfile <<argv[1]<<","<< elapsed<<"\n"; 



        }
    }

   
   //short int arr[9]= {3,5,4,10,5,3,12,2,6};

   
 



}