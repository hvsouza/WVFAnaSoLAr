// ________________________________________ //
// Author: Henrique Souza
// Filename: denoiser.C
// Created: 2024
// ________________________________________ //
#include "WVFAnaSoLAr.h"
class DENOISE{

  public:
    template<class T>
    void TV1D_denoise(T* input, T* output, unsigned int width, const Double_t lambda) {

      unsigned int* indstart_low = (unsigned int*) malloc(sizeof *indstart_low * width);
      unsigned int* indstart_up = (unsigned int*) malloc(sizeof *indstart_up * width);
      unsigned int j_low = 0, j_up = 0, jseg = 0, indjseg = 0, i=1, indjseg2, ind;
      double output_low_first = input[0]-lambda;
      double output_low_curr = output_low_first;
      double output_up_first = input[0]+lambda;
      double output_up_curr = output_up_first;
      const double twolambda=2.0*lambda;
      if (width==1) {output[0]=input[0]; return;}
      indstart_low[0] = 0;
      indstart_up[0] = 0;
      width--;
      for (; i<width; i++) {
        if (input[i]>=output_low_curr) {
          if (input[i]<=output_up_curr) {
            output_up_curr+=(input[i]-output_up_curr)/(i-indstart_up[j_up]+1);
            output[indjseg]=output_up_first;
            while ((j_up>jseg)&&(output_up_curr<=output[ind=indstart_up[j_up-1]]))
              output_up_curr+=(output[ind]-output_up_curr)*
                ((double)(indstart_up[j_up--]-ind)/(i-ind+1));
            if (j_up==jseg) {
              while ((output_up_curr<=output_low_first)&&(jseg<j_low)) {
                indjseg2=indstart_low[++jseg];
                output_up_curr+=(output_up_curr-output_low_first)*
                  ((double)(indjseg2-indjseg)/(i-indjseg2+1));
                while (indjseg<indjseg2) output[indjseg++]=output_low_first;
                output_low_first=output[indjseg];
              }
              output_up_first=output_up_curr;
              indstart_up[j_up=jseg]=indjseg;
            } else output[indstart_up[j_up]]=output_up_curr;
          } else
            output_up_curr=output[i]=input[indstart_up[++j_up]=i];
          output_low_curr+=(input[i]-output_low_curr)/(i-indstart_low[j_low]+1);
          output[indjseg]=output_low_first;
          while ((j_low>jseg)&&(output_low_curr>=output[ind=indstart_low[j_low-1]]))
            output_low_curr+=(output[ind]-output_low_curr)*
              ((double)(indstart_low[j_low--]-ind)/(i-ind+1));
          if (j_low==jseg) {
            while ((output_low_curr>=output_up_first)&&(jseg<j_up)) {
              indjseg2=indstart_up[++jseg];
              output_low_curr+=(output_low_curr-output_up_first)*
                ((double)(indjseg2-indjseg)/(i-indjseg2+1));
              while (indjseg<indjseg2) output[indjseg++]=output_up_first;
              output_up_first=output[indjseg];
            }
            if ((indstart_low[j_low=jseg]=indjseg)==i) output_low_first=output_up_first-twolambda;
            else output_low_first=output_low_curr;
          } else output[indstart_low[j_low]]=output_low_curr;
        } else {
          output_up_curr+=((output_low_curr=output[i]=input[indstart_low[++j_low] = i])-
                           output_up_curr)/(i-indstart_up[j_up]+1);
          output[indjseg]=output_up_first;
          while ((j_up>jseg)&&(output_up_curr<=output[ind=indstart_up[j_up-1]]))
            output_up_curr+=(output[ind]-output_up_curr)*
              ((double)(indstart_up[j_up--]-ind)/(i-ind+1));
          if (j_up==jseg) {
            while ((output_up_curr<=output_low_first)&&(jseg<j_low)) {
              indjseg2=indstart_low[++jseg];
              output_up_curr+=(output_up_curr-output_low_first)*
                ((double)(indjseg2-indjseg)/(i-indjseg2+1));
              while (indjseg<indjseg2) output[indjseg++]=output_low_first;
              output_low_first=output[indjseg];
            }
            if ((indstart_up[j_up=jseg]=indjseg)==i) output_up_first=output_low_first+twolambda;
            else output_up_first=output_up_curr;
          } else output[indstart_up[j_up]]=output_up_curr;
        }
      }
      /* here i==width (with value the actual width minus one) */
      if (input[i]+lambda<=output_low_curr) {
        while (jseg<j_low) {
          indjseg2=indstart_low[++jseg];
          while (indjseg<indjseg2) output[indjseg++]=output_low_first;
          output_low_first=output[indjseg];
        }
        while (indjseg<i) output[indjseg++]=output_low_first;
        output[indjseg]=input[i]+lambda;
      } else if (input[i]-lambda>=output_up_curr) {
        while (jseg<j_up) {
          indjseg2=indstart_up[++jseg];
          while (indjseg<indjseg2) output[indjseg++]=output_up_first;
          output_up_first=output[indjseg];
        }
        while (indjseg<i) output[indjseg++]=output_up_first;
        output[indjseg]=input[i]-lambda;
      } else {
        output_low_curr+=(input[i]+lambda-output_low_curr)/(i-indstart_low[j_low]+1);
        output[indjseg]=output_low_first;
        while ((j_low>jseg)&&(output_low_curr>=output[ind=indstart_low[j_low-1]]))
          output_low_curr+=(output[ind]-output_low_curr)*
            ((double)(indstart_low[j_low--]-ind)/(i-ind+1));
        if (j_low==jseg) {
          if (output_up_first>=output_low_curr)
            while (indjseg<=i) output[indjseg++]=output_low_curr;
          else {
            output_up_curr+=(input[i]-lambda-output_up_curr)/(i-indstart_up[j_up]+1);
            output[indjseg]=output_up_first;
            while ((j_up>jseg)&&(output_up_curr<=output[ind=indstart_up[j_up-1]]))
              output_up_curr+=(output[ind]-output_up_curr)*
                ((double)(indstart_up[j_up--]-ind)/(i-ind+1));
            while (jseg<j_up) {
              indjseg2=indstart_up[++jseg];
              while (indjseg<indjseg2) output[indjseg++]=output_up_first;
              output_up_first=output[indjseg];
            }
            indjseg=indstart_up[j_up];
            while (indjseg<=i) output[indjseg++]=output_up_curr;
          }
        } else {
          while (jseg<j_low) {
            indjseg2=indstart_low[++jseg];
            while (indjseg<indjseg2) output[indjseg++]=output_low_first;
            output_low_first=output[indjseg];
          }
          indjseg=indstart_low[j_low];
          while (indjseg<=i) output[indjseg++]=output_low_curr;
        }
      }
      free(indstart_low);
      free(indstart_up);
    }


    template<class T>
    void cumulativeAverage(T* v, T* res, Int_t myinte = 10, Double_t start = 0, Double_t finish = 0, Bool_t eco = false){

      Double_t sum = 0;
      Double_t avgi = 0;
      Int_t count = 0;
      Double_t normfactor = myinte+1;
      if (eco){
        normfactor = pow(2., myinte); // each point contributes much less (gain some spaces)
      }


      for(Int_t i = 0; i < finish; i++){

        if(i<start){ // make it to start at i = 5
          res[i] = v[i];
        }
        else{
          count += 1;
          if (count > myinte){
            res[i] = avgi + (v[i] - avgi)/normfactor;
            avgi = res[i];
          }
          else if( count == myinte ){
            avgi = sum/myinte;
            res[i] = v[i];
          }
          else{
            res[i] = v[i];
            sum+=v[i];
          }
        }
      }

    }


    template<class T>
    void movingAverage(T* v, T* res, Int_t myinte = 10, Double_t start = 0, Double_t finish = 0){

      Int_t n = finish;
      Int_t midpoint = 0;
      Int_t width = 0;
      Double_t sum = 0;

      if(myinte%2==0){ // if it is even, we insert the middle point, e.g. 8 interactions takes 4 before, mid, 4 later
        myinte+=1;
      }
      midpoint = (myinte-1)/2; // e.g. 9 interactions the midpoint will be 5
      width = myinte;


      for(Int_t i = 0; i < n; i++){

        if(i<midpoint || i>(n-(midpoint+1))){ // make it to start at i = 5 and finish at i = (3000-5) = 2995
          // res[i] = v[i];
          res[i] = 0;
        }
        else if (i > start && i < finish){
          for(Int_t j = (i-midpoint); j < (i+midpoint+1); j++) { //first one: from j = (5-5); j<(5+5)
            sum = sum+v[j];

            //                 cout << sum << endl;
          }
          res[i] = (sum/width);
        }
        else{
          res[i] = 0;
        }
        sum=0;
      }

    }
};



