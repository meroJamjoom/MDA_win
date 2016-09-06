// ==========================================================================
// $Id: PoissonSystem.C 999 2014-05-28 15:07:31Z heidrich $
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2007-, UBC
// 
// Creator: heidrich ()
// Email:   heidrich@cs.ubc.ca
// ==========================================================================

#ifndef IMAGESPACESYSTEMS_POISSONSYSTEM_C
#define IMAGESPACESYSTEMS_POISSONSYSTEM_C

#define MIN_HALVING_DIMENSION 5
#define MAX_GRID_LEVEL 8

#include "PoissonSystem.hh"
#include "MDA/GeometricTransform/Scaling.hh"
#include "MDA/Resampling/Resampler.hh"

namespace MDA {
  // the "using" statements have to be inside the MDA scope so
  // that inclusion of C files as required by gcc does not yield
  // problems with other packages!
  using namespace std;

  template<class T>
  inline
  PoissonSystem<T>*
  PoissonSystem<T>::createHalfSystem( Array<T> &a,
				      ChannelList &gradientCh, 
				      unsigned constraintCh, 
				      unsigned maskCh, 
				      bool computeDivergence, 
				      unsigned sourceCh, 
				      bool cleanMatrix)
  {  
    int i,j,k,t;
    CoordinateVector v=a.getDimension();
    int totalPts=1, mindim=v.vec[0];
    for (i=0; i<v.vec.size(); i++){     
      v.vec[i]=(v.vec[i]+1)/2; 
      totalPts*=v.vec[i];
      if (v.vec[i]<mindim) mindim=v.vec[i];
    }
   
    if (mindim<MIN_HALVING_DIMENSION || gridLevel+1>=MAX_GRID_LEVEL) 
      return NULL; 
      
    PoissonSystem<T>* halfSystem=new PoissonSystem<T>(gridLevel+1);
    unsigned dim= a.getDimension().vec.size();

    a2=Array<T>(v); 
    int uc=halfRescale(a,a2,maskCh,constraintCh,gradientCh);   
    Vector rhs=Vector(totalPts);
    halfSystem->_setup(a2,rhs, gradientCh, constraintCh, maskCh, 
		       computeDivergence, sourceCh, cleanMatrix, true);
    return halfSystem;
  }
  
  template<class T>
  inline
  void 
  PoissonSystem<T>::prolongate(const Vector &half, Vector &v)
  { 
    vector < long int > dimensions = getDimensions().vec;
    long int numDimensions = dimensions.size();
    vector < long int > halfDimensions = halfSystem->getDimensions().vec;
    if(numDimensions == 1){
      int x=dimensions[0];
      int x2=halfDimensions[0];
      int numCoordsUnitBox = 2;
      int ix[] = {0,1};
      for(int i = 0; i < x2; ++i){
        for(int k=0; k<2; ++k){
          int halfSystemIndex = i;
          if (((PoissonSystem<T>*)halfSystem)->maskPtr[halfSystemIndex] == valueMask) 
                continue; 	             
          if ((2*i+ix[k])<x){
            v[(2*i+ix[k])]=half[halfSystemIndex];
          }
        }
      }
    }else if(numDimensions == 2){ 
      int x=dimensions[1];
      int y=dimensions[0];
      int x2=halfDimensions[1];
      int y2=halfDimensions[0];
      errorCond( x2*y2== half.getSize() && x*y== v.getSize(),
	         "  incompatible matrix/vector dimensions" );
      int ix[]={0,0,1,1},iy[]={0,1,0,1};
      for (int i=0; i<x2; i++){
        for (int j=0; j<y2; j++){
	  for (int k=0; k<4; k++){
            int halfSystemIndex = i*y2+j;
	    if(((PoissonSystem<T>*)halfSystem)->maskPtr[halfSystemIndex] == 
              valueMask) continue;     
            if((2*i+ix[k])<x && (2*j+iy[k])<y) v[(2*i+ix[k])*y+(2*j+iy[k])]=
              half[halfSystemIndex];
	  }
        }
      }
    }else if(numDimensions == 3){ 
      int x=dimensions[2];
      int y=dimensions[1];
      int z=dimensions[0];
      int x2=halfDimensions[2];
      int y2=halfDimensions[1];     
      int z2=halfDimensions[0];
      int numCoordsUnitBox = 8;
      int ix[]={0,0,0,0,1,1,1,1},iy[]={0,0,1,1,0,0,1,1},iz[]={0,1,0,1,0,1,0,1};
      for(int i =0; i < x2; ++i){
        for(int j =0; j < y2; ++j){
          for(int k =0; k < z2; ++k){
            for(int l =0; l < numCoordsUnitBox; ++l){
              int halfSystemIndex = i*z2*y2+j*z2+k;
              if (((PoissonSystem<T>*)halfSystem)->maskPtr[halfSystemIndex] == valueMask) 
                continue; 	             
              if ((2*i+ix[l])<x && (2*j+iy[l])<y && 
                (2*k+iz[l])<z){
                v[(2*i+ix[l])*z*y+(2*j+iy[l])*z+
                  (2*k+iz[l])]=half[halfSystemIndex];
              }  
            }
          }
        }
      } 
    }else if(numDimensions == 4){
      int x=dimensions[3];
      int y=dimensions[2];
      int z=dimensions[1];
      int w=dimensions[0];
      int x2=halfDimensions[3];
      int y2=halfDimensions[2];     
      int z2=halfDimensions[1];
      int w2=halfDimensions[0];
      int numCoordsUnitBox = 16;
      int ix[]={0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1};
      int iy[]={0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1};
      int iz[]={0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1};
      int iw[]={0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1};
      for(int i =0; i < x2; ++i){
        for(int j =0; j < y2; ++j){
          for(int k =0; k < z2; ++k){
            for(int l =0; l < w2; ++l){
              for(int m =0; m < numCoordsUnitBox; ++m){
                int halfSystemIndex = i*w2*z2*y2 + j*w2*z2 + k*w2 + l;
                if (((PoissonSystem<T>*)halfSystem)->maskPtr[halfSystemIndex] == valueMask) 
                  continue; 	             
                if ((2*i+ix[m])<x && (2*j+iy[m])<y && 
                  (2*k+iz[m])<z && (2*l+iw[m])<w){
                  v[(2*i+ix[m])*w*z*y+
                      (2*j+iy[m])*w*z+
                        (2*k+iz[m])*w+
                          (2*l+iw[m])] = 
                             half[halfSystemIndex];
                }
              }  
            }
          }
        }
      } 
    }
    return;
  }

  template<class T>
  inline
  void 
  PoissonSystem<T>::restrict(const Vector &v, Vector &half)
  { 
    vector < long int > dimensions = getDimensions().vec;
    long int numDimensions = dimensions.size();
    vector < long int > halfDimensions = halfSystem->getDimensions().vec;
    if(numDimensions == 1){
      int x=dimensions[0];
      int x2=halfDimensions[0];
      int ix[]={0,1};
      int counter;
      int fullSystemIndex;
      int halfSystemIndex;
      int numCoordsUnitBox = 2;
      for(int i= 0; i<x2; ++i){
        halfSystemIndex=i;
	half[halfSystemIndex]=0;     
	for(int k=counter=0; k < numCoordsUnitBox; ++k){
	  fullSystemIndex=(2*i+ix[k]);
	  if((2*i+ix[k])<x && 
            ((PoissonSystem<T>*)halfSystem)->maskPtr[halfSystemIndex] ==
              maskPtr[fullSystemIndex]){
	    counter++;
	    half[halfSystemIndex]+=v[fullSystemIndex];
	  }
	}
	
	half[halfSystemIndex]/=counter;
	half[halfSystemIndex]*=numCoordsUnitBox*numDimensions; 
      } 
    }else if(numDimensions == 2){
      int ix[]={0,0,1,1},iy[]={0,1,0,1};
      int x=dimensions[1];
      int y=dimensions[0];
      int x2=halfDimensions[1];
      int y2=halfDimensions[0];
      int cnt;
      for (int i=0; i<x2; i++){
        for (int j=0; j<y2; j++){
	  int halfSystemIndex=i*y2+j;
	  half[halfSystemIndex]=0;
	  for (int k=cnt=0; k<4; k++){
	    int fullSystemIndex=(2*i+ix[k])*y+(2*j+iy[k]);
	    if ((2*i+ix[k])<x && (2*j+iy[k])<y 
	      && ((PoissonSystem<T>*)halfSystem)->maskPtr[halfSystemIndex]==
                maskPtr[fullSystemIndex]){
	      cnt++;
	      half[halfSystemIndex]+=v[fullSystemIndex];
	    }
	  }	    
	  half[halfSystemIndex]/=cnt;
	  half[halfSystemIndex]*=4*2;
        }
        vector < long int > dimensions = getDimensions().vec;
        long int numDimensions = dimensions.size();
        vector < long int > halfDimensions = halfSystem->getDimensions().vec;
      }
    }else if(numDimensions == 3){ 
      int x=dimensions[2];
      int y=dimensions[1];
      int z=dimensions[0];
      int x2=halfDimensions[2];
      int y2=halfDimensions[1];     
      int z2=halfDimensions[0];
      int ix[]={0,0,0,0,1,1,1,1},iy[]={0,0,1,1,0,0,1,1},iz[]={0,1,0,1,0,1,0,1};
      int counter;
      int fullSystemIndex;
      int halfSystemIndex;
      int numCoordsUnitBox = 8;
      for(int i= 0; i<x2; ++i){
        for(int j= 0; j<y2; ++j){
          for(int k= 0; k<z2; ++k){
            halfSystemIndex=i*z2*y2+j*z2+k;
	    half[halfSystemIndex]=0;     
	    for(int l=counter=0; l < numCoordsUnitBox; ++l){
	      fullSystemIndex=(2*i+ix[l])*z*y+(2*j+iy[l])*z+(2*k+iz[l]);
	      if((2*i+ix[l])<x && (2*j+iy[l])<y && (2*k+iz[l])<z
	        && ((PoissonSystem<T>*)halfSystem)->maskPtr[halfSystemIndex]==
                  maskPtr[fullSystemIndex]){
	        counter++;
	        half[halfSystemIndex]+=v[fullSystemIndex];
	      }
	    }
	  }

	  half[halfSystemIndex]/=counter;
	  half[halfSystemIndex]*=numCoordsUnitBox*numDimensions;
        }
      } 
    }else if(numDimensions == 4){
      int x=dimensions[3];
      int y=dimensions[2];
      int z=dimensions[1];
      int w=dimensions[0];
      int x2=halfDimensions[3];
      int y2=halfDimensions[2];     
      int z2=halfDimensions[1];
      int w2=halfDimensions[0];
      int counter;
      int fullSystemIndex;
      int halfSystemIndex;
      int numCoordsUnitBox = 16;
      int ix[]={0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1};
      int iy[]={0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1};
      int iz[]={0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1};
      int iw[]={0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1};
      for(int i= 0; i<x2; ++i){
        for(int j= 0; j<y2; ++j){
          for(int k= 0; k<z2; ++k){
            for(int l= 0; l<w2; ++l){
              halfSystemIndex=i*w2*z2*y2+j*w2*z2+k*w2+l;
	      half[halfSystemIndex]=0;     
	      for(int m=counter=0; m < numCoordsUnitBox; ++m){
	        fullSystemIndex=(2*i+ix[m])*w*z*y+(2*j+iy[m])*w*z+
                  (2*k+iz[m])*w+(2*l+iw[m]);
	        if((2*i+ix[m])<x && (2*j+iy[m])<y && (2*k+iz[m])<z && (2*l+iw[m])<w
	          && ((PoissonSystem<T>*)halfSystem)->maskPtr[halfSystemIndex]==
                    maskPtr[fullSystemIndex]){
	          counter++;
	          half[halfSystemIndex]+=v[fullSystemIndex];
                }
	      }
	    }
	  }

	  half[halfSystemIndex]/=counter;
	  half[halfSystemIndex]*=numCoordsUnitBox*numDimensions;
        }
      } 
    }

    return;
  }

  template<class T>
  inline
  int
  PoissonSystem<T>::halfRescale(Array<T> &a, 
				Array<T> &a2, 
				int maskCh,
				int constraintCh,
				ChannelList &gradientCh){
    vector < long int > dimensions = getDimensions().vec;
    long int numDimensions = dimensions.size();
    CoordinateVector v=a2.getDimension();
    int unconstrainedCount=0;      
    if(numDimensions == 1){
      int i,k,k1,t;
      int ix[]={0,1};
      int x2=a2.getDimension().vec[0];      
      int x = a.getDimension().vec[0];   
      int numCoordsUnitBox = 2;
      int li[]={maskCh,constraintCh,gradientCh.vec[0]};
      int sizeLi = 3;
      int needed=a.getNumChannels()-a2.getNumChannels(); 
      for (k=0; k<needed; ++k) 
        a2.addChannel();
      for (k1=0; k1<3; ++k1){
        k=li[k1];
        for(i=0; i<v.vec[0]; i++){                            
          a2[k][0][i]=0;      
	  if (k==maskCh){
	    a2[k][0][i]=unconstrainedMask;          
	    for (t=0; t<numCoordsUnitBox; t++){
	      T curr;
	      if (2*i+ix[t]<x)
	        curr=a[k][0][(2*i+ix[t])];
	      else curr=unconstrainedMask;
	      if (curr==valueMask || 
	        (curr==gradientMask && a2[k][0][i]!=valueMask)) 
	        a2[k][0][i]=curr;
	    }
	    if (a2[k][0][i]==unconstrainedMask) 
	      unconstrainedCount++;
	  }else if (k==constraintCh){
	    int cnt=0;
	    for (t=0; t<numCoordsUnitBox; t++){
	      if (2*i+ix[t]<x){
	        if (a[maskCh][0][(2*i+ix[t])]==
	          a2[maskCh][0][i]){
		  cnt++;
		  a2[k][0][i]+=a[k][0][(2*i+ix[t])];
		}
	      }
	    }	        	       
	    a2[k][0][i]/=cnt;
	  }else{ 
	    for (t=0; t<numCoordsUnitBox; t++){
	      T curr=0;
	      if ((2*i+ix[t]<x))
	        curr=a[k][0][(2*i+ix[t])];
	        a2[k][0][i]+=curr;
	    }
	  } 
	}
      }       
    }else if(numDimensions == 2){   
      int i,j,k,k1;
      int ix[]={0,0,1,1},iy[]={0,1,0,1};
      int y2=a2.getDimension().vec[0];
      int y = a.getDimension().vec[0];
      int x = a.getDimension().vec[1];
      int li[]={maskCh,constraintCh,gradientCh.vec[0],gradientCh.vec[1]}; 
      int needed=a.getNumChannels()-a2.getNumChannels(); 
      for (k=0; k<needed; k++) 
        a2.addChannel(); 
      for (k1=0; k1<4; k1++){
        k=li[k1];
        for (i=0; i<v.vec[1]; i++){
	  for (j=0; j<v.vec[0]; j++){
	    a2[k][0][i*y2+j]=0;
	    if (k==maskCh){
	      a2[k][0][i*y2+j]=unconstrainedMask;
	      for (int t=0; t<4; t++){
	        T curr;
	        if (2*i+ix[t]<x && 2*j+iy[t]<y)
		  curr=a[k][0][(2*i+ix[t])*y+(2*j+iy[t])];
	        else curr=unconstrainedMask;
	        if (curr==valueMask || 
		    (curr==gradientMask && a2[k][0][i*y2+j]!=valueMask)) 
	  	  a2[k][0][i*y2+j]=curr;
	      }
	      if (a2[k][0][i*y2+j]==unconstrainedMask) 
	      unconstrainedCount++;
	    }else if (k==constraintCh){
	      int cnt=0;
	      for (int t=0; t<4; t++){
	        if (2*i+ix[t]<x && 2*j+iy[t]<y){
		  if (a[maskCh][0][(2*i+ix[t])*y+(2*j+iy[t])]==
		      a2[maskCh][0][i*y2+j]){
		    cnt++;
		    a2[k][0][i*y2+j]+=a[k][0][(2*i+ix[t])*y+(2*j+iy[t])];
		  }
	        }
	      }   
	      a2[k][0][i*y2+j]/=cnt;
	    }else{
	      for (int t=0; t<4; t++){
	        T curr=0;
	        if (2*i+ix[t]<x && 2*j+iy[t]<y)
		  curr=a[k][0][(2*i+ix[t])*y+(2*j+iy[t])];
	        a2[k][0][i*y2+j]+=curr;
	      }
	    }
	  }
        }
      }    
    }else if(numDimensions == 3){
      int i,j,k,k1,t;
      int ix[]={0,0,0,0,1,1,1,1},iy[]={0,0,1,1,0,0,1,1},iz[]={0,1,0,1,0,1,0,1};
      int z2=a2.getDimension().vec[0];
      int y2=a2.getDimension().vec[1];      
      int x2=a2.getDimension().vec[2];      
      int z = a.getDimension().vec[0];
      int y = a.getDimension().vec[1];
      int x = a.getDimension().vec[2];   
      int numCoordsUnitBox = 8;
      int li[]={maskCh,constraintCh,gradientCh.vec[0],gradientCh.vec[1],gradientCh.vec[2]};
      int sizeLi = 5;
      int needed=a.getNumChannels()-a2.getNumChannels(); 
      for (k=0; k<needed; k++) 
        a2.addChannel();
      for (k1=0; k1<5; k1++){
        k=li[k1];
        for(i=0; i<v.vec[2]; i++){          
          for (j=0; j<v.vec[1]; j++){ 
            for (int l=0; l<v.vec[0]; l++){
              int halfSystemIndex = i*z2*y2+j*z2+l;                   
              a2[k][0][halfSystemIndex]=0;      
	      if (k==maskCh){
	        a2[k][0][halfSystemIndex]=unconstrainedMask;          
	        for (t=0; t<8; t++){
	          T curr;
	          if (2*i+ix[t]<x && 2*j+iy[t]<y && 2*l+iz[t]<z)
		    curr=a[k][0][(2*i+ix[t])*z*y+(2*j+iy[t])*z+(2*l+iz[t])];
	          else curr=unconstrainedMask;
	          if (curr==valueMask || 
		      (curr==gradientMask && a2[k][0][halfSystemIndex]!=valueMask)) 
	  	    a2[k][0][halfSystemIndex]=curr;
	        }
	        if (a2[k][0][halfSystemIndex]==unconstrainedMask) 
	          unconstrainedCount++;
	      }else if (k==constraintCh){
	        int cnt=0;
	        for (t=0; t<8; t++){
	          if (2*i+ix[t]<x && 2*j+iy[t]<y && 2*l+iz[t]<z){
		    if (a[maskCh][0][(2*i+ix[t])*z*y+(2*j+iy[t])*z+(2*l+iz[t])]==
		        a2[maskCh][0][halfSystemIndex]){
		      cnt++;
		      a2[k][0][halfSystemIndex]+=a[k][0][(2*i+ix[t])*z*y+
                        (2*j+iy[t])*z+(2*l+iz[t])];
		    }
	          }
	        }	        	       
	        a2[k][0][halfSystemIndex]/=cnt;
	      }else{ 
	        for (t=0; t<8; t++){
	          T curr=0;
	          if ((2*i+ix[t]<x) && (2*j+iy[t]<y) && (2*l+iz[t]<z))
		    curr=a[k][0][(2*i+ix[t])*z*y+(2*j+iy[t])*z+(2*l+iz[t])];
	          a2[k][0][halfSystemIndex]+=curr;
	        }
	      } 
	    }
          }
        }     
      }
    }else if(numDimensions == 4){
      int i,j,k,k1,t;
      int ix[]={0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1};
      int iy[]={0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1};
      int iz[]={0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1};
      int iw[]={0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1};
      int w2=a2.getDimension().vec[0];
      int z2=a2.getDimension().vec[1];
      int y2=a2.getDimension().vec[2];      
      int x2=a2.getDimension().vec[3];
      int w = a.getDimension().vec[0];
      int z = a.getDimension().vec[1];
      int y = a.getDimension().vec[2];
      int x = a.getDimension().vec[3];   
      int numCoordsUnitBox = 16;
      int li[]={maskCh,constraintCh,gradientCh.vec[0],gradientCh.vec[1],
        gradientCh.vec[2],gradientCh.vec[3]};
      int sizeLi = 6;
      int needed=a.getNumChannels()-a2.getNumChannels(); 
      for (k=0; k<needed; k++) 
        a2.addChannel();
      for (k1=0; k1<6; k1++){
        k=li[k1];
        for(i=0; i<v.vec[3]; i++){          
          for(j=0; j<v.vec[2]; j++){ 
            for(int l=0; l<v.vec[1]; l++){
              for(int m=0; m<v.vec[0]; m++){
                int halfSystemIndex = i*w2*z2*y2+j*w2*z2+l*w2+m;                  
                a2[k][0][halfSystemIndex]=0;      
	        if (k==maskCh){
	          a2[k][0][halfSystemIndex]=unconstrainedMask;          
	          for (t=0; t<numCoordsUnitBox; t++){
	            T curr;
	            if (2*i+ix[t]<x && 2*j+iy[t]<y && 2*l+iz[t]<z  && 2*m+iw[t]<w)
		      curr=a[k][0][(2*i+ix[t])*w*z*y+(2*j+iy[t])*w*z+
                      (2*l+iz[t])*w
                      +(2*m+iw[t])];
	            else curr=unconstrainedMask;
	            if (curr==valueMask || 
		        (curr==gradientMask && a2[k][0][halfSystemIndex]!=valueMask)) 
	  	      a2[k][0][halfSystemIndex]=curr;
	          }
	          if (a2[k][0][halfSystemIndex]==unconstrainedMask) 
	            unconstrainedCount++;
	        }else if (k==constraintCh){
	          int cnt=0;
	          for (t=0; t<numCoordsUnitBox; t++){
	            if (2*i+ix[t]<x && 2*j+iy[t]<y && 2*l+iz[t]<z && 2*m+iw[t]<w){
		      if (a[maskCh][0][(2*i+ix[t])*w*z*y+(2*j+iy[t])*w*z+
                        (2*l+iz[t])*w+(2*m+iw[t])]== a2[maskCh][0][halfSystemIndex]){
		        cnt++;
		        a2[k][0][halfSystemIndex]+=a[k][0][(2*i+ix[t])*w*z*y+
                          (2*j+iy[t])*w*z+(2*l+iz[t])*w+(2*l+iw[t])];
		      }
	            }
	          }	        	       
	          a2[k][0][halfSystemIndex]/=cnt;
	        }else{ 
	          for (t=0; t<numCoordsUnitBox; t++){
	            T curr=0;
	            if ((2*i+ix[t]<x) && (2*j+iy[t]<y) && (2*l+iz[t]<z) && (2*m+iw[t]<w))
		      curr=a[k][0][(2*i+ix[t])*w*z*y+(2*j+iy[t])*w*z+(2*l+iz[t])*w+
                        (2*l+iw[t])];
	            a2[k][0][halfSystemIndex]+=curr;
	          }
	        } 
	      }
            }
          }     
        }
      }
    }
    return unconstrainedCount;
  }

  /** Setup - includes computing the RHS vector, which is returned */
  template<class T>
  bool
  PoissonSystem<T>::_setup( Array<T> &a, 
			    Vector &rhs, 
			    ChannelList &gradientCh, 
			    unsigned constraintCh, 
			    unsigned maskCh, 
			    bool computeDivergence, 
			    unsigned sourceCh, 
			    bool cleanMatrix,
			    bool multigrid)
  {
    halfSystem=NULL;
    unsigned i;
    long j; 
    dim= a.getDimension();
    dimension= gradientCh.vec.size();
    totalPts= 1;
    for( i= 0 ; i< dimension ; i++ )
      totalPts*= dim.vec[i];
    if( gradPtr!= NULL ) delete [] gradPtr;
    gradPtr= new T*[dimension];
    if( !warnCond( dimension== dim.vec.size(), "number of gradients must match array dimension!\n" ) )
      return false;

    for( i= 0 ; i< dimension ; i++ ){
      gradPtr[i]= &((*a[gradientCh.vec[i]])[0]);
      if( !warnCond( gradPtr[i]!= NULL, "Gradient channel out of range\n" ) )
	return false;
    }

    constraints= &((*a[constraintCh])[0]);
    if( !warnCond( constraints!= NULL, "Constraint channel out of range\n" ) )
      return false;
  
    maskPtr= &((*a[maskCh])[0]);
    if( !warnCond( maskPtr!= NULL, "Mask channel out of range\n" ) )
      return false;

    T * source = NULL;
    if( !computeDivergence ){
      source = &((*a[sourceCh])[0]);
      if( !warnCond( source!= NULL, "Potential channel out of range\n" ) )
	return false;
    }
    setupOffsetTable();
    recalculateRHS(source,constraints,rhs,computeDivergence);  
    if (multigrid){
      halfSystem=createHalfSystem(a, gradientCh, constraintCh, 
				  maskCh, computeDivergence, sourceCh, 
				  cleanMatrix);
    }
    
    return true;
  }
  
  template<class T>
  void
  PoissonSystem<T>::recalculateRHS(T* source, 
				   T* constraints_, 
				   Vector& rhs, 
				   bool computeDivergence){
    
    if(constraints_==NULL) constraints_=constraints;

    setupCaseTable();
    computeNeighborCounts( );

    if(!computeDivergence){
      if(!warnCond( source!= NULL, "Potential channel NULL\n"))
	return;
    }
    long j;
   
    for( j= 0; j < totalPts; j++ ){
      if(maskPtr[j]!=valueMask){
	rhs[j]= computeRHS(j, computeDivergence, source) + 
	  computeBCContribution(j, constraints_);
      }else{
	rhs[j] = constraints[j];
      }
    }
  }

  template<typename T>
  bool
  PoissonSystem<T>::setup( Array<T> &a, unsigned maskCh)
  {
    halfSystem=NULL;
    unsigned i;
    long j;
    int newMaskCh = -1;
    dim= a.getDimension();
    dimension= dim.vec.size();
    totalPts= 1;
    for( i= 0 ; i< dimension ; i++ )
      totalPts*= dim.vec[i];
  
    maskPtr= &((*a[maskCh])[0]);
    if( !warnCond( maskPtr!= NULL, 
		   "Mask channel out of range, using default mask\n" ) ){
      newMaskCh = a.addChannel( unconstrainedMask );
      maskPtr= &((*a[newMaskCh])[0]);
    }
    if(gradPtr!=NULL){
      delete [] gradPtr;
      gradPtr = NULL;
    }
    setupOffsetTable();
    setupCaseTable();
    computeNeighborCounts( );
    if(newMaskCh >= 0){
      a.deleteChannel(newMaskCh);
      maskPtr=NULL;
    }
    return true;
  }

  /**Convenience method for computing the contribution of the boundary
     conditions (i.e. constrained pixel values or gradient values) to
     the RHS vector. */
  template<class T>
  inline T
  PoissonSystem<T>::computeBCContribution(long index, T* constraints)
  {
    unsigned currDim = 0;
    long currOffset;
    T contribution = 0;
    unsigned nonZero = neighbourConfig[index];

    for(;nonZero!=0 && currDim < dimension; nonZero >>= 2, currDim++){
      currOffset = offsetTable[currDim];
      if ((nonZero&3)==1)
	contribution += gradPtr[currDim][index];
      else
	if(maskPtr[index + currOffset]==valueMask){ 
	  contribution += constraints[index + currOffset];
	  setBit(neighbourConfig[index],2*currDim,false);
	}
	
      if ((nonZero&3)==2)
	contribution -= gradPtr[currDim][index];
      else
	if(maskPtr[index - currOffset]==valueMask){
	  contribution += constraints[index - currOffset];
	  setBit(neighbourConfig[index],2*currDim+1,false);
	}
    }
    return contribution;
  }

  /**convenience method for computing the RHS contribution from either
     the potential function or the divergence of the gradient*/
  template<class T>
  inline T 
  PoissonSystem<T>::computeRHS(long index, bool computeDivergence, T * source)
  {
    if(!computeDivergence){ 
      errorCond(source!=NULL, "source channel cannot be null");
      T ret=source[index];
      return ret;
    }

    T div = 0;
    long offset;
    unsigned nonZero = neighbourConfig[index];
    for(unsigned dim= 0; nonZero!=0 && dim < dimension; nonZero>>=2, dim++){
      offset = offsetTable[dim]; 
      switch(nonZero & 1) {
      case 0: 
        div += gradPtr[dim][index] - gradPtr[dim][index - offset];
        break;
      case 1: 
        div += gradPtr[dim][index + offset] - gradPtr[dim][index];
        break;      
      }
    }
    return -div;
  }

  template<class T>
  Vector &
  PoissonSystem<T>::rightMultiply(const Vector &v, Vector &res) const
  {
    double h;
    unsigned config, nbrs, num;
    unsigned char j;

    int c1=0,c2=0;
    int nz=0;

    for(unsigned long i=0; i<totalPts; i++){
      if(maskPtr[i] == valueMask){
	res[i] = v[i];
	c1++;
      }else{
	config= ImageSpaceSystem::nonZeroes[i];
	nbrs = neighbourConfig[i];
	num= ImageSpaceSystem::supportSize[nbrs];
	h = ImageSpaceSystem::supportSize[config] * v[i];
	for(j=1; j<= num; j++)
	  h-= v[i + ImageSpaceSystem::elems[nbrs][j].first];
	res[i] = h;
	if (fabs(res[i])<1e-10) nz++;
	c2++;
      }
    }
    
    return res;
  }

  //remove comments and replace template params with specific classes
  template class PoissonSystem<float>;
  template class PoissonSystem<double>;
} /* namespace */

#endif /* IMAGESPACESYSTEMS_POISSONSYSTEM_C */

