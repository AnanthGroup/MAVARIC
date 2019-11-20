/* This file contains somewhat miscellaneous functions that are necessary to run
   Hamiltonian.h. */

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <valarray>
#include <complex>
#include <vector>
#include <fstream>
#include <algorithm>

using namespace std;

/* Given a nested vector of type T, construct_tensor creates a dim1 X dim2 X dim3 tensor with
   each element initialized to zero. */
template <typename T>
void construct_tensor(vector<vector<vector<T> > > &tensor, int dim1, int dim2, int dim3){

  T zero = 0.0;
  tensor.resize(dim1);

  for(int i=0; i<dim1; i++){
    tensor[i].resize(dim2);
    for(int j=0; j<dim2; j++){
      tensor[i][j].resize(dim3,zero);
    }
  } 
}

/* Given a nested vector of type T, construct_tensor creates a dim1 X dim2 tensor with
   each element initialized to zero. */
template <typename T>
void construct_tensor(vector<vector<T> > &tensor, int dim1, int dim2){

  tensor.resize(dim1);
  T one = 0.0;

  for(int i=0; i<dim1; i++){
    tensor[i].resize(dim2,one);
  }
}

/* Standard implementation of kronecker delta function. */
inline double kronecker_delta(int i, int j){

  if(i == j) {return 1.0;}
  else{return 0.0;}
}

/* Given two square matricies of type T, m1 and m2, multiply them and store in result.
   Note: m1 and m2 MUST be square and have the same dimensions. */
template <typename T>
void matrix_mult(const vector<vector<T> > &m1,const vector<vector<T > > &m2,
 vector<vector<T > > &result){

  int dim = m1.size();
  T sum = 0.0;

  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      sum = 0.0;
         
      for(int k=0; k<dim; k++){
        sum = sum + m1[i][k]*m2[k][j];
      } 
      result[i][j] = sum;
    }
  }
} 

/* Given m, a nested vector of type T, initialize m to a dim X dim identity matrix. */
template<typename T>
void set_identity(vector<vector<T> >& m, int dim){

  m.resize(dim);

  for(int i=0; i<dim; i++){
    m[i].resize(dim);
  }

  for(int i=0; i<dim; i++){
    m[i][i] = 1;
  }
}

/* Given m, a dim X dim square matrix of type T, return the trace of m. */
template <typename T>
T trace(const vector<vector<T> > &m, int dim){

  T sum = 0;
  for(int i=0; i<dim; i++){
    sum = sum + m[i][i];
  }

  return sum;
}

/* Given v and w, two dim dimensional vectors, return the dot product. */
template <typename T>
T dot(const valarray<T> &v, const valarray<T> &w, int dim){

  T sum = (v * w).sum();
  return sum;
}

/* Compute the histogram of vector X and store the results in a file
 * callsed "fileName." The number of bins is specified by the parameters
 * bins. */
template <typename T>
void histogram(string fileName, int bins, T X){

  vector<double> ordered_X;
  int X_size = X.size();

  for(int i=0; i<X_size;i++){
    ordered_X.push_back(X[i]);
  }

  sort(ordered_X.begin(),ordered_X.end());

  double min = ordered_X[0];
  double max = ordered_X[X_size-1];
  double width = (max - min)/bins;

  double average = 0;

  for(int i=0; i<X_size; i++){
    average = average + ordered_X[i];
  }

  average = average/(X_size);

  double std = 0;

  for(int i=0; i<X_size; i++){
    std = std + (ordered_X[i] - average)*(ordered_X[i] - average);
  }

  std = sqrt(std/(X_size-1));

  double skew = 0;
  for(int i=0; i<X_size; i++){
    skew = skew + (ordered_X[i] - average)*(ordered_X[i] - average)
      *(ordered_X[i] - average)/(std*std*std);
  }

  skew = skew/(X_size);

  vector<double> finalHistogram (bins,0);

  int j = 1;
  for(int i=0; i<X_size; i++){
    if(ordered_X[i] <= (min + j*width)){
      finalHistogram[j-1] = finalHistogram[j-1] + 1;
    }
    else{j = j+1;}
                                                    
  }
 ofstream histogramOut;
  histogramOut.open(fileName);

  if(!histogramOut.is_open()){
    cout << "Failed to open histogram file." << endl;
  }
  else { /* no statement */ }

  double mode = finalHistogram[0];
  double largestValue = finalHistogram[0];

  for(int i=0; i<bins; i++){
    histogramOut << double(min + i*width) << " " << finalHistogram[i]/(X_size) << endl;

    if(finalHistogram[i] > largestValue){
      largestValue = finalHistogram[i];
      mode = (min + i*width);
    }
    else { /* no statement */ }
  }

  histogramOut.close();

  /* Print Statistics about the distribution. */

  cout << endl;
  cout << "\t--- Histogram Statistics ---" << endl;
  cout << "\t    Mode: " << mode << endl;
  cout << "\t    Average: " << average << endl;
  cout << "\t    Standard Deviation: " << std << endl;
  cout << "\t    Skew: " << skew << endl << endl;
}

/* Return the sign of x. */
inline double sign(double x){
  if(x < 0) {return -1;}
  else {return 1;}
}

#endif
