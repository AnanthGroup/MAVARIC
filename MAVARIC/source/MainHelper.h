/* Author: Elliot Eklund 
 * Contains helper function of main.cpp */

#ifndef MAINHELPER_H
#define MAINHELPER_H

#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <vector>

using namespace std;

/* Read in the parameters from the InputFile "fileName". Return
 * these parameters in a vector.*/
vector<double> get_parameters(string fileName){

  vector<double> params;

  ifstream myFile;
  myFile.open(fileName);

  if(!myFile.is_open()) {
    cout << "Could not open file" << endl;
  }

  string myLine;

  while(getline(myFile, myLine)){

    if(!myLine.empty()){
      size_t foundAt = myLine.find(":"); //location of colon
      string goodStuff = myLine.substr(foundAt + 1, myLine.size());
      params.push_back(atof(goodStuff.c_str()));
    }
    else { /* no statement */ }
  }

  myFile.close();
  return params;
}

/* corruption_check returns true if fileName is corrupt.
 * A file is corrupt if it meets any of the following criteria.
 * 1. The number of lines in fileName does not match correct_num_lines.
 * 2. The format DESCRIPTION:VALUE is altered.
 * 3. Any white space or non-numerical characters are found after the colon.*/
bool corruption_check(string fileName, int correct_num_lines ){

  ifstream myFile;
  myFile.open(fileName);
 
  if(!myFile.is_open()) {
    cout << "Could not open file" << endl;
  }
 
  string myLine;
  int num_lines = 0;
  bool file_corrupt = false;

  while(getline(myFile, myLine)){

    num_lines += 1;

    if(!myLine.empty()){
 
      int foundAt = myLine.find(":"); //location of colon
      string goodStuff = myLine.substr(foundAt + 1, myLine.size());

      if((foundAt <= 0) || (goodStuff.length()==0)){
        file_corrupt = true;
      }

      if(goodStuff.find(' ') != -1){
        file_corrupt = true;
      }

      char* p;
      double converted = strtod(goodStuff.c_str(), &p);
    
      if (*p) {
        file_corrupt = true;
      }
    }
  }

  if(num_lines != correct_num_lines){
    file_corrupt = true;
  }

  myFile.close();

  if(file_corrupt){
    cout << "ERROR: File " << fileName << " is corrupt. Aborting calculation." << endl;
  }

  return file_corrupt;

}

/* Return true if x is a positive number. Otherwise, return false. */
bool check_positive(double x, string fileName, string varName){

  if (x < 0){
    cout << "ERROR: In file " + fileName + ", " + varName + " must be a positive number."
      " Aborting calculation." << endl;
    return false;
  }
  else{
    return true;
  }
}
 
/* Return true if x is a positive integer. Otherwise, return false. */
bool check_positive_int(double x, string fileName, string varName){

  if ((x < 0) || (x != floor(x))){
    cout << "ERROR: In file " + fileName + ", " + varName + " must be a positive integer."
      " Aborting calculation." << endl;
    return false;
  }
  else{
    return true;
  } 
}
 
 /* Return true if x i either 0 or 1. Otherwise, return false. */
bool check_binary(double x, string fileName, string varName){

  if((x!=0) && (x!=1)){
    cout << "ERROR: In file " + fileName + ", " + varName + " must be either 1 or 0."
      " Aborting calculation." << endl;
    return false;
  }

  else{
    return true;
  } 
}

/* Calling this function assumes that SysParameters is not corrupt and passed 
 * the call to corruption_check. */
int SysParams_checker(vector<double> & params){

  int result = 0;
  ifstream myFile;
  string fileName = "InputFiles/SystemParameters";
  myFile.open(fileName);
 
  if(!myFile.is_open()) {
    cout << "Could not open file" << endl;
  }

  string myLine;
    
  while(getline(myFile,myLine)){

    int foundAt = myLine.find(":"); //location of colon
    string goodStuff = myLine.substr(foundAt + 1, myLine.size());
    params.push_back(atof(goodStuff.c_str()));
  
  }     

  myFile.close();

  /* Check Mass */
  if(!check_positive(params[0],fileName,"Mass")){
    result = -1; 
  }
  
  /* Check Beads */
  if(!check_positive_int(params[1],fileName,"Beads")){
      result = -1;
  }

  /* Check Temperature */
  if(!check_positive(params[2],fileName,"Temperature")){
    result = -1;
  }

  /* Check MC Step Size */
  if(!check_positive(params[3],fileName,"MC Step Size")){
    result = -1;
  }

  /* If no errors have been caught, return 0. */
  return result;
}

/* Calling this function assumes that ElecParameters is not corrupt and passed
 * the call to corruption_check. */
int ElecParams_checker(vector<double> &params){

  int result = 0;

  ifstream myFile;
  string fileName = "InputFiles/ElecParameters";
  myFile.open(fileName);
  
  if(!myFile.is_open()) {
    cout << "Could not open file" << endl;
  } 

  string myLine;

  while(getline(myFile,myLine)){

    int foundAt = myLine.find(":"); //location of colon
    string goodStuff = myLine.substr(foundAt + 1, myLine.size());
    params.push_back(atof(goodStuff.c_str()));
  
  } 

  myFile.close();

  /* Check Number of States */
  if(!check_positive_int(params[0],fileName,"States")){
    result = -1;
  } 
  
  /* Check MC Step Size */
  if(!check_positive(params[1],fileName,"MC Step Size")){
    result = -1; 
  };
  
    
  /* If no errors have been caught, return 0. */
  return result;  
}

/* Calling this function assumes that MonteCarlo is not corrupt and passed 
 * the call to corruption check. */
int MCParams_checker(vector<double> &params){

  int result = 0;
  ifstream myFile;
  string fileName = "InputFiles/MonteCarlo";
  myFile.open(fileName);
 
  if(!myFile.is_open()) {
    cout << "Could not open file" << endl;
  }

  string myLine;
  
  while(getline(myFile,myLine)){
        
    int foundAt = myLine.find(":"); //location of colon
    string goodStuff = myLine.substr(foundAt + 1, myLine.size());
    params.push_back(atof(goodStuff.c_str()));
    
  }
  
  myFile.close();
  
  /* Check that Run MC is binary. */
  if(!check_binary(params[0],fileName,"Run MC")){
    result = -1;
  }
  
  /* Check that MC Moves ia positive integer. */
  if(!check_positive_int(params[1],fileName,"MC Moves")){
    result = -1;
  }

  /* Check that Estimator Rate ia positive integer. */
  if(!check_positive_int(params[2],fileName,"Estimator Rate")){
    result = -1;
  }

  /* Check that Estimator Rate is smaller than MC Moves */
  if(params[2] >= params[1]){
    cout << "ERROR: In file MonteCarlo, Estimator Rate must be less than MC Moves."
      " Aborting calculation." << endl;
    return -1; 
  }

  /* Check that MC Moves is divisible by Estimator Rate. */
  if(int(params[1]) % int(params[2]) != 0){
    cout << "ERROR: In file MonteCarlo, MC Moves must be divisible by Estimator Rate."
      " Aborting calculation." << endl;
    return -1;
  }

  /*Check that Save PSV is binary */
  if(!check_binary(params[3],fileName,"Save PSV")){
    result = -1;
  }

  /*Check that Save MC Data is binary */
  if(!check_binary(params[4],fileName,"Save MC Data")){
    result = -1;
  }

  /*Check that Read PSV is binary */
  if(!check_binary(params[5],fileName,"Read PSV")){
    result = -1;
  }

  /*Check that Read MC Data is binary */
  if(!check_binary(params[6],fileName,"Read PSV")){
    result = -1;
  }

  return result;
}

/* Calling this function assumes that Sampling is not corrupt and passed 
 * the call to corruption check. */
int SampParams_checker(vector<double> &params){

  int result = 0;
  ifstream myFile;
  string fileName = "InputFiles/Sampling";
  myFile.open(fileName);
 
  if(!myFile.is_open()) {
    cout << "Could not open file" << endl;
  }
        
  string myLine;
  
  while(getline(myFile,myLine)){
    
    int foundAt = myLine.find(":"); //location of colon
    string goodStuff = myLine.substr(foundAt + 1, myLine.size());
    params.push_back(atof(goodStuff.c_str()));
  
  }
  
  myFile.close();

  
  /* Check that Run Sampling is binary. */
  if(!check_binary(params[0],fileName,"Run Sampling")){
    result = -1;
  }

  /* Check that Number of Trajectories is a positive integer. */
  if(!check_positive_int(params[1],fileName,"Number of Trajectories")){
    result = -1;
  }

  /* Check that Decorrelation Length is a positive integer. */
  if(!check_positive_int(params[2],fileName,"Decorrelation Length")){
    result = -1;
  }


  /* Check that Save Sampled Trajectories is binary. */
  if(!check_binary(params[3],fileName,"Save Sampled Trajectories")){
    result = -1;
  }
/* Check that Histogram Positions is binary. */
  if(!check_binary(params[4],fileName,"Histogram Positions")){
    result = -1;
  }

  /* Check that Number of Bins is a positive integer. */
  if(!check_positive_int(params[5],fileName,"Number of Bins")){
    result = -1;
  }

  /* Check that Read PSV is binary. */
  if(!check_binary(params[6],fileName,"Read PSV")){
    result = -1;
  }


  return result;

}

/* Calling this function assumes that Dynamics is not corrupt and passed 
 * the call to corruption check. */
int DynParams_checker(vector<double> &params){

  int result = 0;
  ifstream myFile;
  string fileName = "InputFiles/Dynamics";
  myFile.open(fileName);
 
  if(!myFile.is_open()) {
    cout << "Could not open file" << endl;
  }
  
  string myLine; 
  
  while(getline(myFile,myLine)){
  
    int foundAt = myLine.find(":"); //location of colon
    string goodStuff = myLine.substr(foundAt + 1, myLine.size());
    params.push_back(atof(goodStuff.c_str()));

  }

  myFile.close();

  /* Check that Run Dynamics is binary. */
  if(!check_binary(params[0],fileName,"Run Dynamics")){
    result = -1;
  }

  /* Check that Run Time is a positive number. */
  if(!check_positive(params[1],fileName,"Run Time")){
    result = -1;
  }

  /* Check that Time Step is a positive number. */
  if(!check_positive(params[2],fileName,"Time Step")){
    result = -1;
  }

  double intpart; //used to check for integer

  /* Check that Run Time is divisible by Time Step */
  if(!(modf(params[1]/params[2],&intpart) == 0.0)){
    cout << "ERROR: In file Dynamics, Run Time must be divisible by Time Step."
      " Aborting calculation." << endl;
    return -1;
  }

  /* Check that Read Trajectories is binary. */
  if(!check_binary(params[3],fileName,"Read Trajectories")){
    result = -1;
  }

  /* Check that Check Energy Conservation is binary. */
  if(!check_binary(params[4],fileName,"Check Energy Conservation")){
    result = -1;
  }

  /* Check that Conservation Tolerance is a positive number. */
  if(!check_positive(params[5],fileName,"Conservation Tolerance")){
    result = -1;
  }

  return result;
}

/* Each of the five input vectors well be filled with their corresponding values after
 * input_file_handler is called. These vectors will only be returned if the function does
 * not find any corrupt files or errors.*/
int input_file_handler(vector<double> &sysParams,vector<double> &elecParams,
      vector<double> &MCParams,vector<double> &sampParams, vector<double> &dynParams){

  string ElecFile = "InputFiles/ElecParameters";
  string SysFile = "InputFiles/SystemParameters";
  string MCFile = "InputFiles/MonteCarlo";
  string SampFile ="InputFiles/Sampling";
  string DynFile = "InputFiles/Dynamics";

  int result = 0;

  /* First, check if any of the input files are corrupt */
  if(corruption_check(ElecFile,2)){
    result = -1;
  }

  if(corruption_check(SysFile,4)){
    result = -1;
  }

  if(corruption_check(MCFile,7)){
    result = -1;
  }

  if(corruption_check(SampFile,7)){
    result = -1;
  }

  if(corruption_check(DynFile,6)){
    result = -1;
  }

  /* If an input file is corrupt, do not proceed exit
   * input_file_handler and return -1 */
  if(result == -1){
    return -1;
  }


  /* Check that the parameters in each file are acceptable. */
  if(SysParams_checker(sysParams) == -1){
    result = -1;
  };

  if(ElecParams_checker(elecParams) == -1){
    result = -1;
  };

  if(MCParams_checker(MCParams) == -1){
    result = -1;
  }

  if(SampParams_checker(sampParams) == -1){
    result = -1;
  }

  if(DynParams_checker(dynParams) == -1){
    result = -1;
  }

  return result;
}

#endif
