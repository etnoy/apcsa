#ifndef _CRS_H_
#define _CRS_H_

#include <omp.h>
#include <map>
#include <iostream>
#include <assert.h>

using namespace std;

class CRSmat {

private:

    map<int,double>* row_;

    int nrows_;
    int ncols_;
    int nnz_;



public:

    CRSmat() : nrows_(0), ncols_(0), nnz_(0), row_(0) { }

    ~CRSmat(){ delete [] row_;   }


    void beginAssembly(int n, int m){
        nrows_ = m;
        ncols_ = n;
        row_ = new map<int, double>[m];      assert( row_ != 0);
    }


    void completeAssembly(){
        nnz_ = 0;
        for (int i = 0; i < nrows_; i++) 
            nnz_ += row_[i].size();
        cout << "(" << 100.0*((double)nnz_/(double)nrows_)/(double)ncols_ 
            << "% non-zero entries in the matrix)" << endl << endl;
    }



    void setAij(int i, int j, double value){
        row_[i][j] = value;
    }


    double getAij(int i, int j) const{
        map<int,double>::iterator row_iter = row_[i].find(j);
        if (row_iter != row_[i].end())
            return (row_iter->second);
        else
            return (0.0);
    }

    
    // ------------------- START IMPLEMENTATION ---------------------------


    // c1)

    void mult_serial(const double u[], double v[]) const{



    }



    // c2)

    void mult_OMP(const double u[], double v[]) const{



    }



    // d)

    void mult_trans_OMP(const double u[], double v[]) const{



    }


    // ------------------- END IMPLEMENTATION ---------------------------


    void transpose(CRSmat & B) const{
        for (int i = 0; i < nrows_; i++){
            map<int,double>::iterator rowi;
            for (rowi = row_[i].begin(); rowi != row_[i].end(); ++rowi)
                B.setAij(rowi->first,i,getAij(i,rowi->first)); 
        } 
    }



    void     MatrixDump() const{
        for (int i = 0; i < nrows_; i++){
            map<int,double>::iterator rowi;
            for (rowi = row_[i].begin(); rowi != row_[i].end(); ++rowi)
                cout << "a(" << i << "," << rowi->first <<") = " <<  rowi->second << endl;
        } 
    }

    /*   void     dumpMatlab(char filename[]) const{ */
    /*     fstream outputFile(filename,ios::out); */
    /*     if (!outputFile){ */
    /*       cerr << "\007Hmmm! unable to create file " << filename << endl; */
    /*       exit(1);  */
    /*     } */
    /*     outputFile.setf(ios::scientific); outputFile.precision(16); */
    /*     outputFile << "D = [\n";   */
    /*     for (int i = 0; i < nrows_; ++i){ */
    /*       map<int,double>::iterator rowi; */
    /*       for (rowi = row_[i].begin(); rowi != row_[i].end(); ++rowi){ */
    /* 	outputFile << i+1 << " " << rowi->first+1 << " "; */
    /* 	outputFile.width(25); outputFile << rowi->second << endl; */
    /*       } */
    /*     } */
    /*     outputFile << "];\n"; */

    /*     char message[80]; */
    /*     sprintf(message,"A = spalloc(%d,%d,%d);\n",nrows_,ncols_,nnz_); */
    /*     outputFile << message; */
    /*     sprintf(message,"A = sparse(D(:,1),D(:,2),D(:,3));\n"); */
    /*     outputFile << message; */

    /*     outputFile.close(); */
    /*   } */

};



#endif



