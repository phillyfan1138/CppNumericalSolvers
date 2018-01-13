#ifndef NELDERMEADSOLVER_H_
#define NELDERMEADSOLVER_H_
#include "FunctionalUtilities.h"

namespace neldermead {

    template <typename T>
    auto getRowAsVector(const std::vector<T>& vM, int col, int rowLength){
        auto first = myVec.begin() + col*rowLength;
        auto last = first+rowLength;
        std::vector<T> newVec(first, last);
        return newVec;
    }
    template <typename T>
    auto getItemAt(int row, int col, int rowLength, const std::vector<T>& vM){
        return vM[col*rowLength+row];
    }
    template <typename T>
    void setItemAt(int row, int col, int rowLength, const T& item, std::vector<T>&& vM){
        vM[col*rowLength+row]=item;
    }

    template<typename T>
    std::vector<T> makeInitialSimplex(const std::vector<T> &x) { //Array is std::vector...please
        int DIM = x.size();
        std::vector<T> s=std::vector<T>(DIM*(DIM+1));
        int numCol=DIM+1;
   // MatrixType s = MatrixType::Zero(DIM, DIM + 1);
        for (int c = 0; c < numCol; ++c) {
            for (int r = 0; r < DIM; ++r) {
                setItemAt(r, c, numCol, x[r], s);
                if (r == c - 1) {
                    if (x[r] == 0) {
                        setItemAt(r, c, numCol, .00025, s);
                    } else {
                        setItemAt(r, c, numCol, (1 + 0.05) * x[r], s);
                    }
                }
            }
        }
        return s;
    }

    template<typename ObjFun, typename T>
    auto optimize(const ObjFun& objFun, const std::vector<T>& x){
        auto simplex=makeInitialSimplex(x);
        
    }

}
