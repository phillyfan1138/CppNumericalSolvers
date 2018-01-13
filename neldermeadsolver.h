// CppNumericalSolver
#ifndef NELDERMEADSOLVER_H_
#define NELDERMEADSOLVER_H_
#include <cmath>
#include "FunctionalUtilities.h"


namespace cppoptlib {
  template <typename T>
  auto getColumnAsVector(const std::vector<T>& vM, int col, int numRows){
    auto first = vM.begin() + col*numRows;
    auto last = first+numRows;
    std::vector<T> newVec(first, last);
    return newVec;
  }
  template <typename T>
  auto getItemAt(int row, int col, int numRows, const std::vector<T>& vM){
    return vM[col*numRows+row];
  }
  template <typename T>
  auto setItemAt(int row, int col, int numRows, const T& item, std::vector<T>&& vM){
    vM[col*numRows+row]=item;
    return std::move(vM);
  }

  template<typename T>
  std::vector<T> makeInitialSimplex(const std::vector<T> &x) { 
    const int numRows = x.size();
    const int numCols=numRows+1;
    std::vector<T> s=std::vector<T>(numRows*numCols, 0.0);
    for (int c = 0; c < numCols; ++c) {
      for (int r = 0; r < numRows; ++r) {
        s=setItemAt(r, c, numRows, x[r], std::move(s));
        if (r == c - 1) {
          if (x[r] == 0) {
            s=setItemAt(r, c, numRows, .00025, std::move(s));
          } else {
            s=setItemAt(r, c, numRows, (1 + 0.05) * x[r], std::move(s));
          }
        }
      }
    }
    return s;
  }
  template<typename T>
  T getMaxCoef(const std::vector<T>& col){
    return futilities::reduce_to_single(col, [](const auto& prev, const auto& val, const auto& index){
      return prev>val?prev:val;
    });
  }
  template<typename T>
  T getMaxCoefOfDiff(std::vector<T>&& col1, std::vector<T>&& col2){
    return getMaxCoef(futilities::for_each_parallel(col1, [&](const auto& val, const auto& index){
      return fabs(val-col2[index]);
    }));
  }

/*
  template<typename ObjFunc, typename T >
  void minimize(ObjFunc &&objFunc, const std::vector<T> &x, const int maxIter) {
    const T rho = 1.;    // rho > 0
    const T xi  = 2.;    // xi  > max(rho, 1)
    const T gam = 0.5;   // 0 < gam < 1
    const int numRows = x.rows();
    const int numCols=numRows+1;
    // create initial simplex
    auto x0 = makeInitialSimplex(x);
    // compute function values
    auto f=futilities::for_each_parallel(0, numCols, [&](const auto& ind){
      return objFunc(getColumnAsVector(x0, ind, numRows));
    });
    auto index=futilities::for_each_parallel(0, numCols, [&](const auto& ind){
      return ind;
    });
    sort(index.begin(), index.end(), [&](int a, int b)-> bool { return f[a] < f[b]; });
    int iter = 0;
    while (
      (iter < maxIter) && 
    ) {
      // conv-check
      T maxFn=futilities::reduce_to_single(f, [&](const auto& prev, const auto& val, const auto& ind){
        auto currDif=fabs(val-f[index[0]]);
        return prev>currDif?prev:currDif;
      });
      T maxCoef=futilities::reduce_to_single(index, [&](const auto& prev, const auto& val, const auto& ind){
        auto currDif=getMaxCoefOfDiff(getColumnAsVector(x0, val, numRows), getColumnAsVector(x0, index[0], numRows));
        return prev>currDif?prev:currDif;
      });

      const T tt1 = std::max(T(1.e-04), 10 * std::nextafter(f[index[0]], std::numeric_limits<T>::epsilon()) - f[index[0]]);
      const T tt2 = std::max(T(1.e-04), 10 * (std::nextafter(getMaxCoef(getColumnAsVector(x0, 0, numRows)) , std::numeric_limits<T>::epsilon())- getMaxCoef(getColumnAsVector(x0, 0, numRows)));


      //maxFn and maxCoef are used for stopping criteria...if gets to certain size, stop
      //compare maxFn to tt1 and maxCoef to tt2....if both are smaller than corresponding tt, then stop

     

      //////////////////////////

      // midpoint of the simplex opposite the worst point
      std::vector<T> x_bar = std::vector<T>(numRows);
      for (int i = 0; i < int(DIM); ++i) {
        x_bar += x0.col(index[i]);
      }
      x_bar /= Scalar(DIM);

      // Compute the reflection point
      const TVector x_r   = ( 1. + rho ) * x_bar - rho   * x0.col(index[DIM]);
      const Scalar f_r = objFunc(x_r);
      lastOp = SimplexOp::Reflect;

      if (f_r < f[index[0]]) {
        // the expansion point
        const TVector x_e = ( 1. + rho * xi ) * x_bar - rho * xi   * x0.col(index[DIM]);
        const Scalar f_e = objFunc(x_e);
        if ( f_e < f_r ) {
          // expand
          lastOp = SimplexOp::Expand;
          x0.col(index[DIM]) = x_e;
          f[index[DIM]] = f_e;
        } else {
          // reflect
          lastOp = SimplexOp::Reflect;
          x0.col(index[DIM]) = x_r;
          f[index[DIM]] = f_r;
        }
      } else {
        if ( f_r < f[index[DIM - 1]] ) {
          x0.col(index[DIM]) = x_r;
          f[index[DIM]] = f_r;
        } else {
          // contraction
          if (f_r < f[index[DIM]]) {
            const TVector x_c = (1 + rho * gam) * x_bar - rho * gam * x0.col(index[DIM]);
            const Scalar f_c = objFunc(x_c);
            if ( f_c <= f_r ) {
              // outside
              x0.col(index[DIM]) = x_c;
              f[index[DIM]] = f_c;
              lastOp = SimplexOp::ContractOut;
            } else {
              shrink(x0, index, f, objFunc);
              lastOp = SimplexOp::Shrink;
            }
          } else {
            // inside
            const TVector x_c = ( 1 - gam ) * x_bar + gam   * x0.col(index[DIM]);
            const Scalar f_c = objFunc(x_c);
            if (f_c < f[index[DIM]]) {
              x0.col(index[DIM]) = x_c;
              f[index[DIM]] = f_c;
              lastOp = SimplexOp::ContractIn;
            } else {
              shrink(x0, index, f, objFunc);
              lastOp = SimplexOp::Shrink;
            }
          }
        }
      }
      sort(index.begin(), index.end(), [&](int a, int b)-> bool { return f[a] < f[b]; });
      iter++;
      if (iter >= maxIter) {
        stop_condition = Status::IterationLimit;
      }
      else {
        stop_condition = Status::UserDefined; // if stopped in the callback in while()
      }
    } // while loop

    // report the last result
    objFunc.detailed_callback(this->m_current, lastOp, index[0], x0, f);
    x = x0.col(index[0]);
  }

  void shrink(MatrixType &x, std::vector<int> &index, std::vector<Scalar> &f, ProblemType &objFunc) {
    const Scalar sig = 0.5;   // 0 < sig < 1
    const int DIM = x.rows();
    f[index[0]] = objFunc(x.col(index[0]));
    for (int i = 1; i < DIM + 1; ++i) {
      x.col(index[i]) = sig * x.col(index[i]) + (1. - sig) * x.col(index[0]);
      f[index[i]] = objFunc(x.col(index[i]));
    }
  }

  // Need our own checker here to get rid of the gradient test used in other solvers
  template<typename T>
  Status checkConvergence(const Criteria<T> &stop, const Criteria<T> &current) {
    if ((stop.iterations > 0) && (current.iterations > stop.iterations)) {
      return Status::IterationLimit;
    }
    if ((stop.xDelta > 0) && (current.xDelta < stop.xDelta)) {
      return Status::XDeltaTolerance;
    }
    if ((stop.fDelta > 0) && (current.fDelta < stop.fDelta)) {
      return Status::FDeltaTolerance;
    }
    return Status::Continue;
  }
*/


}

#endif 
