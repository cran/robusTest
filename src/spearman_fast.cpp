#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>

//' @useDynLib robusTest, .registration=TRUE
//' @importFrom Rcpp sourceCpp

using namespace Rcpp;
using namespace std;

class ArgComparator {
private:
  const NumericVector& _v;
public:
  ArgComparator(const NumericVector& v) : _v(v) {}
  bool operator()(size_t i1, size_t i2) {
    return _v[i1] < _v[i2];
  }
};


vector<size_t> argsort2(const NumericVector& v) {
  vector<size_t> idx(v.size());
  ArgComparator cmp(v);
  for(size_t i = 0; i < v.size(); i++) {
    idx[i] = i;
  }
  stable_sort(idx.begin(), idx.end(), cmp);
  return idx;
}

vector<size_t> invPerm(const vector<size_t>& sigma) {
  vector<size_t> ret(sigma.size());
  for(size_t i = 0; i < sigma.size(); i++) {
    ret[sigma[i]] = i;
  }
  return ret;
}

// [[Rcpp::export]]
List spearmanCore(const NumericVector& x, const NumericVector& y) {
  size_t n = x.size();
  vector<size_t> orderX = argsort2(x);
  vector<size_t> orderY = argsort2(y);
  vector<size_t> rankX = invPerm(orderX);
  vector<size_t> rankY = invPerm(orderY);
  vector<double> H(n);
  long long int sumR = 0;
  for(size_t i = 0; i < n; i++) {
    H[i] += (rankX[i]/double(n) - 1)*(rankY[i]/double(n) - 1) - 1;
    sumR += (2*rankX[i]-n+1)*(2*rankY[i]-n+1);
    int cpt = 0;
    for(size_t j = 0; j < n; j++) {
      H[orderX[j]] += double(cpt)/(n*n);
      H[orderY[i]] += double(cpt)/(n*n);
      if (y[orderX[j]] < y[orderY[i]]) cpt++;
    }
  }
  List ret;
  ret["sumR"] = sumR;
  ret["H"] = H;
  return ret;
}
