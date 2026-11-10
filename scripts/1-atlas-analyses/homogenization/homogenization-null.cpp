#include <Rcpp.h>
#include <random>
#include <vector>
#include <cstdint>
#include <numeric>
using namespace Rcpp;

struct Community {
  int n, s, words;
  std::vector<uint64_t> bits;
  std::vector<int> richness;
  Community(IntegerMatrix x): n(x.nrow()), s(x.ncol()), words((s+63)/64),
    bits(n*words,0), richness(n,0) {
    for(int j=0;j<s;++j) for(int i=0;i<n;++i) {
      if(x(i,j)!=0 && x(i,j)!=1) stop("Nonbinary matrix");
      if(x(i,j)) {bits[i*words+j/64] |= uint64_t(1)<<(j%64); ++richness[i];}
    }
  }
  bool has(int i,int j) const {return bits[i*words+j/64] & (uint64_t(1)<<(j%64));}
  void change(int i,int j,bool gain) {
    if(has(i,j)==gain) stop("Invalid transition or duplicated selection");
    bits[i*words+j/64] ^= uint64_t(1)<<(j%64);
    richness[i] += gain?1:-1;
  }
  double beta() const {
    if(n<2) stop("At least two cells required");
    for(int i=0;i<n;++i) if(!richness[i]) stop("Empty simulated community at 1-based cell %d",i+1);
    double sum=0;
    for(int i=1;i<n;++i) for(int k=0;k<i;++k) {
      int shared=0;
      for(int w=0;w<words;++w) shared+=__builtin_popcountll(bits[i*words+w]&bits[k*words+w]);
      int united=richness[i]+richness[k]-shared;
      sum += double(united-shared)/united;
    }
    return sum/(double(n)*(n-1)/2);
  }
};
// [[Rcpp::export]]
double null_spatial_beta(IntegerMatrix x) {return Community(x).beta();}

// Uniform partial Fisher-Yates sampling without replacement. Pools always
// contain fixed baseline eligible cells; their order may persist between draws.
void transitions(Community &z,std::vector<int> &pool,int count,int species,bool gain,
                 std::mt19937 &rng) {
  if(count<0 || count>(int)pool.size()) stop("Infeasible transition count");
  for(int q=0;q<count;++q) {
    std::uniform_int_distribution<int> pick(q,pool.size()-1);
    int r=pick(rng); std::swap(pool[q],pool[r]); z.change(pool[q],species,gain);
  }
}
// [[Rcpp::export]]
List simulate_spatial_null(IntegerMatrix earlier,IntegerMatrix later,int nsim,
                          int seed,bool debug=false) {
  if(nsim<1 || earlier.nrow()!=later.nrow() || earlier.ncol()!=later.ncol()) stop("Invalid simulation dimensions/count");
  Community base(earlier), endpoint(later);
  base.beta(); endpoint.beta(); // Observed communities must be nonempty too.
  int n=base.n,s=base.s;
  std::vector<std::vector<int> > on(s),off(s);
  std::vector<int> loss(s,0),gain(s,0),before(s,0),after(s,0);
  for(int j=0;j<s;++j) for(int i=0;i<n;++i) {
    if(earlier(i,j)) {on[j].push_back(i);++before[j];} else off[j].push_back(i);
    if(later(i,j)) ++after[j];
    if(earlier(i,j)&&!later(i,j)) ++loss[j];
    if(!earlier(i,j)&&later(i,j)) ++gain[j];
  }
  int expected=std::accumulate(after.begin(),after.end(),0);
  std::mt19937 rng(seed); NumericVector beta(nsim);
  IntegerMatrix last(debug?n:0,debug?s:0);
  for(int b=0;b<nsim;++b) {
    Community z=base;
    for(int j=0;j<s;++j) {
      transitions(z,on[j],loss[j],j,false,rng);
      transitions(z,off[j],gain[j],j,true,rng);
      if(before[j]-loss[j]+gain[j]!=after[j]) stop("Endpoint occupancy invariant failed");
    }
    if(std::accumulate(z.richness.begin(),z.richness.end(),0)!=expected) stop("Richness total invariant failed");
    try { beta[b]=z.beta(); } catch(std::exception &e) {
      stop("Draw %d (chunk seed %d): %s",b+1,seed,e.what());
    }
    if(debug && b==nsim-1) for(int j=0;j<s;++j) for(int i=0;i<n;++i) last(i,j)=z.has(i,j);
  }
  return List::create(_["beta_later"]=beta,_["debug_last"]=last);
}
