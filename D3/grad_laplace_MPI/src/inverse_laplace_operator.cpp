#include "../header/inverse_laplace_operator.hpp"

#ifdef __cplusplus
extern "C" {
#endif

void inverse_laplace_operator_(double *out,double *in,double *sigma,int *L,int *V)
  {inverse_laplace_operator(out,in,*sigma,*L,*V);}

#ifdef __cplusplus
}
#endif
