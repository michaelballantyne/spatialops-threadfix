#include <spatialops/structured/FVStaggeredTypes.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/FieldExpressions.h>
#include <spatialops/FieldExpressionsExtended.h>
#include <spatialops/FieldReductions.h>

#include <iostream>
#include <vector>

using namespace SpatialOps;
using namespace structured;

int print_length = 8;

template<typename Field>
void print(Field const & given) {
  typename Field::const_iterator ig = given.begin();
  
  for(int i = 0; i < print_length; i++) {
    std::cout << *ig << std::endl;
    ++ig;
  };
  
  std::cout << std::endl;
};

template<typename Field>
void print_all(Field const & given) {
  typename Field::const_iterator ig = given.begin();
  
  for(int i = 0; ig != given.end(); i++) {
    std::cout << *ig << std::endl;
    ++ig;
  };
  
  std::cout << std::endl;
};

double add_vec(std::vector<double> const & inputs) {
  double answer = 0;
  
  for(std::vector<double>::const_iterator ii = inputs.begin();
      ii != inputs.end();
      ++ii) {
    answer += *ii;
  };
  
  return answer;
};

double mult_vec(std::vector<double> const & inputs) {
  double answer = 1;
  
  for(std::vector<double>::const_iterator ii = inputs.begin();
      ii != inputs.end();
      ++ii) {
    answer *= *ii;
  };
  
  return answer;
};

int main()
{
  const int nx=10, ny=12, nz=14;

  typedef SVolField Field;
  
  const MemoryWindow window( get_window_with_ghost<Field>(IntVec(nx,ny,nz),true,true,true) );

  Field a( window, NULL );
  Field b( window, NULL );
  Field c( window, NULL );
  
  std::vector<Field> vec = std::vector<Field>();
  
  Field::iterator ia1 = a.begin();
  for(int i = 0; ia1!=a.end(); i++, ++ia1 ){
    *ia1 = i;
  }

  b = 2.0;

  DEFINE_$_ANONYMOUS_ARGUMENTS(Field);
  
  // example of what we would like to do:
  // c = b + a + sin(b);
  
//   c <<= -a;
//   print(c);
//   c <<= (3 + a) - 2;
//   print(c);
//   c <<= 3 + (a + a);
//   print(c);
//   c <<= (3 - (3 * $0))(a);
//   print(c);
//   c <<= (a + 3) - a;
//   print(c);
//   c <<= (a * (a + 1)) / (a - 1);
//   c <<= (a - (a * $0))(a);
//   c <<= ((a + a) * ((a + a) - $0))(a);
//   c <<= (($0 + 3) - 3)(a);
//   c <<= (($0 + a) - a)(3.0);
//   c <<= (($0 + (a + a)) - (a + a))(3.0);
//   c <<= (($0 * $0) - $0)(a + a);
//   c <<= ($0 - ($0 * $0))(a + a);
//   c <<= (($1 * $2) - $0)(3 * a, a, 3.0);
//   c <<= (($1 * $2) - $0)(3 * a, a)(3.0);
//   c <<= (($1 * $2) - $0)(3 * a)(a, 3.0);
//   c <<= (($1 * $2) - $0)(3 * a)(a)(3.0);
  //  c <<= BinOp(+, a, b);
  
  std::cout << fieldMax(a) << std::endl << std::endl;
  std::cout << fieldMin(a) << std::endl << std::endl;
  std::cout << fieldSum(a) << std::endl << std::endl;
  std::cout << fieldNorm(a) << std::endl << std::endl;
  std::cout << fieldNorm(b) << std::endl << std::endl;
  
//   vec.push_back(a);
//   vec.push_back(b);
  
//   print(a);
  
  // what we currently must do:
  Field::const_iterator ia = a.begin();
  for( Field::iterator ic=c.begin(); ic!=c.end(); ++ic, ++ia ){
    *ic = sin(*ia);
  }
  c += a;  // operator overload
  c += b;  // operator overload

  c *= a;
  c /= a;

  return 0;
}
