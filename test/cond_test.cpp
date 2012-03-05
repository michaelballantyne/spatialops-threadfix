#include <spatialops/structured/FVStaggeredFieldTypes.h>
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

int main()
{
  const int nx=10, ny=12, nz=14;

  typedef SVolField Field;
  
  const MemoryWindow window( get_window_with_ghost<Field>(IntVec(nx,ny,nz),true,true,true) );

  Field a( window, NULL );
  Field b( window, NULL );
  Field c( window, NULL );
  Field d( window, NULL );
  
  std::vector<Field> vec = std::vector<Field>();

  int const max = nx * ny * nz;

  Field::iterator ia = a.begin();
  Field::iterator ib = b.begin();
  Field::iterator id = d.begin();
  for(int i = 0; ia!=a.end(); i++, ++ia, ++ib, ++id ){
    *ia = i;
    *ib = max - i;
    *id = i % 2;
  }

  //Cond<Nil> Final cases:
  c <<= cond(1.1);
  c <<= cond(a);
  c <<= cond(a + b);

  //Cond<Simple> Final cases:
  c <<= cond(false, 0.0)
            (1.1);
  c <<= cond(true, 0.0)
            (1.1);
  c <<= cond(false, 0.0)
            (a);
  c <<= cond(true, 0.0)
            (a);
  c <<= cond(false, 0.0)
            (a + b);
  c <<= cond(true, 0.0)
            (a + b);

  //Cond<Nebo> Final cases:
  c <<= cond(true, a)
            (1.1);
  c <<= cond(false, a)
            (1.1);
  c <<= cond(false, a)
            (a);
  c <<= cond(true, a)
            (a);
  c <<= cond(false, a)
            (a + b);
  c <<= cond(true, a)
            (a + b);

  //Cond<Simple> Simple only cases:
  c <<= cond(true, 0.0)
            (1.1);
  c <<= cond(false, 0.0)
            (1.1);
  c <<= cond(false, 0.0)
            (true, 2.0)
            (1.1);
  c <<= cond(false, 0.0)
            (false, 2.0)
            (1.1);

  //Cond<Nil> bool Test cases:
  c <<= cond(true, 0.0)
            (1.1);
  c <<= cond(false, 0.0)
            (1.1);
  c <<= cond(true, a)
            (1.1);
  c <<= cond(false, a)
            (1.1);
  c <<= cond(true, a + 1)
            (1.1);
  c <<= cond(false, a + 1)
            (1.1);

  //Cond<Nil> NBE Test cases:
  c <<= cond(d > 0, 0.0)
            (1.1);
  c <<= cond(d > 0, a)
            (1.1);
  c <<= cond(d > 0, a + 1)
            (1.1);

  //Cond<Simple> bool Test cases:
  c <<= cond(false, 0.01)
            (true, 0.0)
            (1.1);
  c <<= cond(false, 0.01)
            (false, 0.0)
            (1.1);
  c <<= cond(false, 0.01)
            (true, a)
            (1.1);
  c <<= cond(false, 0.01)
            (false, a)
            (1.1);
  c <<= cond(false, 0.01)
            (true, a + 1)
            (1.1);
  c <<= cond(false, 0.01)
            (false, a + 1)
            (1.1);

  //Cond<Simple> NBE Test cases:
  c <<= cond(false, 0.01)
            (d > 0, 0.0)
            (1.1);
  c <<= cond(false, 0.01)
            (d > 0, a)
            (1.1);
  c <<= cond(false, 0.01)
            (d > 0, a + 1)
            (1.1);

  //Cond<Cond> bool Test cases:
  c <<= cond(false, b + 0.01)
            (true, 0.0)
            (1.1);
  c <<= cond(false, b + 0.01)
            (false, 0.0)
            (1.1);
  c <<= cond(false, b + 0.01)
            (true, a)
            (1.1);
  c <<= cond(false, b + 0.01)
            (false, a)
            (1.1);
  c <<= cond(false, b + 0.01)
            (true, a + 1)
            (1.1);
  c <<= cond(false, b + 0.01)
            (false, a + 1)
            (1.1);

  //Cond<Cond> NBE Test cases:
  c <<= cond(false, b + 0.01)
            (d > 0, 0.0)
            (1.1);
  c <<= cond(false, b + 0.01)
            (d > 0, a)
            (1.1);
  c <<= cond(false, b + 0.01)
            (d > 0, a + 1)
            (1.1);

  print(c);
  
  return 0;
}
