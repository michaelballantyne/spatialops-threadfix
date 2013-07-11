#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>

#include <spatialops/Nebo.h>

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
  
  const int nghost = 1;
  const GhostDataRT ghost(nghost);
  const BoundaryCellInfo bc = BoundaryCellInfo::build<Field>(true,true,true);
  const MemoryWindow window( get_window_with_ghost(IntVec(nx,ny,nz),ghost,bc) );

  Field a( window, bc, ghost, NULL );
  Field b( window, bc, ghost, NULL );
  Field c( window, bc, ghost, NULL );
  Field d( window, bc, ghost, NULL );
  
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

  //Nebo Boolean Expression tests:
  c <<= nebo_cond(a == 0, 0.0)
                 (a == a, 0.0)
                 (0 == a, 0.0)
                 (a != 0, 0.0)
                 (a != a, 0.0)
                 (0 != a, 0.0)
                 (a < 0, 0.0)
                 (a < a, 0.0)
                 (0 < a, 0.0)
                 (a <= 0, 0.0)
                 (a <= a, 0.0)
                 (0 <= a, 0.0)
                 (a > 0, 0.0)
                 (a > a, 0.0)
                 (0 > a, 0.0)
                 (a >= 0, 0.0)
                 (a >= a, 0.0)
                 (0 >= a, 0.0)
                 (1.0);


  //Nebo Boolean Expression tests:
  c <<= nebo_cond(true && true, 0.0)
                 (true && a < 0, 0.0)
                 (true && !(a < 0), 0.0)
                 (a < 0 && true, 0.0)
                 (a < 0 && a < 0, 0.0)
                 (a < 0 && !(a < 0), 0.0)
                 (!(a < 0) && true, 0.0)
                 (!(a < 0) && a < 0, 0.0)
                 (!(a < 0) && !(a < 0), 0.0)
                 (true || true, 0.0)
                 (true || a < 0, 0.0)
                 (true || !(a < 0), 0.0)
                 (a < 0 || true, 0.0)
                 (a < 0 || a < 0, 0.0)
                 (a < 0 || !(a < 0), 0.0)
                 (!(a < 0) || true, 0.0)
                 (!(a < 0) || a < 0, 0.0)
                 (!(a < 0) || !(a < 0), 0.0)
                 (!true, 0.0)
                 (!(a < 0), 0.0)
                 (!(a < 0 || a < 0), 0.0)
                 (1.0);


  //Cond<Nil> Final cases:
  c <<= nebo_cond(1.1);
  c <<= nebo_cond(a);
  c <<= nebo_cond(a + b);

  //Cond<Simple> Final cases:
  c <<= nebo_cond(false, 0.0)
                 (1.1);
  c <<= nebo_cond(true, 0.0)
                 (1.1);
  c <<= nebo_cond(false, 0.0)
                 (a);
  c <<= nebo_cond(true, 0.0)
                 (a);
  c <<= nebo_cond(false, 0.0)
                 (a + b);
  c <<= nebo_cond(true, 0.0)
                 (a + b);

  //Cond<Nebo> Final cases:
  c <<= nebo_cond(true, a)
                 (1.1);
  c <<= nebo_cond(false, a)
                 (1.1);
  c <<= nebo_cond(false, a)
                 (a);
  c <<= nebo_cond(true, a)
                 (a);
  c <<= nebo_cond(false, a)
                 (a + b);
  c <<= nebo_cond(true, a)
                 (a + b);

  //Cond<Simple> Simple only cases:
  c <<= nebo_cond(true, 0.0)
                 (1.1);
  c <<= nebo_cond(false, 0.0)
                 (1.1);
  c <<= nebo_cond(false, 0.0)
                 (true, 2.0)
                 (1.1);
  c <<= nebo_cond(false, 0.0)
                 (false, 2.0)
                 (1.1);

  //Cond<Nil> bool Test cases:
  c <<= nebo_cond(true, 0.0)
                 (1.1);
  c <<= nebo_cond(false, 0.0)
                 (1.1);
  c <<= nebo_cond(true, a)
                 (1.1);
  c <<= nebo_cond(false, a)
                 (1.1);
  c <<= nebo_cond(true, a + 1)
                 (1.1);
  c <<= nebo_cond(false, a + 1)
                 (1.1);

  //Cond<Nil> NBE Test cases:
  c <<= nebo_cond(d > 0, 0.0)
                 (1.1);
  c <<= nebo_cond(d > 0, a)
                 (1.1);
  c <<= nebo_cond(d > 0, a + 1)
                 (1.1);

  //Cond<Simple> bool Test cases:
  c <<= nebo_cond(false, 0.01)
                 (true, 0.0)
                 (1.1);
  c <<= nebo_cond(false, 0.01)
                 (false, 0.0)
                 (1.1);
  c <<= nebo_cond(false, 0.01)
                 (true, a)
                 (1.1);
  c <<= nebo_cond(false, 0.01)
                 (false, a)
                 (1.1);
  c <<= nebo_cond(false, 0.01)
                 (true, a + 1)
                 (1.1);
  c <<= nebo_cond(false, 0.01)
                 (false, a + 1)
                 (1.1);

  //Cond<Simple> NBE Test cases:
  c <<= nebo_cond(false, 0.01)
                 (d > 0, 0.0)
                 (1.1);
  c <<= nebo_cond(false, 0.01)
                 (d > 0, a)
                 (1.1);
  c <<= nebo_cond(false, 0.01)
                 (d > 0, a + 1)
                 (1.1);

  //Cond<Nebo> bool Test cases:
  c <<= nebo_cond(false, b + 0.01)
                 (true, 0.0)
                 (1.1);
  c <<= nebo_cond(false, b + 0.01)
                 (false, 0.0)
                 (1.1);
  c <<= nebo_cond(false, b + 0.01)
                 (true, a)
                 (1.1);
  c <<= nebo_cond(false, b + 0.01)
                 (false, a)
                 (1.1);
  c <<= nebo_cond(false, b + 0.01)
                 (true, a + 1)
                 (1.1);
  c <<= nebo_cond(false, b + 0.01)
                 (false, a + 1)
                 (1.1);

  //Cond<Nebo> NBE Test cases:
  c <<= nebo_cond(false, b + 0.01)
                 (d > 0, 0.0)
                 (1.1);
  c <<= nebo_cond(false, b + 0.01)
                 (d > 0, a)
                 (1.1);
  c <<= nebo_cond(false, b + 0.01)
                 (d > 0, a + 1)
                 (1.1);

  //Cond<nebo> isNan test case:
  //c <<= nebo_cond(isnan(sqrt(a)), 1.1)
  //               (std::isnan(sqrt(-1)), 0.0)
  //               (1.0);

  print(c);
  
  return 0;
}
