#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/Nebo.h>
#include <test/TestHelper.h>
#include <test/FieldHelper.h>

#include <spatialops/structured/FieldComparisons.h>

#include <limits>
#include <cmath>
#include <boost/math/special_functions/next.hpp>

#include <sstream>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace SpatialOps;
using namespace structured;
using std::cout;
using std::endl;

//--------------------------------------------------------------------
#define COMPARE_SPRNKLE_FIELDS( F1RANGE, F1SPVAL, F1AMOUNT,                       \
                                F2RANGE, F2SPVAL, F2AMOUNT,                       \
                                EXTREMEVAL, EXTREMEAMOUNT,                        \
                                ERRORMODE, ERROR, MESSAGE, EXPECTED)              \
{                                                                                 \
  FieldT sf1(window, NULL, InternalStorage, memType1);                            \
  FieldT sf2(window, NULL, InternalStorage, memType2);                            \
                                                                                  \
  /*Intialize Field 1*/                                                           \
  fill_field_range(&sf1, 0, F1RANGE);                                             \
  double values1[F1AMOUNT+EXTREMEAMOUNT];                                         \
  size_t inc = F1AMOUNT/(EXTREMEAMOUNT + 1);                                      \
  int i = 0;                                                                      \
  for(int j = 0; j < EXTREMEAMOUNT; j++) {                                        \
    for(int k = 0; k < inc; k++) {                                                \
      values1[i] = F1SPVAL;                                                       \
      i++;                                                                        \
    }                                                                             \
    values1[i] = EXTREMEVAL;                                                      \
    i++;                                                                          \
  }                                                                               \
  for(; i < F1AMOUNT+EXTREMEAMOUNT; i++) {                                        \
    values1[i] = F1SPVAL;                                                         \
  }                                                                               \
  sprinkle_in_field(&sf1, values1, F1AMOUNT+EXTREMEAMOUNT);                       \
                                                                                  \
  /*Intialize Field 2*/                                                           \
  fill_field_range(&sf2, 0, F2RANGE);                                             \
  double values2[F2AMOUNT+EXTREMEAMOUNT];                                         \
  inc = F2AMOUNT/(EXTREMEAMOUNT + 1);                                             \
  i = 0;                                                                          \
  for(int j = 0; j < EXTREMEAMOUNT; j++) {                                        \
    for(int k = 0; k < inc; k++) {                                                \
      values2[i] = F2SPVAL;                                                       \
      i++;                                                                        \
    }                                                                             \
    values2[i] = EXTREMEVAL;                                                      \
    i++;                                                                          \
  }                                                                               \
  for(; i < F2AMOUNT+EXTREMEAMOUNT; i++) {                                        \
    values2[i] = F2SPVAL;                                                         \
  }                                                                               \
  sprinkle_in_field(&sf2, values2, F2AMOUNT+EXTREMEAMOUNT);                       \
                                                                                  \
  COMPARE_FIELDS(sf1, sf2, ERRORMODE, ERROR, MESSAGE, EXPECTED, 0);               \
}

#define COMPARE_FIELDS( F1VAL, F2VAL, ERRORMODE, ERROR, MESSAGE, EXPECTED, ABSERROR)                           \
{                                                                                                              \
  FieldT F1(window, NULL, InternalStorage, memType1);                                                          \
  FieldT F2(window, NULL, InternalStorage, memType2);                                                          \
  F1 <<= F1VAL;                                                                                                \
  F2 <<= F2VAL;                                                                                                \
  status(manual_error_compare(F1, F2, ERROR, ERRORMODE, test_field_not_equal, verboseOutput, EXPECTED, ABSERROR), MESSAGE);   \
}

#define COMPARE_MEMORY_WINDOWS( F1NPTS, F2NPTS, ERRORMODE, MESSAGE)               \
{                                                                                 \
  FieldT F1(MemoryWindow(F1NPTS), NULL, InternalStorage, memType1);               \
  FieldT F2(MemoryWindow(F2NPTS), NULL, InternalStorage, memType2);               \
                                                                                  \
  F1 <<= 1.0;                                                                     \
  F2 <<= 1.0;                                                                     \
                                                                                  \
  bool result = true;                                                             \
  try {                                                                           \
    switch(et) {                                                                  \
      case RELATIVE:                                                              \
        field_equal(F1, F2, 0.0);                                                 \
          break;                                                                  \
      case ABSOLUTE:                                                              \
        field_equal_abs(F1, F2, 0.0);                                             \
          break;                                                                  \
      case ULP:                                                                   \
        field_equal_ulp(F1, F2, 0);                                               \
          break;                                                                  \
    }                                                                             \
    result = false;                                                               \
  }                                                                               \
  catch(std::runtime_error& e) {}                                                 \
                                                                                  \
  status(result, MESSAGE);                                                        \
}


/**
 * Function fill_field_range
 *
 * Parameters:
 *     f1 = Field to fill
 *     start = start of sin period
 *     range = range filled values
 *
 * Return:
 *     void
 *
 * Use initialize_field with specified start and multiply each value by range.
 * This has the effect of creating a 'psuedo' random field of values between
 * [-range, range].
 */
template<typename FieldT>
void fill_field_range(FieldT * f1, double start, double range)
{
  FieldT * f;
  bool created_field = false;
  if(f1->memory_device_type() == EXTERNAL_CUDA_GPU) {
    created_field = true;
    f = new FieldT(f1->window_with_ghost(), NULL, InternalStorage);
  }
  else {
    f = f1;
  }

  initialize_field(*f, start, false, range);

  if(created_field) {
    *f1 = *f;
    delete f;
  }
}

/**
 * Function sprinkle_in_field
 *
 * Parameters:
 *     f1 = field to fill
 *     vals = array contianing values to insert into field
 *     size = size of vals array
 *
 * Return:
 *     void
 *
 * Evenly space the values in vals array within the field.
 */
template<typename FieldT>
void sprinkle_in_field(FieldT * f1, typename FieldT::AtomicT* vals, size_t size)
{
  size_t cells = f1->window_with_ghost().local_npts();
  if(size > cells) {
    throw(std::runtime_error("Too many values to sprinkle"));
  }

  FieldT * f;
  bool created_field = false;
  if(f1->memory_device_type() == EXTERNAL_CUDA_GPU) {
    created_field = true;
    f = new FieldT(f1->window_with_ghost(), NULL, InternalStorage);
    *f = *f1;
  }
  else {
    f = f1;
  }

  typename FieldT::iterator ifd = f->begin();
  size_t inc = cells/size;
  for(size_t i = 0; i < size; i++) {
    ifd[i*inc] = vals[i];
  }


  if(created_field) {
    *f1 = *f;
    delete f;
  }
}



enum ErrorType {RELATIVE, ABSOLUTE, ULP};
template<typename FieldT>
bool manual_error_compare(FieldT& f1,
                    FieldT& f2,
                    const double error,
                    const ErrorType et,
                    const bool test_field_not_equal,
                    const bool verboseOutput,
                    bool expected_equal,
                    const double abs_error)
{
  //copy the fields to local ram if applicable
  if(f1.memory_device_type() == EXTERNAL_CUDA_GPU) {
    f1.add_consumer(LOCAL_RAM, 0);
  }
  if(f2.memory_device_type() == EXTERNAL_CUDA_GPU) {
    f2.add_consumer(LOCAL_RAM, 0);
  }

  //iterate through fields.
  typename FieldT::const_iterator if1 = const_cast<const FieldT&>(f1).begin();
  typename FieldT::const_iterator if1e = const_cast<const FieldT&>(f1).end();
  typename FieldT::const_iterator if2 = const_cast<const FieldT&>(f2).begin();
  const std::numeric_limits<double> nl;

  //manually determine equality
  bool man_equal = true;
  for(; if1 != if1e; ++if1, ++if2) {
    double diff = *if1 - *if2;
    switch(et) {
      case RELATIVE:
        double denom;
        if(abs_error) {
          denom = std::abs(*if1) + abs_error;
        }
        else {
          //Defualt absolute error in SpatialField
          denom = std::abs(*if1) + nebo_norm(f1) * error * .000001;
        }

        if( std::abs(diff)/denom > error ) {
          man_equal = false;
        }
        break;
      case ABSOLUTE:
        if( std::abs(diff) > error ) {
          man_equal = false;
        }
        break;
      case ULP:
        double limit = *if1;
        const double inf = diff > 0 ? -nl.infinity() : nl.infinity();
        for(int i = 0; i < error; i++) {
          try {
            //get floats in direction of *if2
            limit = nextafter(limit, inf);
          } catch(...) {
            limit = inf;
            break;
          }
        }
        if( std::abs(diff) > std::abs(*if1 - limit) ) {
          man_equal = false;
        }
        break;
    }
    if (!man_equal) break;
  }

  std::ostringstream msg;
  msg << "Manual Compare Result: " << (man_equal ? "Equal" : "Not Equal");
  TestHelper tmp(verboseOutput);
  tmp(man_equal == expected_equal, msg.str());

  //switch expected_equal and manual equal result based on test_field_not_equal compare
  if(test_field_not_equal) {man_equal = !man_equal; expected_equal = !expected_equal;}

  //compare manual to function
  switch(et) {
    case RELATIVE:
      if(test_field_not_equal) {
        if(abs_error)
          return (man_equal == field_not_equal(f1, f2, error, abs_error)) && (man_equal == expected_equal);
        else
          return (man_equal == field_not_equal(f1, f2, error)) && (man_equal == expected_equal);
      }
      else {
        if(abs_error)
          return (man_equal == field_equal(f1, f2, error, abs_error)) && (man_equal == expected_equal);
        else
          return (man_equal == field_equal(f1, f2, error)) && (man_equal == expected_equal);
      }
    case ABSOLUTE:
      if(test_field_not_equal) {
        return (man_equal == field_not_equal_abs(f1, f2, error)) && (man_equal == expected_equal);
      }
      else {
        return (man_equal == field_equal_abs(f1, f2, error)) && (man_equal == expected_equal);
      }
    case ULP:
      if(test_field_not_equal) {
        return (man_equal == field_not_equal_ulp(f1, f2, error)) && (man_equal == expected_equal);
      }
      else {
        return (man_equal == field_equal_ulp(f1, f2, error)) && (man_equal == expected_equal);
      }
  }
  return false;
}

template<typename FieldT>
bool test_field_equal( const IntVec npts,
                       const MemoryType memType1,
                       const MemoryType memType2,
                       const ErrorType et,
                       const bool test_field_not_equal,
                       const bool verboseOutput)
{
  TestHelper status(verboseOutput);
  const std::numeric_limits<double> nl;
  const MemoryWindow window(npts);
  const int total = npts[0] * npts[1] * npts[2];

  FieldT* f1;
  FieldT* f2;

  //local fields
  FieldT lf1(window, NULL, InternalStorage);
  FieldT lf2(window, NULL, InternalStorage);
  FieldT lf3(window, NULL, InternalStorage);

  initialize_field(lf1, 0);
  initialize_field(lf2, 0);
  initialize_field(lf3, total);
#ifdef __CUDACC__
  //gpu fields
  FieldT gf1(window, NULL, InternalStorage, EXTERNAL_CUDA_GPU, 0);
  FieldT gf2(window, NULL, InternalStorage, EXTERNAL_CUDA_GPU, 0);
  FieldT gf3(window, NULL, InternalStorage, EXTERNAL_CUDA_GPU, 0);
  //move local initialized fields to gpu if necessary
  if(memType1 == EXTERNAL_CUDA_GPU) {
    lf1.add_consumer(EXTERNAL_CUDA_GPU, 0);
    gf1 <<= lf1;
    f1 = &gf1;
  }
  else {
    f1 = &lf1;
  }
  if(memType2 == EXTERNAL_CUDA_GPU) {
    lf2.add_consumer(EXTERNAL_CUDA_GPU, 0);
    gf2 <<= lf2;
    f2 = &gf2;
  }
  else {
    f2 = &lf2;
  }
#else
  f1 = &lf1;
  f2 = &lf2;
#endif

  //two duplicate fields exactly equal
  status(manual_error_compare(*f1, *f2, 0.0, et, test_field_not_equal, verboseOutput, true, 0), "Duplicate Fields Equal");

  //change second field and compare not equal
#ifdef __CUDACC__
  if(memType2 == EXTERNAL_CUDA_GPU) {
    lf3.add_consumer(EXTERNAL_CUDA_GPU, 0);
    gf3 <<= lf3;
    f2 = &gf3;
  }
  else {
    f2 = &lf3;
  }
#else
  f2 = &lf3;
#endif
  status(manual_error_compare(*f1, *f2, 0.0, et, test_field_not_equal, verboseOutput, false, 0), "Non-Duplicate Fields Not Equal");


  switch(et) {
    case RELATIVE:
      //relative error test
      COMPARE_FIELDS(3.0, 7.49, RELATIVE, 1.5, "Off By 150% (Equal)", true, 0);
      COMPARE_FIELDS(3.0, 7.51, RELATIVE, 1.5, "Off By 150% (Not Equal)", false, 0);
      COMPARE_FIELDS(3.0, 3.98, RELATIVE, .33, "Off By 33% (Equal)", true, 0);
      COMPARE_FIELDS(3.0, 4.10, RELATIVE, .33, "Off By 33% (Not Equal)", false, 0);
      COMPARE_FIELDS(3.0, 2.97, RELATIVE, .01, "Off By 1% (Equal)", true, 0);
      COMPARE_FIELDS(3.0, 2.96, RELATIVE, .01, "Off By 1% (Not Equal)", false, 0);

      //near zero value tests
      COMPARE_SPRNKLE_FIELDS(1, 0, 8, 1, 1e-9, 8, 0, 0, RELATIVE, .01, "Near Zero Field Range [-1, 1] Off By 1% (Equal)", true);
      COMPARE_SPRNKLE_FIELDS(1, 0, 8, 1, 1e-8, 8, 0, 0, RELATIVE, .01, "Near Zero Field Range [-1, 1] Off By 1% (Not Equal)", false);
      COMPARE_SPRNKLE_FIELDS(10, 0, 8, 10, 1e-8, 8, 0, 0, RELATIVE, .01, "Near Zero Field Range [-10, 10] Off By 1% (Equal)", true);
      COMPARE_SPRNKLE_FIELDS(10, 0, 8, 10, 1e-7, 8, 0, 0, RELATIVE, .01, "Near Zero Field Range [-10, 10] Off By 1% (Not Equal)", false);
      COMPARE_SPRNKLE_FIELDS(100000, 0, 8, 100000, 1e-4, 8, 0, 0, RELATIVE, .01, "Near Zero Field Range [-100000, 100000] Off By 1% (Equal)", true);
      COMPARE_SPRNKLE_FIELDS(100000, 0, 8, 100000, 1e-3, 8, 0, 0, RELATIVE, .01, "Near Zero Field Range [-100000, 100000] Off By 1% (Not Equal)", false);
      COMPARE_SPRNKLE_FIELDS(1, 0, 8, 1, 1e-7, 8, 1000, 1, RELATIVE, .01, "Outlier Near Zero Field Range [-1, 1] Off By 1% (Equal)", true);
      COMPARE_SPRNKLE_FIELDS(1, 0, 8, 1, 1e-6, 8, 1000, 1, RELATIVE, .01, "Outlier Near Zero Field Range [-1, 1] Off By 1% (Not Equal)", false);
      COMPARE_SPRNKLE_FIELDS(1, 0, 8, 1, 1e-7, 8, 1000, 5, RELATIVE, .01, "Five Outlier Near Zero Field Range [-1, 1] Off By 1% (Equal)", true);
      COMPARE_SPRNKLE_FIELDS(1, 0, 8, 1, 1e-6, 8, 1000, 5, RELATIVE, .01, "Five Outlier Near Zero Field Range [-1, 1] Off By 1% (Not Equal)", false);

      //Custom Absolute Error
      COMPARE_FIELDS(1.0, 3.0, RELATIVE, 1.0, "Absolute Tolerance 1: 1-3 Off By 100% (Equal)", true, 1.0);
      COMPARE_FIELDS(1.0, 4.0, RELATIVE, 1.0, "Absolute Tolerance 1: 1-4 Off By 100% (Not Equal)", false, 1.0);
      COMPARE_FIELDS(1.0, 5.0, RELATIVE, 2.0, "Absolute Tolerance 1: 1-5 Off By 200% (Equal)", true, 1.0);
      COMPARE_FIELDS(1.0, 5.0, RELATIVE, 1.0, "Absolute Tolerance 1: 1-5 Off By 100% (Not Equal)", false, 1.0);
      break;
    case ABSOLUTE:
      //absolute error test
      COMPARE_FIELDS(1.0, 1.0 + nl.epsilon(), ABSOLUTE, nl.epsilon(), "Off By epsilon (Equal)", true, 0);
      COMPARE_FIELDS(1.0, 1.0 + 2*nl.epsilon(), ABSOLUTE, nl.epsilon(), "Off By epsilon (Not Equal)", false, 0);
      COMPARE_FIELDS(0.033, 0.024, ABSOLUTE, std::pow(10,-2), "Off By 10^-2 (Equal)", true, 0);
      COMPARE_FIELDS(0.033, 0.022, ABSOLUTE, std::pow(10,-2), "Off By 10^-2 (Not Equal)", false, 0);
      COMPARE_FIELDS(4679000.0, 4680000.0, ABSOLUTE, std::pow(10,3), "Off By 10^3 (Equal)", true, 0);
      COMPARE_FIELDS(4679000.0, 4681000.0, ABSOLUTE, std::pow(10,3), "Off By 10^3 (Not Equal)", false, 0);
      COMPARE_FIELDS(4679000.0, 6890330.0, ABSOLUTE, 11569300.0, "Large Number Check, Exact Error (Equal)", true, 0);
      break;
    case ULP:
      using boost::math::float_prior;
      using boost::math::float_next;
      using boost::math::float_advance;

      //near zero value tests
      COMPARE_FIELDS(0.0, float_next(0.0), ULP, 1, "Near Zero Off By 1 ulp (Equal)", true, 0);
      COMPARE_FIELDS(0.0, float_advance(0.0, 2), ULP, 1, "Near Zero Off By 1 ulp (Not Equal)", false, 0);
      COMPARE_FIELDS(0.0, -float_advance(0.0, 2), ULP, 2, "Near Zero Off By 2 ulps (Equal)", true, 0);
      COMPARE_FIELDS(0.0, -float_advance(0.0, 3), ULP, 2, "Near Zero Off By 2 ulps (Not Equal)", false, 0);
      COMPARE_FIELDS(-float_advance(0.0, 1), 0.0, ULP, 1, "Near Zero Reversed Off By 1 ulp (Equal)", true, 0);
      COMPARE_FIELDS(float_advance(0.0, 2), 0.0, ULP, 1, "Near Zero Reversed Off By 1 ulp (Not Equal)", false, 0);

      //machine epsilon
      COMPARE_FIELDS(1.0, 1.0+nl.epsilon(), ULP, 1, "Machine epsilon at 1.0 is 1 ulp (Equal)", true, 0);
      COMPARE_FIELDS(1.0, 1.0+nl.epsilon(), ULP, 0, "Machine epsilon at 1.0 is not 0 ulps (Not Equal)", false, 0);

      //ulps error test
      COMPARE_FIELDS(3.0, float_prior(3.0), ULP, 1, "Off By 1 ulp (Equal)", true, 0);
      COMPARE_FIELDS(3.0, float_advance(3.0, -2), ULP, 1, "Off By 1 ulp (Not Equal)", false, 0);
      COMPARE_FIELDS(3.0, float_advance(3.0, 3), ULP, 3, "Off By 3 ulps (Equal)", true, 0);
      COMPARE_FIELDS(3.0, float_advance(3.0, 4), ULP, 3, "Off By 3 ulps (Not Equal)", false, 0);
      COMPARE_FIELDS(3.0, float_advance(3.0, 20), ULP, 20, "Off By 20 ulps (Equal)", true, 0);
      COMPARE_FIELDS(3.0, float_advance(3.0, -21), ULP, 20, "Off By 20 ulps (Not Equal)", false, 0);

      //limits test
      COMPARE_FIELDS(nl.max(), float_advance(nl.max(), -1), ULP, 1, "Max Double Off By -1 ulp (Equal)", true, 0);
      COMPARE_FIELDS(nl.min(), float_advance(nl.min(), 1), ULP, 1, "Min > 0 Off By 1 ulp (Equal)", true, 0);
      COMPARE_FIELDS(nl.min(), float_advance(nl.min(), -1), ULP, 1, "Min > 0 Off By -1 ulp (Equal)", true, 0);
      COMPARE_FIELDS(-nl.max(), float_advance(-nl.max(), 1), ULP, 1, "Min Off By 1 ulps (Equal)", true, 0);
      try {
        //manual compare useless
        COMPARE_FIELDS(nl.max(), nl.infinity(), ULP, 20, "fail", true, 0);
        status(false, "Max Double = Infinity By 20 ulps (Throws Exception)");
      } catch(std::domain_error) {status(true,  "Max Double = Infinity By 20 ulps (Throws Exception)");}
      try {
        //manual compare useless
        COMPARE_FIELDS(-nl.max(), -nl.infinity(), ULP, 20, "fail", true, 0);
        status(false, "Min Double = Infinity By 20 ulps (Throws Exception)");
      } catch(std::domain_error) {status(true, "Min Double = Infinity By 20 ulps (Throws Exception)");}
      try {
        //manual compare useless
        COMPARE_FIELDS(nan(""), nan(""), ULP, 20, "fail", true, 0);
        status(false, "NAN = NAN By 20 ulps (Throws Exception)");
      } catch(std::domain_error) {status(true, "NAN = NAN By 20 ulps (Throws Exception)");}
      break;
  }

  //unequal memory window tests
  COMPARE_MEMORY_WINDOWS(npts, IntVec(npts[0]+5, npts[1], npts[2]), et,
      "Cannot Compare Windows Of Different Dimensions X");
  COMPARE_MEMORY_WINDOWS(npts, IntVec(npts[0], npts[1]+2, npts[2]), et,
      "Cannot Compare Windows Of Different Dimensions Y");
  COMPARE_MEMORY_WINDOWS(npts, IntVec(npts[0], npts[1], npts[2]+7), et,
      "Cannot Compare Windows Of Different Dimensions Z");
  COMPARE_MEMORY_WINDOWS(npts, IntVec(npts[1], npts[2], npts[0]), et,
      "Cannot Compare Windows Of Different Dimensions, Same Memory Size");

  return status.ok();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////                                          /////////////////////////////
/////////////////////////////          SCALAR IMPLEMENTATION           /////////////////////////////
/////////////////////////////                                          /////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename FieldT>
bool manual_error_compare(double d,
                    FieldT& f1,
                    const double error,
                    const ErrorType et,
                    const bool test_field_not_equal,
                    const bool verboseOutput,
                    bool expected_equal,
                    const double abs_error)
{
  //copy the fields to local ram if applicable
  if(f1.memory_device_type() == EXTERNAL_CUDA_GPU) {
    f1.add_consumer(LOCAL_RAM, 0);
  }

  //iterate through fields.
  typename FieldT::const_iterator if1 = const_cast<const FieldT&>(f1).begin();
  typename FieldT::const_iterator if1e = const_cast<const FieldT&>(f1).end();
  const std::numeric_limits<double> nl;

  //manually determine equality
  bool man_equal = true;
  double denom;
  if(abs_error) {
    denom = std::abs(d) + abs_error;
  }
  else {
    //Defualt absolute error in SpatialField
    denom = std::abs(d) + nebo_norm(f1) * error * .000001;
  }
  for(; if1 != if1e; ++if1) {
    double diff = d - *if1;
    switch(et) {
      case RELATIVE:

        if( std::abs(diff)/denom > error ) {
          man_equal = false;
        }
        break;
      case ABSOLUTE:
        if( std::abs(diff) > error ) {
          man_equal = false;
        }
        break;
      case ULP:
        double limit = d;
        const double inf = diff > 0 ? -nl.infinity() : nl.infinity();
        for(int i = 0; i < error; i++) {
          try {
            //get floats in direction of d
            limit = nextafter(limit, inf);
          } catch(...) {
            limit = inf;
            break;
          }
        }
        if( std::abs(diff) > std::abs(d - limit) ) {
          man_equal = false;
        }
        break;
    }
    if (!man_equal) break;
  }

  std::ostringstream msg;
  msg << "Manual Compare Result: " << (man_equal ? "Equal" : "Not Equal");
  TestHelper tmp(verboseOutput);
  tmp(man_equal == expected_equal, msg.str());

  //switch expected_equal and manual equal result based on test_field_not_equal compare
  if(test_field_not_equal) {man_equal = !man_equal; expected_equal = !expected_equal;}

  //compare manual to function
  switch(et) {
    case RELATIVE:
      if(test_field_not_equal) {
        if(abs_error)
          return (man_equal == field_not_equal(d, f1, error, abs_error)) && (man_equal == expected_equal);
        else
          return (man_equal == field_not_equal(d, f1, error)) && (man_equal == expected_equal);
      }
      else {
        if(abs_error)
          return (man_equal == field_equal(d, f1, error, abs_error)) && (man_equal == expected_equal);
        else
          return (man_equal == field_equal(d, f1, error)) && (man_equal == expected_equal);
      }
    case ABSOLUTE:
      if(test_field_not_equal) {
        return (man_equal == field_not_equal_abs(d, f1, error)) && (man_equal == expected_equal);
      }
      else {
        return (man_equal == field_equal_abs(d, f1, error)) && (man_equal == expected_equal);
      }
    case ULP:
      if(test_field_not_equal) {
        return (man_equal == field_not_equal_ulp(d, f1, error)) && (man_equal == expected_equal);
      }
      else {
        return (man_equal == field_equal_ulp(d, f1, error)) && (man_equal == expected_equal);
      }
  }
  return false;
}

#define COMPARE_FIELDS_SCALAR( VAL, F1VAL, ERRORMODE, ERROR, MESSAGE, EXPECTED, ABSERROR)                      \
{                                                                                                              \
  FieldT F1(window, NULL, InternalStorage, memType1);                                                          \
  F1 <<= F1VAL;                                                                                                \
  status(manual_error_compare(VAL, F1, ERROR, ERRORMODE, test_field_not_equal, verboseOutput, EXPECTED, ABSERROR), MESSAGE);   \
}

template<typename FieldT>
bool test_field_equal_scalar( const IntVec npts,
                              const MemoryType memType1,
                              const ErrorType et,
                              const bool test_field_not_equal,
                              const bool verboseOutput)
{
  TestHelper status(verboseOutput);
  const std::numeric_limits<double> nl;
  const MemoryWindow window(npts);
  const int total = npts[0] * npts[1] * npts[2];

  FieldT* f1;

  //local field
  FieldT lf1(window, NULL, InternalStorage);

  lf1 <<= 42.0;
#ifdef __CUDACC__
  //gpu fields
  FieldT gf1(window, NULL, InternalStorage, EXTERNAL_CUDA_GPU, 0);
  //move local initialized field to gpu if necessary
  if(memType1 == EXTERNAL_CUDA_GPU) {
    lf1.add_consumer(EXTERNAL_CUDA_GPU, 0);
    gf1 <<= lf1;
    f1 = &gf1;
  }
  else {
    f1 = &lf1;
  }
#else
  f1 = &lf1;
#endif

  //field exactly equal
  status(manual_error_compare(42.0, *f1, 0.0, et, test_field_not_equal, verboseOutput, true, 0), "Duplicate Fields Equal");

  //field not equal
  status(manual_error_compare(21.0, *f1, 0.0, et, test_field_not_equal, verboseOutput, false, 0), "Non-Duplicate Fields Not Equal");


  switch(et) {
    case RELATIVE:
      //relative error test
      COMPARE_FIELDS_SCALAR(3.0, 7.49, RELATIVE, 1.5, "Off By 150% (Equal)", true, 0);
      COMPARE_FIELDS_SCALAR(3.0, 7.51, RELATIVE, 1.5, "Off By 150% (Not Equal)", false, 0);
      COMPARE_FIELDS_SCALAR(3.0, 3.98, RELATIVE, .33, "Off By 33% (Equal)", true, 0);
      COMPARE_FIELDS_SCALAR(3.0, 4.10, RELATIVE, .33, "Off By 33% (Not Equal)", false, 0);
      COMPARE_FIELDS_SCALAR(3.0, 2.97, RELATIVE, .01, "Off By 1% (Equal)", true, 0);
      COMPARE_FIELDS_SCALAR(3.0, 2.96, RELATIVE, .01, "Off By 1% (Not Equal)", false, 0);

      //Custom Absolute Error
      COMPARE_FIELDS_SCALAR(1.0, 3.0, RELATIVE, 1.0, "Absolute Tolerance 1: 1-3 Off By 100% (Equal)", true, 1.0);
      COMPARE_FIELDS_SCALAR(1.0, 4.0, RELATIVE, 1.0, "Absolute Tolerance 1: 1-4 Off By 100% (Not Equal)", false, 1.0);
      COMPARE_FIELDS_SCALAR(1.0, 5.0, RELATIVE, 2.0, "Absolute Tolerance 1: 1-5 Off By 200% (Equal)", true, 1.0);
      COMPARE_FIELDS_SCALAR(1.0, 5.0, RELATIVE, 1.0, "Absolute Tolerance 1: 1-5 Off By 100% (Not Equal)", false, 1.0);
      break;
    case ABSOLUTE:
      //absolute error test
      COMPARE_FIELDS_SCALAR(1.0, 1.0 + nl.epsilon(), ABSOLUTE, nl.epsilon(), "Off By epsilon (Equal)", true, 0);
      COMPARE_FIELDS_SCALAR(1.0, 1.0 + 2*nl.epsilon(), ABSOLUTE, nl.epsilon(), "Off By epsilon (Not Equal)", false, 0);
      COMPARE_FIELDS_SCALAR(0.033, 0.024, ABSOLUTE, std::pow(10,-2), "Off By 10^-2 (Equal)", true, 0);
      COMPARE_FIELDS_SCALAR(0.033, 0.022, ABSOLUTE, std::pow(10,-2), "Off By 10^-2 (Not Equal)", false, 0);
      COMPARE_FIELDS_SCALAR(4679000.0, 4680000.0, ABSOLUTE, std::pow(10,3), "Off By 10^3 (Equal)", true, 0);
      COMPARE_FIELDS_SCALAR(4679000.0, 4681000.0, ABSOLUTE, std::pow(10,3), "Off By 10^3 (Not Equal)", false, 0);
      COMPARE_FIELDS_SCALAR(4679000.0, 6890330.0, ABSOLUTE, 11569300.0, "Large Number Check, Exact Error (Equal)", true, 0);
      break;
    case ULP:
      using boost::math::float_prior;
      using boost::math::float_next;
      using boost::math::float_advance;

      //near zero value tests
      COMPARE_FIELDS_SCALAR(0.0, float_next(0.0), ULP, 1, "Near Zero Off By 1 ulp (Equal)", true, 0);
      COMPARE_FIELDS_SCALAR(0.0, float_advance(0.0, 2), ULP, 1, "Near Zero Off By 1 ulp (Not Equal)", false, 0);
      COMPARE_FIELDS_SCALAR(0.0, -float_advance(0.0, 2), ULP, 2, "Near Zero Off By 2 ulps (Equal)", true, 0);
      COMPARE_FIELDS_SCALAR(0.0, -float_advance(0.0, 3), ULP, 2, "Near Zero Off By 2 ulps (Not Equal)", false, 0);
      COMPARE_FIELDS_SCALAR(-float_advance(0.0, 1), 0.0, ULP, 1, "Near Zero Reversed Off By 1 ulp (Equal)", true, 0);
      COMPARE_FIELDS_SCALAR(float_advance(0.0, 2), 0.0, ULP, 1, "Near Zero Reversed Off By 1 ulp (Not Equal)", false, 0);

      //machine epsilon
      COMPARE_FIELDS_SCALAR(1.0, 1.0+nl.epsilon(), ULP, 1, "Machine epsilon at 1.0 is 1 ulp (Equal)", true, 0);
      COMPARE_FIELDS_SCALAR(1.0, 1.0+nl.epsilon(), ULP, 0, "Machine epsilon at 1.0 is not 0 ulps (Not Equal)", false, 0);

      //ulps error test
      COMPARE_FIELDS_SCALAR(3.0, float_prior(3.0), ULP, 1, "Off By 1 ulp (Equal)", true, 0);
      COMPARE_FIELDS_SCALAR(3.0, float_advance(3.0, -2), ULP, 1, "Off By 1 ulp (Not Equal)", false, 0);
      COMPARE_FIELDS_SCALAR(3.0, float_advance(3.0, 3), ULP, 3, "Off By 3 ulps (Equal)", true, 0);
      COMPARE_FIELDS_SCALAR(3.0, float_advance(3.0, 4), ULP, 3, "Off By 3 ulps (Not Equal)", false, 0);
      COMPARE_FIELDS_SCALAR(3.0, float_advance(3.0, 20), ULP, 20, "Off By 20 ulps (Equal)", true, 0);
      COMPARE_FIELDS_SCALAR(3.0, float_advance(3.0, -21), ULP, 20, "Off By 20 ulps (Not Equal)", false, 0);

      //limits test
      COMPARE_FIELDS_SCALAR(nl.max(), float_advance(nl.max(), -1), ULP, 1, "Max Double Off By -1 ulp (Equal)", true, 0);
      COMPARE_FIELDS_SCALAR(nl.min(), float_advance(nl.min(), 1), ULP, 1, "Min > 0 Off By 1 ulp (Equal)", true, 0);
      COMPARE_FIELDS_SCALAR(nl.min(), float_advance(nl.min(), -1), ULP, 1, "Min > 0 Off By -1 ulp (Equal)", true, 0);
      COMPARE_FIELDS_SCALAR(-nl.max(), float_advance(-nl.max(), 1), ULP, 1, "Min Off By 1 ulps (Equal)", true, 0);
      try {
        //manual compare useless
        COMPARE_FIELDS_SCALAR(nl.max(), nl.infinity(), ULP, 20, "fail", true, 0);
        status(false, "Max Double = Infinity By 20 ulps (Throws Exception)");
      } catch(std::domain_error) {status(true,  "Max Double = Infinity By 20 ulps (Throws Exception)");}
      try {
        //manual compare useless
        COMPARE_FIELDS_SCALAR(-nl.max(), -nl.infinity(), ULP, 20, "fail", true, 0);
        status(false, "Min Double = Infinity By 20 ulps (Throws Exception)");
      } catch(std::domain_error) {status(true, "Min Double = Infinity By 20 ulps (Throws Exception)");}
      try {
        //manual compare useless
        COMPARE_FIELDS_SCALAR(nan(""), nan(""), ULP, 20, "fail", true, 0);
        status(false, "NAN = NAN By 20 ulps (Throws Exception)");
      } catch(std::domain_error) {status(true, "NAN = NAN By 20 ulps (Throws Exception)");}
      break;
  }

  return status.ok();
}

int main(int argc, const char *argv[])
{

  int nx, ny, nz;
  po::options_description desc("Supported Options");
  desc.add_options()
    ( "help", "print help message\n" )
    ( "nx",   po::value<int>(&nx)->default_value(10), "number of points in x-dir for base mesh" )
    ( "ny",   po::value<int>(&ny)->default_value(11), "number of points in y-dir for base mesh" )
    ( "nz",   po::value<int>(&nz)->default_value(12), "number of points in z-dir for base mesh" );

  po::variables_map args;
  po::store( po::parse_command_line(argc,argv,desc), args );
  po::notify(args);
  if( args.count("help") ){
    cout << desc << endl
      << "Examples:" << endl
      << " test_field_compare --nx 5 --ny 10 --nz 3" << endl
      << " test_field_compare --nx 50" << endl
      << endl;
    return -1;
  }

  TestHelper overall(true);

  bool fieldEqualVerbose = false;
  IntVec winSize(nx, ny, nz);

  //test field comparison functions
  {
    //test field_equal
    overall(test_field_equal<SVolField>(winSize, LOCAL_RAM, LOCAL_RAM, RELATIVE, false, fieldEqualVerbose), "LOCAL_RAM x LOCAL_RAM Relative Equal Test");
#ifdef __CUDACC__
    overall(test_field_equal<SVolField>(winSize, LOCAL_RAM, EXTERNAL_CUDA_GPU, RELATIVE, false, fieldEqualVerbose), "LOCAL_RAM x EXTERNAL_CUDA_GPU Relative Equal Test");
    overall(test_field_equal<SVolField>(winSize, EXTERNAL_CUDA_GPU, LOCAL_RAM, RELATIVE, false, fieldEqualVerbose), "EXTERNAL_CUDA_GPU x LOCAL_RAM Relative Equal Test");
    overall(test_field_equal<SVolField>(winSize, EXTERNAL_CUDA_GPU, EXTERNAL_CUDA_GPU, RELATIVE, false, fieldEqualVerbose), "EXTERNAL_CUDA_GPU x EXTERNAL_CUDA_GPU Relative Equal Test");
#endif

    //test field_equal_abs
    overall(test_field_equal<SVolField>(winSize, LOCAL_RAM, LOCAL_RAM, ABSOLUTE, false, fieldEqualVerbose), "LOCAL_RAM x LOCAL_RAM Absolute Error Equal Test");
#ifdef __CUDACC__
    overall(test_field_equal<SVolField>(winSize, LOCAL_RAM, EXTERNAL_CUDA_GPU, ABSOLUTE, false, fieldEqualVerbose), "LOCAL_RAM x EXTERNAL_CUDA_GPU Absolute Error Equal Test");
    overall(test_field_equal<SVolField>(winSize, EXTERNAL_CUDA_GPU, LOCAL_RAM, ABSOLUTE, false, fieldEqualVerbose), "EXTERNAL_CUDA_GPU x LOCAL_RAM Absolute Error Equal Test");
    overall(test_field_equal<SVolField>(winSize, EXTERNAL_CUDA_GPU, EXTERNAL_CUDA_GPU, ABSOLUTE, false, fieldEqualVerbose), "EXTERNAL_CUDA_GPU x EXTERNAL_CUDA_GPU Absolute Error Equal Test");
#endif

    //test field_equal_ulp
    overall(test_field_equal<SVolField>(winSize, LOCAL_RAM, LOCAL_RAM, ULP, false, fieldEqualVerbose), "LOCAL_RAM x LOCAL_RAM Ulps Equal Test");
#ifdef __CUDACC__
    overall(test_field_equal<SVolField>(winSize, LOCAL_RAM, EXTERNAL_CUDA_GPU, ULP, false, fieldEqualVerbose), "LOCAL_RAM x EXTERNAL_CUDA_GPU Ulps Equal Test");
    overall(test_field_equal<SVolField>(winSize, EXTERNAL_CUDA_GPU, LOCAL_RAM, ULP, false, fieldEqualVerbose), "EXTERNAL_CUDA_GPU x LOCAL_RAM Ulps Equal Test");
    overall(test_field_equal<SVolField>(winSize, EXTERNAL_CUDA_GPU, EXTERNAL_CUDA_GPU, ULP, false, fieldEqualVerbose), "EXTERNAL_CUDA_GPU x EXTERNAL_CUDA_GPU Ulps Equal Test");
#endif

    //test field_not_equal
    overall(test_field_equal<SVolField>(winSize, LOCAL_RAM, LOCAL_RAM, RELATIVE, true, fieldEqualVerbose), "LOCAL_RAM x LOCAL_RAM Relative Not Equal Test");
#ifdef __CUDACC__
    overall(test_field_equal<SVolField>(winSize, LOCAL_RAM, EXTERNAL_CUDA_GPU, RELATIVE, true, fieldEqualVerbose), "LOCAL_RAM x EXTERNAL_CUDA_GPU Relative Not Equal Test");
    overall(test_field_equal<SVolField>(winSize, EXTERNAL_CUDA_GPU, LOCAL_RAM, RELATIVE, true, fieldEqualVerbose), "EXTERNAL_CUDA_GPU x LOCAL_RAM Relative Not Equal Test");
    overall(test_field_equal<SVolField>(winSize, EXTERNAL_CUDA_GPU, EXTERNAL_CUDA_GPU, RELATIVE, true, fieldEqualVerbose), "EXTERNAL_CUDA_GPU x EXTERNAL_CUDA_GPU Relative Not Equal Test");
#endif

    //test field_not_equal_abs
    overall(test_field_equal<SVolField>(winSize, LOCAL_RAM, LOCAL_RAM, ABSOLUTE, true, fieldEqualVerbose), "LOCAL_RAM x LOCAL_RAM Absolute Error Not Equal Test");
#ifdef __CUDACC__
    overall(test_field_equal<SVolField>(winSize, LOCAL_RAM, EXTERNAL_CUDA_GPU, ABSOLUTE, true, fieldEqualVerbose), "LOCAL_RAM x EXTERNAL_CUDA_GPU Absolute Error Not Equal Test");
    overall(test_field_equal<SVolField>(winSize, EXTERNAL_CUDA_GPU, LOCAL_RAM, ABSOLUTE, true, fieldEqualVerbose), "EXTERNAL_CUDA_GPU x LOCAL_RAM Absolute Error Not Equal Test");
    overall(test_field_equal<SVolField>(winSize, EXTERNAL_CUDA_GPU, EXTERNAL_CUDA_GPU, ABSOLUTE, true, fieldEqualVerbose), "EXTERNAL_CUDA_GPU x EXTERNAL_CUDA_GPU Absolute Error Not Equal Test");
#endif

    //test field_not_equal_ulp
    overall(test_field_equal<SVolField>(winSize, LOCAL_RAM, LOCAL_RAM, ULP, true, fieldEqualVerbose), "LOCAL_RAM x LOCAL_RAM Ulps Not Equal Test");
#ifdef __CUDACC__
    overall(test_field_equal<SVolField>(winSize, LOCAL_RAM, EXTERNAL_CUDA_GPU, ULP, true, fieldEqualVerbose), "LOCAL_RAM x EXTERNAL_CUDA_GPU Ulps Not Equal Test");
    overall(test_field_equal<SVolField>(winSize, EXTERNAL_CUDA_GPU, LOCAL_RAM, ULP, true, fieldEqualVerbose), "EXTERNAL_CUDA_GPU x LOCAL_RAM Ulps Not Equal Test");
    overall(test_field_equal<SVolField>(winSize, EXTERNAL_CUDA_GPU, EXTERNAL_CUDA_GPU, ULP, true, fieldEqualVerbose), "EXTERNAL_CUDA_GPU x EXTERNAL_CUDA_GPU Ulps Not Equal Test");
#endif
  }


  //test field comparison functions with a scalar value
  {
    //test field_equal with scalar
    overall(test_field_equal_scalar<SVolField>(winSize, LOCAL_RAM, RELATIVE, false, fieldEqualVerbose), "SCALAR x LOCAL_RAM Relative Equal Test");
#ifdef __CUDACC__
    overall(test_field_equal_scalar<SVolField>(winSize, EXTERNAL_CUDA_GPU, RELATIVE, false, fieldEqualVerbose), "SCALAR x EXTERNAL_CUDA_GPU Relative Equal Test");
#endif

    //test field_equal_abs with scalar
    overall(test_field_equal_scalar<SVolField>(winSize, LOCAL_RAM, ABSOLUTE, false, fieldEqualVerbose), "SCALAR x LOCAL_RAM Absolute Error Equal Test");
#ifdef __CUDACC__
    overall(test_field_equal_scalar<SVolField>(winSize, EXTERNAL_CUDA_GPU, ABSOLUTE, false, fieldEqualVerbose), "SCALAR x EXTERNAL_CUDA_GPU Absolute Error Equal Test");
#endif

    //test field_equal_ulp with scalar
    overall(test_field_equal_scalar<SVolField>(winSize, LOCAL_RAM, ULP, false, fieldEqualVerbose), "SCALAR x LOCAL_RAM Ulps Equal Test");
#ifdef __CUDACC__
    overall(test_field_equal_scalar<SVolField>(winSize, EXTERNAL_CUDA_GPU, ULP, false, fieldEqualVerbose), "SCALAR x EXTERNAL_CUDA_GPU Ulps Equal Test");
#endif

    //test field_not_equal with scalar
    overall(test_field_equal_scalar<SVolField>(winSize, LOCAL_RAM, RELATIVE, true, fieldEqualVerbose), "SCALAR x LOCAL_RAM Relative Not Equal Test");
#ifdef __CUDACC__
    overall(test_field_equal_scalar<SVolField>(winSize, EXTERNAL_CUDA_GPU, RELATIVE, true, fieldEqualVerbose), "SCALAR x EXTERNAL_CUDA_GPU Relative Not Equal Test");
#endif

    //test field_not_equal_abs with scalar
    overall(test_field_equal_scalar<SVolField>(winSize, LOCAL_RAM, ABSOLUTE, true, fieldEqualVerbose), "SCALAR x LOCAL_RAM Absolute Error Not Equal Test");
#ifdef __CUDACC__
    overall(test_field_equal_scalar<SVolField>(winSize, EXTERNAL_CUDA_GPU, ABSOLUTE, true, fieldEqualVerbose), "SCALAR x EXTERNAL_CUDA_GPU Absolute Error Not Equal Test");
#endif

    //test field_not_equal_ulp with scalar
    overall(test_field_equal_scalar<SVolField>(winSize, LOCAL_RAM, ULP, true, fieldEqualVerbose), "SCALAR x LOCAL_RAM Ulps Not Equal Test");
#ifdef __CUDACC__
    overall(test_field_equal_scalar<SVolField>(winSize, EXTERNAL_CUDA_GPU, ULP, true, fieldEqualVerbose), "SCALAR x EXTERNAL_CUDA_GPU Ulps Not Equal Test");
#endif
  }

  if( overall.isfailed() ){
    std::cout << "FAIL!" << std::endl;
    return -1;
  }
  std::cout << "PASS" << std::endl;
  return 0;
}
