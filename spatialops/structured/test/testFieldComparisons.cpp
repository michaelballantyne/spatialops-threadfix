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
#define COMPARE_MEMORY_WINDOWS( F1NPTS, F2NPTS, ERRORMODE, MESSAGE)               \
{                                                                                 \
  FieldT F1(MemoryWindow(F1NPTS), bc, gd, NULL, InternalStorage, memType1);       \
  FieldT F2(MemoryWindow(F2NPTS), bc, gd, NULL, InternalStorage, memType2);       \
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
    f = new FieldT(f1->window_with_ghost(), f1->boundary_info(), f1->get_valid_ghost_data(), NULL, InternalStorage);
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
void sprinkle_in_field(FieldT * f1, typename FieldT::value_type* vals, size_t size)
{
  size_t cells = f1->window_with_ghost().local_npts();
  if(size > cells) {
    throw(std::runtime_error("Too many values to sprinkle"));
  }

  FieldT * f;
  bool created_field = false;
  if(f1->memory_device_type() == EXTERNAL_CUDA_GPU) {
    created_field = true;
    f = new FieldT(f1->window_with_ghost(), f1->boundary_info(), f1->get_valid_ghost_data(), NULL, InternalStorage);
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
                    const bool testFieldNotEqualFunction,
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
  std::numeric_limits<double> nl;

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

  //switch expected_equal and manual equal result based on testFieldNotEqualFunction compare
  if(testFieldNotEqualFunction) {man_equal = !man_equal; expected_equal = !expected_equal;}

  //compare manual to function
  switch(et) {
    case RELATIVE:
      if(testFieldNotEqualFunction) {
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
      if(testFieldNotEqualFunction) {
        return (man_equal == field_not_equal_abs(f1, f2, error)) && (man_equal == expected_equal);
      }
      else {
        return (man_equal == field_equal_abs(f1, f2, error)) && (man_equal == expected_equal);
      }
    case ULP:
      if(testFieldNotEqualFunction) {
        return (man_equal == field_not_equal_ulp(f1, f2, error)) && (man_equal == expected_equal);
      }
      else {
        return (man_equal == field_equal_ulp(f1, f2, error)) && (man_equal == expected_equal);
      }
  }
  return false;
}

template<typename FieldT>
class TestFieldEqual
{
  public:

    const MemoryWindow window;
    const BoundaryCellInfo bc;
    const GhostData gd;
    const int total;

    MemoryType memType1;
    MemoryType memType2;

    TestFieldEqual(const IntVec npts)
      : memType1(LOCAL_RAM), 
      memType2(LOCAL_RAM), 
      window(npts),
      bc(BoundaryCellInfo::build<FieldT>()),
      gd(1),
      total(npts[0]*npts[1]*npts[2])
  {}

    bool compare_sprinkle_fields(
        double const f1Range, double const  f1SprinkleVal, size_t const f1Amount,
        double const f2Range, double const  f2SprinkleVal, size_t const f2Amount,
        double const extremeVal, size_t const extremeAmount,
        ErrorType const & et, 
        double const error, 
        bool const expectedEqual,
        bool const testFieldNotEqualFunction,
        bool const verboseOutput)
    {
      FieldT sf1(window, bc, gd, NULL, InternalStorage, memType1);
      FieldT sf2(window, bc, gd, NULL, InternalStorage, memType2);

      /* Fill in Field 1 Values to Sprinkle */
      fill_field_range(&sf1, 0, f1Range);
      double values1[f1Amount+extremeAmount];
      size_t inc = f1Amount/(extremeAmount + 1);
      int i = 0;
      for(int j = 0; j < extremeAmount; j++) {
        for(int k = 0; k < inc; k++) {
          values1[i] = f1SprinkleVal;
          i++;
        }
        values1[i] = extremeVal;
        i++;
      }
      for(; i < f1Amount+extremeAmount; i++) {
        values1[i] = f1SprinkleVal;
      }
      /*Intialize Field 1*/
      sprinkle_in_field(&sf1, values1, f1Amount+extremeAmount);

      /* Fill in Field 2 Values to Sprinkle */
      fill_field_range(&sf2, 0, f2Range);
      double values2[f2Amount+extremeAmount];
      inc = f2Amount/(extremeAmount + 1);
      i = 0;
      for(int j = 0; j < extremeAmount; j++) {
        for(int k = 0; k < inc; k++) {
          values2[i] = f2SprinkleVal;
          i++;
        }
        values2[i] = extremeVal;
        i++;
      }
      for(; i < f2Amount+extremeAmount; i++) {
        values2[i] = f2SprinkleVal;
      }
      /*Intialize Field 2*/
      sprinkle_in_field(&sf2, values2, f2Amount+extremeAmount);

      return compare_fields(&sf1, &sf2, et, error, expectedEqual, 0, testFieldNotEqualFunction, verboseOutput);
    }

    inline bool compare_fields(double const f1Val,
        double const f2Val,
        ErrorType const et,
        double const error,
        bool const expectedEqual,
        double const absError,
        bool const testFieldNotEqualFunction,
        bool const verboseOutput)
    {
      FieldT f1(window, bc, gd, NULL, InternalStorage, memType1);
      FieldT f2(window, bc, gd, NULL, InternalStorage, memType2);
      f1 <<= f1Val;
      f2 <<= f2Val;
      return compare_fields(&f1, &f2, et, error, expectedEqual, absError, testFieldNotEqualFunction, verboseOutput);
    }

    inline bool compare_fields(FieldT * f1,
        FieldT * f2,
        ErrorType const et,
        double const error,
        bool const expectedEqual,
        double const absError,
        bool const testFieldNotEqualFunction,
        bool const verboseOutput)
    {
      return manual_error_compare(*f1, *f2, error, et, testFieldNotEqualFunction, verboseOutput, expectedEqual, absError);
    }

    bool test(const MemoryType memType1,
        const MemoryType memType2,
        const ErrorType et,
        const bool testFieldNotEqualFunction,
        const bool verboseOutput)
    {
      TestHelper status(verboseOutput);
      std::numeric_limits<double> nl;
      this->memType1 = memType1;
      this->memType2 = memType2;

      FieldT* f1;
      FieldT* f2;

      //local fields
      FieldT lf1(window, bc, gd, NULL, InternalStorage);
      FieldT lf2(window, bc, gd, NULL, InternalStorage);
      FieldT lf3(window, bc, gd, NULL, InternalStorage);

      initialize_field(lf1, 0);
      initialize_field(lf2, 0);
      initialize_field(lf3, total);
#ifdef __CUDACC__
      //gpu fields
      FieldT gf1(window, bc, gd, NULL, InternalStorage, EXTERNAL_CUDA_GPU, 0);
      FieldT gf2(window, bc, gd, NULL, InternalStorage, EXTERNAL_CUDA_GPU, 0);
      FieldT gf3(window, bc, gd, NULL, InternalStorage, EXTERNAL_CUDA_GPU, 0);
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
      status(manual_error_compare(*f1, *f2, 0.0, et, testFieldNotEqualFunction, verboseOutput, true, 0), "Duplicate Fields Equal");

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
      status(manual_error_compare(*f1, *f2, 0.0, et, testFieldNotEqualFunction, verboseOutput, false, 0), "Non-Duplicate Fields Not Equal");


      switch(et) {
        case RELATIVE:
          //relative error test
          status(compare_fields(3.0, 7.49, RELATIVE, 1.5, true, 0, testFieldNotEqualFunction, verboseOutput), "Off By 150% (Equal)");
          status(compare_fields(3.0, 7.51, RELATIVE, 1.5, false, 0, testFieldNotEqualFunction, verboseOutput), "Off By 150% (Not Equal)");
          status(compare_fields(3.0, 3.98, RELATIVE, .33, true, 0, testFieldNotEqualFunction, verboseOutput), "Off By 33% (Equal)");
          status(compare_fields(3.0, 4.10, RELATIVE, .33, false, 0, testFieldNotEqualFunction, verboseOutput), "Off By 33% (Not Equal)");
          status(compare_fields(3.0, 2.97, RELATIVE, .01, true, 0, testFieldNotEqualFunction, verboseOutput), "Off By 1% (Equal)");
          status(compare_fields(3.0, 2.96, RELATIVE, .01, false, 0, testFieldNotEqualFunction, verboseOutput), "Off By 1% (Not Equal)");

          //near zero value tests
          status(compare_sprinkle_fields(1, 0, 8, 1, 1e-9, 8, 0, 0, RELATIVE, .01, true, testFieldNotEqualFunction, verboseOutput),
              "Near Zero Field Range [-1, 1] Off By 1% (Equal)");
          status(compare_sprinkle_fields(1, 0, 8, 1, 1e-8, 8, 0, 0, RELATIVE, .01, false, testFieldNotEqualFunction, verboseOutput),
              "Near Zero Field Range [-1, 1] Off By 1% (Not Equal)");
          status(compare_sprinkle_fields(10, 0, 8, 10, 1e-8, 8, 0, 0, RELATIVE, .01, true, testFieldNotEqualFunction, verboseOutput),
              "Near Zero Field Range [-10, 10] Off By 1% (Equal)");
          status(compare_sprinkle_fields(10, 0, 8, 10, 1e-7, 8, 0, 0, RELATIVE, .01, false, testFieldNotEqualFunction, verboseOutput),
              "Near Zero Field Range [-10, 10] Off By 1% (Not Equal)");
          status(compare_sprinkle_fields(100000, 0, 8, 100000, 1e-4, 8, 0, 0, RELATIVE, .01, true, testFieldNotEqualFunction, verboseOutput),
              "Near Zero Field Range [-100000, 100000] Off By 1% (Equal)");
          status(compare_sprinkle_fields(100000, 0, 8, 100000, 1e-3, 8, 0, 0, RELATIVE, .01, false, testFieldNotEqualFunction, verboseOutput),
              "Near Zero Field Range [-100000, 100000] Off By 1% (Not Equal)");
          status(compare_sprinkle_fields(1, 0, 8, 1, 1e-7, 8, 1000, 1, RELATIVE, .01, true, testFieldNotEqualFunction, verboseOutput),
              "Outlier Near Zero Field Range [-1, 1] Off By 1% (Equal)");
          status(compare_sprinkle_fields(1, 0, 8, 1, 1e-6, 8, 1000, 1, RELATIVE, .01, false, testFieldNotEqualFunction, verboseOutput),
              "Outlier Near Zero Field Range [-1, 1] Off By 1% (Not Equal)");
          status(compare_sprinkle_fields(1, 0, 8, 1, 1e-7, 8, 1000, 5, RELATIVE, .01, true, testFieldNotEqualFunction, verboseOutput),
              "Five Outlier Near Zero Field Range [-1, 1] Off By 1% (Equal)");
          status(compare_sprinkle_fields(1, 0, 8, 1, 1e-6, 8, 1000, 5, RELATIVE, .01, false, testFieldNotEqualFunction, verboseOutput),
              "Five Outlier Near Zero Field Range [-1, 1] Off By 1% (Not Equal)");

          //Custom Absolute Error
          status(compare_fields(1.0, 3.0, RELATIVE, 1.0, true, 1.0, testFieldNotEqualFunction, verboseOutput), "Absolute Tolerance 1: 1-3 Off By 100% (Equal)");
          status(compare_fields(1.0, 4.0, RELATIVE, 1.0, false, 1.0, testFieldNotEqualFunction, verboseOutput), "Absolute Tolerance 1: 1-4 Off By 100% (Not Equal)");
          status(compare_fields(1.0, 5.0, RELATIVE, 2.0, true, 1.0, testFieldNotEqualFunction, verboseOutput), "Absolute Tolerance 1: 1-5 Off By 200% (Equal)");
          status(compare_fields(1.0, 5.0, RELATIVE, 1.0, false, 1.0, testFieldNotEqualFunction, verboseOutput), "Absolute Tolerance 1: 1-5 Off By 100% (Not Equal)");
          break;
        case ABSOLUTE:
          //absolute error test
          status(compare_fields(1.0, 1.0 + nl.epsilon(), ABSOLUTE, nl.epsilon(), true, 0, testFieldNotEqualFunction, verboseOutput), "Off By epsilon (Equal)");
          status(compare_fields(1.0, 1.0 + 2*nl.epsilon(), ABSOLUTE, nl.epsilon(), false, 0, testFieldNotEqualFunction, verboseOutput), "Off By epsilon (Not Equal)");
          status(compare_fields(0.033, 0.024, ABSOLUTE, std::pow(10.0,-2), true, 0, testFieldNotEqualFunction, verboseOutput), "Off By 10^-2 (Equal)");
          status(compare_fields(0.033, 0.022, ABSOLUTE, std::pow(10.0,-2), false, 0, testFieldNotEqualFunction, verboseOutput), "Off By 10^-2 (Not Equal)");
          status(compare_fields(4679000.0, 4680000.0, ABSOLUTE, std::pow(10.0,3), true, 0, testFieldNotEqualFunction, verboseOutput), "Off By 10^3 (Equal)");
          status(compare_fields(4679000.0, 4681000.0, ABSOLUTE, std::pow(10.0,3), false, 0, testFieldNotEqualFunction, verboseOutput), "Off By 10^3 (Not Equal)");
          status(compare_fields(4679000.0, 6890330.0, ABSOLUTE, 11569300.0, true, 0, testFieldNotEqualFunction, verboseOutput), "Large Number Check, Exact Error (Equal)");
          break;
        case ULP:
          using boost::math::float_prior;
          using boost::math::float_next;
          using boost::math::float_advance;

          //near zero value tests
          status(compare_fields(0.0, float_next(0.0), ULP, 1, true, 0, testFieldNotEqualFunction, verboseOutput), "Near Zero Off By 1 ulp (Equal)");
          status(compare_fields(0.0, float_advance(0.0, 2), ULP, 1, false, 0, testFieldNotEqualFunction, verboseOutput), "Near Zero Off By 1 ulp (Not Equal)");
          status(compare_fields(0.0, -float_advance(0.0, 2), ULP, 2, true, 0, testFieldNotEqualFunction, verboseOutput), "Near Zero Off By 2 ulps (Equal)");
          status(compare_fields(0.0, -float_advance(0.0, 3), ULP, 2, false, 0, testFieldNotEqualFunction, verboseOutput), "Near Zero Off By 2 ulps (Not Equal)");
          status(compare_fields(-float_advance(0.0, 1), 0.0, ULP, 1, true, 0, testFieldNotEqualFunction, verboseOutput), "Near Zero Reversed Off By 1 ulp (Equal)");
          status(compare_fields(float_advance(0.0, 2), 0.0, ULP, 1, false, 0, testFieldNotEqualFunction, verboseOutput), "Near Zero Reversed Off By 1 ulp (Not Equal)");

          //machine epsilon
          status(compare_fields(1.0, 1.0+nl.epsilon(), ULP, 1, true, 0, testFieldNotEqualFunction, verboseOutput), "Machine epsilon at 1.0 is 1 ulp (Equal)");
          status(compare_fields(1.0, 1.0+nl.epsilon(), ULP, 0, false, 0, testFieldNotEqualFunction, verboseOutput), "Machine epsilon at 1.0 is not 0 ulps (Not Equal)");

          //ulps error test
          status(compare_fields(3.0, float_prior(3.0), ULP, 1, true, 0, testFieldNotEqualFunction, verboseOutput), "Off By 1 ulp (Equal)");
          status(compare_fields(3.0, float_advance(3.0, -2), ULP, 1, false, 0, testFieldNotEqualFunction, verboseOutput), "Off By 1 ulp (Not Equal)");
          status(compare_fields(3.0, float_advance(3.0, 3), ULP, 3, true, 0, testFieldNotEqualFunction, verboseOutput), "Off By 3 ulps (Equal)");
          status(compare_fields(3.0, float_advance(3.0, 4), ULP, 3, false, 0, testFieldNotEqualFunction, verboseOutput), "Off By 3 ulps (Not Equal)");
          status(compare_fields(3.0, float_advance(3.0, 20), ULP, 20, true, 0, testFieldNotEqualFunction, verboseOutput), "Off By 20 ulps (Equal)");
          status(compare_fields(3.0, float_advance(3.0, -21), ULP, 20, false, 0, testFieldNotEqualFunction, verboseOutput), "Off By 20 ulps (Not Equal)");

          //limits test
          status(compare_fields(nl.max(), float_advance(nl.max(), -1), ULP, 1, true, 0, testFieldNotEqualFunction, verboseOutput), "Max Double Off By -1 ulp (Equal)");
          status(compare_fields(nl.min(), float_advance(nl.min(), 1), ULP, 1, true, 0, testFieldNotEqualFunction, verboseOutput), "Min > 0 Off By 1 ulp (Equal)");
          status(compare_fields(nl.min(), float_advance(nl.min(), -1), ULP, 1, true, 0, testFieldNotEqualFunction, verboseOutput), "Min > 0 Off By -1 ulp (Equal)");
          status(compare_fields(-nl.max(), float_advance(-nl.max(), 1), ULP, 1, true, 0, testFieldNotEqualFunction, verboseOutput), "Min Off By 1 ulps (Equal)");
          try {
            //manual compare useless
            status(compare_fields(nl.max(), nl.infinity(), ULP, 20, true, 0, testFieldNotEqualFunction, verboseOutput), "fail");
            status(false, "Max Double = Infinity By 20 ulps (Throws Exception)");
          } catch(std::domain_error) {status(true,  "Max Double = Infinity By 20 ulps (Throws Exception)");}
          try {
            //manual compare useless
            status(compare_fields(-nl.max(), -nl.infinity(), ULP, 20, true, 0, testFieldNotEqualFunction, verboseOutput), "fail");
            status(false, "Min Double = Infinity By 20 ulps (Throws Exception)");
          } catch(std::domain_error) {status(true, "Min Double = Infinity By 20 ulps (Throws Exception)");}
          try {
            //manual compare useless
            status(compare_fields(nan(""), nan(""), ULP, 20, true, 0, testFieldNotEqualFunction, verboseOutput), "fail");
            status(false, "NAN = NAN By 20 ulps (Throws Exception)");
          } catch(std::domain_error) {status(true, "NAN = NAN By 20 ulps (Throws Exception)");}
          break;
      }

      //unequal memory window tests
      IntVec npts(window.glob_dim());
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
};

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
                    const bool testFieldNotEqualFunction,
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
  std::numeric_limits<double> nl;

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

  //switch expected_equal and manual equal result based on testFieldNotEqualFunction compare
  if(testFieldNotEqualFunction) {man_equal = !man_equal; expected_equal = !expected_equal;}

  //compare manual to function
  switch(et) {
    case RELATIVE:
      if(testFieldNotEqualFunction) {
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
      if(testFieldNotEqualFunction) {
        return (man_equal == field_not_equal_abs(d, f1, error)) && (man_equal == expected_equal);
      }
      else {
        return (man_equal == field_equal_abs(d, f1, error)) && (man_equal == expected_equal);
      }
    case ULP:
      if(testFieldNotEqualFunction) {
        return (man_equal == field_not_equal_ulp(d, f1, error)) && (man_equal == expected_equal);
      }
      else {
        return (man_equal == field_equal_ulp(d, f1, error)) && (man_equal == expected_equal);
      }
  }
  return false;
}

template<typename FieldT>
class TestFieldEqualScalar
{
  private:
    const MemoryWindow window;
    const BoundaryCellInfo bc;
    const GhostData gd;
    const int total;

    MemoryType memType1;

  public:
    TestFieldEqualScalar(const IntVec npts)
      : memType1(LOCAL_RAM), window(npts), bc(BoundaryCellInfo::build<FieldT>()), 
      gd(1), total(npts[0]*npts[1]*npts[2])
    {}

    inline bool compare_field_scalar(double const val,
        double const f1Val,
        ErrorType const et,
        double const error,
        bool const expectedEqual,
        double const absError,
        bool const testFieldNotEqualFunction,
        bool const verboseOutput)
    {
      FieldT f1(window, bc, gd, NULL, InternalStorage, this->memType1);
      f1 <<= f1Val;
      return manual_error_compare(val, f1, error, et, testFieldNotEqualFunction, verboseOutput, expectedEqual, absError);
    }
    bool test(const MemoryType memType1, const ErrorType et, const bool testFieldNotEqualFunction, const bool verboseOutput)
    {
      TestHelper status(verboseOutput);
      std::numeric_limits<double> nl;
      this->memType1 = memType1;

      //field exactly equal
      status(compare_field_scalar(42.0, 42.0, et, 0.0, true, 0, testFieldNotEqualFunction, verboseOutput), "Duplicate Fields Equal");

      //field not equal
      status(compare_field_scalar(21.0, 42.0, et, 0.0, false, 0, testFieldNotEqualFunction, verboseOutput), "Non-Duplicate Fields Not Equal");


      switch(et) {
        case RELATIVE:
          //relative error test
          status(compare_field_scalar(3.0, 7.49, RELATIVE, 1.5, true, 0, testFieldNotEqualFunction, verboseOutput), "Off By 150% (Equal)");
          status(compare_field_scalar(3.0, 7.51, RELATIVE, 1.5, false, 0, testFieldNotEqualFunction, verboseOutput), "Off By 150% (Not Equal)");
          status(compare_field_scalar(3.0, 3.98, RELATIVE, .33, true, 0, testFieldNotEqualFunction, verboseOutput), "Off By 33% (Equal)");
          status(compare_field_scalar(3.0, 4.10, RELATIVE, .33, false, 0, testFieldNotEqualFunction, verboseOutput), "Off By 33% (Not Equal)");
          status(compare_field_scalar(3.0, 2.97, RELATIVE, .01, true, 0, testFieldNotEqualFunction, verboseOutput), "Off By 1% (Equal)");
          status(compare_field_scalar(3.0, 2.96, RELATIVE, .01, false, 0, testFieldNotEqualFunction, verboseOutput), "Off By 1% (Not Equal)");

          //Custom Absolute Error
          status(compare_field_scalar(1.0, 3.0, RELATIVE, 1.0, true, 1.0, testFieldNotEqualFunction, verboseOutput), "Absolute Tolerance 1: 1-3 Off By 100% (Equal)");
          status(compare_field_scalar(1.0, 4.0, RELATIVE, 1.0, false, 1.0, testFieldNotEqualFunction, verboseOutput), "Absolute Tolerance 1: 1-4 Off By 100% (Not Equal)");
          status(compare_field_scalar(1.0, 5.0, RELATIVE, 2.0, true, 1.0, testFieldNotEqualFunction, verboseOutput), "Absolute Tolerance 1: 1-5 Off By 200% (Equal)");
          status(compare_field_scalar(1.0, 5.0, RELATIVE, 1.0, false, 1.0, testFieldNotEqualFunction, verboseOutput), "Absolute Tolerance 1: 1-5 Off By 100% (Not Equal)");
          break;
        case ABSOLUTE:
          //absolute error test
          status(compare_field_scalar(1.0, 1.0 + nl.epsilon(), ABSOLUTE, nl.epsilon(), true, 0, testFieldNotEqualFunction, verboseOutput), "Off By epsilon (Equal)");
          status(compare_field_scalar(1.0, 1.0 + 2*nl.epsilon(), ABSOLUTE, nl.epsilon(), false, 0, testFieldNotEqualFunction, verboseOutput), "Off By epsilon (Not Equal)");
          status(compare_field_scalar(0.033, 0.024, ABSOLUTE, std::pow(10.0,-2), true, 0, testFieldNotEqualFunction, verboseOutput), "Off By 10^-2 (Equal)");
          status(compare_field_scalar(0.033, 0.022, ABSOLUTE, std::pow(10.0,-2), false, 0, testFieldNotEqualFunction, verboseOutput), "Off By 10^-2 (Not Equal)");
          status(compare_field_scalar(4679000.0, 4680000.0, ABSOLUTE, std::pow(10.0,3), true, 0, testFieldNotEqualFunction, verboseOutput), "Off By 10^3 (Equal)");
          status(compare_field_scalar(4679000.0, 4681000.0, ABSOLUTE, std::pow(10.0,3), false, 0, testFieldNotEqualFunction, verboseOutput), "Off By 10^3 (Not Equal)");
          status(compare_field_scalar(4679000.0, 6890330.0, ABSOLUTE, 11569300.0, true, 0, testFieldNotEqualFunction, verboseOutput), "Large Number Check, Exact Error (Equal)");
          break;
        case ULP:
          using boost::math::float_prior;
          using boost::math::float_next;
          using boost::math::float_advance;

          //near zero value tests
          status(compare_field_scalar(0.0, float_next(0.0), ULP, 1, true, 0, testFieldNotEqualFunction, verboseOutput), "Near Zero Off By 1 ulp (Equal)");
          status(compare_field_scalar(0.0, float_advance(0.0, 2), ULP, 1, false, 0, testFieldNotEqualFunction, verboseOutput), "Near Zero Off By 1 ulp (Not Equal)");
          status(compare_field_scalar(0.0, -float_advance(0.0, 2), ULP, 2, true, 0, testFieldNotEqualFunction, verboseOutput), "Near Zero Off By 2 ulps (Equal)");
          status(compare_field_scalar(0.0, -float_advance(0.0, 3), ULP, 2, false, 0, testFieldNotEqualFunction, verboseOutput), "Near Zero Off By 2 ulps (Not Equal)");
          status(compare_field_scalar(-float_advance(0.0, 1), 0.0, ULP, 1, true, 0, testFieldNotEqualFunction, verboseOutput), "Near Zero Reversed Off By 1 ulp (Equal)");
          status(compare_field_scalar(float_advance(0.0, 2), 0.0, ULP, 1, false, 0, testFieldNotEqualFunction, verboseOutput), "Near Zero Reversed Off By 1 ulp (Not Equal)");

          //machine epsilon
          status(compare_field_scalar(1.0, 1.0+nl.epsilon(), ULP, 1, true, 0, testFieldNotEqualFunction, verboseOutput), "Machine epsilon at 1.0 is 1 ulp (Equal)");
          status(compare_field_scalar(1.0, 1.0+nl.epsilon(), ULP, 0, false, 0, testFieldNotEqualFunction, verboseOutput), "Machine epsilon at 1.0 is not 0 ulps (Not Equal)");

          //ulps error test
          status(compare_field_scalar(3.0, float_prior(3.0), ULP, 1, true, 0, testFieldNotEqualFunction, verboseOutput), "Off By 1 ulp (Equal)");
          status(compare_field_scalar(3.0, float_advance(3.0, -2), ULP, 1, false, 0, testFieldNotEqualFunction, verboseOutput), "Off By 1 ulp (Not Equal)");
          status(compare_field_scalar(3.0, float_advance(3.0, 3), ULP, 3, true, 0, testFieldNotEqualFunction, verboseOutput), "Off By 3 ulps (Equal)");
          status(compare_field_scalar(3.0, float_advance(3.0, 4), ULP, 3, false, 0, testFieldNotEqualFunction, verboseOutput), "Off By 3 ulps (Not Equal)");
          status(compare_field_scalar(3.0, float_advance(3.0, 20), ULP, 20, true, 0, testFieldNotEqualFunction, verboseOutput), "Off By 20 ulps (Equal)");
          status(compare_field_scalar(3.0, float_advance(3.0, -21), ULP, 20, false, 0, testFieldNotEqualFunction, verboseOutput), "Off By 20 ulps (Not Equal)");

          //limits test
          status(compare_field_scalar(nl.max(), float_advance(nl.max(), -1), ULP, 1, true, 0, testFieldNotEqualFunction, verboseOutput), "Max Double Off By -1 ulp (Equal)");
          status(compare_field_scalar(nl.min(), float_advance(nl.min(), 1), ULP, 1, true, 0, testFieldNotEqualFunction, verboseOutput), "Min > 0 Off By 1 ulp (Equal)");
          status(compare_field_scalar(nl.min(), float_advance(nl.min(), -1), ULP, 1, true, 0, testFieldNotEqualFunction, verboseOutput), "Min > 0 Off By -1 ulp (Equal)");
          status(compare_field_scalar(-nl.max(), float_advance(-nl.max(), 1), ULP, 1, true, 0, testFieldNotEqualFunction, verboseOutput), "Min Off By 1 ulps (Equal)");
          try {
            //manual compare useless
            status(compare_field_scalar(nl.max(), nl.infinity(), ULP, 20, true, 0, testFieldNotEqualFunction, verboseOutput), "fail");
            status(false, "Max Double = Infinity By 20 ulps (Throws Exception)");
          } catch(std::domain_error) {status(true,  "Max Double = Infinity By 20 ulps (Throws Exception)");}
          try {
            //manual compare useless
            status(compare_field_scalar(-nl.max(), -nl.infinity(), ULP, 20, true, 0, testFieldNotEqualFunction, verboseOutput), "fail");
            status(false, "Min Double = Infinity By 20 ulps (Throws Exception)");
          } catch(std::domain_error) {status(true, "Min Double = Infinity By 20 ulps (Throws Exception)");}
          try {
            //manual compare useless
            status(compare_field_scalar(nan(""), nan(""), ULP, 20, true, 0, testFieldNotEqualFunction, verboseOutput), "fail");
            status(false, "NAN = NAN By 20 ulps (Throws Exception)");
          } catch(std::domain_error) {status(true, "NAN = NAN By 20 ulps (Throws Exception)");}
          break;
      }

      return status.ok();
    }
};




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
    TestFieldEqual<SVolField> tfe(winSize);
    //test field_equal
    overall(tfe.test(LOCAL_RAM, LOCAL_RAM, RELATIVE, false, fieldEqualVerbose), "LOCAL_RAM x LOCAL_RAM Relative Equal Test");
#ifdef __CUDACC__
    overall(tfe.test(LOCAL_RAM, EXTERNAL_CUDA_GPU, RELATIVE, false, fieldEqualVerbose), "LOCAL_RAM x EXTERNAL_CUDA_GPU Relative Equal Test");
    overall(tfe.test(EXTERNAL_CUDA_GPU, LOCAL_RAM, RELATIVE, false, fieldEqualVerbose), "EXTERNAL_CUDA_GPU x LOCAL_RAM Relative Equal Test");
    overall(tfe.test(EXTERNAL_CUDA_GPU, EXTERNAL_CUDA_GPU, RELATIVE, false, fieldEqualVerbose), "EXTERNAL_CUDA_GPU x EXTERNAL_CUDA_GPU Relative Equal Test");
#endif

    //test field_equal_abs
    overall(tfe.test(LOCAL_RAM, LOCAL_RAM, ABSOLUTE, false, fieldEqualVerbose), "LOCAL_RAM x LOCAL_RAM Absolute Error Equal Test");
#ifdef __CUDACC__
    overall(tfe.test(LOCAL_RAM, EXTERNAL_CUDA_GPU, ABSOLUTE, false, fieldEqualVerbose), "LOCAL_RAM x EXTERNAL_CUDA_GPU Absolute Error Equal Test");
    overall(tfe.test(EXTERNAL_CUDA_GPU, LOCAL_RAM, ABSOLUTE, false, fieldEqualVerbose), "EXTERNAL_CUDA_GPU x LOCAL_RAM Absolute Error Equal Test");
    overall(tfe.test(EXTERNAL_CUDA_GPU, EXTERNAL_CUDA_GPU, ABSOLUTE, false, fieldEqualVerbose), "EXTERNAL_CUDA_GPU x EXTERNAL_CUDA_GPU Absolute Error Equal Test");
#endif

    //test field_equal_ulp
    overall(tfe.test(LOCAL_RAM, LOCAL_RAM, ULP, false, fieldEqualVerbose), "LOCAL_RAM x LOCAL_RAM Ulps Equal Test");
#ifdef __CUDACC__
    overall(tfe.test(LOCAL_RAM, EXTERNAL_CUDA_GPU, ULP, false, fieldEqualVerbose), "LOCAL_RAM x EXTERNAL_CUDA_GPU Ulps Equal Test");
    overall(tfe.test(EXTERNAL_CUDA_GPU, LOCAL_RAM, ULP, false, fieldEqualVerbose), "EXTERNAL_CUDA_GPU x LOCAL_RAM Ulps Equal Test");
    overall(tfe.test(EXTERNAL_CUDA_GPU, EXTERNAL_CUDA_GPU, ULP, false, fieldEqualVerbose), "EXTERNAL_CUDA_GPU x EXTERNAL_CUDA_GPU Ulps Equal Test");
#endif

    //test field_not_equal
    overall(tfe.test(LOCAL_RAM, LOCAL_RAM, RELATIVE, true, fieldEqualVerbose), "LOCAL_RAM x LOCAL_RAM Relative Not Equal Test");
#ifdef __CUDACC__
    overall(tfe.test(LOCAL_RAM, EXTERNAL_CUDA_GPU, RELATIVE, true, fieldEqualVerbose), "LOCAL_RAM x EXTERNAL_CUDA_GPU Relative Not Equal Test");
    overall(tfe.test(EXTERNAL_CUDA_GPU, LOCAL_RAM, RELATIVE, true, fieldEqualVerbose), "EXTERNAL_CUDA_GPU x LOCAL_RAM Relative Not Equal Test");
    overall(tfe.test(EXTERNAL_CUDA_GPU, EXTERNAL_CUDA_GPU, RELATIVE, true, fieldEqualVerbose), "EXTERNAL_CUDA_GPU x EXTERNAL_CUDA_GPU Relative Not Equal Test");
#endif

    //test field_not_equal_abs
    overall(tfe.test(LOCAL_RAM, LOCAL_RAM, ABSOLUTE, true, fieldEqualVerbose), "LOCAL_RAM x LOCAL_RAM Absolute Error Not Equal Test");
#ifdef __CUDACC__
    overall(tfe.test(LOCAL_RAM, EXTERNAL_CUDA_GPU, ABSOLUTE, true, fieldEqualVerbose), "LOCAL_RAM x EXTERNAL_CUDA_GPU Absolute Error Not Equal Test");
    overall(tfe.test(EXTERNAL_CUDA_GPU, LOCAL_RAM, ABSOLUTE, true, fieldEqualVerbose), "EXTERNAL_CUDA_GPU x LOCAL_RAM Absolute Error Not Equal Test");
    overall(tfe.test(EXTERNAL_CUDA_GPU, EXTERNAL_CUDA_GPU, ABSOLUTE, true, fieldEqualVerbose), "EXTERNAL_CUDA_GPU x EXTERNAL_CUDA_GPU Absolute Error Not Equal Test");
#endif

    //test field_not_equal_ulp
    overall(tfe.test(LOCAL_RAM, LOCAL_RAM, ULP, true, fieldEqualVerbose), "LOCAL_RAM x LOCAL_RAM Ulps Not Equal Test");
#ifdef __CUDACC__
    overall(tfe.test(LOCAL_RAM, EXTERNAL_CUDA_GPU, ULP, true, fieldEqualVerbose), "LOCAL_RAM x EXTERNAL_CUDA_GPU Ulps Not Equal Test");
    overall(tfe.test(EXTERNAL_CUDA_GPU, LOCAL_RAM, ULP, true, fieldEqualVerbose), "EXTERNAL_CUDA_GPU x LOCAL_RAM Ulps Not Equal Test");
    overall(tfe.test(EXTERNAL_CUDA_GPU, EXTERNAL_CUDA_GPU, ULP, true, fieldEqualVerbose), "EXTERNAL_CUDA_GPU x EXTERNAL_CUDA_GPU Ulps Not Equal Test");
#endif
  }


  //test field comparison functions with a scalar value
  {
    //create testing object with basic state
    TestFieldEqualScalar<SVolField> tfes(winSize);

    //test field_equal with scalar
    overall(tfes.test(LOCAL_RAM, RELATIVE, false, fieldEqualVerbose), "SCALAR x LOCAL_RAM Relative Equal Test");
#ifdef __CUDACC__
    overall(tfes.test(EXTERNAL_CUDA_GPU, RELATIVE, false, fieldEqualVerbose), "SCALAR x EXTERNAL_CUDA_GPU Relative Equal Test");
#endif

    //test field_equal_abs with scalar
    overall(tfes.test(LOCAL_RAM, ABSOLUTE, false, fieldEqualVerbose), "SCALAR x LOCAL_RAM Absolute Error Equal Test");
#ifdef __CUDACC__
    overall(tfes.test(EXTERNAL_CUDA_GPU, ABSOLUTE, false, fieldEqualVerbose), "SCALAR x EXTERNAL_CUDA_GPU Absolute Error Equal Test");
#endif

    //test field_equal_ulp with scalar
    overall(tfes.test(LOCAL_RAM, ULP, false, fieldEqualVerbose), "SCALAR x LOCAL_RAM Ulps Equal Test");
#ifdef __CUDACC__
    overall(tfes.test(EXTERNAL_CUDA_GPU, ULP, false, fieldEqualVerbose), "SCALAR x EXTERNAL_CUDA_GPU Ulps Equal Test");
#endif

    //test field_not_equal with scalar
    overall(tfes.test(LOCAL_RAM, RELATIVE, true, fieldEqualVerbose), "SCALAR x LOCAL_RAM Relative Not Equal Test");
#ifdef __CUDACC__
    overall(tfes.test(EXTERNAL_CUDA_GPU, RELATIVE, true, fieldEqualVerbose), "SCALAR x EXTERNAL_CUDA_GPU Relative Not Equal Test");
#endif

    //test field_not_equal_abs with scalar
    overall(tfes.test(LOCAL_RAM, ABSOLUTE, true, fieldEqualVerbose), "SCALAR x LOCAL_RAM Absolute Error Not Equal Test");
#ifdef __CUDACC__
    overall(tfes.test(EXTERNAL_CUDA_GPU, ABSOLUTE, true, fieldEqualVerbose), "SCALAR x EXTERNAL_CUDA_GPU Absolute Error Not Equal Test");
#endif

    //test field_not_equal_ulp with scalar
    overall(tfes.test(LOCAL_RAM, ULP, true, fieldEqualVerbose), "SCALAR x LOCAL_RAM Ulps Not Equal Test");
#ifdef __CUDACC__
    overall(tfes.test(EXTERNAL_CUDA_GPU, ULP, true, fieldEqualVerbose), "SCALAR x EXTERNAL_CUDA_GPU Ulps Not Equal Test");
#endif
  }

  if( overall.isfailed() ){
    std::cout << "FAIL!" << std::endl;
    return -1;
  }
  std::cout << "PASS" << std::endl;
  return 0;
}
