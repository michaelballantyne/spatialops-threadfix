/*
 * Copyright (c) 2014 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */
#ifndef SpatialOps_FieldComparisons_h
#define SpatialOps_FieldComparisons_h

#include <spatialops/Nebo.h>
#include <sstream>

#define FIELDCOMPARISONS_ABS_ERROR_CONST .000001

/**
 * @brief Comparison operators
 * WARNING: Slow in general and comparison with external fields will incur copy penalties.
 */
namespace SpatialOps{
namespace structured{

template<typename FieldT>
class FieldComparisonHelper;

/**
 * @brief Returns if f1 is element-wise not equal to f2 within a certain relative
 * tolerance.
 *
 * This function simply calls field_equal and negates it.
 * \c return !field_equal(f1, f2, error, error_abs);
 * error_abs is defined as default to be the L2 norm of \c f1 multiplied by \c error
 * and \c FIELDCOMPARISONS_ABS_ERROR_CONST.
 *
 * WARNING: Undefined behavior if f1 is a field of all 0's.
 *
 * WARNING: Slow in general and comparison with external fields will incur copy penalties.
 *
 * @tparam FieldT -- Any type of SpatialField
 * @param f1 -- Field 1
 * @param f2 -- Field 2
 * @param error -- Allowable percentage of error. i.e. 1% = .01
 */
template<typename FieldT>
bool field_not_equal(const FieldT& f1, const FieldT& f2, double error=0.0) {
  double error_abs = error ? nebo_norm(f1)*error*FIELDCOMPARISONS_ABS_ERROR_CONST : 0;
  return !field_equal(f1, f2, error, error_abs);
}

/**
 * @brief Returns if f1 is element-wise not equal to f2 within a certain relative
 * tolerance.
 *
 * This function simply calls field_equal and negates it.
 * \c return !field_equal(f1, f2, error, error_abs);
 *
 * WARNING: Slow in general and comparison with external fields will incur copy penalties.
 *
 * @tparam FieldT -- Any type of SpatialField
 * @param f1 -- Field 1
 * @param f2 -- Field 2
 * @param error -- Allowable percentage of error. i.e. 1% = .01
 * @param error_abs -- Allowable absolute error passed on to \c field_equal
 */
template<typename FieldT>
bool field_not_equal(const FieldT& f1, const FieldT& f2, double error, const double error_abs) {
  return !field_equal(f1, f2, error, error_abs);
}
//------------------------------------------------------------------

/**
 * @brief Returns if f1 is element-wise equal to f2 within a certain relative
 * tolerance.
 *
 * This function returns the result of |f1 - f2|/(error_abs + |f1|) > error element wise.
 * error_abs is defined as default to be the L2 norm of \c f1 multiplied by \c error
 * and \c FIELDCOMPARISONS_ABS_ERROR_CONST.
 *
 * WARNING: Undefined behavior if f1 is a field of all 0's.
 *
 * WARNING: Slow in general and comparison with external fields will incur copy penalties.
 *
 * @tparam FieldT -- Any type of SpatialField
 * @param f1 -- Field 1
 * @param f2 -- Field 2
 * @param error -- Allowable percentage of error. i.e. 1% = .01
 */
template<typename FieldT>
bool field_equal(const FieldT& f1, const FieldT& f2, double error=0.0)
{
  double error_abs = error ? nebo_norm(f1)*error*FIELDCOMPARISONS_ABS_ERROR_CONST : 0;
  return field_equal(f1, f2, error, error_abs);
}

/**
 * @brief Returns if f1 is element-wise equal to f2 within a certain relative
 * tolerance.
 *
 * This function returns the result of |f1 - f2|/(error_abs + |f1|) > error element wise.
 *
 * WARNING: Slow in general and comparison with external fields will incur copy penalties.
 *
 * @tparam FieldT -- Any type of SpatialField
 * @param f1 -- Field 1
 * @param f2 -- Field 2
 * @param error -- Allowable percentage of error. i.e. 1% = .01
 * @param error_abs -- Allowable absolute error.  This term becomes significant
 * in the calculation as f1 approaches zero.
 */
template<typename FieldT>
bool field_equal(const FieldT& f1, const FieldT& f2, double error, const double error_abs)
{
  MemoryWindow w1(f1.window_with_ghost());
  MemoryWindow w2(f2.window_with_ghost());
  BoundaryCellInfo bc1(f1.boundary_info());
  BoundaryCellInfo bc2(f2.boundary_info());
  GhostData gd1(f1.get_valid_ghost_data());
  GhostData gd2(f2.get_valid_ghost_data());

  if(w1 != w2) {
    throw( std::runtime_error( "Attempted comparison between fields of unequal size." ) );
  }

  error = std::abs(error);
  bool exact_comparison = error == 0.0;
  FieldT *temp1 = 0;
  FieldT *temp2 = 0;
  typename FieldT::const_iterator *if1 = 0;
  typename FieldT::const_iterator *iend = 0;
  typename FieldT::const_iterator *if2 = 0;

  //initialize if1 and iend
  if(!(if1 = FieldComparisonHelper<FieldT>::init_iterator(f1, w1, bc1, gd1, &temp1))) {
    std::ostringstream msg;
    msg << "Attempted comparison operation on unsupported memory type, at "
      << __FILE__
      << " : " << __LINE__ << std::endl;
    msg << "\t - "
      << "Attempted to compare with a field in "
      << SpatialOps::DeviceTypeTools::get_memory_type_description(
          f1.device_index()) << std::endl;
    throw(std::runtime_error(msg.str()));
  }
  if(temp1)
    iend = new typename FieldT::const_iterator(temp1->end());
  else
    iend = new typename FieldT::const_iterator(f1.end());


  //initialize if2
  if(!(if2 = FieldComparisonHelper<FieldT>::init_iterator(f2, w2, bc2, gd2, &temp2))) {
    std::ostringstream msg;
    msg << "Attempted comparison operation on unsupported memory type, at "
      << __FILE__
      << " : " << __LINE__ << std::endl;
    msg << "\t - "
      << "Attempted to compare with a field in "
      << SpatialOps::DeviceTypeTools::get_memory_type_description(
          f2.device_index()) << std::endl;
    throw(std::runtime_error(msg.str()));
  } //memory leak if runtime error thrown

  //do comparison
  bool result = true;
  for (; *if1 != *iend; ++(*if1), ++(*if2)) {
    if(exact_comparison) {
      if (**if1 != **if2) {
        result = false;
        break;
      }
    }
    else {
      const double denom = std::abs(**if1) + error_abs;

      if (std::abs(**if1 - **if2)/denom > error) {
        result =  false;
        break;
      }
    }
  }

  if(temp1) delete temp1;
  if(temp2) delete temp2;
  delete if1;
  delete iend;
  delete if2;

  return result;
}
//------------------------------------------------------------------

/**
 * @brief Returns if f1 is element-wise not equal to f2 within a certain absolute
 * tolerance.
 *
 * This function simply returns the negated result of field_equal_abs
 *
 * WARNING: Slow in general and comparison with external fields will incur copy penalties.
 *
 * @tparam FieldT -- Any type of SpatialField
 * @param f1 -- Field 1
 * @param f2 -- Field 2
 * @param error -- Allowable absolute value of error.
 */
template<typename FieldT>
bool field_not_equal_abs(const FieldT& f1, const FieldT& f2, double error=0.0) {
  return !field_equal_abs(f1, f2, error);
}
//------------------------------------------------------------------

/**
 * @brief Returns if f1 is element-wise equal to f2 within a certain absolute
 * tolerance.
 *
 * WARNING: Slow in general and comparison with external fields will incur copy penalties.
 *
 * @tparam FieldT -- Any type of SpatialField
 * @param f1 -- Field 1
 * @param f2 -- Field 2
 * @param error -- Allowable absolute value of error.
 */
template<typename FieldT>
bool field_equal_abs(const FieldT& f1, const FieldT& f2, double error=0.0)
{
  MemoryWindow w1(f1.window_with_ghost());
  MemoryWindow w2(f2.window_with_ghost());
  BoundaryCellInfo bc1(f1.boundary_info());
  BoundaryCellInfo bc2(f2.boundary_info());
  GhostData gd1(f1.get_valid_ghost_data());
  GhostData gd2(f2.get_valid_ghost_data());

  if(w1 != w2) {
    throw( std::runtime_error( "Attempted comparison between fields of unequal size." ) );
  }

  error = std::abs(error);
  bool exact_comparison = error == 0.0;
  FieldT *temp1 = 0;
  FieldT *temp2 = 0;
  typename FieldT::const_iterator *if1 = 0;
  typename FieldT::const_iterator *iend = 0;
  typename FieldT::const_iterator *if2 = 0;

  //initialize if1 and iend
  if(!(if1 = FieldComparisonHelper<FieldT>::init_iterator(f1, w1, bc1, gd1, &temp1))) {
    std::ostringstream msg;
    msg << "Attempted comparison operation on unsupported memory type, at "
      << __FILE__
      << " : " << __LINE__ << std::endl;
    msg << "\t - "
      << "Attempted to compare with a field in "
      << SpatialOps::DeviceTypeTools::get_memory_type_description(
          f1.device_index()) << std::endl;
    throw(std::runtime_error(msg.str()));
  }
  if(temp1)
    iend = new typename FieldT::const_iterator(temp1->end());
  else
    iend = new typename FieldT::const_iterator(f1.end());

  //initialize if2
  if(!(if2 = FieldComparisonHelper<FieldT>::init_iterator(f2, w2, bc2, gd2, &temp2))) {
    std::ostringstream msg;
    msg << "Attempted comparison operation on unsupported memory type, at "
      << __FILE__
      << " : " << __LINE__ << std::endl;
    msg << "\t - "
      << "Attempted to compare with a field in "
      << SpatialOps::DeviceTypeTools::get_memory_type_description(
          f2.device_index()) << std::endl;
    throw(std::runtime_error(msg.str()));
  } //memory leak if runtime error thrown

  //do comparison
  bool result = true;
  for (; *if1 != *iend; ++(*if1), ++(*if2)) {
    if(exact_comparison) {
      if(**if1 != **if2) {
        result = false;
        break;
      }
    }
    else {
      if(std::abs(**if1 - **if2) > error) {
        result = false;
        break;
      }
    }
  }

  if(temp1) delete temp1;
  if(temp2) delete temp2;
  delete if1;
  delete iend;
  delete if2;

  return result;
}
//------------------------------------------------------------------

/**
 * @brief Returns if f1 is element-wise not equal to f2 within a certain number
 * of ulps.
 *
 * This function simply returns the negated result of field_equal_ulp
 *
 * WARNING: Slow in general and comparison with external fields will incur copy penalties.
 *
 * @tparam FieldT -- Any type of SpatialField
 * @param f1 -- Field 1
 * @param f2 -- Field 2
 * @param ulps -- Allowable difference in ulps
 */
template<typename FieldT>
bool field_not_equal_ulp(const FieldT& f1, const FieldT& f2, const unsigned int ulps) {
  return !field_equal_ulp(f1, f2, ulps);
}
//------------------------------------------------------------------

/**
 * @brief Returns if f1 is element-wise equal to f2 within a certain number
 * of ulps.
 *
 * This function determines the amount of ulps two floating point numbers are
 * off and compares them to the allowed tolerance.  Ulp stands for Unit in the
 * Last Place and is a measure of rounding error in floating point numbers.  A
 * more detailed article can be found at:
 * http://en.wikipedia.org/wiki/Unit_in_the_last_place
 *
 * WARNING: Slow in general and comparison with external fields will incur copy penalties.
 *
 * @tparam FieldT -- Any type of SpatialField
 * @param f1 -- Field 1
 * @param f2 -- Field 2
 * @param ulps -- Allowable difference in ulps
 */
template<typename FieldT>
bool field_equal_ulp(const FieldT& f1, const FieldT& f2, const unsigned int ulps)
{
  MemoryWindow w1(f1.window_with_ghost());
  MemoryWindow w2(f2.window_with_ghost());
  BoundaryCellInfo bc1(f1.boundary_info());
  BoundaryCellInfo bc2(f2.boundary_info());
  GhostData gd1(f1.get_valid_ghost_data());
  GhostData gd2(f2.get_valid_ghost_data());

  if(w1 != w2) {
    throw( std::runtime_error( "Attempted comparison between fields of unequal size." ) );
  }

  bool exact_comparison = ulps == 0;
  FieldT *temp1 = 0;
  FieldT *temp2 = 0;
  typename FieldT::const_iterator *if1 = 0;
  typename FieldT::const_iterator *iend = 0;
  typename FieldT::const_iterator *if2 = 0;

  //initialize if1 and iend
  if(!(if1 = FieldComparisonHelper<FieldT>::init_iterator(f1, w1, bc1, gd1, &temp1))) {
    std::ostringstream msg;
    msg << "Attempted comparison operation on unsupported memory type, at "
      << __FILE__
      << " : " << __LINE__ << std::endl;
    msg << "\t - "
      << "Attempted to compare with a field in "
      << SpatialOps::DeviceTypeTools::get_memory_type_description(
          f1.device_index()) << std::endl;
    throw(std::runtime_error(msg.str()));
  }
  if(temp1)
    iend = new typename FieldT::const_iterator(temp1->end());
  else
    iend = new typename FieldT::const_iterator(f1.end());

  //initialize if2
  if(!(if2 = FieldComparisonHelper<FieldT>::init_iterator(f2, w2, bc2, gd2, &temp2))) {
    std::ostringstream msg;
    msg << "Attempted comparison operation on unsupported memory type, at "
      << __FILE__
      << " : " << __LINE__ << std::endl;
    msg << "\t - "
      << "Attempted to compare with a field in "
      << SpatialOps::DeviceTypeTools::get_memory_type_description(
          f2.device_index()) << std::endl;
    throw(std::runtime_error(msg.str()));
  } //memory leak if runtime error thrown

  //do comparison
  bool result = true;
  for (; *if1 != *iend; ++(*if1), ++(*if2)) {
    if(exact_comparison) {
      if (boost::math::float_distance(**if1, **if2) != 0) {
        result = false;
        break;
      }
    }
    else {
      if (std::abs(boost::math::float_distance(**if1, **if2)) > ulps) {
        result =  false;
        break;
      }
    }
  }

  if(temp1) delete temp1;
  if(temp2) delete temp2;
  delete if1;
  delete iend;
  delete if2;

  return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////                                          /////////////////////////////
/////////////////////////////          SCALAR IMPLEMENTATION           /////////////////////////////
/////////////////////////////                                          /////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Returns if f1 is element-wise not equal to the scalar value d
 * within a certain relative tolerance.
 *
 * This function simply calls field_equal and negates it.
 * \c return !field_equal(d, f1, error, error_abs);
 * error_abs is defined as default to be the L2 norm of \c f1 multiplied by \c error
 * and \c FIELDCOMPARISONS_ABS_ERROR_CONST.
 *
 * WARNING: Undefined behavior if f1 is a field of all 0's.
 *
 * WARNING: Slow in general and comparison with external fields will incur copy penalties.
 *
 * @tparam FieldT -- Any type of SpatialField
 * @param d -- Scalar value
 * @param f1 -- Field 1
 * @param error -- Allowable percentage of error. i.e. 1% = .01
 */
template<typename FieldT>
bool field_not_equal(const double d, const FieldT& f1, double error=0.0) {
  double error_abs = error ? nebo_norm(f1)*error*FIELDCOMPARISONS_ABS_ERROR_CONST : 0;
  return !field_equal(d, f1, error, error_abs);
}

/**
 * @brief Returns if f1 is element-wise not equal to the scalar value d
 * within a certain relative tolerance.
 *
 * This function simply calls field_equal and negates it.
 * \c return !field_equal(d, f1, error, error_abs);
 *
 * WARNING: Slow in general and comparison with external fields will incur copy penalties.
 *
 * @tparam FieldT -- Any type of SpatialField
 * @param d -- Scalar value
 * @param f1 -- Field 1
 * @param error -- Allowable percentage of error. i.e. 1% = .01
 * @param error_abs -- Allowable absolute error passed on to \c field_equal
 */
template<typename FieldT>
bool field_not_equal(const double d, const FieldT& f1, double error, const double error_abs) {
  return !field_equal(d, f1, error, error_abs);
}
//------------------------------------------------------------------

/**
 * @brief Returns if f1 is element-wise equal to the scalar value d
 * within a certain relative tolerance.
 *
 * This function returns the result of |d - f1|/(error_abs + |d|) > error element wise.
 * error_abs is defined as default to be the L2 norm of \c f1 multiplied by \c error
 * and \c FIELDCOMPARISONS_ABS_ERROR_CONST.
 *
 * WARNING: Undefined behavior if f1 is a field of all 0's.
 *
 * WARNING: Slow in general and comparison with external fields will incur copy penalties.
 *
 * @tparam FieldT -- Any type of SpatialField
 * @param d -- Scalar value
 * @param f1 -- Field 1
 * @param error -- Allowable percentage of error. i.e. 1% = .01
 */
template<typename FieldT>
bool field_equal(const double d, const FieldT& f1, double error=0.0)
{
  double error_abs = error ? nebo_norm(f1)*error*FIELDCOMPARISONS_ABS_ERROR_CONST : 0;
  return field_equal(d, f1, error, error_abs);
}

/**
 * @brief Returns if f1 is element-wise equal to the scalar value d
 * within a certain relative tolerance.
 *
 * This function returns the result of |d - f1|/(error_abs + |d|) > error element wise.
 *
 * WARNING: Slow in general and comparison with external fields will incur copy penalties.
 *
 * @tparam FieldT -- Any type of SpatialField
 * @param d -- Scalar value
 * @param f1 -- Field 1
 * @param error -- Allowable percentage of error. i.e. 1% = .01
 * @param error_abs -- Allowable absolute error.  This term becomes significant
 * in the calculation as d approaches zero.
 */
template<typename FieldT>
bool field_equal(const double d, const FieldT& f1, double error, const double error_abs)
{
  MemoryWindow w1(f1.window_with_ghost());
  BoundaryCellInfo bc1(f1.boundary_info());
  GhostData gd1(f1.get_valid_ghost_data());

  error = std::abs(error);
  bool exact_comparison = error == 0.0;
  FieldT *temp1 = 0;
  typename FieldT::const_iterator *if1 = 0;
  typename FieldT::const_iterator *iend = 0;

  //initialize if1 and iend
  if(!(if1 = FieldComparisonHelper<FieldT>::init_iterator(f1, w1, bc1, gd1, &temp1))) {
    std::ostringstream msg;
    msg << "Attempted comparison operation on unsupported memory type, at "
      << __FILE__
      << " : " << __LINE__ << std::endl;
    msg << "\t - "
      << "Attempted to compare with a field in "
      << SpatialOps::DeviceTypeTools::get_memory_type_description(
          f1.device_index()) << std::endl;
    throw(std::runtime_error(msg.str()));
  }
  if(temp1)
    iend = new typename FieldT::const_iterator(temp1->end());
  else
    iend = new typename FieldT::const_iterator(f1.end());

  //do comparison
  bool result = true;
  const double denom = std::abs(d) + error_abs;
  for (; *if1 != *iend; ++(*if1)) {
    if(exact_comparison) {
      if (**if1 != d) {
        result = false;
        break;
      }
    }
    else {

      if (std::abs(d - **if1)/denom > error) {
        result =  false;
        break;
      }
    }
  }

  if(temp1) delete temp1;
  delete if1;
  delete iend;

  return result;
}
//------------------------------------------------------------------

/**
 * @brief Returns if f1 is element-wise not equal to Scalar value d within a
 * certain absolute tolerance.
 *
 * This function simply returns the negated result of field_equal_abs
 *
 * WARNING: Slow in general and comparison with external fields will incur copy penalties.
 *
 * @tparam FieldT -- Any type of SpatialField
 * @param d -- Scalar value
 * @param f1 -- Field 1
 * @param error -- Allowable absolute value of error.
 */
template<typename FieldT>
bool field_not_equal_abs(const double d, const FieldT& f1, double error=0.0) {
  return !field_equal_abs(d, f1, error);
}
//------------------------------------------------------------------

/**
 * @brief Returns if f1 is element-wise equal to Scalar value d within a
 * certain absolute tolerance.
 *
 * WARNING: Slow in general and comparison with external fields will incur copy penalties.
 *
 * @tparam FieldT -- Any type of SpatialField
 * @param d -- Scalar value
 * @param f1 -- Field 1
 * @param error -- Allowable absolute value of error.
 */
template<typename FieldT>
bool field_equal_abs(const double d, const FieldT& f1, double error=0.0)
{
  MemoryWindow w1(f1.window_with_ghost());
  BoundaryCellInfo bc1(f1.boundary_info());
  GhostData gd1(f1.get_valid_ghost_data());

  error = std::abs(error);
  bool exact_comparison = error == 0.0;
  FieldT *temp1 = 0;
  typename FieldT::const_iterator *if1 = 0;
  typename FieldT::const_iterator *iend = 0;

  //initialize if1 and iend
  if(!(if1 = FieldComparisonHelper<FieldT>::init_iterator(f1, w1, bc1, gd1, &temp1))) {
    std::ostringstream msg;
    msg << "Attempted comparison operation on unsupported memory type, at "
      << __FILE__
      << " : " << __LINE__ << std::endl;
    msg << "\t - "
      << "Attempted to compare with a field in "
      << SpatialOps::DeviceTypeTools::get_memory_type_description(
          f1.device_index()) << std::endl;
    throw(std::runtime_error(msg.str()));
  }
  if(temp1)
    iend = new typename FieldT::const_iterator(temp1->end());
  else
    iend = new typename FieldT::const_iterator(f1.end());

  //do comparison
  bool result = true;
  for (; *if1 != *iend; ++(*if1)) {
    if(exact_comparison) {
      if(**if1 != d) {
        result = false;
        break;
      }
    }
    else {
      if(std::abs(d - **if1) > error) {
        result = false;
        break;
      }
    }
  }

  if(temp1) delete temp1;
  delete if1;
  delete iend;

  return result;
}
//------------------------------------------------------------------

/**
 * @brief Returns if f1 is element-wise not equal to Scalar value d within a
 * certain number of ulps.
 *
 * This function simply returns the negated result of field_equal_ulp
 *
 * WARNING: Slow in general and comparison with external fields will incur copy penalties.
 *
 * @tparam FieldT -- Any type of SpatialField
 * @param d -- Scalar value
 * @param f1 -- Field 1
 * @param ulps -- Allowable difference in ulps
 */
template<typename FieldT>
bool field_not_equal_ulp(const double d, const FieldT& f1, const unsigned int ulps) {
  return !field_equal_ulp(d, f1, ulps);
}
//------------------------------------------------------------------

/**
 * @brief Returns if f1 is element-wise equal to Scalar value d within a
 * certain number of ulps.
 *
 * This function determines the amount of ulps two floating point numbers are
 * off and compares them to the allowed tolerance.  Ulp stands for Unit in the
 * Last Place and is a measure of rounding error in floating point numbers.  A
 * more detailed article can be found at:
 * http://en.wikipedia.org/wiki/Unit_in_the_last_place
 *
 * WARNING: Slow in general and comparison with external fields will incur copy penalties.
 *
 * @tparam FieldT -- Any type of SpatialField
 * @param d -- Scalar value
 * @param f1 -- Field 1
 * @param ulps -- Allowable difference in ulps
 */
template<typename FieldT>
bool field_equal_ulp(const double d, const FieldT& f1, const unsigned int ulps)
{
  MemoryWindow w1(f1.window_with_ghost());
  BoundaryCellInfo bc1(f1.boundary_info());
  GhostData gd1(f1.get_valid_ghost_data());

  bool exact_comparison = ulps == 0;
  FieldT *temp1 = 0;
  typename FieldT::const_iterator *if1 = 0;
  typename FieldT::const_iterator *iend = 0;

  //initialize if1 and iend
  if(!(if1 = FieldComparisonHelper<FieldT>::init_iterator(f1, w1, bc1, gd1, &temp1))) {
    std::ostringstream msg;
    msg << "Attempted comparison operation on unsupported memory type, at "
      << __FILE__
      << " : " << __LINE__ << std::endl;
    msg << "\t - "
      << "Attempted to compare with a field in "
      << SpatialOps::DeviceTypeTools::get_memory_type_description(
          f1.device_index()) << std::endl;
    throw(std::runtime_error(msg.str()));
  }
  if(temp1)
    iend = new typename FieldT::const_iterator(temp1->end());
  else
    iend = new typename FieldT::const_iterator(f1.end());

  //do comparison
  bool result = true;
  for (; *if1 != *iend; ++(*if1)) {
    if(exact_comparison) {
      if (boost::math::float_distance(d, **if1) != 0) {
        result = false;
        break;
      }
    }
    else {
      if (std::abs(boost::math::float_distance(d, **if1)) > ulps) {
        result =  false;
        break;
      }
    }
  }

  if(temp1) delete temp1;
  delete if1;
  delete iend;

  return result;
}

/**
 *  @class FieldComparisonHelper
 *  @ingroup structured
 *
 *  @brief static helper class for field comparison functions
 *
 *  This class contains private helper functions for the field_equal functions.
 *
 *  @tparam FieldT -- Any type of SpatialField
 *
 */
template<typename FieldT>
class FieldComparisonHelper
{
  /**
   * @brief Returns a field iterator to data in \c field or null.  This function
   * is only intended to be used by the functions in this file.
   *
   * This function will return a \c const_iterator to the data in \c field or a copy
   * of the data in \c field.  This function will copy the data of \c field if \c
   * field is not on local ram, and assign the pointer to the new memory created
   * for the copy in \c localPtr.  It is up for the caller of this function to
   * delete the iterator returned and the \c *localPtr if they are not null.
   *
   * If this function returns null, it means \c field had an unsupported memory
   * type.
   *
   * @tparam FieldT -- Any type of SpatialField
   * @param field -- Field to return an iterator to.
   * @param window -- Memory window of \c field
   * @param bcinfo -- Boundary cell info of \c field
   * @param ghost -- Ghost data of \c field
   * @param localPtr -- Pointer to a pointer which will be overriden if \c field
   * does not exist on local ram.  This pointer will point to memory created for
   * the copy of \c field to local ram.
   */
  private:
    inline static typename FieldT::const_iterator* init_iterator(FieldT const & field,
        MemoryWindow const & window,
        BoundaryCellInfo const & bcinfo,
        GhostData const & ghost,
        FieldT** const localPtr)
    {
      //we can get the iterator directly if on local ram
      if(field.device_index() == CPU_INDEX) {
        return new typename FieldT::const_iterator(field.begin());
      }
      //we need to transfer the field to local ram to iterate
      else if (IS_GPU_INDEX(field.device_index())) {
        *localPtr = new FieldT(window, bcinfo, ghost, NULL, InternalStorage);
        **localPtr = field;

        return new typename FieldT::const_iterator((*localPtr)->begin());
      }
      else {
        //unsupported memory comparison
        return 0;
      }
    }

    friend bool field_equal<FieldT>(const FieldT& f1, const FieldT& f2, double error, const double error_abs);
    friend bool field_equal_abs<FieldT>(const FieldT& f1, const FieldT& f2, double error);
    friend bool field_equal_ulp<FieldT>(const FieldT& f1, const FieldT& f2, const unsigned int ulps);

    friend bool field_equal<FieldT>(const double d, const FieldT& f1, double error, const double error_abs);
    friend bool field_equal_abs<FieldT>(const double d, const FieldT& f1, double error);
    friend bool field_equal_ulp<FieldT>(const double d, const FieldT& f1, const unsigned int ulps);
};

} // namespace structured
} // namespace SpatialOps

#endif //SpatialOps_FieldComparisons_h
