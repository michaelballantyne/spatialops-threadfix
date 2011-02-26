#ifndef SpatialOps_Structured_Stencil_h
#define SpatialOps_Structured_Stencil_h

namespace SpatialOps{
namespace structured{

  class MemoryWindow;
  class IntVec;

  /**
   *  \struct Stencil2Helper
   *  \author James C. Sutherland
   *
   *  \brief Provides methods to adjust iterator positions when
   *         applying 2-point stencils to fields.  See also Stencil2
   */
  template< typename SrcT, typename DestT >
  struct Stencil2Helper
  {
    Stencil2Helper( const MemoryWindow& wsrc, ///< MemoryWindow for the source field
                    const MemoryWindow& wdest ///< MemoryWindow for the destination field
                    );
    IntVec high() const;           ///< upper bounds for mesh loop over dest field
    IntVec low() const;            ///< lower bounds for mesh loop over dest field
    IntVec src_increment() const;  ///< source field increment count after each loop (x,y,z)
    IntVec dest_increment() const; ///< destination field increment count after each loop (x,y,z)
    unsigned int src_offset_lo() const;  ///< offset for the "low" (minus-side) source field iterator
    unsigned int src_offset_hi() const;  ///< offset for the "high" (plus-side) source field iterator
    unsigned int dest_offset() const;    ///< offset for the destination field iterator (nominally 0)
  };

  /**
   *  \class Stencil2
   *  \author James C. Sutherland
   *
   *  \brief Support for implementing simple two-point stencils in
   *         one-dimension on structured meshes.
   *
   *  \tparam OpT - the type of operator
   *  \tparam SrcT - the type of field the operator is applied to
   *  \tparam DestT - the type of field the operator produces
   *
   *  See also Stencil2Helper
   */
  template< typename OperatorT, typename SrcFieldT, typename DestFieldT >
  class Stencil2
  {
    const double coefLo_, coefHi_;
  public:

    typedef OperatorT  OpT;
    typedef SrcFieldT  SrcFieldType;
    typedef DestFieldT DestFieldType;

    /**
     *  \brief construct a stencil with the specified coefficients
     *  \param coefLo the coefficient to multiply the (-) side field by
     *  \param coefHi the coefficient to multiply the (+) side field by
     */
    Stencil2( const double coefLo, const double coefHi );

    ~Stencil2();

    void apply_to_field( const SrcFieldType& src, DestFieldType& dest ) const;

    double get_minus_coef() const{ return coefLo_; }
    double get_plus_coef () const{ return coefHi_; }
  };

} // namespace structured
} // namespace SpatialOps

#endif // SpatialOps_Structured_Stencil_h
