/**
 *  \file   IntVec.h
 *
 *  \date   Sep 28, 2011
 *  \author James C. Sutherland
 */

#ifndef SpatialOps_IntVec_h
#define SpatialOps_IntVec_h

namespace SpatialOps{
  namespace structured{

    /**
     *  \class IntVec
     *  \author James C. Sutherland
     *  \ingroup structured
     *  \brief provides a lightweight class to deal with a 3D vector of integers.
     */
    class IntVec
    {
      friend std::ostream& operator<<(std::ostream&, const IntVec&);

      int ijk[3];

#   ifdef SOPS_BOOST_SERIALIZATION
      friend class boost::serialization::access;
      template<typename Archive>
      void serialize( Archive& ar, const unsigned int version )
      {
        ar & ijk;
      }
#   endif

    public:

      IntVec(){ ijk[0]=0; ijk[1]=0; ijk[2]=0; }

      inline IntVec( const int i, const int j, const int k ){
        ijk[0]=i; ijk[1]=j; ijk[2]=k;
      }

      inline IntVec( const int vec[3] ){
        ijk[0]=vec[0];  ijk[1]=vec[1];  ijk[2]=vec[2];
      }

      IntVec( const std::vector<int>& vec ){
        ijk[0]=vec[0]; ijk[1]=vec[1]; ijk[2]=vec[2];
      }

      inline IntVec( const IntVec& x ){
        ijk[0]=x.ijk[0];  ijk[1]=x.ijk[1];  ijk[2]=x.ijk[2];
      }

      inline int  operator[](const size_t i) const{ return ijk[i]; }
      inline int& operator[](const size_t i)      { return ijk[i]; }

      IntVec& operator=(const IntVec& x)
      {
        for( size_t i=0; i<3; ++i ) ijk[i] = x.ijk[i];
        return *this;
      }

      inline bool operator==(const IntVec& v) const
        {
        return (ijk[0]==v.ijk[0]) & (ijk[1]==v.ijk[1]) & (ijk[2]==v.ijk[2]);
        }
      inline bool operator!=(const IntVec& v) const
        {
        return (ijk[0]!=v.ijk[0]) | (ijk[1]!=v.ijk[1]) | (ijk[2]!=v.ijk[2]);
        }

      inline IntVec operator+( const IntVec& v ) const{
        return IntVec( ijk[0] + v.ijk[0],
            ijk[1] + v.ijk[1],
            ijk[2] + v.ijk[2] );
      }
      inline IntVec operator-( const IntVec& v ) const{
        return IntVec( ijk[0] - v.ijk[0],
            ijk[1] - v.ijk[1],
            ijk[2] - v.ijk[2] );
      }

      inline IntVec& operator+=( const IntVec& v ){
        ijk[0] += v.ijk[0];
        ijk[1] += v.ijk[1];
        ijk[2] += v.ijk[2];
        return *this;
      }
      inline IntVec& operator-=( const IntVec& v ){
        ijk[0] -= v.ijk[0];
        ijk[1] -= v.ijk[1];
        ijk[2] -= v.ijk[2];
        return *this;
      }
    };

  } // namespace structured
} // namespace SpatialOps

#endif /* SpatialOps_IntVec_h */
