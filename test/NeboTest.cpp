#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/structured/FVTools.h>

#include <spatialops/Nebo.h>
#include <test/TestHelper.h>
#include <test/FieldHelper.h>

#include <spatialops/structured/FieldComparisons.h>

#include <iostream>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace SpatialOps;
using namespace structured;
using std::cout;
using std::endl;

#define MANUAL(EXPR)                                            \
    {                                                           \
        const int xLength = window.extent(0);                   \
        const int yLength = window.extent(1);                   \
        const int zLength = window.extent(2);                   \
                                                                \
        for(int kk = 0; kk < zLength; kk++) {                   \
            for(int jj = 0; jj < yLength; jj++) {               \
                for(int ii = 0; ii < xLength; ii++) {           \
                    ref(ii, jj, kk) = EXPR;                     \
                };                                              \
            };                                                  \
        };                                                      \
    }
#define RUNTEST(NEBOEXPR, EXPR, MESSAGE)                        \
    {                                                           \
        test <<= NEBOEXPR;                                      \
        MANUAL(EXPR);                                           \
        status(field_equal(ref, test, 0.0), MESSAGE);            \
    }                                                           \

#define RUNTESTULP(NEBOEXPR, EXPR, MESSAGE, ULPS)               \
    {                                                           \
        test <<= NEBOEXPR;                                      \
        MANUAL(EXPR);                                           \
        status(field_equal_ulp(ref, test, ULPS), MESSAGE);       \
    }                                                           \

#define RUNTESTIGNORENANULP(NEBOEXPR, EXPR, MESSAGE, ULPS)      \
    {                                                           \
        test <<= cond(NEBOEXPR != NEBOEXPR, 0.0)(NEBOEXPR);     \
        MANUAL(EXPR != EXPR ? 0.0 : EXPR);                      \
        status(field_equal_ulp(ref, test, ULPS), MESSAGE);       \
    }                                                           \

#define RUN_BINARY_OP_TEST(OP, TESTTYPE)                                                                            \
    {                                                                                                               \
        RUNTEST(50.3 OP 4.7, 50.3 OP 4.7, TESTTYPE" (Scalar x Scalar) test");                                       \
        RUNTEST(input1 OP 4.7, INPUT1 OP 4.7, TESTTYPE" (Field x Scalar) test");                                    \
        RUNTEST(50.3 OP input2, 50.3 OP INPUT2, TESTTYPE" (Scalar x Field) test");                                  \
        RUNTEST(input1 OP input2, INPUT1 OP INPUT2, TESTTYPE" (Field x Field) test");                               \
        RUNTEST(input1 OP (input2 OP input3), INPUT1 OP (INPUT2 OP INPUT3), TESTTYPE" (Field x SubExpr) test");     \
        RUNTEST((input1 OP input2) OP input3, (INPUT1 OP INPUT2) OP INPUT3, TESTTYPE" (SubExpr x Field) test");     \
        RUNTEST((input1 OP input2) OP (input3 OP input4),                                                           \
                (INPUT1 OP INPUT2) OP (INPUT3 OP INPUT4),                                                           \
                TESTTYPE" (SubExpr x SubExpr) test");                                                               \
        RUNTEST(50.3 OP (input2 OP input3), 50.3 OP (INPUT2 OP INPUT3), TESTTYPE" (Scalar x SubExpr) test");        \
        RUNTEST((input1 OP input2) OP 4.7, (INPUT1 OP INPUT2) OP 4.7, TESTTYPE" (SubExpr x Scalar) test");          \
    }                                                                                                               \

#define INPUT1 input1(ii, jj, kk)
#define INPUT2 input2(ii, jj, kk)
#define INPUT3 input3(ii, jj, kk)
#define INPUT4 input4(ii, jj, kk)

int main( int iarg, char* carg[] )
{
    int nx, ny, nz;
    bool bcplus[] = { false, false, false };

    po::options_description desc("Supported Options");
    desc.add_options()
        ( "help", "print help message\n" )
        ( "nx",   po::value<int>(&nx)->default_value(11), "number of points in x-dir for base mesh" )
        ( "ny",   po::value<int>(&ny)->default_value(11), "number of points in y-dir for base mesh" )
        ( "nz",   po::value<int>(&nz)->default_value(11), "number of points in z-dir for base mesh" )
        ( "bcx",  "physical boundary on +x side?" )
        ( "bcy",  "physical boundary on +y side?" )
        ( "bcz",  "physical boundary on +z side?" );

    po::variables_map args;
    po::store( po::parse_command_line(iarg,carg,desc), args );
    po::notify(args);

    if( args.count("bcx") ) bcplus[0] = true;
    if( args.count("bcy") ) bcplus[1] = true;
    if( args.count("bcz") ) bcplus[2] = true;

    if( args.count("help") ){
      cout << desc << endl
           << "Examples:" << endl
           << " test_nebo --nx 5 --ny 10 --nz 3 --bcx" << endl
           << " test_nebo --bcx --bcy --bcz" << endl
           << " test_nebo --nx 50 --bcz" << endl
           << endl;
      return -1;
    }

    typedef SVolField Field;

    const int nghost = 1;
    const GhostData ghost(nghost);
    const BoundaryCellInfo bcinfo = BoundaryCellInfo::build<Field>(bcplus[0],bcplus[1],bcplus[2]);
    const MemoryWindow window( get_window_with_ghost(IntVec(nx,ny,nz),ghost,bcinfo) );

    Field input1( window, bcinfo, ghost, NULL );
    Field input2( window, bcinfo, ghost, NULL );
    Field input3( window, bcinfo, ghost, NULL );
    Field input4( window, bcinfo, ghost, NULL );
    Field   test( window, bcinfo, ghost, NULL );
    Field    ref( window, bcinfo, ghost, NULL );

    const int total = nx * ny * nz;

    initialize_field(input1, 0.0);
    initialize_field(input2, total);
    initialize_field(input3, 2 * total);
    initialize_field(input4, 3 * total);

    TestHelper status(true);

    RUNTEST(0.0, 0.0, "scalar assignment test"); cout << "\n";

    RUN_BINARY_OP_TEST(+, "summation"); cout << "\n";
    RUN_BINARY_OP_TEST(-, "difference"); cout << "\n";
    RUN_BINARY_OP_TEST(/, "product"); cout << "\n";
    RUN_BINARY_OP_TEST(*, "division"); cout << "\n";

    RUNTESTULP(sin(input1), std::sin(INPUT1), "sin test", 1);
    RUNTESTULP(cos(input1), std::cos(INPUT1), "cos test", 1);
    RUNTESTULP(tan(input1), std::tan(INPUT1), "tan test", 2);
    RUNTESTULP(exp(input1), std::exp(INPUT1), "exp test", 1);
    //documentation says with 1 ulp, empirically found to be 2
    RUNTESTULP(tanh(input1), std::tanh(INPUT1), "tanh test", 2);
    //display_fields_compare(ref, test, true, true);
    RUNTEST(abs(input1), std::abs(INPUT1), "abs test");
    RUNTEST(-input1, -INPUT1, "negation test");
    RUNTESTIGNORENANULP(pow(input1, input2), std::pow(INPUT1, INPUT2), "power test", 2);
    RUNTESTIGNORENANULP(sqrt(input1), std::sqrt(INPUT1), "square root test", 0);
    RUNTESTIGNORENANULP(log(input1), std::log(INPUT1), "log test", 1);

    RUNTEST(cond(input1 == input2, true)(false), INPUT1 == INPUT2, "equivalence test");
    RUNTEST(cond(input1 != input2, true)(false), INPUT1 != INPUT2, "non-equivalence test");

    RUNTEST(cond(input1 < input2, true)(false), INPUT1 < INPUT2, "less than test");
    RUNTEST(cond(input1 <= input2, true)(false), INPUT1 <= INPUT2, "less than or equal test");
    RUNTEST(cond(input1 > input2, true)(false), INPUT1 > INPUT2, "greater than test");
    RUNTEST(cond(input1 >= input2, true)(false), INPUT1 >= INPUT2, "greater than or equal test");

    //these tests are dependent on
    //nebo less than working
    RUNTEST(cond(input1 < input2 && input3 < input4, true)(false), INPUT1 < INPUT2 && INPUT3 < INPUT4, "boolean and test");
    RUNTEST(cond(input1 < input2 || input3 < input4, true)(false), INPUT1 < INPUT2 || INPUT3 < INPUT4, "boolean or test");
    RUNTEST(cond(!(input1 < input2), true)(false), !(INPUT1 < INPUT2), "boolean not test");

    RUNTEST(max(input1, 0.0), std::max(INPUT1, 0.0), "max test");
    RUNTEST(min(input1, 0.0), std::min(INPUT1, 0.0), "min test");

    if( status.ok() ) {
      cout << "ALL TESTS PASSED :)" << endl;
      return 0;
    }
    else {
        cout << "******************************" << endl
             << " At least one test FAILED! :(" << endl
             << "******************************" << endl;
        return -1;
    };
}
