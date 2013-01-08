#include <iostream>
#include <vector>

//--- SpatialOps includes ---//
#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/FieldExpressions.h>
#include <spatialops/FieldReductions.h>

//-- boost includes ---//
#include <boost/program_options.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <test/FieldHelper.h>

namespace po = boost::program_options;

using namespace SpatialOps;

#define RUN_TESTS(TEST,							\
		  TYPE)							\
    ii = 0;								\
    start = boost::posix_time::microsec_clock::universal_time();        \
                                                                        \
    for(; ii < number_of_runs; ii++) {					\
        TEST;								\
    };									\
                                                                        \
    end = boost::posix_time::microsec_clock::universal_time();		\
    std::cout << TYPE;							\
    std::cout << " runs: ";						\
    std::cout << number_of_runs;                                        \
    std::cout << " result: ";						\
    std::cout << (end - start).total_microseconds()*1e-6;               \
    std::cout << std::endl;

#define RUN_GPU_TESTS(TEST,                                             \
                      TYPE)                                             \
    ii = 0;								\
                                                                        \
    copy_to_GPU();                                                      \
                                                                        \
    start = boost::posix_time::microsec_clock::universal_time();        \
                                                                        \
    for(; ii < number_of_runs; ii++) {					\
        TEST;								\
    };									\
                                                                        \
    end = boost::posix_time::microsec_clock::universal_time();		\
                                                                        \
    copy_from_GPU();                                                    \
                                                                        \
    std::cout << TYPE;							\
    std::cout << " runs: ";						\
    std::cout << number_of_runs;                                        \
    std::cout << " result: ";						\
    std::cout << (end - start).total_microseconds()*1e-6;               \
    std::cout << std::endl;

#define RUN_FULL_GPU_TESTS(TEST,                                        \
                           TYPE)                                        \
    ii = 0;								\
    start = boost::posix_time::microsec_clock::universal_time();        \
                                                                        \
    for(; ii < number_of_runs; ii++) {					\
        copy_to_GPU();                                                  \
        TEST;								\
        copy_from_GPU();                                                \
    };									\
                                                                        \
    end = boost::posix_time::microsec_clock::universal_time();		\
                                                                        \
    std::cout << TYPE;							\
    std::cout << " runs: ";						\
    std::cout << number_of_runs;                                        \
    std::cout << " result: ";						\
    std::cout << (end - start).total_microseconds()*1e-6;               \
    std::cout << std::endl;

int main( int iarg, char* carg[] )
{
    typedef SpatialOps::structured::SVolField Field;

    std::vector<int> npts(3,1);
    int number_of_runs;
#ifdef FIELD_EXPRESSION_THREADS
    int thread_count;
#endif

    // parse the command line options input describing the problem
    {
        po::options_description desc("Supported Options");
        desc.add_options()
            ( "help", "print help message" )
            ( "nx", po::value<int>(&npts[0])->default_value(10), "Grid in x" )
            ( "ny", po::value<int>(&npts[1])->default_value(10), "Grid in y" )
            ( "nz", po::value<int>(&npts[2])->default_value(10), "Grid in z" )
#ifdef FIELD_EXPRESSION_THREADS
            ( "tc", po::value<int>(&thread_count)->default_value(NTHREADS), "Number of threads for Nebo")
#endif
            ( "runs", po::value<int>(&number_of_runs)->default_value(1), "Number of iterations of each test");

        po::variables_map args;
        po::store( po::parse_command_line(iarg,carg,desc), args );
        po::notify(args);

        if (args.count("help")) {
            std::cout << desc << "\n";
            return 1;
        }

#ifdef FIELD_EXPRESSION_THREADS
        set_hard_thread_count(thread_count);
#endif
    }

    const SpatialOps::structured::MemoryWindow window( SpatialOps::structured::get_window_with_ghost<Field>(npts,true,true,true) );
    int rawSize = ((npts[0] > ? npts[0] + 2 : 1) *
                   (npts[1] > ? npts[1] + 2 : 1) *
                   (npts[2] > ? npts[2] + 2 : 1));

    std::vector<Field> fs(106, Field(window, NULL));

    for(int x = 0; x < 36; x++)
        initialize_field(fs[x], x * rawSize, false);

    boost::posix_time::ptime start( boost::posix_time::microsec_clock::universal_time() );
    boost::posix_time::ptime end( boost::posix_time::microsec_clock::universal_time() );
    int ii;

    //this is to warm up the system:
    result <<= ((((sin(fs[0]) + sin(fs[1]) + sin(fs[2])) + (sin(fs[3]) + sin(fs[4]) + sin(fs[5])) + (sin(fs[6]) + sin(fs[7]) + sin(fs[8]))) +
                 ((sin(fs[9]) + sin(fs[10]) + sin(fs[11])) + (sin(fs[12]) + sin(fs[13]) + sin(fs[14])) + (sin(fs[15]) + sin(fs[16]) + sin(fs[17])))) +
                (((sin(fs[18]) + sin(fs[19]) + sin(fs[20])) + (sin(fs[21]) + sin(fs[22]) + sin(fs[23])) + (sin(fs[24]) + sin(fs[25]) + sin(fs[26]))) +
                 ((sin(fs[27]) + sin(fs[28]) + sin(fs[29])) + (sin(fs[30]) + sin(fs[31]) + sin(fs[32])) + (sin(fs[33]) + sin(fs[34]) + sin(fs[35])))));

    //1 Loop
    RUN_TESTS(result <<= ((((sin(fs[0]) + sin(fs[1]) + sin(fs[2])) + (sin(fs[3]) + sin(fs[4]) + sin(fs[5])) + (sin(fs[6]) + sin(fs[7]) + sin(fs[8]))) +
                           ((sin(fs[9]) + sin(fs[10]) + sin(fs[11])) + (sin(fs[12]) + sin(fs[13]) + sin(fs[14])) + (sin(fs[15]) + sin(fs[16]) + sin(fs[17])))) +
                          (((sin(fs[18]) + sin(fs[19]) + sin(fs[20])) + (sin(fs[21]) + sin(fs[22]) + sin(fs[23])) + (sin(fs[24]) + sin(fs[25]) + sin(fs[26]))) +
                           ((sin(fs[27]) + sin(fs[28]) + sin(fs[29])) + (sin(fs[30]) + sin(fs[31]) + sin(fs[32])) + (sin(fs[33]) + sin(fs[34]) + sin(fs[35]))))),
              "1-loop");

    //3 Loop
    RUN_TESTS(fs[36] <<= (((sin(fs[0]) + sin(fs[1]) + sin(fs[2])) + (sin(fs[3]) + sin(fs[4]) + sin(fs[5])) + (sin(fs[6]) + sin(fs[7]) + sin(fs[8]))) +
                          ((sin(fs[9]) + sin(fs[10]) + sin(fs[11])) + (sin(fs[12]) + sin(fs[13]) + sin(fs[14])) + (sin(fs[15]) + sin(fs[16]) + sin(fs[17]))));
              fs[37] <<= (((sin(fs[18]) + sin(fs[19]) + sin(fs[20])) + (sin(fs[21]) + sin(fs[22]) + sin(fs[23])) + (sin(fs[24]) + sin(fs[25]) + sin(fs[26]))) +
                          ((sin(fs[27]) + sin(fs[28]) + sin(fs[29])) + (sin(fs[30]) + sin(fs[31]) + sin(fs[32])) + (sin(fs[33]) + sin(fs[34]) + sin(fs[35]))));
              result <<= fs[36] + fs[37],
              "3-loop");
    //5 Loop
    RUN_TESTS(fs[36] <<= ((sin(fs[0]) + sin(fs[1]) + sin(fs[2])) + (sin(fs[3]) + sin(fs[4]) + sin(fs[5])) + (sin(fs[6]) + sin(fs[7]) + sin(fs[8])));
              fs[37] <<= ((sin(fs[9]) + sin(fs[10]) + sin(fs[11])) + (sin(fs[12]) + sin(fs[13]) + sin(fs[14])) + (sin(fs[15]) + sin(fs[16]) + sin(fs[17])));
              fs[38] <<= ((sin(fs[18]) + sin(fs[19]) + sin(fs[20])) + (sin(fs[21]) + sin(fs[22]) + sin(fs[23])) + (sin(fs[24]) + sin(fs[25]) + sin(fs[26])));
              fs[39] <<= ((sin(fs[27]) + sin(fs[28]) + sin(fs[29])) + (sin(fs[30]) + sin(fs[31]) + sin(fs[32])) + (sin(fs[33]) + sin(fs[34]) + sin(fs[35])));
              result <<= ((fs[36] + fs[37]) +
                          (fs[38] + fs[39])),
              "5-loop");

    //13 Loop
    RUN_TESTS(fs[36] <<= (sin(fs[0]) + sin(fs[1]) + sin(fs[2]));
              fs[37] <<= (sin(fs[3]) + sin(fs[4]) + sin(fs[5]));
              fs[38] <<= (sin(fs[6]) + sin(fs[7]) + sin(fs[8]));
              fs[39] <<= (sin(fs[9]) + sin(fs[10]) + sin(fs[11]));
              fs[40] <<= (sin(fs[12]) + sin(fs[13]) + sin(fs[14]));
              fs[41] <<= (sin(fs[15]) + sin(fs[16]) + sin(fs[17]));
              fs[42] <<= (sin(fs[18]) + sin(fs[19]) + sin(fs[20]));
              fs[43] <<= (sin(fs[21]) + sin(fs[22]) + sin(fs[23]));
              fs[44] <<= (sin(fs[24]) + sin(fs[25]) + sin(fs[26]));
              fs[45] <<= (sin(fs[27]) + sin(fs[28]) + sin(fs[29]));
              fs[46] <<= (sin(fs[30]) + sin(fs[31]) + sin(fs[32]));
              fs[47] <<= (sin(fs[33]) + sin(fs[34]) + sin(fs[35]));
              result <<= (((fs[36] + fs[37] + fs[38]) +
                           (fs[39] + fs[40] + fs[41])) +
                          ((fs[42] + fs[43] + fs[44]) +
                           (fs[45] + fs[46] + fs[47]))),
              "13-loop");

    //35 Loop
    RUN_TESTS(fs[36] <<= (sin(fs[0]) + sin(fs[1]));
              fs[37] <<= (fs[36] + sin(fs[2]));
              fs[38] <<= (sin(fs[3]) + sin(fs[4]));
              fs[39] <<= (fs[38] + sin(fs[5]));
              fs[40] <<= (sin(fs[6]) + sin(fs[7]));
              fs[41] <<= (fs[40] + sin(fs[8]));
              fs[42] <<= (sin(fs[9]) + sin(fs[10]));
              fs[43] <<= (fs[42] + sin(fs[11]));
              fs[44] <<= (sin(fs[12]) + sin(fs[13]));
              fs[45] <<= (fs[44] + sin(fs[14]));
              fs[46] <<= (sin(fs[15]) + sin(fs[16]));
              fs[47] <<= (fs[46] + sin(fs[17]));
              fs[48] <<= (sin(fs[18]) + sin(fs[19]));
              fs[49] <<= (fs[48] + sin(fs[20]));
              fs[50] <<= (sin(fs[21]) + sin(fs[22]));
              fs[51] <<= (fs[50] + sin(fs[23]));
              fs[52] <<= (sin(fs[24]) + sin(fs[25]));
              fs[53] <<= (fs[52] + sin(fs[26]));
              fs[54] <<= (sin(fs[27]) + sin(fs[28]));
              fs[55] <<= (fs[54] + sin(fs[29]));
              fs[56] <<= (sin(fs[30]) + sin(fs[31]));
              fs[57] <<= (fs[56] + sin(fs[32]));
              fs[58] <<= (sin(fs[33]) + sin(fs[34]));
              fs[59] <<= (fs[58] + sin(fs[35]));
              fs[60] <<= (fs[37] + fs[39]);
              fs[61] <<= (fs[60] + fs[41]);
              fs[62] <<= (fs[43] + fs[45]);
              fs[63] <<= (fs[62] + fs[47]);
              fs[64] <<= (fs[49] + fs[51]);
              fs[65] <<= (fs[64] + fs[53]);
              fs[66] <<= (fs[55] + fs[57]);
              fs[67] <<= (fs[66] + fs[59]);
              fs[68] <<= (fs[61] + fs[63]);
              fs[69] <<= (fs[65] + fs[67]);
              result <<= (fs[68] + fs[69]),
              "35-loop");

    //71 Loop
    RUN_TESTS(fs[36] <<= sin(fs[0]);
              fs[37] <<= sin(fs[1]);
              fs[38] <<= (fs[36] + fs[37]);
              fs[39] <<= sin(fs[2]);
              fs[40] <<= (fs[38] + fs[39]);
              fs[41] <<= sin(fs[3]);
              fs[42] <<= sin(fs[4]);
              fs[43] <<= (fs[41] + fs[42]);
              fs[44] <<= sin(fs[5]);
              fs[45] <<= (fs[43] + fs[44]);
              fs[46] <<= (fs[40] + fs[45]);
              fs[47] <<= sin(fs[6]);
              fs[48] <<= sin(fs[7]);
              fs[49] <<= (fs[47] + fs[48]);
              fs[50] <<= sin(fs[8]);
              fs[51] <<= (fs[49] + fs[50]);
              fs[52] <<= (fs[46] + fs[51]);
              fs[53] <<= sin(fs[9]);
              fs[54] <<= sin(fs[10]);
              fs[55] <<= (fs[53] + fs[54]);
              fs[56] <<= sin(fs[11]);
              fs[57] <<= (fs[55] + fs[56]);
              fs[58] <<= sin(fs[12]);
              fs[59] <<= sin(fs[13]);
              fs[60] <<= (fs[58] + fs[59]);
              fs[61] <<= sin(fs[14]);
              fs[62] <<= (fs[60] + fs[61]);
              fs[63] <<= (fs[57] + fs[62]);
              fs[64] <<= sin(fs[15]);
              fs[65] <<= sin(fs[16]);
              fs[66] <<= (fs[64] + fs[65]);
              fs[67] <<= sin(fs[17]);
              fs[68] <<= (fs[66] + fs[67]);
              fs[69] <<= (fs[63] + fs[68]);
              fs[70] <<= (fs[52] + fs[69]); //****
              fs[71] <<= sin(fs[18]);
              fs[72] <<= sin(fs[19]);
              fs[73] <<= (fs[71] + fs[72]);
              fs[74] <<= sin(fs[20]);
              fs[75] <<= (fs[73] + fs[74]);
              fs[76] <<= sin(fs[21]);
              fs[77] <<= sin(fs[22]);
              fs[78] <<= (fs[76] + fs[77]);
              fs[79] <<= sin(fs[23]);
              fs[80] <<= (fs[78] + fs[79]);
              fs[81] <<= (fs[75] + fs[80]);
              fs[82] <<= sin(fs[24]);
              fs[83] <<= sin(fs[25]);
              fs[84] <<= (fs[82] + fs[83]);
              fs[85] <<= sin(fs[26]);
              fs[86] <<= (fs[84] + fs[85]);
              fs[87] <<= (fs[81] + fs[86]);
              fs[88] <<= sin(fs[27]);
              fs[89] <<= sin(fs[28]);
              fs[90] <<= (fs[88] + fs[89]);
              fs[91] <<= sin(fs[29]);
              fs[92] <<= (fs[90] + fs[91]);
              fs[93] <<= sin(fs[30]);
              fs[94] <<= sin(fs[31]);
              fs[95] <<= (fs[93] + fs[94]);
              fs[96] <<= sin(fs[32]);
              fs[97] <<= (fs[95] + fs[96]);
              fs[98] <<= (fs[92] + fs[97]);
              fs[99] <<= sin(fs[33]);
              fs[100] <<= sin(fs[34]);
              fs[101] <<= (fs[99] + fs[100]);
              fs[102] <<= sin(fs[35]);
              fs[103] <<= (fs[101] + fs[102]);
              fs[104] <<= (fs[98] + fs[103]);
              fs[105] <<= (fs[87] + fs[104]); //****
              result <<= (fs[70] + fs[105]),
              "71-loop");

    return 0;
}
