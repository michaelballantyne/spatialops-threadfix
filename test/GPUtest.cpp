#include <iostream>
#include <vector>

//--- SpatialOps includes ---//
#include <spatialops/SpatialOpsConfigure.h>
#include <spatialops/structured/FVTools.h>
#include <spatialops/structured/FVStaggeredFieldTypes.h>
#include <spatialops/FieldExpressions.h>

#include <spatialops/structured/GPUtest.h>

//-- boost includes ---//
#include <boost/program_options.hpp>

#include <test/FieldHelper.h>

namespace po = boost::program_options;

using namespace SpatialOps;

template<typename FieldType>
void nebo_allocate(typename FieldType::memory_window mw,
                   std::vector<FieldType *> & ptr_fields,
                   std::vector<FieldType> & fields,
                   int numFields) {
    ptr_fields.clear();
    fields.clear();
    for(int ii = 0; ii < numFields; ii++) {
        ptr_fields.push_back(new FieldType(mw, NULL));
        fields.push_back(* (ptr_fields.back()));
    };
};

template<typename FieldType>
void nebo_initialize(std::vector<FieldType> & fields,
                     int start) {
    for(int ii = 0; ii < fields.size(); ii++)
        initialize_field(fields[ii], ii * start, false);
};

template<typename FieldType>
void nebo_cuda_allocate(typename FieldType::memory_window mw,
                        std::vector<FieldType *> & ptr_gpu_fields,
                        std::vector<FieldType> & gpu_fields,
                        int numFields) {
    ptr_gpu_fields.clear();
    gpu_fields.clear();
    for(int ii = 0; ii < numFields; ii++) {
#ifdef ENABLE_CUDA
        ptr_gpu_fields.push_back(new FieldType(mw,
                                               NULL,
                                               structured::InternalStorage,
                                               EXTERNAL_CUDA_GPU,
                                               0));
#else
        ptr_gpu_fields.push_back(new FieldType(mw, NULL));
#endif

        gpu_fields.push_back(* (ptr_gpu_fields.back()));
    };
};

template<typename FieldType>
void nebo_copy(std::vector<FieldType> & dest_fields,
               std::vector<FieldType> & src_fields,
               int numFields) {
    for(int ii = 0; ii < numFields; ii++)
        dest_fields[ii] = src_fields[ii];
};

template<typename FieldType>
void nebo_cuda_free(std::vector<FieldType *> & ptr_gpu_fields,
                    std::vector<FieldType> & gpu_fields) {
    ptr_gpu_fields.clear();
    gpu_fields.clear();
};

int main( int iarg, char* carg[] )
{
    typedef SpatialOps::structured::SVolField Field;

    std::vector<int> npts(3,1);
    int number_of_runs;

    // parse the command line options input describing the problem
    {
        po::options_description desc("Supported Options");
        desc.add_options()
            ( "help", "print help message" )
            ( "nx", po::value<int>(&npts[0])->default_value(10), "Grid in x" )
            ( "ny", po::value<int>(&npts[1])->default_value(10), "Grid in y" )
            ( "nz", po::value<int>(&npts[2])->default_value(10), "Grid in z" )
            ( "runs", po::value<int>(&number_of_runs)->default_value(1), "Number of iterations of each test");

        po::variables_map args;
        po::store( po::parse_command_line(iarg,carg,desc), args );
        po::notify(args);

        if (args.count("help")) {
            std::cout << desc << "\n";
            return 1;
        }
    }

    const SpatialOps::structured::MemoryWindow window( SpatialOps::structured::get_window_with_ghost<Field>(npts,true,true,true) );
    int rawSize = ((npts[0] > 1 ? npts[0] + 2 : 1) *
                   (npts[1] > 1 ? npts[1] + 2 : 1) *
                   (npts[2] > 1 ? npts[2] + 2 : 1));

#define COUNT 4

    //CPU fields:
    std::vector<Field *> ptr_fs;
    std::vector<Field> fs;
    Field result(window, NULL);
    Field result_from_gpu(window, NULL);

    //GPU fields:
    std::vector<Field *> ptr_gpu_fs;
    std::vector<Field> gpu_fs;
#ifdef ENABLE_CUDA
    Field gpu_result(window, NULL, structured::InternalStorage, EXTERNAL_CUDA_GPU, 0);
#else
    Field gpu_result(window, NULL);
#endif
    //allocate space on CPU:
    nebo_allocate(window, ptr_fs, fs, COUNT);

    //initialize CPU fields:
    nebo_initialize(fs, rawSize);

    //allocate space on GPU:
    nebo_cuda_allocate(window, ptr_gpu_fs, gpu_fs, COUNT);

    //copy input fields to GPU:
    nebo_copy(gpu_fs, fs, 2);

    //compute result:
    result <<= fs[0] + sin(fs[1]);
    addsin(gpu_result, gpu_fs[0], gpu_fs[1]);

    //copy GPU result back to CPU:
    result_from_gpu = gpu_result;

    //compare result:
    bool check = display_fields_compare(result,
                                        result_from_gpu,
                                        true, false);

    std::cout << (check ? "Perfect match!" : "Not a perfect match.") << std::endl;

    return 0;
}
