/*
 * Copyright (c) 2011 The University of Utah
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

#include<spatialops/structured/MemoryWindow.h>

using std::cout;
using std::endl;
using SpatialOps::structured::MemoryWindow;

template<typename Field>
inline void internal_initialize_field(typename Field::iterator fi,
                                      MemoryWindow mw,
                                      bool print) {
    int xLength = mw.extent(0);
    int yLength = mw.extent(1);
    int zLength = mw.extent(2);

    for(int z = 1; z <= zLength; z++) {
        for(int y = 1; y <= yLength; y++) {
            for(int x = 1; x <= xLength; x++, fi++) {
                *fi = sin(x + y * xLength + z * xLength * yLength);
                if(print) cout << *fi << " ";
            };
            if(print) cout << endl;
        };
        if(print) cout << endl;
    };
 };

template<typename Field>
inline void initialize_field(Field & f,
                             bool print = false) {
    internal_initialize_field<Field>(f.begin(), f.window_with_ghost(), print);
};

template<typename Field>
inline void interior_initialize_field(Field & f,
                                      bool print = false) {
    internal_initialize_field<Field>(f.interior_begin(), f.window_without_ghost(), print);
};

template<typename Field>
inline void internal_print_field(typename Field::const_iterator fi,
                                 MemoryWindow const & mw) {
    int xLength = mw.extent(0);
    int yLength = mw.extent(1);
    int zLength = mw.extent(2);

    for(int z = 1; z <= zLength; z++) {
        for(int y = 1; y <= yLength; y++) {
            for(int x = 1; x <= xLength; x++, fi++) {
                cout << *fi << " ";
            };
            cout << endl;
        };
        cout << endl;
    };
};

template<typename Field>
inline void print_field(Field const & f) {
    internal_print_field<Field>(f.begin(), f.window_with_ghost());
};

template<typename Field>
inline void interior_print_field(Field const & f) { 
    internal_print_field<Field>(f.interior_begin(), f.window_without_ghost());
};

template<typename Field>
inline bool internal_display_fields_compare(typename Field::const_iterator fi1,
                                            typename Field::const_iterator if2,
                                            MemoryWindow const & mw,
                                            bool display,
                                            bool print) {
    bool result = true;
    int xLength = mw.extent(0);
    int yLength = mw.extent(1);
    int zLength = mw.extent(2);

    //copies of iterators for printing
    typename Field::const_iterator cfi1 = fi1;
    typename Field::const_iterator cif2 = if2;

    // end condition for each test: index < axisLength && (result || print)
    //  this ends the loops early if and only if the result has been found to be false
    //                                       AND printing has been turned only
    for(int z = 0; z < zLength && (result || print); z++) {
        for(int y = 0; y < yLength && (result || print); y++) {
            for(int x = 0; x < xLength && (result || print); x++, fi1++, if2++) {
                bool compare = (*fi1 == *if2);
                result = result && compare;
                if(display) cout << compare << " ";
            };
            if(print) {
                cout << "\t\t";
                for(int x = 0; x < xLength; x++, cfi1++) {
                    cout << *cfi1 << " ";
                };
                cout << "\t\t";
                for(int x = 0; x < xLength; x++, cif2++) {
                    cout << *cif2 << " ";
                };
                cout << endl;
            };
        };
        if(print) cout << endl;
    };

    return result;
};

template<typename Field>
inline bool display_fields_compare(Field const & field1,
                                   Field const & field2,
                                   bool display = false,
                                   bool print = false) {
    return internal_display_fields_compare<Field>(field1.begin(),
                                                  field2.begin(),
                                                  field1.window_with_ghost(),
                                                  display,
                                                  print);
};

template<typename Field>
inline bool interior_display_fields_compare(Field const & field1,
                                            Field const & field2,
                                            bool display = false,
                                            bool print = false) {
    return internal_display_fields_compare<Field>(field1.interior_begin(),
                                                  field2.interior_begin(),
                                                  field1.window_without_ghost(),
                                                  display,
                                                  print);
};

