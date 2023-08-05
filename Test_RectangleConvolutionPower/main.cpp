#include <iostream>
#include "CycleSampler.hpp"
#include "../src/RectangleConvolutionPower.hpp"

using namespace Tools;
using namespace Tensors;
using namespace CycleSampler;

using Real = double;
using Int  = int;

int main(int argc, const char * argv[])
{
    // insert code here...
    std::cout << "Hello, World!" << std::endl;
    
    
    RectangleConvolutionPower<Real,Int> f ( 100 );
    
    dump( f.Value( 0.3 ) );
    
    
    constexpr Int m = 5;
    constexpr Int n = 2;
    
    Int p [m] = {0,1,2,3,4};
    
    Int p_dims[1] = {m};
    
    Real A [m][n] = { {100,200},{101,201},{102,202},{103,203},{104,204} };
    
    Int A_dims [2] = { m, n } ;
    
//        print( ArrayToString( &p[0], &p_dims[0], 1 ) );
//        print( ArrayToString( &A[0], &A_dims[0], 2 ) );

    
    
    print( ArrayToString( &p[0], &p_dims[0], 1 ) );
    
    for( Int iter = 0; iter < Factorial(m); ++iter )
    {
//        print("");
//        dump(iter);
        const bool result = NextPermutation( &A[0][0], &p[0], m, n );
        
        dump(result);
        print( ArrayToString( &p[0], &p_dims[0], 1 ) );
//        print( ArrayToString( &A[0][0], &A_dims[0], 2 ) );
    }

    return 0;
}
