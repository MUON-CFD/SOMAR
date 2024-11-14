#pragma once
#include <chrono>
#include "BasicIO.H"
// A simple utility function to time the execution of a procedure. Note that what is timed is f().
// Therefore, to time something that requires input values we can wrap it into a lambda function.
// For example, if we want to time
// T r=MyClass.MyFunc(a,b,c), where a and b are passed by reference, and c by value
// we do the following
// T r=Util::FuncTimer(" MyFunc ", [&MyClass, &a, &b, c]() ->T{ return MyClass.MyFunc(a,b,c);});
// T can be a reference to a type (e.g., int&), or can be void. In the latter case, the lambda function does not
// obviously return anything and the call would be 
// Util::FuncTimer(" MyFunc ", [&MyClass, &a, &b, c]() ->void {MyClass.MyFunc(a,b,c);});
// The interface can be simplified to capture the entire scope by reference
// T r=Util::FuncTimer(" MyFunc ", [&]() ->T{ return MyClass.MyFunc(a,b,c);});

namespace Util
{


template <typename F>
auto FuncTimer(std::string Name, F f) -> typename std::enable_if<!std::is_same<decltype(f()), void>::value, decltype(f())>::type
{
    auto start = std::chrono::high_resolution_clock::now();
    decltype(f()) r= f();
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    BasicIO::tout(0) << "time in function " << Name << ": " << elapsed.count() << " s" << endl;
    return r;
}
template<typename F>
auto FuncTimer(std::string Name, F f) -> typename std::enable_if<std::is_same<decltype(f()), void>::value, decltype(f())>::type
{
    auto start = std::chrono::high_resolution_clock::now();
    f();
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    BasicIO::tout(0) << "time in function " << Name << ": " << elapsed.count() << " s" << endl;
    
}

} // namespace Util

