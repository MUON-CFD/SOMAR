#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TUPLE_H_
#define _TUPLE_H_

#include <cstdlib>
#include "CH_assert.H"

#include "BaseNamespaceHeader.H"

using namespace std;

//
/// Ordered Tuples for Types T
/**

  This class represents ordered tuples of some user-specified concrete
  type T for N > 0. The type T must have a default constructor.  If the
  non-default constructor, copy constructor, or copy assignment operator
  are used, T must also have a copy constructor.
*/

template <class T, size_t N>
class Tuple
{
public:
    /**: The default constructor.  For user-defined types T, the
               default constructor for T will be run on each of the N
               objects in the Tuple.  For builtin (intrinsic) types,
               the values in the Tuple will be garbage.
    */
#if defined(__INTEL_COMPILER) || defined(__clang__)
    template <class S = T, typename std::enable_if<std::is_constructible<S>::value>::type * = nullptr>
#else 
    template <class S = T, typename std::enable_if<std::is_copy_constructible<S>::value>::type * = nullptr>
#endif
    Tuple() : vect{}
    {
    }
#if defined(__INTEL_COMPILER) || defined(__clang__)
    template <class S = T, typename std::enable_if<!std::is_constructible<S>::value>::type * = nullptr>
#else 
    template <class S = T, typename std::enable_if<!std::is_copy_constructible<S>::value>::type * = nullptr>
#endif
    Tuple() 
    {
    }
    /**: Constructs a Tuple, initializing the elements in the Tuple
               with the corresponding elements in the vector v.  This assumes
               that v contains at least N elements of type T -- an assumption
               that is NOT checked.  For user-defined types, T must have a
               well-defined and accessible copy constructor.
    */
    //explicit Tuple (const T* v);
    //
    // The copy constructor.
    //
    template <class S = T, typename std::enable_if<std::is_copy_assignable<S>::value>::type * = nullptr>
    Tuple(const Tuple &rhs)
    {
        for (size_t i = 0; i < N; ++i)
            vect[i] = rhs.vect[i];
    }
    template <class S = T, typename std::enable_if<!std::is_copy_assignable<S>::value>::type * = nullptr>
    Tuple(const Tuple &rhs)
    {
       static_assert(std::is_copy_assignable<S>::value, "Tuple<T> copy constructor called on a class which is not copy assignable");
    }
    //
    // The copy assignment operator.
    //
    template <class S = T, typename std::enable_if<std::is_copy_assignable<S>::value>::type * = nullptr>
    Tuple &operator=(const Tuple &rhs)
    {
        for (size_t i = 0; i < N; ++i)
            vect[i] = rhs.vect[i];
        return *this;
    }
    template <class S = T, typename std::enable_if<!std::is_copy_assignable<S>::value>::type * = nullptr>
    Tuple &operator=(const Tuple &rhs)
    {
        static_assert(std::is_copy_assignable<S>::value, "Tuple<T> assignment operator called on a class which is not copy assignable");
    }

    /**: Returns a modifiable lvalue reference to the i'th
               element in the Tuple, counting from zero.  Performs
               range checking when the library is compiled in debug
               mode.  */
    T &operator[](int i)
    {
        CH_assert(0 <= i && i < (int)N);
        return vect[i];
    }

    /**: Returns a constant reference to the i'th element in the Tuple,
               counting from zero.  Performs range checking when the library
               is compiled in debug mode.
    */
    const T &operator[](int i) const
    {
        CH_assert(0 <= i && i < (int)N);
        return vect[i];
    }

    /**: Returns the address of the underlying vector of T
               representation.  This should ONLY be used when interfacing
               to Fortran as it breaks the encapsulation of the class.
    */
    operator const T *() const
    {
        return &vect[0];
    }

protected:
    //
    // The underlying vector of T representing the Tuple.
    //
    T vect[N];
};

//
// Inlines.
//

//template <class T, size_t N>
//Tuple<T,N>::Tuple (const T* v)
//{
//    CH_assert(v != 0);
//    for (size_t i = 0; i < N; ++i)
//        vect[i] = v[i];
//}

#include "BaseNamespaceFooter.H"
#endif /*CH_TUPLE_H*/
