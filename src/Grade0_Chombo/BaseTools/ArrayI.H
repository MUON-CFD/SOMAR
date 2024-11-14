#ifndef M4b858ba808f4a6a391f8ea74ce227cdf
#   error Do not include ArrayI.H directly. Use Array.H instead.
#endif

#ifndef USE_ARRAYI_H__
#define USE_ARRAYI_H__

// -----------------------------------------------------------------------------
template <typename T>
T**
create2DArray(const size_t a_N, const size_t a_M)
{
    if (a_N == 0) throw std::invalid_argument("number of rows is 0");
    if (a_M == 0) throw std::invalid_argument("number of columns is 0");

    T** ptr  = nullptr;
    T*  pool = nullptr;

    try {
        // Allocate pointers
        ptr  = new T*[a_N];

        // Allocate pool
        pool = new T[a_N * a_M];

        // Assign pointers to the pool.
        ptr[0] = pool;
        for (size_t i = 1; i < a_N; ++i) {
            ptr[i] = &pool[i * a_M];
        }

        // Done.
        return ptr;

    } catch (std::bad_alloc& ex) {
        // memory allocation error
        delete[] ptr;
        throw ex;
    }
}


// -----------------------------------------------------------------------------
template <typename T>
void
delete2DArray(T** a_array)
{
    if (a_array) {
        delete[] a_array[0];  // remove the pool
        delete[] a_array;     // remove the pointers
        a_array = nullptr;
    }
}


// -----------------------------------------------------------------------------
template <typename T>
T*
create1DArray(const size_t a_N)
{
    if (a_N == 0) throw std::invalid_argument("number of elements is 0");

    T* ptr = nullptr;

    try {
        ptr  = new T[a_N];
        return ptr;

    } catch (std::bad_alloc& ex) {
        // memory allocation error
        delete[] ptr;
        throw ex;
    }
}


// -----------------------------------------------------------------------------
template <typename T>
void
delete1DArray(T* a_array)
{
   delete[] a_array;     // remove the pointers
   a_array = nullptr;
}


#endif //!USE_ARRAYI_H__
