/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2024 Thomas Jefferson University and Arizona State University
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 *  USA
 *
 *  For up-to-date contact information, please visit the repository homepage,
 *  https://github.com/somarhub.
 ******************************************************************************/
#pragma once
#include <tuple>
#include <vector>
#include <type_traits>


struct enumerator
{
    size_t b;
    size_t end;
    enumerator() : b(0), end(SIZE_MAX) {}
    enumerator(const size_t _i) : b(0), end(_i) {}
    void operator++() { ++b; }
    bool operator!=(const void *) { return b != end; }
    size_t &operator*() { return b; }
};

struct enumerated
{
    enumerator b;

    enumerated(const size_t size = 1) : b(size) {}
    auto begin() { return b; }
    auto begin() const { return b; }
    void *end() { return nullptr; }
    void *end() const { return nullptr; }
};



namespace Utils {

template<typename _Tp>
using remove_cv_t = typename std::remove_cv<_Tp>::type;

template<bool _Cond, typename _Tp = void>
using enable_if_t = typename std::enable_if<_Cond, _Tp>::type;

template <typename _Tp, typename _Up>
inline constexpr bool is_same_v = std::is_same<_Tp, _Up>::value;


// A C++11 implementation of std::integer_sequence from C++14
template <typename T, T... N>
struct integer_sequence {
    typedef T value_type;
    static_assert(
        std::is_integral<T>::value,
        "std::integer_sequence can only be instantiated with an integral type");

    static inline std::size_t
    size()
    {
        return (sizeof...(N));
    }
};

template <std::size_t... N>
using index_sequence = integer_sequence<std::size_t, N...>;

namespace integer_sequence_detail
{

template <typename T, std::size_t... _Extra>
struct repeat;

template <typename T, T... N, std::size_t... _Extra>
struct repeat<integer_sequence<T, N...>, _Extra...> {
    typedef integer_sequence<T,
                             N...,
                             1 * sizeof...(N) + N...,
                             2 * sizeof...(N) + N...,
                             3 * sizeof...(N) + N...,
                             4 * sizeof...(N) + N...,
                             5 * sizeof...(N) + N...,
                             6 * sizeof...(N) + N...,
                             7 * sizeof...(N) + N...,
                             _Extra...>
        type;
};

template <std::size_t N>
struct parity;
template <std::size_t N>
struct make : parity<N % 8>::template pmake<N> {
};

template <>
struct make<0> {
    typedef integer_sequence<std::size_t> type;
};
template <>
struct make<1> {
    typedef integer_sequence<std::size_t, 0> type;
};
template <>
struct make<2> {
    typedef integer_sequence<std::size_t, 0, 1> type;
};
template <>
struct make<3> {
    typedef integer_sequence<std::size_t, 0, 1, 2> type;
};
template <>
struct make<4> {
    typedef integer_sequence<std::size_t, 0, 1, 2, 3> type;
};
template <>
struct make<5> {
    typedef integer_sequence<std::size_t, 0, 1, 2, 3, 4> type;
};
template <>
struct make<6> {
    typedef integer_sequence<std::size_t, 0, 1, 2, 3, 4, 5> type;
};
template <>
struct make<7> {
    typedef integer_sequence<std::size_t, 0, 1, 2, 3, 4, 5, 6> type;
};

template <>
struct parity<0> {
    template <std::size_t N>
    struct pmake : repeat<typename make<N / 8>::type> {
    };
};
template <>
struct parity<1> {
    template <std::size_t N>
    struct pmake : repeat<typename make<N / 8>::type, N - 1> {
    };
};
template <>
struct parity<2> {
    template <std::size_t N>
    struct pmake : repeat<typename make<N / 8>::type, N - 2, N - 1> {
    };
};
template <>
struct parity<3> {
    template <std::size_t N>
    struct pmake : repeat<typename make<N / 8>::type, N - 3, N - 2, N - 1> {
    };
};
template <>
struct parity<4> {
    template <std::size_t N>
    struct pmake
    : repeat<typename make<N / 8>::type, N - 4, N - 3, N - 2, N - 1> {
    };
};
template <>
struct parity<5> {
    template <std::size_t N>
    struct pmake
    : repeat<typename make<N / 8>::type, N - 5, N - 4, N - 3, N - 2, N - 1> {
    };
};
template <>
struct parity<6> {
    template <std::size_t N>
    struct pmake : repeat<typename make<N / 8>::type,
                          N - 6,
                          N - 5,
                          N - 4,
                          N - 3,
                          N - 2,
                          N - 1> {
    };
};
template <>
struct parity<7> {
    template <std::size_t N>
    struct pmake : repeat<typename make<N / 8>::type,
                          N - 7,
                          N - 6,
                          N - 5,
                          N - 4,
                          N - 3,
                          N - 2,
                          N - 1> {
    };
};

template <typename T, typename U>
struct convert {
    template <typename>
    struct result;

    template <T... N>
    struct result<integer_sequence<T, N...>> {
        typedef integer_sequence<U, N...> type;
    };
};

template <typename T>
struct convert<T, T> {
    template <typename U>
    struct result {
        typedef U type;
    };
};

template <typename T, T N>
using make_integer_sequence_unchecked =
    typename convert<std::size_t,
                     T>::template result<typename make<N>::type>::type;

template <typename T, T N>
struct make_integer_sequence {
    static_assert(std::is_integral<T>::value,
                  "std::make_integer_sequence can only be instantiated with an "
                  "integral type");
    static_assert(0 <= N,
                  "std::make_integer_sequence input shall not be negative");

    typedef make_integer_sequence_unchecked<T, N> type;
};

}  // namespace integer_sequence_detail

template <typename T, T N>
using make_integer_sequence =
    typename integer_sequence_detail::make_integer_sequence<T, N>::type;

template <std::size_t N>
using make_index_sequence = make_integer_sequence<std::size_t, N>;

template <typename... T>
using index_sequence_for = make_index_sequence<sizeof...(T)>;

template <class...>
struct conjunction : std::true_type {
};
template <class B1>
struct conjunction<B1> : B1 {
};
template <class B1, class... Bn>
struct conjunction<B1, Bn...>
: std::conditional<bool(B1::value), conjunction<Bn...>, B1>::type {
};



template <typename T>
struct _is_vector : std::false_type {
};

template <typename... Args>
struct _is_vector<std::vector<Args...>> : std::true_type {
};

template <typename... Args>
struct _is_vector<const std::vector<Args...>> : std::true_type {
};

template <>
struct _is_vector<enumerated> : std::true_type {
};
template <>
struct _is_vector<const enumerated> : std::true_type {
};


template <typename, typename>
struct prepend_to_type_seq {
};

template <typename T, typename... Ts, template <typename...> class TT>
struct prepend_to_type_seq<T, TT<Ts...>> {
    using type = TT<T, Ts...>;
};

template <typename>
struct extract_first_type {
};

template <typename... Ts, template <typename...> class TT>
struct extract_first_type<TT<Ts...>> {
    using type = typename std::tuple_element<0, std::tuple<Ts...>>::type;
};

template <typename... Ts>
constexpr bool are_vectors = conjunction<_is_vector<Ts>...>::value;

template <std::size_t N>
struct num {
    static const constexpr auto value = N;
};

template <class F, std::size_t... Is>
void
for_(F func, index_sequence<Is...>)
{
    (func(num<Is>{}), ...);
}

template <std::size_t N, typename F>
void
for_(F func)
{
    for_(func, make_index_sequence<N>());
}

template <typename Tuple, std::size_t... Is>
constexpr auto
getElementsAtPointersInTuple_(Tuple& tuple, index_sequence<Is...>)
{
    return std::tuple_cat(std::tie(*(std::get<num<Is>::value>(tuple)))...);
}

template <typename Tuple>
constexpr auto
getElementsAtPointersInTuple(Tuple& tuple)
{
    constexpr size_t N = std::tuple_size<Tuple>::value;
    return getElementsAtPointersInTuple_(tuple, make_index_sequence<N>());
}

}  // end namespace Utils





// -----------------------------------------------------------------------------
template <typename Tuple, typename... Ts>
struct zip_iterator {
private:
    std::tuple<decltype(Ts().begin())...> element;
    std::tuple<decltype(Ts().end())...>   ends;

    static constexpr size_t nElements = std::tuple_size<Tuple>::value;

public:
    zip_iterator()                                        = delete;
    zip_iterator(const zip_iterator<Tuple, Ts...>& other) = default;
    zip_iterator(zip_iterator<Tuple, Ts...>&& other)      = default;
    zip_iterator(const Tuple& zipped)
    {
        Utils::for_<nElements>([&](auto i) {
            std::get<i.value>(element) = std::get<i.value>(zipped).begin();
            std::get<i.value>(ends)    = std::get<i.value>(zipped).end();
        });
    };
    void
    operator++()
    {
        Utils::for_<nElements>([&](auto i) { ++std::get<i.value>(element); });
    }

    bool
    operator!=(const zip_iterator<Tuple, Ts...>& other)
    {
        bool dummy = true;
        Utils::for_<nElements>([&](auto i) {
            dummy = dummy && (std::get<i.value>(element) !=
                              std::get<i.value>(other.ends));
        });
        return dummy;
    }

    auto
    operator*()
    {
        return Utils::getElementsAtPointersInTuple(element);
    }
};


template <typename... Ts>
struct _vector_zip {
    std::tuple<Ts&...> my_zip;

    _vector_zip(Ts&... ts) : my_zip(std::tuple_cat(std::tie(ts)...))
    {
    }

    auto
    begin() const
    {
        return zip_iterator<std::tuple<Ts&...>, Ts...>(my_zip);
    }
    auto
    end() const
    {
        return this->begin();
    }
};

template <typename... Ts>
struct zip : public _vector_zip<Ts...> {
    zip(Ts&... args) : _vector_zip<Ts...>(args...)
    {
    }
};

template <typename... Ts>
struct enumerate {
    const enumerated e;
    const typename Utils::prepend_to_type_seq<decltype(e), zip<Ts...>>::type
        zipper;
    enumerate(Ts&... args) : e(SIZE_MAX), zipper(e, args...)
    {
    }
    auto
    begin() const
    {
        return zipper.begin();
    }
    auto
    end() const
    {
        return zipper.end();
    }
};
