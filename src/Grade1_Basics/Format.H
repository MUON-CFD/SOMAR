/*******************************************************************************
 *    FILE: Format.H
 *    DESC: A set of text and stream formatting tools.
 *    DATE: Jan 5, 2016 1:48 PM EST
 *    MAIL: santillie@philau.edu
 ******************************************************************************/
#ifndef ___Format_H__INCLUDED___
#define ___Format_H__INCLUDED___

#include <array>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stack>
#include <vector>
#include "GNUC_Extensions.H"
#include "RealVect.H"
#include "RefCountedPtr.H"
#include "MayDay.H"

namespace Format {

// -----------------------------------------------------------------------------
// A stream buffer that indents each '\n'. The original stream is restored
// once the Indenter goes out of scope. Also, Indenters can be piled on top of
// one another.
//
// Credit goes to James Kanze from stackoverflow.com for his excellent code.
// http://stackoverflow.com/questions/9599807/how-to-add-indention-to-the-stream-operator
// -----------------------------------------------------------------------------
class Indenter: public std::streambuf
{
public:
    explicit Indenter (std::streambuf* dest,
                       const int       indent = 2,
                       const char*     bullet = "");

    explicit Indenter (std::ostream& dest,
                       const int     indent = 2,
                       const char*   bullet = "");

    virtual ~Indenter ();

protected:
    virtual int overflow (int ch);

    std::streambuf*     m_dest;
    bool                m_isAtStartOfLine;
    std::string         m_indent;
    std::ostream*       m_owner;
};


// -----------------------------------------------------------------------------
// Stream manipulator that indents and/or bullets all incoming stream data until
// unindent is used.
// NOTE: The syntax needs to be corrected. Here is a current use case:
// std::cout << indent() << "This is indented." << unindent << std::endl;
// Notice that intent() needs parentheses while unindent does not.
// -----------------------------------------------------------------------------
class indent {
    friend std::ostream& unindent (std::ostream& os);
    friend std::ostream& unindentAll (std::ostream& os);

public:
    inline indent (const int idt = 2, const char* const bullet = "")
    : m_indent(idt), m_bullet(bullet) {}

    void operator() (std::ostream& os) const;

private:
    int         m_indent;
    const char* m_bullet;
    static std::stack<RefCountedPtr<Indenter> > s_indentStack;
};

std::ostream& unindent (std::ostream& os);
std::ostream& unindentAll (std::ostream& os);
std::ostream& operator<< (std::ostream& os, const indent& idt);


// -----------------------------------------------------------------------------
// Stream manipulator that formats numeric data. Use these instead of
// std::fixed and std::scientific to provide uniform formatting.
// -----------------------------------------------------------------------------
std::ostream& fixed (std::ostream& os);
std::ostream& scientific (std::ostream& os);


// -----------------------------------------------------------------------------
// Stream manipulator that chooses fixed or scientific, based on the number.
// This looks like a mess, but the code mimics std::setprecision. So it works.
// -----------------------------------------------------------------------------
struct _Number {Real _M_n;};
inline _Number number(Real __n) {return { __n };} // The manipulator.

template <typename _CharT, typename _Traits>
inline std::basic_ostream<_CharT, _Traits>&
operator<<(std::basic_ostream<_CharT, _Traits>& __os, _Number __f)
{
    if (1.0e-3 < std::abs(__f._M_n) && std::abs(__f._M_n) < 1.0e3) {
        return __os << std::fixed << std::setprecision(6) << __f._M_n;
    } else {
        return __os << std::scientific << std::setprecision(6) << __f._M_n;
    }
}


// -----------------------------------------------------------------------------
// Stream manipulators that stores/restores a streams flags.
// NOTE: This does NOT store or restore the indentation state.
// -----------------------------------------------------------------------------
std::ostream& pushFlags (std::ostream& os);
std::ostream& popFlags (std::ostream& os);


// -----------------------------------------------------------------------------
// Creates a section heading for use by a stream (typically, cout or pout).
// -----------------------------------------------------------------------------
std::string banner(const std::string& a_str);


// -----------------------------------------------------------------------------
// Converts a time in seconds (double) into the more readable days, hours, mins,
// and seconds (std::string).
// -----------------------------------------------------------------------------
std::string textTime (const double inSeconds);


// -----------------------------------------------------------------------------
// Writes a formatted warning message to an ostream.
// -----------------------------------------------------------------------------
class Warning {
public:
    Warning (const char* msg, const char* file = "", const int line = 0)
    : m_msg(msg), m_file(file), m_line(line) {}

    void operator() (std::ostream& os) const;

private:
    static const char* const s_prefix;
    const char* m_msg;
    const char* m_file;
    int         m_line;
};

std::ostream& operator<< (std::ostream& os, const Warning& warning);

#define WARNING(m)                                      \
    do {                                                \
        ostringstream msg;                              \
        msg << ::Format::Warning(m, __FILE__, __LINE__);\
        MayDay::Warning(msg.str().c_str());             \
    } while(0)


// -----------------------------------------------------------------------------
// Writes a formatted error message to an ostream.
// -----------------------------------------------------------------------------
class Error {
public:
    Error (const char* msg, const char* file = "", const int line = 0)
    : m_msg(msg), m_file(file), m_line(line) {}

    void operator() (std::ostream& os) const;

private:
    static const char* const s_prefix;
    const char* m_msg;
    const char* m_file;
    int         m_line;
};

std::ostream& operator<< (std::ostream& os, const Error& error);

#define ERROR(m)                                        \
    do {                                                \
        ostringstream msg;                              \
        msg << ::Format::Error(m, __FILE__, __LINE__);  \
        MayDay::Error(msg.str().c_str());               \
    } while(0)



// -----------------------------------------------------------------------------
// Color stream manipulators
// -----------------------------------------------------------------------------
// #define USE_COLOR_IN_DEBUG
#define USE_COLOR_IN_RELEASE

#if ((!defined NDEBUG) && (defined USE_COLOR_IN_DEBUG)) || ((defined NDEBUG) && (defined USE_COLOR_IN_RELEASE))
    // Return to default color
    static const char __nowarn_unused *none     = "\033[0m";

    // Low intensity colors
    static const char __nowarn_unused *black    = "\033[0;30m";
    static const char __nowarn_unused *blue     = "\033[0;34m";
    static const char __nowarn_unused *green    = "\033[0;32m";
    static const char __nowarn_unused *cyan     = "\033[0;36m";
    static const char __nowarn_unused *red      = "\033[0;31m";
    static const char __nowarn_unused *purple   = "\033[0;35m";
    static const char __nowarn_unused *brown    = "\033[0;33m";
    static const char __nowarn_unused *white    = "\033[0;37m";

    // High intensity colors
    static const char __nowarn_unused *hiblack  = "\033[1;30m";
    static const char __nowarn_unused *hiblue   = "\033[1;34m";
    static const char __nowarn_unused *higreen  = "\033[1;32m";
    static const char __nowarn_unused *hicyan   = "\033[1;36m";
    static const char __nowarn_unused *hired    = "\033[1;31m";
    static const char __nowarn_unused *hipurple = "\033[1;35m";
    static const char __nowarn_unused *hibrown  = "\033[1;33m";
    static const char __nowarn_unused *hiwhite  = "\033[1;37m";

    // Aliases
    static const char __nowarn_unused *yellow   = hibrown;
#else
    // Return to default color
    static const char __nowarn_unused *none     = "";

    // Low intensity colors
    static const char __nowarn_unused *black    = "";
    static const char __nowarn_unused *blue     = "";
    static const char __nowarn_unused *green    = "";
    static const char __nowarn_unused *cyan     = "";
    static const char __nowarn_unused *red      = "";
    static const char __nowarn_unused *purple   = "";
    static const char __nowarn_unused *brown    = "";
    static const char __nowarn_unused *white    = "";

    // High intensity colors
    static const char __nowarn_unused *hiblack  = "";
    static const char __nowarn_unused *hiblue   = "";
    static const char __nowarn_unused *higreen  = "";
    static const char __nowarn_unused *hicyan   = "";
    static const char __nowarn_unused *hired    = "";
    static const char __nowarn_unused *hipurple = "";
    static const char __nowarn_unused *hibrown  = "";
    static const char __nowarn_unused *hiwhite  = "";

    // Aliases
    static const char __nowarn_unused *yellow   = hibrown;
#endif


// -----------------------------------------------------------------------------
// Used to make text more readable.
// -----------------------------------------------------------------------------
std::string stripFuncName(const std::string& a_func);
std::string stripFileName(const std::string& a_file);

} // namespace Format


// These functions must be at global scope...

// -----------------------------------------------------------------------------
// Sends an IntVect to an ostream.
// -----------------------------------------------------------------------------
inline std::ostream&
operator<<(std::ostream& stream, const IntVect& v)
{
    stream << "("
           D_TERM( << v[0],
                   << ", " << v[1],
                   << ", " << v[2] )
           << ")";
    return stream;
}


// // -----------------------------------------------------------------------------
// // Sends a RealVect to an ostream.
// // -----------------------------------------------------------------------------
// inline std::ostream&
// operator<<(std::ostream& stream, const RealVect& v)
// {
//     stream << "("
//            D_TERM( << v[0],
//                    << ", " << v[1],
//                    << ", " << v[2] )
//            << ")";
//     return stream;
// }


// -----------------------------------------------------------------------------
// Sends a std::vector to an ostream.
// -----------------------------------------------------------------------------
template<typename T>
inline std::ostream&
operator<<(std::ostream& stream, const std::vector<T>& v)
{
    if (!v.empty()) {
        stream << "(" << v[0];
        for (size_t idx = 1; idx < v.size(); ++idx)
            stream << ", " << v[idx];
        stream << ")";
    } else {
        stream << "()";
    }
    return stream;
}


// -----------------------------------------------------------------------------
// Sends a std::array to an ostream.
// -----------------------------------------------------------------------------
template<typename T, size_t N>
inline std::ostream&
operator<<(std::ostream& stream, const std::array<T, N>& a)
{
    if constexpr (N > 0) {
        stream << "[" << a[0];
        for (size_t idx = 1; idx < N; ++idx)
            stream << ", " << a[idx];
        stream << "]";
    } else {
        stream << "[]";
    }
    return stream;
}

#endif //!___Format_H__INCLUDED___
