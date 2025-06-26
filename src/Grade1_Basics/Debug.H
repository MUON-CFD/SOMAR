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
 *  https://github.com/MUON-CFD/SOMAR.
 ******************************************************************************/

// This file includes all of the debugging utilities.

#ifndef __DEBUG_H__INCLUDED__
#define __DEBUG_H__INCLUDED__

#include <array>
#include <sstream>
#include <string>
#include "FArrayBox.H"
#include "NodeFArrayBox.H"
#include "FluxBox.H"
#include "Format.H"
#include "LevelData.H"
#include "LayoutTools.H"
#include "MayDay.H"
#include "BasicIO.H"

#ifndef NDEBUG
#include "DebugOut.H"
#endif


// Like CH_assert, but also works in release mode.
#define CH_verify_str(s) #s
#define CH_verify_xstr(s) CH_verify_str(s)
#define CH_verify(cond)                                         \
    if (!(cond)) {                                              \
        MayDay::Error(__FILE__ ":" CH_verify_xstr(              \
            __LINE__) ": Verification of `" #cond "' failed."); \
    }


// Nicer than using #ifndef NDEBUG...
#ifdef NDEBUG
#define DEBUG_RAWCODE(x)
#define DEBUG_CODE(x)
#else
#define DEBUG_RAWCODE(x) x
#define DEBUG_CODE(x) \
    {                 \
        x             \
    }
#endif


// Allows inline assignments that depend on release mode
// For example:
// const char* dbgMsg = IF_DEBUG("Debug mode!", "Release mode!");
#ifndef NDEBUG
#define IF_DEBUG(d, r) d
#else
#define IF_DEBUG(d, r) r
#endif

static bool debugMode __attribute__((unused))  = IF_DEBUG(true, false);
static bool releaseMode __attribute__((unused))= IF_DEBUG(false, true);

// This will exist in release mode as well.
#define POUT(x)                                    \
    do {                                           \
        pout() << #x << " = " << (x) << std::endl; \
    } while (0)
#define TOUT(r, x)                                           \
    do {                                                     \
        BasicIO::tout((r)) << #x << " = " << x << std::endl; \
    } while (0)


// -----------------------------------------------------------------------------
// Puts a checkpoint into pout.*
// -----------------------------------------------------------------------------
#define CHECKPOINT() \
    _WriteCheckpoint(__PRETTY_FUNCTION__, __FILE__, __LINE__)
#define CHECKMARK(mark) \
    _WriteCheckpoint(__PRETTY_FUNCTION__, __FILE__, __LINE__, mark)
void
_WriteCheckpoint(std::string  a_func,
                 std::string  a_file,
                 unsigned int a_line,
                 std::string  a_mark = "");


// -----------------------------------------------------------------------------
// Throws an error and reports location
// -----------------------------------------------------------------------------
#define UNDEFINED_FUNCTION()                                           \
    do {                                                               \
        ostringstream msg1, msg2;                                      \
        msg1 << __PRETTY_FUNCTION__ << " is undefined";                \
        msg2 << Format::Error(msg1.str().c_str(), __FILE__, __LINE__); \
        MayDay::Error(msg2.str().c_str());                             \
    } while (0)


// -----------------------------------------------------------------------------
// Write a note/todo to the terminal.
// -----------------------------------------------------------------------------
void
_alwaysNote(const char*       a_filename,
            const char*       a_funcname,
            const int         a_linenumber,
            const std::string a_msg,
            const std::string a_premsg = "");

#ifndef NDEBUG
inline void
_debugNote(const char*       a_filename,
           const char*       a_funcname,
           const int         a_linenumber,
           const std::string a_msg,
           const std::string a_premsg = "")
{
    _alwaysNote(a_filename, a_funcname, a_linenumber, a_msg, a_premsg);
}

inline void
_releaseNote(const char*       /*a_filename*/,
             const char*       /*a_funcname*/,
             const int         /*a_linenumber*/,
             const std::string /*a_msg*/,
             const std::string /*a_premsg*/ = "")
{}
#else
inline void
_debugNote(const char*       /*a_filename*/,
           const char*       /*a_funcname*/,
           const int         /*a_linenumber*/,
           const std::string /*a_msg*/,
           const std::string /*a_premsg*/ = "")
{}

inline void
_releaseNote(const char*       a_filename,
             const char*       a_funcname,
             const int         a_linenumber,
             const std::string a_msg,
             const std::string a_premsg = "")
{
    _alwaysNote(a_filename, a_funcname, a_linenumber, a_msg, a_premsg);
}
#endif

#define NOTE(msg)                                                          \
    do {                                                                   \
        static bool w = true;                                              \
        if (w) {                                                           \
            ostringstream ss;                                              \
            ss << msg;                                                     \
            _debugNote(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str()); \
            w = false;                                                     \
        }                                                                  \
    } while (0)

#define TODO()                                                                 \
    do {                                                                       \
        static bool w = true;                                                  \
        if (w) {                                                               \
            ostringstream ss;                                                  \
            ss << Format::brown << "TODO: ";                                   \
            _debugNote(__FILE__, __PRETTY_FUNCTION__, __LINE__, "", ss.str()); \
            w = false;                                                         \
        }                                                                      \
    } while (0)

#define TODONOTE(msg)                           \
    do {                                        \
        static bool w = true;                   \
        if (w) {                                \
            ostringstream sspre;                \
            sspre << Format::brown << "TODO: "; \
            ostringstream ss;                   \
            ss << msg;                          \
            _debugNote(__FILE__,                \
                  __PRETTY_FUNCTION__,          \
                  __LINE__,                     \
                  ss.str(),                     \
                  sspre.str());                 \
            w = false;                          \
        }                                       \
    } while(0)

#define BUG(msg)                                \
    do {                                        \
        static bool w = true;                   \
        if (w) {                                \
            ostringstream sspre;                \
            sspre << Format::hired << "*BUG: "; \
            ostringstream ss;                   \
            ss << msg;                          \
            _alwaysNote(__FILE__,               \
                        __PRETTY_FUNCTION__,    \
                        __LINE__,               \
                        ss.str(),               \
                        sspre.str());           \
            w = false;                          \
        }                                       \
    } while (0)

#define SIMPLEWARNING(msg)                             \
    do {                                               \
        static bool w = true;                          \
        if (w) {                                       \
            ostringstream sspre;                       \
            sspre << Format::yellow << "\n*WARNING: "; \
            ostringstream ss;                          \
            ss << msg << '\n';                         \
            _alwaysNote(__FILE__,                      \
                        __PRETTY_FUNCTION__,           \
                        __LINE__,                      \
                        ss.str(),                      \
                        sspre.str());                  \
            w = false;                                 \
        }                                              \
    } while (0)

#define MAYDAYWARNING(usrmsg)                                                  \
    do {                                                                       \
        static bool warned = false;                                            \
        if (warned) break;                                                     \
        pout() << Format::unindentAll;                                         \
        ostringstream msg;                                                     \
        msg << Format::hired << "***** WARNING *****\n"                        \
            << Format::hiwhite << "File: " << Format::none                     \
            << Format::stripFileName(__FILE__) << '\n'                         \
            << Format::hiwhite << "Line: " << Format::none << __LINE__ << '\n' \
            << Format::hiwhite << "Func: " << Format::none                     \
            << Format::stripFuncName(__PRETTY_FUNCTION__) << '\n'              \
            << Format::hiwhite << "Note: " << Format::none << usrmsg           \
            << Format::none << '\n';                                           \
        MayDay::Warning(msg.str().c_str());                                    \
        warned = true;                                                         \
    } while (0)

#define MAYDAYERROR(usrmsg)                                                    \
    do {                                                                       \
        pout() << Format::unindentAll;                                         \
        ostringstream msg;                                                     \
        msg << Format::hired << "***** ERROR *****\n"                          \
            << Format::hiwhite << "File: " << Format::none                     \
            << Format::stripFileName(__FILE__) << '\n'                         \
            << Format::hiwhite << "Line: " << Format::none << __LINE__ << '\n' \
            << Format::hiwhite << "Func: " << Format::none                     \
            << Format::stripFuncName(__PRETTY_FUNCTION__) << '\n'              \
            << Format::hiwhite << "Note: " << Format::none << usrmsg           \
            << Format::none << '\n';                                           \
        MayDay::Error(msg.str().c_str());                                      \
    } while (0)

// -----------------------------------------------------------------------------
// Put this in all functions that need to be tested.
// This will not throw and error, just display a message at runtime.
// -----------------------------------------------------------------------------
#define TEST()                                              \
    {                                                       \
        static bool w = true;                               \
        if (w) {                                            \
            _test(__FILE__, __PRETTY_FUNCTION__, __LINE__); \
            w = false;                                      \
        }                                                   \
    }
    IF_DEBUG(, inline) void _test(const char* a_filename   [[maybe_unused]],
                                  const char* a_funcname   [[maybe_unused]],
                                  const int   a_linenumber [[maybe_unused]])
                                  IF_DEBUG(;, { ; })


// -----------------------------------------------------------------------------
// Performs an enhanced compatibility check on 2 DisjointBoxLayouts
// -----------------------------------------------------------------------------
#ifndef NDEBUG
#define DBLCHECK(dbl1, dbl2)             \
    if (!(dbl1 == dbl2)) {               \
        ostringstream s1, s2;            \
        s1 << dbl1;                      \
        s2 << dbl2;                      \
        CH_assert(s1.str() == s2.str()); \
    }                                    \
    CH_assert(dbl1.compatible(dbl2))
#else
#define DBLCHECK(dbl1, dbl2)
#endif


// -----------------------------------------------------------------------------
// Throws an error if a NAN is found in the testBox / valid regions.
// -----------------------------------------------------------------------------
typedef struct {
    bool    problemFound;
    IntVect pos;
    int     comp;
    Box     box;
    Real    val;
    string  msg;
} naninfo_type;

#ifndef NDEBUG
#define checkForNAN(data,testbox) _checkForNAN(data, testbox, __FILE__, __LINE__, nullptr);
#else
#define checkForNAN(data,testbox)
#endif

void
_checkForNAN(const FArrayBox&  a_data,
             const Box&        a_testBox,
             const std::string a_file,
             const int         a_line,
             naninfo_type*     a_info = nullptr);


#ifndef NDEBUG
#define checkForValidNAN(d) _checkForValidNAN(d, __FILE__, __LINE__);
#else
#define checkForValidNAN(d)
#endif

void
_checkForValidNAN(const LevelData<FArrayBox>& a_data,
                  const std::string           a_file,
                  const int                   a_line);

void
_checkForValidNAN(const LevelData<FluxBox>& a_data,
                  const std::string         a_file,
                  const int                 a_line);

void
_checkForValidNAN(const BoxLayoutData<FluxBox>&  a_data,
                  const std::string a_file,
                  const int         a_line);

// -----------------------------------------------------------------------------
// LevelData<FluxBox> overlap checks
// -----------------------------------------------------------------------------
#define DO_debugCheckValidFaceOverlap_IN_DEBUG_MODE
// #define DO_debugCheckValidFaceOverlap_IN_RELEASE_MODE

#if defined(NDEBUG)
#   if defined(DO_debugCheckValidFaceOverlap_IN_RELEASE_MODE)
#       define DO_debugCheckValidFaceOverlap
#   endif
#else
#   if defined(DO_debugCheckValidFaceOverlap_IN_DEBUG_MODE)
#       define DO_debugCheckValidFaceOverlap
#   endif
#endif

#if defined(DO_debugCheckValidFaceOverlap)
#   define debugCheckValidFaceOverlap(l) LayoutTools::checkValidFaceOverlap(l, __FILE__, __LINE__)
#else
#   define debugCheckValidFaceOverlap(l)
#endif


// -----------------------------------------------------------------------------
// Initializes data holders to NAN
// -----------------------------------------------------------------------------
#define USE_INIT_FUNCTIONS_IN_DEBUG_MODE
// #define USE_INIT_FUNCTIONS_IN_RELEASE_MODE

#if defined(NDEBUG)
#   if defined(USE_INIT_FUNCTIONS_IN_RELEASE_MODE)
#       define USE_INIT_FUNCTIONS
#   endif
#else
#   if defined(USE_INIT_FUNCTIONS_IN_DEBUG_MODE)
#       define USE_INIT_FUNCTIONS
#   endif
#endif

#if defined(USE_INIT_FUNCTIONS)
    // Set a_levels[a_min:a_max] to NAN.
    void debugInitLevels(Vector<LevelData<FArrayBox>*>& a_levels,
                         int                            a_min,
                         int                            a_max);

    // Set a_level to NAN.
    void
    debugInitLevel(LevelData<FArrayBox>& a_level);

    // Set a_level to NAN.
    void
    debugInitLevel(BoxLayoutData<FArrayBox>& a_level);

    // Set a_level to NAN.
    void
    debugInitLevel(LevelData<FluxBox>& a_level);

    // Set a_level to NAN.
    void
    debugInitLevel(BoxLayoutData<FluxBox>& a_level);

    // Set a_level to NAN.
    void
    debugInitLevel(LevelData<NodeFArrayBox>& a_level);

    // Set a_level to NAN.
    template <size_t N>
    void
    debugInitLevel(std::array<LevelData<FArrayBox>, N>& a_level)
    {
        for (size_t i = 0; i < N; ++i) {
            debugInitLevel(a_level[i]);
        }
    }

    // Set a_level to NAN.
    template <size_t N1, size_t N2>
    void
    debugInitLevel(std::array<std::array<LevelData<FArrayBox>, N1>, N2>& a_level)
    {
        for (size_t i = 0; i < N1; ++i) {
            debugInitLevel(a_level[i]);
        }
    }

    // Set FArrayBox to NAN.
    void
    debugInit(FArrayBox& a_fab);

    // Set FluxBox to NAN.
    void
    debugInit(FluxBox& a_flub);

    // Set ghosts of a_level to NAN.
    void
    debugInitLevelGhosts(LevelData<FArrayBox>& a_level);

    // Set ghosts of a_level to NAN.
    void
    debugInitLevelGhosts(LevelData<FluxBox>& a_level);

    // Set data of a_fab that lie outside of a_valid to NAN.
    // a_fab MUST have the same centering as a_valid.
    void
    debugInitGhosts(FArrayBox& a_fab, Box a_valid);

    // Set faces of a_fab that lie outside of a_ccValid to NAN.
    void
    debugInitGhosts(FluxBox& a_flub, const Box& a_ccValid);

    // Set physical and CFI boundary faces of a_level to NAN.
    void
    debugInitLevelBoundaryFaces(LevelData<FluxBox>& a_level);

#else

    inline void debugInitLevels(Vector<LevelData<FArrayBox>*>& /*a_levels*/,
                                int                            /*a_min*/,
                                int                            /*a_max*/)
    {}

    inline void
    debugInitLevel(BoxLayoutData<FArrayBox>& /*a_level*/)
    {}

    inline void
    debugInitLevel(LevelData<FArrayBox>& /*a_level*/)
    {}

    inline void
    debugInitLevel(LevelData<FluxBox>& /*a_level*/)
    {}

    inline void
    debugInitLevel(BoxLayoutData<FluxBox>& /*a_level*/)
    {}

    inline void
    debugInitLevel(LevelData<NodeFArrayBox>& /*a_level*/)
    {}

    // Set a_level to NAN.
    template <size_t N>
    void
    debugInitLevel(std::array<LevelData<FArrayBox>, N>& /*a_level*/)
    {}

    // Set a_level to NAN.
    template <size_t N1, size_t N2>
    void
    debugInitLevel(std::array<std::array<LevelData<FArrayBox>, N1>, N2>& /*a_level*/)
    {}

    inline void
    debugInit(FArrayBox& /*a_fab*/)
    {}

    inline void
    debugInit(FluxBox& /*a_flub*/)
    {}

    inline void
    debugInitLevelGhosts(LevelData<FArrayBox>& /*a_level*/)
    {}

    inline void
    debugInitLevelGhosts(LevelData<FluxBox>& /*a_level*/)
    {}

    inline void
    debugInitGhosts(FArrayBox& /*a_fab*/, Box /*a_valid*/)
    {}

    inline void
    debugInitGhosts(FluxBox& /*a_flub*/, const Box& /*a_ccValid*/)
    {}

    inline void
    debugInitLevelBoundaryFaces(LevelData<FluxBox>& /*a_level*/)
    {}
#endif


#endif  //!__DEBUG_H__INCLUDED__
