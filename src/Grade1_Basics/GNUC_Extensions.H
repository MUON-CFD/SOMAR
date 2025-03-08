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
#ifndef __GNUC_EXTENSIONS_H__INCLUDED__
#define __GNUC_EXTENSIONS_H__INCLUDED__


// GNU C only!
#ifdef __GNUC__

#	ifndef __nowarn_unused
#		define __nowarn_unused	__attribute__((unused))
#	endif //__nowarn_unused

#	ifndef __deprecated
#		define __deprecated		__attribute__((deprecated))
#	endif //__deprecated

// #	ifndef __error
// #		define __error(x)		__attribute__((error(x)))
// #	endif //__error

#	ifndef __warning
#		define __warning(x)		__attribute__((warning(x)))
#	endif //__warning

#	ifndef __always_inline
#		define __always_inline	__attribute__((always_inline))
#	endif //__always_inline

#	ifndef __flatten
#		define __flatten		__attribute__((flatten))
#	endif //__flatten

#	ifndef __hot
#		define __hot			__attribute__((hot))
#	endif //__hot

#	ifndef __cold
#		define __cold			__attribute__((cold))
#	endif //__cold

#	ifndef __fastcall
#		define __fastcall		__attribute__((fastcall))
#	endif //__fastcall

#	ifndef __thiscall
#		define __thiscall		__attribute__((thiscall))
#	endif //__thiscall

#	ifndef __malloc
#		define __malloc			__attribute__((malloc))
#	endif //__malloc

#	ifndef __pure
#		define __pure			__attribute__((pure))
#	endif //__pure

#	ifndef __const
#		define __const			__attribute__((const))
#	endif //__const

#   ifndef __likely
#       define __likely(x)      __builtin_expect (!!(x), 1)
#   endif //__likely

#   ifndef __unlikely
#       define __unlikely(x)    __builtin_expect (!!(x), 0)
#   endif //__unlikely

#	ifndef __fallthrough
#		define __fallthrough	__attribute__((fallthrough))
#	endif //__pure

#else //__GNUC__

#	ifndef __nowarn_unused
#		define __nowarn_unused
#	endif //__nowarn_unused

#	ifndef __deprecated
#		define __deprecated
#	endif //__deprecated

// #	ifndef __error
// #		define __error(x)
// #	endif //__error

#	ifndef __warning
#		define __warning(x)
#	endif //__warning

#	ifndef __always_inline
#		define __always_inline
#	endif //__always_inline

#	ifndef __flatten
#		define __flatten
#	endif //__flatten

#	ifndef __hot
#		define __hot
#	endif //__hot

#	ifndef __cold
#		define __cold
#	endif //__cold

#	ifndef __fastcall
#		define __fastcall
#	endif //__fastcall

#	ifndef __thiscall
#		define __thiscall
#	endif //__thiscall

#	ifndef __malloc
#		define __malloc
#	endif //__malloc

#	ifndef __pure
#		define __pure
#	endif //__pure

#	ifndef __const
#		define __const
#	endif //__const

#   ifndef __likely
#       define __likely(x)
#   endif //__likely

#   ifndef __unlikely
#       define __unlikely(x)
#   endif //__unlikely

#	ifndef __fallthrough
#		define __fallthrough
#	endif //__pure

#endif //__GNUC__


#endif //!__GNUC_EXTENSIONS_H__INCLUDED__

