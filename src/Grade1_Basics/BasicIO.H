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
 *  https://github.com/MUON-CFD/somar.
 ******************************************************************************/
#ifndef ___BasicIO_H__INCLUDED___
#define ___BasicIO_H__INCLUDED___

#include <iostream>


class BasicIO
{
public:
    // Print to terminal from one or all procs
    static std::ostream&
    tout (int a_procid = -1);

protected:
    // -------------------------------------------------------------------------
    // Stream to nowhere...
    // This code came from:
    // https://stackoverflow.com/questions/19200207/redirect-debug-output-to-null-stream-instead-of-stdcerr
    // -------------------------------------------------------------------------
    template <class cT, class traits = std::char_traits<cT> >
    class basic_nullbuf: public std::basic_streambuf<cT, traits> {
        inline typename traits::int_type
        overflow (typename traits::int_type c)
        {
            return traits::not_eof(c); // indicate success
        }
    };

    template <class cT, class traits = std::char_traits<cT> >
    class basic_onullstream: public std::basic_ostream<cT, traits> {
    public:
        inline basic_onullstream ():
        std::basic_ios<cT, traits>(&m_sbuf),
        std::basic_ostream<cT, traits>(&m_sbuf)
        {
            // note: the original code is missing the required this->
            this->init(&m_sbuf);
        }

    private:
        basic_nullbuf<cT, traits> m_sbuf;
    };

    static basic_onullstream<char> s_null_ostream;
};


#endif //!___BasicIO_H__INCLUDED___
