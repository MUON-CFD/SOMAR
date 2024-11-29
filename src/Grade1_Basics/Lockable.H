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
#ifndef ___Lockable_H__INCLUDED___
#define ___Lockable_H__INCLUDED___

#include "Debug.H"

// Typical usage:
// class Foo: public Lockable {
//   Foo() {};
//
//   do(const Key& a_key) {
//     CH_verify(this->accessGranted(a_key));
//   }
// };
//
// int main () {
//   Foo foo;
//   Foo::Key key = foo.lock();
//   foo.do(key);
//   foo.unlock(key);
// }


class Lockable
{
public:
    class Key final
    {
        friend class Lockable;
    private:
        Key(): value(nullptr) {}
        int* value;
    };

    Lockable()
    : m_isLocked(false)
    , m_key()
    {}

    // Copying m_key.value would be bad!
    Lockable& operator=(const Lockable&) = delete;
    Lockable(const Lockable&) = delete;

    virtual
    ~Lockable()
    {
        if (m_isLocked) {
            MAYDAYERROR("Unlock object before calling destructor.");
        }
    }

    inline virtual bool
    isLocked() const final
    {
        return m_isLocked;
    }

    inline virtual bool
    isUnlocked() const final
    {
        return !m_isLocked;
    }

    inline virtual bool
    accessGranted(const Key& a_key) const final
    {
        return (a_key.value == m_key.value);
    }

    inline virtual const Key&
    lock() final
    {
        if (m_isLocked) {
            MAYDAYERROR("Object is already locked.");
        }
        m_isLocked = true;
        m_key.value = new int;
        return m_key;
    }

    inline virtual void
    unlock(const Key& a_key) final
    {
        if (this->accessGranted(a_key)) {
            delete m_key.value;
            m_key.value = nullptr;
            m_isLocked = false;
        } else {
            MAYDAYERROR("Cannot unlock. Wrong key provided.");
        }
    }

private:
    bool m_isLocked;
    Key  m_key;
};


#endif //!___Lockable_H__INCLUDED___
