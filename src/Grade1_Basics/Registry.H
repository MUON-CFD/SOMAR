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
#ifndef ___BaseRegistry_H__INCLUDED___
#define ___BaseRegistry_H__INCLUDED___

#include <map>
#include <sstream>
#include "LevelData.H"
#include "FluxBox.H"
#include "RealVect.H"


template<class Entry>
class BaseRegistry
{
protected:
    // Not for public use!
    friend class Registry;

    // Weak constructor
    BaseRegistry () {
    }

    // Destructor
    virtual
    ~BaseRegistry () {
        _clear();
    }

    virtual void
    _clear () {
        m_regMap.clear();
    }

    virtual int
    _set (Entry& a_entry) {
        // Create a unique key.
        int key;
        do {
            key = rand();
        } while (m_regMap.find(key) != m_regMap.end());

        // Put entry in map.
        m_regMap[key] = a_entry;

        // Return the new key to the user.
        return key;
    }

    virtual Entry&
    _get (int a_key) {
        // Find entry.
        typename RegMapType::iterator it = m_regMap.find(a_key);

        // If key was not found, abort.
        if (it == m_regMap.end()) {
            MayDay::Error("Key provided was not found in m_regMap.");
        }

        // Return the item.
        return it->second;
    }

    virtual const Entry&
    _get (int a_key) const {
        // Find entry.
        typename RegMapType::const_iterator it = m_regMap.find(a_key);

        // If key was not found, abort.
        if (it == m_regMap.end()) {
            MayDay::Error("Key provided was not found in m_regMap.");
        }

        // Return the item.
        return it->second;
    }


protected:
    // The registry.
    typedef std::map<int, Entry> RegMapType;
    RegMapType m_regMap;
};



// Example usage 1:
//  int a = 7;
//  int handle = myRegistryObj.set<int>(a);
//  ...
//  int b = myRegistryObj.get<int>(handle); // b is now 7.
//
// Example useage 2:
//  int a = 7;
//  myRegistryObj.set<int>(a, "My fancy identifier");
//  ...
//  int b = myRegistryObj.get<int>("My fancy identifier"); // b is now 7.
//
// Example usage 3:
//  void* aptr = (void*)(new int(8));
//  int h = m_oldStatePtr->set<void*>(aptr);
//  ...
//  int* bptr = m_oldStatePtr->getCast<int*>(h);
//  delete bptr;

// If you store a pointer, it's up to you to do the cleanup.
// So, you may want to store RefcountedPtrs instead.
//
// I know, the inheritance is a bit ridiculous.
class Registry
: public BaseRegistry<int>,
  public BaseRegistry<Real>,
  public BaseRegistry<IntVect>,
  public BaseRegistry<RealVect>,
  public BaseRegistry<Interval>,
  public BaseRegistry<std::string>,
  public BaseRegistry<ProblemDomain>,
  public BaseRegistry<BoxLayout>,
  public BaseRegistry<DisjointBoxLayout>,
  public BaseRegistry<LevelData<FArrayBox>*>,
  public BaseRegistry<RefCountedPtr<LevelData<FArrayBox> > >,
  public BaseRegistry<LevelData<FluxBox>*>,
  public BaseRegistry<RefCountedPtr<LevelData<FluxBox> > >,

  public BaseRegistry<Vector<int> >,
  public BaseRegistry<Vector<Real> >,
  public BaseRegistry<Vector<IntVect> >,
  public BaseRegistry<Vector<RealVect> >,
  public BaseRegistry<Vector<Interval> >,
  public BaseRegistry<Vector<std::string> >,
  public BaseRegistry<Vector<ProblemDomain> >,
  public BaseRegistry<Vector<BoxLayout> >,
  public BaseRegistry<Vector<DisjointBoxLayout> >,
  public BaseRegistry<Vector<LevelData<FArrayBox>*> >,
  public BaseRegistry<Vector<RefCountedPtr<LevelData<FArrayBox> > > >,
  public BaseRegistry<Vector<LevelData<FluxBox>*> >,
  public BaseRegistry<Vector<RefCountedPtr<LevelData<FluxBox> > > >
{
public:
    // Use this to create a new entry.
    template<class Entry>
    int
    set(Entry a_entry, std::string a_name = "") {
        int key = static_cast<BaseRegistry<Entry>* const>(this)->_set(a_entry);

        if (a_name.length() != 0) {
            if (m_phoneBook.find(a_name) != m_phoneBook.end()) {
                ostringstream msg;
                msg << a_name
                    << " is not a unique name. Cannot put in registry";
                MayDay::Error(msg.str().c_str());
            }
            m_phoneBook[a_name] = key;
        }

        return key;
    }

    // This will retrieve the key if you don't remember.
    int
    getKey(std::string a_name) const
    {
        std::map<std::string, int>::const_iterator it =
            m_phoneBook.find(a_name);

        if (it == m_phoneBook.end()) {
            ostringstream msg;
            msg << a_name << " cannot be found in registry";
            MayDay::Error(msg.str().c_str());
        }

        return it->second;
    }

    // Use this to retrieve a reference to an entry by key.
    template<class Entry>
    Entry&
    get(int a_key) {
        return BaseRegistry<Entry>::_get(a_key);
    }

    // Use this to retrieve a reference to an entry by name.
    template<class Entry>
    Entry&
    get(std::string a_name) {
        return BaseRegistry<Entry>::_get(getKey(a_name));
    }

    // Use this to retrieve a reference to an entry by key. (const version)
    template<class Entry>
    const Entry&
    get(int a_key) const {
        return BaseRegistry<Entry>::_get(a_key);
    }

    // Use this to retrieve a reference to an entry by name. (const version)
    template<class Entry>
    const Entry&
    get(std::string a_name) const {
        return BaseRegistry<Entry>::_get(getKey(a_name));
    }

    // Removes all entries from a single map.
    template<class Entry>
    void
    clear() {
        BaseRegistry<Entry>::_clear();
    }

    // Removes all entries from a all maps.
    virtual void
    clear() {
        BaseRegistry<int>::_clear();
        BaseRegistry<Real>::_clear();
        BaseRegistry<IntVect>::_clear();
        BaseRegistry<RealVect>::_clear();
        BaseRegistry<Interval>::_clear();
        BaseRegistry<std::string>::_clear();
        BaseRegistry<ProblemDomain>::_clear();
        BaseRegistry<BoxLayout>::_clear();
        BaseRegistry<DisjointBoxLayout>::_clear();
        BaseRegistry<LevelData<FArrayBox>*>::_clear();
        BaseRegistry<RefCountedPtr<LevelData<FArrayBox> > >::_clear();
        BaseRegistry<LevelData<FluxBox>*>::_clear();
        BaseRegistry<RefCountedPtr<LevelData<FluxBox> > >::_clear();

        BaseRegistry< Vector<int> >::_clear();
        BaseRegistry< Vector<Real> >::_clear();
        BaseRegistry< Vector<IntVect> >::_clear();
        BaseRegistry< Vector<RealVect> >::_clear();
        BaseRegistry< Vector<Interval> >::_clear();
        BaseRegistry< Vector<std::string> >::_clear();
        BaseRegistry< Vector<ProblemDomain> >::_clear();
        BaseRegistry< Vector<BoxLayout> >::_clear();
        BaseRegistry< Vector<DisjointBoxLayout> >::_clear();
        BaseRegistry< Vector<LevelData<FArrayBox>*> >::_clear();
        BaseRegistry< Vector<RefCountedPtr<LevelData<FArrayBox> > > >::_clear();
        BaseRegistry< Vector<LevelData<FluxBox>*> >::_clear();
        BaseRegistry< Vector<RefCountedPtr<LevelData<FluxBox> > > >::_clear();
    }

protected:
    // If you don't know the int needed to recover your data, you can look
    // it up here.
    // KIDS: A phone book was a 1000 page book delivered to your house. It was
    //  used to find people's phone numbers, hold up the low end of a couch,
    //  or tear in half at wrestling matches.
    std::map<std::string, int> m_phoneBook;
};


#endif //!___BaseRegistry_H__INCLUDED___
