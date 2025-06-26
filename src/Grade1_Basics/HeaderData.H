#pragma once

#include <Box.H>
#include <IntVect.H>
#include <REAL.H>
#include <RealVect.H>
#include <iostream>
#include <map>
#include <ostream>
#include <string>
//#include "IO.H"

#ifdef CH_USE_PYTHON
#include "PyGlue.H"
#endif



class HeaderData
{
public:
    ///
    /**
     * simple constructor, does not populate the fields, normally used to write
     * to the header file.
     */
    HeaderData(){};
    /**
     * This constructor populates the map files with the keys from Root
     */
    HeaderData(const std::string& a_filename, const int a_verb);
    ///
    ///
    /**
     * This constructor populates the map files with the keys from Root/level
     */
    HeaderData(const std::string& a_filename,
               const int          a_level,
               const int          a_verb);
    ///
    /**
     * destructor invokes clear on the data member maps
     */
    ~HeaderData()
    {
        this->clear();
    }
    ///
    /**
       Writes this HDF5HeaderData's current attribute list to the
       current group in 'file.'  Returns 0 on success, returns the
       error code from HDF5 on failure.

   */
    int
    writeToFile(const std::string& a_fileName,
                const std::string& a_groupName) const;




  ///
  /**
     Reads into this HDF5HeaderData's attribute list from file.  Read
     process is add/change, does not remove key-value pairs. Reads
     from current group.  Returns 0 on success, positive number if a
     particular member of group caused an error, negative on general
     error.

 */
  int readFromFile(const std::string &a_fileName,
                   const std::string &a_groupName);


  ///
  void clear() {
    m_real.clear();
    m_int.clear();
    m_string.clear();
    m_intvect.clear();
    m_realvect.clear();
    m_box.clear();
  }

  ///
  std::map<std::string, Real> m_real;

  ///
  std::map<std::string, int> m_int;

  ///
  std::map<std::string, std::string> m_string;

  ///
  std::map<std::string, IntVect> m_intvect;

  ///
  std::map<std::string, Box> m_box;

  ///
  std::map<std::string, RealVect> m_realvect;

  /// This dump the content of the attributes etc. (basically, everything but data)
  void dump() const {
  #ifdef CH_USE_PYTHON
  Py::PythonFunction("IO", "dumpFile");
  #endif
  }

  /// This returns the same information that dump() outputs into a string.
  std::string headerAsString() const {
  #ifdef CH_USE_PYTHON
    return Py::PythonReturnFunction<std::string>("IO", "GetHeaderString");
  #else
    std::string s="";
    return s;
 #endif
  }

  friend ostream & operator << (ostream &out, const HeaderData &header);

 private:
  std::vector<std::string> rootIV{"base.nx", "base.nxOffset"};
  std::vector<std::string> rootRV{"base.L", "rhs.coriolisF"};
  std::vector<std::string> rootIntV{
      "base.blockFactor",   "amr.maxLevel",
      "amr.useSubcycling", "iteration",
      "max_level",         "numQComps",
      "num_levels",        "numScalars",
      "regrid_interval_0", "rhs.eddyViscMethodSize",
      "rhs.eddyPrandtlScalarsSize", "rhs.scalarsKappaSize"};
  std::vector<std::string> rootRealV{ "time.dtMult",       "time.maxDt",
                                      "time.fixedDt",      "rhs.TKappa",
                                      "rhs.SKappa",
                                      "rhs.eddyPrandtlT", "rhs.eddyPrandtlS",
                                      "rhs.nu",           "time" };
  std::vector<std::string> rootStringV;

  std::vector<std::string> levelBox{"domainBox"};
  std::vector<std::string> levelIV{"ref_ratio"};
  std::vector<std::string> levelRV{"vec_dx"};
  std::vector<std::string> levelRealV{"dt", "time"};
  std::vector<std::string> levelIntV{"finestExtantLevel", "isEmpty"};
  std::vector<std::string> levelStringV{};


public:
// -----------------------------------------------------------------------------
// Simple reading/writing utilities
template<typename T>
static void writeHeaderValue (const std::string& a_fileName,
                 const std::string& a_groupName,
                 const std::string& a_varName,
                 const T         a_data)
{
  #ifdef CH_USE_PYTHON
    Py::PythonFunction("IO", "WriteToHeader", a_fileName, a_groupName, a_varName,
                          a_data);
  #endif
}




template<typename T>
static T readHeaderValue (const std::string& a_fileName,
                const std::string& a_groupName,
                const std::string& a_varName)
{
  #ifdef CH_USE_PYTHON
  Py::PythonFunction("IO", "OpenFileForRead", a_fileName);
  T X = Py::PythonReturnFunction<T>("IO", "ReadFromHeader", a_fileName, a_groupName, a_varName);
  Py::PythonFunction("IO", "CloseFile", a_fileName);
  return X;
#endif
}
};
