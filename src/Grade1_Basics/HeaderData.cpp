#include "HeaderData.H"

// constructor that populates root attributes
HeaderData::HeaderData(const std::string& a_fileName, const int a_verb)
{
    for (const auto& value : rootIntV) m_int[value];
    for (const auto& value : rootRealV) m_real[value];
    for (const auto& value : rootRV) m_realvect[value];
    for (const auto& value : rootStringV) m_string[value];
    for (const auto& value : rootIV) m_intvect[value];
    std::string Root = "/";
    this->readFromFile(a_fileName, Root);

    // add extra components
    for (int d = 0; d < SpaceDim; ++d) {
        char D[20];
        sprintf(D, "%d", d);
        std::string label = std::string("base.isPeriodic_") + D;
        m_int[label];
    }
    for (int d = 0; d < m_int["rhs.eddyViscMethodSize"]; ++d) {
        char D[20];
        sprintf(D, "%d", d);
        std::string label = std::string("rhs.eddyViscMethod_") + D;
        m_int[label];
    }
    for (int d = 0; d < m_int["numScalars"]; ++d) {
        char D[20];
        sprintf(D, "%d", d);
        std::string label = std::string("scalarComponent_") + D;
        m_string[label];
    }
    for (int d = 0; d < m_int["rhs.scalarsKappaSize"]; ++d) {
        char D[20];
        sprintf(D, "%d", d);
        std::string label = std::string("rhs.scalarsKappa_") + D;
        m_real[label];
    }
    // reread
    this->readFromFile(a_fileName, Root);

    if (a_verb >= 1) {
        pout() << "hdf5 header(root) data: " << endl;
        for (const auto& value : rootIntV)
            pout() << value << " , " << m_int[value] << endl;
        for (const auto& value : rootRealV)
            pout() << value << " , " << m_real[value] << endl;
        for (const auto& value : rootRV)
            pout() << value << " , " << m_realvect[value] << endl;
        for (const auto& value : rootStringV)
            pout() << value << " , " << m_string[value] << endl;
        for (const auto& value : rootIV)
            pout() << value << " , " << m_intvect[value] << endl;
    }
}
// constructor that populates level_xxx attributes if a_rw!="r"
HeaderData::HeaderData(const std::string& a_fileName,
                       const int          level,
                       const int          a_verb)
{
    for (const auto& value : levelIntV) m_int[value];
    for (const auto& value : levelRealV) m_real[value];
    for (const auto& value : levelRV) m_realvect[value];
    for (const auto& value : levelStringV) m_string[value];
    for (const auto& value : levelIV) m_intvect[value];
    for (const auto& value : levelBox) m_box[value];
    char level_str[20];
    sprintf(level_str, "%d", level);
    const std::string label = std::string("level_") + level_str;
    this->readFromFile(a_fileName, label);
    if (a_verb >= 1) {
        pout() << "hdf5 " + label + " data: " << endl;
        for (const auto& value : rootIntV)
            pout() << value << " , " << m_int[value] << endl;
        for (const auto& value : rootRealV)
            pout() << value << " , " << m_real[value] << endl;
        for (const auto& value : rootRV)
            pout() << value << " , " << m_realvect[value] << endl;
        for (const auto& value : rootStringV)
            pout() << value << " , " << m_string[value] << endl;
        for (const auto& value : rootIV)
            pout() << value << " , " << m_intvect[value] << endl;
    }
}

int
HeaderData::writeToFile(const std::string& a_fileName,
                        const std::string& a_groupName) const
{
#ifdef CH_USE_PYTHON
  // integers
  for (const auto &value : m_int)
    Py::PythonFunction("IO", "WriteToHeader", a_fileName, a_groupName, value.first,
                          value.second);
  // Reals
  for (const auto &value : m_real)
    Py::PythonFunction("IO", "WriteToHeader", a_fileName, a_groupName, value.first,
                          value.second);
  // strings
  for (const auto &value : m_string)
    Py::PythonFunction("IO", "WriteToHeader", a_fileName, a_groupName, value.first,
                          value.second);
  // IntVect
  for (const auto &value : m_intvect)
    Py::PythonFunction("IO", "WriteToHeader", a_fileName, a_groupName, value.first,
                          value.second);
  // RealVect
  for (const auto &value : m_realvect)
    Py::PythonFunction("IO", "WriteToHeader", a_fileName, a_groupName, value.first,
                          value.second);
  // box
  for (const auto &value : m_box)
    Py::PythonFunction("IO", "WriteToHeader", a_fileName, a_groupName, value.first,
                          value.second);
#endif
    return 0;
}

int
HeaderData::readFromFile(const std::string& a_fileName,
                         const std::string& a_groupName)
{
#ifdef CH_USE_PYTHON
  Py::PythonFunction("IO", "OpenFileForRead", a_fileName);
  // integers
  for (auto &value : m_int)
    value.second = Py::PythonReturnFunction<int>("IO", "ReadFromHeader", a_fileName,
                                                    a_groupName, value.first);
  // Reals
  for (auto &value : m_real)
    value.second = Py::PythonReturnFunction<Real>("IO", "ReadFromHeader", a_fileName,
                                                     a_groupName, value.first);
  // strings
  for (auto &value : m_string)
    value.second = Py::PythonReturnFunction<std::string>(
							 "IO", "ReadFromHeader", a_fileName, a_groupName, value.first);
  // RealVect
  for (auto &value : m_realvect)
    value.second = Py::PythonReturnFunction<RealVect>(
						      "IO", "ReadFromHeader", a_fileName, a_groupName, value.first);
  // IntVect
  for (auto &value : m_intvect)
    value.second = Py::PythonReturnFunction<IntVect>(
						     "IO", "ReadFromHeader", a_fileName, a_groupName, value.first);

  // Box
  for (auto &value : m_box)
    value.second = Py::PythonReturnFunction<Box>("IO", "ReadFromHeader", a_fileName,
                                                    a_groupName, value.first);

  Py::PythonFunction("IO", "CloseFile", a_fileName);
#endif
    return 0;
}

ostream&
operator<<(ostream& out, const HeaderData& header)
{
    // integers
    for (const auto& value : header.m_int)
        out << value.first << ", " << value.second << endl;

    // Reals
    for (const auto& value : header.m_real)
        out << value.first << ", " << value.second << endl;

    // strings
    for (const auto& value : header.m_string)
        out << value.first << ", " << value.second << endl;

    // IntVect
    for (const auto& value : header.m_intvect)
        out << value.first << ", " << value.second << endl;

    // RealVect
    for (const auto& value : header.m_realvect)
        out << value.first << ", " << value.second << endl;

    // box
    for (const auto& value : header.m_box)
        out << value.first << ", " << value.second << endl;
    return out;
}


