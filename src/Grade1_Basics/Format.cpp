#include "Format.H"
#include <iomanip>
#include <sstream>
namespace Format {

// Indenter...

// #define TURN_INDENTER_OFF


// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
Indenter::Indenter (std::streambuf* dest,
                    const int       indent,
                    const char*     bullet)
: m_dest(dest),
  m_isAtStartOfLine(true),
  m_indent(indent, ' '),
  m_owner(NULL)
{
    m_indent += bullet;
}


// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
Indenter::Indenter (std::ostream& dest,
                    const int     indent,
                    const char*   bullet)
: m_dest(dest.rdbuf()),
  m_isAtStartOfLine(true),
  m_indent(indent, ' '),
  m_owner(&dest)
{
#ifndef TURN_INDENTER_OFF
    if (m_dest != NULL) {
        m_indent += bullet;
        m_owner->rdbuf(this);
    }
#endif
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
Indenter::~Indenter ()
{
#ifndef TURN_INDENTER_OFF
    if (m_dest != NULL && m_owner != NULL) {
        m_owner->rdbuf(m_dest);
    }
#endif
}


// -----------------------------------------------------------------------------
// Overflow override
// -----------------------------------------------------------------------------
int Indenter::overflow (int ch)
{
    if (m_isAtStartOfLine && ch != '\n') {
        m_dest->sputn(m_indent.data(), m_indent.size());
    }
    m_isAtStartOfLine = ch == '\n';
    return m_dest->sputc(ch);
}


// indent/unindent...

// -----------------------------------------------------------------------------
// Static variables.
// -----------------------------------------------------------------------------
std::stack<RefCountedPtr<Indenter> > indent::s_indentStack;


// -----------------------------------------------------------------------------
// Indent() stream manipulator. All stream data will be indented or bulleted
// until unindent is used.
// -----------------------------------------------------------------------------
void indent::operator() (std::ostream& os) const
{
    if (os.rdbuf() != 0) {
        s_indentStack.push(
            RefCountedPtr<Indenter>(new Indenter(os, m_indent, m_bullet))
        );
    }
}


// -----------------------------------------------------------------------------
// Removes the indentation / restores the stream.
// -----------------------------------------------------------------------------
std::ostream& unindent (std::ostream& os)
{
    if (os.rdbuf() != 0 && !indent::s_indentStack.empty()) {
        indent::s_indentStack.pop();
    }
    return os;
}


// -----------------------------------------------------------------------------
// Removes the indentation / restores the stream.
// -----------------------------------------------------------------------------
std::ostream& unindentAll (std::ostream& os)
{
    while (os.rdbuf() != 0 && !indent::s_indentStack.empty()) {
        indent::s_indentStack.pop();
    }
    return os;
}


// -----------------------------------------------------------------------------
// Allows class indent to be used in a stream.
// -----------------------------------------------------------------------------
std::ostream& operator<< (std::ostream& os, const indent& idt)
{
    idt(os);
    return os;
}


// fixed/scientific...

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
std::ostream& fixed (std::ostream& os)
{
    return (os << std::fixed);
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
std::ostream& scientific (std::ostream& os)
{
    return (os << std::scientific << std::setprecision(6));
}


// pushFlags/popFlags...

// -----------------------------------------------------------------------------
// Private, static stack for pushFlags/popFlags.
// -----------------------------------------------------------------------------
struct ___flagsData {
    friend std::ostream& pushFlags (std::ostream& os);
    friend std::ostream& popFlags (std::ostream& os);
private:
    static std::stack<std::ios_base::fmtflags> s_flagsStack;
};

std::stack<std::ios_base::fmtflags> ___flagsData::s_flagsStack;


// -----------------------------------------------------------------------------
// Stores the flags.
// -----------------------------------------------------------------------------
std::ostream& pushFlags (std::ostream& os)
{
    ___flagsData::s_flagsStack.push(os.flags());
    return os;
}


// -----------------------------------------------------------------------------
// Restores the flags.
// -----------------------------------------------------------------------------
std::ostream& popFlags (std::ostream& os)
{
    os.flags(___flagsData::s_flagsStack.top());
    ___flagsData::s_flagsStack.pop();
    return os;
}


// banner...

// -----------------------------------------------------------------------------
// Creates a section heading for use by a stream (typically, cout or pout).
// -----------------------------------------------------------------------------
std::string
banner(const std::string& a_str)
{
    const char border = '#';

    if (a_str.length() == 0) {
        return "\n" + std::string(80, border) + "\n";
    }
    int leftdashlen  = (80 - a_str.length()) / 2 - 1;
    int rightdashlen = 78 - leftdashlen - a_str.length();
    return "\n" + std::string(leftdashlen, border) + " " + a_str + " " +
           std::string(rightdashlen, border) + "\n";
}


// textTime...

// -----------------------------------------------------------------------------
// Converts a time in seconds (double) into the more readable days, hours, mins,
// and seconds (std::string).
// -----------------------------------------------------------------------------
std::string textTime (const double inSeconds)
{
    std::ostringstream retStr;
    retStr << fixed;

    if (inSeconds < 60) {
        retStr << inSeconds << " secs";
    } else {
        const long longSecs = (long)inSeconds;
        const int seconds = longSecs % 60;
        const int minutes = (longSecs / 60) % 60;
        const int hours = (longSecs / (60 * 60)) % 24;
        const int days = longSecs / (24 * 60 * 60);

        if (days > 0)
            retStr << days << " days, ";
        if (days > 0 || hours > 0)
            retStr << hours << " hrs, ";
        if (days > 0 || hours > 0 || minutes > 0)
            retStr << minutes << " mins, ";
        retStr << seconds << " secs";
    }

    return retStr.str();
}


// Warning...

// -----------------------------------------------------------------------------
// Static variables.
// -----------------------------------------------------------------------------
const char* const Warning::s_prefix = "\n   ";


// -----------------------------------------------------------------------------
// Writes a formatted warning message to an ostream.
// -----------------------------------------------------------------------------
void Warning::operator() (std::ostream& os) const
{
    os << pushFlags << fixed
       << std::endl
       << yellow << "WARNING: "
       << hiwhite << m_msg
       << none << std::flush;

    if (m_file[0] != 0) {
        os << s_prefix << "File: " << m_file;
    }

    if (m_line != 0) {
        os << s_prefix << "Line: " << m_line;
    }

    os << s_prefix << "Note: " << popFlags;
}


// -----------------------------------------------------------------------------
// Allows Warning to be used in an ostream.
// -----------------------------------------------------------------------------
std::ostream& operator<< (std::ostream& os, const Warning& warning) {
    warning(os);
    return os;
}


// Error...

const char* const Error::s_prefix = "\n ";

// -----------------------------------------------------------------------------
// Writes a formatted error message to an ostream.
// -----------------------------------------------------------------------------
void Error::operator() (std::ostream& os) const
{
    os << pushFlags << fixed
       << std::endl
       << hired << "\nERROR: "
       << hiwhite << m_msg
       << none << std::flush;

    if (m_file[0] != 0) {
        os << s_prefix << "File: " << m_file;
    }

    if (m_line != 0) {
        os << s_prefix << "Line: " << m_line;
    }

    os << s_prefix << "Note: " << popFlags;
}

// -----------------------------------------------------------------------------
// Allows Error to be used in an ostream.
// -----------------------------------------------------------------------------
std::ostream& operator<< (std::ostream& os, const Error& error) {
    error(os);
    return os;
}



// -----------------------------------------------------------------------------
// Used to make TODO() output readable.
// -----------------------------------------------------------------------------
std::string stripFuncName(const std::string& a_func) {
    // std::string::size_type  start, end, center;
    // center = a_func.find_last_of("::");
    // start = a_func.rfind(" ", center-1) + 1;
    // end   = a_func.find_first_of("(", center+1) - 1;
    // return a_func.substr(start, end - start + 1) + "()";
    return a_func;
}


// -----------------------------------------------------------------------------
// Used to make TODO() output readable.
// -----------------------------------------------------------------------------
std::string stripFileName(const std::string& a_file) {
    std::string::size_type start = a_file.find_last_of("/") + 1;
    if(start != std::string::npos && start < a_file.length())
        return a_file.substr(start, a_file.length() - start);
    return a_file;
}


} // namespace Format


