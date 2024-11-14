#ifndef Me949b7e6032829e355cb3ebad4d68ff6
#   error LayoutToolsI.H is not for public use! Try including LayoutTools.H instead.
#endif


// *****************************************************************************
// The VertShifter class
//
// This shifts shiftFAB so that its lower vertical end coincides with the lower
// vertical end of target*. The centering will also be changed. You MUST call
// restore() before this object goes out of scope or an error will be thrown.
// *****************************************************************************

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template<typename T>
VertShifter<T>::VertShifter (DataType&       a_shiftFAB,
                             const DataType& a_targetFAB)
: m_dataPtr(NULL)
{
    this->shift(a_shiftFAB, a_targetFAB.box());
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template<typename T>
VertShifter<T>::VertShifter (DataType&  a_shiftFAB,
                             const Box& a_targetBox)
: m_dataPtr(NULL)
{
    this->shift(a_shiftFAB, a_targetBox);
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template<typename T>
VertShifter<T>::~VertShifter ()
{
#ifndef NDEBUG
    if (m_dataPtr != NULL) {
        MayDay::Error("You forgot to call VertShifter<T>::restore");
    }
#endif
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template<typename T>
void
VertShifter<T>::restore ()
{
    CH_assert(m_dataPtr != NULL);
    m_dataPtr->shift(SpaceDim-1, -m_totalShift);
    if (m_typeIsChanged) {
        m_dataPtr->shiftHalf(SpaceDim-1, -1);
    }
    m_dataPtr = NULL;
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template<typename T>
void
VertShifter<T>::shift (DataType&  a_shiftFAB,
                       const Box& a_targetBox)
{
    CH_assert(m_dataPtr == NULL);
    m_dataPtr = &a_shiftFAB;

    // Change type
    if (a_targetBox.type(SpaceDim-1) != a_shiftFAB.box().type(SpaceDim-1)) {
        a_shiftFAB.shiftHalf(SpaceDim-1, 1);
        m_typeIsChanged = true;
    } else {
        m_typeIsChanged = false;
    }

    // Shift
    m_totalShift = a_targetBox.smallEnd(SpaceDim-1)
                 - a_shiftFAB.box().smallEnd(SpaceDim-1);
    a_shiftFAB.shift(SpaceDim-1, m_totalShift);
}
