/*******************************************************************************
 *    FILE: AnisotropicAMRLevelFactory.H
 *    DESC: Almost identical to Chombo's AMRLevelFactory.
 *    DATE: Wed 14 Aug 2013 07:07:25 PM EDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __AnisotropicAMRLevelFactory_H__INCLUDED__
#define __AnisotropicAMRLevelFactory_H__INCLUDED__

#include "AnisotropicAMRLevel.H"

/// Factory class to produce a AnisotropicAMRLevel.
/**
   This is needed for construction of AnisotropicAMRLevel-derived classes by
   the AMR class.
*/
class AnisotropicAMRLevelFactory
{
public:
  ///
  /**
      Factory function.  Creates a new AnisotropicAMRLevel.
      In your instantiation you should include defining
      all the stuff that is not part of the AnisotropicAMRLevel API.
  */
  virtual AnisotropicAMRLevel* new_amrlevel() const = 0;

  // to make gcc v4 happy
  virtual ~AnisotropicAMRLevelFactory()
  {
  }
};


#endif //!__AnisotropicAMRLevelFactory_H__INCLUDED__
