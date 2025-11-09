/*
 *  Copyright (C) 2015 TU Delft. All rights reserved.
 *
 *
 *  Penalty contact method, only for small sliding displacements
 *  
 *  Author: Yaolu Liu
 *  Date: August 2016
 *
 */


#ifndef CONTACT_COHESIVE_MATERIAL_H
#define CONTACT_COHESIVE_MATERIAL_H

#include <jem/base/String.h>
#include <jem/util/Flex.h>

#include "CohesiveMat.h"

using jem::idx_t;
using jem::String;
using jem::util::Flex;


// =======================================================
//  class ContactCohesiveMat
// =======================================================


class ContactCohesiveMat : public virtual CohesiveMat
{
 public:


  static const char*      DUMMY_PROP;
  

  explicit                ContactCohesiveMat

    ( const idx_t           rank,
      const Properties&     globdat );
		


  /*
   *  configure and getConfig (from input file)
   */

  virtual void            configure

    ( const Properties&     props,
      const Properties&     globdat );

  virtual void            getConfig

    ( const Properties&     conf,
      const Properties&     globdat )         const;

  /*
   *  compute the traction (t) and cohesive tangent stiffness matrix
   *  (stiff) at material point mpoint given the displacement jump (jump)
   *   jump[0] = crack opening displacement
   *   jump[1] = crack sliding displacement
   */

  virtual void            update

    ( Vector&               traction,
      Matrix&               stiff,
      const Vector&         jump,
      idx_t                 mpoint );
      
      
  virtual void            elasticUpdate

    ( Vector&               traction,
      Matrix&               stiff,
      const Vector&         jump );

  /*
   *  Called when the Newton Raphson converged
   *  Swap newHist_ to oldHist_ 
   */

  virtual void            commit           ();

  /**
   *  Allocate for history variables if it is not yet
   *  initialized. Otherwise, extend it by appending at the end.
   *  Initial value of damage is optional
   */

  virtual void            allocPoints

    ( const idx_t           count );

  virtual void            deallocPoints

    ( idx_t                 count );		

  virtual double          giveHistory       

    ( idx_t                 point  ) const;

 protected:

  virtual                ~ContactCohesiveMat   ();
  
  void                    allocPoints_

    ( const idx_t           count,
      const idx_t           rank );

 protected:

  // history variable (equivalent crack opening), involve in time

  	class                   Hist_
  	{
   		public:
                            Hist_(const idx_t rank);

    	Vector                  jumpi;
    	Vector                  tractioni;
  	};
	
  	Flex<Hist_>             preHist_;      // history of the previous load step
  	Flex<Hist_>             newHist_;      // history of the current iteration
  	Flex<Hist_>*            latestHist_;   // points to latest hist

 // variables...

     double                	dummy_; 
		 Matrix									PenaltyStiff_;
		 double									minTol_;
};




#endif 
