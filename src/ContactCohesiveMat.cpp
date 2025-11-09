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

#include <jem/base/System.h>
#include <jem/base/Error.h>
#include <jem/base/array/utilities.h>
#include <jem/base/array/select.h>
#include <jem/base/array/operators.h>
#include <jem/base/array/tensor.h>
#include <jem/base/array/intrinsics.h>
#include <jem/io/PrintWriter.h>
#include <jem/io/FileWriter.h>
#include <jem/util/Properties.h>
#include <jem/numeric/utilities.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/MatmulChain.h>
#include <jem/numeric/algebra/LUSolver.h>
#include <iostream>
#include <cmath>


#include <jive/model/Actions.h>

#include "ContactCohesiveMat.h"

using namespace jem;
using namespace jem::io;

using jem::numeric::matmul;
using jem::numeric::norm2;
using jem::numeric::MatmulChain;
using jem::numeric::LUSolver;
using jive::model::Actions;


//=======================================================================
//   class ContactCohesiveMat
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char*  ContactCohesiveMat::DUMMY_PROP      = "dummy";

//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------

ContactCohesiveMat::ContactCohesiveMat

  ( const idx_t        rank,
    const Properties&  globdat )

  : CohesiveMat ( rank, globdat )
  
{
  JEM_PRECHECK ( rank >= 1 && rank <= 3 );
}


ContactCohesiveMat::~ContactCohesiveMat ()
{}



ContactCohesiveMat::Hist_::Hist_ (const idx_t rank )
{

    jumpi.resize ( rank );
    jumpi = 0.;
    tractioni.resize ( rank );
    tractioni = 0.;
    
}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void ContactCohesiveMat::configure

( const Properties&     props,
 const Properties&     globdat )

{
    
    using jem::maxOf;
    
    minTol_ = -20.e-3;
    
    props.get (dummy_, DUMMY_PROP, 0., maxOf( dummy_ ) );
    props.find (minTol_, "contactTol" );    
	

    
    
    PenaltyStiff_.resize ( rank_ , rank_ );
    PenaltyStiff_ = 0.0;
    
    // only the normal direction has stiffness
		PenaltyStiff_( 0, 0 ) = dummy_ ;
		
    
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void ContactCohesiveMat::getConfig

( const Properties& conf ,
 const Properties& globdat ) const

{
    
    conf.set ( DUMMY_PROP, dummy_ );
    conf.set ( "ContactTol", minTol_ );
        
}




//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void  ContactCohesiveMat::update

(		Vector&               traction,
		Matrix&               stiff,
		const Vector&         jump,
		idx_t                 mpoint )

{
		// update traction and stiffness

		if ( jump[0] > minTol_ )
		{
			traction = 0.;
			stiff = 0.;
		}	
		else
		{	
			Vector  gapdisp( rank_ );
			gapdisp = 0.;		
			gapdisp[0] = minTol_;

		
		
			traction = matmul( PenaltyStiff_ , jump ) - matmul( PenaltyStiff_, gapdisp );
			stiff = PenaltyStiff_ ;
		}
		
		
		
		// update history variables	
		
		newHist_[mpoint].jumpi = jump ;
		newHist_[mpoint].tractioni = traction ;
		
		latestHist_ = &newHist_;
}



//-----------------------------------------------------------------------
//   elasticUpdate (regular)
//-----------------------------------------------------------------------

void  ContactCohesiveMat::elasticUpdate

( Vector&               traction,
 Matrix&               stiff,
 const Vector&         jump )
{
    TensorIndex i,j;
    
    traction = dummy_ * jump;
    
    stiff(i,j) = where ( i == j, dummy_, 0. );
}




// --------------------------------------------------------------------
//  commit
// --------------------------------------------------------------------

void  ContactCohesiveMat::commit()

{
    newHist_.swap ( preHist_ );    
    latestHist_ = &preHist_;
}


// --------------------------------------------------------------------
//  allocPoints
// --------------------------------------------------------------------

void  ContactCohesiveMat::allocPoints 

    ( const idx_t  count )

{
  	allocPoints_ ( count, rank_ );
}

// --------------------------------------------------------------------
//  allocPoints
// --------------------------------------------------------------------

void  ContactCohesiveMat::allocPoints_

    ( const idx_t   count,
      const idx_t   rank)

{
    System::out() << "I have already allocated points" << endl;

    for ( idx_t i  = 0; i < count; ++i )
    {
        preHist_.pushBack ( Hist_(rank) );
        newHist_.pushBack ( Hist_(rank) );
    }
  
}

// --------------------------------------------------------------------
//  deallocPoints
// --------------------------------------------------------------------

void ContactCohesiveMat::deallocPoints ( idx_t count )
{
  	preHist_.popBack ( count );
  	newHist_.popBack ( count );
}

// --------------------------------------------------------------------
//  giveHistory
// --------------------------------------------------------------------

double ContactCohesiveMat::giveHistory ( idx_t point ) const
{
  return 0.;
}


