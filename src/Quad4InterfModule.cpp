/*
 *
 *  This module reorders quadrilateral elements for usage 
 *  as shell-shell line interface elements.
 *
 *  Sergio Gustavo Ferreira Cordeiro, July 2025
 *  
 */

#include <jem/base/array/utilities.h>
#include <jem/base/array/operators.h>
#include <jem/base/array/select.h>
#include <jem/base/Error.h>
#include <jem/base/System.h>
#include <jem/io/PrintWriter.h>
#include <jem/numeric/utilities.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/util/SparseArray.h>
#include <jem/util/ArrayBuffer.h>
#include <jive/app/ModuleFactory.h>
#include <jive/fem/NodeGroup.h>

#include "Quad4InterfModule.h"

using jem::System;
using jem::Error;
using jem::ALL;
using jem::Tuple;
using jem::Slice;
using jem::io::endl;
using jem::io::PrintWriter;
using jem::numeric::norm2;
using jive::fem::ElementSet;
using jive::fem::toXElementSet;

//=======================================================================
//   class Quad4InterfModule
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char* Quad4InterfModule::TYPE_NAME       = "Quad4Interf";
const char* Quad4InterfModule::ELEM_GROUPS     = "elemGroups";
const char* Quad4InterfModule::THICK_PROP       = "thickness";

//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------

Quad4InterfModule::Quad4InterfModule

  ( const String&  name ) :

      Super   ( name   )

{}

Quad4InterfModule::~Quad4InterfModule ()
{}


//-----------------------------------------------------------------------
//   init
//-----------------------------------------------------------------------

Module::Status Quad4InterfModule::init

  ( const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  String      contxt  = getContext();
  Properties  myConf  = conf .makeProps ( myName_ );
  Properties  myProps = props.findProps ( myName_ );

  myProps.get ( elemGroups_, ELEM_GROUPS );
  myConf.set  ( ELEM_GROUPS, elemGroups_ );

  elems_ = toXElementSet ( ElementSet::get ( globdat, getContext() ) );
  nodes_ = NodeSet::get ( globdat, contxt );

  rank_ = nodes_.rank();

  thickness_ = 0.;

  myProps.find( thickness_, THICK_PROP );
  myConf.set  ( THICK_PROP, thickness_ );

  JEM_PRECHECK ( rank_ == 3 );

  for ( idx_t i = 0; i < elemGroups_.size(); ++i )
  {
    reorderElements_ ( elemGroups_[i], globdat );
  }

  return DONE;
}

//-----------------------------------------------------------------------
//   getThickness_
//-----------------------------------------------------------------------

void Quad4InterfModule::getThickness_

  ( double thickness ) const 

{
  thickness = thickness_;
}

//-----------------------------------------------------------------------
//   reorderElements_
//-----------------------------------------------------------------------

void Quad4InterfModule::reorderElements_

  ( const String&     name,
    const Properties& globdat )

{
  Ref<PrintWriter> out = newInstance<PrintWriter>
    ( &System::debug("quad4interf") );

  String    contxt = getContext();
  ElemGroup egroup = ElemGroup::get ( name, elems_, globdat, contxt );
  IdxVector ielems = egroup.getIndices();
  IdxVector iperm;

  // find permutation 
  // (assuming all elements in the ElemGroup have the same ordering)

  idx_t nodeCount = findPermutation_ ( iperm, ielems );
  counterClockwise_ ( iperm, ielems );

  // reorder elements

  IdxVector inodes ( elems_.maxElemNodeCountOf ( ielems ) );

  for ( idx_t ie = 0; ie < ielems.size(); ++ie )
  {
    idx_t ielem = ielems[ie];

    elems_.getElemNodes ( inodes, ielem );

    IdxVector permutated ( nodeCount );

    permutated = inodes[iperm];

    *out << inodes << " --> " << permutated << endl;

    elems_.setElemNodes ( ielem, permutated );
  }
}

//-----------------------------------------------------------------------
//   counterClockwise_
//-----------------------------------------------------------------------

void Quad4InterfModule::counterClockwise_

  ( const IdxVector&  iperm,
    const IdxVector&  ielems ) const

{
  Ref<PrintWriter> out = newInstance<PrintWriter>
    ( &System::debug("quad4interf") );

  // after permutation, check wether element that shares first two nodes 
  // neighboring element is on the right hand side
  // if not: update permutation to make numbering counterclockwise

  idx_t maxNodeCount  = elems_.maxElemNodeCount ();
  idx_t nodeCount = iperm.size();

  IdxVector  inodes     ( maxNodeCount );
  IdxVector  jnodes     ( maxNodeCount );
  IdxVector  permutated (    nodeCount );
  
  //  3-------2       2-------3
  //  |       |  -->  |       |
  //  1-------0       0-------1   etc

  Matrix     coordsi ( rank_, nodeCount    );
  Matrix     coordsj ( rank_, maxNodeCount );

  elems_.getElemNodes ( inodes, ielems[0] );

  permutated = inodes[iperm];

  nodes_.getSomeCoords ( coordsi, permutated );

  idx_t i0 = permutated[0];
  idx_t i1 = permutated[1];

  Tuple<double,2>  p0 ( coordsi(0,0), coordsi(1,0) );
  Tuple<double,2>  p1 ( coordsi(0,1), coordsi(1,1) );

  idx_t clockwise = -1;
  for ( idx_t ie = 0; ie < elems_.size(); ++ie )
  {
    idx_t nn = elems_.getElemNodes ( jnodes, ie );

    Slice sl ( 0, nn );

    if ( testany ( jnodes[sl] == i0 ) &&
         testany ( jnodes[sl] == i1 ) &&
        !testany ( ie == ielems ) )
    {
      nodes_.getSomeCoords ( coordsj, jnodes[sl] );

      double X = sum(coordsj(0,sl)) / nn;
      double Y = sum(coordsj(1,sl)) / nn;

      double dist = (p1[0]-p0[0])*(Y-p0[1]) - (p1[1]-p0[1])*(X-p0[0]);

      *out << p0 << p1 << Tuple<double,2>(X,Y) << dist << endl;

      clockwise = dist > 0.;

      break;
    }
  }

  *out << "clockwise: " << clockwise << endl;

  if ( clockwise == -1 )
  {
    System::out() << "Something went wrong in Quad4Interf::counterClockwise_"
      << "\n continuing assuming the order was already counterclockwise";
  }
  else if ( clockwise == 1 )
  {
    IdxVector jperm = iperm.clone();

    for ( idx_t j = 0; j < nodeCount/2; ++j )
    {
      iperm[j] = jperm[nodeCount/2-j-1];
      iperm[j+nodeCount/2] = jperm[nodeCount-j-1];
    }

    *out << "correcting order: " << iperm << " --> " << jperm << endl;
  }
}

//-----------------------------------------------------------------------
//   findPermutation
//-----------------------------------------------------------------------

idx_t Quad4InterfModule::findPermutation_

  (       IdxVector&  iperm,
    const IdxVector&  ielems ) const

{
  idx_t nodeCount  = elems_.maxElemNodeCountOf ( ielems );

  IdxVector inodes ( nodeCount );
  Matrix    coords ( rank_, nodeCount );

  elems_.getElemNodes ( inodes, ielems[0] );
  nodes_.getSomeCoords ( coords, inodes );

  // find node that is closest to node 0

  double tol = 1.e-5;
  idx_t  imin = -1;

  for ( idx_t in = 1; in < nodeCount; ++in )
  {
    double dist = norm2 ( coords(ALL,0) - coords(ALL,in) );
    dist = dist - thickness_/2.0; 
    
    if ( std::abs(dist) < tol )
    {
      // minDist = dist;
      imin = in;
      break;
    }
  }

  if ( nodeCount == 4 )
  {
    // linear line element 

    iperm.resize ( 4 );

    if ( imin == 1 )
    {
      //  2-------1       2-------3
      //  |       |  -->  |       |
      //  3-------0       0-------1

      iperm[0] = 3;
      iperm[1] = 0;
      iperm[2] = 2;
      iperm[3] = 1;

      return 4;
    }
    else if ( imin == 3 )
    {
      //  3-------2       2-------3
      //  |       |  -->  |       |
      //  0-------1       0-------1

      iperm[0] = 0;
      iperm[1] = 1;
      iperm[2] = 3;
      iperm[3] = 2;

      return 4;
    }
  }
  else if ( nodeCount == 8 || nodeCount == 9 )
  {
    // quadratic line element 

    iperm.resize ( 6 );

    if ( imin == 1 || imin == 2 )
    {
      //  4---3---2       3---4---5
      //  5  (8)  1  -->  |       |
      //  6---7---0       0---1---2

      iperm[0] = 6;
      iperm[1] = 7;
      iperm[2] = 0;
      iperm[3] = 4;
      iperm[4] = 3;
      iperm[5] = 2;

      return 6;
    }

    if ( imin == 6 || imin == 7 )
    {
      //  6---5---4       3---4---5
      //  7  (8)  3  -->  |       |
      //  0---1---2       0---1---2

      iperm[0] = 0;
      iperm[1] = 1;
      iperm[2] = 2;
      iperm[3] = 6;
      iperm[4] = 5;
      iperm[5] = 4;

      return 6;
    }
  }
  System::out() << "nodeCount " << nodeCount << " imin " << imin << endl;
  throw Error ( JEM_FUNC, "unexpected element configuration" );
}

//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------

Ref<Module>  Quad4InterfModule::makeNew

  ( const String&           name,
    const Properties&       conf,
    const Properties&       props,
    const Properties&       globdat )

{
  return newInstance<Self> ( name );
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   declareQuad4InterfModule
//-----------------------------------------------------------------------

void declareQuad4InterfModule ()
{
  using jive::app::ModuleFactory;

  ModuleFactory::declare ( Quad4InterfModule::TYPE_NAME,
                         & Quad4InterfModule::makeNew );
}
