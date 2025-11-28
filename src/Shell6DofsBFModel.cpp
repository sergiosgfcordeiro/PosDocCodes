/*
 *
 *  Copyright (C) 2025 TU Delft. All rights reserved.
 *
 *  
 *  This class implements the triangular Kirchhoff-Love shell element. It is an assemble of a plate element developed for orthotropic materials to model
 *  the plies in laminated composites and a membrane element with rotational degrees of freedom:
 *  
 *  Allman 1976 - A simple cubic displacement element for plate bending.
 *  Ai et al 2024 - Structural cohesive element for the modelling of delamination in composite laminates without the cohesive zone limit. 
 *  Bergan & Felippa 1985 - A triangular membrane element with rotational degrees of freedom.      
 *  
 *  Author:  S.G.F. Cordeiro, sferreiracorde@tudelft.nl
 *
 *  27 July 2025:
 *  
 */

#include <jem/base/Array.h>
#include <jem/base/array/operators.h>
#include <jem/base/array/select.h>
#include <jem/base/array/tensor.h>
#include <jem/base/array/utilities.h>
#include <jem/base/Error.h>
#include <jem/base/IllegalInputException.h>
#include <jem/base/System.h>
#include <jem/io/FileWriter.h>
#include <jem/numeric/algebra/MatmulChain.h>
#include <jem/numeric/algebra/Cholesky.h>
#include <jem/numeric/algebra/LUSolver.h>

#include <jive/geom/Geometries.h>
#include <jive/model/ModelFactory.h>

#include "Shell6DofsBFModel.h"
#include "TbFiller.h"

#include <iostream>

using namespace jem;
using jem::io::FileWriter;
using jem::numeric::matmul;
using jem::numeric::MatmulChain;
using jem::TensorIndex;

using jive::geom::Geometries;

typedef MatmulChain<double,3>   MChain3;
typedef MatmulChain<double,2>   MChain2;
typedef MatmulChain<double,1>   MChain1;

//======================================================================
//   definition
//======================================================================

const char* Shell6DofsBFModel::DOF_NAMES[6]     = {"u","v","w","wx","wy","wz"};
const char* Shell6DofsBFModel::SHAPE_PROP       = "shape";
const char* Shell6DofsBFModel::MATERIAL_PROP    = "material";
const char* Shell6DofsBFModel::ALPHA_PROP       = "alpha";
const char* Shell6DofsBFModel::BETA_PROP        = "beta";
const char* Shell6DofsBFModel::THICK_PROP       = "thickness";
const char* Shell6DofsBFModel::ORIENTATION_PROP[3] = {"v1","v2","v3"};    
const char* Shell6DofsBFModel::LARGE_DISP_PROP  = "largeDisp";
      idx_t Shell6DofsBFModel::nodesWritten_    = 0;
Ref<PrintWriter> Shell6DofsBFModel::nodeOut_    = NIL;

//-----------------------------------------------------------------------
//   constructor
//-----------------------------------------------------------------------

Shell6DofsBFModel::Shell6DofsBFModel

   ( const String&       name,
     const Properties&   conf,
     const Properties&   props,
     const Properties&   globdat ) : Super(name)
{
  using jive::util::joinNames;
  using jive::geom::IShapeFactory;

  // create myTag_ (last part of myName_)
  
  StringVector names ( StringUtils::split( myName_, '.' ) );
  myTag_     = names [ names.size() - 1 ];

  Properties  myProps = props.getProps ( myName_ );
  Properties  myConf  = conf.makeProps ( myName_ );

  const String context = getContext();

  egroup_ = ElemGroup::get ( myConf, myProps, globdat, context );

  numElem_   = egroup_.size();
  ielems_    . resize( numElem_ );
  ielems_    = egroup_.getIndices ();
  elems_     = egroup_.getElements ( );
  nodes_     = elems_.getNodes     ( );
  nodeRank_  = nodes_.rank         ( );
  rank_      = nodeRank_ - 1; 
  numNode_   = nodes_.size         ( );
  int numDof_  = sizeof(DOF_NAMES)/sizeof(DOF_NAMES[0]);
  strCount_  = STRAIN_COUNTS[rank_];

  // Make sure that the number of spatial dimensions (the rank of the
  // mesh) is valid.

  if ( nodeRank_ < 2 || nodeRank_ > 3 )
  {
    throw IllegalInputException (
      context,
      String::format (
        "invalid node rank: %d (should be 2 or 3)", nodeRank_
      )
    );
  }

  shape_  = IShapeFactory::newInstance(
    joinNames (myName_, SHAPE_PROP ),
    conf,
    props );

  nodeCount_  = shape_->nodeCount   ();
  ipCount_    = shape_->ipointCount ();
  dofCount_   = numDof_ * nodeCount_;

  // Make sure that the rank of the shape matches the rank of the
  // mesh.

  if ( shape_->globalRank() != rank_ )
  {
    throw IllegalInputException (
      context,
      String::format (
        "shape has invalid rank: %d (should be %d)",
        shape_->globalRank (),
        rank_
      )
    );
  }

  // Make sure that each element has the same number of nodes as the
  // shape object.

  elems_.checkSomeElements (
    context,
    ielems_,
    shape_->nodeCount  ()
  );

  dofs_ = XDofSpace::get ( nodes_.getData(), globdat );
  
  dofTypes_.resize( numDof_ );

  for( idx_t i = 0; i < dofTypes_.size(); i++)
  {
    dofTypes_[i] = dofs_->addType ( DOF_NAMES[i]);
  }

  dofs_->addDofs (
    elems_.getUniqueNodesOf ( ielems_ ),
    dofTypes_
  );

  // Compute the total number of integration points.

  idx_t  ipCount = shape_->ipointCount() * egroup_.size();

  // Create a material model object.

  material_ = newMaterial ( MATERIAL_PROP, myConf, myProps, globdat );

  material_-> allocPoints  ( ipCount );

  softening_ = dynamicCast<Softening> ( material_ );

  getShapeGrads_ = getShapeGradsFunc ( rank_ );

  alpha_ = 1.5;

  myProps.find( alpha_, ALPHA_PROP );
  myConf.set  ( ALPHA_PROP, alpha_ );

  beta_ = 0.;

  myProps.find( beta_, BETA_PROP );
  myConf.set  ( BETA_PROP, beta_ );

  thickness_ = 1.;

  myProps.find( thickness_, THICK_PROP );
  myConf.set  ( THICK_PROP, thickness_ );

  orientVec_[0] = 1.0;
  orientVec_[1] = 0.0;
  orientVec_[2] = 0.0;

  myProps.find( orientVec_[0], ORIENTATION_PROP[0] );
  myConf.set  ( ORIENTATION_PROP[0], orientVec_[0] );

  myProps.find( orientVec_[1], ORIENTATION_PROP[1] );
  myConf.set  ( ORIENTATION_PROP[1], orientVec_[1] );

  myProps.find( orientVec_[2], ORIENTATION_PROP[2] );
  myConf.set  ( ORIENTATION_PROP[2], orientVec_[2] );

  largeDisp_ = false;
  myProps.find(   largeDisp_, LARGE_DISP_PROP );
  myConf.set  (   LARGE_DISP_PROP, largeDisp_ );
  // shape_->setLargeDisp ( largeDisp_ );

  crackBandMethod_ = false;
}

Shell6DofsBFModel::~Shell6DofsBFModel()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void Shell6DofsBFModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  Properties  myProps  = props.findProps ( myName_ );
  Properties  matProps = myProps.findProps ( MATERIAL_PROP );

  material_->configure ( matProps );

  Properties  params;
  Properties  locOutProps;

  crackBandMethod_ = ( softening_ != NIL && !material_->isViscous() );

  if ( crackBandMethod_ )
  {
    initCharLength_();
  }

  // read properties related to writing of local stress/strain history
  doLocalOutput_ = false;

  if ( myProps.find ( locOutProps, "localOutput" ) )
  {
    // get element numbers

    locDoAll_ = false;
    locOutProps.find ( locDoAll_, "doAll" );

    ArrayBuffer<idx_t> iebuf;
    ArrayBuffer<idx_t> ipbuf;

    if ( locDoAll_)
    {
      doLocalOutput_ = true;
      for ( idx_t i = 0; i < numElem_; ++i )
      {
        for ( idx_t j = 0; j < ipCount_; ++j )
        {
          iebuf.pushBack ( i );
          ipbuf.pushBack ( j );
        }
      }
      locIpNumbers_.ref ( ipbuf.toArray() );
    }
    else
    {
      locOutProps.get ( locElems_, "ielems" );

      if ( locElems_.size() > 0 )
      {
        if ( locOutProps.find ( locIpNumbers_, "ipnumbers" ) )
        {
          JEM_PRECHECK ( locIpNumbers_.size() == locElems_.size() );
        }
        else
        {
          // default: use first ip (works for Gauss1 integration)
          locIpNumbers_ . resize ( locElems_.size() );
          locIpNumbers_ = 0;
        }

        // find element indices for lookup from ielems_

        for ( idx_t i = 0; i < locElems_.size(); ++i )
        {
          for ( idx_t j = 0; j < ielems_.size(); ++j )
          {
            if ( ielems_[j] == locElems_[i] )
            {
              doLocalOutput_ = true;
              iebuf.pushBack ( j );
              ipbuf.pushBack ( locIpNumbers_[i] );
            }
          }
        }
      }
    }
    locIes_.ref ( iebuf.toArray() );
    locIps_.ref ( ipbuf.toArray() );
  }
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void Shell6DofsBFModel::getConfig 

  ( const Properties& conf,
    const Properties& globdat ) const

{
  Properties  myConf  = conf.makeProps ( myName_ );
  Properties  matConf = myConf.makeProps ( MATERIAL_PROP );

  material_->getConfig ( matConf );

  Properties  locOutConf = myConf.makeProps ( "localOutput" );

  if ( locDoAll_ )
  {
    locOutConf.set ( "doAll", true );
  }
  else
  {
    locOutConf.set ( "ielems", locElems_ );
    locOutConf.set ( "ipnumbers", locIpNumbers_ );
  }
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool Shell6DofsBFModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;

  if ( action == Actions::GET_MATRIX0 
    || action == Actions::GET_INT_VECTOR )
  {
    Ref<MatrixBuilder>  mbuilder;

    Vector  disp;
    Vector  force;

    // Get the current displacements.

    StateVector::get ( disp, dofs_, globdat );

    // Get the matrix builder and the internal force vector.

    params.find( mbuilder, ActionParams::MATRIX0 );
    params.get ( force,    ActionParams::INT_VECTOR );

    getMatrix_ ( mbuilder, force, disp );

    return true;
  }

  // if ( action == Actions::GET_MATRIX2 )
  // {
  //   Ref<MatrixBuilder> mbuilder;

  //   params.get ( mbuilder, ActionParams::MATRIX2 );
    
  //   getMatrix2_( *mbuilder );

  //   return true;
  // }

  if ( action == "GET_DISSIPATION" )
  {
    getDissipation_ ( params );

    return false;
  }

  if ( action == "DESPAIR" )
  {
    return material_->despair();
  }

  if ( action == "END_DESPAIR" )
  {
    material_->endDespair();

    return true;
  }

  if ( action == Actions::COMMIT )
  {
    material_->commit ();

    if ( doLocalOutput_ )
    {
      Vector  disp;

      StateVector::get ( disp, dofs_, globdat );

      writeLocalOutput_ ( disp, globdat );
    }

    return true;
  }
  
  if ( action == SolverNames::GET_DISS_FORCE )
  {
    Ref<Plasticity> p = dynamicCast<Plasticity> ( material_ );

    if ( p == NIL ) return false;

    Vector disp;
    Vector fDiss;

    StateVector::getOld ( disp, dofs_, globdat );
    globdat.get ( fDiss, SolverNames::DISSIPATION_FORCE );

    getDissForce_ ( p, fDiss, disp );

    return true;
  }

  if ( action == SolverNames::SET_STEP_SIZE )
  {
    double             dt;
    params.get       ( dt, SolverNames::STEP_SIZE );

    material_->setDT ( dt );
    return true;
  }

  if ( action == Actions::GET_TABLE )
  {
    return getTable_ ( params, globdat );
  }

  if ( action == "WRITE_XOUTPUT" )
  {
    if ( dw_.amFirstWriter ( this ) )
    {
      writeDisplacements_ ( params, globdat );
    }
  }

  return false;
}


//-----------------------------------------------------------------------
//   getMatrix_
//-----------------------------------------------------------------------


void Shell6DofsBFModel::getMatrix_

  ( Ref<MatrixBuilder>  mbuilder,
    const Vector&       force,
    const Vector&       disp ) const

{
  Matrix      stiff       ( strCount_, strCount_ );
  Matrix      D_mem       ( strCount_, strCount_ );
  Matrix      D_plate     ( strCount_, strCount_ );
  
  Matrix      coords      ( nodeRank_, nodeCount_ );
  Matrix      transMat    ( nodeRank_, nodeRank_ );
  Matrix      transMatDof ( dofCount_, dofCount_ );
  Matrix      xcoords     ( rank_, nodeCount_ );

  Matrix      Lmat        ( 3*nodeCount_, strCount_ );
  Matrix      Bb_mem      ( strCount_, 2*nodeCount_ );
  Matrix      Kb_mem      ( 3*nodeCount_, 3*nodeCount_ );
  Matrix      Gmat        ( 3*nodeCount_, 3*nodeCount_ );
  Matrix      Ginv        ( 3*nodeCount_, 3*nodeCount_ );
  Matrix      Hh          ( nodeCount_, 3*nodeCount_ );

  Matrix      Bh_mem      ( strCount_, nodeCount_ );
  Matrix      Kqh_mem     ( nodeCount_, nodeCount_ );
  Matrix      Kh_mem      ( 3*nodeCount_, 3*nodeCount_ );

  Matrix      K_mem       ( 3*nodeCount_, 3*nodeCount_ );
  Vector      Disp_mem    ( 2*nodeCount_ );

  Matrix      K_plate     ( 3*nodeCount_, 3*nodeCount_ );
  Matrix      b_plate     ( 7, 4*nodeCount_);
  Matrix      Tmat        ( 4*nodeCount_, 3*nodeCount_ );
  Matrix      bT          ( 7, 3*nodeCount_ );
  Matrix      Hmat        ( 7, 7);
  Matrix      Hinv        ( 7, 7);

  Matrix      K_shell     ( dofCount_, dofCount_  );
  Vector      elemForce   ( dofCount_ );
  Vector      elemDisp    ( dofCount_ );
  Vector      elemxDisp   ( dofCount_ );

  Vector      strain_mem  ( strCount_ );
  Vector      stress_mem  ( strCount_ );

  IdxVector   inodes      ( nodeCount_ );
  IdxVector   idofs       ( dofCount_ );

  Cubix       grads       ( rank_, nodeCount_, ipCount_ );
  Vector      ipWeights   ( ipCount_ );

  double Area = 0.0;
  double det;

  MChain1     mc1;
  MChain3     mc3;

  idx_t       ipoint = 0;

  // Iterate over all elements assigned to this model.

  for ( idx_t ie = 0; ie < numElem_; ie++ )
  {
    // Get the global element index.

    idx_t  ielem = ielems_[ie];

    // Get the element coordinates and DOFs.

    elems_.getElemNodes  ( inodes, ielem    );
    nodes_.getSomeCoords ( coords, inodes );
    dofs_->getDofIndices ( idofs,  inodes, dofTypes_ );

    // Get the orientation vector: v.

    Vector v ( nodeRank_ );
    for ( idx_t i = 0; i < nodeRank_; i++ )
    {
      v[i] = orientVec_[i];
    }

    // Get the transformation matrices (Change of basis) to the shell element local system.

    getTransMatrix_ ( transMatDof, transMat, coords, v );

    // Get the nodal coordinates in the Local coordinate systems.

    get2DLocalcoordinates_(xcoords, transMat, coords);

    // Compute the spatial derivatives of the element shape functions.

    shape_->getShapeGradients ( grads, ipWeights, xcoords );
    getArea_(Area, xcoords);
    ipWeights /= Area;

    // Get the Lumping matrix Lmat.

    getLmat_(Lmat, xcoords, alpha_);

    // Get the G matrix relating generalized Dofs to nodal Dofs: Gmat.

    getGmat_(Gmat, xcoords, Area);

    // Get Ginv and the Hh matrices.

    Ginv = Gmat;
    jem::numeric::LUSolver::invert(Ginv, det);

    Hh = Ginv(slice(2*nodeCount_,3*nodeCount_),slice(0,3*nodeCount_));  

    // Get the displacements at the element nodes in the Global coordinate system.

    elemDisp = disp[idofs];

    // Transform the displacements at the element nodes to the Local coordinate system.

    elemxDisp = mc1.matmul ( transMatDof, elemDisp );
    
    // Get the membrane displacements at the element nodes.

    getMembraneDisp_(Disp_mem, elemxDisp);

    // Assemble the basic membrane stiffness matrix.

    Kb_mem = 0.0;
    Kh_mem = 0.0;
    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    { 
      // Compute the B-matrix for the membrane element at this integration point.

      getShapeGrads_ ( Bb_mem, grads(ALL,ALL,ip) );

      // Compute the membrane strain vector of this integration point.

      matmul ( strain_mem, Bb_mem, Disp_mem );

      // Get the tangent stiffness matrix and the stress vector.

      material_->update ( stress_mem, stiff, strain_mem, ipoint++ );

      // Assemble the membrane constitutive matrix.

      D_mem = stiff*thickness_;

      // Compute the basic membrane stiffness matrix.

      Kb_mem += ipWeights[ip] * (1/Area) * mc3.matmul ( Lmat, D_mem, Lmat.transpose() );

      // Get the generalized stiffness matrix for higher-order modes: Kqh.

      getKqhmat_(Kqh_mem, xcoords, D_mem, Area);

      // Compute the higher-order membrane stiffness matrix.

      Kh_mem += ipWeights[ip] * mc3.matmul ( Hh.transpose(), Kqh_mem, Hh );
    }

    // Compute the membrane stiffness matrix.

    K_mem = Kb_mem + beta_ * Kh_mem;

    // Get the tangent stiffness matrix and the stress vector at the first integration point.

    int ip = 0;
    getShapeGrads_ ( Bb_mem, grads(ALL,ALL,ip) );
    matmul ( strain_mem, Bb_mem, Disp_mem );
    material_->update ( stress_mem, stiff, strain_mem, ipoint++ );
    
    // Assemble the plate stiffness matrix.

    D_plate = stiff*(pow(thickness_, 3))/12.0;

    // Assemble the element matrix: Hmat.

    getHmatAnalytic_(Hmat, D_plate, xcoords);

    // Assemble the B matrix of the plate element.

    getBmat_(b_plate, D_plate, xcoords);

    // Get the element transformation matrix Tmat.

    getTmat_(Tmat, xcoords);

    // Assemble the plate stiffness matrix.

    bT = matmul ( b_plate, Tmat );
    Hinv = Hmat;
    jem::numeric::Cholesky::invert(Hinv);
    K_plate = mc3.matmul ( bT.transpose(), Hinv, bT );

    // // Print the contents of the K_plate matrix for debugging purposes.
    // for ( idx_t i = 0; i < K_plate.size(0); i++ )
    // {
    //   for ( idx_t j = 0; j < K_plate.size(1); j++ )
    //   {
    //     jem::System::out() << K_plate(i,j) << " ";
    //   }
    //   jem::System::out() << "\n";
    // }

    // Assemble the shell stiffness matrix.

    K_shell = 0.0;
    for (idx_t i = 0; i < nodeCount_; ++i) {
      for (idx_t j = 0; j < nodeCount_; ++j) {

        // u,v, wz in-plane (membrane) contribution.
        K_shell(6*i+0, 6*j+0) = K_mem(3*i+0, 3*j+0);     
        K_shell(6*i+0, 6*j+1) = K_mem(3*i+0, 3*j+1);
        K_shell(6*i+0, 6*j+5) = K_mem(3*i+0, 3*j+2);

        K_shell(6*i+1, 6*j+0) = K_mem(3*i+1, 3*j+0);
        K_shell(6*i+1, 6*j+1) = K_mem(3*i+1, 3*j+1);
        K_shell(6*i+1, 6*j+5) = K_mem(3*i+1, 3*j+2);

        K_shell(6*i+5, 6*j+0) = K_mem(3*i+2, 3*j+0);
        K_shell(6*i+5, 6*j+1) = K_mem(3*i+2, 3*j+1);
        K_shell(6*i+5, 6*j+5) = K_mem(3*i+2, 3*j+2);

        // w, wx, wy out-of-plane (plate) contribution.
        K_shell(6*i+2, 6*j+2) = K_plate(3*i+0, 3*j+0);
        K_shell(6*i+2, 6*j+3) = K_plate(3*i+0, 3*j+1);
        K_shell(6*i+2, 6*j+4) = K_plate(3*i+0, 3*j+2);

        K_shell(6*i+3, 6*j+2) = K_plate(3*i+1, 3*j+0);
        K_shell(6*i+3, 6*j+3) = K_plate(3*i+1, 3*j+1);
        K_shell(6*i+3, 6*j+4) = K_plate(3*i+1, 3*j+2);

        K_shell(6*i+4, 6*j+2) = K_plate(3*i+2, 3*j+0);
        K_shell(6*i+4, 6*j+3) = K_plate(3*i+2, 3*j+1);
        K_shell(6*i+4, 6*j+4) = K_plate(3*i+2, 3*j+2);
      }
    }

    // Transforme the shell stiffness matrix back to the global system.

    K_shell = mc3.matmul ( transMatDof.transpose(), K_shell, transMatDof );

    // Add the element matrix to the global stiffness matrix.

    if ( mbuilder != NIL )
    {
      mbuilder->addBlock ( idofs, idofs, K_shell );
    }

    // Compute and add the element force vector to the global force vector.

    elemForce = matmul ( K_shell, elemDisp );
    force[idofs] += elemForce;
  }
}


//-----------------------------------------------------------------------
//   writeLocalOutput_       NOT DONE FOR THE SHELL MODEL YET
//-----------------------------------------------------------------------


void Shell6DofsBFModel::writeLocalOutput_

  ( const Vector&       disp,
    const Properties&   globdat )

{
  idx_t       it;

  Matrix      stiff       ( strCount_, strCount_ );
  Matrix      coords      ( nodeRank_, nodeCount_ );
  Matrix      transMat    ( nodeRank_, nodeRank_ );
  Matrix      transMatDof ( dofCount_, dofCount_ );
  Matrix      xcoords     ( rank_, nodeCount_ );
  Vector      elemDisp    ( dofCount_ );
  Vector      strain      ( strCount_ );
  Vector      stress      ( strCount_ );

  Matrix      b           ( strCount_, dofCount_  );

  IdxVector   inodes      ( nodeCount_ );
  IdxVector   idofs       ( dofCount_  );

  Cubix       grads       ( rank_, nodeCount_, ipCount_ );
  Vector      ipWeights   ( ipCount_   );

  if ( locOut_ == NIL )
  {
    initLocOutWriter_ ();
  }

  globdat.get ( it, Globdat::TIME_STEP );

  *locOut_ << "newXOutput " << it << " (strain, stress, dissipation, history)\n";

  // Iterate over all ips for which output is to be written

  for ( idx_t i = 0; i < locIes_.size(); ++i )
  {
    idx_t ie = locIes_[i];
    idx_t ip = locIps_[i];
    idx_t ielem = ielems_[ie];
    idx_t ipoint = ipCount_ * ie + ip;

    elems_.getElemNodes  ( inodes, ielem    );
    nodes_.getSomeCoords ( coords, inodes );
    dofs_->getDofIndices ( idofs,  inodes, dofTypes_ );

    Vector v ( nodeRank_ );
    for ( idx_t i = 0; i < nodeRank_; i++ )
    {
      v[i] = orientVec_[i];
    }

    getTransMatrix_ ( transMatDof, transMat, coords, v );

    get2DLocalcoordinates_(xcoords, transMat, coords);

    shape_->getShapeGradients ( grads, ipWeights, xcoords );

    elemDisp = select ( disp, idofs );

    getShapeGrads_ ( b, grads(ALL,ALL,ip) );
    matmul ( strain, b, elemDisp );

    if ( crackBandMethod_ )
    {
      double le = charLength_[ielem];
      softening_->update ( stress, stiff, strain, ipoint, le );
    }
    else
    {
      material_->update ( stress, stiff, strain, ipoint );
    }
    for ( idx_t j = 0; j < strCount_; ++j ) 
    { 
      *locOut_ << strain[j] << " ";
    }
    for ( idx_t j = 0; j < strCount_; ++j ) 
    { 
      *locOut_ << stress[j] << " ";
    }
    *locOut_ << material_->giveDissipation ( ipoint );
    *locOut_ << endl;
  }
  locOut_->flush();
}


//-----------------------------------------------------------------------
//   getXOutTable_       NOT DONE FOR THE SHELL MODEL YET
//-----------------------------------------------------------------------


void Shell6DofsBFModel::getXOutTable_

  ( Ref<XTable>        table,
    const Vector&      weights,
    const String&      contents,
    const Vector&      disp )

{
  Vector       ndWeights   ( nodeCount_ );
  StringVector hisNames    = material_->getHistoryNames ();

  Cubix        grads       ( rank_, nodeCount_, ipCount_ );
  Matrix       coords      ( nodeRank_,     nodeCount_ );
  Matrix       transMat    ( nodeRank_, nodeRank_ );
  Matrix       transMatDof ( dofCount_, dofCount_ );
  Matrix       xcoords     ( rank_,     nodeCount_ );
  Matrix       b           ( strCount_, dofCount_  );
  Matrix       stiff       ( strCount_, strCount_  );

  Vector       elemDisp    ( dofCount_ );

  IdxVector    inodes      ( nodeCount_ );
  IdxVector    idofs       ( dofCount_  );

  const bool   tri6 = ( shape_->getGeometry() == Geometries::TRIANGLE
                      && nodeCount_ == 6 );

  // tell TbFiller which types are available to write

  TbFiller   tbFiller   ( rank_ );

  Slice      iistrain   = tbFiller.announce ( "strain.tensor" );
  Slice      iistress   = tbFiller.announce ( "stress.tensor" );
  Slice      iihistory  = tbFiller.announce ( hisNames );

  Vector     ipValues   ( tbFiller.typeCount() );

  Vector     strain     ( ipValues[iistrain]   );
  Vector     stress     ( ipValues[iistress]   );
  Vector     history    ( ipValues[iihistory]  );

  // Let TbFiller find out which columns of ndValues to write to 
  // which columns of the table (based on filter in input file)

  IdxVector  i2table;
  IdxVector  jcols;

  tbFiller . setFilter   ( contents );
  tbFiller . prepareTable( i2table, jcols, table );

  Matrix     ndValuesOut ( nodeCount_, i2table.size() );
  Vector     ipValuesOut ( i2table.size() );

  // fill table in loop over elements

  idx_t      ipCount;
  idx_t      ipoint = 0;

  Vector     ipWeights;

  // Add the columns for the stress components to the table.

  const idx_t nel = ielems_.size();

  Matrix      sfuncs     = shape_->getShapeFunctions ();

  for ( idx_t ie = 0; ie < nel; ++ie )
  {
    idx_t ielem = ielems_[ie];

    ndValuesOut = 0.;
    ndWeights   = 0.;

    elems_.getElemNodes  ( inodes, ielem );
    dofs_->getDofIndices ( idofs,  inodes,  dofTypes_ );

    ipCount  = shape_->ipointCount ();

    ipWeights.resize ( ipCount   );

    nodes_.getSomeCoords ( coords, inodes );

    Vector v ( nodeRank_ );
    for ( idx_t i = 0; i < nodeRank_; i++ )
    {
      v[i] = orientVec_[i];
    }

    getTransMatrix_ ( transMatDof, transMat, coords, v );

    get2DLocalcoordinates_(xcoords, transMat, coords);

    shape_->getShapeGradients ( grads, ipWeights, xcoords );

    elemDisp = select ( disp, idofs );

    // Iterate over the integration points.
    // Gather all data, no matter which is asked, to keep code neat
    // The option to specify output is primarily for disk size, not CPU time

    for ( idx_t ip = 0; ip < ipCount; ip++, ++ipoint )
    {
      getShapeGrads_ ( b, grads(ALL,ALL,ip) );
      matmul ( strain, b, elemDisp );

      if ( crackBandMethod_ )
      {
        double le = charLength_[ielem];
        softening_->update ( stress, stiff, strain, ipoint, le );
      }
      else
      {
        material_->update ( stress, stiff, strain, ipoint );
      }
      material_-> getHistory ( history, ipoint );

      ipValuesOut  = ipValues[i2table];
      ndValuesOut += matmul ( sfuncs(ALL,ip), ipValuesOut ); 
      ndWeights   += sfuncs(ALL,ip);
    }

    if ( tri6 ) TbFiller::permTri6 ( ndWeights, ndValuesOut );

    select ( weights, inodes ) += ndWeights;

    // Add the stresses to the table.

    table->addBlock ( inodes, jcols, ndValuesOut );
  }
}


//-----------------------------------------------------------------------
//   getTransMatrix_
//-----------------------------------------------------------------------


void Shell6DofsBFModel::getTransMatrix_
   
  ( Matrix&             transMatDof,
    Matrix&             transMat,
    const Matrix&       coords,
    const Vector&       v   ) const

{
  Vector      e1            ( nodeRank_ );
  Vector      e2            ( nodeRank_ );
  Vector      e3            ( nodeRank_ );
  Vector      e1_l          ( nodeRank_ );
  Vector      e2_l          ( nodeRank_ );
  Vector      e3_l          ( nodeRank_ );
  Matrix      transMat_WiRi ( nodeRank_, nodeRank_ );
  Matrix      transMatRot   ( nodeRank_, nodeRank_ );
  Matrix      transMatBlock ( 2*nodeRank_, 2*nodeRank_ );

  e1 = {1.0, 0.0, 0.0};
  e2 = {0.0, 1.0, 0.0};
  e3 = {0.0, 0.0, 1.0};

  // Get the local system and the transformation matrix.

  e3_l[0] = (coords(1, 1) - coords(1, 0))*
              (coords(2, 2) - coords(2, 0)) -
              (coords(1, 2) - coords(1, 0))*
              (coords(2, 1) - coords(2, 0));
  e3_l[1] = (coords(2, 1) - coords(2, 0))*
              (coords(0, 2) - coords(0, 0)) -
              (coords(2, 2) - coords(2, 0))*
              (coords(0, 1) - coords(0, 0));
  e3_l[2] = (coords(0, 1) - coords(0, 0))*
              (coords(1, 2) - coords(1, 0)) -
              (coords(0, 2) - coords(0, 0))*
              (coords(1, 1) - coords(1, 0));
  double norm = sqrt ( e3_l[0]*e3_l[0] + 
                       e3_l[1]*e3_l[1] + 
                       e3_l[2]*e3_l[2] );
  e3_l  /= norm;
    
  e1_l = v - (v[0]*e3_l[0] + v[1]*e3_l[1] + v[2]*e3_l[2])*e3_l;

  norm = sqrt ( e1_l[0]*e1_l[0] +
                e1_l[1]*e1_l[1] +
                e1_l[2]*e1_l[2] );
  e1_l /= norm;

  e2_l[0] = e3_l[1]*e1_l[2] - e3_l[2]*e1_l[1];
  e2_l[1] = e3_l[2]*e1_l[0] - e3_l[0]*e1_l[2];
  e2_l[2] = e3_l[0]*e1_l[1] - e3_l[1]*e1_l[0];

  // Get the transformation matrix for translations.

  transMat(0,0) = e1_l[0]*e1[0] + e1_l[1]*e1[1] + e1_l[2]*e1[2];
  transMat(0,1) = e1_l[0]*e2[0] + e1_l[1]*e2[1] + e1_l[2]*e2[2];
  transMat(0,2) = e1_l[0]*e3[0] + e1_l[1]*e3[1] + e1_l[2]*e3[2];
  transMat(1,0) = e2_l[0]*e1[0] + e2_l[1]*e1[1] + e2_l[2]*e1[2];
  transMat(1,1) = e2_l[0]*e2[0] + e2_l[1]*e2[1] + e2_l[2]*e2[2];
  transMat(1,2) = e2_l[0]*e3[0] + e2_l[1]*e3[1] + e2_l[2]*e3[2];
  transMat(2,0) = e3_l[0]*e1[0] + e3_l[1]*e1[1] + e3_l[2]*e1[2];
  transMat(2,1) = e3_l[0]*e2[0] + e3_l[1]*e2[1] + e3_l[2]*e2[2];
  transMat(2,2) = e3_l[0]*e3[0] + e3_l[1]*e3[1] + e3_l[2]*e3[2];

  // Get the transformation matrix for rotations.

  transMat_WiRi = 0.0;
  transMat_WiRi(0,0) = 0.0;
  transMat_WiRi(0,1) = -1.0;
  transMat_WiRi(1,0) = 1.0;
  transMat_WiRi(1,1) = 0.0;
  transMat_WiRi(2,2) = 1.0;

  transMatRot = matmul ( transMat_WiRi, transMat );

  // Get the transformation matrix of the Dofs of a single node.

  transMatBlock = 0.0;
  transMatBlock(slice(0,nodeRank_),slice(0,nodeRank_)) = transMat;
  transMatBlock(slice(nodeRank_,2*nodeRank_),slice(nodeRank_,2*nodeRank_)) = transMatRot;

  // Get the transformation matrix of the nodal element Dofs.

  transMatDof = 0.0;
  transMatDof(slice(0,6),slice(0,6)) = transMatBlock;
  transMatDof(slice(6,12),slice(6,12)) = transMatBlock;
  transMatDof(slice(12,18),slice(12,18)) = transMatBlock;
}


//-----------------------------------------------------------------------
//   get2DLocalcoordinates_
//-----------------------------------------------------------------------


void Shell6DofsBFModel::get2DLocalcoordinates_

  ( Matrix&  xcoords,
    const Matrix&   transMat,
    const Matrix&   coords )   const
    
{
  using jem::ALL;
  using jem::numeric::matmul;
  Matrix  xcoords3D ( coords.size(0), coords.size(1) );

  // Get the cg coordinates of the element.

  double Xcg = 0.0;
  double Ycg = 0.0;
  double Zcg = 0.0;
  for ( idx_t i = 0; i < coords.size(1); i++ )
  {
    Xcg += coords(0, i);
    Ycg += coords(1, i);
    Zcg += coords(2, i);
  }
  Xcg /= coords.size(1);
  Ycg /= coords.size(1);
  Zcg /= coords.size(1);

  // Translation to the baricenter of the element.

  xcoords3D(0, ALL) = coords(0, ALL) - Xcg;
  xcoords3D(1, ALL) = coords(1, ALL) - Ycg;
  xcoords3D(2, ALL) = coords(2, ALL) - Zcg;

  // Apply the transformation to the 3Dxcoords.

  for ( idx_t i = 0; i < xcoords.size(1); i++ )
  {
    xcoords3D(ALL, i) = matmul ( transMat, xcoords3D(ALL, i) );
  }

  // Get the plane xcoords from the xcoords3D.

  xcoords(0, ALL) = xcoords3D(0, ALL);
  xcoords(1, ALL) = xcoords3D(1, ALL);
}


//-----------------------------------------------------------------------
//   getArea_
//-----------------------------------------------------------------------


void Shell6DofsBFModel::getArea_

  ( double&  Area,
    const Matrix&   xcoords) const
  
{
  double x1, x2, x3, y1, y2, y3;

  // Get the Area of the triangle.

  x1 = xcoords(0, 0);
  x2 = xcoords(0, 1);
  x3 = xcoords(0, 2);
  y1 = xcoords(1, 0);
  y2 = xcoords(1, 1);
  y3 = xcoords(1, 2);

  Area = (1.0/2.0) * (x2*y3 - x3*y2 + x3*y1 - x1*y3 + x1*y2 - x2*y1);
}


//-----------------------------------------------------------------------
//   getLmat_
//-----------------------------------------------------------------------


void Shell6DofsBFModel::getLmat_

  ( Matrix&  Lmat,
    const Matrix&   xcoords,
    const double&   alpha_ ) const
  
{
  using jem::ALL;

  Matrix      L1 ( nodeRank_, strCount_ );
  Matrix      L2 ( nodeRank_, strCount_ );
  Matrix      L3 ( nodeRank_, strCount_ );
  double y_ki, y_ji, y_kj, x_ij, x_jk, x_ik;

  // Get the L1 part of the Lumping matrix

  y_ki = (xcoords(1, 1) - xcoords(1, 2));
  y_ji = (xcoords(1, 0) - xcoords(1, 2));
  y_kj = (xcoords(1, 1) - xcoords(1, 0));

  x_ij = (xcoords(0, 2) - xcoords(0, 0));
  x_jk = (xcoords(0, 0) - xcoords(0, 1));
  x_ik = (xcoords(0, 2) - xcoords(0, 1));

  L1 = 0.0;
  L1(0, 0) = y_ki;
  L1(0, 2) = x_ik;
  L1(1,1) = x_ik;
  L1(1, 2) = y_ki;
  L1(2,0) = (1.0/6.0)*alpha_*(y_ji*y_ji - y_kj*y_kj);
  L1(2,1) = (1.0/6.0)*alpha_*(x_ij*x_ij - x_jk*x_jk);
  L1(2,2) = (1.0/3.0)*alpha_*(x_ij*y_ji - x_jk*y_kj);
  L1 = (1.0/2.0) * L1;

  // Get the L2 part of the Lumping matrix

  y_ki = (xcoords(1, 2) - xcoords(1, 0));
  y_ji = (xcoords(1, 1) - xcoords(1, 0));
  y_kj = (xcoords(1, 2) - xcoords(1, 1));

  x_ij = (xcoords(0, 0) - xcoords(0, 1));
  x_jk = (xcoords(0, 1) - xcoords(0, 2));
  x_ik = (xcoords(0, 0) - xcoords(0, 2));

  L2 = 0.0;
  L2(0, 0) = y_ki;
  L2(0, 2) = x_ik;
  L2(1,1) = x_ik;
  L2(1, 2) = y_ki;
  L2(2,0) = (1.0/6.0)*alpha_*(y_ji*y_ji - y_kj*y_kj);
  L2(2,1) = (1.0/6.0)*alpha_*(x_ij*x_ij - x_jk*x_jk);
  L2(2,2) = (1.0/3.0)*alpha_*(x_ij*y_ji - x_jk*y_kj);
  L2 = (1.0/2.0) * L2;

  // Get the L3 part of the Lumping matrix

  y_ki = (xcoords(1, 0) - xcoords(1, 1));
  y_ji = (xcoords(1, 2) - xcoords(1, 1));
  y_kj = (xcoords(1, 0) - xcoords(1, 2));

  x_ij = (xcoords(0, 1) - xcoords(0, 2));
  x_jk = (xcoords(0, 2) - xcoords(0, 0));
  x_ik = (xcoords(0, 1) - xcoords(0, 0));

  L3 = 0.0;
  L3(0, 0) = y_ki;
  L3(0, 2) = x_ik;
  L3(1,1) = x_ik;
  L3(1, 2) = y_ki;
  L3(2,0) = (1.0/6.0)*alpha_*(y_ji*y_ji - y_kj*y_kj);
  L3(2,1) = (1.0/6.0)*alpha_*(x_ij*x_ij - x_jk*x_jk);
  L3(2,2) = (1.0/3.0)*alpha_*(x_ij*y_ji - x_jk*y_kj);
  L3 = (1.0/2.0) * L3;

  // Get the Lumping matrix

  Lmat = 0.0;
  Lmat(slice(0,nodeRank_),ALL) = L1;
  Lmat(slice(nodeRank_,2*nodeRank_),ALL) = L2;
  Lmat(slice(2*nodeRank_,3*nodeRank_),ALL) = L3;
}


//-----------------------------------------------------------------------
//   getGmat_
//-----------------------------------------------------------------------


void Shell6DofsBFModel::getGmat_

  ( Matrix&  Gmat,
    const Matrix&   xcoords,
    const double&   Area ) const
  
{
  double lambda;
  double x_c, y_c;

  Vector Xi  ( nodeRank_ );
  Vector Eta ( nodeRank_ );
  Vector xm  ( nodeCount_ );
  Vector ym  ( nodeCount_ );
  Vector lm  ( nodeCount_ );
  Vector S   ( nodeCount_ );
  Vector C   ( nodeCount_ );
  Matrix a   ( nodeRank_, nodeCount_ );
  Matrix b   ( nodeRank_, nodeCount_ );

  Matrix Grc   ( 3*nodeCount_, 2*nodeCount_ );
  Matrix Gh    ( 3*nodeCount_, nodeCount_ );
  Matrix Gh_11  ( nodeCount_, 1 );
  Matrix Gh_12  ( nodeCount_, 1 );
  Matrix Gh_13  ( nodeCount_, 1 );
  Matrix Gh_21  ( nodeCount_, 1 );
  Matrix Gh_22  ( nodeCount_, 1 );
  Matrix Gh_23  ( nodeCount_, 1 );
  Matrix Gh_31  ( nodeCount_, 1 );
  Matrix Gh_32  ( nodeCount_, 1 );
  Matrix Gh_33  ( nodeCount_, 1 );

  // Compute geometric parameters for Grc matrix.
  lambda = 1.0 / sqrt(Area);
  
  x_c = (xcoords(0, 0) + xcoords(0, 1) + xcoords(0, 2)) / 3.0;
  y_c = (xcoords(1, 0) + xcoords(1, 1) + xcoords(1, 2)) / 3.0;
 
  Xi[0] = lambda*(xcoords(0, 0) - x_c);
  Xi[1] = lambda*(xcoords(0, 1) - x_c);
  Xi[2] = lambda*(xcoords(0, 2) - x_c);
  
  Eta[0] = lambda*(xcoords(1, 0) - y_c);
  Eta[1] = lambda*(xcoords(1, 1) - y_c);
  Eta[2] = lambda*(xcoords(1, 2) - y_c);

  // Compute the Grc matrix
  Grc = 0.0;
  Grc(0, 0) = 1.0;
  Grc(3, 0) = 1.0;
  Grc(6, 0) = 1.0; 
  Grc(1, 1) = 1.0;
  Grc(4, 1) = 1.0;
  Grc(7, 1) = 1.0; 
  Grc(0,2) = -Eta[0];
  Grc(1,2) = Xi[0];
  Grc(2,2) = lambda;
  Grc(3,2) = -Eta[1];
  Grc(4,2) = Xi[1];
  Grc(5,2) = lambda;
  Grc(6,2) = -Eta[2];
  Grc(7,2) = Xi[2];
  Grc(8,2) = lambda;
  Grc(0,3) = Xi[0];
  Grc(3,3) = Xi[1];
  Grc(6,3) = Xi[2];
  Grc(1,4) = Eta[0];
  Grc(4,4) = Eta[1];
  Grc(7,4) = Eta[2];
  Grc(0,5) = Eta[0];
  Grc(1,5) = Xi[0];
  Grc(3,5) = Eta[1];
  Grc(4,5) = Xi[1];
  Grc(6,5) = Eta[2];
  Grc(7,5) = Xi[2];
  
  // Compute more geometric parameters for Gh matrix.
  xm[0] = (xcoords(0, 1) + xcoords(0, 2)) / 2.0;
  ym[0] = (xcoords(1, 1) + xcoords(1, 2)) / 2.0;
  lm[0] = sqrt(pow((xm[0] - xcoords(0, 0)),2) + pow((ym[0] - xcoords(1, 0)),2));
  S[0] = (ym[0]-xcoords(1, 0)) / lm[0];
  C[0] = (xm[0]-xcoords(0, 0)) / lm[0];
  a(0,0) = -(1.0/2.0)*S[0]*pow(C[0],2);
  a(1,0) = pow(C[0],3);
  a(2,0) = (1.0/2.0)*pow(S[0],3) + S[0]*pow(C[0],2);
  b(0,0) = -pow(S[0],2)*C[0] - (1.0/2.0)*pow(C[0],3);
  b(1,0) = -pow(S[0],3);
  b(2,0) = (1.0/2.0)*pow(S[0],2)*C[0];

  xm[1] = (xcoords(0, 2) + xcoords(0, 0)) / 2.0;
  ym[1] = (xcoords(1, 2) + xcoords(1, 0)) / 2.0;
  lm[1] = sqrt(pow((xm[1] - xcoords(0, 1)),2) + pow((ym[1] - xcoords(1, 1)),2));
  S[1] = (ym[1] - xcoords(1, 1)) / lm[1];
  C[1] = (xm[1] - xcoords(0, 1)) / lm[1];
  a(0,1) = -(1.0/2.0)*S[1]*pow(C[1],2);
  a(1,1) = pow(C[1],3);
  a(2,1) = (1.0/2.0)*pow(S[1],3) + S[1]*pow(C[1],2);
  b(0,1) = -pow(S[1],2)*C[1] - (1.0/2.0)*pow(C[1],3);
  b(1,1) = -pow(S[1],3);
  b(2,1) = (1.0/2.0)*pow(S[1],2)*C[1];

  xm[2] = (xcoords(0, 0) + xcoords(0, 1)) / 2.0;
  ym[2] = (xcoords(1, 0) + xcoords(1, 1)) / 2.0;
  lm[2] = sqrt(pow((xm[2] - xcoords(0, 2)),2) + pow((ym[2] - xcoords(1, 2)),2));
  S[2] = (ym[2] - xcoords(1, 2)) / lm[2];
  C[2] = (xm[2] - xcoords(0, 2)) / lm[2];
  a(0,2) = -(1.0/2.0)*S[2]*pow(C[2],2);
  a(1,2) = pow(C[2],3);
  a(2,2) = (1.0/2.0)*pow(S[2],3) + S[2]*pow(C[2],2);
  b(0,2) = -pow(S[2],2)*C[2] - (1.0/2.0)*pow(C[2],3);
  b(1,2) = -pow(S[2],3);
  b(2,2) = (1.0/2.0)*pow(S[2],2)*C[2];

  // Compute the Gh_ij Column  matrix
  Gh_11(0,0) = a(0,0)*Xi[0]*Xi[0] + a(1,0)*Xi[0]*Eta[0] + a(2,0)*Eta[0]*Eta[0];
  Gh_11(1,0) = b(0,0)*Xi[0]*Xi[0] + b(1,0)*Xi[0]*Eta[0] + b(2,0)*Eta[0]*Eta[0];
  Gh_11(2,0) = -lambda*(C[0]*Xi[0] + S[0]*Eta[0]);

  Gh_12(0,0) = a(0,1)*Xi[0]*Xi[0] + a(1,1)*Xi[0]*Eta[0] + a(2,1)*Eta[0]*Eta[0];
  Gh_12(1,0) = b(0,1)*Xi[0]*Xi[0] + b(1,1)*Xi[0]*Eta[0] + b(2,1)*Eta[0]*Eta[0];
  Gh_12(2,0) = -lambda*(C[1]*Xi[0] + S[1]*Eta[0]);

  Gh_13(0,0) = a(0,2)*Xi[0]*Xi[0] + a(1,2)*Xi[0]*Eta[0] + a(2,2)*Eta[0]*Eta[0];
  Gh_13(1,0) = b(0,2)*Xi[0]*Xi[0] + b(1,2)*Xi[0]*Eta[0] + b(2,2)*Eta[0]*Eta[0];
  Gh_13(2,0) = -lambda*(C[2]*Xi[0] + S[2]*Eta[0]);

  Gh_21(0,0) = a(0,0)*Xi[1]*Xi[1] + a(1,0)*Xi[1]*Eta[1] + a(2,0)*Eta[1]*Eta[1];
  Gh_21(1,0) = b(0,0)*Xi[1]*Xi[1] + b(1,0)*Xi[1]*Eta[1] + b(2,0)*Eta[1]*Eta[1];
  Gh_21(2,0) = -lambda*(C[0]*Xi[1] + S[0]*Eta[1]);

  Gh_22(0,0) = a(0,1)*Xi[1]*Xi[1] + a(1,1)*Xi[1]*Eta[1] + a(2,1)*Eta[1]*Eta[1];
  Gh_22(1,0) = b(0,1)*Xi[1]*Xi[1] + b(1,1)*Xi[1]*Eta[1] + b(2,1)*Eta[1]*Eta[1];
  Gh_22(2,0) = -lambda*(C[1]*Xi[1] + S[1]*Eta[1]);

  Gh_23(0,0) = a(0,2)*Xi[1]*Xi[1] + a(1,2)*Xi[1]*Eta[1] + a(2,2)*Eta[1]*Eta[1];
  Gh_23(1,0) = b(0,2)*Xi[1]*Xi[1] + b(1,2)*Xi[1]*Eta[1] + b(2,2)*Eta[1]*Eta[1];
  Gh_23(2,0) = -lambda*(C[2]*Xi[1] + S[2]*Eta[1]);

  Gh_31(0,0) = a(0,0)*Xi[2]*Xi[2] + a(1,0)*Xi[2]*Eta[2] + a(2,0)*Eta[2]*Eta[2];
  Gh_31(1,0) = b(0,0)*Xi[2]*Xi[2] + b(1,0)*Xi[2]*Eta[2] + b(2,0)*Eta[2]*Eta[2];
  Gh_31(2,0) = -lambda*(C[0]*Xi[2] + S[0]*Eta[2]);

  Gh_32(0,0) = a(0,1)*Xi[2]*Xi[2] + a(1,1)*Xi[2]*Eta[2] + a(2,1)*Eta[2]*Eta[2];
  Gh_32(1,0) = b(0,1)*Xi[2]*Xi[2] + b(1,1)*Xi[2]*Eta[2] + b(2,1)*Eta[2]*Eta[2];
  Gh_32(2,0) = -lambda*(C[1]*Xi[2] + S[1]*Eta[2]);

  Gh_33(0,0) = a(0,2)*Xi[2]*Xi[2] + a(1,2)*Xi[2]*Eta[2] + a(2,2)*Eta[2]*Eta[2];
  Gh_33(1,0) = b(0,2)*Xi[2]*Xi[2] + b(1,2)*Xi[2]*Eta[2] + b(2,2)*Eta[2]*Eta[2];
  Gh_33(2,0) = -lambda*(C[2]*Xi[2] + S[2]*Eta[2]);

  // Compute the Gh matrix
  Gh = 0.0;
  Gh(slice(0,nodeCount_),0) = Gh_11(slice(0,nodeCount_),0);  
  Gh(slice(0,nodeCount_),1) = Gh_12(slice(0,nodeCount_),0);  
  Gh(slice(0,nodeCount_),2) = Gh_13(slice(0,nodeCount_),0);  
    
  Gh(slice(nodeCount_,2*nodeCount_),0) = Gh_21(slice(0,nodeCount_),0);  
  Gh(slice(nodeCount_,2*nodeCount_),1) = Gh_22(slice(0,nodeCount_),0);  
  Gh(slice(nodeCount_,2*nodeCount_),2) = Gh_23(slice(0,nodeCount_),0);  

  Gh(slice(2*nodeCount_,3*nodeCount_),0) = Gh_31(slice(0,nodeCount_),0);  
  Gh(slice(2*nodeCount_,3*nodeCount_),1) = Gh_32(slice(0,nodeCount_),0);  
  Gh(slice(2*nodeCount_,3*nodeCount_),2) = Gh_33(slice(0,nodeCount_),0);  

  // Compute the Gmat matrix
  Gmat = 0.0;
  Gmat(slice(0,3*nodeCount_),slice(0,2*nodeCount_)) = Grc;
  Gmat(slice(0,3*nodeCount_),slice(2*nodeCount_,3*nodeCount_)) = Gh;
}


//-----------------------------------------------------------------------
//   getKqhmat_
//-----------------------------------------------------------------------


void Shell6DofsBFModel::getKqhmat_

  ( Matrix&  Kqh_mem,
    const Matrix&   xcoords,
    const Matrix&   D_mem,
    const double&   Area ) const
  
{
  using jem::ALL;

  double lambda;
  double x_c, y_c;
  double J_xixi, J_xieta, J_etaeta;

  MChain3     mc3;

  Vector Xi  ( nodeRank_ );
  Vector Eta ( nodeRank_ );
  Vector xm  ( nodeCount_ );
  Vector ym  ( nodeCount_ );
  Vector lm  ( nodeCount_ );
  Vector S   ( nodeCount_ );
  Vector C   ( nodeCount_ );
  Matrix a   ( nodeRank_, nodeCount_ );
  Matrix b   ( nodeRank_, nodeCount_ );
  Matrix BXi_0 ( strCount_, 1 );
  Matrix BXi_1 ( strCount_, 1 );
  Matrix BXi_2 ( strCount_, 1 );
  Matrix BEta_0 ( strCount_, 1 );
  Matrix BEta_1 ( strCount_, 1 );
  Matrix BEta_2 ( strCount_, 1 );

  Matrix Scalar ( 1, 1 );

  // Compute geometric parameters for Grc matrix.
  lambda = 1.0 / sqrt(Area);
  
  x_c = (xcoords(0, 0) + xcoords(0, 1) + xcoords(0, 2)) / 3.0;
  y_c = (xcoords(1, 0) + xcoords(1, 1) + xcoords(1, 2)) / 3.0;
 
  Xi[0] = lambda*(xcoords(0, 0) - x_c);
  Xi[1] = lambda*(xcoords(0, 1) - x_c);
  Xi[2] = lambda*(xcoords(0, 2) - x_c);
  
  Eta[0] = lambda*(xcoords(1, 0) - y_c);
  Eta[1] = lambda*(xcoords(1, 1) - y_c);
  Eta[2] = lambda*(xcoords(1, 2) - y_c);

  J_xixi = - (1.0/6.0)*Area*(Xi[0]*Xi[1] + Xi[1]*Xi[2] + Xi[2]*Xi[0]);
  J_xieta =  (1.0/12.0)*Area*(Xi[0]*Eta[0] + Xi[1]*Eta[1] + Xi[2]*Eta[2]);
  J_etaeta = - (1.0/6.0)*Area*(Eta[0]*Eta[1] + Eta[1]*Eta[2] + Eta[2]*Eta[0]);
  
  // Compute more geometric parameters for Gh matrix.
  xm[0] = (xcoords(0, 1) + xcoords(0, 2)) / 2.0;
  ym[0] = (xcoords(1, 1) + xcoords(1, 2)) / 2.0;
  lm[0] = sqrt(pow((xm[0] - xcoords(0, 0)),2) + pow((ym[0] - xcoords(1, 0)),2));
  S[0] = (ym[0]-xcoords(1, 0)) / lm[0];
  C[0] = (xm[0]-xcoords(0, 0)) / lm[0];
  a(0,0) = -(1.0/2.0)*S[0]*pow(C[0],2);
  a(1,0) = pow(C[0],3);
  a(2,0) = (1.0/2.0)*pow(S[0],3) + S[0]*pow(C[0],2);
  b(0,0) = -pow(S[0],2)*C[0] - (1.0/2.0)*pow(C[0],3);
  b(1,0) = -pow(S[0],3);
  b(2,0) = (1.0/2.0)*pow(S[0],2)*C[0];

  xm[1] = (xcoords(0, 2) + xcoords(0, 0)) / 2.0;
  ym[1] = (xcoords(1, 2) + xcoords(1, 0)) / 2.0;
  lm[1] = sqrt(pow((xm[1] - xcoords(0, 1)),2) + pow((ym[1] - xcoords(1, 1)),2));
  S[1] = (ym[1] - xcoords(1, 1)) / lm[1];
  C[1] = (xm[1] - xcoords(0, 1)) / lm[1];
  a(0,1) = -(1.0/2.0)*S[1]*pow(C[1],2);
  a(1,1) = pow(C[1],3);
  a(2,1) = (1.0/2.0)*pow(S[1],3) + S[1]*pow(C[1],2);
  b(0,1) = -pow(S[1],2)*C[1] - (1.0/2.0)*pow(C[1],3);
  b(1,1) = -pow(S[1],3);
  b(2,1) = (1.0/2.0)*pow(S[1],2)*C[1];

  xm[2] = (xcoords(0, 0) + xcoords(0, 1)) / 2.0;
  ym[2] = (xcoords(1, 0) + xcoords(1, 1)) / 2.0;
  lm[2] = sqrt(pow((xm[2] - xcoords(0, 2)),2) + pow((ym[2] - xcoords(1, 2)),2));
  S[2] = (ym[2] - xcoords(1, 2)) / lm[2];
  C[2] = (xm[2] - xcoords(0, 2)) / lm[2];
  a(0,2) = -(1.0/2.0)*S[2]*pow(C[2],2);
  a(1,2) = pow(C[2],3);
  a(2,2) = (1.0/2.0)*pow(S[2],3) + S[2]*pow(C[2],2);
  b(0,2) = -pow(S[2],2)*C[2] - (1.0/2.0)*pow(C[2],3);
  b(1,2) = -pow(S[2],3);
  b(2,2) = (1.0/2.0)*pow(S[2],2)*C[2];

  // Compute the B_Xi matrix.

  BXi_0(0,0) = 2*a(0,0)*lambda;
  BXi_0(1,0) = b(1,0)*lambda;
  BXi_0(2,0) = -4*b(2,0)*lambda;

  BXi_1(0,0) = 2*a(0,1)*lambda;
  BXi_1(1,0) = b(1,1)*lambda;
  BXi_1(2,0) = -4*b(2,1)*lambda;

  BXi_2(0,0) = 2*a(0,2)*lambda;
  BXi_2(1,0) = b(1,2)*lambda;
  BXi_2(2,0) = -4*b(2,2)*lambda;

  // Compute the B_Eta matrix.

  BEta_0(0,0) = a(1,0)*lambda;
  BEta_0(1,0) = 2*b(2,0)*lambda;
  BEta_0(2,0) = -4*a(0,0)*lambda;

  BEta_1(0,0) = a(1,1)*lambda;
  BEta_1(1,0) = 2*b(2,1)*lambda;
  BEta_1(2,0) = -4*a(0,1)*lambda;

  BEta_2(0,0) = a(1,2)*lambda;
  BEta_2(1,0) = 2*b(2,2)*lambda;
  BEta_2(2,0) = -4*a(0,2)*lambda;

  // Compute the Kqh_mem matrix.
  Scalar = 0.0;
  Scalar += J_xixi * mc3.matmul ( BXi_0.transpose(), D_mem, BXi_0 );
  Scalar += J_xieta * mc3.matmul ( BXi_0.transpose(), D_mem, BEta_0 );
  Scalar += J_xieta * mc3.matmul ( BEta_0.transpose(), D_mem, BXi_0 );
  Scalar += J_etaeta * mc3.matmul ( BEta_0.transpose(), D_mem, BEta_0 );
  Kqh_mem(0,0) = Scalar(0,0);

  Scalar = 0.0;
  Scalar += J_xixi * mc3.matmul ( BXi_0.transpose(), D_mem, BXi_1 );
  Scalar += J_xieta * mc3.matmul ( BXi_0.transpose(), D_mem, BEta_1 );
  Scalar += J_xieta * mc3.matmul ( BEta_0.transpose(), D_mem, BXi_1 );
  Scalar += J_etaeta * mc3.matmul ( BEta_0.transpose(), D_mem, BEta_1 );
  Kqh_mem(0,1) = Scalar(0,0);

  Scalar = 0.0;
  Scalar += J_xixi * mc3.matmul ( BXi_0.transpose(), D_mem, BXi_2 );
  Scalar += J_xieta * mc3.matmul ( BXi_0.transpose(), D_mem, BEta_2 );
  Scalar += J_xieta * mc3.matmul ( BEta_0.transpose(), D_mem, BXi_2 );
  Scalar += J_etaeta * mc3.matmul ( BEta_0.transpose(), D_mem, BEta_2 );
  Kqh_mem(0,2) = Scalar(0,0);

  Scalar = 0.0;
  Scalar += J_xixi * mc3.matmul ( BXi_1.transpose(), D_mem, BXi_0 );
  Scalar += J_xieta * mc3.matmul ( BXi_1.transpose(), D_mem, BEta_0 );
  Scalar += J_xieta * mc3.matmul ( BEta_1.transpose(), D_mem, BXi_0 );
  Scalar += J_etaeta * mc3.matmul ( BEta_1.transpose(), D_mem, BEta_0 );
  Kqh_mem(1,0) = Scalar(0,0);

  Scalar = 0.0;
  Scalar += J_xixi * mc3.matmul ( BXi_1.transpose(), D_mem, BXi_1 );
  Scalar += J_xieta * mc3.matmul ( BXi_1.transpose(), D_mem, BEta_1 );
  Scalar += J_xieta * mc3.matmul ( BEta_1.transpose(), D_mem, BXi_1 );
  Scalar += J_etaeta * mc3.matmul ( BEta_1.transpose(), D_mem, BEta_1 );
  Kqh_mem(1,1) = Scalar(0,0);

  Scalar = 0.0;
  Scalar += J_xixi * mc3.matmul ( BXi_1.transpose(), D_mem, BXi_2 );
  Scalar += J_xieta * mc3.matmul ( BXi_1.transpose(), D_mem, BEta_2 );
  Scalar += J_xieta * mc3.matmul ( BEta_1.transpose(), D_mem, BXi_2 );
  Scalar += J_etaeta * mc3.matmul ( BEta_1.transpose(), D_mem, BEta_2 );
  Kqh_mem(1,2) = Scalar(0,0);

  Scalar = 0.0;
  Scalar += J_xixi * mc3.matmul ( BXi_2.transpose(), D_mem, BXi_0 );
  Scalar += J_xieta * mc3.matmul ( BXi_2.transpose(), D_mem, BEta_0 );
  Scalar += J_xieta * mc3.matmul ( BEta_2.transpose(), D_mem, BXi_0 );
  Scalar += J_etaeta * mc3.matmul ( BEta_2.transpose(), D_mem, BEta_0 );
  Kqh_mem(2,0) = Scalar(0,0);

  Scalar = 0.0;
  Scalar += J_xixi * mc3.matmul ( BXi_2.transpose(), D_mem, BXi_1 );
  Scalar += J_xieta * mc3.matmul ( BXi_2.transpose(), D_mem, BEta_1 );
  Scalar += J_xieta * mc3.matmul ( BEta_2.transpose(), D_mem, BXi_1 );
  Scalar += J_etaeta * mc3.matmul ( BEta_2.transpose(), D_mem, BEta_1 );
  Kqh_mem(2,1) = Scalar(0,0);

  Scalar = 0.0;
  Scalar += J_xixi * mc3.matmul ( BXi_2.transpose(), D_mem, BXi_2 );
  Scalar += J_xieta * mc3.matmul ( BXi_2.transpose(), D_mem, BEta_2 );
  Scalar += J_xieta * mc3.matmul ( BEta_2.transpose(), D_mem, BXi_2 );
  Scalar += J_etaeta * mc3.matmul ( BEta_2.transpose(), D_mem, BEta_2 );
  Kqh_mem(2,2) = Scalar(0,0);
}


//-----------------------------------------------------------------------
//   getTmat_
//-----------------------------------------------------------------------


void Shell6DofsBFModel::getTmat_

  ( Matrix&  Tmat,
    const Matrix&   xcoords ) const
  
{
  double l_12, l_23, l_31;
  l_12 = sqrt(pow(xcoords(0, 1) - xcoords(0, 0), 2) + pow(xcoords(1, 1) - xcoords(1, 0), 2));
  l_23 = sqrt(pow(xcoords(0, 2) - xcoords(0, 1), 2) + pow(xcoords(1, 2) - xcoords(1, 1), 2));
  l_31 = sqrt(pow(xcoords(0, 0) - xcoords(0, 2), 2) + pow(xcoords(1, 0) - xcoords(1, 2), 2));
  double cos_12, cos_23, cos_31;
  cos_12 = (xcoords(1, 1) - xcoords(1, 0))/l_12;
  cos_23 = (xcoords(1, 2) - xcoords(1, 1))/l_23; 
  cos_31 = (xcoords(1, 0) - xcoords(1, 2))/l_31;
  double sin_12, sin_23, sin_31;
  sin_12 = -(xcoords(0, 1) - xcoords(0, 0))/l_12;
  sin_23 = -(xcoords(0, 2) - xcoords(0, 1))/l_23;
  sin_31 = -(xcoords(0, 0) - xcoords(0, 2))/l_31;

  Tmat = 0.0;
  Tmat(0, 0) = 1.0;
  Tmat(1, 3) = 1.0;
  Tmat(2, 6) = 1.0;
  Tmat(3, 0) = l_12/2.0;
  Tmat(3, 1) = -(pow(l_12, 2)/12.0)*sin_12;
  Tmat(3, 2) = (pow(l_12, 2)/12.0)*cos_12;
  Tmat(3, 3) = l_12/2.0;
  Tmat(3, 4) = (pow(l_12, 2)/12.0)*sin_12;
  Tmat(3, 5) = -(pow(l_12, 2)/12.0)*cos_12;
  Tmat(4, 3) = l_23/2.0;
  Tmat(4, 4) = -(pow(l_23, 2)/12.0)*sin_23;
  Tmat(4, 5) = (pow(l_23, 2)/12.0)*cos_23;
  Tmat(4, 6) = l_23/2.0;
  Tmat(4, 7) = (pow(l_23, 2)/12.0)*sin_23;
  Tmat(4, 8) = -(pow(l_23, 2)/12.0)*cos_23;
  Tmat(5, 0) = l_31/2.0;
  Tmat(5, 1) = (pow(l_31, 2)/12.0)*sin_31;
  Tmat(5, 2) = -(pow(l_31, 2)/12.0)*cos_31;
  Tmat(5, 6) = l_31/2.0;
  Tmat(5, 7) = -(pow(l_31, 2)/12.0)*sin_31;
  Tmat(5, 8) = (pow(l_31, 2)/12.0)*cos_31;
  Tmat(6, 1) = -(l_12/3.0)*cos_12;
  Tmat(6, 2) = -(l_12/3.0)*sin_12;
  Tmat(6, 4) = -(l_12/6.0)*cos_12;
  Tmat(6, 5) = -(l_12/6.0)*sin_12;
  Tmat(7, 1) = -(l_12/6.0)*cos_12;
  Tmat(7, 2) = -(l_12/6.0)*sin_12;
  Tmat(7, 4) = -(l_12/3.0)*cos_12;
  Tmat(7, 5) = -(l_12/3.0)*sin_12;
  Tmat(8, 4) = -(l_23/3.0)*cos_23;
  Tmat(8, 5) = -(l_23/3.0)*sin_23;
  Tmat(8, 7) = -(l_23/6.0)*cos_23;
  Tmat(8, 8) = -(l_23/6.0)*sin_23;
  Tmat(9, 4) = -(l_23/6.0)*cos_23;
  Tmat(9, 5) = -(l_23/6.0)*sin_23;
  Tmat(9, 7) = -(l_23/3.0)*cos_23;
  Tmat(9, 8) = -(l_23/3.0)*sin_23;
  Tmat(10, 1) = -(l_31/6.0)*cos_31;
  Tmat(10, 2) = -(l_31/6.0)*sin_31;
  Tmat(10, 7) = -(l_31/3.0)*cos_31;
  Tmat(10, 8) = -(l_31/3.0)*sin_31;
  Tmat(11, 1) = -(l_31/3.0)*cos_31;
  Tmat(11, 2) = -(l_31/3.0)*sin_31;
  Tmat(11, 7) = -(l_31/6.0)*cos_31;
  Tmat(11, 8) = -(l_31/6.0)*sin_31;
}


//-----------------------------------------------------------------------
//   getMembraneDisp_
//-----------------------------------------------------------------------


void Shell6DofsBFModel::getMembraneDisp_

  ( Vector&         Disp_mem,
    const Vector&   elemxDisp )   const
    
{
  for ( idx_t i = 0; i < nodeCount_; ++i )
  {
    Disp_mem[2*i]   = elemxDisp[6*i];
    Disp_mem[2*i+1] = elemxDisp[6*i+1];
  }
}


//-----------------------------------------------------------------------
//   getHmatAnalytic_
//-----------------------------------------------------------------------


void Shell6DofsBFModel::getHmatAnalytic_

  ( Matrix&  Hmat,
    const Matrix&   D_plate,
    const Matrix&   xcoords ) const

{
  Matrix       H1mat      ( 7, 7);
  Matrix       H2mat      ( 7, 7);
  Matrix       H3mat      ( 7, 7);
  Matrix       H4mat      ( 7, 7);
  Matrix       H5mat      ( 7, 7);
  Matrix       H6mat      ( 7, 7);
  H1mat = 0.0;
  H2mat = 0.0;
  H3mat = 0.0;
  H4mat = 0.0;
  H5mat = 0.0;
  H6mat = 0.0;

  double I_dA, I_xdA, I_ydA, I_xxdA, I_xydA, I_yydA;
  I_dA = (1.0/2.0)*(xcoords(0, 0)*xcoords(1, 1) + xcoords(0, 1)*xcoords(1, 2) + xcoords(0, 2)*xcoords(1, 0) - xcoords(0, 0)*xcoords(1, 2)-xcoords(0, 1)*xcoords(1, 0)-xcoords(0, 2)*xcoords(1, 1));
  I_xdA = 0.0;
  I_ydA = 0.0;
  I_xxdA = (I_dA/12.0)*(pow(xcoords(0, 0), 2) + pow(xcoords(0, 1), 2) + pow(xcoords(0, 2), 2));
  I_xydA = (I_dA/12.0)*(xcoords(0, 0)*xcoords(1, 0) + xcoords(0, 1)*xcoords(1, 1) + xcoords(0, 2)*xcoords(1, 2));
  I_yydA = (I_dA/12.0)*(pow(xcoords(1, 0), 2) + pow(xcoords(1, 1), 2) + pow(xcoords(1, 2), 2));

  H1mat(0, 0) = D_plate(0, 0)*4.0*I_dA;
  H1mat(0, 3) = D_plate(0, 0)*12.0*I_xdA;
  H1mat(0, 4) = D_plate(0, 0)*4.0*I_ydA;
  H1mat(3, 0) = D_plate(0, 0)*12.0*I_xdA;
  H1mat(3, 3) = D_plate(0, 0)*36.0*I_xxdA;
  H1mat(3, 4) = D_plate(0, 0)*12.0*I_xydA;
  H1mat(4, 0) = D_plate(0, 0)*4.0*I_ydA;
  H1mat(4, 3) = D_plate(0, 0)*12.0*I_xydA;
  H1mat(4, 4) = D_plate(0, 0)*4.0*I_yydA;

  H2mat(2, 2) = D_plate(1, 1)*4.0*I_dA;
  H2mat(2, 5) = D_plate(1, 1)*4.0*I_xdA;
  H2mat(2, 6) = D_plate(1, 1)*12.0*I_ydA;
  H2mat(5, 2) = D_plate(1, 1)*4.0*I_xdA;
  H2mat(5, 5) = D_plate(1, 1)*4.0*I_xxdA;
  H2mat(5, 6) = D_plate(1, 1)*12.0*I_xydA;
  H2mat(6, 2) = D_plate(1, 1)*12.0*I_ydA;
  H2mat(6, 5) = D_plate(1, 1)*12.0*I_xydA;
  H2mat(6, 6) = D_plate(1, 1)*36.0*I_yydA;

  H3mat(0, 2) = D_plate(0, 1)*4.0*I_dA;
  H3mat(0, 5) = D_plate(0, 1)*4.0*I_xdA;
  H3mat(0, 6) = D_plate(0, 1)*12.0*I_ydA;
  H3mat(2, 0) = D_plate(0, 1)*4.0*I_dA;
  H3mat(2, 3) = D_plate(0, 1)*12.0*I_xdA;
  H3mat(2, 4) = D_plate(0, 1)*4.0*I_ydA;
  H3mat(3, 2) = D_plate(0, 1)*12.0*I_xdA;
  H3mat(3, 5) = D_plate(0, 1)*12.0*I_xxdA;
  H3mat(3, 6) = D_plate(0, 1)*36.0*I_xydA;
  H3mat(4, 2) = D_plate(0, 1)*4.0*I_ydA;
  H3mat(4, 5) = D_plate(0, 1)*4.0*I_xydA;
  H3mat(4, 6) = D_plate(0, 1)*12.0*I_yydA;
  H3mat(5, 0) = D_plate(0, 1)*4.0*I_xdA;
  H3mat(5, 3) = D_plate(0, 1)*12.0*I_xxdA;
  H3mat(5, 4) = D_plate(0, 1)*4.0*I_xydA;
  H3mat(6, 0) = D_plate(0, 1)*12.0*I_ydA;
  H3mat(6, 3) = D_plate(0, 1)*36.0*I_xydA;
  H3mat(6, 4) = D_plate(0, 1)*12.0*I_yydA;

  H4mat(0, 1) = D_plate(0, 2)*4.0*I_dA;
  H4mat(0, 4) = D_plate(0, 2)*8.0*I_xdA;
  H4mat(0, 5) = D_plate(0, 2)*8.0*I_ydA;
  H4mat(1, 0) = D_plate(0, 2)*4.0*I_dA;
  H4mat(1, 3) = D_plate(0, 2)*12.0*I_xdA;
  H4mat(1, 4) = D_plate(0, 2)*4.0*I_ydA;
  H4mat(3, 1) = D_plate(0, 2)*12.0*I_xdA;
  H4mat(3, 4) = D_plate(0, 2)*24.0*I_xxdA;
  H4mat(3, 5) = D_plate(0, 2)*24.0*I_xydA;
  H4mat(4, 0) = D_plate(0, 2)*8.0*I_xdA;
  H4mat(4, 1) = D_plate(0, 2)*4.0*I_ydA;
  H4mat(4, 3) = D_plate(0, 2)*24.0*I_xxdA;
  H4mat(4, 4) = D_plate(0, 2)*16.0*I_xydA;
  H4mat(4, 5) = D_plate(0, 2)*8.0*I_yydA;
  H4mat(5, 0) = D_plate(0, 2)*8.0*I_ydA;
  H4mat(5, 3) = D_plate(0, 2)*24.0*I_xydA;
  H4mat(5, 4) = D_plate(0, 2)*8.0*I_yydA;

  H5mat(1, 2) = D_plate(1, 2)*4.0*I_dA;
  H5mat(1, 5) = D_plate(1, 2)*4.0*I_xdA;
  H5mat(1, 6) = D_plate(1, 2)*12.0*I_ydA;
  H5mat(2, 1) = D_plate(1, 2)*4.0*I_dA;
  H5mat(2, 4) = D_plate(1, 2)*8.0*I_xdA;
  H5mat(2, 5) = D_plate(1, 2)*8.0*I_ydA;
  H5mat(4, 2) = D_plate(1, 2)*8.0*I_xdA;
  H5mat(4, 5) = D_plate(1, 2)*8.0*I_xxdA;
  H5mat(4, 6) = D_plate(1, 2)*24.0*I_xydA;
  H5mat(5, 1) = D_plate(1, 2)*4.0*I_xdA;
  H5mat(5, 2) = D_plate(1, 2)*8.0*I_ydA;
  H5mat(5, 4) = D_plate(1, 2)*8.0*I_xxdA;
  H5mat(5, 5) = D_plate(1, 2)*16.0*I_xydA;
  H5mat(5, 6) = D_plate(1, 2)*24.0*I_yydA;
  H5mat(6, 1) = D_plate(1, 2)*12.0*I_ydA;
  H5mat(6, 4) = D_plate(1, 2)*24.0*I_xydA;
  H5mat(6, 5) = D_plate(1, 2)*24.0*I_yydA;

  H6mat(1, 1) = D_plate(2, 2)*4.0*I_dA;
  H6mat(1, 4) = D_plate(2, 2)*8.0*I_xdA;
  H6mat(1, 5) = D_plate(2, 2)*8.0*I_ydA;
  H6mat(4, 1) = D_plate(2, 2)*8.0*I_xdA;
  H6mat(4, 4) = D_plate(2, 2)*16.0*I_xxdA;
  H6mat(4, 5) = D_plate(2, 2)*16.0*I_xydA;
  H6mat(5, 1) = D_plate(2, 2)*8.0*I_ydA;
  H6mat(5, 4) = D_plate(2, 2)*16.0*I_xydA;
  H6mat(5, 5) = D_plate(2, 2)*16.0*I_yydA;

  Hmat = H1mat + H2mat + H3mat + H4mat + H5mat + H6mat;
}


//-----------------------------------------------------------------------
//   getHmatNumeric_
//-----------------------------------------------------------------------


void Shell6DofsBFModel::getHmatNumeric_

  ( Matrix&  Hmat,
    const Matrix&   D_plate,
    const Matrix&   coords,
    const Matrix&   ipcoords,
    const Vector&   weights ) const

{
  Matrix       H1mat      ( 7, 7);
  Matrix       H2mat      ( 7, 7);
  Matrix       H3mat      ( 7, 7);
  Matrix       H4mat      ( 7, 7);
  Matrix       H5mat      ( 7, 7);
  Matrix       H6mat      ( 7, 7);
  H1mat = 0.0;
  H2mat = 0.0;
  H3mat = 0.0;
  H4mat = 0.0;
  H5mat = 0.0;
  H6mat = 0.0;

  double Xcg = (coords(0, 0) + coords(0, 1) + coords(0, 2))/3.0;
  double Ycg = (coords(1, 0) + coords(1, 1) + coords(1, 2))/3.0;

  for ( idx_t ip = 0; ip < ipcoords.size(1); ip++ )
  {
    double x = ipcoords(0, ip)-Xcg;
    double y = ipcoords(1, ip)-Ycg;
    double w = weights[ip];
  
    H1mat(0, 0) += w*(D_plate(0, 0)*4.0);
    H1mat(0, 3) += w*(D_plate(0, 0)*12.0*x);
    H1mat(0, 4) += w*(D_plate(0, 0)*4.0*y);
    H1mat(3, 0) += w*(D_plate(0, 0)*12.0*x);
    H1mat(3, 3) += w*(D_plate(0, 0)*36.0*x*x);
    H1mat(3, 4) += w*(D_plate(0, 0)*12.0*x*y);
    H1mat(4, 0) += w*(D_plate(0, 0)*4.0*y);
    H1mat(4, 3) += w*(D_plate(0, 0)*12.0*x*y);
    H1mat(4, 4) += w*(D_plate(0, 0)*4.0*y*y);

    H2mat(2, 2) += w*(D_plate(1, 1)*4.0);
    H2mat(2, 5) += w*(D_plate(1, 1)*4.0*x);
    H2mat(2, 6) += w*(D_plate(1, 1)*12.0*y);
    H2mat(5, 2) += w*(D_plate(1, 1)*4.0*x);
    H2mat(5, 5) += w*(D_plate(1, 1)*4.0*x*x);
    H2mat(5, 6) += w*(D_plate(1, 1)*12.0*x*y);
    H2mat(6, 2) += w*(D_plate(1, 1)*12.0*y);
    H2mat(6, 5) += w*(D_plate(1, 1)*12.0*x*y);
    H2mat(6, 6) += w*(D_plate(1, 1)*36.0*y*y);

    H3mat(0, 2) += w*(D_plate(0, 1)*4.0);
    H3mat(0, 5) += w*(D_plate(0, 1)*4.0*x);
    H3mat(0, 6) += w*(D_plate(0, 1)*12.0*y);
    H3mat(2, 0) += w*(D_plate(0, 1)*4.0);
    H3mat(2, 3) += w*(D_plate(0, 1)*12.0*x);
    H3mat(2, 4) += w*(D_plate(0, 1)*4.0*y);
    H3mat(3, 2) += w*(D_plate(0, 1)*12.0*x);
    H3mat(3, 5) += w*(D_plate(0, 1)*12.0*x*x);
    H3mat(3, 6) += w*(D_plate(0, 1)*36.0*x*y);
    H3mat(4, 2) += w*(D_plate(0, 1)*4.0*y);
    H3mat(4, 5) += w*(D_plate(0, 1)*4.0*x*y);
    H3mat(4, 6) += w*(D_plate(0, 1)*12.0*y*y);
    H3mat(5, 0) += w*(D_plate(0, 1)*4.0*x);
    H3mat(5, 3) += w*(D_plate(0, 1)*12.0*x*x);
    H3mat(5, 4) += w*(D_plate(0, 1)*4.0*x*y);
    H3mat(6, 0) += w*(D_plate(0, 1)*12.0*y);
    H3mat(6, 3) += w*(D_plate(0, 1)*36.0*x*y);
    H3mat(6, 4) += w*(D_plate(0, 1)*12.0*y*y);

    H4mat(0, 1) += w*(D_plate(0, 2)*4.0);
    H4mat(0, 4) += w*(D_plate(0, 2)*8.0*x);
    H4mat(0, 5) += w*(D_plate(0, 2)*8.0*y);
    H4mat(1, 0) += w*(D_plate(0, 2)*4.0);
    H4mat(1, 3) += w*(D_plate(0, 2)*12.0*x);
    H4mat(1, 4) += w*(D_plate(0, 2)*4.0*y);
    H4mat(3, 1) += w*(D_plate(0, 2)*12.0*x);
    H4mat(3, 4) += w*(D_plate(0, 2)*24.0*x*x);
    H4mat(3, 5) += w*(D_plate(0, 2)*24.0*x*y);
    H4mat(4, 0) += w*(D_plate(0, 2)*8.0*x);
    H4mat(4, 1) += w*(D_plate(0, 2)*4.0*y);
    H4mat(4, 3) += w*(D_plate(0, 2)*24.0*x*x);
    H4mat(4, 4) += w*(D_plate(0, 2)*16.0*x*y);
    H4mat(4, 5) += w*(D_plate(0, 2)*8.0*y*y);
    H4mat(5, 0) += w*(D_plate(0, 2)*8.0*y);
    H4mat(5, 3) += w*(D_plate(0, 2)*24.0*x*y);
    H4mat(5, 4) += w*(D_plate(0, 2)*8.0*y*y);

    H5mat(1, 2) += w*(D_plate(1, 2)*4.0);
    H5mat(1, 5) += w*(D_plate(1, 2)*4.0*x);
    H5mat(1, 6) += w*(D_plate(1, 2)*12.0*y);
    H5mat(2, 1) += w*(D_plate(1, 2)*4.0);
    H5mat(2, 4) += w*(D_plate(1, 2)*8.0*x);
    H5mat(2, 5) += w*(D_plate(1, 2)*8.0*y);
    H5mat(4, 2) += w*(D_plate(1, 2)*8.0*x);
    H5mat(4, 5) += w*(D_plate(1, 2)*8.0*x*x);
    H5mat(4, 6) += w*(D_plate(1, 2)*24.0*x*y);
    H5mat(5, 1) += w*(D_plate(1, 2)*4.0*x);
    H5mat(5, 2) += w*(D_plate(1, 2)*8.0*y);
    H5mat(5, 4) += w*(D_plate(1, 2)*8.0*x*x);
    H5mat(5, 5) += w*(D_plate(1, 2)*16.0*x*y);
    H5mat(5, 6) += w*(D_plate(1, 2)*24.0*y*y);
    H5mat(6, 1) += w*(D_plate(1, 2)*12.0*y);
    H5mat(6, 4) += w*(D_plate(1, 2)*24.0*x*y);
    H5mat(6, 5) += w*(D_plate(1, 2)*24.0*y*y);

    H6mat(1, 1) += w*(D_plate(2, 2)*4.0);
    H6mat(1, 4) += w*(D_plate(2, 2)*8.0*x);
    H6mat(1, 5) += w*(D_plate(2, 2)*8.0*y);
    H6mat(4, 1) += w*(D_plate(2, 2)*8.0*x);
    H6mat(4, 4) += w*(D_plate(2, 2)*16.0*x*x);
    H6mat(4, 5) += w*(D_plate(2, 2)*16.0*x*y);
    H6mat(5, 1) += w*(D_plate(2, 2)*8.0*y);
    H6mat(5, 4) += w*(D_plate(2, 2)*16.0*x*y);
    H6mat(5, 5) += w*(D_plate(2, 2)*16.0*y*y);
  }   

  Hmat = H1mat + H2mat + H3mat + H4mat + H5mat + H6mat;
}


//-----------------------------------------------------------------------
//   getBmat_
//-----------------------------------------------------------------------


void Shell6DofsBFModel::getBmat_

  ( Matrix&  Bplate,
    const Matrix&   Dplate,
    const Matrix&   xcoords ) const

{
  using jem::ALL;

  // Computing trigonometric terms.

  double l_12, l_23, l_31;
  l_12 = sqrt(pow(xcoords(0, 1) - xcoords(0, 0), 2) + pow(xcoords(1, 1) - xcoords(1, 0), 2));
  l_23 = sqrt(pow(xcoords(0, 2) - xcoords(0, 1), 2) + pow(xcoords(1, 2) - xcoords(1, 1), 2));
  l_31 = sqrt(pow(xcoords(0, 0) - xcoords(0, 2), 2) + pow(xcoords(1, 0) - xcoords(1, 2), 2));
  double cos_12, cos_23, cos_31;
  cos_12 = (xcoords(1, 1) - xcoords(1, 0))/l_12;
  cos_23 = (xcoords(1, 2) - xcoords(1, 1))/l_23; 
  cos_31 = (xcoords(1, 0) - xcoords(1, 2))/l_31;
  double sin_12, sin_23, sin_31;
  sin_12 = -(xcoords(0, 1) - xcoords(0, 0))/l_12;
  sin_23 = -(xcoords(0, 2) - xcoords(0, 1))/l_23;
  sin_31 = -(xcoords(0, 0) - xcoords(0, 2))/l_31;
  double sin_2_12, sin_2_23, sin_2_31;
  sin_2_12 = 2.0*sin_12*cos_12;
  sin_2_23 = 2.0*sin_23*cos_23;
  sin_2_31 = 2.0*sin_31*cos_31; 
  double cos_2_12, cos_2_23, cos_2_31;
  cos_2_12 = 2.0*cos_12*cos_12 - 1;
  cos_2_23 = 2.0*cos_23*cos_23 - 1; 
  cos_2_31 = 2.0*cos_31*cos_31 - 1;
  double S_1, S_2, S_3;
  S_1 = sin_2_12 - sin_2_31;
  S_2 = sin_2_23 - sin_2_12;
  S_3 = sin_2_31 - sin_2_23;
  double C_1, C_2, C_3;
  C_1 = cos_2_12 - cos_2_31;
  C_2 = cos_2_23 - cos_2_12;
  C_3 = cos_2_31 - cos_2_23;

 // Computing The BMx, BMy and BMxy Matrices.

  Matrix  BMx  ( 7, xcoords.size(1) );
  BMx(0, ALL) = -2.0*Dplate(0, 0);
  BMx(1, ALL) = -2.0*Dplate(0, 2);
  BMx(2, ALL) = -2.0*Dplate(0, 1);
  BMx(3, ALL) = -6.0*xcoords(0, ALL)*Dplate(0, 0);
  BMx(4, ALL) = -(2.0*xcoords(1, ALL)*Dplate(0, 0)+4.0*xcoords(0, ALL)*Dplate(0, 2));
  BMx(5, ALL) = -(2.0*xcoords(0, ALL)*Dplate(0, 1)+4.0*xcoords(1, ALL)*Dplate(0, 2));
  BMx(6, ALL) = -6.0*xcoords(1, ALL)*Dplate(0, 1);

  Matrix  BMy  ( 7, xcoords.size(1) );
  BMy(0, ALL) = -2.0*Dplate(0, 1);
  BMy(1, ALL) = -2.0*Dplate(1, 2);
  BMy(2, ALL) = -2.0*Dplate(1, 1);
  BMy(3, ALL) = -6.0*xcoords(0, ALL)*Dplate(0, 1);
  BMy(4, ALL) = -(2.0*xcoords(1, ALL)*Dplate(0, 1)+4.0*xcoords(0, ALL)*Dplate(1, 2));
  BMy(5, ALL) = -(2.0*xcoords(0, ALL)*Dplate(1, 1)+4.0*xcoords(1, ALL)*Dplate(1, 2));
  BMy(6, ALL) = -6.0*xcoords(1, ALL)*Dplate(1, 1);

  Matrix  BMxy ( 7, xcoords.size(1) );
  BMxy(0, ALL) = -2.0*Dplate(0, 2);
  BMxy(1, ALL) = -2.0*Dplate(2, 2);
  BMxy(2, ALL) = -2.0*Dplate(1, 2);
  BMxy(3, ALL) = -6.0*xcoords(0, ALL)*Dplate(0, 2);
  BMxy(4, ALL) = -(4.0*xcoords(0, ALL)*Dplate(2, 2)+2.0*xcoords(1, ALL)*Dplate(0, 2));
  BMxy(5, ALL) = -(4.0*xcoords(1, ALL)*Dplate(2, 2)+2.0*xcoords(0, ALL)*Dplate(1, 2));
  BMxy(6, ALL) = -6.0*xcoords(1, ALL)*Dplate(1, 2);

  // Computing the derivatives of BMx, BMy and BMxy Matrices.

  Matrix  dBMx_dx  ( 7, xcoords.size(1) );
  Matrix  dBMy_dx  ( 7, xcoords.size(1) );
  Matrix  dBMxy_dx ( 7, xcoords.size(1) );
  dBMx_dx = 0.0;
  dBMy_dx = 0.0;
  dBMxy_dx = 0.0;
  dBMx_dx(3, ALL) = -6.0*Dplate(0, 0);
  dBMx_dx(4, ALL) = -4.0*Dplate(0, 2);
  dBMx_dx(5, ALL) = -2.0*Dplate(0, 1);
  dBMy_dx(3, ALL) = -6.0*Dplate(0, 1);
  dBMy_dx(4, ALL) = -4.0*Dplate(1, 2);
  dBMy_dx(5, ALL) = -2.0*Dplate(1, 1);
  dBMxy_dx(3,ALL) = -6.0*Dplate(0, 2);
  dBMxy_dx(4,ALL) = -4.0*Dplate(2, 2);
  dBMxy_dx(5,ALL) = -2.0*Dplate(1, 2);

  Matrix  dBMx_dy  ( 7, xcoords.size(1) );
  Matrix  dBMy_dy  ( 7, xcoords.size(1) );
  Matrix  dBMxy_dy ( 7, xcoords.size(1) );
  dBMx_dy = 0.0;
  dBMy_dy = 0.0;
  dBMxy_dy = 0.0;
  dBMx_dy(4, ALL) = -2.0*Dplate(0, 0);
  dBMx_dy(5, ALL) = -4.0*Dplate(0, 2);
  dBMx_dy(6, ALL) = -6.0*Dplate(0, 1);
  dBMy_dy(4, ALL) = -2.0*Dplate(0, 1);
  dBMy_dy(5, ALL) = -4.0*Dplate(1, 2);
  dBMy_dy(6, ALL) = -6.0*Dplate(1, 1);
  dBMxy_dy(4,ALL) = -2.0*Dplate(0, 2);
  dBMxy_dy(5,ALL) = -4.0*Dplate(2, 2);
  dBMxy_dy(6,ALL) = -6.0*Dplate(1, 2);

  // Computing the BR_i Matrices.

  Matrix  BR_1  ( 7, 1 );
  Matrix  BR_2  ( 7, 1 );
  Matrix  BR_3  ( 7, 1 );
  BR_1(ALL, 0) = (1.0/2.0)*S_1*(BMy(ALL,0)-BMx(ALL,0)) + C_1*BMxy(ALL,0);
  BR_2(ALL, 0) = (1.0/2.0)*S_2*(BMy(ALL,1)-BMx(ALL,1)) + C_2*BMxy(ALL,1);
  BR_3(ALL, 0) = (1.0/2.0)*S_3*(BMy(ALL,2)-BMx(ALL,2)) + C_3*BMxy(ALL,2);

  // Computing the BMn_ij Matrices.

  Matrix  BMn_12  ( 7, 1 );
  Matrix  BMn_21  ( 7, 1 );
  Matrix  BMn_23  ( 7, 1 );
  Matrix  BMn_32  ( 7, 1 );
  Matrix  BMn_31  ( 7, 1 );
  Matrix  BMn_13  ( 7, 1 );
  BMn_12(ALL, 0) = cos_12*cos_12*BMx(ALL,0)+sin_12*sin_12*BMy(ALL,0)+sin_2_12*BMxy(ALL,0);
  BMn_21(ALL, 0) = cos_12*cos_12*BMx(ALL,1)+sin_12*sin_12*BMy(ALL,1)+sin_2_12*BMxy(ALL,1);
  BMn_23(ALL, 0) = cos_23*cos_23*BMx(ALL,1)+sin_23*sin_23*BMy(ALL,1)+sin_2_23*BMxy(ALL,1);
  BMn_32(ALL, 0) = cos_23*cos_23*BMx(ALL,2)+sin_23*sin_23*BMy(ALL,2)+sin_2_23*BMxy(ALL,2);
  BMn_31(ALL, 0) = cos_31*cos_31*BMx(ALL,2)+sin_31*sin_31*BMy(ALL,2)+sin_2_31*BMxy(ALL,2);
  BMn_13(ALL, 0) = cos_31*cos_31*BMx(ALL,0)+sin_31*sin_31*BMy(ALL,0)+sin_2_31*BMxy(ALL,0);

  // Computing the BMnx_ij, BMny_ij, BMnsx_ij, BMnsy_ij  Matrices.

  Matrix  BMnx_12  ( 7, 1 );
  Matrix  BMnx_23  ( 7, 1 );
  Matrix  BMnx_31  ( 7, 1 );
  BMnx_12(ALL, 0) = cos_12*cos_12*dBMx_dx(ALL,0)+sin_12*sin_12*dBMy_dx(ALL,0)+sin_2_12*dBMxy_dx(ALL,0);
  BMnx_23(ALL, 0) = cos_23*cos_23*dBMx_dx(ALL,1)+sin_23*sin_23*dBMy_dx(ALL,1)+sin_2_23*dBMxy_dx(ALL,1);
  BMnx_31(ALL, 0) = cos_31*cos_31*dBMx_dx(ALL,2)+sin_31*sin_31*dBMy_dx(ALL,2)+sin_2_31*dBMxy_dx(ALL,2);

  Matrix  BMny_12  ( 7, 1 );
  Matrix  BMny_23  ( 7, 1 );
  Matrix  BMny_31  ( 7, 1 );
  BMny_12(ALL, 0) = cos_12*cos_12*dBMx_dy(ALL,0)+sin_12*sin_12*dBMy_dy(ALL,0)+sin_2_12*dBMxy_dy(ALL,0);
  BMny_23(ALL, 0) = cos_23*cos_23*dBMx_dy(ALL,1)+sin_23*sin_23*dBMy_dy(ALL,1)+sin_2_23*dBMxy_dy(ALL,1);
  BMny_31(ALL, 0) = cos_31*cos_31*dBMx_dy(ALL,2)+sin_31*sin_31*dBMy_dy(ALL,2)+sin_2_31*dBMxy_dy(ALL,2);

  Matrix  BMnsx_12  ( 7, 1 );
  Matrix  BMnsx_23  ( 7, 1 );
  Matrix  BMnsx_31  ( 7, 1 );
  BMnsx_12(ALL, 0) = 0.5*sin_2_12*(dBMy_dx(ALL,0)-dBMx_dx(ALL,0))+cos_2_12*dBMxy_dx(ALL,0);
  BMnsx_23(ALL, 0) = 0.5*sin_2_23*(dBMy_dx(ALL,1)-dBMx_dx(ALL,1))+cos_2_23*dBMxy_dx(ALL,1);
  BMnsx_31(ALL, 0) = 0.5*sin_2_31*(dBMy_dx(ALL,2)-dBMx_dx(ALL,2))+cos_2_31*dBMxy_dx(ALL,2);

  Matrix  BMnsy_12  ( 7, 1 );
  Matrix  BMnsy_23  ( 7, 1 );
  Matrix  BMnsy_31  ( 7, 1 );
  BMnsy_12(ALL, 0) = 0.5*sin_2_12*(dBMy_dy(ALL,0)-dBMx_dy(ALL,0))+cos_2_12*dBMxy_dy(ALL,0);
  BMnsy_23(ALL, 0) = 0.5*sin_2_23*(dBMy_dy(ALL,1)-dBMx_dy(ALL,1))+cos_2_23*dBMxy_dy(ALL,1);
  BMnsy_31(ALL, 0) = 0.5*sin_2_31*(dBMy_dy(ALL,2)-dBMx_dy(ALL,2))+cos_2_31*dBMxy_dy(ALL,2);
    
  // Computing the BVn_ij  Matrices.

  Matrix  BVn_12  ( 7, 1 );
  Matrix  BVn_23  ( 7, 1 );
  Matrix  BVn_31  ( 7, 1 );
  BVn_12 = cos_12*BMnx_12+sin_12*BMny_12-2.0*sin_12*BMnsx_12+2.0*cos_12*BMnsy_12;
  BVn_23 = cos_23*BMnx_23+sin_23*BMny_23-2.0*sin_23*BMnsx_23+2.0*cos_23*BMnsy_23;
  BVn_31 = cos_31*BMnx_31+sin_31*BMny_31-2.0*sin_31*BMnsx_31+2.0*cos_31*BMnsy_31;

 // Computing the Bplate Matrix.

  Bplate = 0.0;
  Bplate(ALL,0) =  BR_1(ALL,0);
  Bplate(ALL,1) =  BR_2(ALL,0);
  Bplate(ALL,2) =  BR_3(ALL,0);
  Bplate(ALL,3) =  BVn_12(ALL,0);
  Bplate(ALL,4) =  BVn_23(ALL,0);
  Bplate(ALL,5) =  BVn_31(ALL,0);
  Bplate(ALL,6) =  BMn_12(ALL,0);
  Bplate(ALL,7) =  BMn_21(ALL,0);
  Bplate(ALL,8) =  BMn_23(ALL,0);
  Bplate(ALL,9) =  BMn_32(ALL,0);
  Bplate(ALL,10) = BMn_31(ALL,0);
  Bplate(ALL,11) = BMn_13(ALL,0);
}


//-----------------------------------------------------------------------
//   getDissForce_       NOT DONE FOR THE SHELL MODEL YET
//-----------------------------------------------------------------------


void Shell6DofsBFModel::getDissForce_

  ( const Ref<Plasticity> p,
    const Vector&         fstar,
    const Vector&         disp )   const

{
  Matrix      coords      ( nodeRank_, nodeCount_ );
  Matrix      transMat    ( nodeRank_, nodeRank_ );
  Matrix      transMatDof ( dofCount_, dofCount_ );
  Matrix      xcoords     ( rank_, nodeCount_ );
  Matrix      b           ( strCount_, dofCount_  );
  Matrix      bt          = b.transpose ();
  idx_t       ipoint      = 0;

  Cubix       grads       ( rank_, nodeCount_, ipCount_ );
  Vector      ipWeights   ( ipCount_ );

  IdxVector   inodes      ( nodeCount_ );
  IdxVector   idofs       ( dofCount_  );
  Vector      strain      ( strCount_  );
  Vector      elemForce   ( dofCount_  );
  Vector      elemDisp    ( dofCount_  );
  Vector      sstar       ( strCount_  );

  for ( idx_t ie = 0; ie < numElem_; ++ie )
  {
    idx_t ielem = ielems_[ie];

    elems_.getElemNodes  ( inodes, ielem             );
    nodes_.getSomeCoords ( coords, inodes            );
    dofs_->getDofIndices ( idofs , inodes, dofTypes_ );

    Vector v ( nodeRank_ );
    for ( idx_t i = 0; i < nodeRank_; i++ )
    {
      v[i] = orientVec_[i];
    }

    getTransMatrix_ ( transMatDof, transMat, coords, v );

    get2DLocalcoordinates_(xcoords, transMat, coords);

    shape_->getShapeGradients ( grads, ipWeights, xcoords );

    elemDisp   = select ( disp, idofs );

    ipWeights *= thickness_;

    elemForce = 0.0;

    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    {
      // Compute the B-matrix for this integration point.

      getShapeGrads_ ( b, grads(ALL,ALL,ip) );

      // Compute the strain vector of this integration point

      matmul ( strain, b, elemDisp );

      // get dissipation stress F^T*sigma+D^T*eps^p

      p->getDissipationStress ( sstar, strain, ipoint++ );

      elemForce += ipWeights[ip] * matmul ( bt, sstar );
    }
    // Add the element force vector to the global force vector.

    fstar[idofs] += elemForce;
  }
}


//-----------------------------------------------------------------------
//   initCharLength_       NOT DONE FOR THE SHELL MODEL YET
//-----------------------------------------------------------------------

void Shell6DofsBFModel::initCharLength_ ()

{
  IdxVector   inodes      (            nodeCount_ );
  Matrix      coords      ( nodeRank_,     nodeCount_ );
  Matrix      transMat    ( nodeRank_, nodeRank_ );
  Matrix      transMatDof ( dofCount_, dofCount_ );
  Matrix      xcoords     ( rank_, nodeCount_ );

  // maimi07i for triangles
  // double      fac   = 2. / sqrt ( sqrt(3.) ); 

  // my own expression 6/pi/3^.25
  double      fac   = 1.4512;

  double      maxLe = softening_->maxAllowedLength();

  idx_t       ielem;

  charLength_.resize ( max(ielems_)+1 );

  for ( idx_t ie = 0; ie < numElem_; ++ie )
  {
    ielem   = ielems_[ie];

    Vector    ipWeights ( ipCount_ );

    elems_.getElemNodes  ( inodes, ielem  );
    nodes_.getSomeCoords ( coords, inodes );

    Vector v ( nodeRank_ );
    for ( idx_t i = 0; i < nodeRank_; i++ )
    {
      v[i] = orientVec_[i];
    }

    getTransMatrix_ ( transMatDof, transMat, coords, v );

    get2DLocalcoordinates_(xcoords, transMat, coords);

    shape_->getIntegrationWeights ( ipWeights, xcoords );

    double area = sum ( ipWeights );
    double le   = fac * sqrt ( area );

    if ( le > maxLe )
    {
      System::warn() << "characteristic length of element " <<
        ielem << " results in local snapback!\n" <<
        "changed from " << le << " to " << maxLe  << endl;

      charLength_[ielem] = maxLe;
    }
    else
    {
      charLength_[ielem] = le;
    }
  }
}


//-----------------------------------------------------------------------
//   getDissipation_       NOT DONE FOR THE SHELL MODEL YET
//-----------------------------------------------------------------------

void Shell6DofsBFModel::getDissipation_

  ( const Properties&  params )

{
  const idx_t  nodeCount  = shape_->nodeCount   ();
  const idx_t  ipCount    = shape_->ipointCount ();
  const idx_t  ielemCount = ielems_.size        ();

  IdxVector    inodes      (        nodeCount );
  Matrix       coords      ( nodeRank_, nodeCount );
  Matrix       transMat    ( nodeRank_, nodeRank_ );
  Matrix       transMatDof ( dofCount_, dofCount_ );
  Matrix       xcoords     ( rank_, nodeCount_ );
  Vector       ipWeights   (          ipCount );

  idx_t  ipoint      = 0;
  double dissipation = 0.;

  // bulk damage

  for ( idx_t ie = 0; ie < ielemCount; ++ie )
  {
    idx_t ielem = ielems_[ie];

    elems_.getElemNodes  ( inodes, ielem   );
    nodes_.getSomeCoords ( coords, inodes );

    Vector v ( nodeRank_ );
    for ( idx_t i = 0; i < nodeRank_; i++ )
    {
      v[i] = orientVec_[i];
    }

    getTransMatrix_ ( transMatDof, transMat, coords, v );

    get2DLocalcoordinates_(xcoords, transMat, coords);

    // get the correct shape and then the number of idx_t points

    shape_->getIntegrationWeights ( ipWeights, xcoords );

    for ( idx_t ip = 0; ip < ipCount; ++ip )
    {
      dissipation += ipWeights[ip] * material_->giveDissipation ( ipoint++ );
    }
  }
  params.set ( myTag_, dissipation );
}


//-----------------------------------------------------------------------
//   getIShape_
//-----------------------------------------------------------------------

void Shell6DofsBFModel::getIshape_ 

  ( Ref<IShape> shape ) const

{
  shape = shape_;
}

//-----------------------------------------------------------------------
//   getMaterial_
//-----------------------------------------------------------------------

void Shell6DofsBFModel::getMaterial_

  ( Ref<Material> material ) const 

{
  material = material_;
}

//-----------------------------------------------------------------------
//   getThickness_
//-----------------------------------------------------------------------

void Shell6DofsBFModel::getThickness_

  ( double thickness ) const 

{
  thickness = thickness_;
}

//-----------------------------------------------------------------------
//   getOrientVec_
//-----------------------------------------------------------------------

void Shell6DofsBFModel::getOrientVec_

  ( const double orientVec[3] ) const 

{
  orientVec = orientVec_;
}

//-----------------------------------------------------------------------
//   getStress_       NOT DONE FOR THE SHELL MODEL YET
//-----------------------------------------------------------------------


void Shell6DofsBFModel::getStress_

  ( XTable&        table,
    const Vector&  weights,
    const Vector&  disp )

{
  IdxVector   ielems     = egroup_.getIndices  ();

  Matrix     ndNStress   ( nodeCount_, strCount_ );  // nodal normal stress
  Vector     ndWeights   ( nodeCount_ );
  Matrix     stiff       ( strCount_,  strCount_ );

  Cubix      grads       ( rank_, nodeCount_, ipCount_ );
  Matrix     coords      ( nodeRank_,     nodeCount_ );
  Matrix     transMat    ( nodeRank_, nodeRank_ );
  Matrix     transMatDof ( dofCount_, dofCount_ );
  Matrix     xcoords     ( rank_,     nodeCount_ );
  Matrix     b           ( strCount_, dofCount_  );

  Vector     nStressIp   ( strCount_ );    // normal stress vector at idx_t.pt.
  Vector     strain      ( strCount_ );
  Vector     elemDisp    ( dofCount_ );

  IdxVector  inodes      ( nodeCount_ );
  IdxVector  idofs       ( dofCount_  );
  IdxVector  jcols       ( strCount_  );

  jcols.resize ( strCount_ );

  // Add the columns for the stress components to the table.

  switch ( strCount_ )
  {
  case 1:

    jcols[0] = table.addColumn ( "xx" );

    break;

  case 3:

    jcols[0] = table.addColumn ( "xx" );
    jcols[1] = table.addColumn ( "yy" );
    jcols[2] = table.addColumn ( "xy" );

    break;

  case 6:

    jcols[0] = table.addColumn ( "xx" );
    jcols[1] = table.addColumn ( "yy" );
    jcols[2] = table.addColumn ( "zz" );
    jcols[3] = table.addColumn ( "xy" );
    jcols[4] = table.addColumn ( "yz" );
    jcols[5] = table.addColumn ( "xz" );

    break;

  default:

    throw Error (
      JEM_FUNC,
      "unexpected number of stress components: " +
      String ( strCount_ )
    );
  }

  idx_t         ipoint = 0;

  Vector      ipWeights ( ipCount_ );

  for ( idx_t ie = 0; ie < numElem_; ie++ )
  {
    // Get the global element index.

    idx_t  ielem = ielems[ie];

    ndNStress  = 0.0;
    ndWeights  = 0.0;

    elems_.getElemNodes  ( inodes, ielem );
    dofs_->getDofIndices ( idofs,  inodes,  dofTypes_ );

    nodes_.getSomeCoords ( coords, inodes );

    Vector v ( nodeRank_ );
    for ( idx_t i = 0; i < nodeRank_; i++ )
    {
      v[i] = orientVec_[i];
    }

    getTransMatrix_ ( transMatDof, transMat, coords, v );

    get2DLocalcoordinates_(xcoords, transMat, coords);

    shape_->getShapeGradients ( grads, ipWeights, xcoords );

    elemDisp = select ( disp, idofs );

    Matrix     sfuncs     = shape_->getShapeFunctions ();

    // Iterate over the integration points.

    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    {
      getShapeGrads_ ( b, grads(ALL,ALL,ip) );
      matmul ( strain, b, elemDisp );

      if ( crackBandMethod_ )
      {
        double le = charLength_[ielem];
        softening_->update ( nStressIp, stiff, strain, ipoint, le );
      }
      else
      {
        material_->update ( nStressIp, stiff, strain, ipoint );
      }

      ndNStress += matmul ( sfuncs(ALL,ip), nStressIp );

      ndWeights += sfuncs(ALL,ip);

      ++ipoint; 
    }

    select ( weights, inodes ) += ndWeights;

    // Add the stresses to the table.

    table.addBlock ( inodes, jcols[slice(0,strCount_)],   ndNStress );
  }
}


//-----------------------------------------------------------------------
//    writeElements_       NOT DONE FOR THE SHELL MODEL YET
//-----------------------------------------------------------------------

void Shell6DofsBFModel::writeElements_

  ( const Properties&  params,
    const Properties&  globdat )

{
  IdxVector  inodes ( nodeCount_ );

  if ( elemOut_ == NIL )
  {
    // open file and write connectivity of all elements

    elemOut_ = initWriter_ ( params, "elems" );

    // write connectivity of elements

    for ( idx_t ie = 0; ie < numElem_; ++ie )
    {
      idx_t ielem = ielems_[ie];

      elems_.getElemNodes ( inodes, ielem );

      *elemOut_ << ielem << " ";

      for ( idx_t in = 0; in < nodeCount_; ++in )
      {
        *elemOut_ << inodes[in] << ' ';
      }
      *elemOut_ << '\n';
    }
  }
  elemOut_->flush();
}


//-----------------------------------------------------------------------
//    writeNodes_       NOT DONE FOR THE SHELL MODEL YET
//-----------------------------------------------------------------------

void Shell6DofsBFModel::writeNodes_

  ( const Properties&  params,
    const Properties&  globdat )

{
  if ( nodeOut_ == NIL )
  {
    // all nodes are written in one file, nodeOut_ is static

    String       prepend;
    String       filename;

    if ( params.find( prepend, "prepend" ) )
    {
      filename = prepend + String(".all.nodes");
    }
    else
    {
      filename = "all.nodes";
    }
    nodeOut_ = newInstance<PrintWriter>(
               newInstance<FileWriter> ( filename ) );
  }

  Vector      pCoords    ( nodeRank_ );

  // quasi nodes

  while ( nodesWritten_ < nodes_.size() )
  {
    nodes_.getNodeCoords( pCoords, nodesWritten_ );

    *nodeOut_ << nodesWritten_++ <<  " ";

    for ( idx_t j = 0; j < nodeRank_; ++j )
    {
      *nodeOut_ << pCoords[j] << " ";
    }

    *nodeOut_ << "\n";
  }

  nodeOut_->flush();
}


//-----------------------------------------------------------------------
//    initLocOutWriter_       NOT DONE FOR THE SHELL MODEL YET
//-----------------------------------------------------------------------

void  Shell6DofsBFModel::initLocOutWriter_ ()

{
  Matrix      coords     ( nodeRank_, nodeCount_ );
  Matrix      ipcoords   ( nodeRank_, ipCount_ );
  IdxVector   inodes     ( nodeCount_ );

  // open file

  locOut_ = initWriter_ ( Properties(), "locOut" );

  // write coordinate data

  *locOut_ << "number of points, dim:\n" << locIes_.size() << " " << rank_ 
    << "\ncoords of points for which history is written:\n";

  for ( idx_t i = 0; i < locIes_.size(); i++ )
  {
    idx_t ie = locIes_[i];
    idx_t ip = locIpNumbers_[i];
    idx_t ielem = ielems_[ie];

    elems_.getElemNodes  ( inodes, ielem    );
    nodes_.getSomeCoords ( coords, inodes );
    shape_->getGlobalIntegrationPoints ( ipcoords, coords );

    for ( idx_t j = 0; j < nodeRank_; ++j )
    {
      locOut_->printFloat ( ipcoords(j,ip) ); 
      locOut_->printSpace ();
    }
    locOut_->printLine ();
  }
}


//-----------------------------------------------------------------------
//    initWriter_       NOT DONE FOR THE SHELL MODEL YET
//-----------------------------------------------------------------------

Ref<PrintWriter>  Shell6DofsBFModel::initWriter_

  ( const Properties&  params, 
    const String       name )  const

{
  // Open file for output

  StringVector fileName;
  String       prepend;

  if ( params.find( prepend, "prepend" ) )
  { 
    fileName.resize(3);

    fileName[0] = prepend;
    fileName[1] = myTag_;
    fileName[2] = name;
  }
  else
  {
    fileName.resize(2);

    fileName[0] = myTag_;
    fileName[1] = name;
  }

  return newInstance<PrintWriter>( newInstance<FileWriter> ( 
         StringUtils::join( fileName, "." ) ) );
}


//-----------------------------------------------------------------------
//    writeDisplacements_       NOT DONE FOR THE SHELL MODEL YET
//-----------------------------------------------------------------------


void Shell6DofsBFModel::writeDisplacements_

  ( const Properties&  params,
    const Properties&  globdat )

{
  Vector      dispi      ( 5 );
  IdxVector   idofs      ( 5 );
  Vector      disp;
  idx_t         it;

  if ( dispOut_ == NIL )
  {
    dispOut_ = initWriter_ ( params, "disp" );
  }

  globdat.get ( it, Globdat::TIME_STEP );

  *dispOut_ << "newXOutput " << it << '\n';

  StateVector::get ( disp, dofs_, globdat );

  // regular nodes

  for ( idx_t inode = 0; inode < nodes_.size(); ++inode )
  {
    dofs_->getDofsForItem ( idofs,  dofTypes_,  inode );
    dispi = select ( disp, idofs );

    *dispOut_ << inode << " ";

    for ( idx_t j = 0; j < rank_; ++j )
    {
      *dispOut_ << dispi[j] << " ";
    }
    *dispOut_ << '\n';
  }
  dispOut_->flush();
}


//-----------------------------------------------------------------------
//   getTable_       NOT DONE FOR THE SHELL MODEL YET
//-----------------------------------------------------------------------


bool Shell6DofsBFModel::getTable_

  ( const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;

  String       contents;
  Ref<XTable>  table;
  Vector       weights;
  String       name;

  Vector       disp;

  StateVector::get ( disp, dofs_, globdat );

  // Get the table, the name of the table, and the table row weights
  // from the action parameters.

  params.get ( table,   ActionParams::TABLE );
  params.get ( name,    ActionParams::TABLE_NAME );
  params.get ( weights, ActionParams::TABLE_WEIGHTS );

  // Stress value are computed in the nodes.

  if ( name == "stress" &&
       table->getRowItems() == nodes_.getData() )
  {
    getStress_ ( *table, weights, disp );

    return true;
  }
  else if ( name == "xoutTable" )
  {
    params.get ( contents, "contentString" );

    getXOutTable_ ( table, weights, contents, disp );

    return true;
  }
  return false;
}


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newShell6DofsBFModel
//-----------------------------------------------------------------------


static Ref<Model>     newShell6DofsBFModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<Shell6DofsBFModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareShell6DofsBFModel
//-----------------------------------------------------------------------


void declareShell6DofsBFModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "Shell6DofsBF", & newShell6DofsBFModel );
}