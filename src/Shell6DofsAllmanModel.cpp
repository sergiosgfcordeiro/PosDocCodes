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
 *  Allman 1988 - Evaluation of the constant strain triangle with drilling rotations.      
 *  
 *  Author:  S.G.F. Cordeiro, sferreiracorde@tudelft.nl
 *
 *  06 October 2025:
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

#include "Shell6DofsAllmanModel.h"
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

const char* Shell6DofsAllmanModel::DOF_NAMES[6]     = {"u","v","w","wx","wy","wz"};
const char* Shell6DofsAllmanModel::SHAPE_PROP       = "shape";
const char* Shell6DofsAllmanModel::MATERIAL_PROP    = "material";
const char* Shell6DofsAllmanModel::DRILLING_PROP    = "drilling_stif_factor";
const char* Shell6DofsAllmanModel::THICK_PROP       = "thickness";
const char* Shell6DofsAllmanModel::ORIENTATION_PROP[3] = {"v1","v2","v3"};    
const char* Shell6DofsAllmanModel::LARGE_DISP_PROP  = "largeDisp";
      idx_t Shell6DofsAllmanModel::nodesWritten_    = 0;
Ref<PrintWriter> Shell6DofsAllmanModel::nodeOut_    = NIL;

//-----------------------------------------------------------------------
//   constructor
//-----------------------------------------------------------------------

Shell6DofsAllmanModel::Shell6DofsAllmanModel

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

  thickness_ = 1.;

  myProps.find( thickness_, THICK_PROP );
  myConf.set  ( THICK_PROP, thickness_ );

  // drilling_stif_factor_ = 1.0e-3;
  drilling_stif_factor_ = 0.0;

  myProps.find( drilling_stif_factor_, DRILLING_PROP );
  myConf.set  ( DRILLING_PROP, drilling_stif_factor_ );

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

Shell6DofsAllmanModel::~Shell6DofsAllmanModel()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void Shell6DofsAllmanModel::configure

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


void Shell6DofsAllmanModel::getConfig 

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


bool Shell6DofsAllmanModel::takeAction

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


void Shell6DofsAllmanModel::getMatrix_

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

  Matrix      ipxcoords   ( rank_, ipCount_ );

  Matrix      B_mem       ( strCount_, 3*nodeCount_ );
  Matrix      Bw_mem      ( 1, 3*nodeCount_ );
  Matrix      N_theta     ( 1, 3*nodeCount_ );
  Matrix      N_residue   ( 1, 3*nodeCount_ );
  Matrix      K_mem       ( 3*nodeCount_, 3*nodeCount_ );
  Matrix      K_memstab   ( 3*nodeCount_, 3*nodeCount_ );
  Vector      Disp_mem    ( 3*nodeCount_ );

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

  double Area;

  Properties matProps;
  material_->getConfig(matProps);
  double Ex, Ey, Gxy, Nuxy;
  matProps.find( Ex, "young1" );
  matProps.find( Ey, "young2" );
  matProps.find( Gxy, "shear12" );
  matProps.find( Nuxy, "poisson12" );
  double nuyx = Nuxy*Ey/Ex;

  double mat_factor1 = Gxy * thickness_;
  double mat_factor2 = sqrt(Ex*Ey)*thickness_/(1.0 - Nuxy*nuyx);
  double mat_factor = std::max(mat_factor1, mat_factor2);



  double tau = drilling_stif_factor_ * mat_factor;  // stabilization parameter

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

    // Get the integration points in the local 2D barycentric (area) coordinate system of the T3 element.

    shape_->getGlobalIntegrationPoints ( ipxcoords, xcoords );

    // Get the Area of the T3 element.

    getArea_(Area, xcoords);

    // Compute the ipWeights for numerical integration.

    shape_->getIntegrationWeights ( ipWeights, xcoords );

    // Get the displacements at the element nodes in the Global coordinate system.

    elemDisp = disp[idofs];

    // Transform the displacements at the element nodes to the Local coordinate system.

    elemxDisp = mc1.matmul ( transMatDof, elemDisp );
    
    // Get the membrane displacements at the element nodes.

    getMembraneDisp_(Disp_mem, elemxDisp);

    // Assemble the membrane stiffness matrix.

    K_mem = 0.0;
    K_memstab = 0.0;
    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    { 
      // Get this integration point coordinates.

      double x = ipxcoords(0,ip);
      double y = ipxcoords(1,ip);

      // Compute the original B_matrix for the Allman membrane element at this integration point.

      getBmemAllman_ (B_mem, xcoords, Area, x, y);

      // Compute the Bw_mem matrices for the numerical stabilization of the Allman membrane element.

      getBwmemAllman_ (Bw_mem, xcoords, Area, x, y);

      // Compute the N_theta matrix for the numerical stabilization of the Allman membrane element.

      getNLinear_ (N_theta, xcoords, Area, x, y);

      // Compute the N_residue functions for the stabilization energy residue.

      N_residue = N_theta - Bw_mem;

      // Compute the membrane strain vector of this integration point.

      matmul ( strain_mem, B_mem, Disp_mem );

      // Get the tangent stiffness matrix and the stress vector.

      material_->update ( stress_mem, stiff, strain_mem, ipoint++ );

      // Assemble the membrane constitutive matrix.

      D_mem = stiff*thickness_;

      // Compute the membrane stiffness matrix.

      K_mem += ipWeights[ip] * mc3.matmul ( B_mem.transpose(), D_mem, B_mem );

      // Compute the drilling stabilization stiffness matrix.

      K_memstab += ipWeights[ip] * (2.0*tau) * matmul ( N_residue.transpose(), N_residue );
    }

    // Compute the stabilized membrane stiffness matrix.

    K_mem = K_mem + K_memstab;

    // // Print the contents of the K_mem matrix for debugging purposes.
    // for ( idx_t i = 0; i < K_mem.size(0); i++ )
    // {
    //   for ( idx_t j = 0; j < K_mem.size(1); j++ )
    //   {
    //     jem::System::out() << K_mem(i,j) << " ";
    //   }
    //   jem::System::out() << "\n";
    // }

    // Get the tangent stiffness matrix and the stress vector at the first integration point.

    int ip = 0;
    double  x = ipxcoords(0,ip);
    double  y = ipxcoords(1,ip);
    getBmemAllman_ (B_mem, xcoords, Area, x, y);
    matmul ( strain_mem, B_mem, Disp_mem );
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

    // Transforme the shell stiffness matrix back to the global system.

    K_shell = mc3.matmul ( transMatDof.transpose(), K_shell, transMatDof );

    // Compute and add the element force vector to the global force vector.

    elemForce = matmul ( K_shell, elemDisp );
    force[idofs] += elemForce;
  }
}


//-----------------------------------------------------------------------
//   writeLocalOutput_       NOT DONE FOR THE SHELL MODEL YET
//-----------------------------------------------------------------------


void Shell6DofsAllmanModel::writeLocalOutput_

  ( const Vector&       disp,
    const Properties&   globdat )

{
  idx_t       it;

  Matrix      stiff       ( strCount_, strCount_ );
  Matrix      coords      ( nodeRank_, nodeCount_ );
  Matrix      transMat    ( nodeRank_, nodeRank_ );
  Matrix      transMatDof ( dofCount_, dofCount_ );
  Matrix      xcoords     ( rank_, nodeCount_ );
  Matrix      ipxcoords   ( rank_, ipCount_ );
  Vector      elemDisp    ( dofCount_ );
  Vector      elemxDisp   ( dofCount_ );
  Vector      Disp_mem    ( 3*nodeCount_ );
  Vector      strain      ( strCount_ );
  Vector      stress      ( strCount_ );

  Matrix      B_mem       ( strCount_, 3*nodeCount_ );

  IdxVector   inodes      ( nodeCount_ );
  IdxVector   idofs       ( dofCount_  );

  Cubix       grads       ( rank_, nodeCount_, ipCount_ );
  Vector      ipWeights   ( ipCount_   );

  double Area;

  MChain1     mc1;

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

    shape_->getGlobalIntegrationPoints ( ipxcoords, xcoords );

    getArea_(Area, xcoords);

    elemDisp = disp[idofs];

    elemxDisp = mc1.matmul ( transMatDof, elemDisp );

    getMembraneDisp_(Disp_mem, elemxDisp);

    double x = ipxcoords(0,ip);
    double y = ipxcoords(1,ip);

    getBmemAllman_ (B_mem, xcoords, Area, x, y);

    matmul ( strain, B_mem, Disp_mem );

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


void Shell6DofsAllmanModel::getXOutTable_

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
  Matrix       ipxcoords   ( rank_, ipCount_ );
  Matrix       B_mem       ( strCount_, 3*nodeCount_ );
  Matrix       stiff       ( strCount_, strCount_  );

  Vector       elemDisp    ( dofCount_ );
  Vector       elemxDisp   ( dofCount_ );
  Vector       Disp_mem    ( 3*nodeCount_ );

  IdxVector    inodes      ( nodeCount_ );
  IdxVector    idofs       ( dofCount_  );

  double Area;

  MChain1     mc1;

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

    shape_->getGlobalIntegrationPoints ( ipxcoords, xcoords );

    getArea_(Area, xcoords);

    elemDisp = disp[idofs];

    elemxDisp = mc1.matmul ( transMatDof, elemDisp );

    getMembraneDisp_(Disp_mem, elemxDisp);

    // Iterate over the integration points.
    // Gather all data, no matter which is asked, to keep code neat
    // The option to specify output is primarily for disk size, not CPU time

    for ( idx_t ip = 0; ip < ipCount; ip++, ++ipoint )
    {
      double x = ipxcoords(0,ip);
      double y = ipxcoords(1,ip);
  
      getBmemAllman_ (B_mem, xcoords, Area, x, y);
  
      matmul ( strain, B_mem, Disp_mem );

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


void Shell6DofsAllmanModel::getTransMatrix_
   
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

  // Get the transformation matrix.

  transMat(0,0) = e1_l[0]*e1[0] + e1_l[1]*e1[1] + e1_l[2]*e1[2];
  transMat(0,1) = e1_l[0]*e2[0] + e1_l[1]*e2[1] + e1_l[2]*e2[2];
  transMat(0,2) = e1_l[0]*e3[0] + e1_l[1]*e3[1] + e1_l[2]*e3[2];
  transMat(1,0) = e2_l[0]*e1[0] + e2_l[1]*e1[1] + e2_l[2]*e1[2];
  transMat(1,1) = e2_l[0]*e2[0] + e2_l[1]*e2[1] + e2_l[2]*e2[2];
  transMat(1,2) = e2_l[0]*e3[0] + e2_l[1]*e3[1] + e2_l[2]*e3[2];
  transMat(2,0) = e3_l[0]*e1[0] + e3_l[1]*e1[1] + e3_l[2]*e1[2];
  transMat(2,1) = e3_l[0]*e2[0] + e3_l[1]*e2[1] + e3_l[2]*e2[2];
  transMat(2,2) = e3_l[0]*e3[0] + e3_l[1]*e3[1] + e3_l[2]*e3[2];

  // Get the transformation matrix of the Dofs of a single node.

  transMatBlock = 0.0;
  transMatBlock(slice(0,nodeRank_),slice(0,nodeRank_)) = transMat;
  transMatBlock(slice(nodeRank_,2*nodeRank_),slice(nodeRank_,2*nodeRank_)) = transMat;

  // Get the transformation matrix of the nodal element Dofs.

  transMatDof = 0.0;
  transMatDof(slice(0,6),slice(0,6)) = transMatBlock;
  transMatDof(slice(6,12),slice(6,12)) = transMatBlock;
  transMatDof(slice(12,18),slice(12,18)) = transMatBlock;
}


//-----------------------------------------------------------------------
//   get2DLocalcoordinates_
//-----------------------------------------------------------------------


void Shell6DofsAllmanModel::get2DLocalcoordinates_

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


void Shell6DofsAllmanModel::getArea_

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
//   getTmat_
//-----------------------------------------------------------------------


void Shell6DofsAllmanModel::getTmat_

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


void Shell6DofsAllmanModel::getMembraneDisp_

  ( Vector&         Disp_mem,
    const Vector&   elemxDisp )   const
    
{
  for ( idx_t i = 0; i < nodeCount_; ++i )
  {
    Disp_mem[3*i]   = elemxDisp[6*i];
    Disp_mem[3*i+1] = elemxDisp[6*i+1];
    Disp_mem[3*i+2] = elemxDisp[6*i+5];
  }
}


//-----------------------------------------------------------------------
//   getBmemAllman_
//-----------------------------------------------------------------------


void Shell6DofsAllmanModel::getBmemAllman_

  ( Matrix&           B_mem,
    const Matrix&     xcoords,
    const double&     Area,
    const double&     x,
    const double&     y ) const
{
  
  // Get some geometrical parameters.
  
  double x_1, x_2, x_3;
  double y_1, y_2, y_3;

  x_1 = xcoords(0,0);
  x_2 = xcoords(0,1);
  x_3 = xcoords(0,2);
  y_1 = xcoords(1,0);
  y_2 = xcoords(1,1);
  y_3 = xcoords(1,2);

  double x_12, x_21, x_23, x_32, x_13, x_31;

  x_12 = (xcoords(0, 0) - xcoords(0, 1));
  x_21 = (xcoords(0, 1) - xcoords(0, 0));
  x_23 = (xcoords(0, 1) - xcoords(0, 2));
  x_32 = (xcoords(0, 2) - xcoords(0, 1));
  x_13 = (xcoords(0, 0) - xcoords(0, 2));
  x_31 = (xcoords(0, 2) - xcoords(0, 0));

  double y_12, y_21, y_23, y_32, y_13, y_31;

  y_12 = (xcoords(1, 0) - xcoords(1, 1));
  y_21 = (xcoords(1, 1) - xcoords(1, 0));
  y_23 = (xcoords(1, 1) - xcoords(1, 2));
  y_32 = (xcoords(1, 2) - xcoords(1, 1));
  y_13 = (xcoords(1, 0) - xcoords(1, 2));
  y_31 = (xcoords(1, 2) - xcoords(1, 0));

  double l_12, l_23, l_31;

  l_12 = sqrt(pow(x_12,2) + pow(y_12,2));
  l_23 = sqrt(pow(x_23,2) + pow(y_23,2));
  l_31 = sqrt(pow(x_31,2) + pow(y_31,2));

  double a_12, b_12, a_23, b_23, a_31, b_31;

  a_12 = 2.0*Area/l_12;
  b_12 = (x_12*x_13 + y_12*y_13)/l_12;
  a_23 = 2.0*Area/l_23;
  b_23 = (x_23*x_21 + y_23*y_21)/l_23;
  a_31 = 2.0*Area/l_31;
  b_31 = (x_31*x_32 + y_31*y_32)/l_31;

  double xp_12, yp_12, xp_23, yp_23, xp_31, yp_31;

  xp_12 = x_21*b_12/l_12 + x_1;
  yp_12 = y_21*b_12/l_12 + y_1;
  xp_23 = x_32*b_23/l_23 + x_2;
  yp_23 = y_32*b_23/l_23 + y_2;
  xp_31 = x_13*b_31/l_31 + x_3;
  yp_31 = y_13*b_31/l_31 + y_3;

  double g1, g2, g3, g4, g5, g6;

  g1 = l_12*(xp_12 - x_3)/a_12;
  g2 = l_12*(yp_12 - y_3)/a_12;
  g3 = l_23*(xp_23 - x_1)/a_23;
  g4 = l_23*(yp_23 - y_1)/a_23;
  g5 = l_31*(xp_31 - x_2)/a_31;
  g6 = l_31*(yp_31 - y_2)/a_31;

  // Get the barycentric coordinates: zeta_1, zeta_2, zeta_3, zeta_12, zeta_23, zeta_31.

  double zeta_1, zeta_2, zeta_3;

  zeta_1 = (1.0/(2.0*Area))*((x_2 - x)*(y_3 - y) - (x_3 - x)*(y_2 - y));
  zeta_2 = (1.0/(2.0*Area))*((x_3 - x)*(y_1 - y) - (x_1 - x)*(y_3 - y));
  zeta_3 = (1.0/(2.0*Area))*((x_1 - x)*(y_2 - y) - (x_2 - x)*(y_1 - y));

  double zeta_12, zeta_23, zeta_31;

  zeta_12 = zeta_1*zeta_2;
  zeta_23 = zeta_2*zeta_3;
  zeta_31 = zeta_3*zeta_1;

  // Get the derivatives: dzeta1221_dzetai, dzeta2332_dzetai, dzeta3113_dzetai.

  double dzeta1221_dzeta1, dzeta1221_dzeta2, dzeta1221_dzeta3;

  dzeta1221_dzeta1 = zeta_2*(zeta_2-zeta_1) - zeta_12;
  dzeta1221_dzeta2 = zeta_1*(zeta_2-zeta_1) + zeta_12;
  dzeta1221_dzeta3 = 0.0;

  double dzeta2332_dzeta1, dzeta2332_dzeta2, dzeta2332_dzeta3;

  dzeta2332_dzeta1 = 0.0;
  dzeta2332_dzeta2 = zeta_3*(zeta_3-zeta_2) - zeta_23;
  dzeta2332_dzeta3 = zeta_2*(zeta_3-zeta_2) + zeta_23;

  double dzeta3113_dzeta1, dzeta3113_dzeta2, dzeta3113_dzeta3;

  dzeta3113_dzeta1 = zeta_3*(zeta_1-zeta_3) + zeta_31;
  dzeta3113_dzeta2 = 0.0;
  dzeta3113_dzeta3 = zeta_1*(zeta_1-zeta_3) - zeta_31;

  // Get the derivatives: xr0_dzetai, yro_dzetai.

  double dxr0_dzeta1, dxr0_dzeta2, dxr0_dzeta3;

  dxr0_dzeta1 = (g1*dzeta1221_dzeta1 + g3*dzeta2332_dzeta1 + g5*dzeta3113_dzeta1)/(4.0*Area);
  dxr0_dzeta2 = (g1*dzeta1221_dzeta2 + g3*dzeta2332_dzeta2 + g5*dzeta3113_dzeta2)/(4.0*Area);
  dxr0_dzeta3 = (g1*dzeta1221_dzeta3 + g3*dzeta2332_dzeta3 + g5*dzeta3113_dzeta3)/(4.0*Area);

  double dyr0_dzeta1, dyr0_dzeta2, dyr0_dzeta3;

  dyr0_dzeta1 = (g2*dzeta1221_dzeta1 + g4*dzeta2332_dzeta1 + g6*dzeta3113_dzeta1)/(4.0*Area);
  dyr0_dzeta2 = (g2*dzeta1221_dzeta2 + g4*dzeta2332_dzeta2 + g6*dzeta3113_dzeta2)/(4.0*Area);
  dyr0_dzeta3 = (g2*dzeta1221_dzeta3 + g4*dzeta2332_dzeta3 + g6*dzeta3113_dzeta3)/(4.0*Area);

  // Get the shape functions derivatives with respect to: zeta1, zeta2, zeta3.

  double dNxx1_dzeta1, dNxx1_dzeta2, dNxx1_dzeta3;

  dNxx1_dzeta1 = 1.0 - x_23*dxr0_dzeta1;
  dNxx1_dzeta2 = - x_23*dxr0_dzeta2;
  dNxx1_dzeta3 = - x_23*dxr0_dzeta3;

  double dNxy1_dzeta1, dNxy1_dzeta2, dNxy1_dzeta3;

  dNxy1_dzeta1 = - y_23*dxr0_dzeta1;
  dNxy1_dzeta2 = - y_23*dxr0_dzeta2;
  dNxy1_dzeta3 = - y_23*dxr0_dzeta3;

  double dNxx2_dzeta1, dNxx2_dzeta2, dNxx2_dzeta3;

  dNxx2_dzeta1 = - x_31*dxr0_dzeta1;
  dNxx2_dzeta2 = 1.0 - x_31*dxr0_dzeta2;
  dNxx2_dzeta3 = - x_31*dxr0_dzeta3;

  double dNxy2_dzeta1, dNxy2_dzeta2, dNxy2_dzeta3;

  dNxy2_dzeta1 = - y_31*dxr0_dzeta1;
  dNxy2_dzeta2 = - y_31*dxr0_dzeta2;
  dNxy2_dzeta3 = - y_31*dxr0_dzeta3;

  double dNxx3_dzeta1, dNxx3_dzeta2, dNxx3_dzeta3;

  dNxx3_dzeta1 = - x_12*dxr0_dzeta1;
  dNxx3_dzeta2 = - x_12*dxr0_dzeta2;
  dNxx3_dzeta3 = 1.0 - x_12*dxr0_dzeta3;

  double dNxy3_dzeta1, dNxy3_dzeta2, dNxy3_dzeta3;

  dNxy3_dzeta1 = - y_12*dxr0_dzeta1;
  dNxy3_dzeta2 = - y_12*dxr0_dzeta2;
  dNxy3_dzeta3 = - y_12*dxr0_dzeta3;

  double dNxtheta1_dzeta1, dNxtheta1_dzeta2, dNxtheta1_dzeta3;

  dNxtheta1_dzeta1 = (1.0/2.0)*(-g1*zeta_2 + g5*zeta_3 + g1*dzeta1221_dzeta1 + g5*dzeta3113_dzeta1);
  dNxtheta1_dzeta2 = (1.0/2.0)*(-g1*zeta_1 + g1*dzeta1221_dzeta2);
  dNxtheta1_dzeta3 = (1.0/2.0)*( g5*zeta_1 + g5*dzeta3113_dzeta3);

  double dNxtheta2_dzeta1, dNxtheta2_dzeta2, dNxtheta2_dzeta3;

  dNxtheta2_dzeta1 = (1.0/2.0)*(g1*zeta_2 + g1*dzeta1221_dzeta1);
  dNxtheta2_dzeta2 = (1.0/2.0)*(-g3*zeta_3 + g1*zeta_1 + g3*dzeta2332_dzeta2 + g1*dzeta1221_dzeta2);
  dNxtheta2_dzeta3 = (1.0/2.0)*(-g3*zeta_2 + g3*dzeta2332_dzeta3);

  double dNxtheta3_dzeta1, dNxtheta3_dzeta2, dNxtheta3_dzeta3;

  dNxtheta3_dzeta1 = (1.0/2.0)*(-g5*zeta_3 + g5*dzeta3113_dzeta1);
  dNxtheta3_dzeta2 = (1.0/2.0)*(g3*zeta_3 + g3*dzeta2332_dzeta2);
  dNxtheta3_dzeta3 = (1.0/2.0)*(-g5*zeta_1 + g3*zeta_2 + g5*dzeta3113_dzeta3 + g3*dzeta2332_dzeta3);

  double dNyx1_dzeta1, dNyx1_dzeta2, dNyx1_dzeta3;

  dNyx1_dzeta1 = - x_23*dyr0_dzeta1;
  dNyx1_dzeta2 = - x_23*dyr0_dzeta2;
  dNyx1_dzeta3 = - x_23*dyr0_dzeta3;

  double dNyy1_dzeta1, dNyy1_dzeta2, dNyy1_dzeta3;

  dNyy1_dzeta1 = 1.0 - y_23*dyr0_dzeta1;
  dNyy1_dzeta2 = - y_23*dyr0_dzeta2;
  dNyy1_dzeta3 = - y_23*dyr0_dzeta3;

  double dNyx2_dzeta1, dNyx2_dzeta2, dNyx2_dzeta3;

  dNyx2_dzeta1 = - x_31*dyr0_dzeta1;
  dNyx2_dzeta2 = - x_31*dyr0_dzeta2;
  dNyx2_dzeta3 = - x_31*dyr0_dzeta3;

  double dNyy2_dzeta1, dNyy2_dzeta2, dNyy2_dzeta3;

  dNyy2_dzeta1 = - y_31*dyr0_dzeta1;
  dNyy2_dzeta2 = 1.0 - y_31*dyr0_dzeta2;
  dNyy2_dzeta3 = - y_31*dyr0_dzeta3;

  double dNyx3_dzeta1, dNyx3_dzeta2, dNyx3_dzeta3;

  dNyx3_dzeta1 = - x_12*dyr0_dzeta1;
  dNyx3_dzeta2 = - x_12*dyr0_dzeta2;
  dNyx3_dzeta3 = - x_12*dyr0_dzeta3;

  double dNyy3_dzeta1, dNyy3_dzeta2, dNyy3_dzeta3;

  dNyy3_dzeta1 = - y_12*dyr0_dzeta1;
  dNyy3_dzeta2 = - y_12*dyr0_dzeta2;
  dNyy3_dzeta3 = 1.0 - y_12*dyr0_dzeta3;

  double dNytheta1_dzeta1, dNytheta1_dzeta2, dNytheta1_dzeta3;

  dNytheta1_dzeta1 = (1.0/2.0)*(-g2*zeta_2 + g6*zeta_3 + g2*dzeta1221_dzeta1 + g6*dzeta3113_dzeta1);
  dNytheta1_dzeta2 = (1.0/2.0)*(-g2*zeta_1 + g2*dzeta1221_dzeta2);
  dNytheta1_dzeta3 = (1.0/2.0)*( g6*zeta_1 + g6*dzeta3113_dzeta3);

  double dNytheta2_dzeta1, dNytheta2_dzeta2, dNytheta2_dzeta3;

  dNytheta2_dzeta1 = (1.0/2.0)*(g2*zeta_2 + g2*dzeta1221_dzeta1);
  dNytheta2_dzeta2 = (1.0/2.0)*(-g4*zeta_3 + g2*zeta_1 + g4*dzeta2332_dzeta2 + g2*dzeta1221_dzeta2);
  dNytheta2_dzeta3 = (1.0/2.0)*(-g4*zeta_2 + g4*dzeta2332_dzeta3);

  double dNytheta3_dzeta1, dNytheta3_dzeta2, dNytheta3_dzeta3;

  dNytheta3_dzeta1 = (1.0/2.0)*(-g6*zeta_3 + g6*dzeta3113_dzeta1);
  dNytheta3_dzeta2 = (1.0/2.0)*(g4*zeta_3 + g4*dzeta2332_dzeta2);
  dNytheta3_dzeta3 = (1.0/2.0)*(-g6*zeta_1 + g4*zeta_2 + g6*dzeta3113_dzeta3 + g4*dzeta2332_dzeta3);


  // Get the jacobians: dzetai_dx, dzetai_dy.

  double dzeta1_dx, dzeta2_dx, dzeta3_dx;

  dzeta1_dx = (y_2 - y_3)/(2.0*Area);
  dzeta2_dx = (y_3 - y_1)/(2.0*Area);
  dzeta3_dx = (y_1 - y_2)/(2.0*Area);

  double dzeta1_dy, dzeta2_dy, dzeta3_dy;

  dzeta1_dy = (x_3 - x_2)/(2.0*Area);
  dzeta2_dy = (x_1 - x_3)/(2.0*Area);
  dzeta3_dy = (x_2 - x_1)/(2.0*Area);

  // Get the shape functions derivatives with respect to: x, y.

  // Displacement u shape functions derivatives

  double dNxx1_dx, dNxy1_dx, dNxx2_dx, dNxy2_dx, dNxx3_dx, dNxy3_dx;

  dNxx1_dx = dNxx1_dzeta1*dzeta1_dx + dNxx1_dzeta2*dzeta2_dx + dNxx1_dzeta3*dzeta3_dx;
  dNxy1_dx = dNxy1_dzeta1*dzeta1_dx + dNxy1_dzeta2*dzeta2_dx + dNxy1_dzeta3*dzeta3_dx;
  dNxx2_dx = dNxx2_dzeta1*dzeta1_dx + dNxx2_dzeta2*dzeta2_dx + dNxx2_dzeta3*dzeta3_dx;
  dNxy2_dx = dNxy2_dzeta1*dzeta1_dx + dNxy2_dzeta2*dzeta2_dx + dNxy2_dzeta3*dzeta3_dx;
  dNxx3_dx = dNxx3_dzeta1*dzeta1_dx + dNxx3_dzeta2*dzeta2_dx + dNxx3_dzeta3*dzeta3_dx;
  dNxy3_dx = dNxy3_dzeta1*dzeta1_dx + dNxy3_dzeta2*dzeta2_dx + dNxy3_dzeta3*dzeta3_dx;

  double dNxx1_dy, dNxy1_dy, dNxx2_dy, dNxy2_dy, dNxx3_dy, dNxy3_dy;

  dNxx1_dy = dNxx1_dzeta1*dzeta1_dy + dNxx1_dzeta2*dzeta2_dy + dNxx1_dzeta3*dzeta3_dy;
  dNxy1_dy = dNxy1_dzeta1*dzeta1_dy + dNxy1_dzeta2*dzeta2_dy + dNxy1_dzeta3*dzeta3_dy;
  dNxx2_dy = dNxx2_dzeta1*dzeta1_dy + dNxx2_dzeta2*dzeta2_dy + dNxx2_dzeta3*dzeta3_dy;
  dNxy2_dy = dNxy2_dzeta1*dzeta1_dy + dNxy2_dzeta2*dzeta2_dy + dNxy2_dzeta3*dzeta3_dy;
  dNxx3_dy = dNxx3_dzeta1*dzeta1_dy + dNxx3_dzeta2*dzeta2_dy + dNxx3_dzeta3*dzeta3_dy;
  dNxy3_dy = dNxy3_dzeta1*dzeta1_dy + dNxy3_dzeta2*dzeta2_dy + dNxy3_dzeta3*dzeta3_dy;

  double dNxtheta1_dx, dNxtheta2_dx, dNxtheta3_dx;

  dNxtheta1_dx = dNxtheta1_dzeta1*dzeta1_dx + dNxtheta1_dzeta2*dzeta2_dx + dNxtheta1_dzeta3*dzeta3_dx;
  dNxtheta2_dx = dNxtheta2_dzeta1*dzeta1_dx + dNxtheta2_dzeta2*dzeta2_dx + dNxtheta2_dzeta3*dzeta3_dx;
  dNxtheta3_dx = dNxtheta3_dzeta1*dzeta1_dx + dNxtheta3_dzeta2*dzeta2_dx + dNxtheta3_dzeta3*dzeta3_dx;

  double dNxtheta1_dy, dNxtheta2_dy, dNxtheta3_dy;

  dNxtheta1_dy = dNxtheta1_dzeta1*dzeta1_dy + dNxtheta1_dzeta2*dzeta2_dy + dNxtheta1_dzeta3*dzeta3_dy;
  dNxtheta2_dy = dNxtheta2_dzeta1*dzeta1_dy + dNxtheta2_dzeta2*dzeta2_dy + dNxtheta2_dzeta3*dzeta3_dy;
  dNxtheta3_dy = dNxtheta3_dzeta1*dzeta1_dy + dNxtheta3_dzeta2*dzeta2_dy + dNxtheta3_dzeta3*dzeta3_dy;

  // Displacement v shape functions derivatives

  double dNyx1_dx, dNyy1_dx, dNyx2_dx, dNyy2_dx, dNyx3_dx, dNyy3_dx;

  dNyx1_dx = dNyx1_dzeta1*dzeta1_dx + dNyx1_dzeta2*dzeta2_dx + dNyx1_dzeta3*dzeta3_dx;
  dNyy1_dx = dNyy1_dzeta1*dzeta1_dx + dNyy1_dzeta2*dzeta2_dx + dNyy1_dzeta3*dzeta3_dx;
  dNyx2_dx = dNyx2_dzeta1*dzeta1_dx + dNyx2_dzeta2*dzeta2_dx + dNyx2_dzeta3*dzeta3_dx;
  dNyy2_dx = dNyy2_dzeta1*dzeta1_dx + dNyy2_dzeta2*dzeta2_dx + dNyy2_dzeta3*dzeta3_dx;
  dNyx3_dx = dNyx3_dzeta1*dzeta1_dx + dNyx3_dzeta2*dzeta2_dx + dNyx3_dzeta3*dzeta3_dx;
  dNyy3_dx = dNyy3_dzeta1*dzeta1_dx + dNyy3_dzeta2*dzeta2_dx + dNyy3_dzeta3*dzeta3_dx;

  double dNyx1_dy, dNyy1_dy, dNyx2_dy, dNyy2_dy, dNyx3_dy, dNyy3_dy;

  dNyx1_dy = dNyx1_dzeta1*dzeta1_dy + dNyx1_dzeta2*dzeta2_dy + dNyx1_dzeta3*dzeta3_dy;
  dNyy1_dy = dNyy1_dzeta1*dzeta1_dy + dNyy1_dzeta2*dzeta2_dy + dNyy1_dzeta3*dzeta3_dy;
  dNyx2_dy = dNyx2_dzeta1*dzeta1_dy + dNyx2_dzeta2*dzeta2_dy + dNyx2_dzeta3*dzeta3_dy;
  dNyy2_dy = dNyy2_dzeta1*dzeta1_dy + dNyy2_dzeta2*dzeta2_dy + dNyy2_dzeta3*dzeta3_dy;
  dNyx3_dy = dNyx3_dzeta1*dzeta1_dy + dNyx3_dzeta2*dzeta2_dy + dNyx3_dzeta3*dzeta3_dy;
  dNyy3_dy = dNyy3_dzeta1*dzeta1_dy + dNyy3_dzeta2*dzeta2_dy + dNyy3_dzeta3*dzeta3_dy;

  double dNytheta1_dx, dNytheta2_dx, dNytheta3_dx;

  dNytheta1_dx = dNytheta1_dzeta1*dzeta1_dx + dNytheta1_dzeta2*dzeta2_dx + dNytheta1_dzeta3*dzeta3_dx;
  dNytheta2_dx = dNytheta2_dzeta1*dzeta1_dx + dNytheta2_dzeta2*dzeta2_dx + dNytheta2_dzeta3*dzeta3_dx;
  dNytheta3_dx = dNytheta3_dzeta1*dzeta1_dx + dNytheta3_dzeta2*dzeta2_dx + dNytheta3_dzeta3*dzeta3_dx;

  double dNytheta1_dy, dNytheta2_dy, dNytheta3_dy;

  dNytheta1_dy = dNytheta1_dzeta1*dzeta1_dy + dNytheta1_dzeta2*dzeta2_dy + dNytheta1_dzeta3*dzeta3_dy;
  dNytheta2_dy = dNytheta2_dzeta1*dzeta1_dy + dNytheta2_dzeta2*dzeta2_dy + dNytheta2_dzeta3*dzeta3_dy;
  dNytheta3_dy = dNytheta3_dzeta1*dzeta1_dy + dNytheta3_dzeta2*dzeta2_dy + dNytheta3_dzeta3*dzeta3_dy;

  // Compute the B_mem matrix for the Allman membrane

  B_mem = 0.0;
  B_mem(0,0) = dNxx1_dx;
  B_mem(0,1) = dNxy1_dx; 
  B_mem(0,2) = dNxtheta1_dx; 
  B_mem(0,3) = dNxx2_dx;
  B_mem(0,4) = dNxy2_dx; 
  B_mem(0,5) = dNxtheta2_dx; 
  B_mem(0,6) = dNxx3_dx;
  B_mem(0,7) = dNxy3_dx; 
  B_mem(0,8) = dNxtheta3_dx; 

  B_mem(1,0) = dNyx1_dy;
  B_mem(1,1) = dNyy1_dy; 
  B_mem(1,2) = dNytheta1_dy; 
  B_mem(1,3) = dNyx2_dy;
  B_mem(1,4) = dNyy2_dy; 
  B_mem(1,5) = dNytheta2_dy; 
  B_mem(1,6) = dNyx3_dy;
  B_mem(1,7) = dNyy3_dy; 
  B_mem(1,8) = dNytheta3_dy; 

  B_mem(2,0) = dNxx1_dy + dNyx1_dx;
  B_mem(2,1) = dNxy1_dy + dNyy1_dx; 
  B_mem(2,2) = dNxtheta1_dy + dNytheta1_dx; 
  B_mem(2,3) = dNxx2_dy + dNyx2_dx;
  B_mem(2,4) = dNxy2_dy + dNyy2_dx; 
  B_mem(2,5) = dNxtheta2_dy + dNytheta2_dx; 
  B_mem(2,6) = dNxx3_dy + dNyx3_dx;
  B_mem(2,7) = dNxy3_dy + dNyy3_dx; 
  B_mem(2,8) = dNxtheta3_dy + dNytheta3_dx; 
}

//-----------------------------------------------------------------------
//   getBwmemAllman_
//-----------------------------------------------------------------------


void Shell6DofsAllmanModel::getBwmemAllman_

  ( Matrix&           Bw_mem,
    const Matrix&     xcoords,
    const double&     Area,
    const double&     x,
    const double&     y ) const
{
  
  // Get some geometrical parameters.
  
  double x_1, x_2, x_3;
  double y_1, y_2, y_3;

  x_1 = xcoords(0,0);
  x_2 = xcoords(0,1);
  x_3 = xcoords(0,2);
  y_1 = xcoords(1,0);
  y_2 = xcoords(1,1);
  y_3 = xcoords(1,2);

  double x_12, x_21, x_23, x_32, x_13, x_31;

  x_12 = (xcoords(0, 0) - xcoords(0, 1));
  x_21 = (xcoords(0, 1) - xcoords(0, 0));
  x_23 = (xcoords(0, 1) - xcoords(0, 2));
  x_32 = (xcoords(0, 2) - xcoords(0, 1));
  x_13 = (xcoords(0, 0) - xcoords(0, 2));
  x_31 = (xcoords(0, 2) - xcoords(0, 0));

  double y_12, y_21, y_23, y_32, y_13, y_31;

  y_12 = (xcoords(1, 0) - xcoords(1, 1));
  y_21 = (xcoords(1, 1) - xcoords(1, 0));
  y_23 = (xcoords(1, 1) - xcoords(1, 2));
  y_32 = (xcoords(1, 2) - xcoords(1, 1));
  y_13 = (xcoords(1, 0) - xcoords(1, 2));
  y_31 = (xcoords(1, 2) - xcoords(1, 0));

  double l_12, l_23, l_31;

  l_12 = sqrt(pow(x_12,2) + pow(y_12,2));
  l_23 = sqrt(pow(x_23,2) + pow(y_23,2));
  l_31 = sqrt(pow(x_31,2) + pow(y_31,2));

  double a_12, b_12, a_23, b_23, a_31, b_31;

  a_12 = 2.0*Area/l_12;
  b_12 = (x_12*x_13 + y_12*y_13)/l_12;
  a_23 = 2.0*Area/l_23;
  b_23 = (x_23*x_21 + y_23*y_21)/l_23;
  a_31 = 2.0*Area/l_31;
  b_31 = (x_31*x_32 + y_31*y_32)/l_31;

  double xp_12, yp_12, xp_23, yp_23, xp_31, yp_31;

  xp_12 = x_21*b_12/l_12 + x_1;
  yp_12 = y_21*b_12/l_12 + y_1;
  xp_23 = x_32*b_23/l_23 + x_2;
  yp_23 = y_32*b_23/l_23 + y_2;
  xp_31 = x_13*b_31/l_31 + x_3;
  yp_31 = y_13*b_31/l_31 + y_3;

  double g1, g2, g3, g4, g5, g6;

  g1 = l_12*(xp_12 - x_3)/a_12;
  g2 = l_12*(yp_12 - y_3)/a_12;
  g3 = l_23*(xp_23 - x_1)/a_23;
  g4 = l_23*(yp_23 - y_1)/a_23;
  g5 = l_31*(xp_31 - x_2)/a_31;
  g6 = l_31*(yp_31 - y_2)/a_31;

  // Get the barycentric coordinates: zeta_1, zeta_2, zeta_3, zeta_12, zeta_23, zeta_31.

  double zeta_1, zeta_2, zeta_3;

  zeta_1 = (1.0/(2.0*Area))*((x_2 - x)*(y_3 - y) - (x_3 - x)*(y_2 - y));
  zeta_2 = (1.0/(2.0*Area))*((x_3 - x)*(y_1 - y) - (x_1 - x)*(y_3 - y));
  zeta_3 = (1.0/(2.0*Area))*((x_1 - x)*(y_2 - y) - (x_2 - x)*(y_1 - y));

  double zeta_12, zeta_23, zeta_31;

  zeta_12 = zeta_1*zeta_2;
  zeta_23 = zeta_2*zeta_3;
  zeta_31 = zeta_3*zeta_1;

  // Get the derivatives: dzeta1221_dzetai, dzeta2332_dzetai, dzeta3113_dzetai.

  double dzeta1221_dzeta1, dzeta1221_dzeta2, dzeta1221_dzeta3;

  dzeta1221_dzeta1 = zeta_2*(zeta_2-zeta_1) - zeta_12;
  dzeta1221_dzeta2 = zeta_1*(zeta_2-zeta_1) + zeta_12;
  dzeta1221_dzeta3 = 0.0;

  double dzeta2332_dzeta1, dzeta2332_dzeta2, dzeta2332_dzeta3;

  dzeta2332_dzeta1 = 0.0;
  dzeta2332_dzeta2 = zeta_3*(zeta_3-zeta_2) - zeta_23;
  dzeta2332_dzeta3 = zeta_2*(zeta_3-zeta_2) + zeta_23;

  double dzeta3113_dzeta1, dzeta3113_dzeta2, dzeta3113_dzeta3;

  dzeta3113_dzeta1 = zeta_3*(zeta_1-zeta_3) + zeta_31;
  dzeta3113_dzeta2 = 0.0;
  dzeta3113_dzeta3 = zeta_1*(zeta_1-zeta_3) - zeta_31;

  // Get the derivatives: xr0_dzetai, yro_dzetai.

  double dxr0_dzeta1, dxr0_dzeta2, dxr0_dzeta3;

  dxr0_dzeta1 = (g1*dzeta1221_dzeta1 + g3*dzeta2332_dzeta1 + g5*dzeta3113_dzeta1)/(4.0*Area);
  dxr0_dzeta2 = (g1*dzeta1221_dzeta2 + g3*dzeta2332_dzeta2 + g5*dzeta3113_dzeta2)/(4.0*Area);
  dxr0_dzeta3 = (g1*dzeta1221_dzeta3 + g3*dzeta2332_dzeta3 + g5*dzeta3113_dzeta3)/(4.0*Area);

  double dyr0_dzeta1, dyr0_dzeta2, dyr0_dzeta3;

  dyr0_dzeta1 = (g2*dzeta1221_dzeta1 + g4*dzeta2332_dzeta1 + g6*dzeta3113_dzeta1)/(4.0*Area);
  dyr0_dzeta2 = (g2*dzeta1221_dzeta2 + g4*dzeta2332_dzeta2 + g6*dzeta3113_dzeta2)/(4.0*Area);
  dyr0_dzeta3 = (g2*dzeta1221_dzeta3 + g4*dzeta2332_dzeta3 + g6*dzeta3113_dzeta3)/(4.0*Area);

  // Get the shape functions derivatives with respect to: zeta1, zeta2, zeta3.

  double dNxx1_dzeta1, dNxx1_dzeta2, dNxx1_dzeta3;

  dNxx1_dzeta1 = 1.0 - x_23*dxr0_dzeta1;
  dNxx1_dzeta2 = - x_23*dxr0_dzeta2;
  dNxx1_dzeta3 = - x_23*dxr0_dzeta3;

  double dNxy1_dzeta1, dNxy1_dzeta2, dNxy1_dzeta3;

  dNxy1_dzeta1 = - y_23*dxr0_dzeta1;
  dNxy1_dzeta2 = - y_23*dxr0_dzeta2;
  dNxy1_dzeta3 = - y_23*dxr0_dzeta3;

  double dNxx2_dzeta1, dNxx2_dzeta2, dNxx2_dzeta3;

  dNxx2_dzeta1 = - x_31*dxr0_dzeta1;
  dNxx2_dzeta2 = 1.0 - x_31*dxr0_dzeta2;
  dNxx2_dzeta3 = - x_31*dxr0_dzeta3;

  double dNxy2_dzeta1, dNxy2_dzeta2, dNxy2_dzeta3;

  dNxy2_dzeta1 = - y_31*dxr0_dzeta1;
  dNxy2_dzeta2 = - y_31*dxr0_dzeta2;
  dNxy2_dzeta3 = - y_31*dxr0_dzeta3;

  double dNxx3_dzeta1, dNxx3_dzeta2, dNxx3_dzeta3;

  dNxx3_dzeta1 = - x_12*dxr0_dzeta1;
  dNxx3_dzeta2 = - x_12*dxr0_dzeta2;
  dNxx3_dzeta3 = 1.0 - x_12*dxr0_dzeta3;

  double dNxy3_dzeta1, dNxy3_dzeta2, dNxy3_dzeta3;

  dNxy3_dzeta1 = - y_12*dxr0_dzeta1;
  dNxy3_dzeta2 = - y_12*dxr0_dzeta2;
  dNxy3_dzeta3 = - y_12*dxr0_dzeta3;

  double dNxtheta1_dzeta1, dNxtheta1_dzeta2, dNxtheta1_dzeta3;

  dNxtheta1_dzeta1 = (1.0/2.0)*(-g1*zeta_2 + g5*zeta_3 + g1*dzeta1221_dzeta1 + g5*dzeta3113_dzeta1);
  dNxtheta1_dzeta2 = (1.0/2.0)*(-g1*zeta_1 + g1*dzeta1221_dzeta2);
  dNxtheta1_dzeta3 = (1.0/2.0)*( g5*zeta_1 + g5*dzeta3113_dzeta3);

  double dNxtheta2_dzeta1, dNxtheta2_dzeta2, dNxtheta2_dzeta3;

  dNxtheta2_dzeta1 = (1.0/2.0)*(g1*zeta_2 + g1*dzeta1221_dzeta1);
  dNxtheta2_dzeta2 = (1.0/2.0)*(-g3*zeta_3 + g1*zeta_1 + g3*dzeta2332_dzeta2 + g1*dzeta1221_dzeta2);
  dNxtheta2_dzeta3 = (1.0/2.0)*(-g3*zeta_2 + g3*dzeta2332_dzeta3);

  double dNxtheta3_dzeta1, dNxtheta3_dzeta2, dNxtheta3_dzeta3;

  dNxtheta3_dzeta1 = (1.0/2.0)*(-g5*zeta_3 + g5*dzeta3113_dzeta1);
  dNxtheta3_dzeta2 = (1.0/2.0)*(g3*zeta_3 + g3*dzeta2332_dzeta2);
  dNxtheta3_dzeta3 = (1.0/2.0)*(-g5*zeta_1 + g3*zeta_2 + g5*dzeta3113_dzeta3 + g3*dzeta2332_dzeta3);

  double dNyx1_dzeta1, dNyx1_dzeta2, dNyx1_dzeta3;

  dNyx1_dzeta1 = - x_23*dyr0_dzeta1;
  dNyx1_dzeta2 = - x_23*dyr0_dzeta2;
  dNyx1_dzeta3 = - x_23*dyr0_dzeta3;

  double dNyy1_dzeta1, dNyy1_dzeta2, dNyy1_dzeta3;

  dNyy1_dzeta1 = 1.0 - y_23*dyr0_dzeta1;
  dNyy1_dzeta2 = - y_23*dyr0_dzeta2;
  dNyy1_dzeta3 = - y_23*dyr0_dzeta3;

  double dNyx2_dzeta1, dNyx2_dzeta2, dNyx2_dzeta3;

  dNyx2_dzeta1 = - x_31*dyr0_dzeta1;
  dNyx2_dzeta2 = - x_31*dyr0_dzeta2;
  dNyx2_dzeta3 = - x_31*dyr0_dzeta3;

  double dNyy2_dzeta1, dNyy2_dzeta2, dNyy2_dzeta3;

  dNyy2_dzeta1 = - y_31*dyr0_dzeta1;
  dNyy2_dzeta2 = 1.0 - y_31*dyr0_dzeta2;
  dNyy2_dzeta3 = - y_31*dyr0_dzeta3;

  double dNyx3_dzeta1, dNyx3_dzeta2, dNyx3_dzeta3;

  dNyx3_dzeta1 = - x_12*dyr0_dzeta1;
  dNyx3_dzeta2 = - x_12*dyr0_dzeta2;
  dNyx3_dzeta3 = - x_12*dyr0_dzeta3;

  double dNyy3_dzeta1, dNyy3_dzeta2, dNyy3_dzeta3;

  dNyy3_dzeta1 = - y_12*dyr0_dzeta1;
  dNyy3_dzeta2 = - y_12*dyr0_dzeta2;
  dNyy3_dzeta3 = 1.0 - y_12*dyr0_dzeta3;

  double dNytheta1_dzeta1, dNytheta1_dzeta2, dNytheta1_dzeta3;

  dNytheta1_dzeta1 = (1.0/2.0)*(-g2*zeta_2 + g6*zeta_3 + g2*dzeta1221_dzeta1 + g6*dzeta3113_dzeta1);
  dNytheta1_dzeta2 = (1.0/2.0)*(-g2*zeta_1 + g2*dzeta1221_dzeta2);
  dNytheta1_dzeta3 = (1.0/2.0)*( g6*zeta_1 + g6*dzeta3113_dzeta3);

  double dNytheta2_dzeta1, dNytheta2_dzeta2, dNytheta2_dzeta3;

  dNytheta2_dzeta1 = (1.0/2.0)*(g2*zeta_2 + g2*dzeta1221_dzeta1);
  dNytheta2_dzeta2 = (1.0/2.0)*(-g4*zeta_3 + g2*zeta_1 + g4*dzeta2332_dzeta2 + g2*dzeta1221_dzeta2);
  dNytheta2_dzeta3 = (1.0/2.0)*(-g4*zeta_2 + g4*dzeta2332_dzeta3);

  double dNytheta3_dzeta1, dNytheta3_dzeta2, dNytheta3_dzeta3;

  dNytheta3_dzeta1 = (1.0/2.0)*(-g6*zeta_3 + g6*dzeta3113_dzeta1);
  dNytheta3_dzeta2 = (1.0/2.0)*(g4*zeta_3 + g4*dzeta2332_dzeta2);
  dNytheta3_dzeta3 = (1.0/2.0)*(-g6*zeta_1 + g4*zeta_2 + g6*dzeta3113_dzeta3 + g4*dzeta2332_dzeta3);


  // Get the jacobians: dzetai_dx, dzetai_dy.

  double dzeta1_dx, dzeta2_dx, dzeta3_dx;

  dzeta1_dx = (y_2 - y_3)/(2.0*Area);
  dzeta2_dx = (y_3 - y_1)/(2.0*Area);
  dzeta3_dx = (y_1 - y_2)/(2.0*Area);

  double dzeta1_dy, dzeta2_dy, dzeta3_dy;

  dzeta1_dy = (x_3 - x_2)/(2.0*Area);
  dzeta2_dy = (x_1 - x_3)/(2.0*Area);
  dzeta3_dy = (x_2 - x_1)/(2.0*Area);

  // Get the shape functions derivatives with respect to: x, y.

  // Displacement u shape functions derivatives

  double dNxx1_dy, dNxy1_dy, dNxx2_dy, dNxy2_dy, dNxx3_dy, dNxy3_dy;

  dNxx1_dy = dNxx1_dzeta1*dzeta1_dy + dNxx1_dzeta2*dzeta2_dy + dNxx1_dzeta3*dzeta3_dy;
  dNxy1_dy = dNxy1_dzeta1*dzeta1_dy + dNxy1_dzeta2*dzeta2_dy + dNxy1_dzeta3*dzeta3_dy;
  dNxx2_dy = dNxx2_dzeta1*dzeta1_dy + dNxx2_dzeta2*dzeta2_dy + dNxx2_dzeta3*dzeta3_dy;
  dNxy2_dy = dNxy2_dzeta1*dzeta1_dy + dNxy2_dzeta2*dzeta2_dy + dNxy2_dzeta3*dzeta3_dy;
  dNxx3_dy = dNxx3_dzeta1*dzeta1_dy + dNxx3_dzeta2*dzeta2_dy + dNxx3_dzeta3*dzeta3_dy;
  dNxy3_dy = dNxy3_dzeta1*dzeta1_dy + dNxy3_dzeta2*dzeta2_dy + dNxy3_dzeta3*dzeta3_dy;

  double dNxtheta1_dy, dNxtheta2_dy, dNxtheta3_dy;

  dNxtheta1_dy = dNxtheta1_dzeta1*dzeta1_dy + dNxtheta1_dzeta2*dzeta2_dy + dNxtheta1_dzeta3*dzeta3_dy;
  dNxtheta2_dy = dNxtheta2_dzeta1*dzeta1_dy + dNxtheta2_dzeta2*dzeta2_dy + dNxtheta2_dzeta3*dzeta3_dy;
  dNxtheta3_dy = dNxtheta3_dzeta1*dzeta1_dy + dNxtheta3_dzeta2*dzeta2_dy + dNxtheta3_dzeta3*dzeta3_dy;

  // Displacement v shape functions derivatives

  double dNyx1_dx, dNyy1_dx, dNyx2_dx, dNyy2_dx, dNyx3_dx, dNyy3_dx;

  dNyx1_dx = dNyx1_dzeta1*dzeta1_dx + dNyx1_dzeta2*dzeta2_dx + dNyx1_dzeta3*dzeta3_dx;
  dNyy1_dx = dNyy1_dzeta1*dzeta1_dx + dNyy1_dzeta2*dzeta2_dx + dNyy1_dzeta3*dzeta3_dx;
  dNyx2_dx = dNyx2_dzeta1*dzeta1_dx + dNyx2_dzeta2*dzeta2_dx + dNyx2_dzeta3*dzeta3_dx;
  dNyy2_dx = dNyy2_dzeta1*dzeta1_dx + dNyy2_dzeta2*dzeta2_dx + dNyy2_dzeta3*dzeta3_dx;
  dNyx3_dx = dNyx3_dzeta1*dzeta1_dx + dNyx3_dzeta2*dzeta2_dx + dNyx3_dzeta3*dzeta3_dx;
  dNyy3_dx = dNyy3_dzeta1*dzeta1_dx + dNyy3_dzeta2*dzeta2_dx + dNyy3_dzeta3*dzeta3_dx;

  double dNytheta1_dx, dNytheta2_dx, dNytheta3_dx;

  dNytheta1_dx = dNytheta1_dzeta1*dzeta1_dx + dNytheta1_dzeta2*dzeta2_dx + dNytheta1_dzeta3*dzeta3_dx;
  dNytheta2_dx = dNytheta2_dzeta1*dzeta1_dx + dNytheta2_dzeta2*dzeta2_dx + dNytheta2_dzeta3*dzeta3_dx;
  dNytheta3_dx = dNytheta3_dzeta1*dzeta1_dx + dNytheta3_dzeta2*dzeta2_dx + dNytheta3_dzeta3*dzeta3_dx;

  // Compute the Bw_mem matrix for the Allman membrane

  Bw_mem(0,0) = (1.0/2.0)*(dNyx1_dx - dNxx1_dy);
  Bw_mem(0,1) = (1.0/2.0)*(dNyy1_dx - dNxy1_dy); 
  Bw_mem(0,2) = (1.0/2.0)*(dNytheta1_dx - dNxtheta1_dy); 
  Bw_mem(0,3) = (1.0/2.0)*(dNyx2_dx - dNxx2_dy);
  Bw_mem(0,4) = (1.0/2.0)*(dNyy2_dx - dNxy2_dy); 
  Bw_mem(0,5) = (1.0/2.0)*(dNytheta2_dx - dNxtheta2_dy); 
  Bw_mem(0,6) = (1.0/2.0)*(dNyx3_dx - dNxx3_dy);
  Bw_mem(0,7) = (1.0/2.0)*(dNyy3_dx - dNxy3_dy); 
  Bw_mem(0,8) = (1.0/2.0)*(dNytheta3_dx - dNxtheta3_dy); 
}

//-----------------------------------------------------------------------
//   getNLinear_
//-----------------------------------------------------------------------


void Shell6DofsAllmanModel::getNLinear_

  ( Matrix&           N_theta,
    const Matrix&     xcoords,
    const double&     Area,
    const double&     x,
    const double&     y ) const
{
  
  // Get some geometrical parameters.
  
  double x_1, x_2, x_3;
  double y_1, y_2, y_3;

  x_1 = xcoords(0,0);
  x_2 = xcoords(0,1);
  x_3 = xcoords(0,2);
  y_1 = xcoords(1,0);
  y_2 = xcoords(1,1);
  y_3 = xcoords(1,2);

  // Get the barycentric coordinates: zeta_1, zeta_2, zeta_3, zeta_12, zeta_23, zeta_31.

  double zeta_1, zeta_2, zeta_3;

  zeta_1 = (1.0/(2.0*Area))*((x_2 - x)*(y_3 - y) - (x_3 - x)*(y_2 - y));
  zeta_2 = (1.0/(2.0*Area))*((x_3 - x)*(y_1 - y) - (x_1 - x)*(y_3 - y));
  zeta_3 = (1.0/(2.0*Area))*((x_1 - x)*(y_2 - y) - (x_2 - x)*(y_1 - y));

  // Compute the Bw_mem matrix for the Allman membrane

  N_theta(0,0) = 0.0;
  N_theta(0,1) = 0.0; 
  N_theta(0,2) = zeta_1; 
  N_theta(0,3) = 0.0;
  N_theta(0,4) = 0.0; 
  N_theta(0,5) = zeta_2;
  N_theta(0,6) = 0.0;
  N_theta(0,7) = 0.0; 
  N_theta(0,8) = zeta_3;
}


// //-----------------------------------------------------------------------
// //   getBmemAllmanReduced_
// //-----------------------------------------------------------------------


// void Shell6DofsAllmanModel::getBmemAllmanReduced_

//   ( Matrix&           B_mem,
//     const Matrix&     xcoords,
//     const double&     Area,
//     const double&     x,
//     const double&     y ) const
// {
  
//   // Get some geometrical parameters.
  
//   double x_1, x_2, x_3;
//   double y_1, y_2, y_3;

//   x_1 = xcoords(0,0);
//   x_2 = xcoords(0,1);
//   x_3 = xcoords(0,2);
//   y_1 = xcoords(1,0);
//   y_2 = xcoords(1,1);
//   y_3 = xcoords(1,2);

//   double x_12, x_21, x_23, x_32, x_13, x_31;

//   x_12 = (xcoords(0, 0) - xcoords(0, 1));
//   x_21 = (xcoords(0, 1) - xcoords(0, 0));
//   x_23 = (xcoords(0, 1) - xcoords(0, 2));
//   x_32 = (xcoords(0, 2) - xcoords(0, 1));
//   x_13 = (xcoords(0, 0) - xcoords(0, 2));
//   x_31 = (xcoords(0, 2) - xcoords(0, 0));

//   double y_12, y_21, y_23, y_32, y_13, y_31;

//   y_12 = (xcoords(1, 0) - xcoords(1, 1));
//   y_21 = (xcoords(1, 1) - xcoords(1, 0));
//   y_23 = (xcoords(1, 1) - xcoords(1, 2));
//   y_32 = (xcoords(1, 2) - xcoords(1, 1));
//   y_13 = (xcoords(1, 0) - xcoords(1, 2));
//   y_31 = (xcoords(1, 2) - xcoords(1, 0));

//   double l_12, l_23, l_31;

//   l_12 = sqrt(pow(x_12,2) + pow(y_12,2));
//   l_23 = sqrt(pow(x_23,2) + pow(y_23,2));
//   l_31 = sqrt(pow(x_31,2) + pow(y_31,2));

//   double a_12, b_12, a_23, b_23, a_31, b_31;

//   a_12 = 2.0*Area/l_12;
//   b_12 = (x_12*x_13 + y_12*y_13)/l_12;
//   a_23 = 2.0*Area/l_23;
//   b_23 = (x_23*x_21 + y_23*y_21)/l_23;
//   a_31 = 2.0*Area/l_31;
//   b_31 = (x_31*x_32 + y_31*y_32)/l_31;

//   double xp_12, yp_12, xp_23, yp_23, xp_31, yp_31;

//   xp_12 = x_21*b_12/l_12 + x_1;
//   yp_12 = y_21*b_12/l_12 + y_1;
//   xp_23 = x_32*b_23/l_23 + x_2;
//   yp_23 = y_32*b_23/l_23 + y_2;
//   xp_31 = x_13*b_31/l_31 + x_3;
//   yp_31 = y_13*b_31/l_31 + y_3;

//   double g1, g2, g3, g4, g5, g6;

//   g1 = l_12*(xp_12 - x_3)/a_12;
//   g2 = l_12*(yp_12 - y_3)/a_12;
//   g3 = l_23*(xp_23 - x_1)/a_23;
//   g4 = l_23*(yp_23 - y_1)/a_23;
//   g5 = l_31*(xp_31 - x_2)/a_31;
//   g6 = l_31*(yp_31 - y_2)/a_31;

//   // Get the barycentric coordinates: zeta_1, zeta_2, zeta_3, zeta_12, zeta_23, zeta_31.

//   double zeta_1, zeta_2, zeta_3;

//   zeta_1 = (1.0/(2.0*Area))*((x_2 - x)*(y_3 - y) - (x_3 - x)*(y_2 - y));
//   zeta_2 = (1.0/(2.0*Area))*((x_3 - x)*(y_1 - y) - (x_1 - x)*(y_3 - y));
//   zeta_3 = (1.0/(2.0*Area))*((x_1 - x)*(y_2 - y) - (x_2 - x)*(y_1 - y));

//   // Get the shape functions derivatives with respect to: zeta1, zeta2, zeta3.

//   double dNxx1_dzeta1, dNxx1_dzeta2, dNxx1_dzeta3;

//   dNxx1_dzeta1 = 1.0;
//   dNxx1_dzeta2 = 0.0;
//   dNxx1_dzeta3 = 0.0;

//   double dNxy1_dzeta1, dNxy1_dzeta2, dNxy1_dzeta3;

//   dNxy1_dzeta1 = 0.0;
//   dNxy1_dzeta2 = 0.0;
//   dNxy1_dzeta3 = 0.0;

//   double dNxx2_dzeta1, dNxx2_dzeta2, dNxx2_dzeta3;

//   dNxx2_dzeta1 = 0.0;
//   dNxx2_dzeta2 = 1.0;
//   dNxx2_dzeta3 = 0.0;

//   double dNxy2_dzeta1, dNxy2_dzeta2, dNxy2_dzeta3;

//   dNxy2_dzeta1 = 0.0;
//   dNxy2_dzeta2 = 0.0;
//   dNxy2_dzeta3 = 0.0;

//   double dNxx3_dzeta1, dNxx3_dzeta2, dNxx3_dzeta3;

//   dNxx3_dzeta1 = 0.0;
//   dNxx3_dzeta2 = 0.0;
//   dNxx3_dzeta3 = 1.0;

//   double dNxy3_dzeta1, dNxy3_dzeta2, dNxy3_dzeta3;

//   dNxy3_dzeta1 = 0.0;
//   dNxy3_dzeta2 = 0.0;
//   dNxy3_dzeta3 = 0.0;

//   double dNxtheta1_dzeta1, dNxtheta1_dzeta2, dNxtheta1_dzeta3;

//   dNxtheta1_dzeta1 = (1.0/2.0)*(-g1*zeta_2 + g5*zeta_3);
//   dNxtheta1_dzeta2 = (1.0/2.0)*(-g1*zeta_1);
//   dNxtheta1_dzeta3 = (1.0/2.0)*( g5*zeta_1);

//   double dNxtheta2_dzeta1, dNxtheta2_dzeta2, dNxtheta2_dzeta3;

//   dNxtheta2_dzeta1 = (1.0/2.0)*(g1*zeta_2);
//   dNxtheta2_dzeta2 = (1.0/2.0)*(-g3*zeta_3 + g1*zeta_1);
//   dNxtheta2_dzeta3 = (1.0/2.0)*(-g3*zeta_2);

//   double dNxtheta3_dzeta1, dNxtheta3_dzeta2, dNxtheta3_dzeta3;

//   dNxtheta3_dzeta1 = (1.0/2.0)*(-g5*zeta_3);
//   dNxtheta3_dzeta2 = (1.0/2.0)*(g3*zeta_3);
//   dNxtheta3_dzeta3 = (1.0/2.0)*(-g5*zeta_1 + g3*zeta_2);

//   double dNyx1_dzeta1, dNyx1_dzeta2, dNyx1_dzeta3;

//   dNyx1_dzeta1 = 0.0;
//   dNyx1_dzeta2 = 0.0;
//   dNyx1_dzeta3 = 0.0;

//   double dNyy1_dzeta1, dNyy1_dzeta2, dNyy1_dzeta3;

//   dNyy1_dzeta1 = 1.0;
//   dNyy1_dzeta2 = 0.0;
//   dNyy1_dzeta3 = 0.0;

//   double dNyx2_dzeta1, dNyx2_dzeta2, dNyx2_dzeta3;

//   dNyx2_dzeta1 = 0.0;
//   dNyx2_dzeta2 = 0.0;
//   dNyx2_dzeta3 = 0.0;

//   double dNyy2_dzeta1, dNyy2_dzeta2, dNyy2_dzeta3;

//   dNyy2_dzeta1 = 0.0;
//   dNyy2_dzeta2 = 1.0;
//   dNyy2_dzeta3 = 0.0;

//   double dNyx3_dzeta1, dNyx3_dzeta2, dNyx3_dzeta3;

//   dNyx3_dzeta1 = 0.0;
//   dNyx3_dzeta2 = 0.0;
//   dNyx3_dzeta3 = 0.0;

//   double dNyy3_dzeta1, dNyy3_dzeta2, dNyy3_dzeta3;

//   dNyy3_dzeta1 = 0.0;
//   dNyy3_dzeta2 = 0.0;
//   dNyy3_dzeta3 = 1.0;

//   double dNytheta1_dzeta1, dNytheta1_dzeta2, dNytheta1_dzeta3;

//   dNytheta1_dzeta1 = (1.0/2.0)*(-g2*zeta_2 + g6*zeta_3);
//   dNytheta1_dzeta2 = (1.0/2.0)*(-g2*zeta_1);
//   dNytheta1_dzeta3 = (1.0/2.0)*( g6*zeta_1);

//   double dNytheta2_dzeta1, dNytheta2_dzeta2, dNytheta2_dzeta3;

//   dNytheta2_dzeta1 = (1.0/2.0)*(g2*zeta_2);
//   dNytheta2_dzeta2 = (1.0/2.0)*(-g4*zeta_3 + g2*zeta_1);
//   dNytheta2_dzeta3 = (1.0/2.0)*(-g4*zeta_2);

//   double dNytheta3_dzeta1, dNytheta3_dzeta2, dNytheta3_dzeta3;

//   dNytheta3_dzeta1 = (1.0/2.0)*(-g6*zeta_3);
//   dNytheta3_dzeta2 = (1.0/2.0)*(g4*zeta_3);
//   dNytheta3_dzeta3 = (1.0/2.0)*(-g6*zeta_1 + g4*zeta_2);

//   // Get the jacobians: dzetai_dx, dzetai_dy.

//   double dzeta1_dx, dzeta2_dx, dzeta3_dx;

//   dzeta1_dx = (y_2 - y_3)/(2.0*Area);
//   dzeta2_dx = (y_3 - y_1)/(2.0*Area);
//   dzeta3_dx = (y_1 - y_2)/(2.0*Area);

//   double dzeta1_dy, dzeta2_dy, dzeta3_dy;

//   dzeta1_dy = (x_3 - x_2)/(2.0*Area);
//   dzeta2_dy = (x_1 - x_3)/(2.0*Area);
//   dzeta3_dy = (x_2 - x_1)/(2.0*Area);

//   // Get the shape functions derivatives with respect to: x, y.

//   // Displacement u shape functions derivatives

//   double dNxx1_dx, dNxy1_dx, dNxx2_dx, dNxy2_dx, dNxx3_dx, dNxy3_dx;

//   dNxx1_dx = dNxx1_dzeta1*dzeta1_dx + dNxx1_dzeta2*dzeta2_dx + dNxx1_dzeta3*dzeta3_dx;
//   dNxy1_dx = dNxy1_dzeta1*dzeta1_dx + dNxy1_dzeta2*dzeta2_dx + dNxy1_dzeta3*dzeta3_dx;
//   dNxx2_dx = dNxx2_dzeta1*dzeta1_dx + dNxx2_dzeta2*dzeta2_dx + dNxx2_dzeta3*dzeta3_dx;
//   dNxy2_dx = dNxy2_dzeta1*dzeta1_dx + dNxy2_dzeta2*dzeta2_dx + dNxy2_dzeta3*dzeta3_dx;
//   dNxx3_dx = dNxx3_dzeta1*dzeta1_dx + dNxx3_dzeta2*dzeta2_dx + dNxx3_dzeta3*dzeta3_dx;
//   dNxy3_dx = dNxy3_dzeta1*dzeta1_dx + dNxy3_dzeta2*dzeta2_dx + dNxy3_dzeta3*dzeta3_dx;

//   double dNxx1_dy, dNxy1_dy, dNxx2_dy, dNxy2_dy, dNxx3_dy, dNxy3_dy;

//   dNxx1_dy = dNxx1_dzeta1*dzeta1_dy + dNxx1_dzeta2*dzeta2_dy + dNxx1_dzeta3*dzeta3_dy;
//   dNxy1_dy = dNxy1_dzeta1*dzeta1_dy + dNxy1_dzeta2*dzeta2_dy + dNxy1_dzeta3*dzeta3_dy;
//   dNxx2_dy = dNxx2_dzeta1*dzeta1_dy + dNxx2_dzeta2*dzeta2_dy + dNxx2_dzeta3*dzeta3_dy;
//   dNxy2_dy = dNxy2_dzeta1*dzeta1_dy + dNxy2_dzeta2*dzeta2_dy + dNxy2_dzeta3*dzeta3_dy;
//   dNxx3_dy = dNxx3_dzeta1*dzeta1_dy + dNxx3_dzeta2*dzeta2_dy + dNxx3_dzeta3*dzeta3_dy;
//   dNxy3_dy = dNxy3_dzeta1*dzeta1_dy + dNxy3_dzeta2*dzeta2_dy + dNxy3_dzeta3*dzeta3_dy;

//   double dNxtheta1_dx, dNxtheta2_dx, dNxtheta3_dx;

//   dNxtheta1_dx = dNxtheta1_dzeta1*dzeta1_dx + dNxtheta1_dzeta2*dzeta2_dx + dNxtheta1_dzeta3*dzeta3_dx;
//   dNxtheta2_dx = dNxtheta2_dzeta1*dzeta1_dx + dNxtheta2_dzeta2*dzeta2_dx + dNxtheta2_dzeta3*dzeta3_dx;
//   dNxtheta3_dx = dNxtheta3_dzeta1*dzeta1_dx + dNxtheta3_dzeta2*dzeta2_dx + dNxtheta3_dzeta3*dzeta3_dx;

//   double dNxtheta1_dy, dNxtheta2_dy, dNxtheta3_dy;

//   dNxtheta1_dy = dNxtheta1_dzeta1*dzeta1_dy + dNxtheta1_dzeta2*dzeta2_dy + dNxtheta1_dzeta3*dzeta3_dy;
//   dNxtheta2_dy = dNxtheta2_dzeta1*dzeta1_dy + dNxtheta2_dzeta2*dzeta2_dy + dNxtheta2_dzeta3*dzeta3_dy;
//   dNxtheta3_dy = dNxtheta3_dzeta1*dzeta1_dy + dNxtheta3_dzeta2*dzeta2_dy + dNxtheta3_dzeta3*dzeta3_dy;

//   // Displacement v shape functions derivatives

//   double dNyx1_dx, dNyy1_dx, dNyx2_dx, dNyy2_dx, dNyx3_dx, dNyy3_dx;

//   dNyx1_dx = dNyx1_dzeta1*dzeta1_dx + dNyx1_dzeta2*dzeta2_dx + dNyx1_dzeta3*dzeta3_dx;
//   dNyy1_dx = dNyy1_dzeta1*dzeta1_dx + dNyy1_dzeta2*dzeta2_dx + dNyy1_dzeta3*dzeta3_dx;
//   dNyx2_dx = dNyx2_dzeta1*dzeta1_dx + dNyx2_dzeta2*dzeta2_dx + dNyx2_dzeta3*dzeta3_dx;
//   dNyy2_dx = dNyy2_dzeta1*dzeta1_dx + dNyy2_dzeta2*dzeta2_dx + dNyy2_dzeta3*dzeta3_dx;
//   dNyx3_dx = dNyx3_dzeta1*dzeta1_dx + dNyx3_dzeta2*dzeta2_dx + dNyx3_dzeta3*dzeta3_dx;
//   dNyy3_dx = dNyy3_dzeta1*dzeta1_dx + dNyy3_dzeta2*dzeta2_dx + dNyy3_dzeta3*dzeta3_dx;

//   double dNyx1_dy, dNyy1_dy, dNyx2_dy, dNyy2_dy, dNyx3_dy, dNyy3_dy;

//   dNyx1_dy = dNyx1_dzeta1*dzeta1_dy + dNyx1_dzeta2*dzeta2_dy + dNyx1_dzeta3*dzeta3_dy;
//   dNyy1_dy = dNyy1_dzeta1*dzeta1_dy + dNyy1_dzeta2*dzeta2_dy + dNyy1_dzeta3*dzeta3_dy;
//   dNyx2_dy = dNyx2_dzeta1*dzeta1_dy + dNyx2_dzeta2*dzeta2_dy + dNyx2_dzeta3*dzeta3_dy;
//   dNyy2_dy = dNyy2_dzeta1*dzeta1_dy + dNyy2_dzeta2*dzeta2_dy + dNyy2_dzeta3*dzeta3_dy;
//   dNyx3_dy = dNyx3_dzeta1*dzeta1_dy + dNyx3_dzeta2*dzeta2_dy + dNyx3_dzeta3*dzeta3_dy;
//   dNyy3_dy = dNyy3_dzeta1*dzeta1_dy + dNyy3_dzeta2*dzeta2_dy + dNyy3_dzeta3*dzeta3_dy;

//   double dNytheta1_dx, dNytheta2_dx, dNytheta3_dx;

//   dNytheta1_dx = dNytheta1_dzeta1*dzeta1_dx + dNytheta1_dzeta2*dzeta2_dx + dNytheta1_dzeta3*dzeta3_dx;
//   dNytheta2_dx = dNytheta2_dzeta1*dzeta1_dx + dNytheta2_dzeta2*dzeta2_dx + dNytheta2_dzeta3*dzeta3_dx;
//   dNytheta3_dx = dNytheta3_dzeta1*dzeta1_dx + dNytheta3_dzeta2*dzeta2_dx + dNytheta3_dzeta3*dzeta3_dx;

//   double dNytheta1_dy, dNytheta2_dy, dNytheta3_dy;

//   dNytheta1_dy = dNytheta1_dzeta1*dzeta1_dy + dNytheta1_dzeta2*dzeta2_dy + dNytheta1_dzeta3*dzeta3_dy;
//   dNytheta2_dy = dNytheta2_dzeta1*dzeta1_dy + dNytheta2_dzeta2*dzeta2_dy + dNytheta2_dzeta3*dzeta3_dy;
//   dNytheta3_dy = dNytheta3_dzeta1*dzeta1_dy + dNytheta3_dzeta2*dzeta2_dy + dNytheta3_dzeta3*dzeta3_dy;

//   // Compute the B_mem matrix for the Allman membrane

//   B_mem = 0.0;
//   B_mem(0,0) = dNxx1_dx;
//   B_mem(0,1) = dNxy1_dx; 
//   B_mem(0,2) = dNxtheta1_dx; 
//   B_mem(0,3) = dNxx2_dx;
//   B_mem(0,4) = dNxy2_dx; 
//   B_mem(0,5) = dNxtheta2_dx; 
//   B_mem(0,6) = dNxx3_dx;
//   B_mem(0,7) = dNxy3_dx; 
//   B_mem(0,8) = dNxtheta3_dx; 

//   B_mem(1,0) = dNyx1_dy;
//   B_mem(1,1) = dNyy1_dy; 
//   B_mem(1,2) = dNytheta1_dy; 
//   B_mem(1,3) = dNyx2_dy;
//   B_mem(1,4) = dNyy2_dy; 
//   B_mem(1,5) = dNytheta2_dy; 
//   B_mem(1,6) = dNyx3_dy;
//   B_mem(1,7) = dNyy3_dy; 
//   B_mem(1,8) = dNytheta3_dy; 

//   B_mem(2,0) = dNxx1_dy + dNyx1_dx;
//   B_mem(2,1) = dNxy1_dy + dNyy1_dx; 
//   B_mem(2,2) = dNxtheta1_dy + dNytheta1_dx; 
//   B_mem(2,3) = dNxx2_dy + dNyx2_dx;
//   B_mem(2,4) = dNxy2_dy + dNyy2_dx; 
//   B_mem(2,5) = dNxtheta2_dy + dNytheta2_dx; 
//   B_mem(2,6) = dNxx3_dy + dNyx3_dx;
//   B_mem(2,7) = dNxy3_dy + dNyy3_dx; 
//   B_mem(2,8) = dNxtheta3_dy + dNytheta3_dx; 
// }


//-----------------------------------------------------------------------
//   getHmatAnalytic_
//-----------------------------------------------------------------------


void Shell6DofsAllmanModel::getHmatAnalytic_

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


void Shell6DofsAllmanModel::getHmatNumeric_

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


void Shell6DofsAllmanModel::getBmat_

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


void Shell6DofsAllmanModel::getDissForce_

  ( const Ref<Plasticity> p,
    const Vector&         fstar,
    const Vector&         disp )   const

{
  Matrix      coords      ( nodeRank_, nodeCount_ );
  Matrix      transMat    ( nodeRank_, nodeRank_ );
  Matrix      transMatDof ( dofCount_, dofCount_ );
  Matrix      xcoords     ( rank_, nodeCount_ );
  Matrix      ipxcoords   ( rank_, ipCount_ );
  Matrix      B_mem           ( strCount_, dofCount_  );
  Matrix      Bt          = B_mem.transpose ();
  idx_t       ipoint      = 0;

  Cubix       grads       ( rank_, nodeCount_, ipCount_ );
  Vector      ipWeights   ( ipCount_ );

  IdxVector   inodes      ( nodeCount_ );
  IdxVector   idofs       ( dofCount_  );
  Vector      strain      ( strCount_  );
  Vector      elemForce   ( dofCount_  );
  Vector      elemDisp    ( dofCount_  );
  Vector      elemxDisp   ( dofCount_ );
  Vector      Disp_mem    ( 3*nodeCount_ );
  Vector      sstar       ( strCount_  );

  double      Area;

  MChain1     mc1;

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

    shape_->getGlobalIntegrationPoints ( ipxcoords, xcoords );

    getArea_(Area, xcoords);

    elemDisp = disp[idofs];

    elemxDisp = mc1.matmul ( transMatDof, elemDisp );

    getMembraneDisp_(Disp_mem, elemxDisp);

    ipWeights *= thickness_;

    elemForce = 0.0;

    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    {
      double x = ipxcoords(0,ip);
      double y = ipxcoords(1,ip);
  
      getBmemAllman_ (B_mem, xcoords, Area, x, y);

      // Compute the strain vector of this integration point

      matmul ( strain, B_mem, Disp_mem );

      // get dissipation stress F^T*sigma+D^T*eps^p

      p->getDissipationStress ( sstar, strain, ipoint++ );

      elemForce += ipWeights[ip] * matmul ( Bt, sstar );
    }
    // Add the element force vector to the global force vector.

    fstar[idofs] += elemForce;
  }
}


//-----------------------------------------------------------------------
//   initCharLength_       NOT DONE FOR THE SHELL MODEL YET
//-----------------------------------------------------------------------

void Shell6DofsAllmanModel::initCharLength_ ()

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

void Shell6DofsAllmanModel::getDissipation_

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

void Shell6DofsAllmanModel::getIshape_ 

  ( Ref<IShape> shape ) const

{
  shape = shape_;
}

//-----------------------------------------------------------------------
//   getMaterial_
//-----------------------------------------------------------------------

void Shell6DofsAllmanModel::getMaterial_

  ( Ref<Material> material ) const 

{
  material = material_;
}

//-----------------------------------------------------------------------
//   getThickness_
//-----------------------------------------------------------------------

void Shell6DofsAllmanModel::getThickness_

  ( double thickness ) const 

{
  thickness = thickness_;
}

//-----------------------------------------------------------------------
//   getOrientVec_
//-----------------------------------------------------------------------

void Shell6DofsAllmanModel::getOrientVec_

  ( const double orientVec[3] ) const 

{
  orientVec = orientVec_;
}

//-----------------------------------------------------------------------
//   getStress_       NOT DONE FOR THE SHELL MODEL YET
//-----------------------------------------------------------------------


void Shell6DofsAllmanModel::getStress_

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
  Matrix     ipxcoords   ( rank_, ipCount_ );
  Matrix     B_mem       ( strCount_, dofCount_  );

  Vector     nStressIp   ( strCount_ );    // normal stress vector at idx_t.pt.
  Vector     strain      ( strCount_ );
  Vector     elemDisp    ( dofCount_ );
  Vector     elemxDisp   ( dofCount_ );
  Vector     Disp_mem    ( 3*nodeCount_ );

  IdxVector  inodes      ( nodeCount_ );
  IdxVector  idofs       ( dofCount_  );
  IdxVector  jcols       ( strCount_  );

  double    Area;

  MChain1     mc1;

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

    shape_->getGlobalIntegrationPoints ( ipxcoords, xcoords );

    getArea_(Area, xcoords);

    elemDisp = disp[idofs];

    elemxDisp = mc1.matmul ( transMatDof, elemDisp );

    getMembraneDisp_(Disp_mem, elemxDisp);

    // Matrix     sfuncs     = shape_->getShapeFunctions ();

    // Iterate over the integration points.

    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    {
      double x = ipxcoords(0,ip);
      double y = ipxcoords(1,ip);
  
      getBmemAllman_ (B_mem, xcoords, Area, x, y);
  
      matmul ( strain, B_mem, Disp_mem );

      if ( crackBandMethod_ )
      {
        double le = charLength_[ielem];
        softening_->update ( nStressIp, stiff, strain, ipoint, le );
      }
      else
      {
        material_->update ( nStressIp, stiff, strain, ipoint );
      }

      // ndNStress += matmul ( sfuncs(ALL,ip), nStressIp );

      // ndWeights += sfuncs(ALL,ip);

      ++ipoint; 
    }

    // select ( weights, inodes ) += ndWeights;

    // Add the stresses to the table.

    table.addBlock ( inodes, jcols[slice(0,strCount_)],   ndNStress );
  }
}


//-----------------------------------------------------------------------
//    writeElements_       NOT DONE FOR THE SHELL MODEL YET
//-----------------------------------------------------------------------

void Shell6DofsAllmanModel::writeElements_

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

void Shell6DofsAllmanModel::writeNodes_

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

void  Shell6DofsAllmanModel::initLocOutWriter_ ()

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

Ref<PrintWriter>  Shell6DofsAllmanModel::initWriter_

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


void Shell6DofsAllmanModel::writeDisplacements_

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


bool Shell6DofsAllmanModel::getTable_

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
//   newShell6DofsAllmanModel
//-----------------------------------------------------------------------


static Ref<Model>     newShell6DofsAllmanModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<Shell6DofsAllmanModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareShell6DofsAllmanModel
//-----------------------------------------------------------------------


void declareShell6DofsAllmanModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "Shell6DofsAllman", & newShell6DofsAllmanModel );
}