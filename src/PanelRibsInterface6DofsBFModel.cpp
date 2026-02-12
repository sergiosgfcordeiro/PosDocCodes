/*
 *
 * Model for Panel-Ribs interface elements with 6 DOFs per node: 
 *   - assembly of stiffness and internal force vector
 *   - output
 * 
 * Sergio Gustavo Ferreira Cordeiro, July 2025
 *
 */

#include "PanelRibsInterface6DofsBFModel.h"

#include <jem/base/array/operators.h>
#include <jem/base/array/select.h>
#include <jem/base/Error.h>
#include <jem/base/Float.h>
#include <jem/base/IllegalInputException.h>
#include <jem/base/System.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/numeric/algebra/MatmulChain.h>
#include <jem/numeric/algebra/Cholesky.h>
#include <jem/numeric/algebra/LUSolver.h>
#include <jem/util/StringUtils.h>
#include <jive/model/ModelFactory.h>
#include <jive/model/Actions.h>
#include <jive/model/StateVector.h>
#include <jive/geom/BoundaryShape.h>
#include <jive/geom/BShapeFactory.h>
#include <jive/geom/BoundaryTriangle.h>
#include <jive/geom/BoundaryLine.h>
#include <jive/fem/Globdat.h>

#include "models.h"
#include "SolverNames.h"

using jem::numeric::norm2;
using jem::util::StringUtils;
using jem::numeric::matmul;
using jem::numeric::MatmulChain;
using jive::model::StateVector;
using jive::geom::BShapeFactory;
using jive::geom::BShape;
using jive::fem::Globdat;


//-----------------------------------------------------------------------
//   static constants
//-----------------------------------------------------------------------


const char* PanelRibsInterface6DofsBFModel::DOF_NAMES[6]  = {"u","v","w","wx","wy", "wz"};
const char* PanelRibsInterface6DofsBFModel::SHAPE_PROP    = "shape";
const char* PanelRibsInterface6DofsBFModel::COHEMAT_PROP  = "coheMat";
const char* PanelRibsInterface6DofsBFModel::PANRIB_SHAPE_PROP = "shape";
const char* PanelRibsInterface6DofsBFModel::MATERIAL_PROP[2]  = {"material","material"};
const char* PanelRibsInterface6DofsBFModel::THICK_PROP[2]     = {"thickness","thickness"};
const char* PanelRibsInterface6DofsBFModel::THICK_SHAPE_PROP  = "thick_shape";
const char* PanelRibsInterface6DofsBFModel::PAN_ORIENTATION_PROP[3] = {"v1","v2","v3"};  
const char* PanelRibsInterface6DofsBFModel::RIB_ORIENTATION_PROP[3] = {"v1","v2","v3"}; 
const char* PanelRibsInterface6DofsBFModel::PAN_NAME_PROP = "panel_name";
const char* PanelRibsInterface6DofsBFModel::RIB_NAME_PROP = "rib_name";   


//-----------------------------------------------------------------------
//   constructor
//-----------------------------------------------------------------------


PanelRibsInterface6DofsBFModel::PanelRibsInterface6DofsBFModel

   ( const String&       name,
     const Properties&   conf,
     const Properties&   props,
     const Properties&   globdat ) : Super(name)
{
  Ref<BShape> bshape = NIL;

  using jive::util::joinNames;
  using jive::geom::IShapeFactory;

  // create myTag_ (last part of myName_)
  
  StringVector names ( StringUtils::split( myName_, '.' ) );
  myTag_     = names [ names.size() - 1 ];
  myBase_ = myName_[slice(0,myName_.size() - myTag_.size()-1)];

  Properties  myProps = props.getProps ( myName_ );
  Properties  myConf  = conf.makeProps ( myName_ );
  Properties  shProps = myProps.getProps ( SHAPE_PROP );

  const String context = getContext();

  egroup_ = ElemGroup::get ( myConf, myProps, globdat, context );

  elems_      = egroup_.getElements ( );
  elemCount_  = egroup_.size        ( );

  ielems_.resize ( elemCount_ );
  ielems_ = egroup_.getIndices ( );

  nodes_      = elems_.getNodes     ( );
  nodeRank_   = nodes_.rank         ( );
  rank_       = nodeRank_ - 1;
  thick_rank_ = nodeRank_ - 2;

  myConf.set ( SHAPE_PROP, shProps );

  bshape = BShapeFactory::newInstance
    ( joinNames (myName_, SHAPE_PROP ), conf, props );

  shape_ = newInstance<InterfaceShape> ( "interface", bshape );

  thick_shape_  = IShapeFactory::newInstance(
    joinNames (myName_, THICK_SHAPE_PROP ),
    conf,
    props );

  // NB: Make sure you know what you're doing when 
  // the rank of the shape doesn't match the rank of the mesh

  qRank_   = nodeRank_;
  nodeCount_BLine2_ = shape_->nodeCount   ( );
  ipCount_   = shape_->ipointCount ( );
  nodeCount_Line2_ = thick_shape_->nodeCount   ( );
  jpCount_    = thick_shape_->ipointCount ();

  getShapeFuncs_ = getShapeFuncsFunc ( rank_ );

  // Create a coheMat object.

  coheMat_ = newCohesiveMat 
           ( COHEMAT_PROP, myConf, myProps, globdat );

  frictionMat_ = dynamicCast<AlfanoTuronCoheMat> ( coheMat_ );

  // Get the panel and the ribs Names.

  myProps.find( panel_, PAN_NAME_PROP );
  myConf.set  ( PAN_NAME_PROP, panel_ );
  String panelName = joinNames (myBase_, panel_ );
  
  myProps.find( ribs_, RIB_NAME_PROP );
  myConf.set  ( RIB_NAME_PROP, ribs_ );
  String ribsName = joinNames (myBase_, ribs_ );

  // Get the elements from panel and ribs submodels.

  Properties  panelProps = props.getProps ( panelName );
  Properties  panelConf  = conf.makeProps ( panelName );
  const String panelContext = "model `" + panelName + "'";

  pan_egroup_ = ElemGroup::get ( panelConf, panelProps, globdat, panelContext );
  pan_elems_  = pan_egroup_.getElements ( );
  pan_elemCount_ = pan_egroup_.size ( );
  pan_ielems_.resize ( pan_elemCount_ );
  pan_ielems_ = pan_egroup_.getIndices ( );
  pan_nodes_  = pan_elems_.getNodes ( );

  Properties  ribsProps = props.getProps ( ribsName );
  Properties  ribsConf  = conf.makeProps ( ribsName );
  const String ribsContext = "model `" + ribsName + "'";

  rib_egroup_ = ElemGroup::get ( ribsConf, ribsProps, globdat, ribsContext );
  rib_elems_  = rib_egroup_.getElements ( );
  rib_elemCount_  = rib_egroup_.size ( );
  rib_ielems_.resize ( rib_elemCount_ );
  rib_ielems_ = rib_egroup_.getIndices ( );
  rib_nodes_  = rib_elems_.getNodes ( );

  pan_dofs_ = XDofSpace::get ( pan_nodes_.getData(), globdat );
  rib_dofs_ = XDofSpace::get ( rib_nodes_.getData(), globdat );
  pan_dofTypes_.resize( 6 );
  rib_dofTypes_.resize( 6 );

  for( idx_t i = 0; i < 6; i++)
  {
    pan_dofTypes_[i] = pan_dofs_->addType ( DOF_NAMES[i]);
    rib_dofTypes_[i] = rib_dofs_->addType ( DOF_NAMES[i]);
  }

  pan_dofs_->addDofs ( pan_elems_.getUniqueNodesOf ( pan_ielems_ ), pan_dofTypes_);
  rib_dofs_->addDofs ( rib_elems_.getUniqueNodesOf ( rib_ielems_ ), rib_dofTypes_);

  // Create the InternalShape object for the panel and the ribs.
  
  Properties  panrib_shProps = panelProps.getProps ( PANRIB_SHAPE_PROP );
  panelConf.set ( PANRIB_SHAPE_PROP, panrib_shProps );
  
  panribshape_  = IShapeFactory::newInstance(PANRIB_SHAPE_PROP, panelConf, panelProps );
  
  nodeCount_panrib_ = panribshape_->nodeCount  ( );
  nodeCount_ = 3 * nodeCount_panrib_;
  ipCount_panrib_ = panribshape_->ipointCount ( );

  getShapeGrads_ = getShapeGradsFunc ( rank_ );
  
  // Create two material model objects: panel and ribs.
  
  material_[0] = newMaterial ( MATERIAL_PROP[0], panelConf, panelProps, globdat );
  material_[0]-> allocPoints  ( ipCount_panrib_ );
  softening_[0] = dynamicCast<Softening> ( material_[0] );
  
  material_[1] = newMaterial ( MATERIAL_PROP[1], ribsConf, ribsProps, globdat );
  material_[1]-> allocPoints  ( ipCount_panrib_ );
  softening_[1] = dynamicCast<Softening> ( material_[1] );
  
  // Create two thickness variables: panel and ribs.
  
  panelProps.find( thickness_[0], THICK_PROP[0] );
  panelConf.set  ( THICK_PROP[0], thickness_[0] );
  
  ribsProps.find( thickness_[1], THICK_PROP[1] );
  ribsConf.set  ( THICK_PROP[1], thickness_[1] );
  
  // Create the variables: pan_orientVec_ and rib_orientVec_.
  
  panelProps.find( pan_orientVec_[0], PAN_ORIENTATION_PROP[0] );
  panelConf.set  ( PAN_ORIENTATION_PROP[0], pan_orientVec_[0] );
  
  panelProps.find( pan_orientVec_[1], PAN_ORIENTATION_PROP[1] );
  panelConf.set  ( PAN_ORIENTATION_PROP[1], pan_orientVec_[1] );
  
  panelProps.find( pan_orientVec_[2], PAN_ORIENTATION_PROP[2] );
  panelConf.set  ( PAN_ORIENTATION_PROP[2], pan_orientVec_[2] );
  
  ribsProps.find( rib_orientVec_[0], RIB_ORIENTATION_PROP[0] );
  ribsConf.set  ( RIB_ORIENTATION_PROP[0], rib_orientVec_[0] );
  
  ribsProps.find( rib_orientVec_[1], RIB_ORIENTATION_PROP[1] );
  ribsConf.set  ( RIB_ORIENTATION_PROP[1], rib_orientVec_[1] );
  
  ribsProps.find( rib_orientVec_[2], RIB_ORIENTATION_PROP[2] );
  ribsConf.set  ( RIB_ORIENTATION_PROP[2], rib_orientVec_[2] );
  
  crackBandMethod_[0] = false;
  crackBandMethod_[1] = false;
}

PanelRibsInterface6DofsBFModel::~PanelRibsInterface6DofsBFModel()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  Properties  myProps  = props.findProps ( myName_ );
  Properties  matProps = myProps.findProps ( COHEMAT_PROP );

  coheMat_->configure ( matProps, globdat );

  coheMat_->allocPoints  ( ipCount_ * jpCount_ * elemCount_ );

  // Get the panelName.

  String panelName = joinNames (myBase_, panel_ );

  Properties  panelProps  = props.findProps ( panelName );
  Properties  panelMatProps = panelProps.findProps ( MATERIAL_PROP[0] );

  material_[0]->configure ( panelMatProps );

  // Get the ribsName.

  String ribsName = joinNames (myBase_, ribs_ );

  Properties  ribsProps  = props.findProps ( ribsName );
  Properties  ribsMatProps = ribsProps.findProps ( MATERIAL_PROP[1] );

  material_[1]->configure ( ribsMatProps );

  Properties  params;
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::getConfig 

  ( const Properties& conf,
    const Properties& globdat ) const

{
  Properties  myConf  = conf.makeProps ( myName_ );
  Properties  matConf = myConf.makeProps ( COHEMAT_PROP );

  coheMat_->getConfig ( matConf, globdat );

  // Get the panelName.

  String panelName = joinNames (myBase_, panel_ );

  Properties  panelConf  = conf.makeProps ( panelName );
  Properties  panelMatConf = panelConf.makeProps ( MATERIAL_PROP[0] );

  material_[0]->getConfig ( panelMatConf );

  // Get the ribsName.

  String ribsName = joinNames (myBase_, ribs_ );

  Properties  ribsConf  = conf.makeProps ( ribsName );
  Properties  ribsMatConf = ribsConf.makeProps ( MATERIAL_PROP[1] );

  material_[1]->getConfig ( ribsMatConf );
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool PanelRibsInterface6DofsBFModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;

  if ( action == Actions::GET_MATRIX0 ||
       action == Actions::GET_INT_VECTOR )
  {
    Ref<MatrixBuilder>  mbuilder;

    Vector  pan_disp_;
    Vector  rib_disp_;
    Vector  force;

    // Get the current displacements.

    StateVector::get ( pan_disp_, pan_dofs_, globdat );
    StateVector::get ( rib_disp_, rib_dofs_, globdat );

    // Get the matrix builder and the internal force vector.

    params.find( mbuilder, ActionParams::MATRIX0 );
    params.get ( force,    ActionParams::INT_VECTOR );

    getMatrix_ ( mbuilder, force, pan_disp_, rib_disp_ );

    return true;
  }

  if ( action == Actions::COMMIT )
  {
    coheMat_->commit ();
    return true;
  }

  if ( action == "WRITE_XOUTPUT" )
  {
    if ( xOut_ == NIL )
    {
      initWriter_ ( params );
      writeGeom_();
    }
    writeOutput_ ( globdat );
    return true;
  }

  if ( action == "GET_DISSIPATION" )
  {
    getDissipation_ ( params );
    return true;
  }

  if ( action == SolverNames::GET_DISS_FORCE )
  {
    if ( frictionMat_ != NIL )
    {
      Vector pan_disp;
      Vector rib_disp;
      Vector fDiss;

      StateVector::getOld ( pan_disp, pan_dofs_, globdat );
      StateVector::getOld ( rib_disp, rib_dofs_, globdat );
      globdat.get ( fDiss, SolverNames::DISSIPATION_FORCE );

      getFrictionForce_ ( fDiss, pan_disp, rib_disp );
    }
  }

  if ( action == "DESPAIR" )
  {
    return coheMat_->despair();
  }

  if ( action == "END_DESPAIR" )
  {
    coheMat_->endDespair();

    return true;
  }

  return false;
}


//-----------------------------------------------------------------------
//   getMatrix_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::getMatrix_

  ( Ref<MatrixBuilder>  mbuilder,
    const Vector&       force,
    const Vector&       pan_disp,
    const Vector&       rib_disp ) const

{
  const idx_t dofCount   = 6 * nodeCount_;

  Matrix      stiff            ( qRank_, qRank_ );
  Matrix      BCDmat           ( qRank_, dofCount );
  Matrix      BCTmat           ( qRank_, dofCount );
  Matrix      BCmat            ( qRank_, dofCount );

  Matrix      stiff_panA       ( qRank_, qRank_ );
  Matrix      stiff_panB       ( qRank_, qRank_ );
  Matrix      stiff_rib        ( qRank_, qRank_ );
  Matrix      D_plate_panA     ( qRank_, qRank_ );
  Matrix      D_plate_panB     ( qRank_, qRank_ );
  Matrix      D_plate_rib      ( qRank_, qRank_ );

  Matrix      coords           ( qRank_, nodeCount_BLine2_ );
  Matrix      coords_panA      ( qRank_, nodeCount_panrib_ );
  Matrix      coords_panB      ( qRank_, nodeCount_panrib_ );
  Matrix      coords_rib       ( qRank_, nodeCount_panrib_ );
  Matrix      xcoords          ( rank_, nodeCount_BLine2_ );
  Matrix      thick_xcoords    ( thick_rank_, nodeCount_Line2_ );
  Matrix      coords2D         ( rank_, nodeCount_BLine2_ );
  Matrix      xcoords_panA     ( rank_, nodeCount_panrib_ );
  Matrix      xcoords_panB     ( rank_, nodeCount_panrib_ );
  Matrix      xcoords_rib      ( rank_, nodeCount_panrib_ );

  Matrix      ipxcoords        ( rank_, ipCount_ );
  Matrix      jpxcoords        ( thick_rank_, jpCount_ );
  Matrix      ipcoords2D       ( rank_, ipCount_ );
  Matrix      ipxcoords_panA   ( qRank_, ipCount_ );
  Matrix      ipxcoords_panB   ( qRank_, ipCount_ );
  Matrix      ipxcoords_rib    ( qRank_, ipCount_ );

  Matrix      b_mem_panA       ( qRank_, 2*nodeCount_panrib_ );
  Vector      Disp_mem_panA    ( 2*nodeCount_panrib_ );
  Vector      strain_mem_panA  ( qRank_ );
  Vector      stress_mem_panA  ( qRank_ );
  Matrix      Gmat_panA        ( 3*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Ginv_panA        ( 3*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hr_panA          ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hc_panA          ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hh_panA          ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Nr_panA          ( rank_, nodeCount_panrib_ );
  Matrix      Nc_panA          ( rank_, nodeCount_panrib_ );
  Matrix      Nh_panA          ( rank_, nodeCount_panrib_ );
  Matrix      Nrx_panA         ( rank_, nodeCount_panrib_ );
  Matrix      Ncx_panA         ( rank_, nodeCount_panrib_ );
  Matrix      Nhx_panA         ( rank_, nodeCount_panrib_ );
  Matrix      Nry_panA         ( rank_, nodeCount_panrib_ );
  Matrix      Ncy_panA         ( rank_, nodeCount_panrib_ );
  Matrix      Nhy_panA         ( rank_, nodeCount_panrib_ );

  Matrix      b_mem_panB       ( qRank_, 2*nodeCount_panrib_ );
  Vector      Disp_mem_panB    ( 2*nodeCount_panrib_ );
  Vector      strain_mem_panB  ( qRank_ );
  Vector      stress_mem_panB  ( qRank_ );
  Matrix      Gmat_panB        ( 3*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Ginv_panB        ( 3*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hr_panB          ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hc_panB          ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hh_panB          ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Nr_panB          ( rank_, nodeCount_panrib_ );
  Matrix      Nc_panB          ( rank_, nodeCount_panrib_ );
  Matrix      Nh_panB          ( rank_, nodeCount_panrib_ );
  Matrix      Nrx_panB         ( rank_, nodeCount_panrib_ );
  Matrix      Ncx_panB         ( rank_, nodeCount_panrib_ );
  Matrix      Nhx_panB         ( rank_, nodeCount_panrib_ );
  Matrix      Nry_panB         ( rank_, nodeCount_panrib_ );
  Matrix      Ncy_panB         ( rank_, nodeCount_panrib_ );
  Matrix      Nhy_panB         ( rank_, nodeCount_panrib_ );

  Matrix      b_mem_rib        ( qRank_, 2*nodeCount_panrib_ );
  Vector      Disp_mem_rib     ( 2*nodeCount_panrib_ );
  Vector      strain_mem_rib   ( qRank_ );
  Vector      stress_mem_rib   ( qRank_ );
  Matrix      Gmat_rib         ( 3*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Ginv_rib         ( 3*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hr_rib           ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hc_rib           ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hh_rib           ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Nr_rib           ( rank_, nodeCount_panrib_ );
  Matrix      Nc_rib           ( rank_, nodeCount_panrib_ );
  Matrix      Nh_rib           ( rank_, nodeCount_panrib_ );

  Matrix      Bmat_panA        ( 7, 4*nodeCount_panrib_);
  Matrix      Tmat_panA        ( 4*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hmat_panA        ( 7, 7);
  Matrix      Hinv_panA        ( 7, 7);
  Matrix      Bmat_panB        ( 7, 4*nodeCount_panrib_);
  Matrix      Tmat_panB        ( 4*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hmat_panB        ( 7, 7);
  Matrix      Hinv_panB        ( 7, 7);
  Matrix      Bmat_rib         ( 7, 4*nodeCount_panrib_);
  Matrix      Tmat_rib         ( 4*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hmat_rib         ( 7, 7);
  Matrix      Hinv_rib         ( 7, 7);

  Matrix      BAmat            ( nodeCount_panrib_, 3 * nodeCount_panrib_ );

  Matrix      MAmat_panA       ( nodeCount_panrib_, rank_+1 );
  Matrix      MAinv_panA       ( nodeCount_panrib_, rank_+1 );
  Matrix      Mamat_panA       ( nodeCount_panrib_, 7 );
  Matrix      Cmat_panA        ( 7, 3 * nodeCount_panrib_ );
  Matrix      BA_MaC_panA      ( nodeCount_panrib_, 3 * nodeCount_panrib_ );

  Matrix      MAmat_panB       ( nodeCount_panrib_, rank_+1 );
  Matrix      MAinv_panB       ( nodeCount_panrib_, rank_+1 );
  Matrix      Mamat_panB       ( nodeCount_panrib_, 7 );
  Matrix      Cmat_panB        ( 7, 3 * nodeCount_panrib_ );
  Matrix      BA_MaC_panB      ( nodeCount_panrib_, 3 * nodeCount_panrib_ );

  Matrix      MAmat_rib        ( nodeCount_panrib_, rank_+1 );
  Matrix      MAinv_rib        ( nodeCount_panrib_, rank_+1 );
  Matrix      Mamat_rib        ( nodeCount_panrib_, 7 );
  Matrix      Cmat_rib         ( 7, 3 * nodeCount_panrib_ );
  Matrix      BA_MaC_rib       ( nodeCount_panrib_, 3 * nodeCount_panrib_ );

  Matrix      Bxmat            ( 1, rank_+1 );
  Matrix      Bymat            ( 1, rank_+1 );
  Matrix      Smat_panA        ( rank_+1, 1 );
  Matrix      Rmat_panA        ( 7, 1 );
  Matrix      Rxmat_panA       ( 7, 1 );
  Matrix      Rymat_panA       ( 7, 1 );
  Matrix      Smat_panB        ( rank_+1, 1 );
  Matrix      Rmat_panB        ( 7, 1 );
  Matrix      Rxmat_panB       ( 7, 1 );
  Matrix      Rymat_panB       ( 7, 1 );
  Matrix      Smat_rib         ( rank_+1, 1 );
  Matrix      Rmat_rib         ( 7, 1 );
  Matrix      Rxmat_rib        ( 7, 1 );
  Matrix      Rymat_rib        ( 7, 1 );

  Matrix      Nu_panA          ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nv_panA          ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nuy_panA         ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nvx_panA         ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nw_panA          ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nwx_panA         ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nwy_panA         ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nu_panB          ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nv_panB          ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nuy_panB         ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nvx_panB         ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nw_panB          ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nwx_panB         ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nwy_panB         ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nu_rib           ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nv_rib           ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nw_rib           ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nwx_rib          ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nwy_rib          ( 1, 3 * nodeCount_panrib_ );

  Matrix      elemMat      ( dofCount, dofCount  );
  Vector      elemForce    ( dofCount  );
  Matrix      elemxMat     ( dofCount, dofCount  );
  Vector      elemxForce   ( dofCount  );
  Vector      elemxDisp    ( dofCount  );

  IdxVector   Connect      ( dofCount);
  IdxVector   ConnectInv   ( dofCount);

  Matrix      elemMat0     ( dofCount, dofCount  );
  Vector      elemForce0   ( dofCount  );
  Matrix      elemxMat0    ( dofCount, dofCount  );
  Vector      elemxForce0  ( dofCount  );
  Vector      elemDisp0    ( dofCount  );
  Vector      elemxDisp0   ( dofCount  );

  Vector      jump         ( qRank_     );
  Vector      traction     ( qRank_     );
  Vector      trac         ( qRank_     );

  Matrix      transMat_pan ( qRank_, qRank_ );
  Matrix      transMat_rib ( qRank_, qRank_ );
  Matrix      transMat_BL2 ( qRank_, qRank_ );
  Matrix      transMatDof  ( dofCount, dofCount );

  IdxVector   inodes       ( nodeCount_BLine2_ );
  IdxVector   inodes_panA  ( nodeCount_panrib_ );
  IdxVector   inodes_panB  ( nodeCount_panrib_ );
  IdxVector   inodes_rib   ( nodeCount_panrib_ );
  IdxVector   idofs        ( dofCount );
  Vector      ipWeights    ( ipCount_   );
  Vector      jpWeights    ( jpCount_ );

  Cubix       grads_panA       ( rank_, nodeCount_panrib_, ipCount_panrib_ );
  Cubix       grads_panB       ( rank_, nodeCount_panrib_, ipCount_panrib_ );
  Cubix       grads_rib        ( rank_, nodeCount_panrib_, ipCount_panrib_ );
  Vector      ipWeights_panrib ( ipCount_panrib_ );

  Vector      T_I ( ipCount_*jpCount_ );
  Vector      D_I ( ipCount_*jpCount_ );

  Vector      jump_T   ( qRank_     );
  Vector      R_I      ( ipCount_ );

  double Area_panA, Area_panB, Area_rib;
  double det_panA, det_panB, det_rib;

  MChain1     mc1;
  MChain2     mc2;
  MChain3     mc3;

  idx_t       ipoint = 0;
  idx_t       ipoint_panA = 0;
  idx_t       ipoint_panB = 0;
  idx_t       ipoint_rib = 0;

  // Get panel and rib thickness.

  double h_pan = thickness_[0];
  double h_rib = thickness_[1];

  // Get the thick_xcoords for the through the thickness numerical integration.

  getThick_xcoords_(thick_xcoords, h_rib);

  // Get the integration points for the through the thickness numerical integration.

  thick_shape_->getGlobalIntegrationPoints ( jpxcoords, thick_xcoords );

  // Get the Weights for the through the thickness numerical integration.

  thick_shape_->getIntegrationWeights ( jpWeights, thick_xcoords );

  // Dofs connectivity from the Jive order to the shell-shell cohesive line element order.

  getConnectivity_ ( Connect, ConnectInv );

  // Iterate over all elements assigned to this model.

  for ( idx_t ie = 0; ie < elemCount_; ie++ )
  {
    idx_t  ielem = ielems_[ie];
    idx_t  ielem_panSideA;
    idx_t  ielem_panSideB;
    idx_t  ielem_rib;

    // Get the BLine2 element coordinates.

    elems_.getElemNodes  ( inodes, ielem );
    nodes_.getSomeCoords ( coords, inodes );

    // Get the id of the panel and ribs elements connected to the BLine2 element.

    getPanRibElements_(ielem_panSideA, ielem_panSideB, ielem_rib, pan_elems_, rib_elems_, coords);

    // Get the panel and ribs elements coordinates and the DOFs.

    pan_elems_.getElemNodes  ( inodes_panA, ielem_panSideA );
    pan_nodes_.getSomeCoords ( coords_panA, inodes_panA );

    pan_elems_.getElemNodes  ( inodes_panB, ielem_panSideB );
    pan_nodes_.getSomeCoords ( coords_panB, inodes_panB );

    rib_elems_.getElemNodes  ( inodes_rib, ielem_rib );
    rib_nodes_.getSomeCoords ( coords_rib, inodes_rib );

    pan_dofs_->getDofIndices ( idofs[slice(0,(dofCount/3))], inodes_panA, pan_dofTypes_ );
    pan_dofs_->getDofIndices ( idofs[slice((dofCount/3),(2*dofCount/3))], inodes_panB, pan_dofTypes_ );
    rib_dofs_->getDofIndices ( idofs[slice((2*dofCount/3),dofCount)], inodes_rib, rib_dofTypes_ );

    // Get the orientation vectors: v_panel and v_ribs.

    Vector v_pan ( qRank_ );
    Vector v_rib ( qRank_ );
    for ( idx_t i = 0; i < qRank_; i++ )
    {
      v_pan[i] = rib_orientVec_[i];
      v_rib[i] = rib_orientVec_[i];
    }

    // Get the transformation matrices (change of basis) for the structural CE local systems: panel, ribs, BLine2 and interface.

    getTransMatrix_ ( transMat_pan, coords_panA, v_pan );
    getTransMatrix_ ( transMat_rib, coords_rib, v_rib );
    getTransMatrix_ ( transMat_BL2, coords_panA, v_rib );
    getTransMatrixDof_ ( transMatDof, transMat_pan, transMat_rib );

    // Get the nodal coordinates in the Local coordinate systems: Panel A, Panel B, Ribs and BLine2 nodes.

    get2DLocalcoordinates_(xcoords_panA, transMat_pan, coords_panA);
    get2DLocalcoordinates_(xcoords_panB, transMat_pan, coords_panB);
    get2DLocalcoordinates_(xcoords_rib, transMat_rib, coords_rib);
    get2DLocalcoordinates_(xcoords, transMat_BL2, coords);
    get2DGlobalcoordinates_(coords2D, transMat_BL2, coords);

   // Get the areas of the panel A, Panel B and Ribs triangles.

    getArea_(Area_panA, xcoords_panA);
    getArea_(Area_panB, xcoords_panB);
    getArea_(Area_rib, xcoords_rib);

    // Get the G matrices relating generalized Dofs to nodal Dofs: Gmat_panA, Gmat_panB, Gmat_ribs.

    getGmat_(Gmat_panA, xcoords_panA, Area_panA);
    getGmat_(Gmat_panB, xcoords_panB, Area_panB);
    getGmat_(Gmat_rib, xcoords_rib, Area_rib);

    // Get Ginv and the Hr, Hc and Hh matrices: Panel A.

    Ginv_panA = Gmat_panA;
    jem::numeric::LUSolver::invert(Ginv_panA, det_panA);
    Hr_panA = Ginv_panA(slice(0*nodeCount_panrib_,1*nodeCount_panrib_),slice(0,3*nodeCount_panrib_)); 
    Hc_panA = Ginv_panA(slice(1*nodeCount_panrib_,2*nodeCount_panrib_),slice(0,3*nodeCount_panrib_)); 
    Hh_panA = Ginv_panA(slice(2*nodeCount_panrib_,3*nodeCount_panrib_),slice(0,3*nodeCount_panrib_));  

    // Get Ginv and the Hr, Hc and Hh matrices: Panel B.

    Ginv_panB = Gmat_panB;
    jem::numeric::LUSolver::invert(Ginv_panB, det_panB);
    Hr_panB = Ginv_panB(slice(0*nodeCount_panrib_,1*nodeCount_panrib_),slice(0,3*nodeCount_panrib_)); 
    Hc_panB = Ginv_panB(slice(1*nodeCount_panrib_,2*nodeCount_panrib_),slice(0,3*nodeCount_panrib_)); 
    Hh_panB = Ginv_panB(slice(2*nodeCount_panrib_,3*nodeCount_panrib_),slice(0,3*nodeCount_panrib_));   

    // Get Ginv and the Hr, Hc and Hh matrices: Rib.

    Ginv_rib = Gmat_rib;
    jem::numeric::LUSolver::invert(Ginv_rib, det_rib);
    Hr_rib = Ginv_rib(slice(0*nodeCount_panrib_,1*nodeCount_panrib_),slice(0,3*nodeCount_panrib_)); 
    Hc_rib = Ginv_rib(slice(1*nodeCount_panrib_,2*nodeCount_panrib_),slice(0,3*nodeCount_panrib_)); 
    Hh_rib = Ginv_rib(slice(2*nodeCount_panrib_,3*nodeCount_panrib_),slice(0,3*nodeCount_panrib_));  

    // Get the integration points over the interface in the 2D Local coordinate system of the BLine2 element.

    shape_->getGlobalIntegrationPoints ( ipxcoords, xcoords );
    
    // Get the integration points over the interface in the 2D Global coordinate system of the BLine2 interface.

    shape_->getGlobalIntegrationPoints ( ipcoords2D, coords2D );

    // Get the Weights for numerical integration.

    shape_->getIntegrationWeights ( ipWeights, xcoords );

    // Get the integration points over the interface in the 3D Local coordinate systems of the Panel A, Panel B and Ribs elements.

    get3DLocalcoordinates_(ipxcoords_panA, transMat_pan, ipxcoords, transMat_BL2, coords_panA, coords);
    get3DLocalcoordinates_(ipxcoords_panB, transMat_pan, ipxcoords, transMat_BL2, coords_panB, coords);
    get3DLocalcoordinates_(ipxcoords_rib, transMat_rib, ipxcoords, transMat_BL2, coords_rib, coords);

    // Correction of the y coordinate of the ribs integration points

    for ( idx_t j = 0; j < ipxcoords_rib.size(1); j++ )
    {
      ipxcoords_rib(1,j) = ipxcoords_rib(1,j) + h_pan/2.0;
    }

    // Get the spatial derivatives of the membrane shape functions in the Local element system.

    panribshape_->getShapeGradients ( grads_panA, ipWeights_panrib, xcoords_panA );
    panribshape_->getShapeGradients ( grads_panB, ipWeights_panrib, xcoords_panB );
    panribshape_->getShapeGradients ( grads_rib, ipWeights_panrib, xcoords_rib ); 

    // Get the BA, Bx and By matrices, which are equal for both panel and ribs.

    getBABxBymats_(BAmat, Bxmat, Bymat);

    // Get the MA, Ma and MAinv matrices for the panel A, panel B and ribs.

    getMAMamats_(MAmat_panA, Mamat_panA, xcoords_panA);
    invertMAmat_(MAinv_panA, MAmat_panA);
    getMAMamats_(MAmat_panB, Mamat_panB, xcoords_panB);
    invertMAmat_(MAinv_panB, MAmat_panB);
    getMAMamats_(MAmat_rib, Mamat_rib, xcoords_rib);
    invertMAmat_(MAinv_rib, MAmat_rib);

    // Get the Tmat matrix.

    getTmat_(Tmat_panA, xcoords_panA);
    getTmat_(Tmat_panB, xcoords_panB);
    getTmat_(Tmat_rib, xcoords_rib);

    // Get the displacements at the panel and ribs nodes in the Global coordinate system.

    elemDisp0[slice(0,(dofCount/3))] = pan_disp[idofs[slice(0,(dofCount/3))]];
    elemDisp0[slice((dofCount/3),(2*dofCount/3))] = pan_disp[idofs[slice((dofCount/3),(2*dofCount/3))]];
    elemDisp0[slice((2*dofCount/3),dofCount)] = rib_disp[idofs[slice((2*dofCount/3),dofCount)]];

    // Transform the displacements of the panel and ribs nodes to the Local coordinate system.

    elemxDisp0 = mc1.matmul ( transMatDof, elemDisp0 );

    // Get the panel A, panel B and ribs membrane displacements at the element nodes.

    getMembraneDisp_(Disp_mem_panA, Disp_mem_panB, Disp_mem_rib, elemxDisp0);

    // Reordering the Dofs to match the order in the shell-shell cohesive line element.

    elemxDisp = 0.0;
    for ( idx_t i = 0; i < dofCount; i++ )
    {
      elemxDisp[i] = elemxDisp0[ConnectInv[i]];
    }

    // Assemble the element matrix and force vector in the local coordinate system.

    elemxMat   = 0.0;
    elemxForce = 0.0;

    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    {
      // Get the C matrix for the panel and ribs.

      for ( idx_t ip_panrib = 0; ip_panrib < ipCount_panrib_; ip_panrib++ )
      {
        // Get the plate stiffness matrix for the panel A, panel B and ribs.
        
        getShapeGrads_ ( b_mem_panA, grads_panA(ALL,ALL,ip_panrib) );
        matmul ( strain_mem_panA, b_mem_panA, Disp_mem_panA );
        material_[0]->update ( stress_mem_panA, stiff_panA, strain_mem_panA, ipoint_panA++ );
        D_plate_panA = stiff_panA*(pow(thickness_[0], 3))/12.0;

        getShapeGrads_ ( b_mem_panB, grads_panB(ALL,ALL,ip_panrib) );
        matmul ( strain_mem_panB, b_mem_panB, Disp_mem_panB );
        material_[0]->update ( stress_mem_panB, stiff_panB, strain_mem_panB, ipoint_panB++ );
        D_plate_panB = stiff_panB*(pow(thickness_[0], 3))/12.0;

        getShapeGrads_ ( b_mem_rib, grads_rib(ALL,ALL,ip_panrib) );
        matmul ( strain_mem_rib, b_mem_rib, Disp_mem_rib );
        material_[1]->update ( stress_mem_rib, stiff_rib, strain_mem_rib, ipoint_rib++ );
        D_plate_rib = stiff_rib*(pow(thickness_[1], 3))/12.0;

        // Get the H matrix for the panel A, panel B and ribs.

        getHmatAnalytic_(Hmat_panA, D_plate_panA, xcoords_panA);
        Hinv_panA = Hmat_panA;
        jem::numeric::Cholesky::invert(Hinv_panA);

        getHmatAnalytic_(Hmat_panB, D_plate_panB, xcoords_panB);
        Hinv_panB = Hmat_panB;
        jem::numeric::Cholesky::invert(Hinv_panB);

        getHmatAnalytic_(Hmat_rib, D_plate_rib, xcoords_rib);
        Hinv_rib = Hmat_rib;
        jem::numeric::Cholesky::invert(Hinv_rib);

        // Get the B and C matrices for the panel A, panel B and ribs.

        getBmat_(Bmat_panA, D_plate_panA, xcoords_panA);
        Cmat_panA = mc3.matmul ( Hinv_panA, Bmat_panA, Tmat_panA );

        getBmat_(Bmat_panB, D_plate_panB, xcoords_panB);
        Cmat_panB = mc3.matmul ( Hinv_panB, Bmat_panB, Tmat_panB );

        getBmat_(Bmat_rib, D_plate_rib, xcoords_rib);
        Cmat_rib = mc3.matmul ( Hinv_rib, Bmat_rib, Tmat_rib );
      }

      double x, y, z_rib;

      // Get the integration points in the panel A coordinate system.

      x = ipxcoords_panA(0,ip);
      y = ipxcoords_panA(1,ip);

      // Get the Nr, Nc and Nh matrices for the panel A.

      getNrNcNhmats_(Nr_panA, Nc_panA, Nh_panA, xcoords_panA, Area_panA, x, y);

      // Get the Nu and Nv matrices for the panel A.

      getNuNvmats_(Nu_panA, Nv_panA, Nr_panA, Hr_panA, Nc_panA, Hc_panA, Nh_panA, Hh_panA);

      // Get the Nrx, Ncx and Nhx matrices for the panel A.

      getNrxNcxNhxmats_(Nrx_panA, Ncx_panA, Nhx_panA, xcoords_panA, Area_panA, x, y);

      // Get the Nry, Ncy and Nhy matrices for the panel A.

      getNryNcyNhymats_(Nry_panA, Ncy_panA, Nhy_panA, xcoords_panA, Area_panA, x, y);

      // Get the Nvx and Nuy matrices for the panel A.

      getNvxmat_(Nvx_panA, Nrx_panA, Hr_panA, Ncx_panA, Hc_panA, Nhx_panA, Hh_panA);

      getNuymat_(Nuy_panA, Nry_panA, Hr_panA, Ncy_panA, Hc_panA, Nhy_panA, Hh_panA);

      // Get the S, R, Rx and Ry matrices for the panel A.

      getSRRxRymats_(Smat_panA, Rmat_panA, Rxmat_panA, Rymat_panA, x, y);

      // Get the integration points in the panel B coordinate system.

      x = ipxcoords_panB(0,ip);
      y = ipxcoords_panB(1,ip);

      // Get the Nr, Nc and Nh matrices for the panel B.

      getNrNcNhmats_(Nr_panB, Nc_panB, Nh_panB, xcoords_panB, Area_panB, x, y);

      // Get the Nu and Nv matrices for the panel B.

      getNuNvmats_(Nu_panB, Nv_panB, Nr_panB, Hr_panB, Nc_panB, Hc_panB, Nh_panB, Hh_panB);

      // Get the Nrx, Ncx and Nhx matrices for the panel B.

      getNrxNcxNhxmats_(Nrx_panB, Ncx_panB, Nhx_panB, xcoords_panB, Area_panB, x, y);

      // Get the Nry, Ncy and Nhy matrices for the panel B.

      getNryNcyNhymats_(Nry_panB, Ncy_panB, Nhy_panB, xcoords_panB, Area_panB, x, y);

      // Get the Nvx and Nuy matrices for the panel B.

      getNvxmat_(Nvx_panB, Nrx_panB, Hr_panB, Ncx_panB, Hc_panB, Nhx_panB, Hh_panB);

      getNuymat_(Nuy_panB, Nry_panB, Hr_panB, Ncy_panB, Hc_panB, Nhy_panB, Hh_panB);

      // Get the S, R, Rx and Ry matrices for the panel B.

      getSRRxRymats_(Smat_panB, Rmat_panB, Rxmat_panB, Rymat_panB, x, y);

      // Get the integration points in the ribs coordinate system.

      x = ipxcoords_rib(0,ip);
      y = ipxcoords_rib(1,ip); 

      // Get the Nr, Nc and Nh matrices for the rib.

      getNrNcNhmats_(Nr_rib, Nc_rib, Nh_rib, xcoords_rib, Area_rib, x, y);

      // Get the Nu and Nv matrices for the rib.

      getNuNvmats_(Nu_rib, Nv_rib, Nr_rib, Hr_rib, Nc_rib, Hc_rib, Nh_rib, Hh_rib);

      // Get the S, R, Rx and Ry matrices for the ribs.

      getSRRxRymats_(Smat_rib, Rmat_rib, Rxmat_rib, Rymat_rib, x, y);

      // Get the Nw, Nwx and Nwy matrices for the panel A, panel B and ribs.
        
      BA_MaC_panA = BAmat - matmul(Mamat_panA,Cmat_panA);      
      Nw_panA = matmul ( matmul(Smat_panA.transpose(),MAinv_panA), BA_MaC_panA ) + matmul(Rmat_panA.transpose(),Cmat_panA);
      Nwx_panA = matmul ( matmul(Bxmat,MAinv_panA), BA_MaC_panA ) + matmul(Rxmat_panA.transpose(),Cmat_panA);
      Nwy_panA = matmul ( matmul(Bymat,MAinv_panA), BA_MaC_panA ) + matmul(Rymat_panA.transpose(),Cmat_panA);

      BA_MaC_panB = BAmat - matmul(Mamat_panB,Cmat_panB);      
      Nw_panB = matmul ( matmul(Smat_panB.transpose(),MAinv_panB), BA_MaC_panB ) + matmul(Rmat_panB.transpose(),Cmat_panB);
      Nwx_panB = matmul ( matmul(Bxmat,MAinv_panB), BA_MaC_panB ) + matmul(Rxmat_panB.transpose(),Cmat_panB);
      Nwy_panB = matmul ( matmul(Bymat,MAinv_panB), BA_MaC_panB ) + matmul(Rymat_panB.transpose(),Cmat_panB);

      BA_MaC_rib = BAmat - matmul(Mamat_rib,Cmat_rib);
      Nw_rib = matmul ( matmul(Smat_rib.transpose(),MAinv_rib), BA_MaC_rib ) + matmul(Rmat_rib.transpose(),Cmat_rib);
      Nwx_rib = matmul ( matmul(Bxmat,MAinv_rib), BA_MaC_rib ) + matmul(Rxmat_rib.transpose(),Cmat_rib);
      Nwy_rib = matmul ( matmul(Bymat,MAinv_rib), BA_MaC_rib ) + matmul(Rymat_rib.transpose(),Cmat_rib);

      // Assemble the BCD matrix of the structural cohesive element.

      BCDmat(0,slice(0*nodeCount_panrib_,3*nodeCount_panrib_)) = 0.0;
      BCDmat(0,slice(3*nodeCount_panrib_,6*nodeCount_panrib_)) = -(1.0/2.0)*Nw_panA(0,ALL);
      BCDmat(0,slice(6*nodeCount_panrib_,9*nodeCount_panrib_)) = 0.0;
      BCDmat(0,slice(9*nodeCount_panrib_,12*nodeCount_panrib_)) = -(1.0/2.0)*Nw_panB(0,ALL);
      BCDmat(0,slice(12*nodeCount_panrib_,15*nodeCount_panrib_)) = Nv_rib(0,ALL);
      BCDmat(0,slice(15*nodeCount_panrib_,18*nodeCount_panrib_)) = 0.0;

      BCDmat(1,slice(0*nodeCount_panrib_,3*nodeCount_panrib_)) = -(1.0/2.0)*Nu_panA(0,ALL);
      BCDmat(1,slice(3*nodeCount_panrib_,6*nodeCount_panrib_)) =  (h_pan/4.0)*Nwx_panA(0,ALL);
      BCDmat(1,slice(6*nodeCount_panrib_,9*nodeCount_panrib_)) = -(1.0/2.0)*Nu_panB(0,ALL);
      BCDmat(1,slice(9*nodeCount_panrib_,12*nodeCount_panrib_)) = (h_pan/4.0)*Nwx_panB(0,ALL);
      BCDmat(1,slice(12*nodeCount_panrib_,15*nodeCount_panrib_)) = Nu_rib(0,ALL);
      BCDmat(1,slice(15*nodeCount_panrib_,18*nodeCount_panrib_)) = 0.0;
      
      BCDmat(2,slice(0*nodeCount_panrib_,3*nodeCount_panrib_)) = -(1.0/2.0)*Nv_panA(0,ALL);
      BCDmat(2,slice(3*nodeCount_panrib_,6*nodeCount_panrib_)) =  (h_pan/4.0)*Nwy_panA(0,ALL);
      BCDmat(2,slice(6*nodeCount_panrib_,9*nodeCount_panrib_)) = -(1.0/2.0)*Nv_panB(0,ALL);
      BCDmat(2,slice(9*nodeCount_panrib_,12*nodeCount_panrib_)) = (h_pan/4.0)*Nwy_panB(0,ALL);
      BCDmat(2,slice(12*nodeCount_panrib_,15*nodeCount_panrib_)) = 0.0;
      BCDmat(2,slice(15*nodeCount_panrib_,18*nodeCount_panrib_)) = -Nw_rib(0,ALL);

      // Assemble the BCT matrix of the structural cohesive element. 

      BCTmat(0,slice(0*nodeCount_panrib_,3*nodeCount_panrib_)) = 0.0;
      BCTmat(0,slice(3*nodeCount_panrib_,6*nodeCount_panrib_)) = -(1.0/2.0)*Nwy_panA(0,ALL);
      BCTmat(0,slice(6*nodeCount_panrib_,9*nodeCount_panrib_)) = 0.0;
      BCTmat(0,slice(9*nodeCount_panrib_,12*nodeCount_panrib_)) = -(1.0/2.0)*Nwy_panB(0,ALL);
      BCTmat(0,slice(12*nodeCount_panrib_,15*nodeCount_panrib_)) = 0.0;
      BCTmat(0,slice(15*nodeCount_panrib_,18*nodeCount_panrib_)) = Nwy_rib(0,ALL);

      BCTmat(1,slice(0*nodeCount_panrib_,3*nodeCount_panrib_)) = (1.0/4.0)*(Nvx_panA(0,ALL)-Nuy_panA(0,ALL));
      BCTmat(1,slice(3*nodeCount_panrib_,6*nodeCount_panrib_)) = 0.0;
      BCTmat(1,slice(6*nodeCount_panrib_,9*nodeCount_panrib_)) = (1.0/4.0)*(Nvx_panB(0,ALL)-Nuy_panB(0,ALL));
      BCTmat(1,slice(9*nodeCount_panrib_,12*nodeCount_panrib_)) = 0.0;
      BCTmat(1,slice(12*nodeCount_panrib_,15*nodeCount_panrib_)) = 0.0;
      BCTmat(1,slice(15*nodeCount_panrib_,18*nodeCount_panrib_)) = Nwx_rib(0,ALL);
      
      BCTmat(2,slice(0*nodeCount_panrib_,3*nodeCount_panrib_)) = 0.0;
      BCTmat(2,slice(3*nodeCount_panrib_,6*nodeCount_panrib_)) = 0.0;
      BCTmat(2,slice(6*nodeCount_panrib_,9*nodeCount_panrib_)) = 0.0;
      BCTmat(2,slice(9*nodeCount_panrib_,12*nodeCount_panrib_)) = 0.0;
      BCTmat(2,slice(12*nodeCount_panrib_,15*nodeCount_panrib_)) = 0.0;
      BCTmat(2,slice(15*nodeCount_panrib_,18*nodeCount_panrib_)) = 0.0;

      for ( idx_t jp = 0; jp < jpCount_; jp++ )
      {
        // Get the integration point in the ribs through the thickness direction.

        z_rib = jpxcoords(0,jp);

        // Assemble the BC matrix of the structural cohesive element.

        BCmat = BCDmat + z_rib*BCTmat;

        // Compute the displacement jump (in local {n,s}-frame)

        jump = matmul ( BCmat, elemxDisp );

        // D_I[((ipCount_*jpCount_-1) - (jpCount_*ip + jp))] = jump[0];

        // Get the tangent stiffness matrix and the traction
        coheMat_->update ( traction, stiff, jump, ipoint++ );

        // T_I[((ipCount_*jpCount_-1) - (jpCount_*ip + jp))] = traction[0];

        // Compute the element force vector in the Local coordinate system.

        elemxForce += ipWeights[ip] * jpWeights[jp] * mc1.matmul ( BCmat.transpose(), traction );

        // Compute the stiffness matrix in the Local coordinate system.
        
        elemxMat += ipWeights[ip] * jpWeights[jp] * mc3.matmul ( BCmat.transpose(), stiff, BCmat );
      }
      // // Print the contents of the ipcoords2D(0,:) for debugging purposes.
      // jem::System::out() << ipcoords2D(0,((ipCount_-1)-ip)) << " ";

      // Compute the rotation jump (in local {n,s}-frame)
      // jump_T = matmul ( BCTmat, elemxDisp );
      // R_I[((ipCount_-1)-ip)] = jump_T[2];
    }

    // Print the contents of the R_I vector for debugging purposes.
    // for ( idx_t i = 0; i < R_I.size(0); i++ )
    // {
    //   jem::System::out() << R_I[i] << " ";
    // }
    //   jem::System::out() << "\n";

    // Reordering the elemForce to match the order of the DOfs in Jive.

    elemxForce0 = 0.0;
    for ( idx_t i = 0; i < dofCount; i++ )
    {
      elemxForce0[i]   = elemxForce[Connect[i]];
    }

    // Reordering the elemMat to match the order of the DOfs in Jive.

    elemxMat0 = 0.0;
    for ( idx_t i = 0; i < dofCount; i++ )
    {
      for ( idx_t j = 0; j < dofCount; j++ )
      {
        elemxMat0(i,j) = elemxMat(Connect[i],Connect[j]);
      }
    }

    // Transform the elemxForce and elemxMat to the Global coordinate system.

    elemForce0 = mc1.matmul ( transMatDof.transpose(), elemxForce0 );
    elemMat0 = mc3.matmul ( transMatDof.transpose(), elemxMat0, transMatDof );

    // Add the element matrix to the global stiffness matrix.

    if ( mbuilder != NIL )
    {
      mbuilder->addBlock ( idofs, idofs, elemMat0 );
    }

    // // Print the contents of the elemMat0 matrix for debugging purposes.
    // for ( idx_t i = 0; i < elemMat0.size(0); i++ )
    // {
    //   for ( idx_t j = 0; j < elemMat0.size(1); j++ )
    //   {
    //   jem::System::out() << elemMat0(i,j) << " ";
    //   }
    // }

    // Add the element force vector to the global force vector.

    force[idofs] += elemForce0;

    if ( jem::Float::isNaN( sum( elemForce0 ) ) ||
        jem::Float::isNaN( sum( elemMat0   ) ) )
    { 
      System::out() << "Interf something's wrong in element " << ie << endl;
    }
  }
  jem::System::out() << "\n";
}


//-----------------------------------------------------------------------
//   getFrictionForce_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::getFrictionForce_

  ( const Vector&       fsh0,
    const Vector&       pan_disp,
    const Vector&       rib_disp ) const

{
  const idx_t dofCount   = 6 * nodeCount_;

  Matrix      stiff            ( qRank_, qRank_ );
  Matrix      BCDmat           ( qRank_, dofCount );
  Matrix      BCTmat           ( qRank_, dofCount );
  Matrix      BCmat            ( qRank_, dofCount );

  Matrix      stiff_panA       ( qRank_, qRank_ );
  Matrix      stiff_panB       ( qRank_, qRank_ );
  Matrix      stiff_rib        ( qRank_, qRank_ );
  Matrix      D_plate_panA     ( qRank_, qRank_ );
  Matrix      D_plate_panB     ( qRank_, qRank_ );
  Matrix      D_plate_rib      ( qRank_, qRank_ );

  Matrix      coords           ( qRank_, nodeCount_BLine2_ );
  Matrix      coords_panA      ( qRank_, nodeCount_panrib_ );
  Matrix      coords_panB      ( qRank_, nodeCount_panrib_ );
  Matrix      coords_rib       ( qRank_, nodeCount_panrib_ );
  Matrix      xcoords          ( rank_, nodeCount_BLine2_ );
  Matrix      thick_xcoords    ( thick_rank_, nodeCount_Line2_ );
  Matrix      xcoords_panA     ( rank_, nodeCount_panrib_ );
  Matrix      xcoords_panB     ( rank_, nodeCount_panrib_ );
  Matrix      xcoords_rib      ( rank_, nodeCount_panrib_ );

  Matrix      ipxcoords        ( rank_, ipCount_ );
  Matrix      jpxcoords        ( thick_rank_, jpCount_ );
  Matrix      ipxcoords_panA   ( qRank_, ipCount_ );
  Matrix      ipxcoords_panB   ( qRank_, ipCount_ );
  Matrix      ipxcoords_rib    ( qRank_, ipCount_ );

  Matrix      b_mem_panA       ( qRank_, 2*nodeCount_panrib_ );
  Vector      Disp_mem_panA    ( 2*nodeCount_panrib_ );
  Vector      strain_mem_panA  ( qRank_ );
  Vector      stress_mem_panA  ( qRank_ );
  Matrix      Gmat_panA        ( 3*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Ginv_panA        ( 3*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hr_panA          ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hc_panA          ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hh_panA          ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Nr_panA          ( rank_, nodeCount_panrib_ );
  Matrix      Nc_panA          ( rank_, nodeCount_panrib_ );
  Matrix      Nh_panA          ( rank_, nodeCount_panrib_ );
  Matrix      Nrx_panA         ( rank_, nodeCount_panrib_ );
  Matrix      Ncx_panA         ( rank_, nodeCount_panrib_ );
  Matrix      Nhx_panA         ( rank_, nodeCount_panrib_ );
  Matrix      Nry_panA         ( rank_, nodeCount_panrib_ );
  Matrix      Ncy_panA         ( rank_, nodeCount_panrib_ );
  Matrix      Nhy_panA         ( rank_, nodeCount_panrib_ );

  Matrix      b_mem_panB       ( qRank_, 2*nodeCount_panrib_ );
  Vector      Disp_mem_panB    ( 2*nodeCount_panrib_ );
  Vector      strain_mem_panB  ( qRank_ );
  Vector      stress_mem_panB  ( qRank_ );
  Matrix      Gmat_panB        ( 3*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Ginv_panB        ( 3*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hr_panB          ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hc_panB          ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hh_panB          ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Nr_panB          ( rank_, nodeCount_panrib_ );
  Matrix      Nc_panB          ( rank_, nodeCount_panrib_ );
  Matrix      Nh_panB          ( rank_, nodeCount_panrib_ );
  Matrix      Nrx_panB         ( rank_, nodeCount_panrib_ );
  Matrix      Ncx_panB         ( rank_, nodeCount_panrib_ );
  Matrix      Nhx_panB         ( rank_, nodeCount_panrib_ );
  Matrix      Nry_panB         ( rank_, nodeCount_panrib_ );
  Matrix      Ncy_panB         ( rank_, nodeCount_panrib_ );
  Matrix      Nhy_panB         ( rank_, nodeCount_panrib_ );

  Matrix      b_mem_rib        ( qRank_, 2*nodeCount_panrib_ );
  Vector      Disp_mem_rib     ( 2*nodeCount_panrib_ );
  Vector      strain_mem_rib   ( qRank_ );
  Vector      stress_mem_rib   ( qRank_ );
  Matrix      Gmat_rib         ( 3*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Ginv_rib         ( 3*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hr_rib           ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hc_rib           ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hh_rib           ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Nr_rib           ( rank_, nodeCount_panrib_ );
  Matrix      Nc_rib           ( rank_, nodeCount_panrib_ );
  Matrix      Nh_rib           ( rank_, nodeCount_panrib_ );

  Matrix      Bmat_panA        ( 7, 4*nodeCount_panrib_);
  Matrix      Tmat_panA        ( 4*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hmat_panA        ( 7, 7);
  Matrix      Hinv_panA        ( 7, 7);
  Matrix      Bmat_panB        ( 7, 4*nodeCount_panrib_);
  Matrix      Tmat_panB        ( 4*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hmat_panB        ( 7, 7);
  Matrix      Hinv_panB        ( 7, 7);
  Matrix      Bmat_rib         ( 7, 4*nodeCount_panrib_);
  Matrix      Tmat_rib         ( 4*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hmat_rib         ( 7, 7);
  Matrix      Hinv_rib         ( 7, 7);

  Matrix      BAmat            ( nodeCount_panrib_, 3 * nodeCount_panrib_ );

  Matrix      MAmat_panA       ( nodeCount_panrib_, rank_+1 );
  Matrix      MAinv_panA       ( nodeCount_panrib_, rank_+1 );
  Matrix      Mamat_panA       ( nodeCount_panrib_, 7 );
  Matrix      Cmat_panA        ( 7, 3 * nodeCount_panrib_ );
  Matrix      BA_MaC_panA      ( nodeCount_panrib_, 3 * nodeCount_panrib_ );

  Matrix      MAmat_panB       ( nodeCount_panrib_, rank_+1 );
  Matrix      MAinv_panB       ( nodeCount_panrib_, rank_+1 );
  Matrix      Mamat_panB       ( nodeCount_panrib_, 7 );
  Matrix      Cmat_panB        ( 7, 3 * nodeCount_panrib_ );
  Matrix      BA_MaC_panB      ( nodeCount_panrib_, 3 * nodeCount_panrib_ );

  Matrix      MAmat_rib        ( nodeCount_panrib_, rank_+1 );
  Matrix      MAinv_rib        ( nodeCount_panrib_, rank_+1 );
  Matrix      Mamat_rib        ( nodeCount_panrib_, 7 );
  Matrix      Cmat_rib         ( 7, 3 * nodeCount_panrib_ );
  Matrix      BA_MaC_rib       ( nodeCount_panrib_, 3 * nodeCount_panrib_ );

  Matrix      Bxmat            ( 1, rank_+1 );
  Matrix      Bymat            ( 1, rank_+1 );
  Matrix      Smat_panA        ( rank_+1, 1 );
  Matrix      Rmat_panA        ( 7, 1 );
  Matrix      Rxmat_panA       ( 7, 1 );
  Matrix      Rymat_panA       ( 7, 1 );
  Matrix      Smat_panB        ( rank_+1, 1 );
  Matrix      Rmat_panB        ( 7, 1 );
  Matrix      Rxmat_panB       ( 7, 1 );
  Matrix      Rymat_panB       ( 7, 1 );
  Matrix      Smat_rib         ( rank_+1, 1 );
  Matrix      Rmat_rib         ( 7, 1 );
  Matrix      Rxmat_rib        ( 7, 1 );
  Matrix      Rymat_rib        ( 7, 1 );

  Matrix      Nu_panA          ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nv_panA          ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nuy_panA         ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nvx_panA         ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nw_panA          ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nwx_panA         ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nwy_panA         ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nu_panB          ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nv_panB          ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nuy_panB         ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nvx_panB         ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nw_panB          ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nwx_panB         ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nwy_panB         ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nu_rib           ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nv_rib           ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nw_rib           ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nwx_rib          ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nwy_rib          ( 1, 3 * nodeCount_panrib_ );

  Vector      elemForce    ( dofCount  );
  Vector      elemxForce   ( dofCount  );
  Vector      elemxDisp    ( dofCount  );

  IdxVector   Connect      ( dofCount);
  IdxVector   ConnectInv   ( dofCount);

  Vector      elemForce0   ( dofCount  );
  Vector      elemxForce0  ( dofCount  );
  Vector      elemDisp0    ( dofCount  );
  Vector      elemxDisp0   ( dofCount  );

  Vector      jump         ( qRank_     );
  Vector      traction     ( qRank_     );
  Vector      trac         ( qRank_     );

  Matrix      transMat_pan ( qRank_, qRank_ );
  Matrix      transMat_rib ( qRank_, qRank_ );
  Matrix      transMat_BL2 ( qRank_, qRank_ );
  Matrix      transMatDof  ( dofCount, dofCount );

  IdxVector   inodes       ( nodeCount_BLine2_ );
  IdxVector   inodes_panA  ( nodeCount_panrib_ );
  IdxVector   inodes_panB  ( nodeCount_panrib_ );
  IdxVector   inodes_rib   ( nodeCount_panrib_ );
  IdxVector   idofs        ( dofCount );
  Vector      ipWeights    ( ipCount_ );
  Vector      jpWeights    ( jpCount_ );

  Cubix       grads_panA       ( rank_, nodeCount_panrib_, ipCount_panrib_ );
  Cubix       grads_panB       ( rank_, nodeCount_panrib_, ipCount_panrib_ );
  Cubix       grads_rib        ( rank_, nodeCount_panrib_, ipCount_panrib_ );
  Vector      ipWeights_panrib ( ipCount_panrib_ );

  double Area_panA, Area_panB, Area_rib;
  double det_panA, det_panB, det_rib;

  MChain1     mc1;
  MChain2     mc2;
  MChain3     mc3;

  idx_t       ipoint = 0;
  idx_t       ipoint_panA = 0;
  idx_t       ipoint_panB = 0;
  idx_t       ipoint_rib = 0;

  // Get panel and rib thickness.

  double h_pan = thickness_[0];
  double h_rib = thickness_[1];

  // Get the thick_xcoords for the through the thickness numerical integration.

  getThick_xcoords_(thick_xcoords, h_rib);

  // Get the integration points for the through the thickness numerical integration.

  thick_shape_->getGlobalIntegrationPoints ( jpxcoords, thick_xcoords );

  // Get the Weights for the through the thickness numerical integration.

  thick_shape_->getIntegrationWeights ( jpWeights, thick_xcoords );

  // Dofs connectivity from the Jive order to the shell-shell cohesive line element order.

  getConnectivity_ ( Connect, ConnectInv );

  // Iterate over all elements assigned to this model.

  for ( idx_t ie = 0; ie < elemCount_; ie++ )
  {
    idx_t  ielem = ielems_[ie];
    idx_t  ielem_panSideA;
    idx_t  ielem_panSideB;
    idx_t  ielem_rib;

    // Get the BLine2 element coordinates.

    elems_.getElemNodes  ( inodes, ielem );
    nodes_.getSomeCoords ( coords, inodes );

    // Get the id of the panel and ribs elements connected to the BLine2 element.

    getPanRibElements_(ielem_panSideA, ielem_panSideB, ielem_rib, pan_elems_, rib_elems_, coords);

    // Get the panel and ribs elements coordinates and the DOFs.

    pan_elems_.getElemNodes  ( inodes_panA, ielem_panSideA );
    pan_nodes_.getSomeCoords ( coords_panA, inodes_panA );

    pan_elems_.getElemNodes  ( inodes_panB, ielem_panSideB );
    pan_nodes_.getSomeCoords ( coords_panB, inodes_panB );

    rib_elems_.getElemNodes  ( inodes_rib, ielem_rib );
    rib_nodes_.getSomeCoords ( coords_rib, inodes_rib );

    pan_dofs_->getDofIndices ( idofs[slice(0,(dofCount/3))], inodes_panA, pan_dofTypes_ );
    pan_dofs_->getDofIndices ( idofs[slice((dofCount/3),(2*dofCount/3))], inodes_panB, pan_dofTypes_ );
    rib_dofs_->getDofIndices ( idofs[slice((2*dofCount/3),dofCount)], inodes_rib, rib_dofTypes_ );

    // Get the orientation vectors: v_panel and v_ribs.

    Vector v_pan ( qRank_ );
    Vector v_rib ( qRank_ );
    for ( idx_t i = 0; i < qRank_; i++ )
    {
      v_pan[i] = rib_orientVec_[i];
      v_rib[i] = rib_orientVec_[i];
    }

    // Get the transformation matrices (change of basis) for the structural CE local systems: panel, ribs, BLine2 and interface.

    getTransMatrix_ ( transMat_pan, coords_panA, v_pan );
    getTransMatrix_ ( transMat_rib, coords_rib, v_rib );
    getTransMatrix_ ( transMat_BL2, coords_panA, v_rib );
    getTransMatrixDof_ ( transMatDof, transMat_pan, transMat_rib );

    // Get the nodal coordinates in the Local coordinate systems: Panel A, Panel B, Ribs and BLine2 nodes.

    get2DLocalcoordinates_(xcoords_panA, transMat_pan, coords_panA);
    get2DLocalcoordinates_(xcoords_panB, transMat_pan, coords_panB);
    get2DLocalcoordinates_(xcoords_rib, transMat_rib, coords_rib);
    get2DLocalcoordinates_(xcoords, transMat_BL2, coords);

    // Get the areas of the panel A, Panel B and Ribs triangles.

    getArea_(Area_panA, xcoords_panA);
    getArea_(Area_panB, xcoords_panB);
    getArea_(Area_rib, xcoords_rib);

    // Get the G matrices relating generalized Dofs to nodal Dofs: Gmat_panA, Gmat_panB, Gmat_ribs.

    getGmat_(Gmat_panA, xcoords_panA, Area_panA);
    getGmat_(Gmat_panB, xcoords_panB, Area_panB);
    getGmat_(Gmat_rib, xcoords_rib, Area_rib);

    // Get Ginv and the Hr, Hc and Hh matrices: Panel A.

    Ginv_panA = Gmat_panA;
    jem::numeric::LUSolver::invert(Ginv_panA, det_panA);
    Hr_panA = Ginv_panA(slice(0*nodeCount_panrib_,1*nodeCount_panrib_),slice(0,3*nodeCount_panrib_)); 
    Hc_panA = Ginv_panA(slice(1*nodeCount_panrib_,2*nodeCount_panrib_),slice(0,3*nodeCount_panrib_)); 
    Hh_panA = Ginv_panA(slice(2*nodeCount_panrib_,3*nodeCount_panrib_),slice(0,3*nodeCount_panrib_));  

    // Get Ginv and the Hr, Hc and Hh matrices: Panel B.

    Ginv_panB = Gmat_panB;
    jem::numeric::LUSolver::invert(Ginv_panB, det_panB);
    Hr_panB = Ginv_panB(slice(0*nodeCount_panrib_,1*nodeCount_panrib_),slice(0,3*nodeCount_panrib_)); 
    Hc_panB = Ginv_panB(slice(1*nodeCount_panrib_,2*nodeCount_panrib_),slice(0,3*nodeCount_panrib_)); 
    Hh_panB = Ginv_panB(slice(2*nodeCount_panrib_,3*nodeCount_panrib_),slice(0,3*nodeCount_panrib_));   

    // Get Ginv and the Hr, Hc and Hh matrices: Rib.

    Ginv_rib = Gmat_rib;
    jem::numeric::LUSolver::invert(Ginv_rib, det_rib);
    Hr_rib = Ginv_rib(slice(0*nodeCount_panrib_,1*nodeCount_panrib_),slice(0,3*nodeCount_panrib_)); 
    Hc_rib = Ginv_rib(slice(1*nodeCount_panrib_,2*nodeCount_panrib_),slice(0,3*nodeCount_panrib_)); 
    Hh_rib = Ginv_rib(slice(2*nodeCount_panrib_,3*nodeCount_panrib_),slice(0,3*nodeCount_panrib_));  

    // Get the integration points over the interface in the 2D Local coordinate system of the BLine2 element.

    shape_->getGlobalIntegrationPoints ( ipxcoords, xcoords );

    // Get the Weights for numerical integration.

    shape_->getIntegrationWeights ( ipWeights, xcoords );

    // Get the integration points over the interface in the 3D Local coordinate systems of the Panel A, Panel B and Ribs elements.

    get3DLocalcoordinates_(ipxcoords_panA, transMat_pan, ipxcoords, transMat_BL2, coords_panA, coords);
    get3DLocalcoordinates_(ipxcoords_panB, transMat_pan, ipxcoords, transMat_BL2, coords_panB, coords);
    get3DLocalcoordinates_(ipxcoords_rib, transMat_rib, ipxcoords, transMat_BL2, coords_rib, coords);

    // Correction of the y coordinate of the ribs integration points

    for ( idx_t j = 0; j < ipxcoords_rib.size(1); j++ )
    {
      ipxcoords_rib(1,j) = ipxcoords_rib(1,j) + h_pan/2.0;
    }

    // Get the spatial derivatives of the membrane shape functions in the Local element system.

    panribshape_->getShapeGradients ( grads_panA, ipWeights_panrib, xcoords_panA );
    panribshape_->getShapeGradients ( grads_panB, ipWeights_panrib, xcoords_panB );
    panribshape_->getShapeGradients ( grads_rib, ipWeights_panrib, xcoords_rib ); 

    // Get the BA, Bx and By matrices, which are equal for both panel and ribs.

    getBABxBymats_(BAmat, Bxmat, Bymat);

    // Get the MA, Ma and MAinv matrices for the panel A, panel B and ribs.

    getMAMamats_(MAmat_panA, Mamat_panA, xcoords_panA);
    invertMAmat_(MAinv_panA, MAmat_panA);
    getMAMamats_(MAmat_panB, Mamat_panB, xcoords_panB);
    invertMAmat_(MAinv_panB, MAmat_panB);
    getMAMamats_(MAmat_rib, Mamat_rib, xcoords_rib);
    invertMAmat_(MAinv_rib, MAmat_rib);

    // Get the Tmat matrix.

    getTmat_(Tmat_panA, xcoords_panA);
    getTmat_(Tmat_panB, xcoords_panB);
    getTmat_(Tmat_rib, xcoords_rib);

    // Get the displacements at the panel and ribs nodes in the Global coordinate system.

    elemDisp0[slice(0,(dofCount/3))] = pan_disp[idofs[slice(0,(dofCount/3))]];
    elemDisp0[slice((dofCount/3),(2*dofCount/3))] = pan_disp[idofs[slice((dofCount/3),(2*dofCount/3))]];
    elemDisp0[slice((2*dofCount/3),dofCount)] = rib_disp[idofs[slice((2*dofCount/3),dofCount)]];

    // Transform the displacements of the panel and ribs nodes to the Local coordinate system.

    elemxDisp0 = mc1.matmul ( transMatDof, elemDisp0 );

    // Get the panel A, panel B and ribs membrane displacements at the element nodes.

    getMembraneDisp_(Disp_mem_panA, Disp_mem_panB, Disp_mem_rib, elemxDisp0);

    // Reordering the Dofs to match the order in the shell-shell cohesive line element.

    elemxDisp = 0.0;
    for ( idx_t i = 0; i < dofCount; i++ )
    {
      elemxDisp[i] = elemxDisp0[ConnectInv[i]];
    }

    // Assemble the element force vector.

    elemxForce = 0.0;

    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    {
      // Get the C matrix for the panel and ribs.

      for ( idx_t ip_panrib = 0; ip_panrib < ipCount_panrib_; ip_panrib++ )
      {
        // Get the plate stiffness matrix for the panel A, panel B and ribs.
        
        getShapeGrads_ ( b_mem_panA, grads_panA(ALL,ALL,ip_panrib) );
        matmul ( strain_mem_panA, b_mem_panA, Disp_mem_panA );
        material_[0]->update ( stress_mem_panA, stiff_panA, strain_mem_panA, ipoint_panA++ );
        D_plate_panA = stiff_panA*(pow(thickness_[0], 3))/12.0;

        getShapeGrads_ ( b_mem_panB, grads_panB(ALL,ALL,ip_panrib) );
        matmul ( strain_mem_panB, b_mem_panB, Disp_mem_panB );
        material_[0]->update ( stress_mem_panB, stiff_panB, strain_mem_panB, ipoint_panB++ );
        D_plate_panB = stiff_panB*(pow(thickness_[0], 3))/12.0;

        getShapeGrads_ ( b_mem_rib, grads_rib(ALL,ALL,ip_panrib) );
        matmul ( strain_mem_rib, b_mem_rib, Disp_mem_rib );
        material_[1]->update ( stress_mem_rib, stiff_rib, strain_mem_rib, ipoint_rib++ );
        D_plate_rib = stiff_rib*(pow(thickness_[1], 3))/12.0;

        // Get the H matrix for the panel A, panel B and ribs.

        getHmatAnalytic_(Hmat_panA, D_plate_panA, xcoords_panA);
        Hinv_panA = Hmat_panA;
        jem::numeric::Cholesky::invert(Hinv_panA);

        getHmatAnalytic_(Hmat_panB, D_plate_panB, xcoords_panB);
        Hinv_panB = Hmat_panB;
        jem::numeric::Cholesky::invert(Hinv_panB);

        getHmatAnalytic_(Hmat_rib, D_plate_rib, xcoords_rib);
        Hinv_rib = Hmat_rib;
        jem::numeric::Cholesky::invert(Hinv_rib);

        // Get the B and C matrices for the panel A, panel B and ribs.

        getBmat_(Bmat_panA, D_plate_panA, xcoords_panA);
        Cmat_panA = mc3.matmul ( Hinv_panA, Bmat_panA, Tmat_panA );

        getBmat_(Bmat_panB, D_plate_panB, xcoords_panB);
        Cmat_panB = mc3.matmul ( Hinv_panB, Bmat_panB, Tmat_panB );

        getBmat_(Bmat_rib, D_plate_rib, xcoords_rib);
        Cmat_rib = mc3.matmul ( Hinv_rib, Bmat_rib, Tmat_rib );
      }

      double x, y, z_rib;

      // Get the integration points in the panel A coordinate system.

      x = ipxcoords_panA(0,ip);
      y = ipxcoords_panA(1,ip);

      // Get the Nr, Nc and Nh matrices for the panel A.

      getNrNcNhmats_(Nr_panA, Nc_panA, Nh_panA, xcoords_panA, Area_panA, x, y);

      // Get the Nu and Nv matrices for the panel A.

      getNuNvmats_(Nu_panA, Nv_panA, Nr_panA, Hr_panA, Nc_panA, Hc_panA, Nh_panA, Hh_panA);

      // Get the Nrx, Ncx and Nhx matrices for the panel A.

      getNrxNcxNhxmats_(Nrx_panA, Ncx_panA, Nhx_panA, xcoords_panA, Area_panA, x, y);

      // Get the Nry, Ncy and Nhy matrices for the panel A.

      getNryNcyNhymats_(Nry_panA, Ncy_panA, Nhy_panA, xcoords_panA, Area_panA, x, y);

      // Get the Nvx and Nuy matrices for the panel A.

      getNvxmat_(Nvx_panA, Nrx_panA, Hr_panA, Ncx_panA, Hc_panA, Nhx_panA, Hh_panA);

      getNuymat_(Nuy_panA, Nry_panA, Hr_panA, Ncy_panA, Hc_panA, Nhy_panA, Hh_panA);

      // Get the S, R, Rx and Ry matrices for the panel A.

      getSRRxRymats_(Smat_panA, Rmat_panA, Rxmat_panA, Rymat_panA, x, y);

      // Get the integration points in the panel B coordinate system.

      x = ipxcoords_panB(0,ip);
      y = ipxcoords_panB(1,ip);

      // Get the Nr, Nc and Nh matrices for the panel B.

      getNrNcNhmats_(Nr_panB, Nc_panB, Nh_panB, xcoords_panB, Area_panB, x, y);

      // Get the Nu and Nv matrices for the panel B.

      getNuNvmats_(Nu_panB, Nv_panB, Nr_panB, Hr_panB, Nc_panB, Hc_panB, Nh_panB, Hh_panB);

      // Get the Nrx, Ncx and Nhx matrices for the panel B.

      getNrxNcxNhxmats_(Nrx_panB, Ncx_panB, Nhx_panB, xcoords_panB, Area_panB, x, y);

      // Get the Nry, Ncy and Nhy matrices for the panel B.

      getNryNcyNhymats_(Nry_panB, Ncy_panB, Nhy_panB, xcoords_panB, Area_panB, x, y);

      // Get the Nvx and Nuy matrices for the panel B.

      getNvxmat_(Nvx_panB, Nrx_panB, Hr_panB, Ncx_panB, Hc_panB, Nhx_panB, Hh_panB);

      getNuymat_(Nuy_panB, Nry_panB, Hr_panB, Ncy_panB, Hc_panB, Nhy_panB, Hh_panB);

      // Get the S, R, Rx and Ry matrices for the panel B.

      getSRRxRymats_(Smat_panB, Rmat_panB, Rxmat_panB, Rymat_panB, x, y);

      // Get the integration points in the ribs coordinate system.

      x = ipxcoords_rib(0,ip);
      y = ipxcoords_rib(1,ip); 

      // Get the Nr, Nc and Nh matrices for the rib.

      getNrNcNhmats_(Nr_rib, Nc_rib, Nh_rib, xcoords_rib, Area_rib, x, y);

      // Get the Nu and Nv matrices for the rib.

      getNuNvmats_(Nu_rib, Nv_rib, Nr_rib, Hr_rib, Nc_rib, Hc_rib, Nh_rib, Hh_rib);

      // Get the S, R, Rx and Ry matrices for the ribs.

      getSRRxRymats_(Smat_rib, Rmat_rib, Rxmat_rib, Rymat_rib, x, y);

      // Get the Nw, Nwx and Nwy matrices for the panel A, panel B and ribs.
        
      BA_MaC_panA = BAmat - matmul(Mamat_panA,Cmat_panA);      
      Nw_panA = matmul ( matmul(Smat_panA.transpose(),MAinv_panA), BA_MaC_panA ) + matmul(Rmat_panA.transpose(),Cmat_panA);
      Nwx_panA = matmul ( matmul(Bxmat,MAinv_panA), BA_MaC_panA ) + matmul(Rxmat_panA.transpose(),Cmat_panA);
      Nwy_panA = matmul ( matmul(Bymat,MAinv_panA), BA_MaC_panA ) + matmul(Rymat_panA.transpose(),Cmat_panA);

      BA_MaC_panB = BAmat - matmul(Mamat_panB,Cmat_panB);      
      Nw_panB = matmul ( matmul(Smat_panB.transpose(),MAinv_panB), BA_MaC_panB ) + matmul(Rmat_panB.transpose(),Cmat_panB);
      Nwx_panB = matmul ( matmul(Bxmat,MAinv_panB), BA_MaC_panB ) + matmul(Rxmat_panB.transpose(),Cmat_panB);
      Nwy_panB = matmul ( matmul(Bymat,MAinv_panB), BA_MaC_panB ) + matmul(Rymat_panB.transpose(),Cmat_panB);

      BA_MaC_rib = BAmat - matmul(Mamat_rib,Cmat_rib);
      Nw_rib = matmul ( matmul(Smat_rib.transpose(),MAinv_rib), BA_MaC_rib ) + matmul(Rmat_rib.transpose(),Cmat_rib);
      Nwx_rib = matmul ( matmul(Bxmat,MAinv_rib), BA_MaC_rib ) + matmul(Rxmat_rib.transpose(),Cmat_rib);
      Nwy_rib = matmul ( matmul(Bymat,MAinv_rib), BA_MaC_rib ) + matmul(Rymat_rib.transpose(),Cmat_rib);

      // Assemble the BCE matrix of the structural cohesive element.

      BCDmat(0,slice(0*nodeCount_panrib_,3*nodeCount_panrib_)) = 0.0;
      BCDmat(0,slice(3*nodeCount_panrib_,6*nodeCount_panrib_)) = -(1.0/2.0)*Nw_panA(0,ALL);
      BCDmat(0,slice(6*nodeCount_panrib_,9*nodeCount_panrib_)) = 0.0;
      BCDmat(0,slice(9*nodeCount_panrib_,12*nodeCount_panrib_)) = -(1.0/2.0)*Nw_panB(0,ALL);
      BCDmat(0,slice(12*nodeCount_panrib_,15*nodeCount_panrib_)) = Nv_rib(0,ALL);
      BCDmat(0,slice(15*nodeCount_panrib_,18*nodeCount_panrib_)) = 0.0;

      BCDmat(1,slice(0*nodeCount_panrib_,3*nodeCount_panrib_)) = -(1.0/2.0)*Nu_panA(0,ALL);
      BCDmat(1,slice(3*nodeCount_panrib_,6*nodeCount_panrib_)) = (h_pan/4.0)*Nwx_panA(0,ALL);
      BCDmat(1,slice(6*nodeCount_panrib_,9*nodeCount_panrib_)) = -(1.0/2.0)*Nu_panB(0,ALL);
      BCDmat(1,slice(9*nodeCount_panrib_,12*nodeCount_panrib_)) = (h_pan/4.0)*Nwx_panB(0,ALL);
      BCDmat(1,slice(12*nodeCount_panrib_,15*nodeCount_panrib_)) = Nu_rib(0,ALL);
      BCDmat(1,slice(15*nodeCount_panrib_,18*nodeCount_panrib_)) = 0.0;
      
      BCDmat(2,slice(0*nodeCount_panrib_,3*nodeCount_panrib_)) = -(1.0/2.0)*Nv_panA(0,ALL);
      BCDmat(2,slice(3*nodeCount_panrib_,6*nodeCount_panrib_)) = (h_pan/4.0)*Nwy_panA(0,ALL);
      BCDmat(2,slice(6*nodeCount_panrib_,9*nodeCount_panrib_)) = -(1.0/2.0)*Nv_panB(0,ALL);
      BCDmat(2,slice(9*nodeCount_panrib_,12*nodeCount_panrib_)) = (h_pan/4.0)*Nwy_panB(0,ALL);
      BCDmat(2,slice(12*nodeCount_panrib_,15*nodeCount_panrib_)) = 0.0;
      BCDmat(2,slice(15*nodeCount_panrib_,18*nodeCount_panrib_)) = -Nw_rib(0,ALL);

      // Assemble the BCT matrix of the structural cohesive element.

      BCTmat(0,slice(0*nodeCount_panrib_,3*nodeCount_panrib_)) = 0.0;
      BCTmat(0,slice(3*nodeCount_panrib_,6*nodeCount_panrib_)) = -(1.0/2.0)*Nwy_panA(0,ALL);
      BCTmat(0,slice(6*nodeCount_panrib_,9*nodeCount_panrib_)) = 0.0;
      BCTmat(0,slice(9*nodeCount_panrib_,12*nodeCount_panrib_)) = -(1.0/2.0)*Nwy_panB(0,ALL);
      BCTmat(0,slice(12*nodeCount_panrib_,15*nodeCount_panrib_)) = 0.0;
      BCTmat(0,slice(15*nodeCount_panrib_,18*nodeCount_panrib_)) = Nwy_rib(0,ALL);

      BCTmat(1,slice(0*nodeCount_panrib_,3*nodeCount_panrib_)) = (1.0/4.0)*(Nvx_panA(0,ALL)-Nuy_panA(0,ALL));
      BCTmat(1,slice(3*nodeCount_panrib_,6*nodeCount_panrib_)) = 0.0;
      BCTmat(1,slice(6*nodeCount_panrib_,9*nodeCount_panrib_)) = (1.0/4.0)*(Nvx_panB(0,ALL)-Nuy_panB(0,ALL));
      BCTmat(1,slice(9*nodeCount_panrib_,12*nodeCount_panrib_)) = 0.0;
      BCTmat(1,slice(12*nodeCount_panrib_,15*nodeCount_panrib_)) = 0.0;
      BCTmat(1,slice(15*nodeCount_panrib_,18*nodeCount_panrib_)) = Nwx_rib(0,ALL);
      
      BCTmat(2,slice(0*nodeCount_panrib_,3*nodeCount_panrib_)) = 0.0;
      BCTmat(2,slice(3*nodeCount_panrib_,6*nodeCount_panrib_)) = 0.0;
      BCTmat(2,slice(6*nodeCount_panrib_,9*nodeCount_panrib_)) = 0.0;
      BCTmat(2,slice(9*nodeCount_panrib_,12*nodeCount_panrib_)) = 0.0;
      BCTmat(2,slice(12*nodeCount_panrib_,15*nodeCount_panrib_)) = 0.0;
      BCTmat(2,slice(15*nodeCount_panrib_,18*nodeCount_panrib_)) = 0.0;

      for ( idx_t jp = 0; jp < jpCount_; jp++ )
      {
        // Get the integration point in the ribs through the thickness direction.

        z_rib = jpxcoords(0,jp);

        // Assemble the BC matrix of the structural cohesive element.

        BCmat = BCDmat + z_rib*BCTmat;

        // Compute the displacement jump (in local {n,s}-frame)

        jump = matmul ( BCmat, elemxDisp );

        // Get the tangent stiffness matrix and the traction

        coheMat_->update ( traction, stiff, jump, ipoint++ );

        // Compute the element force vector in the Local coordinate system.

        elemxForce += ipWeights[ip] * jpWeights[jp] * mc1.matmul ( BCmat.transpose(), traction );
      }
    }

    // Reordering the elemForce to match the order of the DOfs in Jive.

    elemxForce0 = 0.0;
    for ( idx_t i = 0; i < dofCount; i++ )
    {
      elemxForce0[i]   = elemxForce[Connect[i]];
    }

    // Transform the elemxForce to the Global coordinate system.

    elemForce0 = mc1.matmul ( transMatDof.transpose(), elemxForce0 );

    // Add the element force vector to the global force vector.

    fsh0[idofs] += elemForce0;
  }
}


//-----------------------------------------------------------------------
//   writeOutput_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::writeOutput_

  ( const Properties&    globdat  ) const

{
  const idx_t dofCount   = 6 * nodeCount_;

  Matrix      stiff            ( qRank_, qRank_ );
  Matrix      BCDmat           ( qRank_, dofCount );
  Matrix      BCTmat           ( qRank_, dofCount );
  Matrix      BCmat            ( qRank_, dofCount );

  Matrix      stiff_panA       ( qRank_, qRank_ );
  Matrix      stiff_panB       ( qRank_, qRank_ );
  Matrix      stiff_rib        ( qRank_, qRank_ );
  Matrix      D_plate_panA     ( qRank_, qRank_ );
  Matrix      D_plate_panB     ( qRank_, qRank_ );
  Matrix      D_plate_rib      ( qRank_, qRank_ );

  Matrix      coords           ( qRank_, nodeCount_BLine2_ );
  Matrix      coords_panA      ( qRank_, nodeCount_panrib_ );
  Matrix      coords_panB      ( qRank_, nodeCount_panrib_ );
  Matrix      coords_rib       ( qRank_, nodeCount_panrib_ );
  Matrix      xcoords          ( rank_, nodeCount_BLine2_ );
  Matrix      thick_xcoords    ( thick_rank_, nodeCount_Line2_ );
  Matrix      xcoords_panA     ( rank_, nodeCount_panrib_ );
  Matrix      xcoords_panB     ( rank_, nodeCount_panrib_ );
  Matrix      xcoords_rib      ( rank_, nodeCount_panrib_ );

  Matrix      ipxcoords        ( rank_, ipCount_ );
  Matrix      jpxcoords        ( thick_rank_, jpCount_ );
  Matrix      ipxcoords_panA   ( qRank_, ipCount_ );
  Matrix      ipxcoords_panB   ( qRank_, ipCount_ );
  Matrix      ipxcoords_rib    ( qRank_, ipCount_ );

  Matrix      b_mem_panA       ( qRank_, 2*nodeCount_panrib_ );
  Vector      Disp_mem_panA    ( 2*nodeCount_panrib_ );
  Vector      strain_mem_panA  ( qRank_ );
  Vector      stress_mem_panA  ( qRank_ );
  Matrix      Gmat_panA        ( 3*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Ginv_panA        ( 3*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hr_panA          ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hc_panA          ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hh_panA          ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Nr_panA          ( rank_, nodeCount_panrib_ );
  Matrix      Nc_panA          ( rank_, nodeCount_panrib_ );
  Matrix      Nh_panA          ( rank_, nodeCount_panrib_ );
  Matrix      Nrx_panA         ( rank_, nodeCount_panrib_ );
  Matrix      Ncx_panA         ( rank_, nodeCount_panrib_ );
  Matrix      Nhx_panA         ( rank_, nodeCount_panrib_ );
  Matrix      Nry_panA         ( rank_, nodeCount_panrib_ );
  Matrix      Ncy_panA         ( rank_, nodeCount_panrib_ );
  Matrix      Nhy_panA         ( rank_, nodeCount_panrib_ );

  Matrix      b_mem_panB       ( qRank_, 2*nodeCount_panrib_ );
  Vector      Disp_mem_panB    ( 2*nodeCount_panrib_ );
  Vector      strain_mem_panB  ( qRank_ );
  Vector      stress_mem_panB  ( qRank_ );
  Matrix      Gmat_panB        ( 3*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Ginv_panB        ( 3*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hr_panB          ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hc_panB          ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hh_panB          ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Nr_panB          ( rank_, nodeCount_panrib_ );
  Matrix      Nc_panB          ( rank_, nodeCount_panrib_ );
  Matrix      Nh_panB          ( rank_, nodeCount_panrib_ );
  Matrix      Nrx_panB         ( rank_, nodeCount_panrib_ );
  Matrix      Ncx_panB         ( rank_, nodeCount_panrib_ );
  Matrix      Nhx_panB         ( rank_, nodeCount_panrib_ );
  Matrix      Nry_panB         ( rank_, nodeCount_panrib_ );
  Matrix      Ncy_panB         ( rank_, nodeCount_panrib_ );
  Matrix      Nhy_panB         ( rank_, nodeCount_panrib_ );

  Matrix      b_mem_rib        ( qRank_, 2*nodeCount_panrib_ );
  Vector      Disp_mem_rib     ( 2*nodeCount_panrib_ );
  Vector      strain_mem_rib   ( qRank_ );
  Vector      stress_mem_rib   ( qRank_ );
  Matrix      Gmat_rib         ( 3*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Ginv_rib         ( 3*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hr_rib           ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hc_rib           ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hh_rib           ( nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Nr_rib           ( rank_, nodeCount_panrib_ );
  Matrix      Nc_rib           ( rank_, nodeCount_panrib_ );
  Matrix      Nh_rib           ( rank_, nodeCount_panrib_ );

  Matrix      Bmat_panA        ( 7, 4*nodeCount_panrib_);
  Matrix      Tmat_panA        ( 4*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hmat_panA        ( 7, 7);
  Matrix      Hinv_panA        ( 7, 7);
  Matrix      Bmat_panB        ( 7, 4*nodeCount_panrib_);
  Matrix      Tmat_panB        ( 4*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hmat_panB        ( 7, 7);
  Matrix      Hinv_panB        ( 7, 7);
  Matrix      Bmat_rib         ( 7, 4*nodeCount_panrib_);
  Matrix      Tmat_rib         ( 4*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Matrix      Hmat_rib         ( 7, 7);
  Matrix      Hinv_rib         ( 7, 7);

  Matrix      BAmat            ( nodeCount_panrib_, 3 * nodeCount_panrib_ );

  Matrix      MAmat_panA       ( nodeCount_panrib_, rank_+1 );
  Matrix      MAinv_panA       ( nodeCount_panrib_, rank_+1 );
  Matrix      Mamat_panA       ( nodeCount_panrib_, 7 );
  Matrix      Cmat_panA        ( 7, 3 * nodeCount_panrib_ );
  Matrix      BA_MaC_panA      ( nodeCount_panrib_, 3 * nodeCount_panrib_ );

  Matrix      MAmat_panB       ( nodeCount_panrib_, rank_+1 );
  Matrix      MAinv_panB       ( nodeCount_panrib_, rank_+1 );
  Matrix      Mamat_panB       ( nodeCount_panrib_, 7 );
  Matrix      Cmat_panB        ( 7, 3 * nodeCount_panrib_ );
  Matrix      BA_MaC_panB      ( nodeCount_panrib_, 3 * nodeCount_panrib_ );

  Matrix      MAmat_rib        ( nodeCount_panrib_, rank_+1 );
  Matrix      MAinv_rib        ( nodeCount_panrib_, rank_+1 );
  Matrix      Mamat_rib        ( nodeCount_panrib_, 7 );
  Matrix      Cmat_rib         ( 7, 3 * nodeCount_panrib_ );
  Matrix      BA_MaC_rib       ( nodeCount_panrib_, 3 * nodeCount_panrib_ );

  Matrix      Bxmat            ( 1, rank_+1 );
  Matrix      Bymat            ( 1, rank_+1 );
  Matrix      Smat_panA        ( rank_+1, 1 );
  Matrix      Rmat_panA        ( 7, 1 );
  Matrix      Rxmat_panA       ( 7, 1 );
  Matrix      Rymat_panA       ( 7, 1 );
  Matrix      Smat_panB        ( rank_+1, 1 );
  Matrix      Rmat_panB        ( 7, 1 );
  Matrix      Rxmat_panB       ( 7, 1 );
  Matrix      Rymat_panB       ( 7, 1 );
  Matrix      Smat_rib         ( rank_+1, 1 );
  Matrix      Rmat_rib         ( 7, 1 );
  Matrix      Rxmat_rib        ( 7, 1 );
  Matrix      Rymat_rib        ( 7, 1 );

  Matrix      Nu_panA          ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nv_panA          ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nuy_panA         ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nvx_panA         ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nw_panA          ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nwx_panA         ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nwy_panA         ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nu_panB          ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nv_panB          ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nuy_panB         ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nvx_panB         ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nw_panB          ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nwx_panB         ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nwy_panB         ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nu_rib           ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nv_rib           ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nw_rib           ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nwx_rib          ( 1, 3 * nodeCount_panrib_ );
  Matrix      Nwy_rib          ( 1, 3 * nodeCount_panrib_ );

  Vector      elemxDisp    ( dofCount  );

  IdxVector   Connect      ( dofCount);
  IdxVector   ConnectInv   ( dofCount);

  Vector      elemDisp0    ( dofCount  );
  Vector      elemxDisp0   ( dofCount  );

  Vector      jump         ( qRank_     );
  Vector      traction     ( qRank_     );
  Vector      trac         ( qRank_     );

  Matrix      transMat_pan ( qRank_, qRank_ );
  Matrix      transMat_rib ( qRank_, qRank_ );
  Matrix      transMat_BL2 ( qRank_, qRank_ );
  Matrix      transMatDof  ( dofCount, dofCount );

  IdxVector   inodes       ( nodeCount_BLine2_ );
  IdxVector   inodes_panA  ( nodeCount_panrib_ );
  IdxVector   inodes_panB  ( nodeCount_panrib_ );
  IdxVector   inodes_rib   ( nodeCount_panrib_ );
  IdxVector   idofs        ( dofCount );
  Vector      ipWeights    ( ipCount_   );
  Vector      jpWeights    ( jpCount_ );

  Cubix       grads_panA       ( rank_, nodeCount_panrib_, ipCount_panrib_ );
  Cubix       grads_panB       ( rank_, nodeCount_panrib_, ipCount_panrib_ );
  Cubix       grads_rib        ( rank_, nodeCount_panrib_, ipCount_panrib_ );
  Vector      ipWeights_panrib ( ipCount_panrib_ );

  double Area_panA, Area_panB, Area_rib;
  double det_panA, det_panB, det_rib;

  MChain1     mc1;
  MChain2     mc2;
  MChain3     mc3;

  Vector      pan_disp;
  Vector      rib_disp;

  StateVector::get ( pan_disp, pan_dofs_, globdat );
  StateVector::get ( rib_disp, rib_dofs_, globdat );

  idx_t       ipoint = 0;
  idx_t       ipoint_panA = 0;
  idx_t       ipoint_panB = 0;
  idx_t       ipoint_rib = 0;

  // Get panel and rib thickness.

  double h_pan = thickness_[0];
  double h_rib = thickness_[1];

  // Get the thick_xcoords for the through the thickness numerical integration.

  getThick_xcoords_(thick_xcoords, h_rib);

  // Get the integration points for the through the thickness numerical integration.

  thick_shape_->getGlobalIntegrationPoints ( jpxcoords, thick_xcoords );

  // Get the Weights for the through the thickness numerical integration.

  thick_shape_->getIntegrationWeights ( jpWeights, thick_xcoords );

  // Dofs connectivity from the Jive order to the shell-shell cohesive line element order.

  getConnectivity_ ( Connect, ConnectInv );

  // Iterate over all elements assigned to this model.

  for ( idx_t ie = 0; ie < elemCount_; ie++ )
  {
    idx_t  ielem = ielems_[ie];
    idx_t  ielem_panSideA;
    idx_t  ielem_panSideB;
    idx_t  ielem_rib;

    // Get the BLine2 element coordinates.

    elems_.getElemNodes  ( inodes, ielem );
    nodes_.getSomeCoords ( coords, inodes );

    // Get the id of the panel and ribs elements connected to the BLine2 element.

    getPanRibElements_(ielem_panSideA, ielem_panSideB, ielem_rib, pan_elems_, rib_elems_, coords);

    // Get the panel and ribs elements coordinates and the DOFs.

    pan_elems_.getElemNodes  ( inodes_panA, ielem_panSideA );
    pan_nodes_.getSomeCoords ( coords_panA, inodes_panA );

    pan_elems_.getElemNodes  ( inodes_panB, ielem_panSideB );
    pan_nodes_.getSomeCoords ( coords_panB, inodes_panB );

    rib_elems_.getElemNodes  ( inodes_rib, ielem_rib );
    rib_nodes_.getSomeCoords ( coords_rib, inodes_rib );

    pan_dofs_->getDofIndices ( idofs[slice(0,(dofCount/3))], inodes_panA, pan_dofTypes_ );
    pan_dofs_->getDofIndices ( idofs[slice((dofCount/3),(2*dofCount/3))], inodes_panB, pan_dofTypes_ );
    rib_dofs_->getDofIndices ( idofs[slice((2*dofCount/3),dofCount)], inodes_rib, rib_dofTypes_ );

    // Get the orientation vectors: v_panel and v_ribs.

    Vector v_pan ( qRank_ );
    Vector v_rib ( qRank_ );
    for ( idx_t i = 0; i < qRank_; i++ )
    {
      v_pan[i] = rib_orientVec_[i];
      v_rib[i] = rib_orientVec_[i];
    }

    // Get the transformation matrices (change of basis) for the structural CE local systems: panel, ribs, BLine2 and interface.

    getTransMatrix_ ( transMat_pan, coords_panA, v_pan );
    getTransMatrix_ ( transMat_rib, coords_rib, v_rib );
    getTransMatrix_ ( transMat_BL2, coords_panA, v_rib );
    getTransMatrixDof_ ( transMatDof, transMat_pan, transMat_rib );

    // Get the nodal coordinates in the Local coordinate systems: Panel A, Panel B, Ribs and BLine2 nodes.

    get2DLocalcoordinates_(xcoords_panA, transMat_pan, coords_panA);
    get2DLocalcoordinates_(xcoords_panB, transMat_pan, coords_panB);
    get2DLocalcoordinates_(xcoords_rib, transMat_rib, coords_rib);
    get2DLocalcoordinates_(xcoords, transMat_BL2, coords);

    // Get the areas of the panel A, Panel B and Ribs triangles.

    getArea_(Area_panA, xcoords_panA);
    getArea_(Area_panB, xcoords_panB);
    getArea_(Area_rib, xcoords_rib);

    // Get the G matrices relating generalized Dofs to nodal Dofs: Gmat_panA, Gmat_panB, Gmat_ribs.

    getGmat_(Gmat_panA, xcoords_panA, Area_panA);
    getGmat_(Gmat_panB, xcoords_panB, Area_panB);
    getGmat_(Gmat_rib, xcoords_rib, Area_rib);

    // Get Ginv and the Hr, Hc and Hh matrices: Panel A.

    Ginv_panA = Gmat_panA;
    jem::numeric::LUSolver::invert(Ginv_panA, det_panA);
    Hr_panA = Ginv_panA(slice(0*nodeCount_panrib_,1*nodeCount_panrib_),slice(0,3*nodeCount_panrib_)); 
    Hc_panA = Ginv_panA(slice(1*nodeCount_panrib_,2*nodeCount_panrib_),slice(0,3*nodeCount_panrib_)); 
    Hh_panA = Ginv_panA(slice(2*nodeCount_panrib_,3*nodeCount_panrib_),slice(0,3*nodeCount_panrib_));  

    // Get Ginv and the Hr, Hc and Hh matrices: Panel B.

    Ginv_panB = Gmat_panB;
    jem::numeric::LUSolver::invert(Ginv_panB, det_panB);
    Hr_panB = Ginv_panB(slice(0*nodeCount_panrib_,1*nodeCount_panrib_),slice(0,3*nodeCount_panrib_)); 
    Hc_panB = Ginv_panB(slice(1*nodeCount_panrib_,2*nodeCount_panrib_),slice(0,3*nodeCount_panrib_)); 
    Hh_panB = Ginv_panB(slice(2*nodeCount_panrib_,3*nodeCount_panrib_),slice(0,3*nodeCount_panrib_));   

    // Get Ginv and the Hr, Hc and Hh matrices: Rib.

    Ginv_rib = Gmat_rib;
    jem::numeric::LUSolver::invert(Ginv_rib, det_rib);
    Hr_rib = Ginv_rib(slice(0*nodeCount_panrib_,1*nodeCount_panrib_),slice(0,3*nodeCount_panrib_)); 
    Hc_rib = Ginv_rib(slice(1*nodeCount_panrib_,2*nodeCount_panrib_),slice(0,3*nodeCount_panrib_)); 
    Hh_rib = Ginv_rib(slice(2*nodeCount_panrib_,3*nodeCount_panrib_),slice(0,3*nodeCount_panrib_)); 

    // Get the integration points over the interface in the 2D Local coordinate system of the BLine2 element.

    shape_->getGlobalIntegrationPoints ( ipxcoords, xcoords );

    // Get the Weights for numerical integration.

    shape_->getIntegrationWeights ( ipWeights, xcoords );

    // Get the integration points over the interface in the 3D Local coordinate systems of the Panel A, Panel B and Ribs elements.

    get3DLocalcoordinates_(ipxcoords_panA, transMat_pan, ipxcoords, transMat_BL2, coords_panA, coords);
    get3DLocalcoordinates_(ipxcoords_panB, transMat_pan, ipxcoords, transMat_BL2, coords_panB, coords);
    get3DLocalcoordinates_(ipxcoords_rib, transMat_rib, ipxcoords, transMat_BL2, coords_rib, coords);

    // Correction of the y coordinate of the ribs integration points

    for ( idx_t j = 0; j < ipxcoords_rib.size(1); j++ )
    {
      ipxcoords_rib(1,j) = ipxcoords_rib(1,j) + h_pan/2.0;
    }

    // Get the spatial derivatives of the membrane shape functions in the Local element system.

    panribshape_->getShapeGradients ( grads_panA, ipWeights_panrib, xcoords_panA );
    panribshape_->getShapeGradients ( grads_panB, ipWeights_panrib, xcoords_panB );
    panribshape_->getShapeGradients ( grads_rib, ipWeights_panrib, xcoords_rib ); 

    // Get the BA, Bx and By matrices, which are equal for both panel and ribs.

    getBABxBymats_(BAmat, Bxmat, Bymat);

    // Get the MA, Ma and MAinv matrices for the panel A, panel B and ribs.

    getMAMamats_(MAmat_panA, Mamat_panA, xcoords_panA);
    invertMAmat_(MAinv_panA, MAmat_panA);
    getMAMamats_(MAmat_panB, Mamat_panB, xcoords_panB);
    invertMAmat_(MAinv_panB, MAmat_panB);
    getMAMamats_(MAmat_rib, Mamat_rib, xcoords_rib);
    invertMAmat_(MAinv_rib, MAmat_rib);

    // Get the Tmat matrix.

    getTmat_(Tmat_panA, xcoords_panA);
    getTmat_(Tmat_panB, xcoords_panB);
    getTmat_(Tmat_rib, xcoords_rib);

    // Get the displacements at the panel and ribs nodes in the Global coordinate system.

    elemDisp0[slice(0,(dofCount/3))] = pan_disp[idofs[slice(0,(dofCount/3))]];
    elemDisp0[slice((dofCount/3),(2*dofCount/3))] = pan_disp[idofs[slice((dofCount/3),(2*dofCount/3))]];
    elemDisp0[slice((2*dofCount/3),dofCount)] = rib_disp[idofs[slice((2*dofCount/3),dofCount)]];

    // Transform the displacements of the panel and ribs nodes to the Local coordinate system.

    elemxDisp0 = mc1.matmul ( transMatDof, elemDisp0 );

    // Get the panel A, panel B and ribs membrane displacements at the element nodes.

    getMembraneDisp_(Disp_mem_panA, Disp_mem_panB, Disp_mem_rib, elemxDisp0);

    // Reordering the Dofs to match the order in the shell-shell cohesive line element.

    elemxDisp = 0.0;
    for ( idx_t i = 0; i < dofCount; i++ )
    {
      elemxDisp[i] = elemxDisp0[ConnectInv[i]];
    }

    // Evaluate and write jump and traction

    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    {
      // Get the C matrix for the panel and ribs.

      for ( idx_t ip_panrib = 0; ip_panrib < ipCount_panrib_; ip_panrib++ )
      {
        // Get the plate stiffness matrix for the panel A, panel B and ribs.
        
        getShapeGrads_ ( b_mem_panA, grads_panA(ALL,ALL,ip_panrib) );
        matmul ( strain_mem_panA, b_mem_panA, Disp_mem_panA );
        material_[0]->update ( stress_mem_panA, stiff_panA, strain_mem_panA, ipoint_panA++ );
        D_plate_panA = stiff_panA*(pow(thickness_[0], 3))/12.0;

        getShapeGrads_ ( b_mem_panB, grads_panB(ALL,ALL,ip_panrib) );
        matmul ( strain_mem_panB, b_mem_panB, Disp_mem_panB );
        material_[0]->update ( stress_mem_panB, stiff_panB, strain_mem_panB, ipoint_panB++ );
        D_plate_panB = stiff_panB*(pow(thickness_[0], 3))/12.0;

        getShapeGrads_ ( b_mem_rib, grads_rib(ALL,ALL,ip_panrib) );
        matmul ( strain_mem_rib, b_mem_rib, Disp_mem_rib );
        material_[1]->update ( stress_mem_rib, stiff_rib, strain_mem_rib, ipoint_rib++ );
        D_plate_rib = stiff_rib*(pow(thickness_[1], 3))/12.0;

        // Get the H matrix for the panel A, panel B and ribs.

        getHmatAnalytic_(Hmat_panA, D_plate_panA, xcoords_panA);
        Hinv_panA = Hmat_panA;
        jem::numeric::Cholesky::invert(Hinv_panA);

        getHmatAnalytic_(Hmat_panB, D_plate_panB, xcoords_panB);
        Hinv_panB = Hmat_panB;
        jem::numeric::Cholesky::invert(Hinv_panB);

        getHmatAnalytic_(Hmat_rib, D_plate_rib, xcoords_rib);
        Hinv_rib = Hmat_rib;
        jem::numeric::Cholesky::invert(Hinv_rib);

        // Get the B and C matrices for the panel A, panel B and ribs.

        getBmat_(Bmat_panA, D_plate_panA, xcoords_panA);
        Cmat_panA = mc3.matmul ( Hinv_panA, Bmat_panA, Tmat_panA );

        getBmat_(Bmat_panB, D_plate_panB, xcoords_panB);
        Cmat_panB = mc3.matmul ( Hinv_panB, Bmat_panB, Tmat_panB );

        getBmat_(Bmat_rib, D_plate_rib, xcoords_rib);
        Cmat_rib = mc3.matmul ( Hinv_rib, Bmat_rib, Tmat_rib );
      }

      double x, y, z_rib;

      // Get the integration points in the panel A coordinate system.

      x = ipxcoords_panA(0,ip);
      y = ipxcoords_panA(1,ip);

      // Get the Nr, Nc and Nh matrices for the panel A.

      getNrNcNhmats_(Nr_panA, Nc_panA, Nh_panA, xcoords_panA, Area_panA, x, y);

      // Get the Nu and Nv matrices for the panel A.

      getNuNvmats_(Nu_panA, Nv_panA, Nr_panA, Hr_panA, Nc_panA, Hc_panA, Nh_panA, Hh_panA);

      // Get the Nrx, Ncx and Nhx matrices for the panel A.

      getNrxNcxNhxmats_(Nrx_panA, Ncx_panA, Nhx_panA, xcoords_panA, Area_panA, x, y);

      // Get the Nry, Ncy and Nhy matrices for the panel A.

      getNryNcyNhymats_(Nry_panA, Ncy_panA, Nhy_panA, xcoords_panA, Area_panA, x, y);

      // Get the Nvx and Nuy matrices for the panel A.

      getNvxmat_(Nvx_panA, Nrx_panA, Hr_panA, Ncx_panA, Hc_panA, Nhx_panA, Hh_panA);

      getNuymat_(Nuy_panA, Nry_panA, Hr_panA, Ncy_panA, Hc_panA, Nhy_panA, Hh_panA);

      // Get the S, R, Rx and Ry matrices for the panel A.

      getSRRxRymats_(Smat_panA, Rmat_panA, Rxmat_panA, Rymat_panA, x, y);

      // Get the integration points in the panel B coordinate system.

      x = ipxcoords_panB(0,ip);
      y = ipxcoords_panB(1,ip);

      // Get the Nr, Nc and Nh matrices for the panel B.

      getNrNcNhmats_(Nr_panB, Nc_panB, Nh_panB, xcoords_panB, Area_panB, x, y);

      // Get the Nu and Nv matrices for the panel B.

      getNuNvmats_(Nu_panB, Nv_panB, Nr_panB, Hr_panB, Nc_panB, Hc_panB, Nh_panB, Hh_panB);

      // Get the Nrx, Ncx and Nhx matrices for the panel B.

      getNrxNcxNhxmats_(Nrx_panB, Ncx_panB, Nhx_panB, xcoords_panB, Area_panB, x, y);

      // Get the Nry, Ncy and Nhy matrices for the panel B.

      getNryNcyNhymats_(Nry_panB, Ncy_panB, Nhy_panB, xcoords_panB, Area_panB, x, y);

      // Get the Nvx and Nuy matrices for the panel B.

      getNvxmat_(Nvx_panB, Nrx_panB, Hr_panB, Ncx_panB, Hc_panB, Nhx_panB, Hh_panB);

      getNuymat_(Nuy_panB, Nry_panB, Hr_panB, Ncy_panB, Hc_panB, Nhy_panB, Hh_panB);

      // Get the S, R, Rx and Ry matrices for the panel B.

      getSRRxRymats_(Smat_panB, Rmat_panB, Rxmat_panB, Rymat_panB, x, y);

      // Get the integration points in the ribs coordinate system.

      x = ipxcoords_rib(0,ip);
      y = ipxcoords_rib(1,ip); 

      // Get the Nr, Nc and Nh matrices for the rib.

      getNrNcNhmats_(Nr_rib, Nc_rib, Nh_rib, xcoords_rib, Area_rib, x, y);

      // Get the Nu and Nv matrices for the rib.

      getNuNvmats_(Nu_rib, Nv_rib, Nr_rib, Hr_rib, Nc_rib, Hc_rib, Nh_rib, Hh_rib);

      // Get the S, R, Rx and Ry matrices for the ribs.

      getSRRxRymats_(Smat_rib, Rmat_rib, Rxmat_rib, Rymat_rib, x, y);

      // Get the Nw, Nwx and Nwy matrices for the panel A, panel B and ribs.
        
      BA_MaC_panA = BAmat - matmul(Mamat_panA,Cmat_panA);      
      Nw_panA = matmul ( matmul(Smat_panA.transpose(),MAinv_panA), BA_MaC_panA ) + matmul(Rmat_panA.transpose(),Cmat_panA);
      Nwx_panA = matmul ( matmul(Bxmat,MAinv_panA), BA_MaC_panA ) + matmul(Rxmat_panA.transpose(),Cmat_panA);
      Nwy_panA = matmul ( matmul(Bymat,MAinv_panA), BA_MaC_panA ) + matmul(Rymat_panA.transpose(),Cmat_panA);

      BA_MaC_panB = BAmat - matmul(Mamat_panB,Cmat_panB);      
      Nw_panB = matmul ( matmul(Smat_panB.transpose(),MAinv_panB), BA_MaC_panB ) + matmul(Rmat_panB.transpose(),Cmat_panB);
      Nwx_panB = matmul ( matmul(Bxmat,MAinv_panB), BA_MaC_panB ) + matmul(Rxmat_panB.transpose(),Cmat_panB);
      Nwy_panB = matmul ( matmul(Bymat,MAinv_panB), BA_MaC_panB ) + matmul(Rymat_panB.transpose(),Cmat_panB);

      BA_MaC_rib = BAmat - matmul(Mamat_rib,Cmat_rib);
      Nw_rib = matmul ( matmul(Smat_rib.transpose(),MAinv_rib), BA_MaC_rib ) + matmul(Rmat_rib.transpose(),Cmat_rib);
      Nwx_rib = matmul ( matmul(Bxmat,MAinv_rib), BA_MaC_rib ) + matmul(Rxmat_rib.transpose(),Cmat_rib);
      Nwy_rib = matmul ( matmul(Bymat,MAinv_rib), BA_MaC_rib ) + matmul(Rymat_rib.transpose(),Cmat_rib);

      // Assemble the BCD matrix of the structural cohesive element.

      BCDmat(0,slice(0*nodeCount_panrib_,3*nodeCount_panrib_)) = 0.0;
      BCDmat(0,slice(3*nodeCount_panrib_,6*nodeCount_panrib_)) = -(1.0/2.0)*Nw_panA(0,ALL);
      BCDmat(0,slice(6*nodeCount_panrib_,9*nodeCount_panrib_)) = 0.0;
      BCDmat(0,slice(9*nodeCount_panrib_,12*nodeCount_panrib_)) = -(1.0/2.0)*Nw_panB(0,ALL);
      BCDmat(0,slice(12*nodeCount_panrib_,15*nodeCount_panrib_)) = Nv_rib(0,ALL);
      BCDmat(0,slice(15*nodeCount_panrib_,18*nodeCount_panrib_)) = 0.0;

      BCDmat(1,slice(0*nodeCount_panrib_,3*nodeCount_panrib_)) = -(1.0/2.0)*Nu_panA(0,ALL);
      BCDmat(1,slice(3*nodeCount_panrib_,6*nodeCount_panrib_)) = (h_pan/4.0)*Nwx_panA(0,ALL);
      BCDmat(1,slice(6*nodeCount_panrib_,9*nodeCount_panrib_)) = -(1.0/2.0)*Nu_panB(0,ALL);
      BCDmat(1,slice(9*nodeCount_panrib_,12*nodeCount_panrib_)) = (h_pan/4.0)*Nwx_panB(0,ALL);
      BCDmat(1,slice(12*nodeCount_panrib_,15*nodeCount_panrib_)) = Nu_rib(0,ALL);
      BCDmat(1,slice(15*nodeCount_panrib_,18*nodeCount_panrib_)) = 0.0;
      
      BCDmat(2,slice(0*nodeCount_panrib_,3*nodeCount_panrib_)) = -(1.0/2.0)*Nv_panA(0,ALL);
      BCDmat(2,slice(3*nodeCount_panrib_,6*nodeCount_panrib_)) = (h_pan/4.0)*Nwy_panA(0,ALL);
      BCDmat(2,slice(6*nodeCount_panrib_,9*nodeCount_panrib_)) = -(1.0/2.0)*Nv_panB(0,ALL);
      BCDmat(2,slice(9*nodeCount_panrib_,12*nodeCount_panrib_)) = (h_pan/4.0)*Nwy_panB(0,ALL);
      BCDmat(2,slice(12*nodeCount_panrib_,15*nodeCount_panrib_)) = 0.0;
      BCDmat(2,slice(15*nodeCount_panrib_,18*nodeCount_panrib_)) = -Nw_rib(0,ALL);
    
      // Assemble the BCT matrix of the structural cohesive element.

      BCTmat(0,slice(0*nodeCount_panrib_,3*nodeCount_panrib_)) = 0.0;
      BCTmat(0,slice(3*nodeCount_panrib_,6*nodeCount_panrib_)) = -(1.0/2.0)*Nwy_panA(0,ALL);
      BCTmat(0,slice(6*nodeCount_panrib_,9*nodeCount_panrib_)) = 0.0;
      BCTmat(0,slice(9*nodeCount_panrib_,12*nodeCount_panrib_)) = -(1.0/2.0)*Nwy_panB(0,ALL);
      BCTmat(0,slice(12*nodeCount_panrib_,15*nodeCount_panrib_)) = 0.0;
      BCTmat(0,slice(15*nodeCount_panrib_,18*nodeCount_panrib_)) = Nwy_rib(0,ALL);

      BCTmat(1,slice(0*nodeCount_panrib_,3*nodeCount_panrib_)) = (1.0/4.0)*(Nvx_panA(0,ALL)-Nuy_panA(0,ALL));
      BCTmat(1,slice(3*nodeCount_panrib_,6*nodeCount_panrib_)) = 0.0;
      BCTmat(1,slice(6*nodeCount_panrib_,9*nodeCount_panrib_)) = (1.0/4.0)*(Nvx_panB(0,ALL)-Nuy_panB(0,ALL));
      BCTmat(1,slice(9*nodeCount_panrib_,12*nodeCount_panrib_)) = 0.0;
      BCTmat(1,slice(12*nodeCount_panrib_,15*nodeCount_panrib_)) = 0.0;
      BCTmat(1,slice(15*nodeCount_panrib_,18*nodeCount_panrib_)) = Nwx_rib(0,ALL);
      
      BCTmat(2,slice(0*nodeCount_panrib_,3*nodeCount_panrib_)) = 0.0;
      BCTmat(2,slice(3*nodeCount_panrib_,6*nodeCount_panrib_)) = 0.0;
      BCTmat(2,slice(6*nodeCount_panrib_,9*nodeCount_panrib_)) = 0.0;
      BCTmat(2,slice(9*nodeCount_panrib_,12*nodeCount_panrib_)) = 0.0;
      BCTmat(2,slice(12*nodeCount_panrib_,15*nodeCount_panrib_)) = 0.0;
      BCTmat(2,slice(15*nodeCount_panrib_,18*nodeCount_panrib_)) = 0.0;

      for ( idx_t jp = 0; jp < jpCount_; jp++ )
      {
        // Get the integration point in the ribs through the thickness direction.

        z_rib = jpxcoords(0,jp);

        // Assemble the BC matrix of the structural cohesive element.

        BCmat = BCDmat + z_rib*BCTmat;

        // Compute the displacement jump (in local {n,s}-frame)

        jump = matmul ( BCmat, elemxDisp );

        // Get the tangent stiffness matrix and the traction

        coheMat_->update ( traction, stiff, jump, ipoint );

        for ( idx_t j = 0; j < rank_; j++ )
        {
          *xOut_ << jump[j] << " " << traction[j] << " ";
        }
        *xOut_ << coheMat_->giveHistory( ipoint ) 
                      << " " << norm2( traction ) 
                      << " " << coheMat_->isLoading( ipoint )
                      << '\n';
        ipoint++;                
      }
    }
  }
  xOut_->flush();
}


//-----------------------------------------------------------------------
//   getDissipation_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::getDissipation_

  ( const Properties& params ) const

{
  Matrix      coords        ( qRank_, nodeCount_BLine2_  );
  Matrix      coords_panA   ( qRank_, nodeCount_panrib_ );
  Matrix      coords_rib    ( qRank_, nodeCount_panrib_ );
  Matrix      xcoords       ( rank_, nodeCount_BLine2_ );
  Matrix      thick_xcoords ( thick_rank_, nodeCount_Line2_ );
  Matrix      transMat_BL2  ( qRank_, qRank_ );
  Matrix      ipxcoords     ( rank_, ipCount_ );
  Vector      ipWeights     ( ipCount_ );
  Matrix      jpxcoords     ( thick_rank_, jpCount_ );
  Vector      jpWeights     ( jpCount_ );
  IdxVector   inodes        ( nodeCount_BLine2_ );
  IdxVector   inodes_panA   ( nodeCount_panrib_ );
  IdxVector   inodes_rib    ( nodeCount_panrib_ );

  double dissipation = 0.;
  idx_t  ipoint = 0;

  // Get panel and rib thickness.

  double h_rib = thickness_[1];

  // Get the thick_xcoords for the through the thickness numerical integration.

  getThick_xcoords_(thick_xcoords, h_rib);

  // Get the integration points for the through the thickness numerical integration.

  thick_shape_->getGlobalIntegrationPoints ( jpxcoords, thick_xcoords );

  // Get the Weights for the through the thickness numerical integration.

  thick_shape_->getIntegrationWeights ( jpWeights, thick_xcoords );

  // loop over integration points

  for ( idx_t ie = 0; ie < elemCount_; ie++ )
  {
    idx_t  ielem = ielems_[ie];
    idx_t  ielem_panSideA;
    idx_t  ielem_panSideB;
    idx_t  ielem_rib;

    // Get the BLine2 element coordinates.

    elems_.getElemNodes  ( inodes, ielem );
    nodes_.getSomeCoords ( coords, inodes );

    // Get the id of the panel and ribs elements connected to the BLine2 element.

    getPanRibElements_(ielem_panSideA, ielem_panSideB, ielem_rib, pan_elems_, rib_elems_, coords);

    // Get the panel and ribs elements coordinates.

    pan_elems_.getElemNodes  ( inodes_panA, ielem_panSideA );
    pan_nodes_.getSomeCoords ( coords_panA, inodes_panA );

    rib_elems_.getElemNodes  ( inodes_rib, ielem_rib );
    rib_nodes_.getSomeCoords ( coords_rib, inodes_rib );

    // Get the orientation vector: v_ribs.

    Vector v_rib ( qRank_ );
    for ( idx_t i = 0; i < qRank_; i++ )
    {
      v_rib[i] = rib_orientVec_[i];
    }

    // Get the transformation matrices (change of basis) for the structural CE local systems: panel, ribs, BLine2 and interface.

    getTransMatrix_ ( transMat_BL2, coords_panA, v_rib );

    // Get the nodal coordinates in the Local coordinate systems: Panel A, Panel B, Ribs and BLine2 nodes.

    get2DLocalcoordinates_(xcoords, transMat_BL2, coords);

    // Get the Weights for numerical integration.

    shape_->getIntegrationWeights ( ipWeights, xcoords );

    // get dissipation

    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    {
      for ( idx_t jp = 0; jp < jpCount_; jp++ )
      {

        double    G  = coheMat_->giveDissipation ( ipoint++ );

        dissipation += ipWeights[ ip ] * jpWeights[ jp ] * G;
      }
    }
  }
  params.set ( myTag_, dissipation );
}

//-----------------------------------------------------------------------
//   getThick_xcoords_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::getThick_xcoords_

  ( Matrix&  thick_xcoords,
    const double&     h_rib )    const
    
{
  thick_xcoords(0,0) = -(h_rib/2.0);
  thick_xcoords(0,1) =  (h_rib/2.0);
}


//-----------------------------------------------------------------------
//   getConnectivity_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::getConnectivity_

  ( IdxVector&      Connect,
    IdxVector&      ConnectInv )   const
    
{
  // DOFs connectivity from Jive style to the U,W style.

  for ( idx_t i = 0; i < 3; i++ )
  {
    Connect[18*i+0] = 18*i+0;
    Connect[18*i+1] = 18*i+1;
    Connect[18*i+2] = 18*i+9;
    Connect[18*i+3] = 18*i+10;
    Connect[18*i+4] = 18*i+11;
    Connect[18*i+5] = 18*i+2;
    Connect[18*i+6] = 18*i+3;
    Connect[18*i+7] = 18*i+4;
    Connect[18*i+8] = 18*i+12;
    Connect[18*i+9] = 18*i+13;
    Connect[18*i+10] = 18*i+14;
    Connect[18*i+11] = 18*i+5;
    Connect[18*i+12] = 18*i+6;
    Connect[18*i+13] = 18*i+7;
    Connect[18*i+14] = 18*i+15;
    Connect[18*i+15] = 18*i+16;
    Connect[18*i+16] = 18*i+17;
    Connect[18*i+17] = 18*i+8;  
  }

  // DOFs connectivity from the U,W style to Jive style - ConnectInv.
  
  for ( idx_t i = 0; i < 3; i++ )
  {
    ConnectInv[18*i+0] = 18*i+0;
    ConnectInv[18*i+1] = 18*i+1;
    ConnectInv[18*i+2] = 18*i+5;
    ConnectInv[18*i+3] = 18*i+6;
    ConnectInv[18*i+4] = 18*i+7;
    ConnectInv[18*i+5] = 18*i+11;
    ConnectInv[18*i+6] = 18*i+12;
    ConnectInv[18*i+7] = 18*i+13;
    ConnectInv[18*i+8] = 18*i+17;
    ConnectInv[18*i+9] = 18*i+2;
    ConnectInv[18*i+10] = 18*i+3;
    ConnectInv[18*i+11] = 18*i+4;
    ConnectInv[18*i+12] = 18*i+8;
    ConnectInv[18*i+13] = 18*i+9;
    ConnectInv[18*i+14] = 18*i+10;
    ConnectInv[18*i+15] = 18*i+14;
    ConnectInv[18*i+16] = 18*i+15;
    ConnectInv[18*i+17] = 18*i+16;
  }
}


//-----------------------------------------------------------------------
//   getPanRibElements_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::getPanRibElements_

  ( idx_t&              ielem_panSideA,
    idx_t&              ielem_panSideB,
    idx_t&              ielem_rib,
    const ElemSet&      pan_elems_,
    const ElemSet&      rib_elems_,
    const Matrix&       coords ) const
{
  IdxVector   inodes_pan       ( nodeCount_panrib_ );
  IdxVector   inodes_rib       ( nodeCount_panrib_ );
  Matrix      coords_pan       ( qRank_, nodeCount_panrib_ );
  Matrix      coords_rib       ( qRank_, nodeCount_panrib_ );
  int count;

  // Find the panel element that is connected to the BLine2 element from side A.

  for ( idx_t pan_ie = 0; pan_ie < pan_elemCount_; pan_ie++ )
  {
    idx_t  pan_ielem = pan_ielems_[pan_ie];

    // Get the panel element coordinates.

    pan_elems_.getElemNodes  ( inodes_pan, pan_ielem );
    pan_nodes_.getSomeCoords ( coords_pan, inodes_pan );

    // Check if the panel element is connected to the BLine2 element.

    count = 0;
    for (idx_t i = 0; i < coords_pan.size(1); ++i) {
      for (idx_t j = 0; j < coords.size(1); ++j) {
        if (coords_pan(0, i) == coords(0, j) && 
            coords_pan(1, i) == coords(1, j) && 
            coords_pan(2, i) == coords(2, j)) {
          count++;
        }
      }
    }
    if (count == 2) {
      // Found the first panel element connected to the BLine2 element.
      ielem_panSideA = pan_ielem;
      break;
    }
  }

  // Find the panel element that is connected to the BLine2 element from side B.

  for ( idx_t pan_ie = 0; pan_ie < pan_elemCount_; pan_ie++ )
  {
    idx_t  pan_ielem = pan_ielems_[pan_ie];

    // Get the panel element coordinates.

    pan_elems_.getElemNodes  ( inodes_pan, pan_ielem );
    pan_nodes_.getSomeCoords ( coords_pan, inodes_pan );

    // Check if the panel element is connected to the BLine2 element and is not the ielem_panSideA.

    count = 0;
    for (idx_t i = 0; i < coords_pan.size(1); ++i) {
      for (idx_t j = 0; j < coords.size(1); ++j) {
        if (coords_pan(0, i) == coords(0, j) && 
            coords_pan(1, i) == coords(1, j) && 
            coords_pan(2, i) == coords(2, j)) {
          count++;
        }
      }
    }
    if (count == 2 && pan_ielem != ielem_panSideA) {
      // Found the second panel element connected to the BLine2 element.
      ielem_panSideB = pan_ielem;
      break;
    }
  }

  // Find the rib element that is connected to the BLine2 element.

  for ( idx_t rib_ie = 0; rib_ie < rib_elemCount_; rib_ie++ )
  {
    idx_t  rib_ielem = rib_ielems_[rib_ie];

    // Get the ribs element coordinates.

    rib_elems_.getElemNodes  ( inodes_rib, rib_ielem );
    rib_nodes_.getSomeCoords ( coords_rib, inodes_rib );

    // Check if the ribs element is connected to the BLine2 element.

    count = 0;
    for (idx_t i = 0; i < coords_rib.size(1); ++i) {
      for (idx_t j = 0; j < coords.size(1); ++j) {
        if (coords_rib(0, i) == coords(0, j) && 
            coords_rib(1, i) == coords(1, j) && 
            coords_rib(2, i) == coords(2, j)) {
          count++;
        }
      }
    }
    if (count == 2) {
      ielem_rib = rib_ielem;
      break;
    }
  }
}


//-----------------------------------------------------------------------
//   getTransMatrix_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::getTransMatrix_
   
  ( Matrix&             transMat,
    const Matrix&       coords,
    const Vector&       v   ) const

{
  using jem::ALL;
  using jem::numeric::matmul;

  Vector      e1       ( qRank_ );
  Vector      e2       ( qRank_ );
  Vector      e3       ( qRank_ );
  Vector      e1_l     ( qRank_ );
  Vector      e2_l     ( qRank_ );
  Vector      e3_l     ( qRank_ );

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
}


//-----------------------------------------------------------------------
//   getTransMatrixDof_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::getTransMatrixDof_
   
  ( Matrix&           transMatDof,
    const Matrix&     transMat_pan,
    const Matrix&     transMat_rib   ) const

{
  Matrix      transMat_WiRi   ( nodeRank_, nodeRank_ );
  Matrix      transMatRot_pan ( nodeRank_, nodeRank_ );
  Matrix      transMatRot_rib ( nodeRank_, nodeRank_ );
  Matrix      transMatBlock   ( 2*nodeRank_, 2*nodeRank_ );

  // Get the transformation matrix for rotational Dofs.

  transMat_WiRi = 0.0;
  transMat_WiRi(0,0) = 0.0;
  transMat_WiRi(0,1) = -1.0;
  transMat_WiRi(1,0) = 1.0;
  transMat_WiRi(1,1) = 0.0;
  transMat_WiRi(2,2) = 1.0;

  transMatRot_pan = matmul ( transMat_WiRi, transMat_pan );
  transMatRot_rib = matmul ( transMat_WiRi, transMat_rib );

  // Get the transformation matrix of the Dofs of a single node for the panel local system.

  transMatBlock = 0.0;
  transMatBlock(slice(0,nodeRank_),slice(0,nodeRank_)) = transMat_pan;
  transMatBlock(slice(nodeRank_,2*nodeRank_),slice(nodeRank_,2*nodeRank_)) = transMatRot_pan;

  transMatDof = 0.0;

  // Get the panelA and Panel B part of the 54x54 transformation matrix of the nodal element Dofs.

  transMatDof(slice(0,2*nodeRank_),slice(0,2*nodeRank_)) = transMatBlock;
  transMatDof(slice(2*nodeRank_,4*nodeRank_),slice(2*nodeRank_,4*nodeRank_)) = transMatBlock;
  transMatDof(slice(4*nodeRank_,6*nodeRank_),slice(4*nodeRank_,6*nodeRank_)) = transMatBlock;

  transMatDof(slice(6*nodeRank_,8*nodeRank_),slice(6*nodeRank_,8*nodeRank_)) = transMatBlock;
  transMatDof(slice(8*nodeRank_,10*nodeRank_),slice(8*nodeRank_,10*nodeRank_)) = transMatBlock;
  transMatDof(slice(10*nodeRank_,12*nodeRank_),slice(10*nodeRank_,12*nodeRank_)) = transMatBlock;

  // Get the transformation matrix of the Dofs of a single node for the ribs local system.

  transMatBlock = 0.0;
  transMatBlock(slice(0,nodeRank_),slice(0,nodeRank_)) = transMat_rib;
  transMatBlock(slice(nodeRank_,2*nodeRank_),slice(nodeRank_,2*nodeRank_)) = transMatRot_rib;

  // Get the ribs part of the 54x54 transformation matrix of the nodal element Dofs.

  transMatDof(slice(12*nodeRank_,14*nodeRank_),slice(12*nodeRank_,14*nodeRank_)) = transMatBlock;
  transMatDof(slice(14*nodeRank_,16*nodeRank_),slice(14*nodeRank_,16*nodeRank_)) = transMatBlock;
  transMatDof(slice(16*nodeRank_,18*nodeRank_),slice(16*nodeRank_,18*nodeRank_)) = transMatBlock;
}


//-----------------------------------------------------------------------
//   get2DLocalcoordinates_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::get2DLocalcoordinates_

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

  // Geting the plane xcoords from the xcoords3D.

  xcoords(0, ALL) = xcoords3D(0, ALL);
  xcoords(1, ALL) = xcoords3D(1, ALL);
}


//-----------------------------------------------------------------------
//   get2DGlobalcoordinates_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::get2DGlobalcoordinates_

  ( Matrix&  coords2D,
    const Matrix&   transMat,
    const Matrix&   coords )   const
    
{
  using jem::ALL;
  using jem::numeric::matmul;
  Matrix  coords3D ( coords.size(0), coords.size(1) );

  coords3D(0, ALL) = coords(0, ALL);
  coords3D(1, ALL) = coords(1, ALL);
  coords3D(2, ALL) = coords(2, ALL);

  // Apply the transformation to the coords3D.

  for ( idx_t i = 0; i < coords2D.size(1); i++ )
  {
    coords3D(ALL, i) = matmul ( transMat, coords3D(ALL, i) );
  }

  // Geting the plane coords2D from the coords3D.

  coords2D(0, ALL) = coords3D(0, ALL);
  coords2D(1, ALL) = coords3D(1, ALL);
}


//-----------------------------------------------------------------------
//   getArea_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::getArea_

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
//   getGmat_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::getGmat_

  ( Matrix&  Gmat,
    const Matrix&   xcoords,
    const double&   Area ) const
  
{
  double lambda;
  double x_c, y_c;

  Vector Xi  ( nodeRank_ );
  Vector Eta ( nodeRank_ );
  Vector xm  ( nodeCount_panrib_ );
  Vector ym  ( nodeCount_panrib_ );
  Vector lm  ( nodeCount_panrib_ );
  Vector S   ( nodeCount_panrib_ );
  Vector C   ( nodeCount_panrib_ );
  Matrix a   ( nodeRank_, nodeCount_panrib_ );
  Matrix b   ( nodeRank_, nodeCount_panrib_ );

  Matrix Grc   ( 3*nodeCount_panrib_, 2*nodeCount_panrib_ );
  Matrix Gh    ( 3*nodeCount_panrib_, nodeCount_panrib_ );
  Matrix Gh_11  ( nodeCount_panrib_, 1 );
  Matrix Gh_12  ( nodeCount_panrib_, 1 );
  Matrix Gh_13  ( nodeCount_panrib_, 1 );
  Matrix Gh_21  ( nodeCount_panrib_, 1 );
  Matrix Gh_22  ( nodeCount_panrib_, 1 );
  Matrix Gh_23  ( nodeCount_panrib_, 1 );
  Matrix Gh_31  ( nodeCount_panrib_, 1 );
  Matrix Gh_32  ( nodeCount_panrib_, 1 );
  Matrix Gh_33  ( nodeCount_panrib_, 1 );

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
  Gh(slice(0,nodeCount_panrib_),0) = Gh_11(slice(0,nodeCount_panrib_),0);  
  Gh(slice(0,nodeCount_panrib_),1) = Gh_12(slice(0,nodeCount_panrib_),0);  
  Gh(slice(0,nodeCount_panrib_),2) = Gh_13(slice(0,nodeCount_panrib_),0);  
    
  Gh(slice(nodeCount_panrib_,2*nodeCount_panrib_),0) = Gh_21(slice(0,nodeCount_panrib_),0);  
  Gh(slice(nodeCount_panrib_,2*nodeCount_panrib_),1) = Gh_22(slice(0,nodeCount_panrib_),0);  
  Gh(slice(nodeCount_panrib_,2*nodeCount_panrib_),2) = Gh_23(slice(0,nodeCount_panrib_),0);  

  Gh(slice(2*nodeCount_panrib_,3*nodeCount_panrib_),0) = Gh_31(slice(0,nodeCount_panrib_),0);  
  Gh(slice(2*nodeCount_panrib_,3*nodeCount_panrib_),1) = Gh_32(slice(0,nodeCount_panrib_),0);  
  Gh(slice(2*nodeCount_panrib_,3*nodeCount_panrib_),2) = Gh_33(slice(0,nodeCount_panrib_),0);  

  // Compute the Gmat matrix
  Gmat = 0.0;
  Gmat(slice(0,3*nodeCount_panrib_),slice(0,2*nodeCount_panrib_)) = Grc;
  Gmat(slice(0,3*nodeCount_panrib_),slice(2*nodeCount_panrib_,3*nodeCount_panrib_)) = Gh;
}


//-----------------------------------------------------------------------
//   GmatInverseError_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::GmatInverseError_

  ( double&  error,
    const Matrix&   Gmat,
    const Matrix&   Ginv ) const
  
{
  // Get the analytical Identity.
  Matrix Identity   ( 3*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Identity = 0.0;
  for ( idx_t i = 0; i < 3*nodeCount_panrib_; i++ )
  {
    Identity(i,i) = 1.0;
  }

  // Get the numerical Identity.
  Matrix Identity_num   ( 3*nodeCount_panrib_, 3*nodeCount_panrib_ );
  Identity_num = matmul ( Gmat, Ginv );

  // Compute the error.
  error = 0.0;
  for ( idx_t i = 0; i < 3*nodeCount_panrib_; i++ )
  {
    for ( idx_t j = 0; j < 3*nodeCount_panrib_; j++ )
    {
      error += pow((Identity(i,j) - Identity_num(i,j)),2);      
    }
  }
  error = sqrt (error);
}


//-----------------------------------------------------------------------
//   get3DLocalcoordinates_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::get3DLocalcoordinates_

  ( Matrix&         ipxcoords_pr,
    const Matrix&   transMat_pr,
    const Matrix&   ipxcoords,
    const Matrix&   transMat_int,
    const Matrix&   coords_pr,
    const Matrix&   coords )   const
    
{
  using jem::ALL;
  using jem::numeric::matmul;
  Matrix  xcoords ( qRank_, nodeCount_BLine2_ );
  Matrix  ipxcoords3D ( qRank_, ipCount_ );
  Matrix  ipcoords3D  ( qRank_, ipCount_ );
  Matrix  transMatT_int = transMat_int.transpose ();

  // Get the cg of the BLine2 element.

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

  // Translation to the baricenter of the BLine2 element.

  xcoords(0, ALL) = coords(0, ALL) - Xcg;
  xcoords(1, ALL) = coords(1, ALL) - Ycg;
  xcoords(2, ALL) = coords(2, ALL) - Zcg;

  // Apply the transformation transMat_int to the 3Dxcoords.

  for ( idx_t i = 0; i < xcoords.size(1); i++ )
  {
    xcoords(ALL, i) = matmul ( transMat_int, xcoords(ALL, i) );
  }

  // Geting the 3D ipxcoords from the 2D ipxcoords and the constant third coordinate from xcoords.

  ipxcoords3D(0,ALL) = ipxcoords(0,ALL);
  ipxcoords3D(1,ALL) = ipxcoords(1,ALL);
  ipxcoords3D(2,ALL) = xcoords(2, 0);   // Assuming the third coordinate is constant for the panel-rib interface.


  // Apply the inverse of the transformation transMat_int to ipxcoords3D.

  for ( idx_t i = 0; i < ipxcoords3D.size(1); i++ )
  {
    ipcoords3D(ALL, i) = matmul ( transMatT_int, ipxcoords3D(ALL, i) );
  }

  // Translation back to the origin of the Global system.

  ipcoords3D(0, ALL) += Xcg;
  ipcoords3D(1, ALL) += Ycg;
  ipcoords3D(2, ALL) += Zcg;

  // Get the cg of the pr element.

  Xcg = 0.0;
  Ycg = 0.0;
  Zcg = 0.0;
  for ( idx_t i = 0; i < coords_pr.size(1); i++ )
  {
    Xcg += coords_pr(0, i);
    Ycg += coords_pr(1, i);
    Zcg += coords_pr(2, i);
  }
  Xcg /= coords_pr.size(1);
  Ycg /= coords_pr.size(1);
  Zcg /= coords_pr.size(1);

  // Translation to the baricenter of the pr element.

  ipcoords3D(0, ALL) -= Xcg;
  ipcoords3D(1, ALL) -= Ycg;
  ipcoords3D(2, ALL) -= Zcg;

  // Apply the transformation transMat_pr to the ipcoords3D.

  for ( idx_t i = 0; i < ipcoords3D.size(1); i++ )
  {
    ipxcoords_pr(ALL, i) = matmul ( transMat_pr, ipcoords3D(ALL, i) );
  }
}


//-----------------------------------------------------------------------
//   getBABxBymats_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::getBABxBymats_

  ( Matrix&  BAmat,
    Matrix&  Bxmat,
    Matrix&  Bymat )   const
    
{
  BAmat = 0.0;
  BAmat(0,0) = 1.0;
  BAmat(1,3) = 1.0;
  BAmat(2,6) = 1.0;

  Bxmat = 0.0;
  Bxmat(0,1) = 1.0;

  Bymat = 0.0;
  Bymat(0,2) = 1.0;
}


//-----------------------------------------------------------------------
//   getMAMamats_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::getMAMamats_

  ( Matrix&  MAmat,
    Matrix&  Mamat,
    const Matrix&   xcoords ) const

{
  using jem::ALL;
  Vector x ( nodeCount_panrib_);
  Vector y ( nodeCount_panrib_);

  x[ALL] = xcoords(0, ALL);
  y[ALL] = xcoords(1, ALL);

  // get MA matrix.
  MAmat = 0.0;
  MAmat(ALL, 0) = 1.0;
  MAmat(ALL, 1) = x[ALL];
  MAmat(ALL, 2) = y[ALL];

  // get Ma matrix.
  Mamat = 0.0;
  for ( idx_t i = 0; i < xcoords.size(1); ++i )
  {
    Mamat(i, 0) = pow(x[i],2);
    Mamat(i, 1) = x[i]*y[i];
    Mamat(i, 2) = pow(y[i],2);
    Mamat(i, 3) = pow(x[i],3);
    Mamat(i, 4) = pow(x[i],2)*y[i];
    Mamat(i, 5) = x[i]*pow(y[i],2);
    Mamat(i, 6) = pow(y[i],3);
  }
}


//-----------------------------------------------------------------------
//   invertMAmat_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::invertMAmat_

  ( Matrix&  MAinv,
    const Matrix&   MAmat_ ) const
{
  double det = MAmat_(0, 0)*MAmat_(1, 1)*MAmat_(2, 2) + MAmat_(0, 1)*MAmat_(1, 2)*MAmat_(2, 0) +
               MAmat_(0, 2)*MAmat_(1, 0)*MAmat_(2, 1) - MAmat_(0, 2)*MAmat_(1, 1)*MAmat_(2, 0) -
               MAmat_(0, 1)*MAmat_(1, 0)*MAmat_(2, 2) - MAmat_(0, 0)*MAmat_(1, 2)*MAmat_(2, 1);

  if ( fabs(det) < 1.0e-10 )
  {
    throw jem::Exception ( "StructuralInterfaceModel::invertMAmat_",
                             "Matrix is singular" );
  }

  MAinv(0, 0) = (MAmat_(1, 1)*MAmat_(2, 2) - MAmat_(1, 2)*MAmat_(2, 1)) / det;
  MAinv(0, 1) = (MAmat_(0, 2)*MAmat_(2, 1) - MAmat_(0, 1)*MAmat_(2, 2)) / det;
  MAinv(0, 2) = (MAmat_(0, 1)*MAmat_(1, 2) - MAmat_(0, 2)*MAmat_(1, 1)) / det;
  MAinv(1, 0) = (MAmat_(1, 2)*MAmat_(2, 0) - MAmat_(1, 0)*MAmat_(2, 2)) / det;
  MAinv(1, 1) = (MAmat_(0, 0)*MAmat_(2, 2) - MAmat_(0, 2)*MAmat_(2, 0)) / det;
  MAinv(1, 2) = (MAmat_(0, 2)*MAmat_(1, 0) - MAmat_(0, 0)*MAmat_(1, 2)) / det;
  MAinv(2, 0) = (MAmat_(1, 0)*MAmat_(2, 1) - MAmat_(1, 1)*MAmat_(2, 0)) / det;
  MAinv(2, 1) = (MAmat_(0, 1)*MAmat_(2, 0) - MAmat_(0, 0)*MAmat_(2, 1)) / det;
  MAinv(2, 2) = (MAmat_(0, 0)*MAmat_(1, 1) - MAmat_(0, 1)*MAmat_(1, 0)) / det;
}


//-----------------------------------------------------------------------
//   getTmat_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::getTmat_

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


void PanelRibsInterface6DofsBFModel::getMembraneDisp_

  ( Vector&         Disp_mem_panA,
    Vector&         Disp_mem_panB,
    Vector&         Disp_mem_rib,
    const Vector&   elemxDisp0 )   const
    
{
  for ( idx_t i = 0; i < nodeCount_; i++ )
  {
    if ( i < nodeCount_panrib_)
    // Panel A membrane displacements.
    {
      Disp_mem_panA[2*i]   = elemxDisp0[5*i];
      Disp_mem_panA[2*i+1] = elemxDisp0[5*i+1];
    }
    else if ( i < 2*nodeCount_panrib_ )
    // Panel B membrane displacements.
    {
      idx_t ii = i - nodeCount_panrib_;
      Disp_mem_panB[2*ii]   = elemxDisp0[5*i];
      Disp_mem_panB[2*ii+1] = elemxDisp0[5*i+1];
    }
    else
    // Ribs membrane displacements.
    {
      idx_t ii = i - 2*nodeCount_panrib_;
      Disp_mem_rib[2*ii]   = elemxDisp0[5*i];
      Disp_mem_rib[2*ii+1] = elemxDisp0[5*i+1];
    }
  }
}


//-----------------------------------------------------------------------
//   getHmatAnalytic_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::getHmatAnalytic_

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
//   getBmat_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::getBmat_

  ( Matrix&  Bmat,
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

 // Computing the Bmat Matrix.

  Bmat = 0.0;
  Bmat(ALL,0) =  BR_1(ALL,0);
  Bmat(ALL,1) =  BR_2(ALL,0);
  Bmat(ALL,2) =  BR_3(ALL,0);
  Bmat(ALL,3) =  BVn_12(ALL,0);
  Bmat(ALL,4) =  BVn_23(ALL,0);
  Bmat(ALL,5) =  BVn_31(ALL,0);
  Bmat(ALL,6) =  BMn_12(ALL,0);
  Bmat(ALL,7) =  BMn_21(ALL,0);
  Bmat(ALL,8) =  BMn_23(ALL,0);
  Bmat(ALL,9) =  BMn_32(ALL,0);
  Bmat(ALL,10) = BMn_31(ALL,0);
  Bmat(ALL,11) = BMn_13(ALL,0);
}


//-----------------------------------------------------------------------
//   getNrNcNhmats_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::getNrNcNhmats_

  ( Matrix&  Nr,
    Matrix&  Nc,
    Matrix&  Nh,
    const Matrix&   xcoords,
    const double&   Area,
    const double&   x,
    const double&   y ) const
  
{
  double lambda;
  double x_c, y_c;

  double xi;
  double eta;

  Vector Xi  ( nodeRank_ );
  Vector Eta ( nodeRank_ );
  Vector xm  ( nodeCount_panrib_ );
  Vector ym  ( nodeCount_panrib_ );
  Vector lm  ( nodeCount_panrib_ );
  Vector S   ( nodeCount_panrib_ );
  Vector C   ( nodeCount_panrib_ );
  Matrix a   ( nodeRank_, nodeCount_panrib_ );
  Matrix b   ( nodeRank_, nodeCount_panrib_ );

  // Compute geometric parameters for the Nr and Nc matrices.
  lambda = 1.0 / sqrt(Area);
  x_c = (xcoords(0, 0) + xcoords(0, 1) + xcoords(0, 2)) / 3.0;
  y_c = (xcoords(1, 0) + xcoords(1, 1) + xcoords(1, 2)) / 3.0;
  xi = lambda*(x - x_c);
  eta = lambda*(y - y_c);

  // Compute the Nr matrix
  Nr = 0.0;
  Nr(0, 0) = 1.0;
  Nr(0, 2) = -eta;
  Nr(1, 1) = 1.0;
  Nr(1, 2) = xi;

  // Compute the Nc matrix
  Nc = 0.0;
  Nc(0, 0) = xi;
  Nc(0, 2) = eta;
  Nc(1, 1) = eta;
  Nc(1, 2) = xi;

  // Compute geometric parameters for the Nh matrix.
  Xi[0] = lambda*(xcoords(0, 0) - x_c);
  Xi[1] = lambda*(xcoords(0, 1) - x_c);
  Xi[2] = lambda*(xcoords(0, 2) - x_c);
  
  Eta[0] = lambda*(xcoords(1, 0) - y_c);
  Eta[1] = lambda*(xcoords(1, 1) - y_c);
  Eta[2] = lambda*(xcoords(1, 2) - y_c);

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

  // Compute the Nh matrix
  Nh = 0.0;
  Nh(0, 0) = (a(0,0)*xi*xi+a(1,0)*xi*eta+a(2,0)*eta*eta);
  Nh(0, 1) = (a(0,1)*xi*xi+a(1,1)*xi*eta+a(2,1)*eta*eta);
  Nh(0, 2) = (a(0,2)*xi*xi+a(1,2)*xi*eta+a(2,2)*eta*eta);
  Nh(1, 0) = (b(0,0)*xi*xi+b(1,0)*xi*eta+b(2,0)*eta*eta);
  Nh(1, 1) = (b(0,1)*xi*xi+b(1,1)*xi*eta+b(2,1)*eta*eta);
  Nh(1, 2) = (b(0,2)*xi*xi+b(1,2)*xi*eta+b(2,2)*eta*eta);
}


//-----------------------------------------------------------------------
//   getNuNvmats_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::getNuNvmats_

  ( Matrix&           Nu,
    Matrix&           Nv,
    const Matrix&     Nr,
    const Matrix&     Hr,
    const Matrix&     Nc,
    const Matrix&     Hc,
    const Matrix&     Nh,
    const Matrix&     Hh ) const
{
  Matrix Nuv   ( rank_, 3*nodeCount_panrib_ );

  // Compute the Nuv matrix
  Nuv = 0.0;
  Nuv += matmul( Nr, Hr );
  Nuv += matmul( Nc, Hc );
  Nuv += matmul( Nh, Hh );

  // Compute the Nu and Nv matrices
  Nu = Nuv(slice(0,1),slice(0,3*nodeCount_panrib_)); 
  Nv = Nuv(slice(1,2),slice(0,3*nodeCount_panrib_)); 
}


//-----------------------------------------------------------------------
//   getNrxNcxNhxmats_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::getNrxNcxNhxmats_

  ( Matrix&  Nrx,
    Matrix&  Ncx,
    Matrix&  Nhx,
    const Matrix&   xcoords,
    const double&   Area,
    const double&   x,
    const double&   y ) const
  
{
  double lambda;
  double x_c, y_c;

  double xi;
  double eta;

  Vector Xi  ( nodeRank_ );
  Vector Eta ( nodeRank_ );
  Vector xm  ( nodeCount_panrib_ );
  Vector ym  ( nodeCount_panrib_ );
  Vector lm  ( nodeCount_panrib_ );
  Vector S   ( nodeCount_panrib_ );
  Vector C   ( nodeCount_panrib_ );
  Matrix a   ( nodeRank_, nodeCount_panrib_ );
  Matrix b   ( nodeRank_, nodeCount_panrib_ );

  // Compute geometric parameters for the Nr and Nc matrices.
  lambda = 1.0 / sqrt(Area);
  x_c = (xcoords(0, 0) + xcoords(0, 1) + xcoords(0, 2)) / 3.0;
  y_c = (xcoords(1, 0) + xcoords(1, 1) + xcoords(1, 2)) / 3.0;
  xi = lambda*(x - x_c);
  eta = lambda*(y - y_c);

  // Compute the Nrx matrix
  Nrx = 0.0;
  Nrx(1, 2) = lambda;

  // Compute the Ncx matrix
  Ncx = 0.0;
  Ncx(0, 0) = lambda;
  Ncx(1, 2) = lambda;

  // Compute geometric parameters for the Nh matrix.
  Xi[0] = lambda*(xcoords(0, 0) - x_c);
  Xi[1] = lambda*(xcoords(0, 1) - x_c);
  Xi[2] = lambda*(xcoords(0, 2) - x_c);
  
  Eta[0] = lambda*(xcoords(1, 0) - y_c);
  Eta[1] = lambda*(xcoords(1, 1) - y_c);
  Eta[2] = lambda*(xcoords(1, 2) - y_c);

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

  // Compute the Nhx matrix
  Nhx = 0.0;
  Nhx(0, 0) = (a(0,0)*2.0*xi*lambda+a(1,0)*eta*lambda);
  Nhx(0, 1) = (a(0,1)*2.0*xi*lambda+a(1,1)*eta*lambda);
  Nhx(0, 2) = (a(0,2)*2.0*xi*lambda+a(1,2)*eta*lambda);
  Nhx(1, 0) = (b(0,0)*2.0*xi*lambda+b(1,0)*eta*lambda);
  Nhx(1, 1) = (b(0,1)*2.0*xi*lambda+b(1,1)*eta*lambda);
  Nhx(1, 2) = (b(0,2)*2.0*xi*lambda+b(1,2)*eta*lambda);
}


//-----------------------------------------------------------------------
//   getNryNcyNhymats_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::getNryNcyNhymats_

  ( Matrix&  Nry,
    Matrix&  Ncy,
    Matrix&  Nhy,
    const Matrix&   xcoords,
    const double&   Area,
    const double&   x,
    const double&   y ) const
  
{
  double lambda;
  double x_c, y_c;

  double xi;
  double eta;

  Vector Xi  ( nodeRank_ );
  Vector Eta ( nodeRank_ );
  Vector xm  ( nodeCount_panrib_ );
  Vector ym  ( nodeCount_panrib_ );
  Vector lm  ( nodeCount_panrib_ );
  Vector S   ( nodeCount_panrib_ );
  Vector C   ( nodeCount_panrib_ );
  Matrix a   ( nodeRank_, nodeCount_panrib_ );
  Matrix b   ( nodeRank_, nodeCount_panrib_ );

  // Compute geometric parameters for the Nr and Nc matrices.
  lambda = 1.0 / sqrt(Area);
  x_c = (xcoords(0, 0) + xcoords(0, 1) + xcoords(0, 2)) / 3.0;
  y_c = (xcoords(1, 0) + xcoords(1, 1) + xcoords(1, 2)) / 3.0;
  xi = lambda*(x - x_c);
  eta = lambda*(y - y_c);

  // Compute the Nry matrix
  Nry = 0.0;
  Nry(0, 2) = -lambda;

  // Compute the Ncy matrix
  Ncy = 0.0;
  Ncy(0, 2) = lambda;
  Ncy(1, 1) = lambda;

  // Compute geometric parameters for the Nh matrix.
  Xi[0] = lambda*(xcoords(0, 0) - x_c);
  Xi[1] = lambda*(xcoords(0, 1) - x_c);
  Xi[2] = lambda*(xcoords(0, 2) - x_c);
  
  Eta[0] = lambda*(xcoords(1, 0) - y_c);
  Eta[1] = lambda*(xcoords(1, 1) - y_c);
  Eta[2] = lambda*(xcoords(1, 2) - y_c);

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

  // Compute the Nhy matrix
  Nhy = 0.0;
  Nhy(0, 0) = (a(1,0)*xi*lambda+a(2,0)*2.0*eta*lambda);
  Nhy(0, 1) = (a(1,1)*xi*lambda+a(2,1)*2.0*eta*lambda);
  Nhy(0, 2) = (a(1,2)*xi*lambda+a(2,2)*2.0*eta*lambda);
  Nhy(1, 0) = (b(1,0)*xi*lambda+b(2,0)*2.0*eta*lambda);
  Nhy(1, 1) = (b(1,1)*xi*lambda+b(2,1)*2.0*eta*lambda);
  Nhy(1, 2) = (b(1,2)*xi*lambda+b(2,2)*2.0*eta*lambda);
}


//-----------------------------------------------------------------------
//   getNvxmat_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::getNvxmat_

  ( Matrix&           Nvx,
    const Matrix&     Nrx,
    const Matrix&     Hr,
    const Matrix&     Ncx,
    const Matrix&     Hc,
    const Matrix&     Nhx,
    const Matrix&     Hh ) const
{
  Matrix Nuvx   ( rank_, 3*nodeCount_panrib_ );

  // Compute the Nuvx matrix
  Nuvx = 0.0;
  Nuvx += matmul( Nrx, Hr );
  Nuvx += matmul( Ncx, Hc );
  Nuvx += matmul( Nhx, Hh );

  // Compute the Nvx matrix
  Nvx = Nuvx(slice(1,2),slice(0,3*nodeCount_panrib_)); 
}


//-----------------------------------------------------------------------
//   getNuymat_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::getNuymat_

  ( Matrix&           Nuy,
    const Matrix&     Nry,
    const Matrix&     Hr,
    const Matrix&     Ncy,
    const Matrix&     Hc,
    const Matrix&     Nhy,
    const Matrix&     Hh ) const
{
  Matrix Nuvy   ( rank_, 3*nodeCount_panrib_ );

  // Compute the Nuvy matrix
  Nuvy = 0.0;
  Nuvy += matmul( Nry, Hr );
  Nuvy += matmul( Ncy, Hc );
  Nuvy += matmul( Nhy, Hh );

  // Compute the Nuy matrix
  Nuy = Nuvy(slice(0,1),slice(0,3*nodeCount_panrib_)); 
}


//-----------------------------------------------------------------------
//   getSRRxRymats_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::getSRRxRymats_

  ( Matrix&           Smat,
    Matrix&           Rmat,
    Matrix&           Rxmat,
    Matrix&           Rymat,
    const double&     x,
    const double&     y ) const
{
  Smat(0,0) = 1.0;
  Smat(1,0) = x;
  Smat(2,0) = y;

  Rmat(0,0) = pow(x,2);
  Rmat(1,0) = x*y;
  Rmat(2,0) = pow(y,2);
  Rmat(3,0) = pow(x,3);
  Rmat(4,0) = pow(x,2)*y;
  Rmat(5,0) = x*pow(y,2);
  Rmat(6,0) = pow(y,3);

  Rxmat(0,0) = 2.0*x;
  Rxmat(1,0) = y;
  Rxmat(2,0) = 0.0;
  Rxmat(3,0) = 3*pow(x,2);
  Rxmat(4,0) = 2.0*x*y;
  Rxmat(5,0) = pow(y,2);
  Rxmat(6,0) = 0.0;

  Rymat(0,0) = 0.0;
  Rymat(1,0) = x;
  Rymat(2,0) = 2.0*y;
  Rymat(3,0) = 0.0;
  Rymat(4,0) = pow(x,2);
  Rymat(5,0) = 2.0*x*y;
  Rymat(6,0) = 3.0*pow(y,2);
}


//-----------------------------------------------------------------------
//   initWriter_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::initWriter_  

  ( const Properties&  params )

{
  // Open file for writeXOutput

  StringVector fileName;
  String       prepend;
  
  if ( params.find( prepend, "prepend" ) )
  { 
    fileName.resize(3);

    fileName[0] = prepend;
    fileName[1] = myTag_;
    fileName[2] = "dat";
  }
  else
  {
    fileName.resize(2);

    fileName[0] = myTag_;
    fileName[1] = "dat";
  }

  xOut_ = newInstance<PrintWriter>(
          newInstance<FileWriter> ( StringUtils::join( fileName, "." ) ) );
}


//-----------------------------------------------------------------------
//   writeGeom_
//-----------------------------------------------------------------------


void PanelRibsInterface6DofsBFModel::writeGeom_  ()  const

{

  // write integration points coordinates

  Matrix      coords     ( qRank_, nodeCount_BLine2_ );
              coords     = 0.;

  Matrix      ipCoords   ( qRank_, ipCount_ );

  IdxVector   inodes     ( nodeCount_BLine2_ );

  *xOut_ << "ipCoords" << '\n';
  *xOut_ << rank_ << '\n';

  for ( idx_t ie = 0; ie < elemCount_; ie++ )
  {
    // Get the element coordinates.

    elems_.getElemNodes  ( inodes, ielems_[ie] );
    nodes_.getSomeCoords ( coords( slice(0,qRank_) , ALL ), 
                           inodes );

    // Get the coordinates of the integration points of the BLine2 element

    shape_->getGlobalIntegrationPoints ( ipCoords, coords( slice(0,rank_) , ALL ) );

    for ( idx_t ip = 0; ip < ipCoords.size( 1 ); ++ip )
    {
      *xOut_ << ie << " " << ip << " ";
      for ( idx_t j = 0; j < qRank_; ++j )
      {
        *xOut_ << ipCoords( j, ip ) << " ";
      }
      *xOut_ << '\n';
    }
  }
  xOut_->flush();
}


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newPanelRibsInterface6DofsBFModel
//-----------------------------------------------------------------------


static Ref<Model>     newPanelRibsInterface6DofsBFModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<PanelRibsInterface6DofsBFModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declarePanelRibsInterface6DofsBFModel
//-----------------------------------------------------------------------


void declarePanelRibsInterface6DofsBFModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "PanelRibsInterface6DofsBFModel", & newPanelRibsInterface6DofsBFModel );
}
