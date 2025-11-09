/*
 *
 * Model for structural interface elements: 
 *   - assembly of stiffness and internal force vector
 *   - output
 * 
 * Sergio Gustavo Ferreira Cordeiro, May 2025
 *
 */

#include "StructuralInterfaceModel.h"

#include <jem/base/array/operators.h>
#include <jem/base/array/select.h>
#include <jem/base/Error.h>
#include <jem/base/Float.h>
#include <jem/base/IllegalInputException.h>
#include <jem/base/System.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/numeric/algebra/MatmulChain.h>
#include <jem/numeric/algebra/Cholesky.h>
#include <jem/util/StringUtils.h>
#include <jive/model/ModelFactory.h>
#include <jive/model/Actions.h>
#include <jive/model/StateVector.h>
#include <jive/geom/BoundaryShape.h>
#include <jive/geom/BShapeFactory.h>
#include <jive/geom/BoundaryTriangle.h>
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


const char* StructuralInterfaceModel::DOF_NAMES[5]  = {"u","v","w","wx","wy"};
const char* StructuralInterfaceModel::SHAPE_PROP    = "shape";
const char* StructuralInterfaceModel::COHEMAT_PROP  = "coheMat";
const char* StructuralInterfaceModel::BOTTOP_SHAPE_PROP   = "shape";
const char* StructuralInterfaceModel::MATERIAL_PROP[2]    = {"material","material"};
const char* StructuralInterfaceModel::THICK_PROP[2]       = {"thickness","thickness"};
const char* StructuralInterfaceModel::BOT_ORIENTATION_PROP[3] = {"v1","v2","v3"};  
const char* StructuralInterfaceModel::TOP_ORIENTATION_PROP[3] = {"v1","v2","v3"};  
const char* StructuralInterfaceModel::BOT_NAME_PROP = "bot_name";
const char* StructuralInterfaceModel::TOP_NAME_PROP = "top_name";  


//-----------------------------------------------------------------------
//   constructor
//-----------------------------------------------------------------------


StructuralInterfaceModel::StructuralInterfaceModel

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

  myConf.set ( SHAPE_PROP, shProps );

  bshape = makeBTriangleNC ( shProps );

  shape_ = newInstance<InterfaceShape> ( "interface", bshape );

  // NB: Make sure you know what you're doing when 
  // the rank of the shape doesn't match the rank of the mesh

  qRank_     = shape_->globalRank  ( );
  nodeCount_ = shape_->nodeCount   ( );
  ipCount_   = shape_->ipointCount ( );

  // Make sure that each element has the same number of nodes as the
  // shape object.

  elems_.checkSomeElements 
    ( context,
      ielems_,
      nodeCount_ );

  dofs_ = XDofSpace::get ( nodes_.getData(), globdat );
  
  dofTypes_.resize( 5 );

  for( idx_t i = 0; i < 5; i++)
  {
    dofTypes_[i] = dofs_->addType ( DOF_NAMES[i]);
  }

  dofs_->addDofs (
    elems_.getUniqueNodesOf ( ielems_ ),
    dofTypes_
  );

  getShapeFuncs_ = getShapeFuncsFunc ( rank_ );

  // Create a coheMat object.

  coheMat_ = newCohesiveMat 
           ( COHEMAT_PROP, myConf, myProps, globdat );

  frictionMat_ = dynamicCast<AlfanoTuronCoheMat> ( coheMat_ );

  // Get the bottom and the top Names.

  myProps.find( bot_, BOT_NAME_PROP );
  myConf.set  ( BOT_NAME_PROP, bot_ );
  String botName = joinNames (myBase_, bot_ );
  
  myProps.find( top_, TOP_NAME_PROP );
  myConf.set  ( TOP_NAME_PROP, top_ );
  String topName = joinNames (myBase_, top_ );

  // Get the Pros and Conf from bottom and top submodels.

  Properties  botProps = props.getProps ( botName );
  Properties  botConf  = conf.makeProps ( botName );
  const String botContext = "model `" + botName + "'";

  Properties  topProps = props.getProps ( topName );
  Properties  topConf  = conf.makeProps ( topName );
  const String topContext = "model `" + topName + "'";

  // Create the InternalShape object for the bottom and top plies.

  Properties  bottop_shProps = botProps.getProps ( BOTTOP_SHAPE_PROP );
  botConf.set ( BOTTOP_SHAPE_PROP, bottop_shProps );

  bottopshape_  = IShapeFactory::newInstance(BOTTOP_SHAPE_PROP, botConf, botProps );

  nodeCount_bottop_ = bottopshape_->nodeCount  ( );
  ipCount_bottop_ = bottopshape_->ipointCount ( );

  getShapeGrads_ = getShapeGradsFunc ( rank_ );

  // Create two material model objects: bottom and top plies.

  material_[0] = newMaterial ( MATERIAL_PROP[0], botConf, botProps, globdat );
  material_[0]-> allocPoints  ( ipCount_bottop_ );
  softening_[0] = dynamicCast<Softening> ( material_[0] );
  
  material_[1] = newMaterial ( MATERIAL_PROP[1], topConf, topProps, globdat );
  material_[1]-> allocPoints  ( ipCount_bottop_ );
  softening_[1] = dynamicCast<Softening> ( material_[1] );

  // Create two thickness variables: bottom and top plies.

  botProps.find( thickness_[0], THICK_PROP[0] );
  botConf.set  ( THICK_PROP[0], thickness_[0] );
  
  topProps.find( thickness_[1], THICK_PROP[1] );
  topConf.set  ( THICK_PROP[1], thickness_[1] );

  // Create the variables: bot_orientVec_ and bot_orientVec_.

  botProps.find( bot_orientVec_[0], BOT_ORIENTATION_PROP[0] );
  botConf.set  ( BOT_ORIENTATION_PROP[0], bot_orientVec_[0] );

  botProps.find( bot_orientVec_[1], BOT_ORIENTATION_PROP[1] );
  botConf.set  ( BOT_ORIENTATION_PROP[1], bot_orientVec_[1] );

  botProps.find( bot_orientVec_[2], BOT_ORIENTATION_PROP[2] );
  botConf.set  ( BOT_ORIENTATION_PROP[2], bot_orientVec_[2] );

  topProps.find( top_orientVec_[0], TOP_ORIENTATION_PROP[0] );
  topConf.set  ( TOP_ORIENTATION_PROP[0], top_orientVec_[0] );

  topProps.find( top_orientVec_[1], TOP_ORIENTATION_PROP[1] );
  topConf.set  ( TOP_ORIENTATION_PROP[1], top_orientVec_[1] );

  topProps.find( top_orientVec_[2], TOP_ORIENTATION_PROP[2] );
  topConf.set  ( TOP_ORIENTATION_PROP[2], top_orientVec_[2] );

  crackBandMethod_[0] = false;
  crackBandMethod_[1] = false;
}

StructuralInterfaceModel::~StructuralInterfaceModel()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void StructuralInterfaceModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  Properties  myProps  = props.findProps ( myName_ );
  Properties  matProps = myProps.findProps ( COHEMAT_PROP );

  coheMat_->configure ( matProps, globdat );

  coheMat_->allocPoints  ( ipCount_ * elemCount_ );

  // Get the botName.

  String botName = joinNames (myBase_, bot_ );

  Properties  botProps  = props.findProps ( botName );
  Properties  botMatProps = botProps.findProps ( MATERIAL_PROP[0] );

  material_[0]->configure ( botMatProps );

  // Get the topName.

  String topName = joinNames (myBase_, top_ );

  Properties  topProps  = props.findProps ( topName );
  Properties  topMatProps = topProps.findProps ( MATERIAL_PROP[1] );

  material_[1]->configure ( topMatProps );

  Properties  params;
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void StructuralInterfaceModel::getConfig 

  ( const Properties& conf,
    const Properties& globdat ) const

{
  Properties  myConf  = conf.makeProps ( myName_ );
  Properties  matConf = myConf.makeProps ( COHEMAT_PROP );

  coheMat_->getConfig ( matConf, globdat );

  // Get the botName.

  String botName = joinNames (myBase_, bot_ );

  Properties  botConf  = conf.makeProps ( botName );
  Properties  botMatConf = botConf.makeProps ( MATERIAL_PROP[0] );

  material_[0]->getConfig ( botMatConf );

  // Get the topName.

  String topName = joinNames (myBase_, top_ );

  Properties  topConf  = conf.makeProps ( topName );
  Properties  topMatConf = topConf.makeProps ( MATERIAL_PROP[1] );

  material_[1]->getConfig ( topMatConf );
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool StructuralInterfaceModel::takeAction

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
      Vector disp;
      Vector fDiss;

      StateVector::getOld ( disp, dofs_, globdat );
      globdat.get ( fDiss, SolverNames::DISSIPATION_FORCE );

      getFrictionForce_ ( fDiss, disp );
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


void StructuralInterfaceModel::getMatrix_

  ( Ref<MatrixBuilder>  mbuilder,
    const Vector&       force,
    const Vector&       disp ) const

{
  const idx_t dofCount   = 5 * nodeCount_;

  Matrix      stiff   ( qRank_, qRank_ );
  Matrix      BCEmat  ( qRank_, dofCount );

  Matrix      stiff_bot   ( qRank_, qRank_ );
  Matrix      stiff_top   ( qRank_, qRank_ );
  Matrix      D_plate_bot ( qRank_, qRank_ );
  Matrix      D_plate_top ( qRank_, qRank_ );

  Matrix      coords      ( qRank_, nodeCount_ );
              coords     = 0.;
  Matrix      coords_bot  ( qRank_, nodeCount_bottop_ );
  Matrix      coords_top  ( qRank_, nodeCount_bottop_ );
  Matrix      xcoords_bot ( rank_, nodeCount_bottop_ );
  Matrix      xcoords_top ( rank_, nodeCount_bottop_ );

  Matrix      ipcoords    ( qRank_, ipCount_ );
  Matrix      ipxcoords_bot ( rank_, ipCount_ );
  Matrix      ipxcoords_top ( rank_, ipCount_ );

  Matrix      b_mem_bot       ( qRank_, 2*nodeCount_bottop_ );
  Vector      Disp_mem_bot    ( 2*nodeCount_bottop_ );
  Vector      strain_mem_bot  ( qRank_ );
  Vector      stress_mem_bot  ( qRank_ );
  Matrix      b_mem_top       ( qRank_, 2*nodeCount_bottop_ );
  Vector      Disp_mem_top    ( 2*nodeCount_bottop_ );
  Vector      strain_mem_top  ( qRank_ );
  Vector      stress_mem_top  ( qRank_ );

  Matrix      Bmat_bot        ( 7, 4*nodeCount_bottop_);
  Matrix      Tmat_bot        ( 4*nodeCount_bottop_, 3*nodeCount_bottop_ );
  Matrix      Hmat_bot        ( 7, 7);
  Matrix      Hinv_bot        ( 7, 7);
  Matrix      Bmat_top        ( 7, 4*nodeCount_bottop_);
  Matrix      Tmat_top        ( 4*nodeCount_bottop_, 3*nodeCount_bottop_ );
  Matrix      Hmat_top        ( 7, 7);
  Matrix      Hinv_top        ( 7, 7);

  Matrix      BAmat           ( nodeCount_bottop_, 3 * nodeCount_bottop_ );

  Matrix      MAmat_bot       ( nodeCount_bottop_, rank_+1 );
  Matrix      MAinv_bot       ( nodeCount_bottop_, rank_+1 );
  Matrix      Mamat_bot       ( nodeCount_bottop_, 7 );
  Matrix      Cmat_bot        ( 7, 3 * nodeCount_bottop_ );
  Matrix      BA_MaC_bot      ( nodeCount_bottop_, 3 * nodeCount_bottop_ );

  Matrix      MAmat_top       ( nodeCount_bottop_, rank_+1 );
  Matrix      MAinv_top       ( nodeCount_bottop_, rank_+1 );
  Matrix      Mamat_top       ( nodeCount_bottop_, 7 );
  Matrix      Cmat_top        ( 7, 3 * nodeCount_bottop_ );
  Matrix      BA_MaC_top      ( nodeCount_bottop_, 3 * nodeCount_bottop_ );

  Matrix      Bxmat           ( 1, rank_+1 );
  Matrix      Bymat           ( 1, rank_+1 );
  Matrix      Smat_bot        ( rank_+1, 1 );
  Matrix      Rmat_bot        ( 7, 1 );
  Matrix      Rxmat_bot       ( 7, 1 );
  Matrix      Rymat_bot       ( 7, 1 );
  Matrix      Smat_top        ( rank_+1, 1 );
  Matrix      Rmat_top        ( 7, 1 );
  Matrix      Rxmat_top       ( 7, 1 );
  Matrix      Rymat_top       ( 7, 1 );
  
  Matrix      sfuncs      = shape_->getShapeFunctions ();

  Matrix      Nu_bot      ( 1, 2 * nodeCount_bottop_ );
  Matrix      Nv_bot      ( 1, 2 * nodeCount_bottop_ );
  Matrix      Nw_bot      ( 1, 3 * nodeCount_bottop_ );
  Matrix      Nthetax_bot ( 1, 3 * nodeCount_bottop_ );
  Matrix      Nthetay_bot ( 1, 3 * nodeCount_bottop_ );
  Matrix      Nu_top      ( 1, 2 * nodeCount_bottop_ );
  Matrix      Nv_top      ( 1, 2 * nodeCount_bottop_ );
  Matrix      Nw_top      ( 1, 3 * nodeCount_bottop_ );
  Matrix      Nthetax_top ( 1, 3 * nodeCount_bottop_ );
  Matrix      Nthetay_top ( 1, 3 * nodeCount_bottop_ );

  Matrix      elemMat     ( dofCount, dofCount  );
  Vector      elemForce   ( dofCount  );
  Vector      elemDisp    ( dofCount  );

  IdxVector   Connect      ( 5 * nodeCount_ );
  IdxVector   ConnectInv   ( 5 * nodeCount_ );

  Matrix      elemMat0    ( dofCount, dofCount  );
  Vector      elemForce0  ( dofCount  );
  Vector      elemDisp0   ( dofCount  );

  Vector      jump        ( qRank_     );
  Vector      traction    ( qRank_     );
  Vector      trac        ( qRank_     );

  Matrix      transMat   ( qRank_, qRank_ );
  Matrix      transMatT  = transMat.transpose ();

  IdxVector   inodes     ( nodeCount_ );
  IdxVector   idofs      ( dofCount   );
  Vector      ipWeights  ( ipCount_   );

  Cubix       grads_bot  ( rank_, nodeCount_bottop_, ipCount_bottop_ );
  Cubix       grads_top  ( rank_, nodeCount_bottop_, ipCount_bottop_ );
  Vector      ipWeights_bottop  ( ipCount_bottop_ );

  MChain1     mc1;
  MChain2     mc2;
  MChain3     mc3;

  idx_t       ipoint = 0;
  idx_t       ipoint_bot = 0;
  idx_t       ipoint_top = 0;

  // DOFs connectivity from Jive to the paper from Ai et al. 2024.

  getConnectivity_ ( Connect, ConnectInv );

  // Iterate over all elements assigned to this model.

  for ( idx_t ie = 0; ie < elemCount_; ie++ )
  {
    idx_t  ielem = ielems_[ie];

    // Get the element coordinates and the DOFs.

    elems_.getElemNodes  ( inodes, ielem  );
    nodes_.getSomeCoords ( coords, inodes );
    coords_bot = coords( ALL, slice(0,nodeCount_/2) );
    coords_top = coords( ALL, slice(nodeCount_/2,nodeCount_) );
    dofs_->getDofIndices ( idofs,  inodes, dofTypes_ );
    shape_->getIntegrationWeights ( ipWeights, coords );

    // Get the orientation vectors: v_bot and v_top.

    Vector v_bot ( qRank_ );
    Vector v_top ( qRank_ );
    for ( idx_t i = 0; i < qRank_; i++ )
    {
      v_bot[i] = bot_orientVec_[i];
      v_top[i] = top_orientVec_[i];
    }

    // Get the transformation matrix (rotation tensor) for the structural CE model.

    getTransMatrix_ ( transMat, coords_bot, coords_top, v_bot, v_top );

    // Get the integration point coordinates in the Global system.

    shape_->getGlobalIntegrationPoints(ipcoords, coords);

    // Change of basis to the Local element system.

    getLocalSystem_(ipxcoords_bot, xcoords_bot, ipcoords, coords_bot, v_bot);
    getLocalSystem_(ipxcoords_top, xcoords_top, ipcoords, coords_top, v_top);

    // Get the spatial derivatives of the membrane shape functions in the Local element system.

    bottopshape_->getShapeGradients ( grads_bot, ipWeights_bottop, xcoords_bot );
    bottopshape_->getShapeGradients ( grads_top, ipWeights_bottop, xcoords_top ); 

    // Get the BA, Bx and By matrices, which are equal for both plies.

    getBABxBymats_(BAmat, Bxmat, Bymat);

    // Get the MA, Ma and MAinv matrices for the bottom and top plies.

    getMAMamat_(MAmat_bot, Mamat_bot, xcoords_bot);
    invertMAmat_(MAinv_bot, MAmat_bot);
    getMAMamat_(MAmat_top, Mamat_top, xcoords_top);
    invertMAmat_(MAinv_top, MAmat_top);

    // Get the Tmat matrix.

    getTmat_(Tmat_bot, xcoords_bot);
    getTmat_(Tmat_top, xcoords_bot);

    // Get the displacements at the element nodes, and the bottom and top membrane displacements at the element nodes.

    elemDisp0 = disp[idofs];

    // Get the bottom and top membrane displacements at the element nodes.

    getMembraneDisp_(Disp_mem_bot, Disp_mem_top, elemDisp0);

    // Reordering the DOFs to match the order in the paper from Ai et al. 2024.

    elemDisp = 0.0;
    for ( idx_t i = 0; i < dofCount; i++ )
    {
      elemDisp[i]   = elemDisp0[ConnectInv[i]];
    }

    // Get laminate thickness.

    double h_bot = thickness_[0];
    double h_top = thickness_[1];

    // Assemble the element matrix and force vector.

    elemMat   = 0.0;
    elemForce = 0.0;

    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    {
      // Get the C matrix for the bottom and top plies.

      for ( idx_t ip_botop = 0; ip_botop < ipCount_bottop_; ip_botop++ )
      {
        // Get the plate stiffness matrix for the bottom and top plies.
        
        getShapeGrads_ ( b_mem_bot, grads_bot(ALL,ALL,ip_botop) );
        matmul ( strain_mem_bot, b_mem_bot, Disp_mem_bot );
        material_[0]->update ( stress_mem_bot, stiff_bot, strain_mem_bot, ipoint_bot++ );
        D_plate_bot = stiff_bot*(pow(thickness_[0], 3))/12.0;

        getShapeGrads_ ( b_mem_top, grads_top(ALL,ALL,ip_botop) );
        matmul ( strain_mem_top, b_mem_top, Disp_mem_top );
        material_[1]->update ( stress_mem_top, stiff_top, strain_mem_top, ipoint_top++ );
        D_plate_top = stiff_top*(pow(thickness_[1], 3))/12.0;

        // Get the H matrix for the bottom and top plies.

        getHmatAnalytic_(Hmat_bot, D_plate_bot, xcoords_bot);
        Hinv_bot = Hmat_bot;
        jem::numeric::Cholesky::invert(Hinv_bot);

        getHmatAnalytic_(Hmat_top, D_plate_top, xcoords_top);
        Hinv_top = Hmat_top;
        jem::numeric::Cholesky::invert(Hinv_top);

        // Get the B and C matrices for the bottom and top plies.

        getBmat_(Bmat_bot, D_plate_bot, xcoords_bot);
        Cmat_bot = mc3.matmul ( Hinv_bot, Bmat_bot, Tmat_bot );
        getBmat_(Bmat_top, D_plate_top, xcoords_top);
        Cmat_top = mc3.matmul ( Hinv_top, Bmat_top, Tmat_top );
      }

      // Get the Nu and Nv matrices for the bottom and top plies.

      getNuNvmats_(Nu_bot, Nv_bot, sfuncs, ip);
      getNuNvmats_(Nu_top, Nv_top, sfuncs, ip);

      // Get the integration points for the bottom ply.

      double x, y;

      x = ipxcoords_bot(0,ip);
      y = ipxcoords_bot(1,ip);

      // Get the S, R, Rx and Ry matrices for the bottom ply.

      getSRRxRymats_(Smat_bot, Rmat_bot, Rxmat_bot, Rymat_bot, x, y);

      // Get the integration points for the top ply.

      x = ipxcoords_top(0,ip);
      y = ipxcoords_top(1,ip);

      // Get the S, R, Rx and Ry matrices for the top ply.

      getSRRxRymats_(Smat_top, Rmat_top, Rxmat_top, Rymat_top, x, y);

      // Get the Nw, Nthetax and Nthetay matrices for the bottom and top plies.
        
      BA_MaC_bot = BAmat - matmul(Mamat_bot,Cmat_bot);      
      Nw_bot = matmul ( matmul(Smat_bot.transpose(),MAinv_bot), BA_MaC_bot ) + matmul(Rmat_bot.transpose(),Cmat_bot);
      Nthetax_bot = matmul ( matmul(Bxmat,MAinv_bot), BA_MaC_bot ) + matmul(Rxmat_bot.transpose(),Cmat_bot);
      Nthetay_bot = matmul ( matmul(Bymat,MAinv_bot), BA_MaC_bot ) + matmul(Rymat_bot.transpose(),Cmat_bot);

      BA_MaC_top = BAmat - matmul(Mamat_top,Cmat_top);
      Nw_top = matmul ( matmul(Smat_top.transpose(),MAinv_top), BA_MaC_top ) + matmul(Rmat_top.transpose(),Cmat_top);
      Nthetax_top = matmul ( matmul(Bxmat,MAinv_top), BA_MaC_top ) + matmul(Rxmat_top.transpose(),Cmat_top);
      Nthetay_top = matmul ( matmul(Bymat,MAinv_top), BA_MaC_top ) + matmul(Rymat_top.transpose(),Cmat_top);

      // Assemble the BCE matrix of the structural cohesive element.

      BCEmat =  0.0;
      BCEmat(0,slice(2*nodeCount_bottop_,(2*nodeCount_bottop_+Nw_bot.size(1)))) = -Nw_bot(0,ALL);
      BCEmat(0,slice((4*nodeCount_bottop_+Nw_bot.size(1)),(4*nodeCount_bottop_+2*Nw_bot.size(1)))) = Nw_top(0,ALL);

      BCEmat(1,slice(0,2*nodeCount_bottop_)) = -Nu_bot(0,ALL);
      BCEmat(1,slice(2*nodeCount_bottop_,(2*nodeCount_bottop_+Nw_bot.size(1)))) = (h_bot/2.0)*Nthetax_bot(0,ALL);
      BCEmat(1,slice((2*nodeCount_bottop_+Nw_bot.size(1)),(4*nodeCount_bottop_+Nw_bot.size(1)))) = Nu_top(0,ALL);
      BCEmat(1,slice((4*nodeCount_bottop_+Nw_bot.size(1)),(4*nodeCount_bottop_+2*Nw_bot.size(1)))) = (h_top/2.0)*Nthetax_top(0,ALL);

      BCEmat(2,slice(0,2*nodeCount_bottop_)) = -Nv_bot(0,ALL);
      BCEmat(2,slice(2*nodeCount_bottop_,(2*nodeCount_bottop_+Nw_bot.size(1)))) = (h_bot/2.0)*Nthetay_bot(0,ALL);
      BCEmat(2,slice((2*nodeCount_bottop_+Nw_bot.size(1)),(4*nodeCount_bottop_+Nw_bot.size(1)))) = Nv_top(0,ALL);
      BCEmat(2,slice((4*nodeCount_bottop_+Nw_bot.size(1)),(4*nodeCount_bottop_+2*Nw_bot.size(1)))) = (h_top/2.0)*Nthetay_top(0,ALL);

      // Compute the displacement jump (in local {n,s}-frame)

      jump = mc2.matmul ( transMat, BCEmat, elemDisp );

      // Get the tangent stiffness matrix and the traction

      coheMat_->update ( traction, stiff, jump, ipoint++ );
      
      // transform traction and stiffness to global {x,y}-frame

      matmul ( trac, transMatT, traction );
      stiff    = mc3.matmul ( transMatT, stiff, transMat );

      // Compute the element force vector

      elemForce += ipWeights[ip] * mc1.matmul ( BCEmat.transpose(), trac );

      // Compute the stiffness matrix
      
      elemMat   += ipWeights[ip] * mc3.matmul ( BCEmat.transpose(), stiff, BCEmat );
    }

    // Reordering the elemForce to match the order in the DOfs in Jive.

    elemForce0 = 0.0;
    for ( idx_t i = 0; i < dofCount; i++ )
    {
      elemForce0[i]   = elemForce[Connect[i]];
    }

    // Reordering the elemMat to match the order in the DOfs in Jive.

    elemMat0 = 0.0;
    for ( idx_t i = 0; i < dofCount; i++ )
    {
      for ( idx_t j = 0; j < dofCount; j++ )
      {
        elemMat0(i,j) = elemMat(Connect[i],Connect[j]);
      }
    }

    // Add the element matrix to the global stiffness matrix.

    if ( mbuilder != NIL )
    {
      mbuilder->addBlock ( idofs, idofs, elemMat0 );
    }

    // Add the element force vector to the global force vector.

    force[idofs] += elemForce0;

    if ( jem::Float::isNaN( sum( elemForce0 ) ) ||
        jem::Float::isNaN( sum( elemMat0   ) ) )
    { 
      System::out() << "Interf something's wrong in element " << ie << endl;
    }
  }
}


//-----------------------------------------------------------------------
//   getFrictionForce_
//-----------------------------------------------------------------------


void StructuralInterfaceModel::getFrictionForce_

  ( const Vector&       fsh0,
    const Vector&       disp ) const

{
  const idx_t dofCount   = 5 * nodeCount_;

  Matrix      stiff   ( qRank_, qRank_ );
  Matrix      BCEmat  ( qRank_, dofCount );
  Matrix      QB      ( qRank_, dofCount  );
  Matrix      QBt   = QB.transpose ();

  Matrix      stiff_bot   ( qRank_, qRank_ );
  Matrix      stiff_top   ( qRank_, qRank_ );
  Matrix      D_plate_bot ( qRank_, qRank_ );
  Matrix      D_plate_top ( qRank_, qRank_ );

  Matrix      coords     ( qRank_, nodeCount_ );
              coords     = 0.;
  Matrix      coords_bot  ( qRank_, nodeCount_bottop_ );
  Matrix      coords_top  ( qRank_, nodeCount_bottop_ );
  Matrix      xcoords_bot ( rank_, nodeCount_bottop_ );
  Matrix      xcoords_top ( rank_, nodeCount_bottop_ );

  Matrix      ipcoords    ( qRank_, ipCount_ );
  Matrix      ipxcoords_bot ( rank_, ipCount_ );
  Matrix      ipxcoords_top ( rank_, ipCount_ );

  Matrix      b_mem_bot       ( qRank_, 2*nodeCount_bottop_ );
  Vector      Disp_mem_bot    ( 2*nodeCount_bottop_ );
  Vector      strain_mem_bot  ( qRank_ );
  Vector      stress_mem_bot  ( qRank_ );
  Matrix      b_mem_top       ( qRank_, 2*nodeCount_bottop_ );
  Vector      Disp_mem_top    ( 2*nodeCount_bottop_ );
  Vector      strain_mem_top  ( qRank_ );
  Vector      stress_mem_top  ( qRank_ );

  Matrix      Bmat_bot        ( 7, 4*nodeCount_bottop_);
  Matrix      Tmat_bot        ( 4*nodeCount_bottop_, 3*nodeCount_bottop_ );
  Matrix      Hmat_bot        ( 7, 7);
  Matrix      Hinv_bot        ( 7, 7);
  Matrix      Bmat_top        ( 7, 4*nodeCount_bottop_);
  Matrix      Tmat_top        ( 4*nodeCount_bottop_, 3*nodeCount_bottop_ );
  Matrix      Hmat_top        ( 7, 7);
  Matrix      Hinv_top        ( 7, 7);

  Matrix      BAmat           ( nodeCount_bottop_, 3 * nodeCount_bottop_ );

  Matrix      MAmat_bot       ( nodeCount_bottop_, rank_+1 );
  Matrix      MAinv_bot       ( nodeCount_bottop_, rank_+1 );
  Matrix      Mamat_bot       ( nodeCount_bottop_, 7 );
  Matrix      Cmat_bot        ( 7, 3 * nodeCount_bottop_ );
  Matrix      BA_MaC_bot      ( nodeCount_bottop_, 3 * nodeCount_bottop_ );

  Matrix      MAmat_top       ( nodeCount_bottop_, rank_+1 );
  Matrix      MAinv_top       ( nodeCount_bottop_, rank_+1 );
  Matrix      Mamat_top       ( nodeCount_bottop_, 7 );
  Matrix      Cmat_top        ( 7, 3 * nodeCount_bottop_ );
  Matrix      BA_MaC_top      ( nodeCount_bottop_, 3 * nodeCount_bottop_ );

  Matrix      Bxmat           ( 1, rank_+1 );
  Matrix      Bymat           ( 1, rank_+1 );
  Matrix      Smat_bot        ( rank_+1, 1 );
  Matrix      Rmat_bot        ( 7, 1 );
  Matrix      Rxmat_bot       ( 7, 1 );
  Matrix      Rymat_bot       ( 7, 1 );
  Matrix      Smat_top        ( rank_+1, 1 );
  Matrix      Rmat_top        ( 7, 1 );
  Matrix      Rxmat_top       ( 7, 1 );
  Matrix      Rymat_top       ( 7, 1 );
  
  Matrix      sfuncs      = shape_->getShapeFunctions ();

  Matrix      Nu_bot      ( 1, 2 * nodeCount_bottop_ );
  Matrix      Nv_bot      ( 1, 2 * nodeCount_bottop_ );
  Matrix      Nw_bot      ( 1, 3 * nodeCount_bottop_ );
  Matrix      Nthetax_bot ( 1, 3 * nodeCount_bottop_ );
  Matrix      Nthetay_bot ( 1, 3 * nodeCount_bottop_ );
  Matrix      Nu_top      ( 1, 2 * nodeCount_bottop_ );
  Matrix      Nv_top      ( 1, 2 * nodeCount_bottop_ );
  Matrix      Nw_top      ( 1, 3 * nodeCount_bottop_ );
  Matrix      Nthetax_top ( 1, 3 * nodeCount_bottop_ );
  Matrix      Nthetay_top ( 1, 3 * nodeCount_bottop_ );

  Vector      elemForce  ( dofCount  );
  Vector      elemDisp   ( dofCount  );

  IdxVector   Connect      ( 5 * nodeCount_ );
  IdxVector   ConnectInv   ( 5 * nodeCount_ );

  Vector      elemForce0  ( dofCount  );
  Vector      elemDisp0   ( dofCount  );

  Vector      jump        ( qRank_     );
  Vector      floc        ( qRank_     );
  Vector      traction    ( qRank_     );
  Vector      trac        ( qRank_     );

  Matrix      transMat   ( qRank_, qRank_ );
  Matrix      transMatT  = transMat.transpose ();

  IdxVector   inodes     ( nodeCount_ );
  IdxVector   idofs      ( dofCount   );
  Vector      ipWeights  ( ipCount_   );

  Cubix       grads_bot  ( rank_, nodeCount_bottop_, ipCount_bottop_ );
  Cubix       grads_top  ( rank_, nodeCount_bottop_, ipCount_bottop_ );
  Vector      ipWeights_bottop  ( ipCount_bottop_ );

  MChain1     mc1;
  MChain2     mc2;
  MChain3     mc3;

  idx_t       ipoint = 0;
  idx_t       ipoint_bot = 0;
  idx_t       ipoint_top = 0;

  // DOFs connectivity from Jive to the paper from Ai et al. 2024.

  getConnectivity_ ( Connect, ConnectInv );

  // Iterate over all elements assigned to this model.

  for ( idx_t ie = 0; ie < elemCount_; ie++ )
  {
    idx_t  ielem = ielems_[ie];

    // Get the element coordinates and the DOFs.

    elems_.getElemNodes  ( inodes, ielem  );
    nodes_.getSomeCoords ( coords, inodes );
    coords_bot = coords( ALL, slice(0,nodeCount_/2) );
    coords_top = coords( ALL, slice(nodeCount_/2,nodeCount_) );
    dofs_->getDofIndices ( idofs,  inodes, dofTypes_ );
    shape_->getIntegrationWeights ( ipWeights, coords );

    // Get the orientation vectors: v_bot and v_top.

    Vector v_bot ( qRank_ );
    Vector v_top ( qRank_ );
    for ( idx_t i = 0; i < qRank_; i++ )
    {
      v_bot[i] = bot_orientVec_[i];
      v_top[i] = top_orientVec_[i];
    }

    // Get the transformation matrix (rotation tensor) for the structural CE model.

    getTransMatrix_ ( transMat, coords_bot, coords_top, v_bot, v_top );

    // Get the integration point coordinates in the Global system.

    shape_->getGlobalIntegrationPoints(ipcoords, coords);

    // Change of basis to the Local element system.

    getLocalSystem_(ipxcoords_bot, xcoords_bot, ipcoords, coords_bot, v_bot);
    getLocalSystem_(ipxcoords_top, xcoords_top, ipcoords, coords_top, v_top);

    // Get the spatial derivatives of the membrane shape functions in the Local element system.

    bottopshape_->getShapeGradients ( grads_bot, ipWeights_bottop, xcoords_bot );
    bottopshape_->getShapeGradients ( grads_top, ipWeights_bottop, xcoords_top ); 

    // Get the BA, Bx and By matrices, hich are equal for both plies.

    getBABxBymats_(BAmat, Bxmat, Bymat);

    // Get the MA, Ma and MAinv matrices for the bottom and top plies.

    getMAMamat_(MAmat_bot, Mamat_bot, xcoords_bot);
    invertMAmat_(MAinv_bot, MAmat_bot);
    getMAMamat_(MAmat_top, Mamat_top, xcoords_top);
    invertMAmat_(MAinv_top, MAmat_top);

    // Get the Tmat matrix.

    getTmat_(Tmat_bot, xcoords_bot);
    getTmat_(Tmat_top, xcoords_bot);

    // Get the displacements at the element nodes, and the bottom and top membrane displacements at the element nodes.

    elemDisp0 = disp[idofs];

    // Reordering the DOFs to match the order in the paper from Ai et al. 2024.

    elemDisp = 0.0;
    for ( idx_t i = 0; i < dofCount; i++ )
    {
      elemDisp[i]   = elemDisp0[ConnectInv[i]];
    }

    // Get the bottom and top membrane displacements at the element nodes.

    getMembraneDisp_(Disp_mem_bot, Disp_mem_top, elemDisp0);

    // Get laminate thickness.

    double h_bot = thickness_[0];
    double h_top = thickness_[1];

    // Assemble the element force vector.

    elemForce = 0.0;

    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    {
      // Get the C matrix for the bottom and top plies.

      for ( idx_t ip_botop = 0; ip_botop < ipCount_bottop_; ip_botop++ )
      {
        // Get the plate stiffness matrix for the bottom and top plies.
        
        getShapeGrads_ ( b_mem_bot, grads_bot(ALL,ALL,ip_botop) );
        matmul ( strain_mem_bot, b_mem_bot, Disp_mem_bot );
        material_[0]->update ( stress_mem_bot, stiff_bot, strain_mem_bot, ipoint_bot++ );
        D_plate_bot = stiff_bot*(pow(thickness_[0], 3))/12.0;

        getShapeGrads_ ( b_mem_top, grads_top(ALL,ALL,ip_botop) );
        matmul ( strain_mem_top, b_mem_top, Disp_mem_top );
        material_[1]->update ( stress_mem_top, stiff_top, strain_mem_top, ipoint_top++ );
        D_plate_top = stiff_top*(pow(thickness_[1], 3))/12.0;

        // Get the H matrix for the bottom and top plies.

        getHmatAnalytic_(Hmat_bot, D_plate_bot, xcoords_bot);
        Hinv_bot = Hmat_bot;
        jem::numeric::Cholesky::invert(Hinv_bot);

        getHmatAnalytic_(Hmat_top, D_plate_top, xcoords_top);
        Hinv_top = Hmat_top;
        jem::numeric::Cholesky::invert(Hinv_top);

        // Get the B and C matrices for the bottom and top plies.

        getBmat_(Bmat_bot, D_plate_bot, xcoords_bot);
        Cmat_bot = mc3.matmul ( Hinv_bot, Bmat_bot, Tmat_bot );
        getBmat_(Bmat_top, D_plate_top, xcoords_top);
        Cmat_top = mc3.matmul ( Hinv_top, Bmat_top, Tmat_top );
      }

      // Get the Nu and Nv matrices for the bottom and top plies.

      getNuNvmats_(Nu_bot, Nv_bot, sfuncs, ip);
      getNuNvmats_(Nu_top, Nv_top, sfuncs, ip);

      // Get the integration points for the bottom ply.

      double x, y;

      x = ipxcoords_bot(0,ip);
      y = ipxcoords_bot(1,ip);

      // Get the S, R, Rx and Ry matrices for the bottom ply.

      getSRRxRymats_(Smat_bot, Rmat_bot, Rxmat_bot, Rymat_bot, x, y);

      // Get the integration points for the top ply.

      x = ipxcoords_top(0,ip);
      y = ipxcoords_top(1,ip);

      // Get the S, R, Rx and Ry matrices for the top ply.

      getSRRxRymats_(Smat_top, Rmat_top, Rxmat_top, Rymat_top, x, y);

      // Get the Nw, Nthetax and Nthetay matrices for the bottom and top plies.
        
      BA_MaC_bot = BAmat - matmul(Mamat_bot,Cmat_bot);      
      Nw_bot = matmul ( matmul(Smat_bot.transpose(),MAinv_bot), BA_MaC_bot ) + matmul(Rmat_bot.transpose(),Cmat_bot);
      Nthetax_bot = matmul ( matmul(Bxmat,MAinv_bot), BA_MaC_bot ) + matmul(Rxmat_bot.transpose(),Cmat_bot);
      Nthetay_bot = matmul ( matmul(Bymat,MAinv_bot), BA_MaC_bot ) + matmul(Rymat_bot.transpose(),Cmat_bot);

      BA_MaC_top = BAmat - matmul(Mamat_top,Cmat_top);
      Nw_top = matmul ( matmul(Smat_top.transpose(),MAinv_top), BA_MaC_top ) + matmul(Rmat_top.transpose(),Cmat_top);
      Nthetax_top = matmul ( matmul(Bxmat,MAinv_top), BA_MaC_top ) + matmul(Rxmat_top.transpose(),Cmat_top);
      Nthetay_top = matmul ( matmul(Bymat,MAinv_top), BA_MaC_top ) + matmul(Rymat_top.transpose(),Cmat_top);

      // Assemble the BCE matrix of the structural cohesive element.

      BCEmat =  0.0;
      BCEmat(0,slice(2*nodeCount_bottop_,(2*nodeCount_bottop_+Nw_bot.size(1)))) = -Nw_bot(0,ALL);
      BCEmat(0,slice((4*nodeCount_bottop_+Nw_bot.size(1)),(4*nodeCount_bottop_+2*Nw_bot.size(1)))) = Nw_top(0,ALL);

      BCEmat(1,slice(0,2*nodeCount_bottop_)) = -Nu_bot(0,ALL);
      BCEmat(1,slice(2*nodeCount_bottop_,(2*nodeCount_bottop_+Nw_bot.size(1)))) = (h_bot/2.0)*Nthetax_bot(0,ALL);
      BCEmat(1,slice((2*nodeCount_bottop_+Nw_bot.size(1)),(4*nodeCount_bottop_+Nw_bot.size(1)))) = Nu_top(0,ALL);
      BCEmat(1,slice((4*nodeCount_bottop_+Nw_bot.size(1)),(4*nodeCount_bottop_+2*Nw_bot.size(1)))) = (h_top/2.0)*Nthetax_top(0,ALL);

      BCEmat(2,slice(0,2*nodeCount_bottop_)) = -Nv_bot(0,ALL);
      BCEmat(2,slice(2*nodeCount_bottop_,(2*nodeCount_bottop_+Nw_bot.size(1)))) = (h_bot/2.0)*Nthetay_bot(0,ALL);
      BCEmat(2,slice((2*nodeCount_bottop_+Nw_bot.size(1)),(4*nodeCount_bottop_+Nw_bot.size(1)))) = Nv_top(0,ALL);
      BCEmat(2,slice((4*nodeCount_bottop_+Nw_bot.size(1)),(4*nodeCount_bottop_+2*Nw_bot.size(1)))) = (h_top/2.0)*Nthetay_top(0,ALL);

      // Compute the displacement jump (in local {n,s}-frame)

      matmul ( QB, transMat, BCEmat );
      matmul ( jump, QB, elemDisp );

      // Get the local contribution

      frictionMat_->getDissForce ( floc, jump, ipoint++ );

      // Compute the element force vector

      elemForce += ipWeights[ip] * mc1.matmul ( QBt, floc );
    }

    // Reordering the elemForce to match the order in the DOfs in Jive.

    elemForce0 = 0.0;
    for ( idx_t i = 0; i < dofCount; i++ )
    {
      elemForce0[i]   = elemForce[Connect[i]];
    }

    // Add the element force vector to the global force vector.

    fsh0[idofs] += elemForce0;
  }
}


//-----------------------------------------------------------------------
//   writeOutput_
//-----------------------------------------------------------------------


void StructuralInterfaceModel::writeOutput_

  ( const Properties&    globdat  ) const

{
  const idx_t dofCount   = 5 * nodeCount_;

  Matrix      stiff   ( qRank_, qRank_ );
  Matrix      BCEmat  ( qRank_, dofCount );

  Matrix      stiff_bot   ( qRank_, qRank_ );
  Matrix      stiff_top   ( qRank_, qRank_ );
  Matrix      D_plate_bot ( qRank_, qRank_ );
  Matrix      D_plate_top ( qRank_, qRank_ );

  Matrix      coords      ( qRank_, nodeCount_ );
              coords     = 0.;
  Matrix      coords_bot  ( qRank_, nodeCount_bottop_ );
  Matrix      coords_top  ( qRank_, nodeCount_bottop_ );
  Matrix      xcoords_bot ( rank_, nodeCount_bottop_ );
  Matrix      xcoords_top ( rank_, nodeCount_bottop_ );

  Matrix      ipcoords    ( qRank_, ipCount_ );
  Matrix      ipxcoords_bot ( rank_, ipCount_ );
  Matrix      ipxcoords_top ( rank_, ipCount_ );

  Matrix      b_mem_bot       ( qRank_, 2*nodeCount_bottop_ );
  Vector      Disp_mem_bot    ( 2*nodeCount_bottop_ );
  Vector      strain_mem_bot  ( qRank_ );
  Vector      stress_mem_bot  ( qRank_ );
  Matrix      b_mem_top       ( qRank_, 2*nodeCount_bottop_ );
  Vector      Disp_mem_top    ( 2*nodeCount_bottop_ );
  Vector      strain_mem_top  ( qRank_ );
  Vector      stress_mem_top  ( qRank_ );

  Matrix      Bmat_bot        ( 7, 4*nodeCount_bottop_);
  Matrix      Tmat_bot        ( 4*nodeCount_bottop_, 3*nodeCount_bottop_ );
  Matrix      Hmat_bot        ( 7, 7);
  Matrix      Hinv_bot        ( 7, 7);
  Matrix      Bmat_top        ( 7, 4*nodeCount_bottop_);
  Matrix      Tmat_top        ( 4*nodeCount_bottop_, 3*nodeCount_bottop_ );
  Matrix      Hmat_top        ( 7, 7);
  Matrix      Hinv_top        ( 7, 7);

  Matrix      BAmat           ( nodeCount_bottop_, 3 * nodeCount_bottop_ );

  Matrix      MAmat_bot       ( nodeCount_bottop_, rank_+1 );
  Matrix      MAinv_bot       ( nodeCount_bottop_, rank_+1 );
  Matrix      Mamat_bot       ( nodeCount_bottop_, 7 );
  Matrix      Cmat_bot        ( 7, 3 * nodeCount_bottop_ );
  Matrix      BA_MaC_bot      ( nodeCount_bottop_, 3 * nodeCount_bottop_ );

  Matrix      MAmat_top       ( nodeCount_bottop_, rank_+1 );
  Matrix      MAinv_top       ( nodeCount_bottop_, rank_+1 );
  Matrix      Mamat_top       ( nodeCount_bottop_, 7 );
  Matrix      Cmat_top        ( 7, 3 * nodeCount_bottop_ );
  Matrix      BA_MaC_top      ( nodeCount_bottop_, 3 * nodeCount_bottop_ );

  Matrix      Bxmat           ( 1, rank_+1 );
  Matrix      Bymat           ( 1, rank_+1 );
  Matrix      Smat_bot        ( rank_+1, 1 );
  Matrix      Rmat_bot        ( 7, 1 );
  Matrix      Rxmat_bot       ( 7, 1 );
  Matrix      Rymat_bot       ( 7, 1 );
  Matrix      Smat_top        ( rank_+1, 1 );
  Matrix      Rmat_top        ( 7, 1 );
  Matrix      Rxmat_top       ( 7, 1 );
  Matrix      Rymat_top       ( 7, 1 );
  
  Matrix      sfuncs      = shape_->getShapeFunctions ();

  Matrix      Nu_bot      ( 1, 2 * nodeCount_bottop_ );
  Matrix      Nv_bot      ( 1, 2 * nodeCount_bottop_ );
  Matrix      Nw_bot      ( 1, 3 * nodeCount_bottop_ );
  Matrix      Nthetax_bot ( 1, 3 * nodeCount_bottop_ );
  Matrix      Nthetay_bot ( 1, 3 * nodeCount_bottop_ );
  Matrix      Nu_top      ( 1, 2 * nodeCount_bottop_ );
  Matrix      Nv_top      ( 1, 2 * nodeCount_bottop_ );
  Matrix      Nw_top      ( 1, 3 * nodeCount_bottop_ );
  Matrix      Nthetax_top ( 1, 3 * nodeCount_bottop_ );
  Matrix      Nthetay_top ( 1, 3 * nodeCount_bottop_ );

  Vector      elemDisp    ( dofCount  );

  IdxVector   Connect      ( 5 * nodeCount_ );
  IdxVector   ConnectInv   ( 5 * nodeCount_ );

  Vector      elemDisp0   ( dofCount  );

  Vector      jump        ( qRank_     );
  Vector      traction    ( qRank_     );
  Vector      trac        ( qRank_     );

  Matrix      transMat   ( qRank_, qRank_ );
  Matrix      transMatT  = transMat.transpose ();

  IdxVector   inodes     ( nodeCount_ );
  IdxVector   idofs      ( dofCount   );
  Vector      ipWeights  ( ipCount_   );

  Cubix       grads_bot  ( rank_, nodeCount_bottop_, ipCount_bottop_ );
  Cubix       grads_top  ( rank_, nodeCount_bottop_, ipCount_bottop_ );
  Vector      ipWeights_bottop  ( ipCount_bottop_ );

  MChain2     mc2;
  MChain3     mc3;

  Vector      disp;

  StateVector::get ( disp, dofs_, globdat );

  idx_t       ipoint = 0;
  idx_t       ipoint_bot = 0;
  idx_t       ipoint_top = 0;

  // DOFs connectivity from Jive to the paper from Ai et al. 2024.

  getConnectivity_ ( Connect, ConnectInv );

  idx_t       it;

  globdat.get ( it, Globdat::TIME_STEP );

  *xOut_ << "newXOutput " << it <<  '\n';

  for ( idx_t ie = 0; ie < elemCount_; ie++ )
  {
    idx_t  ielem = ielems_[ie];

    // Get the element coordinates and the DOFs.

    elems_.getElemNodes  ( inodes, ielem  );
    nodes_.getSomeCoords ( coords, inodes );
    coords_bot = coords( ALL, slice(0,nodeCount_/2) );
    coords_top = coords( ALL, slice(nodeCount_/2,nodeCount_) );
    dofs_->getDofIndices ( idofs,  inodes, dofTypes_ );
    shape_->getIntegrationWeights ( ipWeights, coords );

    // Get the orientation vectors: v_bot and v_top.

    Vector v_bot ( qRank_ );
    Vector v_top ( qRank_ );
    for ( idx_t i = 0; i < qRank_; i++ )
    {
      v_bot[i] = bot_orientVec_[i];
      v_top[i] = top_orientVec_[i];
    }

    // Get the transformation matrix (rotation tensor) for the structural CE model.

    getTransMatrix_ ( transMat, coords_bot, coords_top, v_bot, v_top );

    // Get the integration point coordinates in the Global system.

    shape_->getGlobalIntegrationPoints(ipcoords, coords);

    // Change of basis to the Local element system.

    getLocalSystem_(ipxcoords_bot, xcoords_bot, ipcoords, coords_bot, v_bot);
    getLocalSystem_(ipxcoords_top, xcoords_top, ipcoords, coords_top, v_top);

    // Get the spatial derivatives of the membrane shape functions in the Local element system.

    bottopshape_->getShapeGradients ( grads_bot, ipWeights_bottop, xcoords_bot );
    bottopshape_->getShapeGradients ( grads_top, ipWeights_bottop, xcoords_top ); 

    // Get the BA, Bx and By matrices, hich are equal for both plies.

    getBABxBymats_(BAmat, Bxmat, Bymat);

    // Get the MA, Ma and MAinv matrices for the bottom and top plies.

    getMAMamat_(MAmat_bot, Mamat_bot, xcoords_bot);
    invertMAmat_(MAinv_bot, MAmat_bot);
    getMAMamat_(MAmat_top, Mamat_top, xcoords_top);
    invertMAmat_(MAinv_top, MAmat_top);

    // Get the Tmat matrix.

    getTmat_(Tmat_bot, xcoords_bot);
    getTmat_(Tmat_top, xcoords_bot);

    // Get the displacements at the element nodes, and the bottom and top membrane displacements at the element nodes.

    elemDisp0 = disp[idofs];

    // Reordering the DOFs to match the order in the paper from Ai et al. 2024.

    elemDisp = 0.0;
    for ( idx_t i = 0; i < dofCount; i++ )
    {
      elemDisp[i]   = elemDisp0[ConnectInv[i]];
    }

    // Get the bottom and top membrane displacements at the element nodes.

    getMembraneDisp_(Disp_mem_bot, Disp_mem_top, elemDisp0);

    // Get laminate thickness.

    double h_bot = thickness_[0];
    double h_top = thickness_[1];

    // Evaluate and write jump and traction

    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    {
      // Get the C matrix for the bottom and top plies.

      for ( idx_t ip_botop = 0; ip_botop < ipCount_bottop_; ip_botop++ )
      {
        // Get the plate stiffness matrix for the bottom and top plies.
        
        getShapeGrads_ ( b_mem_bot, grads_bot(ALL,ALL,ip_botop) );
        matmul ( strain_mem_bot, b_mem_bot, Disp_mem_bot );
        material_[0]->update ( stress_mem_bot, stiff_bot, strain_mem_bot, ipoint_bot++ );
        D_plate_bot = stiff_bot*(pow(thickness_[0], 3))/12.0;

        getShapeGrads_ ( b_mem_top, grads_top(ALL,ALL,ip_botop) );
        matmul ( strain_mem_top, b_mem_top, Disp_mem_top );
        material_[1]->update ( stress_mem_top, stiff_top, strain_mem_top, ipoint_top++ );
        D_plate_top = stiff_top*(pow(thickness_[1], 3))/12.0;

        // Get the H matrix for the bottom and top plies.

        getHmatAnalytic_(Hmat_bot, D_plate_bot, xcoords_bot);
        Hinv_bot = Hmat_bot;
        jem::numeric::Cholesky::invert(Hinv_bot);

        getHmatAnalytic_(Hmat_top, D_plate_top, xcoords_top);
        Hinv_top = Hmat_top;
        jem::numeric::Cholesky::invert(Hinv_top);

        // Get the B and C matrices for the bottom and top plies.

        getBmat_(Bmat_bot, D_plate_bot, xcoords_bot);
        Cmat_bot = mc3.matmul ( Hinv_bot, Bmat_bot, Tmat_bot );
        getBmat_(Bmat_top, D_plate_top, xcoords_top);
        Cmat_top = mc3.matmul ( Hinv_top, Bmat_top, Tmat_top );
      }

      // Get the Nu and Nv matrices for the bottom and top plies.

      getNuNvmats_(Nu_bot, Nv_bot, sfuncs, ip);
      getNuNvmats_(Nu_top, Nv_top, sfuncs, ip);

      // Get the integration points for the bottom ply.

      double x, y;

      x = ipxcoords_bot(0,ip);
      y = ipxcoords_bot(1,ip);

      // Get the S, R, Rx and Ry matrices for the bottom ply.

      getSRRxRymats_(Smat_bot, Rmat_bot, Rxmat_bot, Rymat_bot, x, y);

      // Get the integration points for the top ply.

      x = ipxcoords_top(0,ip);
      y = ipxcoords_top(1,ip);

      // Get the S, R, Rx and Ry matrices for the top ply.

      getSRRxRymats_(Smat_top, Rmat_top, Rxmat_top, Rymat_top, x, y);

      // Get the Nw, Nthetax and Nthetay matrices for the bottom and top plies.
        
      BA_MaC_bot = BAmat - matmul(Mamat_bot,Cmat_bot);      
      Nw_bot = matmul ( matmul(Smat_bot.transpose(),MAinv_bot), BA_MaC_bot ) + matmul(Rmat_bot.transpose(),Cmat_bot);
      Nthetax_bot = matmul ( matmul(Bxmat,MAinv_bot), BA_MaC_bot ) + matmul(Rxmat_bot.transpose(),Cmat_bot);
      Nthetay_bot = matmul ( matmul(Bymat,MAinv_bot), BA_MaC_bot ) + matmul(Rymat_bot.transpose(),Cmat_bot);

      BA_MaC_top = BAmat - matmul(Mamat_top,Cmat_top);
      Nw_top = matmul ( matmul(Smat_top.transpose(),MAinv_top), BA_MaC_top ) + matmul(Rmat_top.transpose(),Cmat_top);
      Nthetax_top = matmul ( matmul(Bxmat,MAinv_top), BA_MaC_top ) + matmul(Rxmat_top.transpose(),Cmat_top);
      Nthetay_top = matmul ( matmul(Bymat,MAinv_top), BA_MaC_top ) + matmul(Rymat_top.transpose(),Cmat_top);

      // Assemble the BCE matrix of the structural cohesive element.

      BCEmat =  0.0;
      BCEmat(0,slice(2*nodeCount_bottop_,(2*nodeCount_bottop_+Nw_bot.size(1)))) = -Nw_bot(0,ALL);
      BCEmat(0,slice((4*nodeCount_bottop_+Nw_bot.size(1)),(4*nodeCount_bottop_+2*Nw_bot.size(1)))) = Nw_top(0,ALL);

      BCEmat(1,slice(0,2*nodeCount_bottop_)) = -Nu_bot(0,ALL);
      BCEmat(1,slice(2*nodeCount_bottop_,(2*nodeCount_bottop_+Nw_bot.size(1)))) = (h_bot/2.0)*Nthetax_bot(0,ALL);
      BCEmat(1,slice((2*nodeCount_bottop_+Nw_bot.size(1)),(4*nodeCount_bottop_+Nw_bot.size(1)))) = Nu_top(0,ALL);
      BCEmat(1,slice((4*nodeCount_bottop_+Nw_bot.size(1)),(4*nodeCount_bottop_+2*Nw_bot.size(1)))) = (h_top/2.0)*Nthetax_top(0,ALL);

      BCEmat(2,slice(0,2*nodeCount_bottop_)) = -Nv_bot(0,ALL);
      BCEmat(2,slice(2*nodeCount_bottop_,(2*nodeCount_bottop_+Nw_bot.size(1)))) = (h_bot/2.0)*Nthetay_bot(0,ALL);
      BCEmat(2,slice((2*nodeCount_bottop_+Nw_bot.size(1)),(4*nodeCount_bottop_+Nw_bot.size(1)))) = Nv_top(0,ALL);
      BCEmat(2,slice((4*nodeCount_bottop_+Nw_bot.size(1)),(4*nodeCount_bottop_+2*Nw_bot.size(1)))) = (h_top/2.0)*Nthetay_top(0,ALL);

      // Compute the displacement jump (in local {n,s}-frame)

      jump = mc2.matmul ( transMat, BCEmat, elemDisp );

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
  xOut_->flush();
}


//-----------------------------------------------------------------------
//   getDissipation_
//-----------------------------------------------------------------------


void StructuralInterfaceModel::getDissipation_

  ( const Properties& params ) const

{
  Matrix      coords     ( qRank_, nodeCount_ );
              coords     = 0.;
  Matrix      ipcoords   ( 2, nodeCount_ );
  Vector      ipWeights  (           ipCount_ );
  IdxVector   inodes     (         nodeCount_ );

  double dissipation = 0.;
  idx_t       ipoint = 0;
  double thickness   = 0.;

  // loop over integration points

  for ( idx_t ie = 0; ie < elemCount_; ie++ )
  {
    // get the integrationpoint weights

    elems_.getElemNodes  ( inodes, ielems_[ie] );
    nodes_.getSomeCoords ( coords( slice(0,qRank_) , ALL ), inodes );
    shape_->getIntegrationWeights ( ipWeights, coords );

    // get dissipation

    for ( idx_t ip = 0; ip < ipCount_; ip++ )
    {

      double    G  = coheMat_->giveDissipation ( ipoint++ );

      dissipation += thickness * ipWeights[ ip ] * G;
    }
  }
  params.set ( myTag_, dissipation );
}


//-----------------------------------------------------------------------
//   getConnectivity_
//-----------------------------------------------------------------------


void StructuralInterfaceModel::getConnectivity_

  ( IdxVector&      Connect,
    IdxVector&      ConnectInv )   const
    
{
  // DOFs connectivity from Jive to the paper from Ai et al. 2024 - Connect.

  for ( idx_t i = 0; i < 2; i++ )
  {
    Connect[15*i] = 15*i;
    Connect[15*i+1] = 15*i+1;
    Connect[15*i+2] = 15*i+6;
    Connect[15*i+3] = 15*i+7;
    Connect[15*i+4] = 15*i+8;
    Connect[15*i+5] = 15*i+2;
    Connect[15*i+6] = 15*i+3;
    Connect[15*i+7] = 15*i+9;
    Connect[15*i+8] = 15*i+10;
    Connect[15*i+9] = 15*i+11;
    Connect[15*i+10] = 15*i+4;
    Connect[15*i+11] = 15*i+5;
    Connect[15*i+12] = 15*i+12;
    Connect[15*i+13] = 15*i+13;
    Connect[15*i+14] = 15*i+14;
  }

  // DOFs connectivity from the paper Ai et al. 2024 to Jive - ConnectInv.
  
  for ( idx_t i = 0; i < 2; i++ )
  {
    ConnectInv[15*i] = 15*i;
    ConnectInv[15*i+1] = 15*i+1;
    ConnectInv[15*i+2] = 15*i+5;
    ConnectInv[15*i+3] = 15*i+6;
    ConnectInv[15*i+4] = 15*i+10;
    ConnectInv[15*i+5] = 15*i+11;
    ConnectInv[15*i+6] = 15*i+2;
    ConnectInv[15*i+7] = 15*i+3;
    ConnectInv[15*i+8] = 15*i+4;
    ConnectInv[15*i+9] = 15*i+7;
    ConnectInv[15*i+10] = 15*i+8;
    ConnectInv[15*i+11] = 15*i+9;
    ConnectInv[15*i+12] = 15*i+12;
    ConnectInv[15*i+13] = 15*i+13;
    ConnectInv[15*i+14] = 15*i+14;
  }
}


//-----------------------------------------------------------------------
//   getTransMatrix_
//-----------------------------------------------------------------------


void StructuralInterfaceModel::getTransMatrix_
   
  ( Matrix&       transMat,
    const Matrix&       coords_bot,
    const Matrix&       coords_top,
    const Vector&       v_bot,
    const Vector&       v_top   ) const

{
  using jem::ALL;
  using jem::numeric::matmul;

  Matrix coords ( qRank_, nodeCount_bottop_ );
  Vector v ( qRank_ );

  coords = (coords_bot + coords_top) / 2.0;
  v = (v_bot + v_top) / 2.0;

  Vector      e1       ( 3 );
  Vector      e2       ( 3 );
  Vector      e3       ( 3 );
  Vector      e1_l     ( 3 );
  Vector      e2_l     ( 3 );
  Vector      e3_l     ( 3 );

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
//   getLocalSystem_
//-----------------------------------------------------------------------


void StructuralInterfaceModel::getLocalSystem_

  ( Matrix&  ipxcoords,
    Matrix&  xcoords,
    const Matrix&   ipcoords,
    const Matrix&   coords,
    const Vector&        v )   const
    
{
  using jem::ALL;
  using jem::numeric::matmul;

  Vector      e1       ( qRank_ );
  Vector      e2       ( qRank_ );
  Vector      e3       ( qRank_ );
  Vector      e1_l     ( qRank_ );
  Vector      e2_l     ( qRank_ );
  Vector      e3_l     ( qRank_ );
  Matrix      Transmat ( rank_, rank_ );

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

  Transmat(0,0) = e1_l[0]*e1[0] + e1_l[1]*e1[1] + e1_l[2]*e1[2];
  Transmat(0,1) = e1_l[0]*e2[0] + e1_l[1]*e2[1] + e1_l[2]*e2[2];
  Transmat(1,0) = e2_l[0]*e1[0] + e2_l[1]*e1[1] + e2_l[2]*e1[2];
  Transmat(1,1) = e2_l[0]*e2[0] + e2_l[1]*e2[1] + e2_l[2]*e2[2];

  // Get the local coordinates of the element nodes and integration points.

  double Xcg = (coords(0, 0) + coords(0, 1) + coords(0, 2))/3.0;
  double Ycg = (coords(1, 0) + coords(1, 1) + coords(1, 2))/3.0;

  // Translation to the baricenter of the element.

  xcoords(0, ALL) = coords(0, ALL) - Xcg;
  xcoords(1, ALL) = coords(1, ALL) - Ycg;

  ipxcoords(0, ALL) = ipcoords(0, ALL) - Xcg;
  ipxcoords(1, ALL) = ipcoords(1, ALL) - Ycg;

  // Apply the transformation to the xcoords.

  for ( idx_t i = 0; i < xcoords.size(1); i++ )
  {
    xcoords(ALL, i) = matmul ( Transmat, xcoords(ALL, i) );
  }

  // Apply the transformation to the ipxcoords.

  for ( idx_t i = 0; i < ipxcoords.size(1); i++ )
  {
    ipxcoords(ALL, i) = matmul ( Transmat, ipxcoords(ALL, i) );
  }
}


//-----------------------------------------------------------------------
//   getBABxBymats_
//-----------------------------------------------------------------------


void StructuralInterfaceModel::getBABxBymats_

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
//   getMAMamat_
//-----------------------------------------------------------------------


void StructuralInterfaceModel::getMAMamat_

  ( Matrix&  MAmat,
    Matrix&  Mamat,
    const Matrix&   xcoords ) const

{
  using jem::ALL;
  Vector x ( nodeCount_bottop_);
  Vector y ( nodeCount_bottop_);

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


void StructuralInterfaceModel::invertMAmat_

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


void StructuralInterfaceModel::getTmat_

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


void StructuralInterfaceModel::getMembraneDisp_

  ( Vector&         Disp_mem_bot,
    Vector&         Disp_mem_top,
    const Vector&   elemDisp0 )   const
    
{
  for ( idx_t i = 0; i < nodeCount_; i++ )
  {
    if ( i < (nodeCount_/2))
    {
      Disp_mem_bot[2*i]   = elemDisp0[5*i];
      Disp_mem_bot[2*i+1] = elemDisp0[5*i+1];
    }
    else
    {
      idx_t ii = i - (nodeCount_/2);
      Disp_mem_top[2*ii]   = elemDisp0[5*i];
      Disp_mem_top[2*ii+1] = elemDisp0[5*i+1];
    }
  }
}


//-----------------------------------------------------------------------
//   getHmatAnalytic_
//-----------------------------------------------------------------------


void StructuralInterfaceModel::getHmatAnalytic_

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


void StructuralInterfaceModel::getBmat_

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
//   getNuNvmats_
//-----------------------------------------------------------------------


void StructuralInterfaceModel::getNuNvmats_

  ( Matrix&           Nu,
    Matrix&           Nv,
    const Matrix&     sfuncs,
    const idx_t&      ip ) const
{
  Nu = 0.0;
  Nv = 0.0;
  for ( idx_t i = 0; i < nodeCount_bottop_; i++ )
  {
    Nu(0,2*i)   = sfuncs(i,ip);

    Nv(0,2*i+1) = sfuncs(i,ip);
  }
}


//-----------------------------------------------------------------------
//   getSRRxRymats_
//-----------------------------------------------------------------------


void StructuralInterfaceModel::getSRRxRymats_

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


void StructuralInterfaceModel::initWriter_  

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


void StructuralInterfaceModel::writeGeom_  ()  const

{
  bool        gauss1     = false;
  bool        gauss3     = false;

  // write nodal coordinates when 1pt integration

  if ( ipCount_ == 1 )
  {
    gauss1 = true;
  }

  Matrix      coords     ( qRank_, nodeCount_ );
              coords     = 0.;

  Matrix      ipCoords   ( qRank_, ipCount_ );

  IdxVector   inodes     ( nodeCount_ );

  *xOut_ << "ipCoords" << '\n';
  *xOut_ << rank_ << '\n';

  if ( ( nodeCount_ == 6 || nodeCount_ == 12 ) && ipCount_ == 3 
       && shape_->getBShape()->getName() != "BTriangleNC3" )
  {
    gauss3 = true;
    System::warn() << myName_ << " Writing gauss point data instead of " 
      << "nodal values." << endl;
  }

  for ( idx_t ie = 0; ie < elemCount_; ie++ )
  {
    // Get the element coordinates and DOFs.

    elems_.getElemNodes  ( inodes, ielems_[ie] );
    nodes_.getSomeCoords ( coords( slice(0,qRank_) , ALL ), 
                           inodes );

    // take node coordinates in case of triangle with Gauss3 integration

    if ( gauss3 || gauss1 )
    {
      idx_t i = nodeCount_ == 3 ? 1 : 2;
      ipCoords.ref ( coords ( ALL, slice(0,nodeCount_/2,i) ) );
    }
    else
    {
      shape_->getGlobalIntegrationPoints ( ipCoords, coords );
    }

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


//-----------------------------------------------------------------------
//   checkGeom_
//-----------------------------------------------------------------------


void StructuralInterfaceModel::checkGeom_ () const

{
  Matrix      normals    ( qRank_, nodeCount_ / 2 );
  Matrix      coords     ( qRank_, nodeCount_ );
  IdxVector   inodes     ( nodeCount_ );

  for ( idx_t ie = 0; ie < elemCount_; ie++ )
  {
    elems_.getElemNodes  ( inodes, ielems_[ie] );
    nodes_.getSomeCoords ( coords( slice(0,qRank_) , ALL ), 
                           inodes );

    // check assumption used in getTransMatrix_ 

    if ( qRank_ == 3 )
    {
      shape_->getVertexNormals ( normals, coords );
      
      for ( idx_t in = 0; in < nodeCount_/2; ++in )
      {
        if ( normals ( 2, in ) < 1.-1.e-10 )
        {
          System::warn() << "Interface element normal not in positive"
            "z-direction! " << ielems_[ie] << endl;
        }
      }
    }
  }
}


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newStructuralInterfaceModel
//-----------------------------------------------------------------------


static Ref<Model>     newStructuralInterfaceModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<StructuralInterfaceModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareStructuralInterfaceModel
//-----------------------------------------------------------------------


void declareStructuralInterfaceModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "StructuralInterface", & newStructuralInterfaceModel );
}


//-----------------------------------------------------------------------
//   makeBTriangleNC
//-----------------------------------------------------------------------


Ref<BShape> makeBTriangleNC ( const Properties& props )
{
  using jive::Matrix;
  using jive::geom::BoundaryTriangle3;
  using jive::geom::BoundaryTriangle6;

  String     type;
  String     intScheme;
  Matrix     ischeme;

  props.find (      type,      "type" );
  props.find ( intScheme, "intScheme" );

  // Set up the Newton-Cotes scheme for the interior of the
  // triangle. The first matrix row contains the weights and the other
  // two rows contain the coordinates. Both the weights and
  // coordinates must be specified in the local coordinate system of a
  // standard triangle. The weights must add up to 0.5 (the area of a
  // standard triangle).

  if ( intScheme == "NewtonCotes3" )
  {
    ischeme.resize ( 3, 3 );

    ischeme(0,ALL) = 1.0 / 6.0;

    ischeme(1,0) = 0.0;
    ischeme(2,0) = 0.0;

    ischeme(1,1) = 1.0;
    ischeme(2,1) = 0.0;

    ischeme(1,2) = 0.0;
    ischeme(2,2) = 1.0;
  }
  else if ( intScheme == "NewtonCotes6" )
  {
    ischeme.resize ( 3, 6 );

    ischeme(0,slice(0,END,2)) = 1.0 / 24.0;
    ischeme(0,slice(1,END,2)) = 3.0 / 24.0;

    ischeme(1,0) = 0.0;
    ischeme(2,0) = 0.0;

    ischeme(1,1) = 0.5;
    ischeme(2,1) = 0.0;

    ischeme(1,2) = 1.0;
    ischeme(2,2) = 0.0;

    ischeme(1,3) = 0.5;
    ischeme(2,3) = 0.5;

    ischeme(1,4) = 0.0;
    ischeme(2,4) = 1.0;

    ischeme(1,5) = 0.0;
    ischeme(2,5) = 0.5;
  }
  else if ( intScheme == "NewtonCotes13" )
  {
    double area = 0.5;
    ischeme.resize ( 3, 13 );

    ischeme(0,0) = -0.149570044467670*area;
    ischeme(0,slice(1,4)) = 0.175615257433204*area;
    ischeme(0,slice(4,7)) = 0.053347235608839*area;
    ischeme(0,slice(7,END)) = 0.077113760890257*area;

    ischeme(1,0) = 1.0 / 3.0;
    ischeme(2,0) = 1.0 / 3.0;

    ischeme(1,1) = 0.479308067841923;
    ischeme(2,1) = 0.260345966079038;

    ischeme(1,2) = 0.260345966079038;
    ischeme(2,2) = 0.479308067841923;

    ischeme(1,3) = 0.260345966079038;
    ischeme(2,3) = 0.260345966079038;

    ischeme(1,4) = 0.869739794195568;
    ischeme(2,4) = 0.065130102902216;

    ischeme(1,5) = 0.065130102902216;
    ischeme(2,5) = 0.869739794195568;

    ischeme(1,6) = 0.065130102902216;
    ischeme(2,6) = 0.065130102902216;

    ischeme(1,7) = 0.638444188569809;
    ischeme(2,7) = 0.312865496004875;

    ischeme(1,8) = 0.638444188569809;
    ischeme(2,8) = 0.048690315425316;

    ischeme(1,9) = 0.312865496004875;
    ischeme(2,9) = 0.638444188569809;

    ischeme(1,10) = 0.312865496004875;
    ischeme(2,10) = 0.048690315425316;

    ischeme(1,11) = 0.048690315425316;
    ischeme(2,11) = 0.638444188569809;

    ischeme(1,12) = 0.048690315425316;
    ischeme(2,12) = 0.312865496004875;
  }
  else
  {
    return NIL;
  }

  if ( type == "BTriangle3" ) 
  {
    return BoundaryTriangle3::getShape ( "BTriangleNC3", ischeme );
  }
  else if ( type == "BTriangle6" )
  {
    return BoundaryTriangle6::getShape ( "BTriangleNC6", ischeme );
  }

  return NIL;
}
