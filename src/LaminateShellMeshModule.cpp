#include <jem/base/array/utilities.h>
#include <jem/base/array/operators.h>
#include <jem/base/array/select.h>
#include <jem/base/Exception.h>
#include <jem/base/IllegalInputException.h>
#include <jem/base/System.h>
#include <jem/util/ArrayBuffer.h>
#include <jem/util/Properties.h>
#include <jem/util/Dictionary.h>
#include <jive/app/ModuleFactory.h>
#include <jive/fem/NodeGroup.h>
#include <jive/util/DofSpace.h>
#include <jive/util/Globdat.h>
#include <jive/util/XItemGroup.h>
#include <jive/util/ItemSet.h>
#include <jive/util/ConstraintsParser.h>

#include "LaminateShellMeshModule.h"

using jem::dynamicCast;
using jem::IllegalInputException;
using jem::io::endl;
using jem::util::ArrayBuffer;
using jem::util::Dict;
using jem::util::DictEnum;
using jive::util::Globdat;
using jive::IdxVector;
using jive::IdxMatrix;
using jive::Matrix;
using jive::fem::toXElementSet;
using jive::fem::toXNodeSet;
using jive::fem::newXNodeSet;
using jive::fem::NodeGroup;
using jive::util::Constraints;
using jive::util::DofSpace;
using jive::util::ItemGroup;
using jive::util::XItemGroup;
using jive::util::ItemSet;


//=======================================================================
//   class LaminateShellMeshModule
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char* LaminateShellMeshModule::DOF_NAMES[5]    = {"u","v","w","wx","wy"};
const char* LaminateShellMeshModule::TYPE_NAME       = "LaminateShellMesh";
const char* LaminateShellMeshModule::DIM             = "dim";
const char* LaminateShellMeshModule::N_LAYERS        = "nLayers";
const char* LaminateShellMeshModule::NOTCHED_INTERF  = "notchedInterface";
const char* LaminateShellMeshModule::THICKNESS       = "thickness";
const char* LaminateShellMeshModule::INTERFACES      = "interfaces";
const char* LaminateShellMeshModule::LAYER_NAMES     = "layerNames";
const char* LaminateShellMeshModule::INTERFACE_NAMES = "interfaceNames";
const char* LaminateShellMeshModule::NODE_GROUPS     = "nodeGroups";
const char* LaminateShellMeshModule::TIE_GROUPS      = "tieGroups";
const char* LaminateShellMeshModule::SKIP_NGROUPS    = "skipNGroups";
const char* LaminateShellMeshModule::SKIP_EGROUPS    = "skipEGroups";
const char* LaminateShellMeshModule::NOTCH_NGROUPS    = "notchNGroups";
const char* LaminateShellMeshModule::NOTCH_EGROUPS    = "notchEGroups";
const char* LaminateShellMeshModule::ELAS_EGROUPS    = "elasEGroups";
const char* LaminateShellMeshModule::NONE            = "none";
const char* LaminateShellMeshModule::NORMAL          = "normal";
const char* LaminateShellMeshModule::PERIODIC        = "periodic";
const char* LaminateShellMeshModule::TIE_ALL         = "tieAll";
const char* LaminateShellMeshModule::IELEMS_N        = "var.lamMsh.ielemsN";
const char* LaminateShellMeshModule::IELEMS_E        = "var.lamMsh.ielemsE";
const char* LaminateShellMeshModule::TRANSITION      = "transition";

//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------

LaminateShellMeshModule::LaminateShellMeshModule

  ( const String&  name ) :

      Super   ( name   )

{
  active_      = false;
  tieAll_      = false;
  interfaces_  = NORMAL;
  nElGroups_   = 0;
  rank_        = 0;
  numNodes_    = 0;
  numElems_    = 0;
  dz_          = 0.;
  doTieGroups_ = false;
  groupInput_  = newInstance<GroupInputModule> ( name );
}

LaminateShellMeshModule::~LaminateShellMeshModule ()
{}


//-----------------------------------------------------------------------
//   init
//-----------------------------------------------------------------------


Module::Status LaminateShellMeshModule::init

  ( const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  configure_ ( conf, props, globdat );

  if ( active_ )
  {
    System::out() << "Mesh extrusion with LaminateShellMeshModule \n";

    prepareMesh_ ( globdat );

    createShellMesh_ ( globdat );

    if ( ! tieAll_ ) 
    {
      updateConstraints_ ( globdat );

      updateNodeGroups_ ( globdat );
    }
    groupInput_->init ( conf, props, globdat );

    doTieGroups_ &= tieGroupsVert_ ( globdat );
  }

  // return OK when still need to tie nodegroups

  return doTieGroups_ ? OK : DONE;
}

//-----------------------------------------------------------------------
//   run
//-----------------------------------------------------------------------

Module::Status LaminateShellMeshModule::run ( const Properties& globdat )

{
  // Keep trying to tie the specified nodegroups until a DofSpace is found
  // Generally, that will be in the first 'run' 

  doTieGroups_ &= tieGroupsVert_ ( globdat );

  return doTieGroups_ ? OK : DONE;
}

//-----------------------------------------------------------------------
//   shutdown
//-----------------------------------------------------------------------


void LaminateShellMeshModule::shutdown ( const Properties& globdat )
{
}


//-----------------------------------------------------------------------
//   configure_
//-----------------------------------------------------------------------


void LaminateShellMeshModule::configure_

  ( const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  Properties  myConf  = conf .makeProps ( myName_ );
  Properties  myProps = props.findProps ( myName_ );

  // number of layers (if not specified: module inactive)

  if ( ! myProps.find( nLayers_, N_LAYERS ) )
  {
    active_ = false;
    myConf.set( "active", active_ );
    return;
  }
  myConf.set( N_LAYERS, nLayers_ );

  active_ = true;

  // shell elements in 3D space

  nodeRank_ = 3;
  rank_ = nodeRank_ - 1;

  if ( myProps.find( rank_, DIM ) )
  {
    JEM_ASSERT( rank_ == 2 || rank_ == 3 );
  }

  myConf.set( DIM, rank_ );

  // for Shells, rankFac_==1, because the Shell mesh of the plie has the same number of nodes as the 2D initial mesh, therefore often multiplications with
  // factor 1; use private member rankFac_ for that.

  rankFac_ = 1;

  // get ply thickness 

  myProps.get ( dz_, THICKNESS );
  myConf.set  ( THICKNESS, dz_ );

  // option: generate interface elements between layers

  myProps.find( interfaces_, INTERFACES );

  myConf.set  ( INTERFACES, interfaces_ );

  if      ( interfaces_ == NORMAL )
  {
    nInterf_ = nLayers_ - 1;
  }
  else if ( interfaces_ == PERIODIC )
  {
    nInterf_ = nLayers_;
  }
  else if ( interfaces_ == NONE )
  {
    nInterf_ = 0;
  }
  else if ( interfaces_ == TIE_ALL )
  {
    nInterf_ = 0;
    tieAll_  = true;
  }
  else
  {
    throw IllegalInputException (
      JEM_FUNC,
      String (
        "Invalid interfaces type: `" + interfaces_ + "',\nshould be `" + 
        NORMAL  + "', `" + PERIODIC + "' or `" + NONE + "')\n"
      )
    );
  }

  nElGroups_ = nLayers_ + nInterf_;

  // names that will be given to element groups

  layerNames_    .resize( nLayers_ );
  interfaceNames_.resize( nInterf_ );

  if ( ! myProps.find( layerNames_, LAYER_NAMES ) )
  {
    for ( idx_t i = 0; i < nLayers_; ++i )
    {
      layerNames_[i] = String::format( "layer%i", i );
    }
  }

  if ( ! myProps.find( interfaceNames_, INTERFACE_NAMES ) )
  {
    for ( idx_t i = 0; i < nInterf_; ++i )
    {
      interfaceNames_[i] = String::format( "interface%i", i );
    }
  }

  myConf.set(     LAYER_NAMES,     layerNames_ );
  myConf.set( INTERFACE_NAMES, interfaceNames_ );

  // get element group "all"

  if ( ! myProps.contains ( "elements" ) )
  {
    myProps.set ( "elements", "all" );
  }

  const String context = getContext();

  egroup_ = ElemGroup::get ( myConf, myProps, globdat, context );

  // read names of NodeGroups that are to be created

  myProps.find ( groupNames_, NODE_GROUPS );
  myConf.set   ( NODE_GROUPS, groupNames_ );

  nng_ = groupNames_.size();

  // read names of NodeGroups that are to be tied vertically
  // (for uniform load boundary)
  
  if ( myProps.find ( tieGroups_, TIE_GROUPS ) )
  {
    doTieGroups_ = ( tieGroups_.size() > 0 );
    System::out() << "doTieGroups " << doTieGroups_ << endl;
    myConf.set   ( TIE_GROUPS, tieGroups_ );
  }

  // read names of ElemGroups which should not be duplicated

  myProps.find ( skipEGroups_, SKIP_EGROUPS );
  myConf.set   ( SKIP_EGROUPS, skipEGroups_ );

  // read names of ElemGroups which are going to be duplicated to become notches

  myProps.find ( notchEGroups_, NOTCH_EGROUPS );
  myConf.set   ( NOTCH_EGROUPS, notchEGroups_ );

  // number of the notched interface layer

  myProps.find( notchedInterf_, NOTCHED_INTERF );
  myConf.set( NOTCHED_INTERF, notchedInterf_ );

  // read names of ElemGroups which should be treated elastically by models
  
  myProps.find ( elasticEGroups_, ELAS_EGROUPS );
  myConf.set   ( ELAS_EGROUPS, elasticEGroups_ );

  if ( skipEGroups_.size() > 0 && elasticEGroups_.size() > 0 )
  {
    throw IllegalInputException ( JEM_FUNC,
        "combination of skipEGroups and elasticEGroups not allowed!\n" );
  }

  if ( notchEGroups_.size() > 0 && elasticEGroups_.size() > 0 )
  {
    throw IllegalInputException ( JEM_FUNC,
        "combination of notchEGroups and elasticEGroups not allowed!\n" );
  }

  if ( skipEGroups_.size() > 0 && notchEGroups_.size() > 0 )
  {
    throw IllegalInputException ( JEM_FUNC,
        "combination of skipEGroups and notchEGroups not allowed!\n" );
  }

  // read names of NodeGroups which should not be updated
  // (for partial boundary condition)

  myProps.find ( skipNGroups_, SKIP_NGROUPS );
  myConf.set   ( SKIP_NGROUPS, skipNGroups_ );

  // read names of NodeGroups which should ...
  // ()

  myProps.find ( notchNGroups_, NOTCH_NGROUPS );
  myConf.set   ( NOTCH_NGROUPS, notchNGroups_ );
  
  // add nodeGroup "all" to groups that shouldn't be updated

  idx_t nInput = skipNGroups_.size();
  skipNGroups_.reshape ( nInput + 1 );
  skipNGroups_[nInput] = "all";
}

//-----------------------------------------------------------------------
//   prepareMesh_
//-----------------------------------------------------------------------

void LaminateShellMeshModule::prepareMesh_

  ( const Properties&   globdat )

{
  // get mesh data (as generated by InputModule)

  elems_ = toXElementSet ( egroup_.getElements() );
  nodes_ = toXNodeSet    ( elems_.getNodes()     );

  JEM_ASSERT( nodes_.rank() == nodeRank_ );

  numElems_ = elems_.size();
  numNodes_ = nodes_.size();

  doElem_.resize ( numElems_ );
  doNode_.resize ( numNodes_ );

  // prepare for skipping or notching groups

  if ( skipEGroups_.size() > 0 )
  {
    prepareSkip_ ( globdat );
  }
  else if ( notchEGroups_.size() > 0 )
  {
    prepareNotch_ ( globdat );
  }
  else
  {
    // default: do all elements

    doElem_ = true;
    doNode_ = true;

    ielemsY_.resize ( numElems_ );
    ielemsY_ = ( iarray ( numElems_ ) );

    inodesY_.resize ( numNodes_ );
    inodesY_ = ( iarray ( numNodes_ ) );

    if ( elasticEGroups_.size() > 0 ) 
    {
      // store element numbers for which elastic behavior is desired

      prepareElastic_ ( globdat );
    }
  }

  // tie all option overrules skipEGroups or notchEGroups

  if ( tieAll_ ) doNode_ = false;

  numElemY_ = ielemsY_.size();
  numNodeY_ = inodesY_.size();

  numElemN_ = ielemsN_.size();
  numNodeN_ = inodesN_.size();

  eMatrix_.resize (            nLayers_, numElems_ );
  nMatrix_.resize ( rankFac_ * nLayers_, numNodes_ );

  eMatrix_(0,ALL) = iarray ( numElems_ );
  nMatrix_(0,ALL) = iarray ( numNodes_ );

  nodes_.reserve ( rankFac_ * nLayers_ * numNodes_ );
  elems_.reserve ( nElGroups_ * numElems_ );
}

//-----------------------------------------------------------------------
//   prepareSkip_
//-----------------------------------------------------------------------

void LaminateShellMeshModule::prepareSkip_

  ( const Properties&   globdat )

{
  ArrayBuffer<idx_t>  ieN;
  ArrayBuffer<idx_t>  ieY;

  IdxVector         ieThis;
  IdxVector         inodesT;

  Ref<Dict>         groups;
  Ref<DictEnum>     e;
  Ref<ItemGroup>    group;
  String            name;

  doElem_.resize ( numElems_ );
  doNode_.resize ( numNodes_ );

  groups = ItemGroup::getFor( elems_.getData(), globdat );

  // collect data in loop over ElemGroups

  for ( e = groups->enumerate(); ! e->atEnd(); e->toNext() )
  {
    group = dynamicCast<ItemGroup> ( e->getValue() );
    name  = e->getKey();

    ieThis.ref ( group->getIndices() );

    if ( name == "all" || name == "none" || name == "gmshElems" )
    {
      // do nothing
    }
    else if ( testany ( skipEGroups_ == name ) )
    {
      System::out() << "     " << name << " is skipped\n";

      ieN.pushBack ( ieThis.begin(), ieThis.end() );

      select ( doElem_, ieThis ) = false;
    }
    else
    {
      System::out() << "     " << name << " is extruded\n";

      ieY.pushBack ( ieThis.begin(), ieThis.end() );

      select ( doElem_, ieThis ) = true;
    }
  }

  ielemsY_.ref ( ieY.toArray() );
  inodesY_.ref ( elems_.getUniqueNodesOf ( ielemsY_ ) );

  ielemsN_.ref ( ieN.toArray() );
  inodesN_.ref ( elems_.getUniqueNodesOf ( ielemsN_ ) );

  doNode_ = true;
  select ( doNode_, inodesN_ ) = false;

  // store in globdat (used in XFEMModel)

  globdat.set ( IELEMS_N, ielemsN_ );

  // check if nodes are specified
  // (gmshInput.doElemGroups must be set to true!)

  System::out() << "checking " << ielemsN_.size() << " + " << ielemsY_.size() << " =? "
    << elems_.size() << endl;

  if ( ielemsN_.size() + ielemsY_.size() != elems_.size() )
  {
    throw IllegalInputException ( JEM_FUNC, 
        String( "Not all elements are in an ElementGroup, ") +
        String("set gmshInput.doElemGroups to true") );
  }
}

//-----------------------------------------------------------------------
//   prepareNotch_
//-----------------------------------------------------------------------

void LaminateShellMeshModule::prepareNotch_

  ( const Properties&   globdat )

{
  ArrayBuffer<idx_t>  ieN;
  ArrayBuffer<idx_t>  ieY;

  IdxVector         ieThis;
  IdxVector         inodesT;

  Ref<Dict>         groups;
  Ref<DictEnum>     e;
  Ref<ItemGroup>    group;
  String            name;

  doElem_.resize ( numElems_ );
  doNode_.resize ( numNodes_ );

  groups = ItemGroup::getFor( elems_.getData(), globdat );

  // collect data in loop over ElemGroups

  for ( e = groups->enumerate(); ! e->atEnd(); e->toNext() )
  {
    group = dynamicCast<ItemGroup> ( e->getValue() );
    name  = e->getKey();

    ieThis.ref ( group->getIndices() );

    if ( name == "all" || name == "none" || name == "gmshElems" )
    {
      // do nothing
    }
    else if ( testany ( notchEGroups_ == name ) )
    {
      System::out() << "     " << name << " is notched\n";

      ieN.pushBack ( ieThis.begin(), ieThis.end() );

      select ( doElem_, ieThis ) = true;
    }
    else
    {
      System::out() << "     " << name << " is extruded\n";

      ieY.pushBack ( ieThis.begin(), ieThis.end() );

      select ( doElem_, ieThis ) = true;
    }
  }

  ielemsN_.ref ( ieN.toArray() );
  inodesN_.ref ( elems_.getUniqueNodesOf ( ielemsN_ ) );

  ielemsY_.ref ( ieY.toArray() );
  inodesY_.ref ( elems_.getUniqueNodesOf ( ielemsY_ ) );

  doNode_ = true;

  // store in globdat (used in XFEMModel)

  globdat.set ( IELEMS_N, ielemsN_ );

  // check if nodes are specified
  // (gmshInput.doElemGroups must be set to true!)

  System::out() << "checking " << ielemsN_.size() << " + " << ielemsY_.size() << " =? "
    << elems_.size() << endl;

  if ( ielemsN_.size() + ielemsY_.size() != elems_.size() )
  {
    throw IllegalInputException ( JEM_FUNC, 
        String( "Not all elements are in an ElementGroup, ") +
        String("set gmshInput.doElemGroups to true") );
  }
}
  
//-----------------------------------------------------------------------
//   prepareElastic_
//-----------------------------------------------------------------------

void LaminateShellMeshModule::prepareElastic_

  ( const Properties&   globdat )

{
  ArrayBuffer<idx_t>  ieN;
  IdxVector         ieThis;

  Ref<Dict>         groups;
  Ref<DictEnum>     e;
  Ref<ItemGroup>    group;
  String            name;

  groups = ItemGroup::getFor( elems_.getData(), globdat );

  // collect data in loop over ElemGroups

  for ( e = groups->enumerate(); ! e->atEnd(); e->toNext() )
  {
    group = dynamicCast<ItemGroup> ( e->getValue() );
    name  = e->getKey();

    ieThis.ref ( group->getIndices() );

    if ( testany ( elasticEGroups_ == name ) )
    {
      System::out() << "     " << name << 
        " will have elastic interfaces\n";

      ieN.pushBack ( ieThis.begin(), ieThis.end() );
    }
  }

  ielemsN_.ref ( ieN.toArray() );

  // store in globdat (used in XFEMModel and NCInterfaceModel)

  globdat.set ( IELEMS_N, ielemsN_ );

  globdat.set ( IELEMS_E, ielemsN_ );
}

//-----------------------------------------------------------------------
//   createShellMesh_
//-----------------------------------------------------------------------

void LaminateShellMeshModule::createShellMesh_

  ( const Properties&   globdat )

{
  JEM_ASSERT ( nodeRank_ == 3 );
  JEM_ASSERT ( rank_ == 2 );

  Assignable<ElemGroup> newGroup;

  idx_t nodeCount = elems_.maxElemNodeCount();

  IdxVector     iinodes0;
  IdxVector     jjnodes0;
  IdxVector     jjnodes1;

  Matrix        coords ( nodeRank_ , numNodes_ );
  IdxMatrix     ielems ( nElGroups_, numElems_ );
  IdxMatrix     jelemsY( nInterf_  , numElemY_ );
  IdxMatrix     jelems ( nInterf_  , numElems_ );
  IdxVector     inodes ( nodeCount );
  IdxVector     jnodes ( nodeCount );
  IdxVector     iinodes( nodeCount * 2 );
  IdxVector     jjnodes( nodeCount * 2 );

  iinodes0.ref  ( iinodes [ slice(  0,nodeCount) ] );
  jjnodes0.ref  ( jjnodes [ slice(  0,nodeCount) ] );
  jjnodes1.ref  ( jjnodes [ slice(nodeCount,END) ] );

  // get mesh data (as generated by InputModule)

  nodes_.getCoords ( coords );

  System::out() << " ...Generating " << nLayers_ << " layers of "
    << "shell elements, with element group names \n    "
    << layerNames_ << "\n";

  // generate nodes
  // for doNode==true: one node per layer
  // for doNode==false: skip duplicate nodes (related to interface)

  for ( idx_t iLayer = 1; iLayer < nLayers_; ++iLayer )
  {
    // set z-coords

    idx_t iz = iLayer ;

    coords( 2, ALL ) = iz * dz_;

    // add node to set

    for ( idx_t in = 0; in < numNodes_; ++in )
    {
      if ( doNode_[in] )
      {
        nMatrix_(iLayer,in) = nodes_.addNode ( coords(ALL,in) );
      }
      else
      {
        nMatrix_(iLayer,in) = nMatrix_(iLayer-1,in);
      }
    }
  }

  // generate the first layer of shell elements  
  // (irrespective of doElem)

  for ( idx_t ie = 0; ie < numElems_; ++ie )
  {
    elems_.getElemNodes ( inodes, ie );

    elems_.setElemNodes ( ie, inodes );
  }

  // generate additional layers of shell elements
  // (irrespective of doElem)

  for ( idx_t iLayer = 1; iLayer < nLayers_; ++iLayer )
  {
    idx_t izNode0 = iLayer;

    for ( idx_t ie = 0; ie < numElems_; ++ie )
    {
      elems_.getElemNodes ( inodes, ie );
      
      jnodes = select ( nMatrix_(izNode0,ALL), inodes );

      elems_.addElement ( jnodes );
    }
  }

  // store shell element groups

  for ( idx_t iLayer = 0; iLayer < nLayers_; ++iLayer )
  {
    ielems( iLayer, ALL ) = iarray( numElems_ ) + iLayer * numElems_;

    newGroup = newElementGroup( ielems(iLayer,ALL), elems_ );

    newGroup.store( layerNames_[iLayer], globdat );
  }

  // generate interface elements (if doElem==true)

  if ( nInterf_ > 0 )
  {
    System::out() << " ...Generating " << nInterf_ << " layers of "
      << "interface elements (" << interfaces_ << "), with element group "
      << "names \n    " << interfaceNames_ << "\n";

    for ( idx_t iLayer = 0; iLayer < nInterf_; ++iLayer )
    {
      idx_t izNode0 = iLayer  * 1;
      idx_t izNode1 = izNode0 + 1;

      if ( iLayer == notchedInterf_ - 1 )
      {
        for ( idx_t ieY = 0; ieY < numElemY_; ++ieY )
        {
          idx_t ie = ielemsY_[ieY];
  
          JEM_ASSERT ( doElem_[ie] );
  
          elems_.getElemNodes ( iinodes, ie );
  
          jjnodes0 = select ( nMatrix_(izNode0,ALL), iinodes0 );
  
          jjnodes1 = select ( nMatrix_(izNode1,ALL), iinodes0 );
  
          jelemsY(iLayer,ieY) = elems_.addElement ( jjnodes );
        }

        // store group

        newGroup = newElementGroup( jelemsY(iLayer,ALL), elems_ );

        newGroup.store( interfaceNames_[iLayer], globdat );
      }
      else 
      {
        for ( idx_t ie = 0; ie < numElems_; ++ie )
        {
          // idx_t ie = ielemsY_[ieY];
  
          // JEM_ASSERT ( doElem_[ie] );
  
          elems_.getElemNodes ( iinodes, ie );
  
          jjnodes0 = select ( nMatrix_(izNode0,ALL), iinodes0 );
  
          jjnodes1 = select ( nMatrix_(izNode1,ALL), iinodes0 );
  
          jelems(iLayer,ie) = elems_.addElement ( jjnodes );
        }

        // store group
      
        newGroup = newElementGroup( jelems(iLayer,ALL), elems_ );

        newGroup.store( interfaceNames_[iLayer], globdat );
      }
    }
  }
  else
  {
    System::out() << " ...Generating no interface elements.\n";
  }
}


//-----------------------------------------------------------------------
//   updateConstraints_
//-----------------------------------------------------------------------


void LaminateShellMeshModule::updateConstraints_

  ( const Properties&   globdat )

{
  // This function doesn't do anything, except that it gives a warning
  // when constraints have been applied before the mesh is extracted.
  
  using jive::util::ConstraintsParser;

  StringVector            list = ItemSet::listAll ( globdat );

  Ref<ConstraintsParser>  conParser;
  Ref<Constraints>        cons;
  Ref<ItemSet>            items;

  idx_t                   i, n;

  for ( i = 0, n = list.size(); i < n; i++ )
  {
    items = ItemSet::find ( list[i], globdat );

    if ( items == NIL )
    {
      continue;
    }

    conParser = ConstraintsParser::extract ( items, globdat );

    if ( conParser != NIL )
    {
      if ( conParser->slaveDofCount() > 0 )
      {
        System::warn() << "Constraints are not copied to other layers\n";
      }
    }
  }
}

//-----------------------------------------------------------------------
//   updateNodeGroups_
//-----------------------------------------------------------------------


void LaminateShellMeshModule::updateNodeGroups_

  ( const Properties&   globdat )

{
  Ref<Dict>         groups;
  Ref<DictEnum>     e;
  Ref<XItemGroup>   group;
  String            name;

  groups = ItemGroup::getFor( nodes_.getData(), globdat );

  for ( e = groups->enumerate(); ! e->atEnd(); e->toNext() )
  {
    group = dynamicCast<XItemGroup> ( e->getValue() );
    name  = e->getKey();

    if ( ! ( group == NIL || testany ( skipNGroups_ == name )  || testany ( notchNGroups_ == name ) ) )
    {
      IdxVector inGr ( group->getIndices() );

      ArrayBuffer<idx_t> jnodeBuf;

      for ( idx_t i = 0; i < inGr.size(); ++i )
      {
        idx_t in = inGr[i];

        for ( idx_t iLayer = 1; iLayer < nMatrix_.size(0); ++iLayer )
        {
          idx_t candidate = nMatrix_(iLayer,in);

          if ( candidate != nMatrix_(iLayer-1,in) )
          {
            jnodeBuf.pushBack ( candidate );
          }
        }
      }
      group->append ( jnodeBuf.toArray() );

      System::out() << " ...Nodegroup `" << name <<
        "' now contains " << group->size() << " nodes.\n";
    }
  }
}

//-----------------------------------------------------------------------
//   tieGroupsVert_
//-----------------------------------------------------------------------

bool LaminateShellMeshModule::tieGroupsVert_

  ( const Properties&   globdat )

{
  Ref<DofSpace>     dofs = DofSpace::find ( nodes_.getData(), globdat );

  if ( dofs == NIL ) return true;

  Ref<Constraints>  cons = Constraints::get ( dofs, globdat );

  idx_t             ntype  ( dofs->typeCount() );

  IdxVector         types  ( iarray ( ntype )  );
  IdxVector         idofs  ( ntype             );
  IdxVector         jdofs  ( ntype             );

  Ref<Dict>         groups;
  Ref<DictEnum>     e;
  Ref<ItemGroup>    group;

  groups = ItemGroup::getFor( nodes_.getData(), globdat );

  for ( e = groups->enumerate(); ! e->atEnd(); e->toNext() )
  {
    if ( testany ( e->getKey() == tieGroups_ ) )
    {
      System::out() << " ...Vertically tying nodes from nodeGroup `" 
        << e->getKey() << "'.\n";

      group = dynamicCast<ItemGroup> ( e->getValue() );

      JEM_PRECHECK ( group->size() > 0 );

      IdxVector inodes ( group->getIndices() );

      for ( idx_t in = 0; inodes[in] < numNodes_; ++in )
      {
        idx_t inode = inodes[in];

        dofs->getDofsForItem ( idofs, types, inode );

        for ( idx_t iLayer = 1; iLayer < nLayers_*rankFac_; ++iLayer )
        {
          idx_t jnode = nMatrix_ ( iLayer, inode );

          dofs->getDofsForItem ( jdofs, types, jnode );

          for ( idx_t id = 0; id < types.size(); ++ id )
          {
            try 
            {
              cons->addConstraint ( jdofs[id], idofs[id], 1. );
            }
            catch ( const jem::Exception& ex )
            {
              System::out() << "no constraint added for nodes " <<
                inode << " and " << jnode << " : " << id << endl;
            }
          }
        }
      }
    }
  }
  System::out() << endl;

  return false;
}

//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------

Ref<Module>  LaminateShellMeshModule::makeNew

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
//   declareLaminateShellMeshModule
//-----------------------------------------------------------------------

void declareLaminateShellMeshModule ()
{
  using jive::app::ModuleFactory;

  ModuleFactory::declare ( LaminateShellMeshModule::TYPE_NAME,
                         & LaminateShellMeshModule::makeNew );
}
