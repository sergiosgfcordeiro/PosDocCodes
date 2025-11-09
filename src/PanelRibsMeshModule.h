/*
 *  Copyright (C) 2025 TU Delft. All rights reserved.
 *  
 *  Sérgio Gustavo Ferreira Cordeiro, July 2025. 
 *  
 *  This module inserts BLine2 PanelRibs structural cohesive elements for debonding analysis
 *  in a preexisting PanelRibs mesh.
 *
 *  - nodes on the panel-ribs interface are multiplied (new nodes goes for the ribs and old nodes for the panels)
 *  - An element group with BLine2 interface elements is defined between the panel and the ribs
 *  - NB: constraints are not copied to new plies
 *
 */

#ifndef PANEL_RIBS_MESH_MODULE_H
#define PANEL_RIBS_MESH_MODULE_H

#include <jive/app/Module.h>
#include <jive/fem/ElementGroup.h>
#include <jive/fem/XElementSet.h>
#include <jive/fem/XNodeSet.h>
#include <jive/util/Constraints.h>
#include <jive/util/XDofSpace.h>
#include <jive/util/Assignable.h>

#include "GroupInputModule.h"

using namespace jem;

using jem::util::Properties;
using jive::BoolVector;
using jive::IdxVector;
using jive::IdxMatrix;
using jive::StringVector;
using jive::app::Module;
using jive::fem::ElementGroup;
using jive::fem::XElementSet;
using jive::fem::XNodeSet;
using jive::util::Constraints;
using jive::util::XDofSpace;
using jive::util::Assignable;

typedef ElementGroup        ElemGroup;
typedef XElementSet         XElemSet;

//-----------------------------------------------------------------------
//   class PanelRibsMeshModule
//-----------------------------------------------------------------------


class PanelRibsMeshModule : public Module
{
 public:

  typedef PanelRibsMeshModule Self;
  typedef Module             Super;

  static const char*        DOF_NAMES[5];
  static const char*        TYPE_NAME;
  static const char*        DIM;
  static const char*        INTERFACES;
  static const char*        PANEL_NAMES;
  static const char*        RIBS_NAMES;
  static const char*        INTERFACE_NAMES;
  static const char*        PANEL_EGROUPS;
  static const char*        RIBS_EGROUPS;
  static const char*        NORMAL;
  static const char*        NONE;
  static const char*        IELEMS_PANEL;

  explicit                  PanelRibsMeshModule

    ( const String&           name   = "PanelRibsMesh" );

  virtual Status            init

    ( const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  virtual Status            run

    ( const Properties&       globdat );

  virtual void              shutdown

    ( const Properties&       globdat );

  static Ref<Module>        makeNew

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

 protected:

  virtual                  ~PanelRibsMeshModule  ();

 private:

  idx_t                    nPanel_;
  idx_t                    nRibs_;
  idx_t                    nInterf_;
  idx_t                    notchedInterf_;
  idx_t                    nElGroups_;
  idx_t                    nodeRank_;
  idx_t                    rank_;

  idx_t                    numNodes_;
  idx_t                    numElems_;

  idx_t                    numInterfNodes_;
  idx_t                    numInterfElems_;

  String                   interfaces_;
  bool                     active_;

  StringVector             panelNames_;
  StringVector             ribsNames_;
  StringVector             interfaceNames_;
  // StringVector             skipNGroups_;
  // StringVector             notchNGroups_;

  Assignable<ElemGroup>    egroup_;
  Assignable<XElemSet>     elems_;
  Assignable<XNodeSet>     nodes_;

  Ref<Constraints>         cons_;
  Ref<GroupInputModule>    groupInput_;

  // some quantities for panel and ribs element groups

  StringVector             panelEGroups_;
  StringVector             ribsEGroups_;

  IdxVector                ipanelelems_; // 
  IdxVector                iribselems_;  // 
  IdxVector                iribselemsY_; // Yes: will be duplicated
  IdxVector                iribselemsN_; // No: will not be duplicated
  BoolVector               doElem_;   // same information, stored differently
  BoolVector               doNode_;  

  IdxVector                inodesY_;  // nodes from ielemsY_
  IdxVector                inodesN_;  // nodes from ielemsN_

  IdxVector                nVector_;
  IdxVector                eVector_;

  idx_t                    numElemY_;
  idx_t                    numNodeY_;
  idx_t                    numElemN_;
  idx_t                    numNodeN_;

  // functions

  void                     configure_ 

    ( const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  void                     prepareMesh_

    ( const Properties&   globdat );

  void                     prepareSkipPanel_

  ( const Properties&   globdat );

  void                     createPanelRibsMesh_ 

    ( const Properties&   globdat );

  void                     getInterfaceElems_
  
  ( idx_t&             numInterfElems_,
    const Matrix&      coords,
    const Vector&       v )   const;
};


#endif
