/*
 *
 *  This module reorders quadrilateral elements for usage 
 *  as interface elements.
 *
 *  Frans van der Meer, January 2015
 *  
 */

#ifndef QUAD4INTERF_MODULE_H
#define QUAD4INTERF_MODULE_H

#include <jem/util/Properties.h>
#include <jive/app/Module.h>
#include <jive/fem/NodeSet.h>
#include <jive/fem/XElementSet.h>
#include <jive/fem/ElementGroup.h>
#include <jive/util/Assignable.h>

using jem::String;
using jem::Ref;
using jem::idx_t;
using jem::newInstance;
using jem::util::Properties;
using jive::StringVector;
using jive::app::Module;
using jive::fem::NodeSet;
using jive::fem::XElementSet;
using jive::fem::ElementGroup;
using jive::util::Assignable;
using jive::IdxVector;
using jive::Matrix;

typedef XElementSet         XElemSet;
typedef ElementGroup        ElemGroup;

//-----------------------------------------------------------------------
//   class Quad4InterfModule
//-----------------------------------------------------------------------

class Quad4InterfModule : public Module
{

 public:

  typedef Quad4InterfModule  Self;
  typedef Module             Super;

  static const char*        TYPE_NAME;
  static const char*        ELEM_GROUPS;
  static const char*        THICK_PROP;

  explicit                  Quad4InterfModule

    ( const String&           name   = "quad4interf" );

  virtual Status            init

    ( const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  static Ref<Module>        makeNew

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

 protected:

  virtual                  ~Quad4InterfModule  ();

  void                        getThickness_

    ( double                  thickness ) const;

  void                      reorderElements_

    ( const String&           name,
      const Properties&       globdat );

  void                      counterClockwise_

    ( const IdxVector&        iperm,
      const IdxVector&        ielems ) const;

  idx_t                     findPermutation_

    (       IdxVector&        iperm,
      const IdxVector&        ielems ) const;

 private:

  StringVector             elemGroups_;

  Assignable<XElemSet>     elems_;
  Assignable<NodeSet>      nodes_;

  idx_t                    rank_;
  double                  thickness_;
};

#endif
