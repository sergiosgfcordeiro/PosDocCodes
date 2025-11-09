/*
 *
 * Model for structural interface elements: 
 *   - assembly of stiffness and internal force vector
 *   - output
 * 
 * Sergio Gustavo Ferreira Cordeiro, May 2025
 *
 */

#include <jem/io/PrintWriter.h>
#include <jem/io/FileWriter.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/MatmulChain.h>
#include <jem/util/Properties.h>

#include <jive/util/utilities.h>
#include <jive/util/XDofSpace.h>
#include <jive/util/Assignable.h>
#include <jive/algebra/MatrixBuilder.h>
#include <jive/model/Model.h>
#include <jive/geom/IntegrationSchemes.h>
#include <jive/geom/InterfaceShape.h>
#include <jive/geom/IShapeFactory.h>
#include <jive/geom/InternalShape.h>
#include <jive/fem/ElementGroup.h>

#include "utilities.h"
#include "CohesiveMat.h"
#include "AlfanoTuronCoheMat.h"
#include "Material.h"
#include "Softening.h"
#include "Plasticity.h"

using namespace jem;

using jem::io::PrintWriter;
using jem::io::FileWriter;
using jem::numeric::matmul;
using jem::numeric::MatmulChain;
using jem::util::Properties;
using jive::Cubix;
using jive::IdxVector;
using jive::Vector;
using jive::StringVector;
using jive::util::joinNames;
using jive::util::XDofSpace;
using jive::util::Assignable;
using jive::algebra::MatrixBuilder;
using jive::model::Model;
using jive::geom::InterfaceShape;
using jive::geom::BShape;
using jive::geom::IShape;
using jive::fem::NodeSet;
using jive::fem::ElementSet;
using jive::fem::ElementGroup;

// some typedef to avoid typing

typedef ElementSet              ElemSet;
typedef ElementGroup            ElemGroup;
typedef MatmulChain<double,3>   MChain3;
typedef MatmulChain<double,2>   MChain2;
typedef MatmulChain<double,1>   MChain1;

class StructuralInterfaceModel : public Model
{
 public:

  typedef StructuralInterfaceModel     Self;
  typedef Model              Super;

  static const char*         DOF_NAMES[5];
  static const char*         SHAPE_PROP;
  static const char*         COHEMAT_PROP;
  static const char*         BOTTOP_SHAPE_PROP;
  static const char*         MATERIAL_PROP[2];
  static const char*         THICK_PROP[2];
  static const char*         BOT_ORIENTATION_PROP[3];
  static const char*         TOP_ORIENTATION_PROP[3];
  static const char*         BOT_NAME_PROP;
  static const char*         TOP_NAME_PROP;

                       StructuralInterfaceModel
			 
    ( const String&       name,
      const Properties&   conf,
      const Properties&   props,
      const Properties&   globdat );

  virtual void         configure

    ( const Properties&   props,
      const Properties&   globdat );

  virtual void         getConfig

    ( const Properties&   conf,
      const Properties&   globdat )      const;

  virtual bool         takeAction

    ( const String&       action,
      const Properties&   params,
      const Properties&   globdat );

 protected:

  virtual              ~StructuralInterfaceModel ();

 protected:

  virtual void         getMatrix_

  ( Ref<MatrixBuilder>  mbuilder,
    const Vector&       force,
    const Vector&       disp )       const;

  virtual void         getFrictionForce_

  ( const Vector&       fsh0,
  const Vector&       disp )       const;

  virtual void         writeOutput_

  ( const Properties&   globdat )    const;

  virtual void         getDissipation_ 

  ( const Properties&   params ) const;

  void                 getConnectivity_

  ( IdxVector&           Connect,
    IdxVector&           ConnectInv )    const;

  void                 getTransMatrix_
  
  ( Matrix&             transMat,
    const Matrix&       coords_bot,
    const Matrix&       coords_top,
    const Vector&       v_bot,
      const Vector&     v_top )   const;

  void                 getLocalSystem_

  ( Matrix&             ipxcoords,
    Matrix&             xcoords,
    const Matrix&       ipcoords,
    const Matrix&       coords,
    const Vector&       v )   const;

  void                 getBABxBymats_

  ( Matrix&             BAmat,
    Matrix&             Bxmat,
    Matrix&             Bymat )   const;

  void                 getMAMamat_

  ( Matrix&           MAmat_,
    Matrix&           Mamat_,
      const Matrix&   xcoords )    const;

  void                 invertMAmat_

    ( Matrix&           MAinv,
      const Matrix&     MAmat_ )     const;

  void                 getTmat_

  ( Matrix&             Tmat,
    const Matrix&       xcoords )    const;

  void                 getMembraneDisp_

  ( Vector&             Disp_mem_bot,
    Vector&             Disp_mem_top,
    const Vector&       elemDisp0 )   const;

  void                 getHmatAnalytic_

    ( Matrix&           Hmat,
      const Matrix&     Dplate,
      const Matrix&     xcoords )    const;

  void                 getBmat_

    ( Matrix&           Bmat,
      const Matrix&     Dplate,
      const Matrix&     xcoords )    const;

  void                 getNuNvmats_

  ( Matrix&           Nu,
    Matrix&           Nv,
    const Matrix&     sfuncs,
    const idx_t&      ip )    const;

  void                 getSRRxRymats_

  ( Matrix&           Smat,
    Matrix&           Rmat,
    Matrix&           Rxmat,
    Matrix&           Rymat,
    const double&     x,
    const double&     y )    const;

  void                 initWriter_

    ( const Properties&   params );

  void                 writeGeom_()    const;

  void                 checkGeom_()   const;
    
protected: 

  String                  myTag_;
  String                  myBase_;

  // model rank
  
  idx_t                   nodeRank_;
  idx_t                   rank_;

  // interface shape rank

  idx_t                   qRank_;

  idx_t                   nodeCount_;
  idx_t                   nodeCount_bottop_;
  idx_t                   elemCount_;
  idx_t                   ipCount_;
  idx_t                   ipCount_bottop_;
  double                  thickness_[2];
  double                  bot_orientVec_[3];
  double                  top_orientVec_[3];

  IdxVector               ielems_;

  Assignable<ElemGroup>   egroup_;
  Assignable<ElemSet>     elems_;
  Assignable<NodeSet>     nodes_;

  Ref<InterfaceShape>     shape_;

  // String                  ischeme_;
  String                  bot_;
  String                  top_;

  Ref<XDofSpace>          dofs_;
  IdxVector               dofTypes_;

  Ref<Material>           material_[2];
  Ref<Softening>          softening_[2];
  bool                    crackBandMethod_[2];
  Vector                  charLength_[2];

  Ref<CohesiveMat>        coheMat_;

  Ref<AlfanoTuronCoheMat> frictionMat_;


  Ref<IShape>             bottopshape_;

  ShapeGradsFunc          getShapeGrads_;
  ShapeFuncsFunc          getShapeFuncs_;
  
  Ref<PrintWriter>        xOut_;
};

Ref<BShape> makeBTriangleNC  ( const Properties& props );

Ref<BShape> makeBTriangleNC3 ( const String& name = "BTriangleNC3" );

Ref<BShape> makeBTriangleNC6 ( const String& name = "BTriangleNC6" );

Ref<BShape> makeBTriangleNC13 ( const String& name = "BTriangleNC13" );

