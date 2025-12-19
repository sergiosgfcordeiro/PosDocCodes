/*
 *
 * Model for Panel-Ribs interface elements with 6 DOFs per node: 
 *   - assembly of stiffness and internal force vector
 *   - output
 * 
 * Sergio Gustavo Ferreira Cordeiro, July 2025
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
#include "Shell6DofsAllmanModel.h"
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

class PanelRibsInterface6DofsAllmanModel : public Model
{
 public:

  typedef PanelRibsInterface6DofsAllmanModel     Self;
  typedef Model              Super;

  static const char*         DOF_NAMES[6];
  static const char*         SHAPE_PROP;
  static const char*         COHEMAT_PROP;
  static const char*         PANRIB_SHAPE_PROP;
  static const char*         MATERIAL_PROP[2];
  static const char*         THICK_PROP[2];
  static const char*         THICK_SHAPE_PROP;
  static const char*         PAN_ORIENTATION_PROP[3];
  static const char*         RIB_ORIENTATION_PROP[3];
  static const char*         PAN_NAME_PROP;
  static const char*         RIB_NAME_PROP;

                       PanelRibsInterface6DofsAllmanModel
			 
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

  virtual              ~PanelRibsInterface6DofsAllmanModel ();

 protected:

  virtual void         getMatrix_

  ( Ref<MatrixBuilder>  mbuilder,
    const Vector&       force,
    const Vector&       pan_disp,
    const Vector&       rib_disp )       const;

  virtual void         getFrictionForce_

  ( const Vector&         fsh0,
  const Vector&         pan_disp,
  const Vector&         rib_disp )       const;

  virtual void         writeOutput_

  ( const Properties&   globdat )    const;

  virtual void         getDissipation_ 

  ( const Properties&   params ) const;

  void                 getThick_xcoords_

  ( Matrix&           thick_xcoords,
    const double&     h_rib )    const;

  void                 getConnectivity_

  ( IdxVector&           Connect,
    IdxVector&           ConnectInv )    const;

  void                 getPanRibElements_

  ( idx_t&               ielem_panSideA,
    idx_t&               ielem_panSideB,
    idx_t&               ielem_rib,
    const ElemSet&       pan_elems_,
    const ElemSet&       rib_elems_ ,
    const Matrix&        coords )     const;

  void                 getTransMatrix_
  
  ( Matrix&             transMat,
    const Matrix&       coords,
    const Vector&       v )   const;
  
  void                 getTransMatrixDof_
    
    ( Matrix&           transMatDof,
      const Matrix&     transMat_pan,
      const Matrix&     transMat_rib )   const;
  
  void                 get2DLocalcoordinates_

  ( Matrix&             xcoords,
    const Matrix&       transMat,
    const Matrix&       coords )   const;

  void                 get2DGlobalcoordinates_

  ( Matrix&             coords2D,
    const Matrix&       transMat,
    const Matrix&       coords )   const;

  void                 getArea_

    ( double&             Area,
      const Matrix&       xcoords)    const;
  
  void                 get3DLocalcoordinates_

  ( Matrix&             ipxcoords_pr,
    const Matrix&       transMat_pr,
    const Matrix&       ipxcoords,
    const Matrix&       transMat_int,
    const Matrix&       coords_pr,
    const Matrix&       coords )   const;

  void                 getBABxBymats_

  ( Matrix&             BAmat,
    Matrix&             Bxmat,
    Matrix&             Bymat )   const;

  void                 getMAMamats_

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

  ( Vector&             Disp_mem_panA,
    Vector&             Disp_mem_panB,
    Vector&             Disp_mem_rib,
    const Vector&       elemxDisp0 )   const;

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
    const Matrix&     xcoords,
    const double&     Area,
    const double&     x,
    const double&     y )    const;

  void                 getNuyNvxmats_

  ( Matrix&           Nuy,
    Matrix&           Nvx,
    const Matrix&     xcoords,
    const double&     Area,
    const double&     x,
    const double&     y )    const;

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
    
protected: 

  String                  myTag_;
  String                  myBase_;

  // model rank
  
  idx_t                   nodeRank_;
  idx_t                   rank_;
  idx_t                   thick_rank_;

  // interface shape rank

  idx_t                   qRank_;

  idx_t                   nodeCount_;
  idx_t                   nodeCount_BLine2_;
  idx_t                   nodeCount_Line2_;
  idx_t                   nodeCount_panrib_;
  idx_t                   elemCount_;
  idx_t                   pan_elemCount_;
  idx_t                   rib_elemCount_;
  idx_t                   ipCount_;
  idx_t                   jpCount_;
  idx_t                   ipCount_panrib_;
  double                  thickness_[2];
  double                  pan_orientVec_[3];
  double                  rib_orientVec_[3];

  IdxVector               ielems_;
  IdxVector               pan_ielems_;
  IdxVector               rib_ielems_;
  IdxVector               panrib_ielems_;

  Assignable<ElemGroup>   egroup_;
  Assignable<ElemSet>     elems_;
  Assignable<NodeSet>     nodes_;

  Assignable<ElemGroup>   pan_egroup_;
  Assignable<ElemSet>     pan_elems_;
  Assignable<NodeSet>     pan_nodes_;

  Assignable<ElemGroup>   rib_egroup_;
  Assignable<ElemSet>     rib_elems_;
  Assignable<NodeSet>     rib_nodes_;

  Assignable<ElemSet>     panrib_elems_;
  Assignable<NodeSet>     panrib_nodes_;

  Ref<InterfaceShape>     shape_;
  Ref<IShape>             thick_shape_;

  String                  panel_;
  String                  ribs_;

  Ref<XDofSpace>          pan_dofs_;
  Ref<XDofSpace>          rib_dofs_;
  IdxVector               pan_dofTypes_;
  IdxVector               rib_dofTypes_;

  Ref<Shell6DofsAllmanModel>  shell_[2];

  Ref<Material>           material_[2];
  Ref<Softening>          softening_[2];
  bool                    crackBandMethod_[2];
  Vector                  charLength_[2];

  Ref<CohesiveMat>        coheMat_;

  Ref<AlfanoTuronCoheMat> frictionMat_;

  Ref<IShape>             panribshape_;

  ShapeGradsFunc          getShapeGrads_;
  ShapeFuncsFunc          getShapeFuncs_;
  
  Ref<PrintWriter>        xOut_;
};

