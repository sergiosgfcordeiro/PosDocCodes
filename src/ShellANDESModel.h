/*
 *
 *  Copyright (C) 2025 TU Delft. All rights reserved.
 *
 *  
 *  This class implements the ANDES High Performance triangular Thin (Kirchhoff-Love) shell element by Felippa and Militello. It is an assemble of the AQR ANDES 
 *  template for plate bending triangle and the optimal ANDES membrane triangle with drilling degrees of freedom:
 *  
 *  C. Militello and C. A Felippa, The first ANDES elements: 9-dof plate bending triangles, Computer Meth. Applied Mechanics and Engineering., 93, 217–246, 1991.
 *  C. A. Felippa and C. Militello, Membrane triangles with corner drilling freedoms: II. The ANDES element, Finite Elem. Anal. Des., 12, 189–201, 1992.
 *  C. A Felippa, Advanced Finite Element Methods, Lecture Notes, University of Colorado at Boulder, 1994. 
 *  C. A. Felippa and S. Alexander, Membrane triangles with corner drilling freedoms: III. Implementation and performance evaluation, Finite Elem. Anal. Des., 12, 203–239, 1992.
 *  
 *  Author:  S.G.F. Cordeiro, sferreiracorde@tudelft.nl
 *
 *  25 August 2025:
 *    + triangle3 elements
 */

#ifndef SHELLANDES_MODEL_H
#define SHELLANDES_MODEL_H

#include <jem/io/PrintWriter.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/util/Properties.h>
#include <jem/util/SparseArray.h>
#include <jem/util/StringUtils.h>

#include <jive/util/utilities.h>
#include <jive/util/XDofSpace.h>
#include <jive/util/XTable.h>
#include <jive/util/Assignable.h>
#include <jive/algebra/MatrixBuilder.h>
#include <jive/model/Model.h>
#include <jive/model/Actions.h>
#include <jive/model/StateVector.h>
#include <jive/geom/IShapeFactory.h>
#include <jive/geom/InternalShape.h>
#include <jive/geom/Triangle.h>
#include <jive/fem/ElementGroup.h>
#include <jive/fem/Globdat.h>

#include "models.h"
#include "utilities.h"
#include "Material.h"
#include "SolverNames.h"
#include "Softening.h"
#include "Plasticity.h"
#include "DispWriter.h"

using namespace jem::io;

using namespace jem;

using jem::io::PrintWriter;
using jem::util::Properties;
using jem::util::SparseArray;
using jem::util::StringUtils;
using jive::Vector;
using jive::IdxVector;
using jive::StringVector;
using jive::Cubix;
using jive::util::XDofSpace;
using jive::util::XTable;
using jive::util::Assignable;
using jive::algebra::MatrixBuilder;
using jive::model::Model;
using jive::model::StateVector;
using jive::geom::IShape;
using jive::fem::NodeSet;
using jive::fem::ElementSet;
using jive::fem::ElementGroup;
using jive::fem::Globdat;

// some typedef to avoid typing

typedef ElementSet              ElemSet;
typedef ElementGroup            ElemGroup;


//=======================================================================
//   class ShellANDESModel
//=======================================================================

// Forward declaration
class StructuralInterfaceModel;
class PanelRibsInterface6DofsAlmanModel;


class ShellANDESModel : public Model
{
 friend class StructuralInterfaceModel; 
 friend class PanelRibsInterfaceAllmanModel; 

 public:

  typedef ShellANDESModel    Self;
  typedef Model              Super;

  static const char*         DOF_NAMES[6];
  static const char*         SHAPE_PROP;
  static const char*         MATERIAL_PROP;
  static const char*         ALPHA_PROP;
  static const char*         BETA_PROP[10];
  static const char*         THICK_PROP;
  static const char*         ORIENTATION_PROP[3];
  static const char*	       LARGE_DISP_PROP;

                       ShellANDESModel
			 
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

  virtual              ~ShellANDESModel ();

  virtual void         getMatrix_

    ( Ref<MatrixBuilder>  mbuilder,
      const Vector&       force,
      const Vector&       disp )       const;

  virtual void         writeLocalOutput_

      ( const Vector&       disp,
        const Properties&   globdat );

  virtual void         getXOutTable_

      ( Ref<XTable>             table,
        const Vector&           weights,
        const String&           contents,
        const Vector&           disp );

  void                  getTransMatrix_

    ( Matrix&             transMatDof,
      Matrix&             transMat,
      const Matrix&       coords,
      const Vector&       v )   const;

  void                 get2DLocalcoordinates_

    ( Matrix&             xcoords,
      const Matrix&       transMat,
      const Matrix&       coords )     const;

  void                 getArea_

  ( double&             Area,
    const Matrix&       xcoords)    const;

  void                 getLmat_

  ( Matrix&             Lmat,
    const Matrix&       xcoords,
    const double&       alpha_ )    const;

  void                 getTthetau_

  ( Matrix&             Ttheta_u,
    const Matrix&       xcoords,
    const double&       Area )    const;

  void                 getGamma_

  ( double&             Gamma,
    const Matrix&       stiff)    const;

  void                 getKtheta_

  ( Matrix&             Ktheta,
    const Matrix&       xcoords,
    const Matrix&       D_mem,
    const Vector&       beta, 
    const double&       Area )    const;

  void                 getDmemNat_

    ( Matrix&             D_mem_nat,
      const Matrix&       xcoords,
      const Matrix&       D_mem,
      const double&       Area )    const;

  void                 getTmat_

  ( Matrix&             Tmat,
    const Matrix&       xcoords )    const;

  void                 getMembraneDisp_

  ( Vector&             Disp_mem,
    const Vector&       elemxDisp )   const;

  void                 getHmatAnalytic_

    ( Matrix&             Hmat,
      const Matrix&       Dplate,
      const Matrix&       xcoords )    const;

  void                 getHmatNumeric_

  ( Matrix&             Hmat,
    const Matrix&       Dplate,
    const Matrix&       xcoords,
    const Matrix&       ipcoords,
    const Vector&       weights )    const;

  void                 getBmat_

    ( Matrix&             Bmat,
        const Matrix&     Dplate,
        const Matrix&     xcoords )    const;

  void                 getDissForce_

  ( const Ref<Plasticity>  p,
    const Vector&       fstar,
    const Vector&       disp )    const;

  void                 initCharLength_ ();

  void                 getDissipation_

    ( const Properties&   params );

  void                 getIshape_

  ( Ref<IShape> shape ) const;

  void                 getMaterial_

  ( Ref<Material> material ) const;

  void                 getThickness_

  ( double thickness ) const;

  void                 getOrientVec_

  ( const double orientVec[3] ) const;

  void                 getStress_

    ( XTable&             table,
      const Vector&       weights,
      const Vector&       disp );

  // void                  getStrain_

  //   ( Vector&             strain,
  //     Matrix&             b,
  //     const Vector&       disp )    const;

  void                 writeElements_

  ( const Properties&   params,
    const Properties&   globdat );

  void                 writeNodes_

    ( const Properties&   params,
      const Properties&   globdat );
 
  void                 initLocOutWriter_ ();

  Ref<PrintWriter>     initWriter_

    ( const Properties&   params,
      const String        name ) const;

  void                 writeDisplacements_

    ( const Properties&   params,
      const Properties&   globdat );

  bool                 getTable_

  ( const Properties&   params,
    const Properties&   globdat );

 protected:

  Assignable<ElemGroup>   egroup_;
  Assignable<ElemSet>     elems_;
  Assignable<NodeSet>     nodes_;

  IdxVector               ielems_;

  idx_t                   nodeRank_;
  idx_t                   rank_;
  idx_t                   nodeCount_;
  idx_t                   numElem_;
  idx_t                   numNode_;
  idx_t                   strCount_;
  idx_t                   dofCount_;
  idx_t                   ipCount_;
  double                  alpha_;
  double                  beta_[10];
  double                  thickness_;
  double                  orientVec_[3];
  bool			              largeDisp_;

  Ref<IShape>             shape_;

  Ref<XDofSpace>          dofs_;
  IdxVector               dofTypes_;

  Ref<Material>           material_;
  Ref<Softening>          softening_;
  bool                    crackBandMethod_;
  
  Vector                  charLength_;
    
  Ref<PrintWriter>        dispOut_;
  Ref<PrintWriter>        elemOut_;
  Ref<PrintWriter>        damFOut_;
  Ref<PrintWriter>        damMOut_;

  ShapeGradsFunc          getShapeGrads_;
  ShapeFuncsFunc          getShapeFuncs_;

  String                  myTag_;

  static Ref<PrintWriter> nodeOut_;
  static idx_t            nodesWritten_;

  bool                    doLocalOutput_;
  IdxVector               locElems_;     // input
  IdxVector               locIpNumbers_; // input
  IdxVector               locIes_;       // relevant selection
  IdxVector               locIps_;       // relevant selection
  bool                    locDoAll_;
  Ref<PrintWriter>        locOut_;

  DispWriter              dw_;
};

#endif
