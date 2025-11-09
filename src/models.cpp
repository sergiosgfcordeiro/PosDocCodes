
#include "models.h"


//-----------------------------------------------------------------------
//   declareModels
//-----------------------------------------------------------------------


void declareModels ()
{
  declareDirichletModel();
  declareDispArclenModel();
  declareStructuralInterfaceModel();
  declarePanelRibsInterfaceModel();
  declarePanelRibsInterface6DofsBFModel();
  declarePanelRibsInterface6DofsAllmanModel();
  declareLoadDispModel();
  declareNeumannModel();
  declareShellModel();
  declareShell6DofsBFModel();
  declareShell6DofsAllmanModel();
  declareShellANDESModel();
}


