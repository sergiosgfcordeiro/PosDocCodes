
#include "modules.h"

#include <jive/fem/InputModule.h>
#include <jive/app/ModuleFactory.h>

#include "ParaViewModule.h"


//-----------------------------------------------------------------------
//   declareModules
//-----------------------------------------------------------------------


void declareModules ()
{
  declareAdaptiveStepModule      ();
  declareFlexArclenModule        ();
  declareGmshInputModule         ();
  declareLaminateShellMeshModule ();
  // declarePanelRibsMeshModule     ();
  declareGroupInputModule        ();
  declareInputModule             ();
  declareQuad4InterfModule       ();
  declareXOutputModule           ();
  ParaViewModule::declare        ();
}

void declareInputModule ()
{
  using jive::app::ModuleFactory;
  using jive::fem::InputModule;

  ModuleFactory::declare ( "Input",
                         & InputModule::makeNew );
  ModuleFactory::declare ( "Output",
                          & InputModule::makeNew );
}

