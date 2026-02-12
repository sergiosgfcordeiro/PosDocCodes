log =
{
  pattern = "*.info";
  file    = "-$(CASE_NAME).log";
};

control = 
{
  pause = 0.;
  runWhile = "i < 2";
};

userInput =
{
  modules = [ "GmshInput", "ngroups", "loads" ];

  GmshInput =
  {
    type = "GmshInput";
    file = "cooksmembrane.msh";

    doElemGroups = true; // don't make element groups from physical surfaces
  };

  // define NodeGroups for boundary conditions

  ngroups = 
  {
    type = "GroupInput";
    nodeGroups = [ "left", "right", "bottom",  "top"];
    left = 
    {
      xtype = "min";
    };
    right = 
    {
      xtype = "max";
    };
    bottom = 
    {
      ytype = "min";
    };
    top = 
    {
      ytype = "max";
    };
  };
  
  // define load tables
  
  loads = 
  {
    type = "Input";
    file = "cooks.data";
  };
};

model =
{
  type        = "Matrix";
  matrix.type = "Sparse"; 

  model       =
  {
    type   = "Multi";

    models = [ "cooksmembrane" ,        // first ShellFEMModel
               "arclen" ,       // DispArclenModel
               "lodi" ];        // LoadDispModel

    // first ShellFEMModel: Cantilever Beam

    cooksmembrane =
    {
      type      = "Shell6DofsBF";
      elements  = "gmsh1";  // "gmsh0" is the default element group for panels

      shape =
      {
        type      = "Triangle3";
        intScheme = "Gauss1"   ;
        shapeFuncs =
        {
          type = "Linear"; 
        };
      };

      // Elasticity parameters

      material =
      {
        type    = "Orthotropic";
        dim     = 2;
        state   = "PLANE_STRESS";
        young1  = 1.0;
        young2  = 1.0;
        poisson12 = 0.333333;
        shear12 = 0.375;
        theta      = 0.;
      };

      // Orientation vector and thickness

      v1 = 1.0;
      v2 = 0.0;
      v3 = 0.0;

      alpha = 1.5;
      beta = 0.5;

      thickness = 1.0;
    };

    arclen = 
    {
      type = "DispArclen";

      // thresholds for switching to dissipation increments

      swtIter = 4;   
      swtEnergy = 1.;

      // optimal number of iterations (for adaptive increments)
      
      optIter = 6;   

      // factor for increment reduction

      reduction = .2;

      // initial displacement increment:

      dispIncr = 1.0;

      // bounds for increment size
      
      minDispIncr = 0.00001;
      maxDisp	= 100.0;

      minIncr = 0.0001; // (energy)
      maxIncr = 2.;

      // define constraints for DispArclenBCs

      constraints =
      {
        nodeGroups = [ "left", "left", "left" ];
        dofs = [ "u", "v", "rz" ];

        loads = 
        {
          loadTable = "load";
        };        
      };
    };

    lodi =
    {
      type = "LoadDisp";
      group = "right";
    };
  };
};

sample =
{
    file = "lodi.dat";

    // only write final solution for time step (no further crack growth)

    sampleWhen = "accepted";

    // the data that is written (see 'lodi' above)

    dataSets = [ "i","model.model.lodi.disp[1]",
                 "abs(model.model.lodi.load[1])" ];
};
 
userModules = 
{
  modules = [ "solver"  ];

  solver =
  {
    type = "FlexArclen";

    nonLin = 
    {
      lineSearch = false;
      precision = 1.0e-6;
      solver.type = "Skyline";
      solver.useThreads = true;
      maxIter = 10;
    };

    // XArclenModule for dissipation increments

    arcLen =
    {
      allowDespair = true;  // allow modified Newton-Raphson scheme
      precision = 1.0e-6;
      maxIter = 10;
    };

    arcLen.solver = nonLin.solver;
  };
};

// Parameters for the userOutputModule.

userOutput =
{
modules = [ "paraview" ];

// Parameters controlling the paraview output.

  paraview =
  {
    type = "ParaView";
    output_format = "$(CASE_NAME)/vis%i";
    output_file = "vtu";
    sampleWhen = "(i-1)%1<1";
    groups = [ "gmsh1" ];
   
    gmsh1 =
    {
      shape = "Triangle3";
      disps = [ "u", "v", "w" ];
      // otherDofs = [ "rx", "ry", "rz" ];
    };
  };
};
