log =
{
  pattern = "*.info";
  file    = "-$(CASE_NAME).log";
};

control = 
{
  pause = 0.;
  runWhile = "i < 100";
};

userInput =
{
  modules = [ "GmshInput", "ngroups", "loads" ];

  GmshInput =
  {
    type = "GmshInput";
    file = "twistingmoments.msh";

    doElemGroups = true; // don't make element groups from physical surfaces
  };

  // define NodeGroups for boundary conditions

  ngroups = 
  {
    type = "GroupInput";
    nodeGroups = [ "left", "RT", "RB"];
    left = 
    {
      xtype = "min";
    };
    RT = 
    {
      xtype = "max";
      ytype = "max";
    };
    RB = 
    {
      xtype = "max";
      ytype = "min";
    };
  };
  
  // define load tables
  
  loads = 
  {
    type = "Input";
    file = "twistingmoments.data";
  };
};

model =
{
  type        = "Matrix";
  matrix.type = "Sparse"; 

  model       =
  {
    type   = "Multi";

    models = [ "plate" ,        // first ShellFEMModel
               "arclen" ,       // DispArclenModel
               "lodi" ];        // LoadDispModel

    // first ShellFEMModel: Cantilever Beam

    plate =
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
        young1  = 10000000.0;
        young2  = 10000000.0;
        poisson12 = 0.25;
        shear12 = 4000000.0;
        theta      = 0.;
      };

      // Orientation vector and thickness

      v1 = 1.0;
      v2 = 0.0;
      v3 = 0.0;

      alpha = 1.5;
      beta = 0.5;

      thickness = 0.05;
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

      dispIncr = 0.1;

      // bounds for increment size
      
      minDispIncr = 0.00001;
      maxDisp	= 100.0;

      minIncr = 0.0001; // (energy)
      maxIncr = 2.;

      // define constraints for DispArclenBCs

      constraints =
      {
        nodeGroups = [ "left", "left", "left", "left", "left" ];
        dofs = [ "u", "v", "w", "rx", "ry" ];

        loads = 
        {
          loadTable = "load";
        };        
      };
    };

    lodi =
    {
      type = "LoadDisp";
      group = "RT";
    };
  };
};

sample =
{
    file = "lodi.dat";

    // only write final solution for time step (no further crack growth)

    sampleWhen = "accepted";

    // the data that is written (see 'lodi' above)

    dataSets = [ "i","model.model.lodi.disp[2]",
                 "abs(model.model.lodi.load[2])" ];
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
      precision = 1.0e-1;
      solver.type = "Skyline";
      solver.useThreads = true;
      maxIter = 10;
    };

    // XArclenModule for dissipation increments

    arcLen =
    {
      allowDespair = true;  // allow modified Newton-Raphson scheme
      precision = 1.0e-1;
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
