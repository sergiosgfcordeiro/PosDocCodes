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
    file = "cantileverbeam.msh";

    doElemGroups = true; // don't make element groups from physical surfaces
  };

  // define NodeGroups for boundary conditions

  ngroups = 
  {
    type = "GroupInput";
    nodeGroups = [ "left", "R0", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10", "R11", "R12", "R13", "R14", "R15", "R16", "bottom",  "top"];
    left = 
    {
      xtype = "min";
    };
    R0 = 
    {
      xtype = "max";
      ybounds = [-0.01, 0.01];
    };
    R1 = 
    {
      xtype = "max";
      ybounds = [0.74, 0.76];
    };
    R2 = 
    {
      xtype = "max";
      ybounds = [1.49, 1.51];
    };
    R3 = 
    {
      xtype = "max";
      ybounds = [2.24, 2.26];
    };
    R4 = 
    {
      xtype = "max";
      ybounds = [2.99, 3.01];
    };
    R5 = 
    {
      xtype = "max";
      ybounds = [3.74, 3.76];
    };
    R6 = 
    {
      xtype = "max";
      ybounds = [4.49, 4.51];
    };
    R7 = 
    {
      xtype = "max";
      ybounds = [5.24, 5.26];
    };
    R8 = 
    {
      xtype = "max";
      ybounds = [5.99, 6.01];
    };
    R9 = 
    {
      xtype = "max";
      ybounds = [6.74, 6.76];
    };
    R10 = 
    {
      xtype = "max";
      ybounds = [7.49, 7.51];
    };
    R11 = 
    {
      xtype = "max";
      ybounds = [8.24, 8.26];
    };
    R12 = 
    {
      xtype = "max";
      ybounds = [8.99, 9.01];
    };
    R13 = 
    {
      xtype = "max";
      ybounds = [9.74, 9.76];
    };
    R14 = 
    {
      xtype = "max";
      ybounds = [10.49, 10.51];
    };
    R15 = 
    {
      xtype = "max";
      ybounds = [11.24, 11.26];
    };
    R16 = 
    {
      xtype = "max";
      ybounds = [11.99, 12.01];
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
    file = "cantilever.data";
  };
};

model =
{
  type        = "Matrix";
  matrix.type = "Sparse"; 

  model       =
  {
    type   = "Multi";

    models = [ "cantileverbeam" ,        // first ShellFEMModel
               "arclen" ,       // DispArclenModel
               "lodi" ];        // LoadDispModel

    // first ShellFEMModel: Cantilever Beam

    cantileverbeam =
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
        young1  = 30000.0;
        young2  = 30000.0;
        poisson12 = 0.25;
        shear12 = 12000.0;
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
        nodeGroups = [ "left", "left", "left", "left", "left", "left" ];
        dofs = [ "u", "v", "w", "wx", "wy", "wz" ];

        loads = 
        {
          loadTable = "load";
        };        
      };
    };

    lodi =
    {
      type = "LoadDisp";
      group = "R8";
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
      // otherDofs = [ "wx", "wy" ];
    };
  };
};
