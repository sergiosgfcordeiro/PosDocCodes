log =
{
  pattern = "*.info";
  file    = "-$(CASE_NAME).log";
};

control = 
{
  pause = 0.;
  runWhile = "i < 150";
};

userInput =
{
  modules = [ "GmshInput", "reorder", "ngroups" ];

  GmshInput =
  {
    type = "GmshInput";
    file = "panelribs.msh";

    doElemGroups = true; // don't make element groups from physical surfaces
  };

  // reorderes element numbering from that for a Quad4 element to that for a BLine2 element 

  reorder =
  {
    type = "Quad4Interf";
  
    elemGroups = ["gmsh2"];
    thickness = 1.0;
  };

  // define NodeGroups for boundary conditions

  ngroups = 
  {
    type = "GroupInput";
    nodeGroups = [ "panel", "ribs",  "rib_top"];
    panel = 
    {
      restrictToGroups = "gmsh0";
      xbounds = [0.0, 60.0];
      ybounds = [0.0, 60.0];
    };
    ribs = 
    {
      restrictToGroups = "gmsh1";
      zbounds = [0.0, 30.5];
    };
    rib_top = 
    {
      restrictToGroups = "gmsh1";
      ztype = "max";
    };
  };
};

model =
{
  type        = "Matrix";
  matrix.type = "Sparse"; 

  model       =
  {
    type   = "Multi";

    models = [ "panel" ,        // first ShellFEMModel
               "ribs" ,         // second ShellFEMModel
               "interface1" ,   // firts PanelRibsInterfaceModel
               "diri" ,       // DispArclenModel
               "lodi" ];        // LoadDispModel

    // first ShellFEMModel: Panel-interface 1a

    panel =
    {
      type      = "Shell6DofsBF";
      elements  = "gmsh0";  // "gmsh0" is the default element group for panels

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
        young1  = 139.4e3;
        young2  = 10.16e3;
        poisson12 = 0.30;
        shear12 = 4.6e3;
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

    // second ShellFEMModel: Ribs
    
    ribs =
    {
      type      = "Shell6DofsBF";
      elements  = "gmsh1";  // "gmsh4" is the defaultelement group for ribsfree

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
        young1  = 15.5e3;
        young2  = 15.5e3;
        poisson12 = 0.3;
        shear12 = 5.96154e3;
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

    // first PanelRibsInterfaceModel: interface1

    interface1 = 
    {
      type = "PanelRibsInterface6DofsBFModel";
      elements  = "gmsh2";  // this element group should be created in the future by a PanelRibsShellMeshModule
     
      shape =
      {
        type      = "BLine2";
        intScheme = "Gauss13";
        shapeFuncs =
        {
          type = "Linear"; 
        };
      };

      thick_shape =
      {
        type      = "Line2";
        intScheme = "Gauss7";
      };

      panel_name = "panel";
      rib_name = "ribs";
    
      // Cohesive zone properties
    
      coheMat =
      {
        type   = "Turon";
        dim    = 3;
        // dummy  = 10500.0e3;
        dummy  = 161.78e3;
    
        f2t = 60.;
        f6  = 60.;
        gI  = 0.170;
        gII = 0.494;
        eta = 1.62;
      };
    };

    diri = 
    {
      type = "Dirichlet";

      dispIncr = 0.05;

      nodeGroups = [ "rib_top", "panel", "panel", "panel", "panel", "panel", "panel" ];
      dofs = [ "v", "u", "v", "w", "wx", "wy", "wz" ];
      factors = [1., 0., 0., 0., 0., 0., 0.]; 
    };

    lodi =
    {
      type = "LoadDisp";
      group = "rib_top";
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
  modules = [ "solver", "graph"  ];

  solver =
  {
    type = "AdaptiveStep";

    nonlin = 
    {
      lineSearch = false;
      precision = 1.0e-6;
      solver.type = "Skyline";
      solver.useThreads = true;
      tiny = 1.e-10;
      maxIter = 10;
    };

  };

  graph =
  {
    type = "Graph";

    sampleWhen = "accepted";
    
    dataSets = "loadDisp";
    loadDisp =
    {
      xData = "model.model.lodi.disp[1]";
      yData = "abs(model.model.lodi.load[1])";
    };
    graph.keyPos  = [ 0.060000, 45.00000 ];
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
    groups = [ "gmsh0", "gmsh1", "gmsh2" ];
   
    gmsh0 =
    {
      shape = "Triangle3";
      disps = [ "u", "v", "w" ];
      // otherDofs = [ "wx", "wy" ];
    };
    gmsh1 =
    {
      shape = "Triangle3";
      disps = [ "u", "v", "w" ];
      // otherDofs = [ "wx", "wy" ];
    };
    gmsh2 =
    {
      shape = "BLine2";
      disps = [ "u", "v", "w" ];
      // otherDofs = [ "wx", "wy" ];
    };
  };
};
