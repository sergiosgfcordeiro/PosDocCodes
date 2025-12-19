log =
{
  pattern = "*.info";
  file    = "-$(CASE_NAME).log";
};

control = 
{
  pause = 0.;
  runWhile = "i < 3000";
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
  
    elemGroups = ["gmsh3"];
    thickness = 1.0;
  };

  // define NodeGroups for boundary conditions

  ngroups = 
  {
    type = "GroupInput";
    nodeGroups = [ "left", "RB", "RT"];
    left = 
    {
      restrictToGroups = "gmsh1";
      xtype = "min";
    };
    RB = 
    {
      restrictToGroups = "gmsh1";
      xtype = "max";
    };
    RT = 
    {
      restrictToGroups = "gmsh2";
      xtype = "max";
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
               "arclen" ,       // DispArclenModel
               "lodi" ];        // LoadDispModel

    // first ShellFEMModel: Panel-interface 1a

    panel =
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
      elements  = "gmsh2";  // "gmsh4" is the defaultelement group for ribsfree

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
      elements  = "gmsh3";  // this element group should be created in the future by a PanelRibsShellMeshModule
     
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
        intScheme = "Gauss5";
      };

      panel_name = "panel";
      rib_name = "ribs";
    
      // Cohesive zone properties
    
      coheMat =
      {
        type   = "Turon";
        dim    = 3;
        // dummy  = 10500.0e3;
        dummy  = 161.78e2;
    
        f2t = 60.;
        f6  = 60.;
        gI  = 0.170;
        gII = 0.494;
        eta = 1.62;
      };
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

      dispIncr = 0.001;

      // bounds for increment size
      
      minDispIncr = 0.00001;
      maxDisp	= 100.0;

      minIncr = 0.0001; // (energy)
      maxIncr = 10.;

      // define constraints for DispArclenBCs

      constraints =
      {
        nodeGroups = [ "RT", "RT", "left", "RB", "RB" ];
        dofs = [ "w", "u", "u", "u", "w" ];
        loaded = 0;  // RT.w is the loaded boundary

        // instead of identifying nonzero dirichlet bc's with loaded, 
        // neumann bc's could be indicated here with a loadTable
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
  modules = [ "solver", "graph"  ];

  solver =
  {
    type = "FlexArclen";

    nonLin = 
    {
      lineSearch = false;
      precision = 1.0e-4;
      solver.type = "Skyline";
      solver.useThreads = true;
      maxIter = 20;
    };

    // XArclenModule for dissipation increments

    arcLen =
    {
      allowDespair = true;  // allow modified Newton-Raphson scheme
      precision = 1.0e-4;
      maxIter = 20;
    };

    arcLen.solver = nonLin.solver;
  };

  graph =
  {
    type = "Graph";

    sampleWhen = "accepted";
    
    dataSets = "loadDisp";
    loadDisp =
    {
      xData = "model.model.lodi.disp[2]";
      yData = "abs(model.model.lodi.load[2])";
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
    groups = [ "gmsh1", "gmsh2", "gmsh3" ];
   
    gmsh1 =
    {
      shape = "Triangle3";
      disps = [ "u", "v", "w" ];
      // otherDofs = [ "wx", "wy" ];
    };
    gmsh2 =
    {
      shape = "Triangle3";
      disps = [ "u", "v", "w" ];
      // otherDofs = [ "wx", "wy" ];
    };
    gmsh3 =
    {
      shape = "BLine2";
      disps = [ "u", "v", "w" ];
      // otherDofs = [ "wx", "wy" ];
    };
  };
};
