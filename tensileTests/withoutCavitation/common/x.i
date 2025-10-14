RVE_length =  50.0
strain = 0.01
erate = 0.29988 # in hr. = 8.33e-5/s
dispRate = '${fparse erate*RVE_length}'
endTime = '${fparse strain/erate}'

[Mesh]
  [base]
    type = FileMeshGenerator
  []
  [rename]
    type = RenameBoundaryGenerator
    input = base
    old_boundary = '1 2 3 4 5 6' 
    new_boundary = 'x0 x1 y0 y1 z0 z1'
  []
  use_displaced_mesh = false
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Modules]
  [TensorMechanics]
    [Master]
      [all]
        strain = FINITE
        new_system = true
        formulation = TOTAL
        add_variables = true
        volumetric_locking_correction = true
        generate_output = 'cauchy_stress_xx cauchy_stress_yy cauchy_stress_zz cauchy_stress_yz cauchy_stress_xz cauchy_stress_xy '
                          'mechanical_strain_xx mechanical_strain_yy mechanical_strain_zz mechanical_strain_yz mechanical_strain_xz mechanical_strain_xy'
      []
    []
  []
[]

[UserObjects]
  [./euler_angle_file]
    type = PropertyReadFile
    nprop = 3
    # read_type = block
    # nblock = 385
    use_zero_based_block_indexing = false
  [../]
[]

[Materials]
  [stress]
  # define the bulk material model, euler angles for each grain come from the `euler_angle_file` UserObjects
    type = NEMLCrystalPlasticity
    model = "cpdeformation"
    large_kinematics = true
    euler_angle_reader = euler_angle_file
  []
[]

[BCs]
  [x0]
    type = DirichletBC
    variable = disp_x
    boundary = x0
    value = 0.0
  []
  [y0]
    type = DirichletBC
    variable = disp_y
    boundary = y0
    value = 0.0
  []
  [z0]
    type = DirichletBC
    variable = disp_z
    boundary = z0
    value = 0.0
  []
  [x1]
    type = FunctionDirichletBC
    boundary = x1
    function = tdisp
    variable = disp_x
  []
[]

[Functions]
  [./tdisp]
    type = ParsedFunction
    expression = '${dispRate} * t'
  [../]
[]

[Constraints]
  [x1]
    type = EqualValueBoundaryConstraint
    variable = disp_x
    secondary = 'x1'
    penalty = 1e7
  []
  [y1]
    type = EqualValueBoundaryConstraint
    variable = disp_y
    secondary = 'y1'
    penalty = 1e7
  []
  [z1]
    type = EqualValueBoundaryConstraint
    variable = disp_z
    secondary = 'z1'
    penalty = 1e7
  []
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [cauchy_stress_xx]
    type = ElementAverageValue
    variable = cauchy_stress_xx
  []
  [avg_disp_x]
    type = SideAverageValue
    variable = disp_x
    boundary = x1
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [strain]
    type = ParsedPostprocessor
    pp_names = 'avg_disp_x'
    function = 'avg_disp_x / ${RVE_length}'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [delta_strain]
    type = ChangeOverTimePostprocessor
    postprocessor = strain
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [dt]
    type = TimestepSize
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []
  [strain_rate]
    type = ParsedPostprocessor
    pp_names = 'delta_strain dt'
    function = 'delta_strain / dt'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient

  solve_type = 'newton'

  petsc_options = '-snes_converged_reason -ksp_converged_reason'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_interp_type -pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_agg_nl -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_truncfactor'
  petsc_options_value = 'hypre boomeramg 301 0.7 ext+i HMIS 4 2 0.4'
  

  line_search = none
  automatic_scaling = true
  l_max_its = 300
  # l_tol = 1e-7
  nl_max_its = 15
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  nl_forced_its = 1
  n_max_nonlinear_pingpong = 1
  dtmin = 1e-8 #${dtmin}
  dtmax = 1e2 #${dtmax}
  end_time = '${endTime}'
  
  [./Predictor]
    type = SimplePredictor
    scale = 1.0
    skip_after_failed_timestep = true
  [../]
  
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-4
    growth_factor = 2
    cutback_factor = 0.5
    cutback_factor_at_failure = 0.1
    optimal_iterations = 8
    iteration_window = 1
    linear_iteration_ratio = 1000000000
  []
[]

[Outputs]
  print_linear_residuals = false
  checkpoint = false
  [./out_csv]
    type = CSV
    file_base = amTensile_out_x
  [../]
[]


