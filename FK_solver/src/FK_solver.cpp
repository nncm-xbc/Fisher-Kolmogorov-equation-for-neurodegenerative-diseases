#include "FK_solver.hpp"

void FK_solver::setup()
{
  // Create the mesh.
  {
    pcout << "Initializing the mesh" << std::endl;

    Triangulation<dim> mesh_serial;

    GridIn<dim> grid_in;
    grid_in.attach_triangulation(mesh_serial);

    std::ifstream grid_in_file(mesh_file_name);
    grid_in.read_msh(grid_in_file);

    GridTools::partition_triangulation(mpi_size, mesh_serial);
    const auto construction_data = TriangulationDescription::Utilities::
        create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
    mesh.create_triangulation(construction_data);

    pcout << "  Number of elements = " << mesh.n_global_active_cells()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the finite element space.
  {
    pcout << "Initializing the finite element space" << std::endl;

    fe = std::make_unique<FE_SimplexP<dim>>(r);

    pcout << "  Degree                     = " << fe->degree << std::endl;
    pcout << "  DoFs per cell              = " << fe->dofs_per_cell
          << std::endl;

    quadrature = std::make_unique<QGaussSimplex<dim>>(r + 1);

    pcout << "  Quadrature points per cell = " << quadrature->size()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the DoF handler.
  {
    pcout << "Initializing the DoF handler" << std::endl;

    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    pcout << "  Number of DoFs = " << dof_handler.n_dofs() << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the linear system.
  {
    pcout << "Initializing the linear system" << std::endl;

    pcout << "  Initializing the sparsity pattern" << std::endl;

    TrilinosWrappers::SparsityPattern sparsity(locally_owned_dofs,
                                               MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler, sparsity);
    sparsity.compress();

    pcout << "  Initializing the matrices" << std::endl;
    jacobian_matrix.reinit(sparsity);

    pcout << "  Initializing the system right-hand side" << std::endl;
    residual_vector.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    pcout << "  Initializing the solution vector" << std::endl;
    solution_owned.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    delta_owned.reinit(locally_owned_dofs, MPI_COMM_WORLD);

    solution.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    solution_old = solution;
  }
}

void FK_solver::assemble_system()
{
  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q = quadrature->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_residual(dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  jacobian_matrix = 0.0;
  residual_vector = 0.0;

  // Value and gradient of the solution on current cell.
  std::vector<double> solution_loc(n_q);
  std::vector<Tensor<1, dim>> solution_gradient_loc(n_q);

  // Value of the solution at previous timestep (un) on current cell.
  std::vector<double> solution_old_loc(n_q);

  forcing_term.set_time(time);

  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (!cell->is_locally_owned())
      continue;

    fe_values.reinit(cell);

    cell_matrix = 0.0;
    cell_residual = 0.0;

    fe_values.get_function_values(solution, solution_loc);
    fe_values.get_function_gradients(solution, solution_gradient_loc);
    fe_values.get_function_values(solution_old, solution_old_loc);

    for (unsigned int q = 0; q < n_q; ++q)
    {
      // Evaluate coefficients on this quadrature node.
      const double alpha_loc = alpha.value(fe_values.quadrature_point(q));

      // Define the matrix to store D values and populate it.
      // const FullMatrix<double> D_matrix = D.matrix_value(fe_values.quadrature_point(q));
      const Tensor<2, dim> D_matrix = D.matrix_value(fe_values.quadrature_point(q));

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
        {

          cell_matrix(i, j) += fe_values.shape_value(i, q) *
                               fe_values.shape_value(j, q) / deltat *
                               fe_values.JxW(q);

          cell_matrix(i, j) -= alpha_loc *
                               fe_values.shape_value(j, q) *
                               fe_values.shape_value(i, q) *
                               fe_values.JxW(q);

          cell_matrix(i, j) += 2.0 * alpha_loc *
                               solution_loc[q] *
                               fe_values.shape_value(j, q) *
                               fe_values.shape_value(i, q) *
                               fe_values.JxW(q);

          const Tensor<1, dim> &grad_phi_i = fe_values.shape_grad(i, q);
          const Tensor<1, dim> &grad_phi_j = fe_values.shape_grad(j, q);
          cell_matrix(i, j) += scalar_product(D_matrix * grad_phi_j, grad_phi_i) * fe_values.JxW(q);
        }

        // Time derivative term
        cell_residual(i) -= (solution_loc[q] - solution_old_loc[q]) /
                            deltat * fe_values.shape_value(i, q) *
                            fe_values.JxW(q);

        cell_residual(i) -= scalar_product(D_matrix * solution_gradient_loc[q], fe_values.shape_grad(i, q)) * fe_values.JxW(q);

        // Non-linear term
        cell_residual(i) += alpha_loc *
                            solution_loc[q] *
                            (1.0 - solution_loc[q]) *
                            fe_values.shape_value(i, q) *
                            fe_values.JxW(q);
      }
    }

    cell->get_dof_indices(dof_indices);

    jacobian_matrix.add(dof_indices, cell_matrix);
    residual_vector.add(dof_indices, cell_residual);
  }

  jacobian_matrix.compress(VectorOperation::add);
  residual_vector.compress(VectorOperation::add);
}

void FK_solver::solve_linear_system()
{
  SolverControl solver_control(100000, 1e-12 * residual_vector.l2_norm());

  SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);
  TrilinosWrappers::PreconditionSSOR preconditioner;
  preconditioner.initialize(
      jacobian_matrix, TrilinosWrappers::PreconditionSSOR::AdditionalData(1.0));

  solver.solve(jacobian_matrix, delta_owned, residual_vector, preconditioner);
  pcout << "  " << solver_control.last_step() << " CG iterations" << std::endl;
}

void FK_solver::solve_newton()
{
  const unsigned int n_max_iters = 1000;
  const double residual_tolerance = 1e-6;

  unsigned int n_iter = 0;
  double residual_norm = residual_tolerance + 1;

  while (n_iter < n_max_iters && residual_norm > residual_tolerance)
  {
    assemble_system();
    residual_norm = residual_vector.l2_norm();

    pcout << "  Newton iteration " << n_iter << "/" << n_max_iters
          << " - ||r|| = " << std::scientific << std::setprecision(6)
          << residual_norm << std::flush;

    // We actually solve the system only if the residual is larger than the
    // tolerance.
    if (residual_norm > residual_tolerance)
    {
      solve_linear_system();
      solution_owned += delta_owned;
      solution = solution_owned;
    }
    else
    {
      pcout << " < tolerance" << std::endl;
    }

    ++n_iter;
  }
}

void FK_solver::output(const unsigned int &time_step) const
{
  DataOut<dim> data_out;
  data_out.add_data_vector(dof_handler, solution, "u");

  std::vector<unsigned int> partition_int(mesh.n_active_cells());
  GridTools::get_subdomain_association(mesh, partition_int);
  const Vector<double> partitioning(partition_int.begin(), partition_int.end());
  data_out.add_data_vector(partitioning, "partitioning");

  data_out.build_patches();

  data_out.write_vtu_with_pvtu_record(
      "./", "output", time_step, MPI_COMM_WORLD, 3);
}

void FK_solver::solve()
{
  pcout << "===============================================" << std::endl;

  time = 0.0;

  // Apply the initial condition.
  {
    pcout << "Applying the initial condition" << std::endl;

    VectorTools::interpolate(dof_handler, u_0, solution_owned);
    solution = solution_owned;

    // Output the initial solution.
    output(0);
    pcout << "-----------------------------------------------" << std::endl;
  }

  unsigned int time_step = 0;

  while (time < T - 0.5 * deltat)
  {
    time += deltat;
    ++time_step;

    // Store the old solution, so that it is available for assembly.
    solution_old = solution;

    pcout << "n = " << std::setw(3) << time_step << ", t = " << std::setw(5)
          << std::fixed << time << std::endl;

    // At every time step, we invoke Newton's method to solve the non-linear
    // problem.
    solve_newton();

    output(time_step);

    pcout << std::endl;
  }
}

double
FK_solver::compute_error(const VectorTools::NormType &norm_type)
{
  FE_SimplexP<dim> fe_linear(1);
  MappingFE mapping(fe_linear);

  const QGaussSimplex<dim> quadrature_error = QGaussSimplex<dim>(r + 2);

  exact_solution.set_time(time);

  Vector<double> error_per_cell;
  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    solution,
                                    exact_solution,
                                    error_per_cell,
                                    quadrature_error,
                                    norm_type);

  const double error =
      VectorTools::compute_global_error(mesh, error_per_cell, norm_type);

  return error;
}