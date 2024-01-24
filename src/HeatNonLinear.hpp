#ifndef HEAT_NON_LINEAR_HPP
#define HEAT_NON_LINEAR_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_in.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;

// Class representing the non-linear diffusion problem.
class HeatNonLinear
{
public:
  // Physical dimension (1D, 2D, 3D)
  static constexpr unsigned int dim = 3;

  // Function for the mu_0 coefficient.
  // class FunctionD : public Function<dim>
  // {
  // public:
  //   virtual double
  //   value(const Point<dim> & /*p*/,
  //         const unsigned int /*component*/ = 0) const override
  //   {
  //     return 1.0;
  //   }
  // };

  // Function for the mu_1 coefficient.
  class FunctionAlpha : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 1.0;
    }
  };

  // Function for the forcing term.
  class ForcingTerm : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
          {
          return 0.0;

  }

  };

  // // Function for Dirichlet boundary conditions.
  // class FunctionG : public Function<dim>
  // {
  // public:
  //   virtual double
  //   value(const Point<dim> & /*p*/,
  //         const unsigned int /*component*/ = 0) const override
  //   {
  //     return 0.0;
  //   }
  // };

  // Function for initial conditions.
  class FunctionU0 : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> & p,
          const unsigned int /*component*/ = 0) const override
    {
      double x = std::sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
      return (x*x)/(x*x+1);
    }
  };
    //  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

      // Class representing the matrix function D.
  //      class ForcingTerm : public Function<dim>
  // {
  // public:
  //   virtual void
  //   vector_value(const Point<dim> & /*p*/,
  //                Vector<double> &values) const override
  //   {
  //     for (unsigned int i = 0; i < dim - 1; ++i)
  //       values[i] = 0.0;

  //     values[dim - 1] = -g;
  //   }

  //   virtual double
  //   value(const Point<dim> & /*p*/,
  //         const unsigned int component = 0) const override
  //   {
  //     if (component == dim - 1)
  //       return -g;
  //     else
  //       return 0.0;
  //   }

  // protected:
  //   const double g = 0.0;
  // };

  class FunctionD : public Function<dim>{
    public:
     virtual void
     matrix_value(const Point<dim> & /*p*/,
                  FullMatrix<double> &values) const //override
     {
       for (unsigned int i = 0; i < dim; ++i){
        for (unsigned int j = 0; j < dim; ++j){
          if(i==j)
            values[i][j] = 1.0;
          else
            values[i][j] = 0.0;
        }
       }
     }
     virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int component_1 = 0, const unsigned int component_2 = 1) const //override
    {
      if (component_1 == component_2)
        return 1.0;
      else
        return 0.0;
    }
  };
  

  // Constructor. We provide the final time, time step Delta t and theta method
  // parameter as constructor arguments.
  HeatNonLinear(const std::string  &mesh_file_name_,
                const unsigned int &r_,
                const double       &T_,
                const double       &deltat_)
    : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , pcout(std::cout, mpi_rank == 0)
    , T(T_)
    , mesh_file_name(mesh_file_name_)
    , r(r_)
    , deltat(deltat_)
    , mesh(MPI_COMM_WORLD)
  {}

  // Initialization.
  void
  setup();

  // Solve the problem.
  void
  solve();

protected:
  // Assemble the tangent problem.
  void
  assemble_system();

  // Solve the linear system associated to the tangent problem.
  void
  solve_linear_system();

  // Solve the problem for one time step using Newton's method.
  void
  solve_newton();

  // Output.
  void
  output(const unsigned int &time_step) const;

  // MPI parallel. /////////////////////////////////////////////////////////////

  // Number of MPI processes.
  const unsigned int mpi_size;

  // This MPI process.
  const unsigned int mpi_rank;

  // Parallel output stream.
  ConditionalOStream pcout;

  // Problem definition. ///////////////////////////////////////////////////////

  // mu_0 coefficient.
  FunctionD D;

  // mu_1 coefficient.
  FunctionAlpha alpha;

  // Forcing term.
  ForcingTerm forcing_term;

  // // Dirichlet boundary conditions.
  // FunctionG function_g;

  // Initial conditions.
  FunctionU0 u_0;

  // Current time.
  double time;

  // Final time.
  const double T;

  // Discretization. ///////////////////////////////////////////////////////////

  // Mesh file name.
  const std::string mesh_file_name;

  // Polynomial degree.
  const unsigned int r;

  // Time step.
  const double deltat;

  // Mesh.
  parallel::fullydistributed::Triangulation<dim> mesh;

  // Finite element space.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> quadrature;

  // DoF handler.
  DoFHandler<dim> dof_handler;

  // DoFs owned by current process.
  IndexSet locally_owned_dofs;

  // DoFs relevant to the current process (including ghost DoFs).
  IndexSet locally_relevant_dofs;

  // Jacobian matrix.
  TrilinosWrappers::SparseMatrix jacobian_matrix;

  // Residual vector.
  TrilinosWrappers::MPI::Vector residual_vector;

  // Increment of the solution between Newton iterations.
  TrilinosWrappers::MPI::Vector delta_owned;

  // System solution (without ghost elements).
  TrilinosWrappers::MPI::Vector solution_owned;

  // System solution (including ghost elements).
  TrilinosWrappers::MPI::Vector solution;

  // System solution at previous time step.
  TrilinosWrappers::MPI::Vector solution_old;
};

#endif