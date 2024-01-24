#include "HeatNonLinear.hpp"

// Main function.
int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const unsigned int degree = 2;

  const double T      = 1;
  const double deltat = 0.01;

  HeatNonLinear problem("../mesh/brain-h3.0.msh", degree, T, deltat);

  problem.setup();
  problem.solve();

  return 0;
}