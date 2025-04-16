// Implement member functions for ScalarPoissonProblem class.
#include "scalar-poisson-problem.h"
#include "geometrycentral/numerical/linear_solvers.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
ScalarPoissonProblem::ScalarPoissonProblem(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;

    // TODO: Build member variables A (Laplace matrix), M (mass matrix), total area
    this->A = inputGeo->laplaceMatrix(); 
    this->M = inputGeo->massMatrix(); 
    this->totalArea = inputGeo->totalArea(); 
}

/*
 * Computes the solution of the poisson problem Ax = -M(rho - rhoBar), where A is the POSITIVE DEFINITE Laplace matrix
 * and M is the mass matrix.
 *
 * Input: <rho>, the density of vertices in the mesh.
 * Returns: The solution vector.
 */
Vector<double> ScalarPoissonProblem::solve(const Vector<double>& rho) const {

    // TODO
    // Note: Geometry Central has linear solvers: https://geometry-central.net/numerical/linear_solvers/
    SparseMatrix<double> L = this->A ;
    // solve Lx = -M(rho - rhoBar)

    // rhoBar = /int_M (rho) / totalArea
    double rhoBar_value = 0;
    for (int i = 0; i < rho.rows(); i++) {
        rhoBar_value += rho[i] * geometry->barycentricDualArea(this->mesh->vertex(i));
    }
    rhoBar_value /= this->totalArea;

    // solve Lx = -M(rho - rhoBar)

    Vector<double> rhoBar = Vector<double>::Constant(rho.rows(), rhoBar_value);

    Vector<double> rhs = -M * (rho - rhoBar);

    Vector<double> x = geometrycentral::solvePositiveDefinite(L, rhs);

    return x;
}