#include "heat-method.h"
#include "GeometryCentral/numerical/linear_solvers.h"



using namespace geometrycentral;
using namespace geometrycentral::surface;

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
HeatMethod::HeatMethod(ManifoldSurfaceMesh* surfaceMesh, VertexPositionGeometry* geo) {

    this->mesh = surfaceMesh;
    this->geometry = geo;

    // TODO: Build Laplace and flow matrices.
    // Note: core/geometry.cpp has meanEdgeLength() function
    double t = geometry->meanEdgeLength();
    t = t * t;


    this->A = geometry->laplaceMatrix();
    this->F = geometry->massMatrix() + (t * A);
}

/*
 * Computes the vector field X = -∇u / |∇u|.
 *
 * Input: <u>, a dense vector representing the heat that is allowed to diffuse on the input mesh for a brief period of
 * time.
 * Returns: A MeshData container that stores a Vector3 per face.
 */
FaceData<Vector3> HeatMethod::computeVectorField(const Vector<double>& u) const {

    FaceData<Vector3> vectorField(*mesh, Vector3{0, 0, 0});
    

    for (Face f : mesh->faces()) {
        Vector3 gradient = Vector3::zero();
        Vector3 normal = geometry->faceNormal(f);
        
        for (Halfedge he : f.adjacentHalfedges()) {
            size_t vertexIndex = he.vertex().getIndex();
            double heatValue = u[vertexIndex];
            
            // get Mathcal{J} e_i
            Vector3 edgeVector = geometry->halfedgeVector(he.next());
            Vector3 perpVector = cross(normal, edgeVector);
            
            gradient += perpVector * heatValue;
        }
        vectorField[f] = -gradient.normalizeCutoff();
    }
    
    return vectorField;
}

/*
 * Computes the integrated divergence ∇.X.
 *
 * 
 * 
 * Input: <X>, the vector field -∇u / |∇u| represented as a FaceData container
 * Returns: A dense vector containing the divergence at each vertex
 */
Vector<double> HeatMethod::computeDivergence(const FaceData<Vector3>& X) const {
    // Initialize divergence vector with zeros for each vertex
    Vector<double> divergence = Vector<double>::Zero(mesh->nVertices());
    
    for (Face f : mesh->faces()) {
        Vector3 faceVector = X[f];
        
        // For each halfedge in the face, compute its contribution to the divergence
        for (Halfedge he : f.adjacentHalfedges()) {
            Vector3 edgeVector = geometry->halfedgeVector(he);
            
            double cotWeight = geometry->cotan(he);
            double contribution = 0.5 * cotWeight * dot(edgeVector, faceVector);
            
            size_t tailIndex = he.tailVertex().getIndex();
            size_t tipIndex  = he.tipVertex().getIndex();
            
            divergence[tailIndex] += contribution;
            divergence[tipIndex] -= contribution;
        }
    }
    return divergence;
}

/*
 * Computes the geodesic distances φ using the heat method.
 *
 * Input: <delta>, a dense vector representing the heat sources, i.e., u0 = δ(x). Returns: A dense vector containing the
 * geodesic distances per vertex.
 */
Vector<double> HeatMethod::compute(const Vector<double>& delta) const {
    SparseMatrix<double> F = this->F;
    Eigen::SimplicialLLT<SparseMatrix<double>> llt(F);
    Vector<double> ut = llt.solve(delta);


    FaceData<Vector3> X = computeVectorField(ut);
    Vector<double> divX = computeDivergence(X);


    SparseMatrix<double> A = this->A;
    geometrycentral::PositiveDefiniteSolver<double> solver(A);
    Vector<double> phi = - solver.solve(divX);

    // Since φ is unique up to an additive constant, it should be shifted such that the smallest distance is zero
    this->subtractMinimumDistance(phi);

    return phi;
}