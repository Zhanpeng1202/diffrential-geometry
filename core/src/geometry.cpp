// PLEASE READ:
//
// This file implements additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because
// we are "inside" the class, we no longer have to call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using geometry->cotan(he),
// geometry->barycentricDualArea(v), etc. where "geometry" is a pointer to a VertexPositionGeometry. This avoids having
// to declare a GeometryRoutines object in every project, and also mimics the way that geometry routines are normally
// called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.


#include "geometrycentral/surface/vertex_position_geometry.h"
#include <complex>

namespace geometrycentral {
namespace surface {

/*
 * Compute the Euler characteristic of the mesh.
 */
int VertexPositionGeometry::eulerCharacteristic() const {
    return (int)mesh.nVertices() - (int)mesh.nEdges() + (int)mesh.nFaces();
}

/*
 * Compute the mean length of all the edges in the mesh.
 *
 * Input:
 * Returns: The mean edge length.
 */
double VertexPositionGeometry::meanEdgeLength() const {

    double total = 0.0;
    for (Edge e : mesh.edges()) {
        total += edgeLength(e);
    }
    return total / mesh.nEdges();
}

/*
 * Compute the total surface area of the mesh.
 *
 * Input:
 * Returns: The surface area of the mesh.
 */
double VertexPositionGeometry::totalArea() const {

    double total = 0.0;
    for (Face f : mesh.faces()) {
        total += faceArea(f);
    }
    return total;
}

/*
 * Computes the cotangent of the angle opposite to a halfedge. (Do NOT use built-in function for this)
 *
 * Input: The halfedge whose cotan weight is to be computed.
 * Returns: The cotan of the angle opposite the given halfedge.
 */
double VertexPositionGeometry::cotan(Halfedge he) const {
    // cotant  =  dot product / cross product
    // (Hint: how are the dot and cross product of two vectors related to the cosine and sine of the angle between them?)

    //       v1
    //      /  \ half
    //     /    \ edge
    //    v3 --- v2
    //    cotant for he is cot (v3) 
    //   is (e2 dot e3) / (e2 cross e3)

    Vertex v1 = he.vertex();
    Vertex v2 = he.next().vertex();
    Vertex v3 = he.next().next().vertex();

    Vector3 e2 = inputVertexPositions[v1] - inputVertexPositions[v3];
    Vector3 e1 = inputVertexPositions[v2] - inputVertexPositions[v3];

    Vector3 cross;
    cross.x = e1.y * e2.z - e1.z * e2.y;
    cross.y = e1.z * e2.x - e1.x * e2.z;
    cross.z = e1.x * e2.y - e1.y * e2.x;

    double dot = e1.x * e2.x + e1.y * e2.y + e1.z * e2.z;
    // vector length
    double det = cross.norm();

    // there are cases where e1 parallel to e2
    // why this would happen????? Boundary???

    if (det == 0) {
        printf("--------------------cross product of cotan is 0\n");
        // print what is e1 and e2
        printf("--------------------e1: %f %f %f\n", e1.x, e1.y, e1.z);
        printf("--------------------e2: %f %f %f\n", e2.x, e2.y, e2.z);
    }

    double cot = dot / det;

    return cot;
}

/*
 * Computes the barycentric dual area of a vertex.
 *
 * Input: The vertex whose barycentric dual area is to be computed.
 * Returns: The barycentric dual area of the given vertex.
 */
double VertexPositionGeometry::barycentricDualArea(Vertex v) const {
    //
    //Note that the barycentric dual area 
    // associated with a vertex 
    // is equal to one-third the area of all triangles 
    // touching 
    // , i.e.,
    // 
    // A_i = 1/3 * sum(area(f_j))
    //      for all f_j in F(v)
    // 
    // where F(v) is the set of all faces touching vertex v.

    double area = 0.0;
    for (Face f : v.adjacentFaces()) {
        area += faceArea(f);
    }
    return area / 3.0;
}

/*
 * Computes the angle (in radians) at a given corner. (Do NOT use built-in function for this)
 *
 *
 * Input: The corner at which the angle needs to be computed.
 * Returns: The angle clamped between 0 and π.
 */
double VertexPositionGeometry::angle(Corner c) const {
    // // Get the vertex at this corner
    // Vertex v = c.vertex();
    
    // // Get the next and previous vertices in the face
    // Halfedge he = c.halfedge();
    // Vertex vNext = he.next().vertex();
    // Vertex vPrev = he.next().next().vertex();
    
    // // Compute vectors from the corner vertex to its adjacent vertices
    // Vector3 vec1 = inputVertexPositions[vNext] - inputVertexPositions[v];
    // Vector3 vec2 = inputVertexPositions[vPrev] - inputVertexPositions[v];
    
    // // Normalize the vectors
    // vec1 = vec1.normalize();
    // vec2 = vec2.normalize();
    
    // // Compute the dot product
    // double dotProduct = dot(vec1, vec2);
    
    // // Compute the angle using the arc cosine of the dot product
    // double angleValue = std::acos(dotProduct);
    
    // return angleValue;
    // WARNING: Logic duplicated between cached and immediate version

    Halfedge he = c.halfedge();
    Vector3 pA = inputVertexPositions[he.vertex()];
    he = he.next();
    Vector3 pB = inputVertexPositions[he.vertex()];
    he = he.next();
    Vector3 pC = inputVertexPositions[he.vertex()];

    GC_SAFETY_ASSERT(he.next() == c.halfedge(), "faces mush be triangular");

    double q = dot(unit(pB - pA), unit(pC - pA));
    q = clamp(q, -1.0, 1.0);
    double angle = std::acos(q);
    return angle;
}

/*
 * Computes the signed angle (in radians) between two adjacent faces. (Do NOT use built-in function for this)
 *
 * Input: The halfedge (shared by the two adjacent faces) on which the dihedral angle is computed.
 * Returns: The dihedral angle.
 */
double VertexPositionGeometry::dihedralAngle(Halfedge he) const {
    if (he.edge().isBoundary() || !he.edge().isManifold() ) {
       return 0.;
    }

    Vector3 N1 = faceNormal(he.face());
    Vector3 N2 = faceNormal(he.sibling().face());
    Vector3 pTail = inputVertexPositions[he.vertex()];
    Vector3 pTip = inputVertexPositions[he.next().vertex()];
    Vector3 edgeDir = unit(pTip - pTail);

    double dihedralAngle = atan2(dot(edgeDir, cross(N1, N2)), dot(N1, N2));
    return dihedralAngle;
}

/*
 * Computes the normal at a vertex using the "equally weighted" method.
 * i.e. the normal is the average of the normals of the faces that the vertex is adjacent to
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "equally weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalEquallyWeighted(Vertex v) const {

    Vector3 vertexNormal = Vector3::zero();
    for (Face f : v.adjacentFaces()) {
        vertexNormal += faceNormal(f);
    }

    return unit(vertexNormal);
}

/*
 * Computes the normal at a vertex using the "tip angle weights" method.
 # i.e. N = Sum of (N_i * Phi_i)
 * Noted that Phi_i is the each angle of the corner at all faces that the vertex is adjacent to

 * Input: The vertex on which the normal is to be computed.
 * Returns: The "tip angle weights" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAngleWeighted(Vertex v) const {
    Vector3 vertexNormal = Vector3::zero();
    for (Corner c : v.adjacentCorners()) {
        // vertexNormal += faceNormal(c.halfedge().face()) * angle(c);
        if (c.vertex() == v) {
            vertexNormal += faceNormal(c.halfedge().face()) * angle(c);
        }
    }
    return unit(vertexNormal);
}

/*
 * Computes the normal at a vertex using the "inscribed sphere" method.
 * 
 * i.e. N = Sum of Cross(e_ij, e_ik)
 * Cross  = e_ij x e_ik
 * i is the vertex on which the normal is to be computed
 * j and k are the two vertices that are adjacent to i
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "inscribed sphere" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalSphereInscribed(Vertex v) const {

    Vector3 vertexNormal = Vector3::zero();
    for (Corner c : v.adjacentCorners()) {
        Vector3 cross_product = cross(halfedgeVector(c.halfedge()), halfedgeVector(c.halfedge().next()));

        float eij = edgeLength(c.halfedge().edge());
        float eik = edgeLength(c.halfedge().next().edge());

        vertexNormal += cross_product / (eij * eik);
    }
    return unit(vertexNormal);
}

/*
 * Computes the normal at a vertex using the "face area weights" method.
 * i.e. N = Sum of (N_i * A_i)
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "face area weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAreaWeighted(Vertex v) const {

    Vector3 vertexNormal = Vector3::zero();
    for (Face f : v.adjacentFaces()) {
        vertexNormal += faceNormal(f) * faceArea(f);
    }
    return unit(vertexNormal);
}

/*
 * Computes the normal at a vertex using the "Gauss curvature" method.
 * KN_i = 1/2 * Sum of (theta_{ij} * unit(E_{ij}))
 * theta_{ij} is the dihedral angle for edge ij
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "Gauss curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalGaussianCurvature(Vertex v) const {

    Vector3 vertexNormal = Vector3::zero();
    for (Edge e : v.adjacentEdges()) {
        vertexNormal += dihedralAngle(e.halfedge()) * unit(halfedgeVector(e.halfedge()));
    }
    return unit(vertexNormal);
}

/*
 * Computes the normal at a vertex using the "mean curvature" method (equivalent to the "area gradient" method).
 *
 * Input: The vertex on which the normal is to be computed.
 * i.e. N = 1/2 Sum of (cot(alpha) + cot(beta)) * e_ij
 * alpha and beta are the angles at the vertex opposite to edge e_ij
 * e_ij is the edge between vertex i and j
 * Returns: The "mean curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalMeanCurvature(Vertex v) const {

    Vector3 vertexNormal = Vector3::zero();
    for (Edge e : v.adjacentEdges()) {
        vertexNormal += (cotan(e.halfedge()) + cotan(e.halfedge().twin())) * halfedgeVector(e.halfedge());
    }
    return unit(vertexNormal);
}

/*
 * Computes the angle defect at a vertex.
 *
 * Input: The vertex whose angle defect is to be computed.
 * i.e. angle defect = 2 * pi - sum of angles of all corners at the vertex
 * Returns: The angle defect of the given vertex.
 */
double VertexPositionGeometry::angleDefect(Vertex v) const {

    double angleDefect = 0.0;
    for (Corner c : v.adjacentCorners()) {
        angleDefect += angle(c);
    }
    return 2 * PI - angleDefect;
}

/*
 * Computes the total angle defect of the mesh.
 *
 * Input:
 * Returns: The total angle defect
 * from the gauss-bonnet theorem, it is equal to 2 * pi * (1 - \chi(M))
 */
double VertexPositionGeometry::totalAngleDefect() const {
    int chi = eulerCharacteristic();
    return 2 * PI * (1 - chi);
}

/*
 * Computes the (integrated) scalar mean curvature at a vertex.
 *
 * Input: The vertex whose mean curvature is to be computed.
 * Returns: The mean curvature at the given vertex.
 */
double VertexPositionGeometry::scalarMeanCurvature(Vertex v) const {
   double meanCurvature = 0.;
   for (Halfedge he : v.outgoingHalfedges()) {
      double len = edgeLength(he.edge());
      double alpha = edgeDihedralAngle(he.edge());
      meanCurvature += alpha * len / 2.;
   }
   return meanCurvature/2.;
}

/*
 * Computes the circumcentric dual area of a vertex.
 *
 * Input: The vertex whose circumcentric dual area is to be computed.
 * A_i = 1/8 Sum of (e_{ij}^2 * cot(alpha_{ij}) + e_{ik}^2 * cot(beta_{ik}))
 * we sum over all face, j and k are the two vertices that are adjacent to i
 * Returns: The circumcentric dual area of the given vertex.
 */
double VertexPositionGeometry::circumcentricDualArea(Vertex v) const {
    double area = 0.0;
    for (Face f : v.adjacentFaces()) {
        area += faceArea(f);
    }
    return area / 8.0;
}

/*
 * Computes the (pointwise) minimum and maximum principal curvature values at a vertex.
 *
 * Input: The vertex on which the principal curvatures need to be computed.
 * Returns: A std::pair containing the minimum and maximum principal curvature values at a vertex.
 */
std::pair<double, double> VertexPositionGeometry::principalCurvatures(Vertex v) const {

    double A = vertexDualArea(v);
    double H = vertexMeanCurvature(v) / A;
    double K = vertexGaussianCurvature(v) / A;

    double k1 = H - std::sqrt(H*H - K);
    double k2 = H + std::sqrt(H*H - K);

    return std::make_pair(k1, k2);
}


/*
 * Builds the sparse POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace matrix,
 * multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse positive definite Laplace matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::laplaceMatrix() const {

    // TODO
    return identityMatrix<double>(1); // placeholder
}

/*
 * Builds the sparse diagonal mass matrix containing the barycentric dual area of each vertex.
 *
 * Input:
 * Returns: Sparse mass matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::massMatrix() const {

    // TODO
    return identityMatrix<double>(1); // placeholder
}

/*
 * Builds the sparse complex POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace
 * matrix, multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse complex positive definite Laplace matrix for the mesh.
 */
SparseMatrix<std::complex<double>> VertexPositionGeometry::complexLaplaceMatrix() const {

    // TODO
    return identityMatrix<std::complex<double>>(1); // placeholder
}

/*
 * Compute the center of mass of a mesh.
 */
Vector3 VertexPositionGeometry::centerOfMass() const {

    // Compute center of mass.
    Vector3 center = {0.0, 0.0, 0.0};
    for (Vertex v : mesh.vertices()) {
        center += inputVertexPositions[v];
    }
    center /= mesh.nVertices();

    return center;
}

/*
 * Centers a mesh about the origin.
 * Also rescales the mesh to unit radius if <rescale> == true.
 */
void VertexPositionGeometry::normalize(const Vector3& origin, bool rescale) {

    // Compute center of mass.
    Vector3 center = centerOfMass();

    // Translate to origin [of original mesh].
    double radius = 0;
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] -= center;
        radius = std::max(radius, inputVertexPositions[v].norm());
    }

    // Rescale.
    if (rescale) {
        for (Vertex v : mesh.vertices()) {
            inputVertexPositions[v] /= radius;
        }
    }

    // Translate to origin [of original mesh].
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] += origin;
    }
}

} // namespace surface
} // namespace geometrycentral