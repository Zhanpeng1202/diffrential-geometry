// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because we are
// "inside" the class, we no longer have to call
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
// Functions in this file can be called from other projects simply by using geometry->buildHodgeStar0Form(), etc. where
// "geometry" is a pointer to a VertexPositionGeometry. This avoids having to declare a GeometryRoutines object in every
// project, and also mimics the way that geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {


/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {

    SparseMatrix<double> res = identityMatrix<double>(mesh.nVertices());

    for (int i = 0; i < mesh.nVertices(); i++) {
        Vertex v = mesh.vertex(i);
        double area = this->barycentricDualArea(v);
        res.coeffRef(i, i) = area;
    }

    return res;
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {

    SparseMatrix<double> res = identityMatrix<double>(mesh.nEdges());

    for (int i = 0; i < mesh.nEdges(); i++) {
        Edge e = mesh.edge(i);
        double cotana = this->cotan(e.halfedge());

        double cotanb = this->cotan(e.halfedge().twin());
        res.coeffRef(i, i) = (cotana + cotanb) / 2;
    }

    return res;
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {

    SparseMatrix<double> res = identityMatrix<double>(mesh.nFaces());

    for (int i = 0; i < mesh.nFaces(); i++) {
        Face f = mesh.face(i);
        double area = this->faceArea(f);
        res.coeffRef(i, i) = 1.0 / area;
    }

    return res;
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {

    SparseMatrix<double> res = SparseMatrix<double>(mesh.nEdges(), mesh.nVertices());

    for (int i = 0; i < mesh.nEdges(); i++) {
        Edge e = mesh.edge(i);
        res.coeffRef(i, e.firstVertex().getIndex()) = -1;
        res.coeffRef(i, e.secondVertex().getIndex()) = 1;
    }

    return res;
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {

    SparseMatrix<double> res = SparseMatrix<double>(mesh.nFaces(), mesh.nEdges());

    for (int i = 0; i < mesh.nFaces(); i++) {
        Face f = mesh.face(i);
        for (Halfedge he : f.adjacentHalfedges()) {
            res.coeffRef(i, he.edge().getIndex()) = 
            (he.tailVertex().getIndex() == he.edge().firstVertex().getIndex()) ? 1 : -1;
        }
    }

    return res;
}

} // namespace surface
} // namespace geometrycentral