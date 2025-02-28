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



    SparseMatrix<double> star0 = SparseMatrix<double> (mesh.nVertices(), mesh.nVertices());

    for (Vertex v : mesh.vertices()) {
        double area = barycentricDualArea(v);

        star0.insert(v.getIndex(), v.getIndex()) = area;
    }

    return star0;
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {




    SparseMatrix<double> star1 = SparseMatrix<double> (mesh.nEdges(), mesh.nEdges());

    for (Edge e : mesh.edges()) {
        double cot = 0.0;
        cot += cotan(e.halfedge());
        cot += cotan(e.halfedge().twin());
        star1.insert(e.getIndex(), e.getIndex()) = cot / 2.0;
    }

    return star1;
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {
    SparseMatrix<double> star2 = SparseMatrix<double> (mesh.nFaces(), mesh.nFaces());

    for (Face f : mesh.faces()) {
        double area = faceArea(f);
        star2.insert(f.getIndex(), f.getIndex()) = 1.0 / area;
    }
    return star2;
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {

    SparseMatrix<double> d0 = SparseMatrix<double>(mesh.nEdges(), mesh.nVertices());

    // can we do with halfedge?
    
    // print if edge is oriented
    // for (Edge e : mesh.edges()) {
    //     std::cout << "Edge " << e.getIndex() << " is oriented: " << e.isOriented() << std::endl;
    // }
    for (Edge e : mesh.edges()) {
        // if the edge is on the boundary, reverse the direction
        if (e.isBoundary()) {
            // res.insert(e.getIndex(), e.firstVertex().getIndex()) = -1;
            // res.insert(e.getIndex(), e.secondVertex().getIndex()) = 1;
        } else {
            d0.insert(e.getIndex(), e.firstVertex().getIndex()) = +1;
            d0.insert(e.getIndex(), e.secondVertex().getIndex()) = -1;
        }
    }


    // Halfedge Navigators
//   Halfedge twin() const;
//   Halfedge sibling() const;
//   Halfedge nextOutgoingNeighbor() const; // next halfedge which has the same tail vertex as this, form a cycle
//   Halfedge nextIncomingNeighbor() const; // next halfedge which has the same tip vertex as this, form a cycle
//   Halfedge next() const;
//   Corner corner() const;
//   Vertex vertex() const;
//   Vertex tipVertex() const;
//   Vertex tailVertex() const;
//   Edge edge() const;
//   Face face() const;
//   bool isDead() const;
    // for (Edge e : mesh.edges()) {
    // for (Halfedge he : mesh.halfedges()) {
    //     res.insert(he.edge().getIndex(), he.vertex().getIndex()) = -1;
    //     res.insert(he.edge().getIndex(), he.twin().vertex().getIndex()) = 1;
    // }

    return d0;

}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {

    SparseMatrix<double> d1 = SparseMatrix<double>(mesh.nFaces(), mesh.nEdges());

    // for (Face f : mesh.faces()) {
    //     for (Halfedge he : f.adjacentHalfedges()) {
    //         res.insert(f.getIndex(), he.edge().getIndex()) = (he == he.edge().halfedge()) ? 1 : -1;
    //     }
    // }

    for (Face f : mesh.faces()) {
        for (Halfedge he : f.adjacentHalfedges()) {
            if (he == he.edge().halfedge()) {
                d1.insert(f.getIndex(), he.edge().getIndex()) = 1;
            } else {
                d1.insert(f.getIndex(), he.edge().getIndex()) = -1;
            }

        }
    }
    return d1;
}

} // namespace surface
} // namespace geometrycentral