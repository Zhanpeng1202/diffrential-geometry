// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v]; 
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
        // geometry->faceIndices[f] = idx;
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    // TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.

    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();    

    SparseMatrix<size_t> result(mesh->nEdges(), mesh->nVertices());
    
    for (Edge e : mesh->edges()) {
        for (Vertex v : e.adjacentVertices()) {
            result.insert(geometry->edgeIndices[e], geometry->vertexIndices[v]) = 1;
        }
    }

    return result;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();    

    SparseMatrix<size_t> result(mesh->nFaces(), mesh->nEdges());
    
    for (Face f : mesh->faces()) {
        for (Edge e : f.adjacentEdges()) {
            result.insert(geometry->faceIndices[f], geometry->edgeIndices[e]) = 1;
        }
    }
    return result;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 * 
 * Implement the methods buildVertexVector, buildEdgeVector, buildFaceVector,
 * which each take a subset S of simplices as input, and construct a column vector encoding the
 * vertices, edges, or faces (respectively) in that subset. For instance, in buildVertexVector you
 * should build a column vector with |V| entries that has a “1” for each vertex in the subset, and “0”
 * for all other vertices.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();    

    Vector<size_t> result(mesh->nVertices());

    for (auto v : subset.vertices) {
        result[geometry->vertexIndices[v]] = 1;
    }
    return result;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();    

    Vector<size_t> result(mesh->nEdges());

    for (auto e : subset.edges) {
        result[geometry->edgeIndices[e]] = 1;
    }
    return result;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();    

    Vector<size_t> result(mesh->nFaces());

    for (auto f : subset.faces) {
        result[geometry->faceIndices[f]] = 1;
    }
    return result;
}


/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 * 
 * Implement the method star, which takes as input a subset S of simplices, and
 * computes the simplicial star St(S) of this subset. Hint: What happens if you apply the two unsigned
 * adjacency matrices in sequence? How do you get all the simplices in the star?
 * 
 * 
 * Star: union of simplices containing a given subset of simplices
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    // MeshSubset class:
    //   public:
    //     std::set<size_t> vertices;
    //     std::set<size_t> edges;
    //     std::set<size_t> faces;

    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();    

    // MeshSubset result(subset.vertices, subset.edges, subset.faces);
    MeshSubset result = subset.deepCopy();
    // result.vertices.clear();
    // result.edges.clear();
    // result.faces.clear();

    

    Vector<size_t> v = buildVertexVector(subset);
    Vector<size_t> e = buildEdgeVector(subset);
    Vector<size_t> f = buildFaceVector(subset);

    // A0, A1 are already computed
    // A0 is a vertex-edge adjacency matrix
    // A1 is a face-edge adjacency matrix

    // Vector<size_t> v_star = A;
    // Vector<size_t> e_star = A1 * e;
    // Vector<size_t> f_star = A0 * f + A1 * e;

    Vector<size_t> v_star = v;
    Vector<size_t> e_star = A0 * v;
    Vector<size_t> f_star = A1 * e_star +A1 * e;


    for (size_t i = 0; i < mesh->nVertices(); i++) {
        if (v_star[i] >0) {
            result.vertices.insert(i);
        }
    }

    for (size_t i = 0; i < mesh->nEdges(); i++) {
        if (e_star[i] > 0) {
            result.edges.insert(i);
        }
    }

    for (size_t i = 0; i < mesh->nFaces(); i++) {
        if (f_star[i] > 0) {
            result.faces.insert(i);
        }
    }

    return result;
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    // Closure: smallest simplicial complex containing a given set of simplices
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();   

    MeshSubset result = subset.deepCopy();

    Vector<size_t> v = buildVertexVector(subset);
    Vector<size_t> e = buildEdgeVector(subset);
    Vector<size_t> f = buildFaceVector(subset);

    Vector<size_t> f_cl = f;
    Vector<size_t> e_cl = A1.transpose()* f_cl;
    Vector<size_t> v_cl = A0.transpose() * e_cl + A0.transpose() * e;
    // Vector<size_t> v_cl = v;

    for (size_t i = 0; i < mesh->nVertices(); i++) {
        if (v_cl[i] > 0) {
            result.vertices.insert(i);
        }
    }

    for (size_t i = 0; i < mesh->nEdges(); i++) {
        if (e_cl[i] > 0) {
            result.edges.insert(i);
        }
    }

    for (size_t i = 0; i < mesh->nFaces(); i++) {
        if (f_cl[i] > 0) {
            result.faces.insert(i);
        }
    }



    return result; // placeholder
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    // Link: closure of the star minus the star of the closure
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();   

    // MeshSubset result = subset.deepCopy();
    MeshSubset closure_subset = star(closure(subset));
    MeshSubset star_subset = closure(star(subset));

    // Here is a naming mistake: star_subset is the star of the closure, not the closure of the star
    // Leaving it as is 
    MeshSubset result = star_subset.deepCopy();

    result.deleteVertices(closure_subset.vertices);
    result.deleteEdges(closure_subset.edges);
    result.deleteFaces(closure_subset.faces);

    return result; // placeholder
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 * 
 * Implement the methods isComplex and isPureComplex, which check whether
 * a given subset S is a simplicial complex, and a pure simpicial complex, resp. The latter method
 * should return the degree of the complex if it’s pure, and -1 otherwise. Hint: use the closure method
 * for the first part, plus the adjacency matrices for the second part.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    MeshSubset closure_subset = closure(subset);

    return closure_subset.equals(subset); // placeholder
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *  * Implement the methods isComplex and isPureComplex, which check whether
 * a given subset S is a simplicial complex, and a pure simpicial complex, resp. The latter method
 * should return the degree of the complex if it’s pure, and -1 otherwise. Hint: use the closure method
 * for the first part, plus the adjacency matrices for the second part.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();  
    if (isComplex(subset)) {
        return -1;
    }


    Vector<size_t> v = buildVertexVector(subset);
    Vector<size_t> e = buildEdgeVector(subset);
    Vector<size_t> f = buildFaceVector(subset);

    Vector<size_t> f_cl = f;
    Vector<size_t> e_cl = A1.transpose()* f_cl;
    Vector<size_t> v_cl = A0.transpose() * e_cl;

    std::set<size_t> v_cl_set;
    std::set<size_t> e_cl_set;
    std::set<size_t> f_cl_set;

    for (size_t i = 0; i < mesh->nVertices(); i++) {
        if (v_cl[i] > 0) {
            v_cl_set.insert(i);
        }
    }

    for (size_t i = 0; i < mesh->nEdges(); i++) {
        if (e_cl[i] > 0) {
            e_cl_set.insert(i);
        }
    }

    for (size_t i = 0; i < mesh->nFaces(); i++) {
        if (f_cl[i] > 0) {
            f_cl_set.insert(i);
        }
    }

    MeshSubset pure_complex_subest = MeshSubset(v_cl_set, e_cl_set, f_cl_set);

    if (pure_complex_subest.equals(subset)) {
        return 0;
    }

    return -1; // placeholder
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();  


    Vector<size_t> v = buildVertexVector(subset);
    Vector<size_t> e = buildEdgeVector(subset);
    Vector<size_t> f = buildFaceVector(subset);

    Vector<size_t> f_cl = f;
    Vector<size_t> e_cl = A1.transpose()* f_cl;
    // Vector<size_t> v_cl = A0.transpose() * e_cl;

    std::set<size_t> e_cl_set;
    for (size_t i = 0; i < mesh->nEdges(); i++) {
        if (e_cl[i] == 1) {
            e_cl_set.insert(i);
        }
    }
    std::set<size_t> v_cl_set;
    v_cl_set.clear();
    std::set<size_t> f_cl_set;
    f_cl_set.clear();
    MeshSubset result = MeshSubset(v_cl_set, e_cl_set, f_cl_set);
    MeshSubset closure_subset = closure(result);


    return closure_subset; // placeholder
}
