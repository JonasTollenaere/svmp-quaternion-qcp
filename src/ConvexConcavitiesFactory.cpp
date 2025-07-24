//
// Created by Jonas on 6/05/2024.
//

#include "ConvexConcavitiesFactory.h"
#include "meshcore/core/VertexTriangle.h"
#include "meshcore/acceleration/AABBVolumeHierarchy.h"
#include "meshcore/acceleration/CachingBoundsTreeFactory.h"

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/convex_decomposition_3.h>

//typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>                        SurfaceMesh;
typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor        halfedge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::edge_descriptor            edge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::face_descriptor            face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

typedef CGAL::Polyhedron_3<K> Polyhedron_3;
typedef CGAL::Nef_polyhedron_3<K, CGAL::SNC_indexed_items> Nef_polyhedron_3;
typedef Nef_polyhedron_3::Volume_const_iterator Volume_const_iterator;

std::shared_ptr<ModelSpaceMesh> getBoundingBoxAsModelMesh(const std::shared_ptr<ModelSpaceMesh>& inputMesh, float padding=0.0f);
SurfaceMesh getSurfaceMesh(const std::shared_ptr<ModelSpaceMesh>& inputMesh);
unsigned int findOrEmplaceVertex(std::vector<Vertex>& vertices, const Vertex& vertex);
std::shared_ptr<ModelSpaceMesh> tryConvexMerge(const std::shared_ptr<ModelSpaceMesh>& meshA, const std::shared_ptr<ModelSpaceMesh>& meshB);


std::unordered_map<std::shared_ptr<ModelSpaceMesh>, std::vector<std::shared_ptr<ModelSpaceMesh>>> ConvexConcavitiesFactory::decompositionCache{};
std::unordered_map<std::shared_ptr<ModelSpaceMesh>, std::vector<std::shared_ptr<ModelSpaceMesh>>> ConvexConcavitiesFactory::concavitiesCache{};
std::unordered_map<std::shared_ptr<ModelSpaceMesh>, std::vector<IndexTriangle>> ConvexConcavitiesFactory::concaveTrianglesCache{};

std::vector<std::shared_ptr<ModelSpaceMesh>> ConvexConcavitiesFactory::getConvexDecomposition(const std::shared_ptr<ModelSpaceMesh>& inputMesh){

    if(inputMesh->isConvex()) return {inputMesh};

    // 0. Check if we have already computed the convex concavities
    if(decompositionCache.find(inputMesh) != decompositionCache.end()){
        return decompositionCache[inputMesh];
    }

    // 1. Convert both to CGAL Meshes
    SurfaceMesh itemMesh = getSurfaceMesh(inputMesh);
    auto boundsMesh = getBoundingBoxAsModelMesh(inputMesh, 0.0f);

    SurfaceMesh convexHullMesh = getSurfaceMesh(boundsMesh);
    SurfaceMesh outputMesh;

    // Nef polyhedron of the item
    Polyhedron_3 inputPolyhedron;
    CGAL::copy_face_graph(itemMesh, inputPolyhedron);
    Nef_polyhedron_3 itemNefPolyhedron(inputPolyhedron);

    CGAL::convex_decomposition_3(itemNefPolyhedron);

    std::vector<std::shared_ptr<ModelSpaceMesh>> result;
    Volume_const_iterator ci = ++itemNefPolyhedron.volumes_begin();
    for( ; ci != itemNefPolyhedron.volumes_end(); ++ci) {
        if(ci->mark()) {
            Polyhedron_3 P;
            itemNefPolyhedron.convert_inner_shell_to_polyhedron(ci->shells_begin(), P);

            std::vector<Vertex> convexPartVertices;
            for(auto it = P.vertices_begin(); it != P.vertices_end(); ++it){
                auto vertex = it->point();
                glm::dvec3 itemVertex = {CGAL::to_double(vertex.x()), CGAL::to_double(vertex.y()), CGAL::to_double(vertex.z())};
                convexPartVertices.emplace_back(itemVertex);
            }

            auto emptyModelMesh = std::make_shared<ModelSpaceMesh>(convexPartVertices, std::vector<IndexTriangle>());
            auto convexModelMesh = emptyModelMesh->getConvexHull();
            result.emplace_back(convexModelMesh);
        }
    }
    std::cout << "Decomposition into " << result.size() << " convex parts " << std::endl;

    // Filter out parts that only have coplanar vertices
    {
        std::vector<std::shared_ptr<ModelSpaceMesh>> filteredResult;
        for (const auto &part: result){

            // Derive plane in which first triangle lies
            Vertex vertex0 = part->getVertices().at(part->getTriangles().at(0).vertexIndex0);
            Vertex vertex1 = part->getVertices().at(part->getTriangles().at(0).vertexIndex1);
            Vertex vertex2 = part->getVertices().at(part->getTriangles().at(0).vertexIndex2);
            VertexTriangle triangle({vertex0, vertex1, vertex2});
            auto normal = glm::normalize(triangle.normal);
            auto d = -glm::dot(normal, triangle.vertices[0]);

            // Check if all vertices are coplanar
            bool coplanar = true;
            for (const auto &vertex: part->getVertices()){
                auto delta = glm::normalize(vertex - triangle.vertices[0]);
                auto dot = glm::dot(normal, delta);
                if(glm::epsilonNotEqual(dot, 0.0f, 1e-3f) && glm::epsilonNotEqual(d, -glm::dot(normal, vertex), 1e-3f)){
                    coplanar = false;
                    break;
                }
            }

            // If not coplanar, add to filtered result
            if(!coplanar){
                filteredResult.emplace_back(part);
            }
        }
        result = filteredResult;
        std::cout << "Filtered into " << result.size() << " non-planar parts " << std::endl;
    }

    // Iteratively merge parts that form a convex volume together
    {
        bool mergeFound = true;
        while(mergeFound){
            mergeFound = false;
            for (int i = 0; i < result.size(); ++i){
                for (int j = i + 1; j < result.size(); ++j){

                    auto merged = tryConvexMerge(result[i], result[j]);
                    if(merged){
                        result[i] = merged;
                        result.erase(result.begin() + j);
                        mergeFound = true;
                        break;
                    }
                }
                if(mergeFound){
                    break;
                }
            }
        }
        std::cout << "Merged into " << result.size() << " convex parts " << std::endl;
    }

    // Cache the result
    decompositionCache[inputMesh] = result;
    return result;
}

std::vector<std::shared_ptr<ModelSpaceMesh>> ConvexConcavitiesFactory::getConvexConcavities(const std::shared_ptr<ModelSpaceMesh>& inputMesh){

    // 0. Check if we have already computed the convex concavities
    if(concavitiesCache.find(inputMesh) != concavitiesCache.end()){
        return concavitiesCache[inputMesh];
    }

    // 1. Convert both to CGAL Meshes
    SurfaceMesh itemMesh = getSurfaceMesh(inputMesh);
    auto boundsMesh = getBoundingBoxAsModelMesh(inputMesh, 0.0f);

//    FileParser::saveFile("boundsMesh.obj", boundsMesh);
    SurfaceMesh convexHullMesh = getSurfaceMesh(boundsMesh);
    SurfaceMesh outputMesh;

    // Nef polyhedron of the item
    Polyhedron_3 inputPolyhedron;
    CGAL::copy_face_graph(itemMesh, inputPolyhedron);
    Nef_polyhedron_3 itemNefPolyhedron(inputPolyhedron);

    // Nef polyhedron of bounding box
    Polyhedron_3 boundingBoxPolyhedron;
    auto boundingBox = getBoundingBoxAsModelMesh(inputMesh, 0.0f);
    CGAL::copy_face_graph(getSurfaceMesh(boundingBox), boundingBoxPolyhedron);
    Nef_polyhedron_3 boundingBoxNefPolyhedron(boundingBoxPolyhedron);

    // Nef polyhedron of convex hull
    Polyhedron_3 convexHullPolyhedron;
    CGAL::copy_face_graph(getSurfaceMesh(inputMesh->getConvexHull()), convexHullPolyhedron);
    Nef_polyhedron_3 convexHullNefPolyhedron(convexHullPolyhedron);

    // Get the difference
    auto diff = convexHullNefPolyhedron - itemNefPolyhedron;
//    auto diff = boundingBoxNefPolyhedron - itemNefPolyhedron;

    CGAL::convex_decomposition_3(diff);

    std::vector<std::shared_ptr<ModelSpaceMesh>> result;
    Volume_const_iterator ci = ++diff.volumes_begin();
    for( ; ci != diff.volumes_end(); ++ci) {
        if(ci->mark()) {
            Polyhedron_3 P;
            diff.convert_inner_shell_to_polyhedron(ci->shells_begin(), P);

            std::vector<Vertex> convexPartVertices;
            for(auto it = P.vertices_begin(); it != P.vertices_end(); ++it){
                auto vertex = it->point();
                glm::dvec3 itemVertex = {CGAL::to_double(vertex.x()), CGAL::to_double(vertex.y()), CGAL::to_double(vertex.z())};
                convexPartVertices.emplace_back(itemVertex);
            }

            auto emptyModelMesh = std::make_shared<ModelSpaceMesh>(convexPartVertices, std::vector<IndexTriangle>());
            auto convexModelMesh = emptyModelMesh->getConvexHull();
            result.emplace_back(convexModelMesh);
        }
    }
    std::cout << "Decomposition into " << result.size() << " convex parts " << std::endl;

    // Filter out parts that only have coplanar vertices
    {
        std::vector<std::shared_ptr<ModelSpaceMesh>> filteredResult;
        for (const auto &part: result){

            // Derive plane in which first triangle lies
            Vertex vertex0 = part->getVertices().at(part->getTriangles().at(0).vertexIndex0);
            Vertex vertex1 = part->getVertices().at(part->getTriangles().at(0).vertexIndex1);
            Vertex vertex2 = part->getVertices().at(part->getTriangles().at(0).vertexIndex2);
            VertexTriangle triangle({vertex0, vertex1, vertex2});
            auto normal = glm::normalize(triangle.normal);
            auto d = -glm::dot(normal, triangle.vertices[0]);

            // Check if all vertices are coplanar
            bool coplanar = true;
            for (const auto &vertex: part->getVertices()){
                auto delta = glm::normalize(vertex - triangle.vertices[0]);
                auto dot = glm::dot(normal, delta);
                if(glm::epsilonNotEqual(dot, 0.0f, 1e-3f) && glm::epsilonNotEqual(d, -glm::dot(normal, vertex), 1e-3f)){
                    coplanar = false;
                    break;
                }
            }

            // If not coplanar, add to filtered result
            if(!coplanar){
                filteredResult.emplace_back(part);
            }
        }
        result = filteredResult;
        std::cout << "Filtered into " << result.size() << " non-planar parts " << std::endl;
    }

    // Filter out parts that have no volume
    {
        std::vector<std::shared_ptr<ModelSpaceMesh>> filteredResult;
        for (const auto &part: result){
            if(part->getVolume() > 1e-6f){
                filteredResult.emplace_back(part);
            }
        }
        result = filteredResult;
        std::cout << "Filtered into " << result.size() << " parts with volume " << std::endl;
    }

    // Select parts that have at least one face in common with a concave triangle
//    {
//        std::vector<IndexTriangle> concaveTriangles = getConcaveTriangles(inputMesh);
//
//        std::vector<std::shared_ptr<ModelSpaceMesh>> filteredResult;
//        for (const auto &part: result){
//
//            bool partAdded = false;
//
//            for (const auto &concaveTriangle: concaveTriangles){
//
//                Vertex vertex0 = inputMesh->getVertices().at(concaveTriangle.vertexIndex0);
//                Vertex vertex1 = inputMesh->getVertices().at(concaveTriangle.vertexIndex1);
//                Vertex vertex2 = inputMesh->getVertices().at(concaveTriangle.vertexIndex2);
//
//                VertexTriangle triangle({vertex0, vertex1, vertex2});
//
//                auto normal = glm::normalize(triangle.normal);
//                auto d = -glm::dot(normal, triangle.vertices[0]);
//
//                for (const auto &partTriangle: part->getTriangles()){
//                    Vertex partVertex0 = part->getVertices().at(partTriangle.vertexIndex0);
//                    Vertex partVertex1 = part->getVertices().at(partTriangle.vertexIndex1);
//                    Vertex partVertex2 = part->getVertices().at(partTriangle.vertexIndex2);
//
//                    VertexTriangle partVertexTriangle({partVertex0, partVertex1, partVertex2});
//
//                    auto partNormal = -glm::normalize(partVertexTriangle.normal); // negated!
//                    auto partD = -glm::dot(partNormal, partVertexTriangle.vertices[0]);
//
//                    if(glm::all(glm::epsilonEqual(normal, partNormal, 1e-3f)) && glm::epsilonEqual(d, partD, 1e-3f)){
//                        filteredResult.emplace_back(part);
//                        partAdded = true;
//                        break;
//                    }
//                }
//
//                if(partAdded){
//                    break;
//                }
//            }
//
//            if(partAdded){
//                continue;
//            }
//        }
//        result = filteredResult;
//        std::cout << "Filtered into " << result.size() << " convex parts that contain concave item facets " << std::endl;
//    }

    // Iteratively merge parts that form a convex volume together
    {
        bool mergeFound = true;
        while(mergeFound){
            mergeFound = false;
            for (int i = 0; i < result.size(); ++i){
                for (int j = i + 1; j < result.size(); ++j){

                    auto merged = tryConvexMerge(result[i], result[j]);
                    if(merged){
                        result[i] = merged;
                        result.erase(result.begin() + j);
                        mergeFound = true;
                        break;
                    }
                }
                if(mergeFound){
                    break;
                }
            }
        }
        std::cout << "Merged into " << result.size() << " convex parts " << std::endl;
    }

    // can they be simplified? Only convex hulls of parts that touch the actual container? (distance to container or it's convex hull <= 0 + eps)
//    auto containerTree = CachingBoundsTreeFactory<AABBVolumeHierarchy>::getBoundsTree(inputMesh);
//    for (int i = 0; i < result.size(); ++i){
//
//        // Filter out vertices far from the container
//        std::vector<Vertex> filteredVertices;
//        for (const auto &vertex: result[i]->getVertices()){
//            auto d = containerTree->getShortestDistanceSquared(vertex);
//
//            if(d < 1e-6f) {
//                filteredVertices.emplace_back(vertex);
//            }
//        }
//
//        result[i] = ModelSpaceMesh(filteredVertices, {}).getConvexHull();
//    }

    // Rename the components
    for (int i = 0; i < result.size(); ++i){
        result[i]->setName("Concavity-" + std::to_string(i));
    }


    // Cache the result
    concavitiesCache[inputMesh] = result;
    return result;
}

std::vector<IndexTriangle> ConvexConcavitiesFactory::getConcaveTriangles(const std::shared_ptr<ModelSpaceMesh> &inputMesh) {

    // 0. Check if we have already computed the concave triangles
    if(concaveTrianglesCache.find(inputMesh) != concaveTrianglesCache.end()){
        return concaveTrianglesCache[inputMesh];
    }

    // 1. Find the concave triangles
    std::vector<IndexTriangle> concaveTriangles;
    for (const auto &indexTriangle: inputMesh->getTriangles()){

        Vertex vertex0 = inputMesh->getVertices().at(indexTriangle.vertexIndex0);
        Vertex vertex1 = inputMesh->getVertices().at(indexTriangle.vertexIndex1);
        Vertex vertex2 = inputMesh->getVertices().at(indexTriangle.vertexIndex2);

        VertexTriangle triangle({vertex0, vertex1, vertex2});

        auto normal = glm::normalize(triangle.normal);

        for (const auto &vertex : inputMesh->getVertices()){

            // The dot product between the triangle's normal and the vector from the triangle to the vertex should be negative
            auto delta = glm::normalize(vertex - vertex0);
            auto dot = glm::dot(normal, delta);
            if(dot > 1e-4){
                concaveTriangles.emplace_back(indexTriangle);
                break;
            }
        }
    }

    // Cache the result
//    concaveTrianglesCache[inputMesh] = concaveTriangles;
    return concaveTriangles;
}

std::vector<std::vector<Plane>> ConvexConcavitiesFactory::getGroupedConcavityPlanes(std::shared_ptr<ModelSpaceMesh> inputMesh) {
    std::vector<std::vector<Plane>> result;

    auto convexConcavities = getConvexConcavities(inputMesh);
    auto convexHull = inputMesh->getConvexHull();

    for (const auto &concavity: convexConcavities){
        std::vector<Plane> planes;
        for (const auto &triangle: concavity->getTriangles()){
            Vertex vertex0 = concavity->getVertices().at(triangle.vertexIndex0);
            Vertex vertex1 = concavity->getVertices().at(triangle.vertexIndex1);
            Vertex vertex2 = concavity->getVertices().at(triangle.vertexIndex2);
            auto partNormal = glm::normalize(VertexTriangle({vertex0, vertex1, vertex2}).normal);
            auto partD = -glm::dot(partNormal, vertex0);

            // Shouldn't be added if it's coplanar with a facet on the convex hull
            bool coplanar = false;
            for (const auto &convexTriangle: convexHull->getTriangles()){
                Vertex convexVertex0 = convexHull->getVertices().at(convexTriangle.vertexIndex0);
                Vertex convexVertex1 = convexHull->getVertices().at(convexTriangle.vertexIndex1);
                Vertex convexVertex2 = convexHull->getVertices().at(convexTriangle.vertexIndex2);

                VertexTriangle convexVertexTriangle({convexVertex0, convexVertex1, convexVertex2});
                auto normal = glm::normalize(convexVertexTriangle.normal);
                auto d = -glm::dot(normal, convexVertexTriangle.vertices[0]);


                if(glm::all(glm::epsilonEqual(normal, partNormal, 5e-3f)) && glm::epsilonEqual(d, partD, 5e-3f)){
                    coplanar = true;
                    break;
                }
            }
            if(coplanar){
                continue;
            }

            // But it should align with at least one facet of the original mesh (Inverted normal!)
            bool aligned = false;
            for (const auto &originalTriangle: inputMesh->getTriangles()){
                Vertex originalVertex0 = inputMesh->getVertices().at(originalTriangle.vertexIndex0);
                Vertex originalVertex1 = inputMesh->getVertices().at(originalTriangle.vertexIndex1);
                Vertex originalVertex2 = inputMesh->getVertices().at(originalTriangle.vertexIndex2);

                VertexTriangle originalVertexTriangle({originalVertex0, originalVertex1, originalVertex2});
                auto normal = glm::normalize(originalVertexTriangle.normal);

                normal = - normal; // The facet will be oriented in the opposite direction

                auto d = -glm::dot(normal, originalVertexTriangle.vertices[0]);

                if(glm::all(glm::epsilonEqual(normal, partNormal, 5e-3f)) && glm::epsilonEqual(d, partD, 5e-3f)){
                    aligned = true;
                    break;
                }
            }
            if(!aligned){
                continue;
            }

            // Shouldn't be added if already added
            bool alreadyAdded = false;
            for (const auto &plane: planes){
                if(glm::all(glm::epsilonEqual(plane.getNormal(), partNormal, 5e-3f)) && glm::epsilonEqual(plane.getD(), partD, 5e-3f)){
                    alreadyAdded = true;
                    break;
                }
            }
            if(alreadyAdded){
                continue;
            }

            // Add the plane
            planes.emplace_back(partNormal, partD);
        }
        result.emplace_back(planes);
    }
    return result;
}

SurfaceMesh getSurfaceMesh(const std::shared_ptr<ModelSpaceMesh>& inputMesh){

    // Convert this to a CGAL Mesh
    std::vector<Vertex> vertices;
    std::vector<std::vector<std::size_t>> faces;
    for (const auto &triangle : inputMesh->getTriangles()){
        std::vector<std::size_t> triangle_indices;
        triangle_indices.emplace_back(findOrEmplaceVertex(vertices, inputMesh->getVertices()[triangle.vertexIndex0]));
        triangle_indices.emplace_back(findOrEmplaceVertex(vertices, inputMesh->getVertices()[triangle.vertexIndex1]));
        triangle_indices.emplace_back(findOrEmplaceVertex(vertices, inputMesh->getVertices()[triangle.vertexIndex2]));
        faces.emplace_back(triangle_indices);
    }

    std::vector<K::Point_3> points;
    points.reserve(vertices.size());
    for (const auto &vertex: vertices){
        points.emplace_back(vertex.x, vertex.y, vertex.z);
    }

    SurfaceMesh mesh;
    CGAL::Polygon_mesh_processing::orient_polygon_soup(points, faces);
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, faces, mesh);

    return mesh;
}

unsigned int findOrEmplaceVertex(std::vector<Vertex>& vertices, const Vertex& vertex){
    for(unsigned int i = 0; i < vertices.size(); i++){
        if(vertices[i] == vertex){
            return i;
        }
    }
    vertices.push_back(vertex);
    return vertices.size() - 1;
}

std::shared_ptr<ModelSpaceMesh> getBoundingBoxAsModelMesh(const std::shared_ptr<ModelSpaceMesh>& inputMesh, float padding){
    auto min = inputMesh->getBounds().getMinimum() - glm::vec3(padding);
    auto max = inputMesh->getBounds().getMaximum() + glm::vec3(padding);

    std::vector<Vertex> vertices;
    vertices.emplace_back(min.x, min.y, min.z);
    vertices.emplace_back(max.x, min.y, min.z);
    vertices.emplace_back(max.x, max.y, min.z);
    vertices.emplace_back(min.x, max.y, min.z);
    vertices.emplace_back(min.x, min.y, max.z);
    vertices.emplace_back(max.x, min.y, max.z);
    vertices.emplace_back(max.x, max.y, max.z);
    vertices.emplace_back(min.x, max.y, max.z);

    std::vector<IndexTriangle> triangles;
    triangles.emplace_back(IndexTriangle{1, 0, 2});
    triangles.emplace_back(IndexTriangle{2, 0, 3});
    triangles.emplace_back(IndexTriangle{5, 1, 6});
    triangles.emplace_back(IndexTriangle{6, 1, 2});
    triangles.emplace_back(IndexTriangle{4, 5, 7});
    triangles.emplace_back(IndexTriangle{7, 5, 6});
    triangles.emplace_back(IndexTriangle{0, 4, 3});
    triangles.emplace_back(IndexTriangle{3, 4, 7});
    triangles.emplace_back(IndexTriangle{2, 3, 6});
    triangles.emplace_back(IndexTriangle{6, 3, 7});
    triangles.emplace_back(IndexTriangle{5, 4, 1});
    triangles.emplace_back(IndexTriangle{1, 4, 0});

    return std::make_shared<ModelSpaceMesh>(vertices, triangles);
}

std::shared_ptr<ModelSpaceMesh> tryConvexMerge(const std::shared_ptr<ModelSpaceMesh>& meshA, const std::shared_ptr<ModelSpaceMesh>& meshB){

    // Merge the vertices
    std::vector<Vertex> vertices = meshA->getVertices();
    for (const auto &vertex: meshB->getVertices()){
        vertices.emplace_back(vertex);
    }

    // Create the convex hull around the merged vertices
    auto mergedMesh = std::make_shared<ModelSpaceMesh>(vertices, std::vector<IndexTriangle>());
    auto convexHull = mergedMesh->getConvexHull();

    // Check their volumes
    float convexVolume = convexHull->getVolume();
    float inputVolume = meshA->getVolume() + meshB->getVolume();

    if(glm::epsilonEqual(convexVolume/inputVolume, 1.0f, 1e-3f)){
        return convexHull;
    }
    return nullptr;
}