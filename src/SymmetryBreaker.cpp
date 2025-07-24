//
// Created by Jonas on 7/05/2024.
//

#include "SymmetryBreaker.h"
#include "glm/gtc/epsilon.hpp"
#include "meshcore/factories/OBBFactory.h"
#include "meshcore/core/WorldSpaceMesh.h"
#include "meshcore/acceleration/CachingBoundsTreeFactory.h"
#include "meshcore/acceleration/AABBVolumeHierarchy.h"

std::vector<std::vector<Plane>> SymmetryBreaker::deriveSymmetryBreakingPlanes(const std::shared_ptr<ModelSpaceMesh> &mesh, int maxRotationFold) {

    std::vector<std::vector<Plane>> result;

    // Temporary bypass for tetrahedron
    if(mesh->getName() == "tetrahedron.stl"){
        return result;
    }

    // Derive primary axis aligned with OBB
    OBB obb = OBBFactory::createOBB(mesh);
    auto primaryX = obb.getRotation().rotateVertex(glm::vec3(1, 0, 0));
    auto primaryY = obb.getRotation().rotateVertex(glm::vec3(0, 1, 0));
    auto primaryZ = obb.getRotation().rotateVertex(glm::vec3(0, 0, 1));

    Plane planeYZ(primaryX, obb.getCenter());
    Plane planeXZ(primaryY, obb.getCenter());
    Plane planeXY(primaryZ, obb.getCenter());

    if(isReflectionSymmetricAroundPlane(mesh, planeYZ)){
        result.push_back({planeYZ});
        std::cout << "Found reflection symmetry around OBB YZ plane" << std::endl;
    }
    if(isReflectionSymmetricAroundPlane(mesh, planeXZ)){
        result.push_back({planeXZ});
        std::cout << "Found reflection symmetry around OBB XZ plane" << std::endl;
    }
    if(isReflectionSymmetricAroundPlane(mesh, planeXY)){
        result.push_back({planeXY});
        std::cout << "Found reflection symmetry around OBB XY plane" << std::endl;
    }

    // Center around volume centroid
    auto centroid = mesh->getVolumeCentroid();
    auto worldSpaceMesh = std::make_shared<WorldSpaceMesh>(mesh);
    worldSpaceMesh->getModelTransformation().setPosition(-centroid);
    auto centeredModelSpaceMesh = worldSpaceMesh->getTransformedModelSpaceMesh();

    // Check for N-fold rotation symmetry (We start with the highest fold, as the resulting planes are most restrictive)
    for(auto n=maxRotationFold; n>=2; n--) {
        if (isNFoldRotationSymmetricAroundAxis(centeredModelSpaceMesh, primaryX, n)) {
            Plane plane(primaryY, centroid);
            auto angle = 2 * glm::pi<float>() / float(n);
            auto rotation = glm::angleAxis(angle, primaryX);
            auto rotatedNormal = rotation * primaryY;
            Plane plane2(-rotatedNormal, centroid);

            result.push_back({plane, plane2});
            std::cout << "Found " << n << "-fold rotation symmetry around OBB X axis" << std::endl;
            break;
        }
    }
    for(auto n=maxRotationFold; n>=2; n--) {
        if (isNFoldRotationSymmetricAroundAxis(centeredModelSpaceMesh, primaryY, n)) {
            Plane plane(primaryZ, centroid);
            auto angle = 2 * glm::pi<float>() / float(n);
            auto rotation = glm::angleAxis(angle, primaryY);
            auto rotatedNormal = rotation * primaryZ;
            Plane plane2(-rotatedNormal, centroid);

            result.push_back({plane, plane2});
            std::cout << "Found " << n << "-fold rotation symmetry around OBB Y axis" << std::endl;
            break;
        }
    }
    for(auto n=maxRotationFold; n>=2; n--) {
        if(isNFoldRotationSymmetricAroundAxis(centeredModelSpaceMesh, primaryZ, n)){
            Plane plane(primaryX, centroid);
            auto angle = 2 * glm::pi<float>() / float(n);
            auto rotation = glm::angleAxis(angle, -primaryZ);
            auto rotatedNormal = rotation * primaryX;
            Plane plane2(-rotatedNormal, centroid);

            result.push_back({plane, plane2});
            std::cout << "Found " << n << "-fold rotation symmetry around OBB Z axis" << std::endl;
            break;
        }
    }

    return result;
}

bool SymmetryBreaker::isReflectionSymmetricAroundPlane(const std::shared_ptr<ModelSpaceMesh> &mesh, const Plane &plane) {

    // In a more complete approach, we should check if each facet has a mirrored equivalent and not just the vertices
    // Check if all vertices map to a symmetric vertex
    for (const auto &vertex: mesh->getVertices()){

        // Mirror vertex around the plane
        auto signedDistance = plane.signedDistance(vertex);
        auto mirroredVertex = vertex - 2 * signedDistance * plane.getNormal();

        // Find matching vertex
        bool mirroredVertexFound = false;
        for (const auto &other: mesh->getVertices()){
            if(glm::epsilonEqual(glm::distance(mirroredVertex, other), 0.0f, 1e-3f * glm::length(mesh->getBounds().getHalf()))){
                mirroredVertexFound = true;
                break;
            }
        }
        if(!mirroredVertexFound){
            return false;
        }
    }
    return true;
}

bool SymmetryBreaker::isNFoldRotationSymmetricAroundAxis(const std::shared_ptr<ModelSpaceMesh> &mesh, const glm::vec3 &axis, int N) {

    assert(N > 1);

    auto bvh = CachingBoundsTreeFactory<AABBVolumeHierarchy>::getBoundsTree(mesh);

    // Check if all vertices map to a symmetric vertex
    // In a more complete approach, we should check if each facet has a rotated equivalent and not just the vertices
    for (const auto &vertex: mesh->getVertices()){

        // Rotate vertex around the axis
        auto angle = 2 * glm::pi<float>() / float(N);
        auto rotation = glm::angleAxis(angle, axis);
        auto rotatedVertex = vertex;

        // Each vertex should have an equivalent in each of the N segments
        for(auto n=1; n<N; n++){
            rotatedVertex = rotation * rotatedVertex;

            // Find other vertex
            bool rotatedVertexFound = false;
            for (const auto &other: mesh->getVertices()){

                if(glm::epsilonEqual(glm::distance(rotatedVertex, other), 0.0f, 1e-3f * glm::length(mesh->getBounds().getHalf()))){
                    rotatedVertexFound = true;
                    break;
                }
            }
            if(!rotatedVertexFound){
                return false;
            }
        }
    }
    return true;
}