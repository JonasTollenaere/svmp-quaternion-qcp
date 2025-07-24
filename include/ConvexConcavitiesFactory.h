//
// Created by Jonas on 6/05/2024.
//

#ifndef EXTENDEDMESHCORE_CONVEXCONCAVITIESFACTORY_H
#define EXTENDEDMESHCORE_CONVEXCONCAVITIESFACTORY_H

#include "meshcore/core/ModelSpaceMesh.h"
#include "meshcore/core/Plane.h"
#include <unordered_map>

class ConvexConcavitiesFactory {

    static std::unordered_map<std::shared_ptr<ModelSpaceMesh>, std::vector<std::shared_ptr<ModelSpaceMesh>>> decompositionCache;
    static std::unordered_map<std::shared_ptr<ModelSpaceMesh>, std::vector<std::shared_ptr<ModelSpaceMesh>>> concavitiesCache;
    static std::unordered_map<std::shared_ptr<ModelSpaceMesh>, std::vector<IndexTriangle>> concaveTrianglesCache;

public:
    static std::vector<std::shared_ptr<ModelSpaceMesh>> getConvexDecomposition(const std::shared_ptr<ModelSpaceMesh>& inputMesh);
    static std::vector<std::shared_ptr<ModelSpaceMesh>> getConvexConcavities(const std::shared_ptr<ModelSpaceMesh>& inputMesh);
    static std::vector<IndexTriangle> getConcaveTriangles(const std::shared_ptr<ModelSpaceMesh> &inputMesh);

    static std::vector<std::vector<Plane>> getGroupedConcavityPlanes(std::shared_ptr<ModelSpaceMesh> inputMesh);
};


#endif //EXTENDEDMESHCORE_CONVEXCONCAVITIESFACTORY_H
