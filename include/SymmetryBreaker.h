//
// Created by Jonas on 7/05/2024.
//

#ifndef EXTENDEDMESHCORE_SYMMETRYBREAKER_H
#define EXTENDEDMESHCORE_SYMMETRYBREAKER_H

#include "meshcore/core/ModelSpaceMesh.h"
#include "meshcore/core/Plane.h"

class SymmetryBreaker {
public:
    static std::vector<std::vector<Plane>>  deriveSymmetryBreakingPlanes(const std::shared_ptr<ModelSpaceMesh> &mesh, int maxRotationFold=8);
private:
    static bool isReflectionSymmetricAroundPlane(const std::shared_ptr<ModelSpaceMesh> &mesh, const Plane &plane);
    static bool isNFoldRotationSymmetricAroundAxis(const std::shared_ptr<ModelSpaceMesh> &mesh, const glm::vec3 &axis, int N);
};


#endif //EXTENDEDMESHCORE_SYMMETRYBREAKER_H
