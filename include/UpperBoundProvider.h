//
// Created by Jonas on 7/05/2024.
//

#ifndef EXTENDEDMESHCORE_UPPERBOUNDPROVIDER_H
#define EXTENDEDMESHCORE_UPPERBOUNDPROVIDER_H


#include "meshcore/core/ModelSpaceMesh.h"

class UpperBoundProvider {
public:
    static float getUpperScaleBoundConvexHullItem(const std::shared_ptr<ModelSpaceMesh> &item, const std::shared_ptr<ModelSpaceMesh> &container);
    static float getUpperScaleBound(const std::shared_ptr<ModelSpaceMesh> &item, const std::shared_ptr<ModelSpaceMesh> &container);
};


#endif //EXTENDEDMESHCORE_UPPERBOUNDPROVIDER_H
