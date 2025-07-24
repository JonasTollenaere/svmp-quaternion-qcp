//
// Created by Jonas on 7/05/2024.
//

#include <iostream>
#include "UpperBoundProvider.h"
#include "meshcore/factories/SphereFactory.h"


float UpperBoundProvider::getUpperScaleBoundConvexHullItem(const std::shared_ptr<ModelSpaceMesh> &item, const std::shared_ptr<ModelSpaceMesh> &container) {

    float upperBoundBoundingBox = glm::length(container->getBounds().getHalf())/glm::length(item->getBounds().getHalf());
    std::cout << "Upper bound bounding box: " << upperBoundBoundingBox << std::endl;
    float upperScaleBoundVolume = std::pow(container->getVolume()/item->getConvexHull()->getVolume(),1.0f/3.0f); // Convex hull used here! Not valid when we consider non-convex items
    std::cout << "Upper bound volume: " << upperScaleBoundVolume << std::endl;
    float upperBoundBoundingSphere = SphereFactory::createMinimumBoundingSphere(container->getVertices()).getRadius()/SphereFactory::createMinimumBoundingSphere(item->getVertices()).getRadius();
    std::cout << "Upper bound bounding sphere: " << upperBoundBoundingSphere << std::endl;

    float upperScaleBound = std::min(upperBoundBoundingSphere,std::min(upperBoundBoundingBox, upperScaleBoundVolume));

    return upperScaleBound;
}

float UpperBoundProvider::getUpperScaleBound(const std::shared_ptr<ModelSpaceMesh> &item, const std::shared_ptr<ModelSpaceMesh> &container) {

    float upperBoundBoundingBox = glm::length(container->getBounds().getHalf())/glm::length(item->getBounds().getHalf());
    std::cout << "Upper bound bounding box: " << upperBoundBoundingBox << std::endl;
    float upperScaleBoundVolume = std::pow(container->getVolume()/item->getVolume(),1.0f/3.0f);
    std::cout << "Upper bound volume: " << upperScaleBoundVolume << std::endl;
    float upperBoundBoundingSphere = SphereFactory::createMinimumBoundingSphere(container->getVertices()).getRadius()/SphereFactory::createMinimumBoundingSphere(item->getVertices()).getRadius();
    std::cout << "Upper bound bounding sphere: " << upperBoundBoundingSphere << std::endl;

    float upperScaleBound = std::min(upperBoundBoundingSphere,std::min(upperBoundBoundingBox, upperScaleBoundVolume));

    return upperScaleBound;
}