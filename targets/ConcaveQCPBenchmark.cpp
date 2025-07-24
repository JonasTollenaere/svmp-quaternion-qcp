//
// Created by Jonas on 17/05/2024.
//

#include "gurobi_c++.h"

#include <meshcore/utility/FileParser.h>
#include <meshcore/utility/io.h>
#include <meshcore/optimization/SingleVolumeMaximisationSolution.h>

#include <iostream>
#include <string>
#include <array>
#include <glm/gtx/component_wise.hpp>
#include <fstream>

#include "SymmetryBreaker.h"
#include "UpperBoundProvider.h"
#include "ConvexConcavitiesFactory.h"

static std::string dataFolder = MESHCORE_DATA_DIR + std::string("Tollenaere, J. et al/");

auto containerFiles = {
        "Items/arrow.obj",
        "Items/ring.obj",
        "Items/star.obj",
        "Items/star.stl",
        "Items/torus.stl",
        "Items/bobbin.stl",
        "Items/stone_1.obj",
        "Items/tree.stl",
};
auto itemFiles = {
        "Items/convexstone.stl",
        "Items/cube.obj",
        "Items/dodecahedron.stl",
        "Items/frustum.stl",
        "Items/gem.stl",
        "Items/helmet.stl",
        "Items/oloid.stl",
        "Items/tetrahedron.stl",
};

int main() {

    // Open a .csv writer append
    bool append = true;
    std::ofstream csvFile;
    csvFile.open("QuaternionQCPConcaveResults.csv", append ? std::ios_base::app : std::ios_base::out);

    if(!csvFile.is_open()){
        std::cerr << "Could not open file for writing." << std::endl;
        return 1;
    }

    // Write the header
    if(!append) csvFile << "Convex relaxation bound?,Limit separation planes?,Container,Item,Time(s),Scale,Gap,Explored,Unexplored,TranslationX,TranslationY,TranslationZ,RotationW,RotationX,RotationY,RotationZ" << std::endl;
    csvFile.flush();

    auto useConvexRelaxations = {true};
    auto useLimitSeparationPlanes = {false};

    for (const auto &containerFile: containerFiles) {
        for (const auto &itemFile: itemFiles) {

            for (const auto &useConvexRelaxation: useConvexRelaxations) {
                for(const auto& limitSeparationPlanes: useLimitSeparationPlanes) {

                    auto containerPath = dataFolder + containerFile;
                    auto itemPath = dataFolder + itemFile;

                    auto containerModelSpaceMesh = FileParser::loadMeshFile(containerPath);
                    auto itemModelMesh = FileParser::loadMeshFile(itemPath);

                    if(!containerModelSpaceMesh || !itemModelMesh){
                        std::cout << "Could not load mesh files" << std::endl;
                        continue;
                    }

                    auto containerMesh = std::make_shared<WorldSpaceMesh>(containerModelSpaceMesh);
                    auto itemMesh = std::make_shared<WorldSpaceMesh>(itemModelMesh);
                    auto solution = std::make_shared<SingleVolumeMaximisationSolution>(containerMesh, itemMesh);
                    solution->getItemWorldSpaceMesh()->getModelTransformation().setScale(0.0f);

                    std::cout << "Processing: " << itemFile << " in " << containerFile << ". Convex relaxation bound: " << useConvexRelaxation << ". Limit separation planes: " << limitSeparationPlanes << std::endl;

                    auto upperScaleBound = UpperBoundProvider::getUpperScaleBoundConvexHullItem(itemModelMesh, containerModelSpaceMesh);
                    double timeInFirstModel = 0.0;
                    if(useConvexRelaxation){

                        auto convexContainerHull = containerModelSpaceMesh->getConvexHull();
                        try {
                            // Create an environment
                            GRBEnv env = GRBEnv(true);
                            env.set("LogFile", "quat.log");
                            env.start();

                            // Create an empty model
                            GRBModel model = GRBModel(env);
                            model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
                            model.set(GRB_IntParam_Threads, 64);
                            model.set(GRB_IntParam_NonConvex, 2);

                            // Limit execution time to 3600s
                            model.set(GRB_DoubleParam_TimeLimit, 10 * 3600);

                            // Initialize decision variables
                            auto relaxedScaleBound = UpperBoundProvider::getUpperScaleBound(itemModelMesh, convexContainerHull);
                            GRBVar S = model.addVar(0, relaxedScaleBound, 1.0, GRB_CONTINUOUS, "S");

                            GRBVar Tx = model.addVar(-std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), 0, GRB_CONTINUOUS, "Tx");
                            GRBVar Ty = model.addVar(-std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), 0, GRB_CONTINUOUS, "Ty");
                            GRBVar Tz = model.addVar(-std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), 0, GRB_CONTINUOUS, "Tz");

                            GRBVar Qw = model.addVar(-1, 1, 0, GRB_CONTINUOUS, "Qw");
                            GRBVar Qx = model.addVar(-1, 1, 0, GRB_CONTINUOUS, "Qx");
                            GRBVar Qy = model.addVar(-1, 1, 0, GRB_CONTINUOUS, "Qy");
                            GRBVar Qz = model.addVar(-1, 1, 0, GRB_CONTINUOUS, "Qz");

                            // Constrain the quaternion to be normalized
                            model.addQConstr(Qw * Qw + Qx * Qx + Qy * Qy + Qz * Qz == 1, "|Q| = 1");

                            // Support variables
                            std::array<GRBVar, 9> R; // Rotation matrix
                            std::array<GRBVar, 9> SR; // Scaling + Rotation matrix
                            for (int i = 0; i < 9; i++) {
                                R[i] = model.addVar(-1, 1, 0, GRB_CONTINUOUS, "R[" + std::to_string(i) + "]");
                                SR[i] = model.addVar(-relaxedScaleBound, relaxedScaleBound, 0,  GRB_CONTINUOUS, "SR[" + std::to_string(i) + "]");
                                model.addQConstr(SR[i] == S * R[i], "SR[" + std::to_string(i) + "] = S * R[" + std::to_string(i) + "]");
                            }

                            model.addQConstr(R[0] == 1 - 2 * (Qy * Qy + Qz * Qz), "mat[0]  = 1 - 2 * ( yy + zz )");
                            model.addQConstr(R[1] == 2 * (Qx * Qy - Qz * Qw), "mat[1]  = 2 * ( xy - zw )");
                            model.addQConstr(R[2] == 2 * (Qx * Qz + Qy * Qw), "mat[2]  = 2 * ( xz + yw )");
                            model.addQConstr(R[3] == 2 * (Qx * Qy + Qz * Qw), "mat[3]  = 2 * ( xy + zw )");
                            model.addQConstr(R[4] == 1 - 2 * (Qx * Qx + Qz * Qz), "mat[4]  = 1 - 2 * ( xx + zz )");
                            model.addQConstr(R[5] == 2 * (Qy * Qz - Qx * Qw), "mat[5]  = 2 * ( yz - xw )");
                            model.addQConstr(R[6] == 2 * (Qx * Qz - Qy * Qw), "mat[6]  = 2 * ( xz - yw )");
                            model.addQConstr(R[7] == 2 * (Qy * Qz + Qx * Qw), "mat[7]  = 2 * ( yz + xw )");
                            model.addQConstr(R[8] == 1 - 2 * (Qx * Qx + Qy * Qy), "mat[8]  = 1 - 2 * ( xx + yy )");

                            // Symmetry breaking constraints
                            {
                                auto v = itemModelMesh->getVolumeCentroid();

                                for (const auto &symmetryBreakingPlaneSet: SymmetryBreaker::deriveSymmetryBreakingPlanes(containerModelSpaceMesh, 256)){
                                    for (const auto &symmetryBreakingPlane: symmetryBreakingPlaneSet){

                                        // Force position to be at positive side of the planes
                                        auto n = symmetryBreakingPlane.getNormal();
                                        auto d = symmetryBreakingPlane.getD();
                                        model.addConstr(n.x * (SR[0] * v.x + SR[1] * v.y + SR[2] * v.z + Tx) +
                                                        n.y * (SR[3] * v.x + SR[4] * v.y + SR[5] * v.z + Ty) +
                                                        n.z * (SR[6] * v.x + SR[7] * v.y + SR[8] * v.z + Tz) + d >= 0, "Container symmetry breaking constraint");
                                    }
                                }
                            }
                            {
                                for (const auto &symmetryPlaneSet: SymmetryBreaker::deriveSymmetryBreakingPlanes(itemModelMesh, 256)){

                                    if(symmetryPlaneSet.size() == 1){
                                        auto symmetryPlane = symmetryPlaneSet[0];

                                        std::cout << "Symmetry plane normal: " << symmetryPlane.getNormal().x << ", " << symmetryPlane.getNormal().y << ", " << symmetryPlane.getNormal().z << std::endl;
                                        auto n = symmetryPlane.getNormal();

                                        auto v = n;
                                        model.addConstr(n.x * (R[0] * v.x + R[1] * v.y + R[2] * v.z) +
                                                        n.y * (R[3] * v.x + R[4] * v.y + R[5] * v.z) +
                                                        n.z * (R[6] * v.x + R[7] * v.y + R[8] * v.z) >= 0, "Item symmetry breaking constraint");
                                    }
                                    else if(symmetryPlaneSet.size() == 2){

                                        auto symmetryPlaneA = symmetryPlaneSet[0];
                                        auto symmetryPlaneB = symmetryPlaneSet[1];

                                        std::cout << "Symmetry plane normal A: " << symmetryPlaneA.getNormal().x << ", " << symmetryPlaneA.getNormal().y << ", " << symmetryPlaneA.getNormal().z << std::endl;
                                        std::cout << "Symmetry plane normal B: " << symmetryPlaneB.getNormal().x << ", " << symmetryPlaneB.getNormal().y << ", " << symmetryPlaneB.getNormal().z << std::endl;

                                        auto nA = symmetryPlaneA.getNormal();
                                        auto nB = symmetryPlaneB.getNormal();

                                        auto v = nA + nB;

                                        model.addConstr(nA.x * (R[0] * v.x + R[1] * v.y + R[2] * v.z) +
                                                        nA.y * (R[3] * v.x + R[4] * v.y + R[5] * v.z) +
                                                        nA.z * (R[6] * v.x + R[7] * v.y + R[8] * v.z) >= 0, "Item symmetry breaking constraint");

                                        model.addConstr(nB.x * (R[0] * v.x + R[1] * v.y + R[2] * v.z) +
                                                        nB.y * (R[3] * v.x + R[4] * v.y + R[5] * v.z) +
                                                        nB.z * (R[6] * v.x + R[7] * v.y + R[8] * v.z) >= 0, "Item symmetry breaking constraint");
                                    }
                                }
                            }

                            assert(itemModelMesh->getConvexHull());

                            // Constrain vertices to be at positive side of each facet using the plane equation and transformation
                            for (const auto &itemVertex: itemModelMesh->getConvexHull()->getVertices()){
                                // Because we are solely considering convex containers, just the vertices that are part of the convex hull are sufficient
                                // This way we can (greatly) reduce the amount of constraints
                                const auto& containerVertices = convexContainerHull->getVertices();
                                for (const auto &containerTriangle: convexContainerHull->getTriangles()){
                                    VertexTriangle triangle{containerVertices[containerTriangle.vertexIndex0], containerVertices[containerTriangle.vertexIndex1], containerVertices[containerTriangle.vertexIndex2]};

                                    // Plane equation with transformed coordinates
                                    auto n = triangle.normal;
                                    auto d = -glm::dot(n, triangle.vertices[0]);
                                    auto v = itemVertex;
                                    model.addConstr(n.x * (SR[0] * v.x + SR[1] * v.y + SR[2] * v.z + Tx) +
                                                    n.y * (SR[3] * v.x + SR[4] * v.y + SR[5] * v.z + Ty) +
                                                    n.z * (SR[6] * v.x + SR[7] * v.y + SR[8] * v.z + Tz) + d <= -1e-8, "Vertex on negative side of facet.");
                                }
                            }

                            // Optimize model
                            model.optimize();

                            // Extract the achieved scale and update the upper bound
                            float achievedScale = S.get(GRB_DoubleAttr_X);
                            std::cout << "Original upper bound was: " << upperScaleBound << std::endl;
                            std::cout << "Upper bound found in relaxation: " << achievedScale << std::endl;
                            upperScaleBound = std::min(upperScaleBound, achievedScale);

                            // Extract the time spent in the first model
                            timeInFirstModel = model.get(GRB_DoubleAttr_Runtime);


                        } catch (GRBException &e) {
                            std::cerr << "Error code = " << e.getErrorCode() << std::endl;
                            std::cerr << e.getMessage() << std::endl;
                        }
                    }

                    try {
                        // Create an environment
                        GRBEnv env = GRBEnv(true);
                        env.set("LogFile", "concave.log");
                        env.start();

                        // Create an empty model
                        GRBModel model = GRBModel(env);
                        model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
                        model.set(GRB_IntParam_Threads, 64);
                        model.set(GRB_IntParam_NonConvex, 2);

                        // Limit execution time to 3600s
                        model.set(GRB_DoubleParam_TimeLimit, 10 * 3600 - timeInFirstModel);

                        // Initialize decision variables
                        GRBVar S = model.addVar(0.0, upperScaleBound, 1.0, GRB_CONTINUOUS, "S");

                        GRBVar Tx = model.addVar(-std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), 0, GRB_CONTINUOUS, "Tx");
                        GRBVar Ty = model.addVar(-std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), 0, GRB_CONTINUOUS, "Ty");
                        GRBVar Tz = model.addVar(-std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), 0, GRB_CONTINUOUS, "Tz");

                        GRBVar Qw = model.addVar(-1, 1, 0, GRB_CONTINUOUS, "Qw");
                        GRBVar Qx = model.addVar(-1, 1, 0, GRB_CONTINUOUS, "Qx");
                        GRBVar Qy = model.addVar(-1, 1, 0, GRB_CONTINUOUS, "Qy");
                        GRBVar Qz = model.addVar(-1, 1, 0, GRB_CONTINUOUS, "Qz");

                        // Constrain the quaternion to be normalized
                        model.addQConstr(Qw * Qw + Qx * Qx + Qy * Qy + Qz * Qz == 1, "|Q| = 1");

                        // Support variables
                        std::array<GRBVar, 9> R; // Rotation matrix
                        std::array<GRBVar, 9> SR; // Scaling + Rotation matrix
                        for (int i = 0; i < 9; i++) {
                            R[i] = model.addVar(-1, 1, 0, GRB_CONTINUOUS, "R[" + std::to_string(i) + "]");
                            SR[i] = model.addVar(-upperScaleBound, upperScaleBound, 0,  GRB_CONTINUOUS, "SR[" + std::to_string(i) + "]");
                            model.addQConstr(SR[i] == S * R[i], "SR[" + std::to_string(i) + "] = S * R[" + std::to_string(i) + "]");
                        }

                        model.addQConstr(R[0] == 1 - 2 * (Qy * Qy + Qz * Qz), "mat[0]  = 1 - 2 * ( yy + zz )");
                        model.addQConstr(R[1] == 2 * (Qx * Qy - Qz * Qw), "mat[1]  = 2 * ( xy - zw )");
                        model.addQConstr(R[2] == 2 * (Qx * Qz + Qy * Qw), "mat[2]  = 2 * ( xz + yw )");
                        model.addQConstr(R[3] == 2 * (Qx * Qy + Qz * Qw), "mat[3]  = 2 * ( xy + zw )");
                        model.addQConstr(R[4] == 1 - 2 * (Qx * Qx + Qz * Qz), "mat[4]  = 1 - 2 * ( xx + zz )");
                        model.addQConstr(R[5] == 2 * (Qy * Qz - Qx * Qw), "mat[5]  = 2 * ( yz - xw )");
                        model.addQConstr(R[6] == 2 * (Qx * Qz - Qy * Qw), "mat[6]  = 2 * ( xz - yw )");
                        model.addQConstr(R[7] == 2 * (Qy * Qz + Qx * Qw), "mat[7]  = 2 * ( yz + xw )");
                        model.addQConstr(R[8] == 1 - 2 * (Qx * Qx + Qy * Qy), "mat[8]  = 1 - 2 * ( xx + yy )");

                        assert(itemModelMesh->getConvexHull() && itemModelMesh->getConvexHull()->isConvex());
                        assert(containerModelSpaceMesh->getConvexHull() && containerModelSpaceMesh->getConvexHull()->isConvex());

                        // Symmetry breaking constraints
                        {

                            auto v = itemModelMesh->getVolumeCentroid();

                            for (const auto &symmetryBreakingPlaneSet: SymmetryBreaker::deriveSymmetryBreakingPlanes(containerModelSpaceMesh, 256)){
                                for (const auto &symmetryBreakingPlane: symmetryBreakingPlaneSet){

                                    // Force position to be at positive side of the planes
                                    auto n = symmetryBreakingPlane.getNormal();
                                    auto d = symmetryBreakingPlane.getD();
                                    model.addConstr(n.x * (SR[0] * v.x + SR[1] * v.y + SR[2] * v.z + Tx) +
                                                    n.y * (SR[3] * v.x + SR[4] * v.y + SR[5] * v.z + Ty) +
                                                    n.z * (SR[6] * v.x + SR[7] * v.y + SR[8] * v.z + Tz) + d >= 0, "Container symmetry breaking constraint");
                                }
                            }
                        }

                        {
                            for (const auto &symmetryPlaneSet: SymmetryBreaker::deriveSymmetryBreakingPlanes(itemModelMesh, 256)){

                                if(symmetryPlaneSet.size() == 1){
                                    auto symmetryPlane = symmetryPlaneSet[0];


                                    std::cout << "Symmetry plane normal: " << symmetryPlane.getNormal().x << ", " << symmetryPlane.getNormal().y << ", " << symmetryPlane.getNormal().z << std::endl;
                                    auto n = symmetryPlane.getNormal();

                                    auto v = n;
                                    model.addConstr(n.x * (R[0] * v.x + R[1] * v.y + R[2] * v.z) +
                                                    n.y * (R[3] * v.x + R[4] * v.y + R[5] * v.z) +
                                                    n.z * (R[6] * v.x + R[7] * v.y + R[8] * v.z) >= 0, "Item symmetry breaking constraint");
                                }
                                else if(symmetryPlaneSet.size() == 2){

                                    auto symmetryPlaneA = symmetryPlaneSet[0];
                                    auto symmetryPlaneB = symmetryPlaneSet[1];

                                    auto nA = symmetryPlaneA.getNormal();
                                    auto nB = symmetryPlaneB.getNormal();

                                    auto v = nA + nB;

                                    model.addConstr(nA.x * (R[0] * v.x + R[1] * v.y + R[2] * v.z) +
                                                    nA.y * (R[3] * v.x + R[4] * v.y + R[5] * v.z) +
                                                    nA.z * (R[6] * v.x + R[7] * v.y + R[8] * v.z) >= 0, "Item symmetry breaking constraint");

                                    model.addConstr(nB.x * (R[0] * v.x + R[1] * v.y + R[2] * v.z) +
                                                    nB.y * (R[3] * v.x + R[4] * v.y + R[5] * v.z) +
                                                    nB.z * (R[6] * v.x + R[7] * v.y + R[8] * v.z) >= 0, "Item symmetry breaking constraint");
                                }
                            }
                        }

                        for (const auto &convexContainerTriangle: containerModelSpaceMesh->getConvexHull()->getTriangles()){
                            const auto& containerVertices = containerModelSpaceMesh->getConvexHull()->getVertices();

                            VertexTriangle triangle{containerVertices[convexContainerTriangle.vertexIndex0], containerVertices[convexContainerTriangle.vertexIndex1], containerVertices[convexContainerTriangle.vertexIndex2]};


                            // Plane equation with transformed coordinates
                            auto n = triangle.normal;

                            if(glm::length(n) <= 1e-8f){
                                continue;
                            }

                            n = glm::normalize(n);
                            auto d = -glm::dot(n, triangle.vertices[0]);


                            // Constrain vertices to be at positive side of each facet using the plane equation and transformation
                            for (const auto &v: itemModelMesh->getConvexHull()->getVertices()){
                                model.addConstr(n.x * (SR[0] * v.x + SR[1] * v.y + SR[2] * v.z + Tx) +
                                                n.y * (SR[3] * v.x + SR[4] * v.y + SR[5] * v.z + Ty) +
                                                n.z * (SR[6] * v.x + SR[7] * v.y + SR[8] * v.z + Tz) + d <= -1e-8, "Vertices should be on the negative side of a facet that is part of the convex hull.");
                            }
                        }

                        std::vector<std::shared_ptr<ModelSpaceMesh>> concavities = ConvexConcavitiesFactory::getConvexConcavities(containerModelSpaceMesh);
                        std::vector<std::vector<Plane>> groupedConcavityPlanes = ConvexConcavitiesFactory::getGroupedConcavityPlanes(containerModelSpaceMesh);
                        if(!concavities.empty()){

                            for (int g = 0; g < concavities.size(); ++g){

                                const auto& concavity = concavities[g];

                                float NxMax;
                                float NxMin;
                                float NyMax;
                                float NyMin;
                                float NzMax;
                                float NzMin;
                                float DMax;
                                float DMin;

                                if(limitSeparationPlanes){

                                    const auto& concavityPlanes = groupedConcavityPlanes[g];

                                    assert(concavityPlanes.size() > 1);

                                    if(concavityPlanes.empty()){
                                        continue;
                                    }

                                    // Check if normalized
                                    assert(glm::epsilonEqual(glm::length(concavityPlanes[0].getNormal()), 1.0f, 1e-6f));

                                    NxMax = concavityPlanes[0].getNormal().x;
                                    NxMin = concavityPlanes[0].getNormal().x;
                                    NyMax = concavityPlanes[0].getNormal().y;
                                    NyMin = concavityPlanes[0].getNormal().y;
                                    NzMax = concavityPlanes[0].getNormal().z;
                                    NzMin = concavityPlanes[0].getNormal().z;
                                    for (int i = 1; i < concavityPlanes.size(); ++i){
                                        NxMax = std::max(NxMax, concavityPlanes[i].getNormal().x);
                                        NxMin = std::min(NxMin, concavityPlanes[i].getNormal().x);
                                        NyMax = std::max(NyMax, concavityPlanes[i].getNormal().y);
                                        NyMin = std::min(NyMin, concavityPlanes[i].getNormal().y);
                                        NzMax = std::max(NzMax, concavityPlanes[i].getNormal().z);
                                        NzMin = std::min(NzMin, concavityPlanes[i].getNormal().z);
                                    }

                                    DMax = 0.0f;
                                    DMin = 0.0f;
                                    for (const auto& concavityPlane: concavityPlanes){
                                        DMax += std::max(0.0f, concavityPlane.getD());
                                        DMin += std::min(0.0f, concavityPlane.getD());
                                    }

                                    std::cout << "Nx-" + std::to_string(g) << " \t[" << NxMin << ", " << NxMax << "]" << std::endl;
                                    std::cout << "Ny-" + std::to_string(g) << "\t[" << NyMin << ", " << NyMax << "]" << std::endl;
                                    std::cout << "Nz-" + std::to_string(g) << "\t[" << NzMin << ", " << NzMax << "]" << std::endl;
                                    std::cout << "D-" + std::to_string(g) << "\t[" << DMin << ", " << DMax << "]" << std::endl;
                                }
                                else {

                                    NxMax = 1;
                                    NxMin = -1;
                                    NyMax = 1;
                                    NyMin = -1;
                                    NzMax = 1;
                                    NzMin = -1;
                                    DMax = std::numeric_limits<float>::max();
                                    DMin = -std::numeric_limits<float>::max();

                                }

                                // Define a separating plane between the concavity and the convex item
                                // The normal direction of this plane is expected to be a convex composition of the concavity planes, so we derive bounds from these extreme values
                                GRBVar Nx = model.addVar(NxMin, NxMax, 0, GRB_CONTINUOUS, "Nx-" + std::to_string(g));
                                GRBVar Ny = model.addVar(NyMin, NyMax, 0, GRB_CONTINUOUS, "Ny-" + std::to_string(g));
                                GRBVar Nz = model.addVar(NzMin, NzMax, 0, GRB_CONTINUOUS, "Nz-" + std::to_string(g));
                                GRBVar D = model.addVar(DMin, DMax, 0, GRB_CONTINUOUS, "D-" + std::to_string(g));

                                if(limitSeparationPlanes){
                                    // Force N to be non zero
                                    model.addQConstr(Nx * Nx + Ny * Ny + Nz * Nz <= 1, "|N| <= 1");
                                    model.addQConstr(Nx * Nx + Ny * Ny + Nz * Nz >= 1e-4, "|N| >= 1e-4");
                                }
                                else{
                                    // Force N to be normalized
                                    model.addQConstr(Nx * Nx + Ny * Ny + Nz * Nz == 1, "|N| == 1");
                                }

                                // Force concavity vertices to be on the negative side of the separating plane
                                for (const auto &v: concavity->getVertices()){
                                    model.addConstr(Nx * v.x + Ny * v.y + Nz * v.z + D <= -1e-8, "Vertices should be on the negative side of a separating plane between the concavity and the convex item.");
                                }

                                // Force all transformed item vertices to be on the positive side

                                // Precompute products of (Nx, Ny, Nz) and transformation matrix in separate variables,
                                // they are identical for all vertices and will allow linear plane equations
                                // This also allows us to manually put sensible bounds on them
                                auto NxSR0 = model.addVar(-upperScaleBound, upperScaleBound, 0, GRB_CONTINUOUS, "NxSR0-" + std::to_string(g));
                                auto NxSR1 = model.addVar(-upperScaleBound, upperScaleBound, 0, GRB_CONTINUOUS, "NxSR1-" + std::to_string(g));
                                auto NxSR2 = model.addVar(-upperScaleBound, upperScaleBound, 0, GRB_CONTINUOUS, "NxSR2-" + std::to_string(g));
                                auto NySR3 = model.addVar(-upperScaleBound, upperScaleBound, 0, GRB_CONTINUOUS, "NySR3-" + std::to_string(g));
                                auto NySR4 = model.addVar(-upperScaleBound, upperScaleBound, 0, GRB_CONTINUOUS, "NySR4-" + std::to_string(g));
                                auto NySR5 = model.addVar(-upperScaleBound, upperScaleBound, 0, GRB_CONTINUOUS, "NySR5-" + std::to_string(g));
                                auto NzSR6 = model.addVar(-upperScaleBound, upperScaleBound, 0, GRB_CONTINUOUS, "NzSR6-" + std::to_string(g));
                                auto NzSR7 = model.addVar(-upperScaleBound, upperScaleBound, 0, GRB_CONTINUOUS, "NzSR7-" + std::to_string(g));
                                auto NzSR8 = model.addVar(-upperScaleBound, upperScaleBound, 0, GRB_CONTINUOUS, "NzSR8-" + std::to_string(g));
                                auto NxTx = model.addVar(-std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), 0, GRB_CONTINUOUS, "NxTx-" + std::to_string(g));
                                auto NyTy = model.addVar(-std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), 0, GRB_CONTINUOUS, "NyTy-" + std::to_string(g));
                                auto NzTz = model.addVar(-std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), 0, GRB_CONTINUOUS, "NzTz-" + std::to_string(g));

                                model.addQConstr(NxSR0 == SR[0] * Nx, "NxSR0-" + std::to_string(g));
                                model.addQConstr(NxSR1 == SR[1] * Nx, "NxSR1-" + std::to_string(g));
                                model.addQConstr(NxSR2 == SR[2] * Nx, "NxSR2-" + std::to_string(g));
                                model.addQConstr(NxTx  == Tx    * Nx, "NxTx-"  + std::to_string(g));
                                model.addQConstr(NySR3 == SR[3] * Ny, "NySR3-" + std::to_string(g));
                                model.addQConstr(NySR4 == SR[4] * Ny, "NySR4-" + std::to_string(g));
                                model.addQConstr(NySR5 == SR[5] * Ny, "NySR5-" + std::to_string(g));
                                model.addQConstr(NyTy  == Ty    * Ny, "NyTy-"  + std::to_string(g));
                                model.addQConstr(NzSR6 == SR[6] * Nz, "NzSR6-" + std::to_string(g));
                                model.addQConstr(NzSR7 == SR[7] * Nz, "NzSR7-" + std::to_string(g));
                                model.addQConstr(NzSR8 == SR[8] * Nz, "NzSR8-" + std::to_string(g));
                                model.addQConstr(NzTz  == Tz    * Nz, "NzTz-"  + std::to_string(g));

                                for (const auto &v: itemModelMesh->getConvexHull()->getVertices()){
                                    model.addConstr(NxSR0 * v.x + NxSR1 * v.y + NxSR2 * v.z + NxTx +
                                                    NySR3 * v.x + NySR4 * v.y + NySR5 * v.z + NyTy +
                                                    NzSR6 * v.x + NzSR7 * v.y + NzSR8 * v.z + NzTz + D >= 1e-8, "Vertices should be on the positive side of a separating plane between the concavity and the convex item.");
                                }
                            }
                        }
                        // Optimize model
                        model.optimize();


                        auto totalRunTime = model.get(GRB_DoubleAttr_Runtime) + timeInFirstModel;
                        std::cout << "Total Gurobi runtime: " << totalRunTime << "s" << std::endl;


                        // Write the results to the .csv file
                        auto s = std::to_string(S.get(GRB_DoubleAttr_X));
                        auto gap = std::to_string(model.get(GRB_DoubleAttr_MIPGap));
                        auto explored = std::to_string(model.get(GRB_DoubleAttr_NodeCount));
                        auto unexplored = std::to_string(model.get(GRB_DoubleAttr_NodeCount));
                        auto tx = std::to_string(Tx.get(GRB_DoubleAttr_X));
                        auto ty = std::to_string(Ty.get(GRB_DoubleAttr_X));
                        auto tz = std::to_string(Tz.get(GRB_DoubleAttr_X));
                        auto qw = std::to_string(Qw.get(GRB_DoubleAttr_X));
                        auto qx = std::to_string(Qx.get(GRB_DoubleAttr_X));
                        auto qy = std::to_string(Qy.get(GRB_DoubleAttr_X));
                        auto qz = std::to_string(Qz.get(GRB_DoubleAttr_X));

                        csvFile << useConvexRelaxation << "," << limitSeparationPlanes << "," << containerModelSpaceMesh->getName() << "," << itemModelMesh->getName() << "," << s << "," << totalRunTime << "," << gap << "," << explored << "," << unexplored << "," << tx << "," << ty << "," << tz << "," << qw << "," << qx << "," << qy << "," << qz << std::endl;

                        // Flush the output
                        csvFile.flush();

                    }
                    catch (GRBException &e) {
                        std::cout << "Error code = " << e.getErrorCode() << std::endl;
                        std::cout << e.getMessage() << std::endl;
                    }
                    catch (...) {
                        std::cout << "Exception during optimization" << std::endl;
                    }
                }
            }
        }
    }
}
