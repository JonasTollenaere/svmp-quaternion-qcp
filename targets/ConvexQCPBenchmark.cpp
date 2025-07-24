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

static std::string dataFolder = MESHCORE_DATA_DIR + std::string("Tollenaere, J. et al/");

auto containerFiles = {
        "Containers/convexstone.stl",
        "Containers/dodecahedron.stl",
        "Containers/frustum.stl",
        "Containers/tetrahedron.stl",
};
auto itemFiles = {
        "Items/apple.stl",
        "Items/arrow.obj",
        "Items/banana.stl",
        "Items/bobbin.stl",
        "Items/convexstone.stl",
        "Items/cube.obj",
        "Items/dagger.stl",
        "Items/dodecahedron.stl",
        "Items/dragon.obj",
        "Items/frustum.stl",
        "Items/gem.stl",
        "Items/helmet.stl",
        "Items/oloid.stl",
        "Items/pig.obj",
        "Items/ring.obj",
        "Items/star.obj",
        "Items/star.stl",
        "Items/tetrahedron.stl",
        "Items/torus.stl",
        "Items/tree.stl"
};

int main() {

    // Open a .csv writerhr
    std::ofstream csvFile;
    csvFile.open("QuaternionQCPResultsInverseScalingVolume.csv");

    // Write the header
    csvFile << "Convex hull?,Break symmetry?,Container,Item,Time(s),Scale,Gap,Explored,Unexplored,TranslationX,TranslationY,TranslationZ,RotationW,RotationX,RotationY,RotationZ,Volume" << std::endl;
    csvFile.flush();

    auto useConvexHulls = {true};
    auto breakSymmetries = {true};

    for (const auto &useConvexHull: useConvexHulls) {
        for(const auto& breakSymmetry: breakSymmetries) {
            for (const auto &containerFile: containerFiles) {
                for (const auto &itemFile: itemFiles) {

                    std::cout << "Processing: " << itemFile << " in " << containerFile << ". Convex hull: " << useConvexHull << ". Break symmetry: " << breakSymmetry << std::endl;

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
                    double upperScaleBound = UpperBoundProvider::getUpperScaleBoundConvexHullItem(itemModelMesh, containerModelSpaceMesh);

                    std::cout << "Processing: " << itemModelMesh->getName() << " in " << containerModelSpaceMesh->getName() << std::endl;

                    try {
                        // Create an environment
                        GRBEnv env = GRBEnv(true);
                        env.set("LogFile", "quat.log");
                        env.start();

                        // Create an empty model
                        GRBModel model = GRBModel(env);
                        model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
                        model.set(GRB_IntParam_Threads, 64); // This seems to benefit our model, even on the 12-core M2 pro (In that case, defaults to 12 threads. Processor usage is also higher with 64) (e.g. 75 seconds instead of 120 seconds)
                        model.set(GRB_IntParam_NonConvex, 2);

                        // Limit execution time to 3600s
                        model.set(GRB_DoubleParam_TimeLimit, 3600);

                        // Initialize decision variables
                        GRBVar inverseS = model.addVar(1.0/upperScaleBound, 10000, 1.0, GRB_CONTINUOUS, "S");

                        GRBVar Tx = model.addVar(-std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
                                                 0, GRB_CONTINUOUS, "Tx");
                        GRBVar Ty = model.addVar(-std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
                                                 0, GRB_CONTINUOUS, "Ty");
                        GRBVar Tz = model.addVar(-std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
                                                 0, GRB_CONTINUOUS, "Tz");

                        GRBVar Qw = model.addVar(-1, 1, 0, GRB_CONTINUOUS, "Qw");
                        GRBVar Qx = model.addVar(-1, 1, 0, GRB_CONTINUOUS, "Qx");
                        GRBVar Qy = model.addVar(-1, 1, 0, GRB_CONTINUOUS, "Qy");
                        GRBVar Qz = model.addVar(-1, 1, 0, GRB_CONTINUOUS, "Qz");

                        // Constrain the quaternion to be normalized
                        model.addQConstr(Qw * Qw + Qx * Qx + Qy * Qy + Qz * Qz == 1, "|Q| = 1");

                        // Support variables
                        std::array<GRBVar, 9> R; // Rotation matrix
                        for (int i = 0; i < 9; i++) {
                            R[i] = model.addVar(-1, 1, 0, GRB_CONTINUOUS, "R[" + std::to_string(i) + "]");
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
                        if(breakSymmetry){

                            auto v = itemModelMesh->getVolumeCentroid();
                            //auto v = itemModelMesh->getVertices()[0];

                            for (const auto &symmetryBreakingPlaneSet: SymmetryBreaker::deriveSymmetryBreakingPlanes(containerModelSpaceMesh, 256)){
                                for (const auto &symmetryBreakingPlane: symmetryBreakingPlaneSet){

                                    // Force position to be at positive side of the planes
                                    auto n = symmetryBreakingPlane.getNormal();
                                    auto d = symmetryBreakingPlane.getD();
                                    model.addConstr(n.x * (R[0] * v.x + R[1] * v.y + R[2] * v.z + Tx) +
                                                    n.y * (R[3] * v.x + R[4] * v.y + R[5] * v.z + Ty) +
                                                    n.z * (R[6] * v.x + R[7] * v.y + R[8] * v.z + Tz) + inverseS * d >= 0, "Container symmetry breaking constraint");
                                }
                            }
                        }

                        if(breakSymmetry){

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

                        assert(itemModelMesh->getConvexHull());
                        assert(containerModelSpaceMesh->isConvex());

                        auto vertices = useConvexHull ? itemModelMesh->getConvexHull()->getVertices() : itemModelMesh->getVertices();

                        // Constrain vertices to be at positive side of each facet using the plane equation and transformation
                        for (const auto &itemVertex: vertices) {
                            // Because we are solely considering convex containers, just the vertices that are part of the convex hull are sufficient
                            // This way we can (greatly) reduce the amount of constraints
                            const auto &containerVertices = containerModelSpaceMesh->getVertices();
                            for (const auto &containerTriangle: containerModelSpaceMesh->getTriangles()) {
                                VertexTriangle triangle{containerVertices[containerTriangle.vertexIndex0],
                                                        containerVertices[containerTriangle.vertexIndex1],
                                                        containerVertices[containerTriangle.vertexIndex2]};

                                // Plane equation with transformed coordinates
                                auto n = triangle.normal;
                                auto d = -glm::dot(n, triangle.vertices[0]);
                                auto v = itemVertex;
                                model.addConstr(n.x * (R[0] * v.x + R[1] * v.y + R[2] * v.z + Tx) +
                                                n.y * (R[3] * v.x + R[4] * v.y + R[5] * v.z + Ty) +
                                                n.z * (R[6] * v.x + R[7] * v.y + R[8] * v.z + Tz) + inverseS*d <= -1e-8,
                                                "Vertex on negative side of facet.");
                            }
                        }

                        // Optimize model
                        model.optimize();

                        // Write the results to the .csv file
                        auto volume = itemModelMesh->getVolume() * std::pow(1.0/inverseS.get(GRB_DoubleAttr_X),3.0);
                        csvFile << useConvexHull << "," << breakSymmetry << "," << containerModelSpaceMesh->getName() << "," << itemModelMesh->getName() << "," << model.get(GRB_DoubleAttr_Runtime) << "," << 1.0/inverseS.get(GRB_DoubleAttr_X) << "," << model.get(GRB_DoubleAttr_MIPGap) << "," << model.get(GRB_DoubleAttr_NodeCount) << "," << model.get(GRB_DoubleAttr_NodeCount) << "," << Tx.get(GRB_DoubleAttr_X) / inverseS.get(GRB_DoubleAttr_X) << "," << Ty.get(GRB_DoubleAttr_X) / inverseS.get(GRB_DoubleAttr_X) << "," << Tz.get(GRB_DoubleAttr_X) / inverseS.get(GRB_DoubleAttr_X) << "," << Qw.get(GRB_DoubleAttr_X) << "," << Qx.get(GRB_DoubleAttr_X) << "," << Qy.get(GRB_DoubleAttr_X) << "," << Qz.get(GRB_DoubleAttr_X) << "," << volume << std::endl;

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
