//
// Created by Jonas on 19/03/2024.
//

#include "gurobi_c++.h"

#include <meshcore/tasks/AbstractTask.h>
#include <meshcore/utility/FileParser.h>
#include <meshcore/utility/io.h>
#include <meshcore/optimization/SingleVolumeMaximisationSolution.h>
#include <meshcore/rendering/ApplicationWindow.h>

#include <nlohmann/json.hpp>
#include <iostream>
#include <string>
#include <glm/gtx/component_wise.hpp>

#include "SymmetryBreaker.h"
#include "UpperBoundProvider.h"

#define SYMMETRY_BREAKING true

class Task: public AbstractTask {

    class Callback: public GRBCallback {
    public:

        std::shared_ptr<SingleVolumeMaximisationSolution> solution;

        GRBVar* S;
        GRBVar* Tx;
        GRBVar* Ty;
        GRBVar* Tz;
        GRBVar* Qw;
        GRBVar* Qx;
        GRBVar* Qy;
        GRBVar* Qz;

        Task* task;

        Callback(Task* task, const std::shared_ptr<SingleVolumeMaximisationSolution>& callbackSolution, GRBVar* S, GRBVar* Tx, GRBVar* Ty, GRBVar* Tz, GRBVar* Qw, GRBVar* Qx, GRBVar* Qy, GRBVar* Qz) : S(S), Tx(Tx), Ty(Ty), Tz(Tz), Qw(Qw), Qx(Qx), Qy(Qy), Qz(Qz), solution(callbackSolution), task(task) {}

    protected:
        void callback() override {
            if (where == GRB_CB_MIPSOL) {

                Quaternion rotation = Quaternion(glm::fquat(getSolution(*Qw), getSolution(*Qx), getSolution(*Qy), getSolution(*Qz)));
                Transformation transformation(rotation);
                transformation.setPositionX(getSolution(*Tx));
                transformation.setPositionY(getSolution(*Ty));
                transformation.setPositionZ(getSolution(*Tz));
                transformation.setScale(getSolution(*S));
                solution->getItemWorldSpaceMesh()->setModelTransformation(transformation);

                task->notifyObserversSolution(solution);
            }
        }
    };

    void run() override {

        this->notifyObserversStatus("Initialising");

        auto containerPath = MESHCORE_DATA_DIR + std::string("Tollenaere, J. et al/Items/convexstone.stl");
        auto itemPath = MESHCORE_DATA_DIR + std::string("Tollenaere, J. et al/Items/star.stl");

        auto containerModelSpaceMesh = FileParser::loadMeshFile(containerPath);
        auto containerMesh = std::make_shared<WorldSpaceMesh>(containerModelSpaceMesh);
        auto itemMesh = std::make_shared<WorldSpaceMesh>(FileParser::loadMeshFile(itemPath));

        // Center around volume centroid
        itemMesh->getModelTransformation().setPosition(-itemMesh->getModelSpaceMesh()->getVolumeCentroid());
        itemMesh = std::make_shared<WorldSpaceMesh>(itemMesh->getTransformedModelSpaceMesh());
        auto itemModelMesh = itemMesh->getModelSpaceMesh();

        auto solution = std::make_shared<SingleVolumeMaximisationSolution>(itemMesh, containerMesh);
        solution->getItemWorldSpaceMesh()->getModelTransformation().setScale(0.0f);
        notifyObserversSolution(solution);
        double upperScaleBound = UpperBoundProvider::getUpperScaleBoundConvexHullItem(itemModelMesh, containerModelSpaceMesh);

        std::cout << "Processing: " << itemModelMesh->getName() << " in " << containerModelSpaceMesh->getName() << std::endl;

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

            // Initialize decision variables
            GRBVar S = model.addVar(0, upperScaleBound, 1.0, GRB_CONTINUOUS, "S");

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

#if SYMMETRY_BREAKING
            // Container symmetry breaking constraints
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

            // Item symmetry breaking constraints
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
#endif

            assert(itemModelMesh->getConvexHull());
            assert(containerModelSpaceMesh->isConvex());

            // Constrain vertices to be at positive side of each facet using the plane equation and transformation
            for (const auto &itemVertex: itemModelMesh->getConvexHull()->getVertices()){
                // Because we are solely considering convex containers, just the vertices that are part of the convex hull are sufficient
                // This way we can (greatly) reduce the amount of constraints
                const auto& containerVertices = containerModelSpaceMesh->getVertices();
                for (const auto &containerTriangle: containerModelSpaceMesh->getTriangles()){
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

            // Configure callback
            Callback cb(this, solution, &S, &Tx, &Ty, &Tz, &Qw, &Qx, &Qy, &Qz);
            model.setCallback(&cb);

            // Optimize model
            notifyObserversStatus("Optimizing");
            model.optimize();

            Quaternion rotation = Quaternion(glm::fquat(Qw.get(GRB_DoubleAttr_X), Qx.get(GRB_DoubleAttr_X), Qy.get(GRB_DoubleAttr_X), Qz.get(GRB_DoubleAttr_X)));

            Transformation transformation(rotation);
            transformation.setScale(S.get(GRB_DoubleAttr_X));
            transformation.setPositionX(Tx.get(GRB_DoubleAttr_X));
            transformation.setPositionY(Ty.get(GRB_DoubleAttr_X));
            transformation.setPositionZ(Tz.get(GRB_DoubleAttr_X));
            solution->getItemWorldSpaceMesh()->setModelTransformation(transformation);

            std::cout << transformation << std::endl;

            notifyObserversSolution(solution);

        } catch(GRBException& e) {
            std::cout << "Error code = " << e.getErrorCode() << std::endl;
            std::cout << e.getMessage() << std::endl;
        } catch(...) {
            std::cout << "Exception during optimization" << std::endl;
        }
    }
};

int main(int argc, char *argv[]) {

    Task task;

    QApplication app(argc, argv);
    ApplicationWindow window;
    window.show();
    RenderWidget *renderWidget = window.getRenderWidget();

    renderWidget->observeTask(&task);
    task.start();
    int returnCode = QApplication::exec();
    task.stop();
    task.join();
    return returnCode;
}

