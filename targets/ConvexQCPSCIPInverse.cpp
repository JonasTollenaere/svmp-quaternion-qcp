//
// Created by Jonas on 22/04/2024.
//

#include "scip/scip.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_quadratic.h"
#include "scip/scipdefplugins.h"

#include "meshcore/tasks/AbstractTask.h"

#include <nlohmann/json.hpp>
#include "meshcore/utility/FileParser.h"
#include "meshcore/utility/io.h"
#include <iostream>
#include <string>
#include <glm/gtx/component_wise.hpp>
#include "InflatingLateAcceptance.h"
#include <QtWidgets>
#include "meshcore/rendering/ApplicationWindow.h"
#include "../solutions/SingleItemSolution.h"
#include "SymmetryBreaker.h"
#include "UpperBoundProvider.h"

class Task: public AbstractTask {

    int exec() {

        this->notifyObserversStatus("Initialising");

        auto containerPath = "../../../meshcore/datasets/Tollenaere, J. et al/Items/gem.stl";
        auto itemPath = "../../../meshcore/datasets/Tollenaere, J. et al/Items/dragon.obj";

        auto containerModelSpaceMesh = FileParser::loadMeshFile(containerPath);
        auto containerMesh = std::make_shared<WorldSpaceMesh>(containerModelSpaceMesh);
        auto itemMesh = std::make_shared<WorldSpaceMesh>(FileParser::loadMeshFile(itemPath));

        // Center around volume centroid
        itemMesh->getModelTransformation().setPosition(-itemMesh->getModelSpaceMesh()->getVolumeCentroid());
        std::cout << "Normalizing transformation on item" << itemMesh->getModelTransformation() << std::endl;
        itemMesh = std::make_shared<WorldSpaceMesh>(itemMesh->getTransformedModelSpaceMesh());
        auto itemModelMesh = itemMesh->getModelSpaceMesh();

        auto solution = std::make_shared<SingleItemSolution>(containerMesh, itemMesh);
        solution->getItemModelTransformation().setScale(0.0f);
        notifyObserversSolution(solution);
        double upperScaleBound = UpperBoundProvider::getUpperScaleBoundConvexHullItem(itemModelMesh, containerModelSpaceMesh);


        std::cout << "Processing: " << itemModelMesh->getName() << " in " << containerModelSpaceMesh->getName() << std::endl;
        std::cout << "Upper scale bound: " << upperScaleBound << std::endl;

        try {
            // Create an environment
            SCIP* scip = nullptr;

            /* initialize SCIP environment */
            SCIP_CALL(SCIPcreate(&scip));
            SCIP_CALL(SCIPincludeDefaultPlugins(scip));

            // Create an empty model
            SCIP_CALL(SCIPcreateProbBasic(scip, "QuadraticExample"));
            SCIP_CALL(SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE));

            // Enable concurrency
//            SCIP_CALL(SCIPsetIntParam(scip, "parallel/minnthreads", 1) );
//            SCIP_CALL(SCIPsetIntParam(scip, "parallel/maxnthreads", 20) );

            // Initialize scaling variable
            SCIP_VAR* inverseS = nullptr;
            SCIP_CALL(SCIPcreateVarBasic(scip, &inverseS, "inverseS", 1.0/upperScaleBound, 1000, 1.0, SCIP_VARTYPE_CONTINUOUS));
            SCIP_CALL(SCIPaddVar(scip, inverseS));

            // Initialize translation variables
            SCIP_VAR* Tx = nullptr;
            SCIP_VAR* Ty = nullptr;
            SCIP_VAR* Tz = nullptr;
            SCIP_CALL(SCIPcreateVarBasic(scip, &Tx, "Tx", -SCIP_DEFAULT_INFINITY, SCIP_DEFAULT_INFINITY, 0.0, SCIP_VARTYPE_CONTINUOUS));
            SCIP_CALL(SCIPcreateVarBasic(scip, &Ty, "Ty", -SCIP_DEFAULT_INFINITY, SCIP_DEFAULT_INFINITY, 0.0, SCIP_VARTYPE_CONTINUOUS));
            SCIP_CALL(SCIPcreateVarBasic(scip, &Tz, "Tz", -SCIP_DEFAULT_INFINITY, SCIP_DEFAULT_INFINITY, 0.0, SCIP_VARTYPE_CONTINUOUS));
            SCIP_CALL(SCIPaddVar(scip, Tx));
            SCIP_CALL(SCIPaddVar(scip, Ty));
            SCIP_CALL(SCIPaddVar(scip, Tz));

            // Initialize quaternion variables
            SCIP_VAR* Qw = nullptr;
            SCIP_VAR* Qx = nullptr;
            SCIP_VAR* Qy = nullptr;
            SCIP_VAR* Qz = nullptr;
            SCIP_CALL(SCIPcreateVarBasic(scip, &Qw, "Qw", -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS));
            SCIP_CALL(SCIPcreateVarBasic(scip, &Qx, "Qx", -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS));
            SCIP_CALL(SCIPcreateVarBasic(scip, &Qy, "Qy", -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS));
            SCIP_CALL(SCIPcreateVarBasic(scip, &Qz, "Qz", -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS));
            SCIP_CALL(SCIPaddVar(scip, Qw));
            SCIP_CALL(SCIPaddVar(scip, Qx));
            SCIP_CALL(SCIPaddVar(scip, Qy));
            SCIP_CALL(SCIPaddVar(scip, Qz));

            // Constrain the quaternion to be normalized: Qw * Qw + Qx * Qx + Qy * Qy + Qz * Qz == 1, "|Q| = 1");
            SCIP_CONS* normalizationConstraint = nullptr;
            {
                SCIP_VAR* quad_vars_1[4] = {Qw, Qx, Qy, Qz};
                SCIP_VAR* quad_vars_2[4] = {Qw, Qx, Qy, Qz};
                SCIP_Real quad_coef[4] = {1.0, 1.0, 1.0, 1.0};
                SCIP_CALL(SCIPcreateConsBasicQuadraticNonlinear(scip, &normalizationConstraint, "|Q| = 1", 0, nullptr, nullptr, 4, quad_vars_1, quad_vars_2, quad_coef, 1.0, 1.0)); // lhs = rhs = 1.0
            }
            SCIP_CALL(SCIPaddCons(scip, normalizationConstraint));
            SCIP_CALL( SCIPreleaseCons(scip, &normalizationConstraint));

            // Support variables
            std::array<SCIP_VAR*, 9> R; // Rotation matrix
            for (int i = 0; i < 9; i++) {
                SCIP_CALL(SCIPcreateVarBasic(scip, &R[i], "Ri", -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS));
                SCIP_CALL(SCIPaddVar(scip, R[i]));
            }

            SCIP_CONS* R0Constraint;
            {
                SCIP_VAR* lin_vars[] = {R[0]};
                SCIP_Real lin_coef[] = {1.0};
                SCIP_VAR* quad_vars_1[] = {Qy, Qz};
                SCIP_VAR* quad_vars_2[] = {Qy, Qz};
                SCIP_Real quad_coef[] = {2.0, 2.0};
                SCIP_CALL(SCIPcreateConsBasicQuadraticNonlinear(scip, &R0Constraint, "mat[0]  = 1 - 2 * ( yy + zz )", 1, lin_vars, lin_coef, 2, quad_vars_1, quad_vars_2, quad_coef, 1.0, 1.0));
            }
            SCIP_CALL(SCIPaddCons(scip, R0Constraint));

            SCIP_CONS* R1Constraint;
            {
                SCIP_VAR* lin_vars[] = {R[1]};
                SCIP_Real lin_coef[] = {1.0};
                SCIP_VAR* quad_vars_1[] = {Qx, Qz};
                SCIP_VAR* quad_vars_2[] = {Qy, Qw};
                SCIP_Real quad_coef[] = {-2.0, 2.0};
                SCIP_CALL(SCIPcreateConsBasicQuadraticNonlinear(scip, &R1Constraint, "mat[1]  = 2 * ( xy - zw )", 1, lin_vars, lin_coef, 2, quad_vars_1, quad_vars_2, quad_coef, 0.0, 0.0));
            }
            SCIP_CALL(SCIPaddCons(scip, R1Constraint));

            SCIP_CONS* R2Constraint;
            {
                SCIP_VAR* lin_vars[] = {R[2]};
                SCIP_Real lin_coef[] = {1.0};
                SCIP_VAR* quad_vars_1[] = {Qx, Qy};
                SCIP_VAR* quad_vars_2[] = {Qz, Qw};
                SCIP_Real quad_coef[] = {-2.0, -2.0};
                SCIP_CALL(SCIPcreateConsBasicQuadraticNonlinear(scip, &R2Constraint, "mat[2]  = 2 * ( xz + yw )", 1, lin_vars, lin_coef, 2, quad_vars_1, quad_vars_2, quad_coef, 0.0, 0.0));
            }
            SCIP_CALL(SCIPaddCons(scip, R2Constraint));

            SCIP_CONS* R3Constraint;
            {
                SCIP_VAR* lin_vars[] = {R[3]};
                SCIP_Real lin_coef[] = {1.0};
                SCIP_VAR* quad_vars_1[] = {Qx, Qz};
                SCIP_VAR* quad_vars_2[] = {Qy, Qw};
                SCIP_Real quad_coef[] = {-2.0, -2.0};
                SCIP_CALL(SCIPcreateConsBasicQuadraticNonlinear(scip, &R3Constraint, "mat[3]  = 2 * ( xy + zw )", 1, lin_vars, lin_coef, 2, quad_vars_1, quad_vars_2, quad_coef, 0.0, 0.0));
            }
            SCIP_CALL(SCIPaddCons(scip, R3Constraint));

            SCIP_CONS* R4Constraint;
            {
                SCIP_VAR* lin_vars[] = {R[4]};
                SCIP_Real lin_coef[] = {1.0};
                SCIP_VAR* quad_vars_1[] = {Qx, Qz};
                SCIP_VAR* quad_vars_2[] = {Qx, Qz};
                SCIP_Real quad_coef[] = {2.0, 2.0};
                SCIP_CALL(SCIPcreateConsBasicQuadraticNonlinear(scip, &R4Constraint, "mat[4]  = 1 - 2 * ( xx + zz )", 1, lin_vars, lin_coef, 2, quad_vars_1, quad_vars_2, quad_coef, 1.0, 1.0));
            }
            SCIP_CALL(SCIPaddCons(scip, R4Constraint));

            SCIP_CONS* R5Constraint;
            {
                SCIP_VAR* lin_vars[] = {R[5]};
                SCIP_Real lin_coef[] = {1.0};
                SCIP_VAR* quad_vars_1[] = {Qy, Qx};
                SCIP_VAR* quad_vars_2[] = {Qz, Qw};
                SCIP_Real quad_coef[] = {-2.0, 2.0};
                SCIP_CALL(SCIPcreateConsBasicQuadraticNonlinear(scip, &R5Constraint, "mat[5]  = 2 * ( yz - xw )", 1, lin_vars, lin_coef, 2, quad_vars_1, quad_vars_2, quad_coef, 0.0, 0.0));
            }
            SCIP_CALL(SCIPaddCons(scip, R5Constraint));

            SCIP_CONS* R6Constraint;
            {
                SCIP_VAR* lin_vars[] = {R[6]};
                SCIP_Real lin_coef[] = {1.0};
                SCIP_VAR* quad_vars_1[] = {Qx, Qy};
                SCIP_VAR* quad_vars_2[] = {Qz, Qw};
                SCIP_Real quad_coef[] = {-2.0, 2.0};
                SCIP_CALL(SCIPcreateConsBasicQuadraticNonlinear(scip, &R6Constraint, "mat[6]  = 2 * ( xz - yw )", 1, lin_vars, lin_coef, 2, quad_vars_1, quad_vars_2, quad_coef, 0.0, 0.0));
            }
            SCIP_CALL(SCIPaddCons(scip, R6Constraint));

            SCIP_CONS* R7Constraint;
            {
                SCIP_VAR* lin_vars[] = {R[7]};
                SCIP_Real lin_coef[] = {1.0};
                SCIP_VAR* quad_vars_1[] = {Qy, Qx};
                SCIP_VAR* quad_vars_2[] = {Qz, Qw};
                SCIP_Real quad_coef[] = {-2.0, -2.0};
                SCIP_CALL(SCIPcreateConsBasicQuadraticNonlinear(scip, &R7Constraint, "mat[7]  = 2 * ( yz + xw )", 1, lin_vars, lin_coef, 2, quad_vars_1, quad_vars_2, quad_coef, 0.0, 0.0));
            }
            SCIP_CALL(SCIPaddCons(scip, R7Constraint));

            SCIP_CONS* R8Constraint;
            {
                SCIP_VAR* lin_vars[] = {R[8]};
                SCIP_Real lin_coef[] = {1.0};
                SCIP_VAR* quad_vars_1[] = {Qx, Qy};
                SCIP_VAR* quad_vars_2[] = {Qx, Qy};
                SCIP_Real quad_coef[] = {2.0, 2.0};
                SCIP_CALL(SCIPcreateConsBasicQuadraticNonlinear(scip, &R8Constraint, "mat[8]  = 1 - 2 * ( xx + yy )", 1, lin_vars, lin_coef, 2, quad_vars_1, quad_vars_2, quad_coef, 1.0, 1.0));
            }
            SCIP_CALL(SCIPaddCons(scip, R8Constraint));

            // Symmetry breaking constraints
            auto breakSymmetry = true;
            if(breakSymmetry) {
                // Break container symmetry
                {

                    auto v = itemModelMesh->getVolumeCentroid();
                    //auto v = itemModelMesh->getVertices()[0];

                    for (const auto &symmetryBreakingPlaneSet: SymmetryBreaker::deriveSymmetryBreakingPlanes(containerModelSpaceMesh, 256)){
                        for (const auto &symmetryBreakingPlane: symmetryBreakingPlaneSet){

                            // Force position to be at positive side of the planes
                            auto n = symmetryBreakingPlane.getNormal();
                            auto d = symmetryBreakingPlane.getD();

                            SCIP_CONS* linearPlaneConstraint;
                            {
                                SCIP_VAR* lin_vars[] = {R[0], R[1], R[2],
                                                        R[3], R[4], R[5],
                                                        R[6], R[7], R[8],
                                                        Tx, Ty, Tz,
                                                        inverseS};
                                SCIP_Real lin_coef[] = {n.x * v.x, n.x * v.y, n.x * v.z,
                                                        n.y * v.x, n.y * v.y, n.y * v.z,
                                                        n.z * v.x, n.z * v.y, n.z * v.z,
                                                        n.x, n.y, n.z,
                                                        d};
                                SCIP_CALL(SCIPcreateConsBasicLinear(scip, &linearPlaneConstraint, "Container symmetry breaking constraint", 13, lin_vars, lin_coef, 0, SCIP_DEFAULT_INFINITY));
                            }
                            SCIP_CALL(SCIPaddCons(scip, linearPlaneConstraint));
                        }
                    }
                }

                // Break item symmetry
                {
                    for (const auto &symmetryPlaneSet: SymmetryBreaker::deriveSymmetryBreakingPlanes(itemModelMesh, 256)){

                        if(symmetryPlaneSet.size() == 1){
                            auto symmetryPlane = symmetryPlaneSet[0];


                            std::cout << "Symmetry plane normal: " << symmetryPlane.getNormal().x << ", " << symmetryPlane.getNormal().y << ", " << symmetryPlane.getNormal().z << std::endl;
                            auto n = symmetryPlane.getNormal();
                            auto v = n; // This works for halfspaces, not smaller rotational segments


                            SCIP_CONS* linearPlaneConstraint;
                            {
                                SCIP_VAR* lin_vars[] = {R[0], R[1], R[2],
                                                        R[3], R[4], R[5],
                                                        R[6], R[7], R[8]};
                                SCIP_Real lin_coef[] = {n.x * v.x, n.x * v.y, n.x * v.z,
                                                        n.y * v.x, n.y * v.y, n.y * v.z,
                                                        n.z * v.x, n.z * v.y, n.z * v.z};

                                SCIP_CALL(SCIPcreateConsBasicLinear(scip, &linearPlaneConstraint, "Item symmetry breaking constraint", 12, lin_vars, lin_coef, 0, SCIP_DEFAULT_INFINITY));
                            }
                            SCIP_CALL(SCIPaddCons(scip, linearPlaneConstraint));
                        }
                        else if(symmetryPlaneSet.size() == 2){

                            auto symmetryPlaneA = symmetryPlaneSet[0];
                            auto symmetryPlaneB = symmetryPlaneSet[1];

                            auto nA = symmetryPlaneA.getNormal();
                            auto nB = symmetryPlaneB.getNormal();

                            auto v = nA + nB;

                            SCIP_CONS* linearPlaneConstraintA;
                            {
                                SCIP_VAR* lin_vars[] = {R[0], R[1], R[2],
                                                        R[3], R[4], R[5],
                                                        R[6], R[7], R[8]};
                                SCIP_Real lin_coef[] = {nA.x * v.x, nA.x * v.y, nA.x * v.z,
                                                        nA.y * v.x, nA.y * v.y, nA.y * v.z,
                                                        nA.z * v.x, nA.z * v.y, nA.z * v.z};

                                SCIP_CALL(SCIPcreateConsBasicLinear(scip, &linearPlaneConstraintA, "Item symmetry breaking constraint", 12, lin_vars, lin_coef, 0, SCIP_DEFAULT_INFINITY));
                            }
                            SCIP_CALL(SCIPaddCons(scip, linearPlaneConstraintA));

                            SCIP_CONS* linearPlaneConstraintB;
                            {
                                SCIP_VAR* lin_vars[] = {R[0], R[1], R[2],
                                                        R[3], R[4], R[5],
                                                        R[6], R[7], R[8]};
                                SCIP_Real lin_coef[] = {nB.x * v.x, nB.x * v.y, nB.x * v.z,
                                                        nB.y * v.x, nB.y * v.y, nB.y * v.z,
                                                        nB.z * v.x, nB.z * v.y, nB.z * v.z};

                                SCIP_CALL(SCIPcreateConsBasicLinear(scip, &linearPlaneConstraintB, "Item symmetry breaking constraint", 12, lin_vars, lin_coef, 0, SCIP_DEFAULT_INFINITY));
                            }
                            SCIP_CALL(SCIPaddCons(scip, linearPlaneConstraintB));
                        }
                    }
                }
            }

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

                    SCIP_CONS* linearPlaneConstraint;
                    {
                        SCIP_VAR* lin_vars[] = {R[0], R[1], R[2],
                                                R[3], R[4], R[5],
                                                R[6], R[7], R[8],
                                                Tx, Ty, Tz,
                                                inverseS};
                        SCIP_Real lin_coef[] = {n.x * v.x, n.x * v.y, n.x * v.z,
                                                n.y * v.x, n.y * v.y, n.y * v.z,
                                                n.z * v.x, n.z * v.y, n.z * v.z,
                                                n.x, n.y, n.z,
                                                d};
                        SCIP_CALL(SCIPcreateConsBasicLinear(scip, &linearPlaneConstraint, "Vertex on negative side of facet.", 13, lin_vars, lin_coef, -SCIP_DEFAULT_INFINITY, -0));
                    }
                    SCIP_CALL(SCIPaddCons(scip, linearPlaneConstraint));
                    SCIP_CALL(SCIPreleaseCons(scip, &linearPlaneConstraint));

                }
            }

            // Optimize
            SCIP_CALL(SCIPsolve(scip));

            SCIP_CALL(SCIPprintSol(scip, SCIPgetBestSol(scip), NULL, FALSE) );

            // Extract variable values
            SCIP_Real inverseSValue = SCIPgetSolVal(scip, SCIPgetBestSol(scip), inverseS);
            SCIP_Real TxValue = SCIPgetSolVal(scip, SCIPgetBestSol(scip), Tx);
            SCIP_Real TyValue = SCIPgetSolVal(scip, SCIPgetBestSol(scip), Ty);
            SCIP_Real TzValue = SCIPgetSolVal(scip, SCIPgetBestSol(scip), Tz);
            SCIP_Real QwValue = SCIPgetSolVal(scip, SCIPgetBestSol(scip), Qw);
            SCIP_Real QxValue = SCIPgetSolVal(scip, SCIPgetBestSol(scip), Qx);
            SCIP_Real QyValue = SCIPgetSolVal(scip, SCIPgetBestSol(scip), Qy);
            SCIP_Real QzValue = SCIPgetSolVal(scip, SCIPgetBestSol(scip), Qz);


            Quaternion rotation = Quaternion(glm::fquat(QwValue, QxValue, QyValue, QzValue));

            std::cout << "|Q|: " << glm::length(glm::fquat(rotation)) << std::endl;

            Transformation transformation(rotation);
            transformation.setScale(1.0/inverseSValue);
            transformation.setPositionX(TxValue/inverseSValue);
            transformation.setPositionY(TyValue/inverseSValue);
            transformation.setPositionZ(TzValue/inverseSValue);
            solution->setItemModelTransformation(transformation);

            std::cout << transformation << std::endl;

            notifyObserversSolution(solution);

        } catch(...) {
            std::cout << "Exception during optimization" << std::endl;
        }

        return 0;
    }

    void run() override {
        exec();
    }

} task;

int main(int argc, char *argv[]) {

    QApplication app(argc, argv);
    ApplicationWindow window;
    window.show();
    RenderWidget *renderWidget = window.getRenderWidget();

    std::function<void(RenderWidget* renderWidget, std::shared_ptr<const AbstractSolution> solution)> onSolutionNotified = [](RenderWidget* renderWidget, const std::shared_ptr<const AbstractSolution>& solution){

        renderWidget->clear();

        const auto& meshSolution = dynamic_cast<const AbstractMeshSolution&>(*solution);

        for (int i = 0; i < meshSolution.getItemWorldSpaceMeshes().size(); i++){
            renderWidget->renderWorldSpaceMesh("Items", meshSolution.getItemWorldSpaceMeshes()[i], Color::Red());

            // Render the convex hull as well
            if(meshSolution.getItemWorldSpaceMeshes()[i]->getModelSpaceMesh()->getConvexHull()){
                auto convexhull = std::make_shared<WorldSpaceMesh>(meshSolution.getItemWorldSpaceMeshes()[i]->getModelSpaceMesh()->getConvexHull());
                convexhull->setModelTransformation(meshSolution.getItemWorldSpaceMeshes()[i]->getModelTransformation());
                auto color = Color::Red();
                color.a = 0.5;
                renderWidget->renderWorldSpaceMesh("ConvexHull", convexhull, color);
            }
        }

        for(const auto& worldSpaceMesh: meshSolution.getInclusionWorldSpaceMeshes()){
            renderWidget->renderWorldSpaceMesh("Inclusions", worldSpaceMesh, Color(0,0.5,0.5,1));
        }

        renderWidget->renderWorldSpaceMesh("Container", meshSolution.getOuterWorldSpaceMesh(), Color(1, 1, 1, 0.4));
    };


    renderWidget->observeTask(&task, onSolutionNotified);
    task.start();
    int returnCode = QApplication::exec();
    task.stop();
    task.join();
    return returnCode;
}