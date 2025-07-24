//
// Created by Jonas on 19/04/2024.
//

#include "gurobi_c++.h"

#include <meshcore/tasks/AbstractTask.h>
#include <meshcore/utility/FileParser.h>
#include <meshcore/utility/io.h>
#include <meshcore/utility/random.h>
#include <meshcore/rendering/ApplicationWindow.h>
#include <meshcore/optimization/SingleVolumeMaximisationSolution.h>

#include <nlohmann/json.hpp>
#include <iostream>
#include <string>
#include <glm/gtx/component_wise.hpp>

#include "UpperBoundProvider.h"
#include "ConvexConcavitiesFactory.h"
#include "SymmetryBreaker.h"

class Task: public AbstractTask {

    class Callback: public GRBCallback {
    public:

        std::shared_ptr<SingleVolumeMaximisationSolution> solution;

        GRBVar* inverseS;
        GRBVar* Tx;
        GRBVar* Ty;
        GRBVar* Tz;
        GRBVar* Qw;
        GRBVar* Qx;
        GRBVar* Qy;
        GRBVar* Qz;

        Task* task;

        Callback(Task* task, const std::shared_ptr<SingleVolumeMaximisationSolution>& callbackSolution, GRBVar* inverseS, GRBVar* Tx, GRBVar* Ty, GRBVar* Tz, GRBVar* Qw, GRBVar* Qx, GRBVar* Qy, GRBVar* Qz) : inverseS(inverseS), Tx(Tx), Ty(Ty), Tz(Tz), Qw(Qw), Qx(Qx), Qy(Qy), Qz(Qz), solution(callbackSolution), task(task) {}

    protected:
        void callback() override {
            if (where == GRB_CB_MIPSOL) {

                Quaternion rotation = Quaternion(glm::fquat(getSolution(*Qw), getSolution(*Qx), getSolution(*Qy), getSolution(*Qz)));
                Transformation transformation(rotation);
                transformation.setPositionX(getSolution(*Tx)/getSolution(*inverseS));
                transformation.setPositionY(getSolution(*Ty)/getSolution(*inverseS));
                transformation.setPositionZ(getSolution(*Tz)/getSolution(*inverseS));
                transformation.setScale(1.0/getSolution(*inverseS));
                solution->getItemWorldSpaceMesh()->setModelTransformation(transformation);

                task->notifyObserversSolution(solution);
            }
        }
    };

    void run() override {

        Random random;

        this->notifyObserversStatus("Initialising");

//        auto containerPath = MESHCORE_DATA_DIR + std::string("Tollenaere, J. et al/Items/arrow.obj");
//        auto containerPath = MESHCORE_DATA_DIR + std::string("Tollenaere, J. et al/Items/ring.obj");
//        auto containerPath = MESHCORE_DATA_DIR + std::string("Tollenaere, J. et al/Items/star.obj");
//        auto containerPath = MESHCORE_DATA_DIR + std::string("Tollenaere, J. et al/Items/star.stl");
        auto containerPath = MESHCORE_DATA_DIR + std::string("Tollenaere, J. et al/Items/torus.stl");
//        auto containerPath = MESHCORE_DATA_DIR + std::string("Tollenaere, J. et al/Items/bobbin.stl");
//        auto containerPath = MESHCORE_DATA_DIR + std::string("Tollenaere, J. et al/Items/stone_1.obj");
//        auto containerPath = MESHCORE_DATA_DIR + std::string("Tollenaere, J. et al/Items/tree.stl");

//        auto itemPath = MESHCORE_DATA_DIR + std::string("Tollenaere, J. et al/Items/convexstone.stl");
//        auto itemPath = MESHCORE_DATA_DIR + std::string("Tollenaere, J. et al/Items/cube.obj");
//        auto itemPath = MESHCORE_DATA_DIR + std::string("Tollenaere, J. et al/Items/dodecahedron.stl");
        auto itemPath = MESHCORE_DATA_DIR + std::string("Tollenaere, J. et al/Items/frustum.stl");
//        auto itemPath = MESHCORE_DATA_DIR + std::string("Tollenaere, J. et al/Items/gem.stl");
//        auto itemPath = MESHCORE_DATA_DIR + std::string("Tollenaere, J. et al/Items/helmet.stl");
//        auto itemPath = MESHCORE_DATA_DIR + std::string("Tollenaere, J. et al/Items/oloid.stl");
//        auto itemPath = MESHCORE_DATA_DIR + std::string("Tollenaere, J. et al/Items/tetrahedron.stl");

        auto containerModelSpaceMesh = FileParser::loadMeshFile(containerPath);
        auto itemModelMesh = FileParser::loadMeshFile(itemPath);
        auto containerMesh = std::make_shared<WorldSpaceMesh>(containerModelSpaceMesh);
        auto itemMesh = std::make_shared<WorldSpaceMesh>(itemModelMesh);
        auto solution = std::make_shared<SingleVolumeMaximisationSolution>(itemMesh, containerMesh);
        solution->getItemWorldSpaceMesh()->getModelTransformation().setScale(0.0f);
        notifyObserversSolution(solution);

        std::cout << "Processing: " << itemModelMesh->getName() << " in " << containerModelSpaceMesh->getName() << std::endl;

        auto upperScaleBound = UpperBoundProvider::getUpperScaleBoundConvexHullItem(itemModelMesh, containerModelSpaceMesh);
        double timeInFirstModel = 0.0;
        bool convexRelaxationBound = true;
        if(convexRelaxationBound){

            auto convexContainerHull = containerModelSpaceMesh->getConvexHull();
            try {
                // Create an environment
                GRBEnv env = GRBEnv(true);
                env.set("LogFile", "quat.log");
                env.start();

                // Create an empty model
                GRBModel model = GRBModel(env);
                model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
                model.set(GRB_IntParam_Threads, 64);
                model.set(GRB_IntParam_NonConvex, 2);

                // Initialize decision variables
                auto relaxedScaleBound = UpperBoundProvider::getUpperScaleBound(itemModelMesh, convexContainerHull);
                GRBVar inverseS = model.addVar(1.0/relaxedScaleBound, 1000, 1.0, GRB_CONTINUOUS, "S");

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
                {
                    auto v = itemModelMesh->getVolumeCentroid();

                    for (const auto &symmetryBreakingPlaneSet: SymmetryBreaker::deriveSymmetryBreakingPlanes(containerModelSpaceMesh, 256)){
                        for (const auto &symmetryBreakingPlane: symmetryBreakingPlaneSet){

                            // Force position to be at positive side of the planes
                            auto n = symmetryBreakingPlane.getNormal();
                            auto d = symmetryBreakingPlane.getD();
                            model.addConstr(n.x * (R[0] * v.x + R[1] * v.y + R[2] * v.z + Tx) +
                                            n.y * (R[3] * v.x + R[4] * v.y + R[5] * v.z + Ty) +
                                            n.z * (R[6] * v.x + R[7] * v.y + R[8] * v.z + Tz) + inverseS*d >= 0, "Container symmetry breaking constraint");
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
                        model.addConstr(n.x * (R[0] * v.x + R[1] * v.y + R[2] * v.z + Tx) +
                                        n.y * (R[3] * v.x + R[4] * v.y + R[5] * v.z + Ty) +
                                        n.z * (R[6] * v.x + R[7] * v.y + R[8] * v.z + Tz) + inverseS*d <= -1e-8, "Vertex on negative side of facet.");
                    }
                }

                // Optimize model
                notifyObserversStatus("Solving for upper bound");
                model.optimize();

                // Extract the achieved scale and update the upper bound
                float achievedScale = 1.0f/inverseS.get(GRB_DoubleAttr_X);
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
            model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
            model.set(GRB_IntParam_Threads, 64);
            model.set(GRB_IntParam_NonConvex, 2);

            // Limit execution time to 3600s
            model.set(GRB_DoubleParam_TimeLimit, 3600 - timeInFirstModel);

            // Initialize decision variables
            GRBVar inverseS = model.addVar(1.0/upperScaleBound, 1000, 1.0, GRB_CONTINUOUS, "S");

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

            assert(itemModelMesh->getConvexHull()->isConvex());
            assert(containerModelSpaceMesh->getConvexHull()->isConvex());

            // Symmetry breaking constraints
            {

                auto v = itemModelMesh->getVolumeCentroid();

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
                    model.addConstr(n.x * (R[0] * v.x + R[1] * v.y + R[2] * v.z + Tx) +
                                    n.y * (R[3] * v.x + R[4] * v.y + R[5] * v.z + Ty) +
                                    n.z * (R[6] * v.x + R[7] * v.y + R[8] * v.z + Tz) + inverseS * d <= -1e-8, "Vertices should be on the negative side of a facet that is part of the convex hull.");
                }
            }

            std::vector<std::shared_ptr<ModelSpaceMesh>> concavities = ConvexConcavitiesFactory::getConvexConcavities(containerModelSpaceMesh);
            std::vector<std::vector<Plane>> groupedConcavityPlanes = ConvexConcavitiesFactory::getGroupedConcavityPlanes(containerModelSpaceMesh);
            if(!concavities.empty()){

                for (int g = 0; g < concavities.size(); ++g){

                    const auto& concavity = concavities[g];
                    bool limitSeparationPlanes = false;

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
                    auto precompute = false;
                    if(precompute){
                        // Precompute products of (Nx, Ny, Nz) and transformation matrix in separate variables,
                        // they are identical for all vertices and will allow linear plane equations
                        // This also allows us to manually put sensible bounds on them
                        auto NxR0 = model.addVar(-upperScaleBound, upperScaleBound, 0, GRB_CONTINUOUS, "NxSR0-" + std::to_string(g));
                        auto NxR1 = model.addVar(-upperScaleBound, upperScaleBound, 0, GRB_CONTINUOUS, "NxSR1-" + std::to_string(g));
                        auto NxR2 = model.addVar(-upperScaleBound, upperScaleBound, 0, GRB_CONTINUOUS, "NxSR2-" + std::to_string(g));
                        auto NyR3 = model.addVar(-upperScaleBound, upperScaleBound, 0, GRB_CONTINUOUS, "NySR3-" + std::to_string(g));
                        auto NyR4 = model.addVar(-upperScaleBound, upperScaleBound, 0, GRB_CONTINUOUS, "NySR4-" + std::to_string(g));
                        auto NyR5 = model.addVar(-upperScaleBound, upperScaleBound, 0, GRB_CONTINUOUS, "NySR5-" + std::to_string(g));
                        auto NzR6 = model.addVar(-upperScaleBound, upperScaleBound, 0, GRB_CONTINUOUS, "NzSR6-" + std::to_string(g));
                        auto NzR7 = model.addVar(-upperScaleBound, upperScaleBound, 0, GRB_CONTINUOUS, "NzSR7-" + std::to_string(g));
                        auto NzR8 = model.addVar(-upperScaleBound, upperScaleBound, 0, GRB_CONTINUOUS, "NzSR8-" + std::to_string(g));
                        auto NxTx = model.addVar(-std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), 0, GRB_CONTINUOUS, "NxTx-" + std::to_string(g));
                        auto NyTy = model.addVar(-std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), 0, GRB_CONTINUOUS, "NyTy-" + std::to_string(g));
                        auto NzTz = model.addVar(-std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), 0, GRB_CONTINUOUS, "NzTz-" + std::to_string(g));
                        auto inverseSD = model.addVar(-std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), 0, GRB_CONTINUOUS, "inverseSD-" + std::to_string(g));

                        model.addQConstr(NxR0 == R[0] * Nx, "NxR0-" + std::to_string(g));
                        model.addQConstr(NxR1 == R[1] * Nx, "NxR1-" + std::to_string(g));
                        model.addQConstr(NxR2 == R[2] * Nx, "NxR2-" + std::to_string(g));
                        model.addQConstr(NxTx  == Tx    * Nx, "NxTx-"  + std::to_string(g));
                        model.addQConstr(NyR3 == R[3] * Ny, "NyR3-" + std::to_string(g));
                        model.addQConstr(NyR4 == R[4] * Ny, "NyR4-" + std::to_string(g));
                        model.addQConstr(NyR5 == R[5] * Ny, "NyR5-" + std::to_string(g));
                        model.addQConstr(NyTy  == Ty    * Ny, "NyTy-"  + std::to_string(g));
                        model.addQConstr(NzR6 == R[6] * Nz, "NzR6-" + std::to_string(g));
                        model.addQConstr(NzR7 == R[7] * Nz, "NzR7-" + std::to_string(g));
                        model.addQConstr(NzR8 == R[8] * Nz, "NzR8-" + std::to_string(g));
                        model.addQConstr(NzTz  == Tz    * Nz, "NzTz-"  + std::to_string(g));
                        model.addQConstr(inverseSD == inverseS * D, "inverseSD-" + std::to_string(g));

                        // cube in torus:

                        for (const auto &v: itemModelMesh->getConvexHull()->getVertices()){
                            model.addConstr(NxR0 * v.x + NxR1 * v.y + NxR2 * v.z + NxTx +
                                            NyR3 * v.x + NyR4 * v.y + NyR5 * v.z + NyTy +
                                            NzR6 * v.x + NzR7 * v.y + NzR8 * v.z + NzTz + inverseSD >= 1e-8, "Vertices should be on the positive side of a separating plane between the concavity and the convex item.");
                        }
                    }
                    else{

                        // cube in torus:

                         for (const auto &v: itemModelMesh->getConvexHull()->getVertices()){
                            model.addQConstr(Nx * (R[0] * v.x + R[1] * v.y + R[2] * v.z + Tx) +
                                             Ny * (R[3] * v.x + R[4] * v.y + R[5] * v.z + Ty) +
                                             Nz * (R[6] * v.x + R[7] * v.y + R[8] * v.z + Tz) + inverseS * D >= 1e-8, "Vertices should be on the positive side of a separating plane between the concavity and the convex item.");
                         }
                    }



                }
            }

            // Configure callback
            Callback cb(this, solution, &inverseS, &Tx, &Ty, &Tz, &Qw, &Qx, &Qy, &Qz);
            model.setCallback(&cb);

            // Optimize model
            notifyObserversStatus("Optimizing");
            model.optimize();

            // Print all variables and their values
            for (int i = 0; i < model.get(GRB_IntAttr_NumVars); i++) {
                std::cout << model.getVar(i).get(GRB_StringAttr_VarName) << " = " << model.getVar(i).get(GRB_DoubleAttr_X) << std::endl;
            }

            Quaternion rotation = Quaternion(glm::fquat(Qw.get(GRB_DoubleAttr_X), Qx.get(GRB_DoubleAttr_X), Qy.get(GRB_DoubleAttr_X), Qz.get(GRB_DoubleAttr_X)));
            Transformation transformation(rotation);
            transformation.setScale(1.0/inverseS.get(GRB_DoubleAttr_X));
            transformation.setPositionX(Tx.get(GRB_DoubleAttr_X)/inverseS.get(GRB_DoubleAttr_X));
            transformation.setPositionY(Ty.get(GRB_DoubleAttr_X)/inverseS.get(GRB_DoubleAttr_X));
            transformation.setPositionZ(Tz.get(GRB_DoubleAttr_X)/inverseS.get(GRB_DoubleAttr_X));
            solution->getItemWorldSpaceMesh()->setModelTransformation(transformation);

            std::cout << transformation << std::endl;

            notifyObserversSolution(solution);

            auto totalRunTime = model.get(GRB_DoubleAttr_Runtime) + timeInFirstModel;
            std::cout << "Total Gurobi runtime: " << totalRunTime << "s" << std::endl;

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

