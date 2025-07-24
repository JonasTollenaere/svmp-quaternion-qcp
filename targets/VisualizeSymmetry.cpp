//
// Created by Jonas on 15/05/2024.
//

#include "QApplication"
#include "meshcore/rendering/ApplicationWindow.h"
#include "meshcore/utility/FileParser.h"
#include "SymmetryBreaker.h"
#include "meshcore/utility/random.h"
#include "meshcore/factories/SphereFactory.h"

void run(RenderWidget* renderWidget);

int main(int argc, char *argv[]){

    QApplication app(argc, argv);
    ApplicationWindow window;
    window.setFixedWidth(1300);
    window.setFixedHeight(760);
    window.show();

    auto renderWidget = window.getRenderWidget();
    auto openGLWidget = renderWidget->getOpenGLWidget();
    openGLWidget->setLightMode(true);
    openGLWidget->setFixedWidth(720);
    openGLWidget->setFixedHeight(720);

    std::thread workerThread(run, window.getRenderWidget());

    int returnCode = QApplication::exec();
    workerThread.join();
    return returnCode;
}

void run(RenderWidget* renderWidget){

    {
        // Load a mesh
        auto filePath = MESHCORE_DATA_DIR + std::string("Tollenaere, J. et al/Items/arrow.obj");
//    auto filePath = MESHCORE_DATA_DIR + "Tollenaere, J. et al/Items/star.stl";
//    auto filePath = MESHCORE_DATA_DIR + "Tollenaere, J. et al/Items/bobbin.stl";
        auto modelSpaceMesh = FileParser::loadMeshFile(filePath);
        auto worldSpaceMesh = std::make_shared<WorldSpaceMesh>(modelSpaceMesh);

//    Quaternion frontViewRotation;
//    Quaternion frontViewRotation = Quaternion(glm::half_pi<float>(), 0, 0);
        Quaternion frontViewRotation = Quaternion(0, glm::half_pi<float>(), glm::half_pi<float>());
        Transformation preprocessTransform(frontViewRotation);
        preprocessTransform.setPosition(frontViewRotation.rotateVertex(-modelSpaceMesh->getVolumeCentroid()));
        worldSpaceMesh->setModelTransformation(preprocessTransform);
        modelSpaceMesh = worldSpaceMesh->getTransformedModelSpaceMesh();
        worldSpaceMesh = std::make_shared<WorldSpaceMesh>(modelSpaceMesh);

        renderWidget->renderWorldSpaceMesh("Group-1", worldSpaceMesh, Color(0.8, 0.8, 0.8, 1.0));

        auto symmetryBreakingPlanes = SymmetryBreaker::deriveSymmetryBreakingPlanes(modelSpaceMesh, 64);

        auto modelAABB = modelSpaceMesh->getBounds();
        auto size = glm::length(modelAABB.getHalf());
        auto centroid = modelSpaceMesh->getVolumeCentroid();

        // Visualize the symmetry breaking constraints
        Random random(1);
        for (int i = 0; i < symmetryBreakingPlanes.size(); ++i){

            if(i!=1) continue; // Only render the one used in the paper

            auto color = Color::Red();

            for (int j = 0; j < symmetryBreakingPlanes[i].size(); ++j){
                auto &plane = symmetryBreakingPlanes[i][j];
                std::cout << "Plane " << j << " normal: " << plane.getNormal().x << ", " << plane.getNormal().y << ", " << plane.getNormal().z << std::endl;
                std::cout << "Plane " << j << " d: " << plane.getD() << std::endl;
                //renderWidget->renderPlane("Group-" + std::to_string(i), "Plane-" + std::to_string(j), plane, color);

                renderWidget->renderRay("Group-" + std::to_string(i), "Ray-" + std::to_string(j), Ray(centroid + glm::vec3(5,0,0), plane.getNormal() * 5.0f), color);

                // AABB aligned with this plane

                glm::vec3 min;
                if(symmetryBreakingPlanes[i].size() == 1){
                    min = glm::vec3(-size, -size, -size/100); // Reflectional
                }
                else {
                    min = glm::vec3(-size, -size/200, -size/100); // Rotational symmetry
                }
                glm::vec3 max(size, size, size/100);
                AABB aabb(min, max);

                auto axis = glm::cross(plane.getNormal(), glm::vec3(0,0,1));
                float angle = -glm::acos(glm::dot(plane.getNormal(), glm::vec3(0,0,1)));

                if(glm::length(axis) < 1e-6){
                    axis = glm::vec3(1,0,0);
                }

                Quaternion normalAlignmentRotation(axis, angle);

                Transformation transform;
                transform.setRotation(normalAlignmentRotation);
                transform.setPosition(centroid);

                {
                    std::vector<Vertex> vertices;
                    std::vector<IndexTriangle> triangles;

                    vertices.emplace_back(min.x, min.y, min.z);
                    vertices.emplace_back(max.x, min.y, min.z);
                    vertices.emplace_back(max.x, max.y, min.z);
                    vertices.emplace_back(min.x, max.y, min.z);
                    vertices.emplace_back(min.x, min.y, max.z);
                    vertices.emplace_back(max.x, min.y, max.z);
                    vertices.emplace_back(max.x, max.y, max.z);
                    vertices.emplace_back(min.x, max.y, max.z);

                    // Add the triangles
                    unsigned int numberOfVertices = vertices.size();
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 7, numberOfVertices - 8, numberOfVertices - 6});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 6, numberOfVertices - 8, numberOfVertices - 5});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 1, numberOfVertices - 4, numberOfVertices - 2});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 2, numberOfVertices - 4, numberOfVertices - 3});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 5, numberOfVertices - 8, numberOfVertices - 4});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 1, numberOfVertices - 5, numberOfVertices - 4});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 2, numberOfVertices - 6, numberOfVertices - 5});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 2, numberOfVertices - 5, numberOfVertices - 1});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 4, numberOfVertices - 7, numberOfVertices - 3});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 4, numberOfVertices - 8, numberOfVertices - 7});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 3, numberOfVertices - 7, numberOfVertices - 6});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 3, numberOfVertices - 6, numberOfVertices - 2});

                    auto mesh = std::make_shared<ModelSpaceMesh>(vertices, triangles);
                    mesh->setName("Box-" + std::to_string(j));

                    auto worldMesh = std::make_shared<WorldSpaceMesh>(mesh);

                    worldMesh->setModelTransformation(transform);

                    renderWidget->renderWorldSpaceMesh("Group-" + std::to_string(i), worldMesh, color);

                }



                AABB orthogonalAABB(glm::vec3(-size/100, -size, 0), glm::vec3(size/100, size, size));
//                renderWidget->renderBox("Group-" + std::to_string(i), "Ortho-" + std::to_string(j), orthogonalAABB, transform, color);

                {
                    std::vector<Vertex> vertices;
                    std::vector<IndexTriangle> triangles;

                    vertices.emplace_back(orthogonalAABB.getMinimum().x, orthogonalAABB.getMinimum().y, orthogonalAABB.getMinimum().z);
                    vertices.emplace_back(orthogonalAABB.getMaximum().x, orthogonalAABB.getMinimum().y, orthogonalAABB.getMinimum().z);
                    vertices.emplace_back(orthogonalAABB.getMaximum().x, orthogonalAABB.getMaximum().y, orthogonalAABB.getMinimum().z);
                    vertices.emplace_back(orthogonalAABB.getMinimum().x, orthogonalAABB.getMaximum().y, orthogonalAABB.getMinimum().z);
                    vertices.emplace_back(orthogonalAABB.getMinimum().x, orthogonalAABB.getMinimum().y, orthogonalAABB.getMaximum().z);
                    vertices.emplace_back(orthogonalAABB.getMaximum().x, orthogonalAABB.getMinimum().y, orthogonalAABB.getMaximum().z);
                    vertices.emplace_back(orthogonalAABB.getMaximum().x, orthogonalAABB.getMaximum().y, orthogonalAABB.getMaximum().z);
                    vertices.emplace_back(orthogonalAABB.getMinimum().x, orthogonalAABB.getMaximum().y, orthogonalAABB.getMaximum().z);

                    // Add the triangles
                    unsigned int numberOfVertices = vertices.size();
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 7, numberOfVertices - 8, numberOfVertices - 6});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 6, numberOfVertices - 8, numberOfVertices - 5});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 1, numberOfVertices - 4, numberOfVertices - 2});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 2, numberOfVertices - 4, numberOfVertices - 3});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 5, numberOfVertices - 8, numberOfVertices - 4});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 1, numberOfVertices - 5, numberOfVertices - 4});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 2, numberOfVertices - 6, numberOfVertices - 5});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 2, numberOfVertices - 5, numberOfVertices - 1});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 4, numberOfVertices - 7, numberOfVertices - 3});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 4, numberOfVertices - 8, numberOfVertices - 7});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 3, numberOfVertices - 7, numberOfVertices - 6});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 3, numberOfVertices - 6, numberOfVertices - 2});

                    auto mesh = std::make_shared<ModelSpaceMesh>(vertices, triangles);
                    mesh->setName("Back-" + std::to_string(j));

                    auto worldMesh = std::make_shared<WorldSpaceMesh>(mesh);

                    worldMesh->setModelTransformation(transform);

                    renderWidget->renderWorldSpaceMesh("Group-" + std::to_string(i), worldMesh, Color(85.f/255.f, 170.f/255.f, 0, 155.f/255.f));
                }
            }
        }
    }

    {
        auto filePath = MESHCORE_DATA_DIR + std::string("Tollenaere, J. et al/Items/frustum.stl");
        auto modelSpaceMesh = FileParser::loadMeshFile(filePath);
        auto worldSpaceMesh = std::make_shared<WorldSpaceMesh>(modelSpaceMesh);

        Quaternion frontViewRotation = Quaternion(glm::half_pi<float>(), 0, 0);
        Transformation preprocessTransform(frontViewRotation);
        preprocessTransform.setPosition(frontViewRotation.rotateVertex(-modelSpaceMesh->getVolumeCentroid()));
        worldSpaceMesh->setModelTransformation(preprocessTransform);
        modelSpaceMesh = worldSpaceMesh->getTransformedModelSpaceMesh();
        worldSpaceMesh = std::make_shared<WorldSpaceMesh>(modelSpaceMesh);

        renderWidget->renderWorldSpaceMesh("DGroup-0", worldSpaceMesh, Color(0.8, 0.8, 0.8, 1.0));

        auto symmetryBreakingPlanes = SymmetryBreaker::deriveSymmetryBreakingPlanes(modelSpaceMesh, 64);

        auto modelAABB = modelSpaceMesh->getBounds();
        auto size = glm::length(modelAABB.getHalf());
        auto centroid = modelSpaceMesh->getVolumeCentroid();

        // Visualize the symmetry breaking constraints
        for (int i = 0; i < symmetryBreakingPlanes.size(); ++i){

            auto color = Color::Red();

            std::vector<Vertex> outerpoints;
            for (int j = 0; j < symmetryBreakingPlanes[i].size(); ++j){
                auto &plane = symmetryBreakingPlanes[i][j];
                std::cout << "Plane " << j << " normal: " << plane.getNormal().x << ", " << plane.getNormal().y << ", " << plane.getNormal().z << std::endl;
                std::cout << "Plane " << j << " d: " << plane.getD() << std::endl;
                //renderWidget->renderPlane("Group-" + std::to_string(i), "Plane-" + std::to_string(j), plane, color);

                renderWidget->renderRay("DGroup-" + std::to_string(i), "Ray-" + std::to_string(j), Ray(centroid + glm::vec3(5,0,0), plane.getNormal() * 1.f), color);

                // AABB aligned with this plane

                glm::vec3 min;
                if(symmetryBreakingPlanes[i].size() == 1){
                    min = glm::vec3(-size, -size, -size/100); // Reflectional
                }
                else {
                    min = glm::vec3(-size, -size/200, -size/100); // Rotational symmetry
                }
                glm::vec3 max(size, size, size/100);
                AABB aabb(min, max);

                auto axis = glm::cross(plane.getNormal(), glm::vec3(0,0,1));
                float angle = -glm::acos(glm::dot(plane.getNormal(), glm::vec3(0,0,1)));

                if(j!=0) angle += glm::pi<float>();

                if(glm::length(axis) < 1e-6){
                    axis = glm::vec3(1,0,0);
                }

                Quaternion normalAlignmentRotation(axis, angle);

                outerpoints.emplace_back(Transformation(normalAlignmentRotation).transformVertex(glm::vec3(0,size,0)) + centroid);

                Ray ray(aabb.getCenter() + glm::vec3(2,0,0), glm::vec3(0,0,j!=0?-1:1));
                ray = ray.getTransformed(Transformation(normalAlignmentRotation));
                renderWidget->renderRay("DGroup-" + std::to_string(i), "Ray-" + std::to_string(j), ray, color);

                Transformation transform;
                transform.setRotation(normalAlignmentRotation);
                transform.setPosition(centroid);

                {
                    std::vector<Vertex> vertices;
                    std::vector<IndexTriangle> triangles;

                    vertices.emplace_back(min.x, min.y, min.z);
                    vertices.emplace_back(max.x, min.y, min.z);
                    vertices.emplace_back(max.x, max.y, min.z);
                    vertices.emplace_back(min.x, max.y, min.z);
                    vertices.emplace_back(min.x, min.y, max.z);
                    vertices.emplace_back(max.x, min.y, max.z);
                    vertices.emplace_back(max.x, max.y, max.z);
                    vertices.emplace_back(min.x, max.y, max.z);

                    // Add the triangles
                    unsigned int numberOfVertices = vertices.size();
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 7, numberOfVertices - 8, numberOfVertices - 6});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 6, numberOfVertices - 8, numberOfVertices - 5});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 1, numberOfVertices - 4, numberOfVertices - 2});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 2, numberOfVertices - 4, numberOfVertices - 3});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 5, numberOfVertices - 8, numberOfVertices - 4});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 1, numberOfVertices - 5, numberOfVertices - 4});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 2, numberOfVertices - 6, numberOfVertices - 5});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 2, numberOfVertices - 5, numberOfVertices - 1});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 4, numberOfVertices - 7, numberOfVertices - 3});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 4, numberOfVertices - 8, numberOfVertices - 7});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 3, numberOfVertices - 7, numberOfVertices - 6});
                    triangles.emplace_back(IndexTriangle{numberOfVertices - 3, numberOfVertices - 6, numberOfVertices - 2});

                    auto mesh = std::make_shared<ModelSpaceMesh>(vertices, triangles);
                    mesh->setName("DBox-" + std::to_string(j));

                    auto worldMesh = std::make_shared<WorldSpaceMesh>(mesh);

                    worldMesh->setModelTransformation(transform);

                    renderWidget->renderWorldSpaceMesh("DGroup-" + std::to_string(i), worldMesh, color);
                }
            }

            outerpoints.emplace_back(centroid);

            if (outerpoints.size() >= 3) {
                std::vector<IndexTriangle> triangles = {{0,1,2}};
                auto mesh =  std::make_shared<ModelSpaceMesh>(outerpoints, triangles);
                renderWidget->renderWorldSpaceMesh("DGroup-" + std::to_string(i), std::make_shared<WorldSpaceMesh>(mesh), Color(85.f/255.f, 170.f/255.f, 0, 155.f/255.f));
            }
        }
    }
}
