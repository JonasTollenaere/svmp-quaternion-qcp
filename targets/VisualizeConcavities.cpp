//
// Created by Jonas on 15/05/2024.
//

#include "QApplication"
#include "meshcore/rendering/ApplicationWindow.h"
#include "meshcore/utility/FileParser.h"
#include "ConvexConcavitiesFactory.h"
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

    // Load a mesh
//    auto filePath = MESHCORE_DATA_DIR + "Tollenaere, J. et al/Items/arrow.obj";
//    auto filePath = MESHCORE_DATA_DIR + "Tollenaere, J. et al/Items/bobbin.stl";
//    auto filePath = MESHCORE_DATA_DIR + "Tollenaere, J. et al/Items/dagger.stl";
    auto filePath = MESHCORE_DATA_DIR + std::string("Tollenaere, J. et al/Items/star.obj");
//    auto filePath = MESHCORE_DATA_DIR + "Tollenaere, J. et al/Items/star.stl";
//    auto filePath = MESHCORE_DATA_DIR + "Tollenaere, J. et al/Items/ring.obj";
//    auto filePath = MESHCORE_DATA_DIR + "Tollenaere, J. et al/Items/bobbin.stl";
//    auto filePath = MESHCORE_DATA_DIR + "Tollenaere, J. et al/Items/bobbin.stl";
//    auto filePath = MESHCORE_DATA_DIR + "Tollenaere, J. et al/Items/torus.stl";
//    auto filePath = MESHCORE_DATA_DIR + "Tollenaere, J. et al/Items/tree.stl";
//    auto filePath = MESHCORE_DATA_DIR + "Tollenaere, J. et al/Items/stone_1.obj";
    auto modelSpaceMesh = FileParser::loadMeshFile(filePath);

    if(!modelSpaceMesh){
        return;
    }

    auto worldSpaceMesh = std::make_shared<WorldSpaceMesh>(modelSpaceMesh);

//    Quaternion frontViewRotation;
//    Quaternion frontViewRotation = Quaternion(glm::half_pi<float>(), 0, 0);
//    Quaternion frontViewRotation = Quaternion(0, glm::half_pi<float>(), glm::half_pi<float>());
//    Transformation preprocessTransform(frontViewRotation);
//    preprocessTransform.setPosition(frontViewRotation.rotateVertex(-modelSpaceMesh->getVolumeCentroid()));
//    worldSpaceMesh->setModelTransformation(preprocessTransform);
    modelSpaceMesh = worldSpaceMesh->getTransformedModelSpaceMesh();
    worldSpaceMesh = std::make_shared<WorldSpaceMesh>(modelSpaceMesh);

    renderWidget->renderWorldSpaceMesh("Meshes", worldSpaceMesh, Color(0.8, 0.8, 0.8, 1));

    auto convexConcavities = ConvexConcavitiesFactory::getConvexConcavities(modelSpaceMesh);

    auto concaveContainerColor = Color(0.7, 0.7, 0.65, 1);
    auto transparentContainerColor = Color(0.7, 0.7, 0.65, .35);
    renderWidget->renderWorldSpaceMesh("Hulls", std::make_shared<WorldSpaceMesh>(modelSpaceMesh->getConvexHull()), transparentContainerColor);
    {
        Random random;
        for (const auto &concavity: convexConcavities) {
            auto concavityMesh = std::make_shared<WorldSpaceMesh>(concavity);
            concavityMesh->setModelTransformation(worldSpaceMesh->getModelTransformation());
            random.nextFloat();
            const Color color(random.nextFloat(), random.nextFloat(), random.nextFloat(), 1.0);
            renderWidget->renderWorldSpaceMesh("Concavities", concavityMesh, color);
        }
    }

    {
        Random random(0);
        const Color color(random.nextFloat(), random.nextFloat(), random.nextFloat(), 1.0);
        for (const auto &concavity: convexConcavities) {
            auto concavityMesh = std::make_shared<WorldSpaceMesh>(concavity);
            concavityMesh->setModelTransformation(worldSpaceMesh->getModelTransformation());
            random.nextFloat();
            renderWidget->renderWorldSpaceMesh("ConcavitiesFixedColor", concavityMesh, color);
        }
    }

}
