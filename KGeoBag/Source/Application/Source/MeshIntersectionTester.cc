#include "KGCoreMessage.hh"
#include "KGInterfaceBuilder.hh"
#include "KGMesher.hh"
#include "KGVTKMeshIntersectionTester.hh"
#include "KGVTKMeshPainter.hh"
#include "KXMLInitializer.hh"
#include "KXMLTokenizer.hh"

using namespace KGeoBag;
using namespace katrin;
using namespace std;

int main(int argc, char** argv)
{
    if (argc < 3) {
        cout << "usage: ./MeshViewer <config_file_name.xml> <geometry_path> [...]" << endl;
        return -1;
    }

    coremsg(eNormal) << "starting initialization..." << eom;

    auto& tXML = KXMLInitializer::GetInstance();
    tXML.AddDefaultIncludePath(CONFIG_DEFAULT_DIR);
    tXML.Configure(argc, argv);

    deque<string> tPathList = tXML.GetArguments().ParameterList();
    tPathList.pop_front();  // strip off config file name

    coremsg(eNormal) << "...initialization finished" << eom;

    KVTKWindow tWindow;
    tWindow.SetName("KGeoBag Mesh Viewer");
    tWindow.SetFrameColorRed(0.);
    tWindow.SetFrameColorGreen(0.);
    tWindow.SetFrameColorBlue(0.);
    tWindow.SetDisplayMode(true);
    tWindow.SetWriteMode(true);

    KGMesher tMesher;

    KGVTKMeshPainter tPainter;
    tPainter.SetName("MeshPainter");
    tPainter.SetDisplayMode(true);
    tPainter.SetWriteMode(true);
    tPainter.SetColorMode(KGVTKMeshPainter::sModulo);

    KGVTKMeshIntersectionTester tTester;
    tTester.SetName("MeshIntersectionTester");
    tTester.SetDisplayMode(true);
    tTester.SetWriteMode(true);
    tTester.SetSampleCount(500);
    tTester.SetSampleColor(KGRGBColor(127, 127, 127));
    tTester.SetUnintersectedLineColor(KGRGBColor(127, 127, 127));
    tTester.SetIntersectedLineColor(KGRGBColor(0, 0, 255));
    tTester.SetVertexSize(0.001);
    tTester.SetLineSize(0.001);

    for (auto& tPath : tPathList) {
        for (auto& tSurface : KGInterface::GetInstance()->RetrieveSurfaces(tPath)) {
            tSurface->MakeExtension<KGMesh>();
            tSurface->AcceptNode(&tMesher);
            tSurface->AcceptNode(&tPainter);
            tTester.AddSurface(tSurface);
        }
        for (auto& tSpace : KGInterface::GetInstance()->RetrieveSpaces(tPath)) {
            tSpace->MakeExtension<KGMesh>();
            tSpace->AcceptNode(&tMesher);
            tSpace->AcceptNode(&tPainter);
            tTester.AddSpace(tSpace);
        }
    }

    tWindow.AddPainter(&tPainter);
    tWindow.AddPainter(&tTester);
    tWindow.Render();
    tWindow.Write();
    tWindow.Display();
    tWindow.RemovePainter(&tPainter);
    tWindow.RemovePainter(&tTester);

    return 0;
}
