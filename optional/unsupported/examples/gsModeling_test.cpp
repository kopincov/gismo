#include <iostream>

#include <gismo.h>
#include <gismo_dev.h>

#include <fstream>


using namespace gismo;


int main(int argc, char *argv[])
{
    bool plot = 0; // If user gives --plot as argiment, paraview file is generated and launched on exit
    bool toxml =0;
    std::string filename("off/mushroom_triangulated.off");
    gsCmdLine cmd("Recover the features of a triangulated surface.");
    cmd.addPlainString("filename", "File containing the input mesh", filename);
    cmd.addSwitch("xml", "Output solid to xml file", toxml);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    std::string baseName = gsFileManager::getBasename(filename);

    gsInfo << "Processing " << baseName << "\n";

    gsMesh<>::uPtr m = gsReadFile<>(filename);
    if (m)
      gsInfo<< "Got "<< *m <<"\n";
    else
    {
      gsInfo<< "Problem encountered in file "<<filename<<", quitting." <<"\n";
      return 0;
    }
    std::vector< gsMesh<> *> paraMeshes;
    std::vector< gsMesh<> *> fitMeshes;
    std::vector< gsMesh<> *> patchMeshes;


    std::vector<std::vector<gsVertex<>*> > iPoints;
    std::vector<std::vector<gsVertex<>*> > oPoints;
    std::vector< std::vector<std::vector<gsVertex<>*> > >  innerBdrys;
    std::vector< std::vector<gsVertex<> > > innerBdrysMassP;
    std::vector<std::vector<bool> > oPointsConvexFlag;
    gsInfo<<"finding patches..."<<"\n";
    gsTriMeshToSolid<> tmts(m.get());
    tmts.getPatchData(12,12,0.2,2,iPoints,oPoints,innerBdrys,innerBdrysMassP,oPointsConvexFlag,"unused_file_name",0);
    gsInfo<<"patches found"<<"\n";
    gsSolid<>::uPtr sl(new gsSolid<>());
    gsInfo<<"fitting..."<<"\n";
    tmts.toSolid(*sl,iPoints,oPoints,innerBdrys,innerBdrysMassP,oPointsConvexFlag,paraMeshes,fitMeshes,patchMeshes,4,12,1,500,1,100,1,1,1);
    gsInfo<<"fitting done"<<"\n";

    gsInfo<<*sl<<'\n';
    if (toxml)
    {
        gsInfo << "Writing xml file..." << "\n";

        gsFileData<> newdata;
        newdata << *sl;
        newdata.dump(baseName);
    }

    //write sharpness information to a file:
    // first line: #patches
    // next #patches lines: #outer points of patch i, sharpness of first point ... sharpness of last point of patch i
    // next #patches lines: #holes of patch i
    // next #holes of patch i lines: #points of hole j of patch i, sharpness of first point ... sharpness of last point of hole j of patch i
    std::string outputFilename = baseName + "_corner_data.txt";
    std::ofstream myfile (outputFilename.c_str());
    if (myfile.is_open())
    {
      myfile << oPoints.size()<<"\n";

      for(size_t i=0;i<oPoints.size();i++)
      {
          myfile << oPoints[i].size()<<" ";
          for(size_t j=0;j<oPoints[i].size();j++)
             myfile<<(oPoints[i][j]->numEdges>2)<<" ";
          myfile << "\n";

      }
      myfile << "\n";

      for(size_t i=0;i<innerBdrys.size();i++)
      {
          myfile << innerBdrys[i].size()<<" ";
          for(size_t j=0;j<innerBdrys[i].size();j++)
          {
              myfile << innerBdrys[i][j].size()<<" ";
              for(size_t k=0;k<innerBdrys[i][j].size();k++)
                  myfile<<(innerBdrys[i][j][k]->numEdges>2)<<" ";
          }
          myfile << "\n";
      }


      myfile.close();
    }

    int exitCommand = 0; // Command to execute on program exit
    if (plot)
    {
        // Write a paraview file
        gsInfo<<"Writing paraview file..." << "\n";
        for(size_t i=0;i<fitMeshes.size();i++)
        {
            std::stringstream fitss;
            fitss << baseName<<"_fit"<<i;
            std::string fitstr = fitss.str();
            const char * fitc = fitstr.c_str();
            gsWriteParaview( *fitMeshes[i], fitc);
            std::stringstream parass;
            parass << baseName<<"_para"<<i;
            std::string parastr = parass.str();
            const char * parac = parastr.c_str();
            gsWriteParaview( *paraMeshes[i], parac);
            std::stringstream patchss;
            patchss << baseName<<"_patch"<<i;
            std::string patchstr = patchss.str();
            const char * patchc = patchstr.c_str();
            gsWriteParaview( *patchMeshes[i], patchc);
        }
        gsWriteParaview( *m, "output");

        exitCommand = system("paraview output.vtp &");

    }
    // free meshes
    freeAll(fitMeshes);
    freeAll(paraMeshes);
    freeAll(patchMeshes);

    return exitCommand;
}
