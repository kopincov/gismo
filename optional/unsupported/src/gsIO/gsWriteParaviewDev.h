

#define PLOT_PRECISION 10


namespace gismo
{

/// Visualizing an edge graph of a 3D solid structured as gsHeMesh
//template <class T>
//void gsWriteParaview(gsHeMesh<T> const& sl, std::string const & fn);


/// Visualizing an edge graph of a 3D solid structured as gsHeMesh
template <class T>
void gsWriteParaview(gsHeMesh<T> const& sl, std::string const & fn)
{
    gsHeMesh<T> mesh = sl;
    std::string mfn(fn);
    mfn.append(".vtp");
    std::ofstream file(mfn.c_str());
    if ( ! file.is_open() )
        std::cout<<"Problem opening "<<fn<<std::endl;
    file << std::fixed; // no exponents
    file << std::setprecision (PLOT_PRECISION);

    file <<"<?xml version=\"1.0\"?>\n";
    file <<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file <<"<PolyData>\n";

    /// Number of vertices and number of faces
    file <<"<Piece NumberOfPoints=\""<< mesh.getnumVertices() <<"\" NumberOfVerts=\"0\" NumberOfLines=\""
         << mesh.getnumHalfEdges()<<"\" NumberOfStrips=\"0\" NumberOfPolys=\""<< mesh.getnumHalfFaces() << "\">\n";

    /// Coordinates of vertices
    file <<"<Points>\n";
    file <<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    std::vector< gsHeVertex<T>* > ve = mesh.getVertices();
    for (typename std::vector< gsHeVertex<T>* >::const_iterator it=ve.begin(); it!=ve.end(); ++it)
    {
        const gsHeVertex<T>& vertex = **it;
        file << vertex[0] << " " << vertex[1] << " " << vertex[2] << " \n";
    }

    file << "\n";
    file <<"</DataArray>\n";
    file <<"</Points>\n";

    // Scalar field attached to each face
    // file << "<PointData Scalars=\"point_scalars\">\n";
    // file << "<DataArray type=\"Int32\" Name=\"point_scalars\" format=\"ascii\">\n";
    // for (typename std::vector< gsVertex<T>* >::const_iterator it=sl.vertex.begin();
    //      it!=sl.vertex.end(); ++it)
    // {
    //     file << 0 << " ";
    // }
    // file << "\n";
    // file << "</DataArray>\n";
    // file << "</PointData>\n";

    // Write out edges
    file << "<Lines>\n";
    file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    typename std::vector< gsHalfEdge<T>* > edge = mesh.getHalfEdges();
    for (typename std::vector< gsHalfEdge<T>* >::const_iterator it=edge.begin();
         it!=edge.end(); ++it)
    {
            file << (*it)->getSource()->getId() << " " << (*it)->getNext()->getSource()->getId() << "\n";
    }
    file << "</DataArray>\n";
    file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int count=0;
    edge = mesh.getHalfEdges();
    for (typename std::vector< gsHalfEdge<T>* >::const_iterator it=edge.begin();
         it!=edge.end(); ++it)
    {
        count+=2;
        file << count << " ";
    }
    file << "\n";
    file << "</DataArray>\n";
    file << "</Lines>\n";

    // Scalar field attached to each face (* if edges exists, this has a problem)
    // file << "<CellData Scalars=\"cell_scalars\">\n";
    // file << "<DataArray type=\"Int32\" Name=\"cell_scalars\" format=\"ascii\">\n";
    // for (typename std::vector< gsFace<T>* >::const_iterator it=sl.face.begin();
    //      it!=sl.face.end(); ++it)
    // {
    //     file << 1 << " ";
    // }
    // file << "\n";
    // file << "</DataArray>\n";
    // file << "</CellData>\n";

    /// Which vertices belong to which faces
    file << "<Polys>\n";
    file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    std::vector< gsHalfFace<T>* > face = mesh.getHalfFaces();
    for (typename std::vector< gsHalfFace<T>* >::const_iterator it=face.begin();
         it!=face.end(); ++it)
    {
        std::vector< gsHeVertex<T>* > ve2 = (*it)->getVertices();
        for (typename std::vector< gsHeVertex<T>* >::const_iterator vit= ve2.begin();
             vit!=ve2.end(); ++vit)
        {
            file << (*vit)->getId() << " ";
        }
        file << "\n";
    }
    file << "</DataArray>\n";
    file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    count=0;
    face = mesh.getHalfFaces();
    for (typename std::vector< gsHalfFace<T>* >::const_iterator it=face.begin();
         it!=face.end(); ++it)
    {
        count += (*it)->getVertices().size();
        file << count << " ";
    }
    file << "\n";
    file << "</DataArray>\n";
    file << "</Polys>\n";

    file << "</Piece>\n";
    file <<"</PolyData>\n";
    file <<"</VTKFile>\n";
    file.close();

    makeCollection(fn, ".vtp");

}




}// namespace gismo


#undef PLOT_PRECISION
