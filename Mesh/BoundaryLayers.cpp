// Gmsh - Copyright (C) 1997-2010 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <gmsh@geuz.org>.

#include "GModel.h"
#include "MLine.h"
#include "MTriangle.h"
#include "MQuadrangle.h"
#include "BoundaryLayers.h"
#include "ExtrudeParams.h"
#include "meshGEdge.h"
#include "meshGFace.h"
#include "GmshMessage.h"

#if defined(HAVE_POST)
#include "PView.h"
#include "OctreePost.h"
#else
class OctreePost{ int dummy; };
#endif

template<class T>
static void addExtrudeNormals(std::vector<T*> &elements, int invert, 
                              OctreePost *octree, int index)
{
  // FIXME: generalize this
  if(index < 0 || index > 1){
    Msg::Error("Boundary layer index should be 0 or 1");
    return;
  }

  if(octree){ // get extrusion direction from post-processing view
    std::set<MVertex*> verts;
    for(unsigned int i = 0; i < elements.size(); i++)
      for(int j = 0; j < elements[i]->getNumVertices(); j++)
        verts.insert(elements[i]->getVertex(j));
    for(std::set<MVertex*>::iterator it = verts.begin(); it != verts.end(); it++){
      MVertex *v = *it;
      double nn[3];
#if defined(HAVE_POST)
      octree->searchVector(v->x(), v->y(), v->z(), nn, 0);
#endif
      ExtrudeParams::normals[index]->add(v->x(), v->y(), v->z(), 3, nn);
    }
  }
  else{ // get extrusion data from Gouraud-shaded element normals
    for(unsigned int i = 0; i < elements.size(); i++){
      MElement *ele = elements[i];
      SVector3 n(0, 0, 0);
      if(ele->getDim() == 2)
        n = ele->getFace(0).normal();
      else if(ele->getDim() == 1) // FIXME only valid in XY-plane
        n = crossprod(ele->getEdge(0).tangent(), SVector3(0, 0, 1));
      if(invert) n *= -1.;
      if(n[0] || n[1] || n[2]){
        double nn[3] = {n[0], n[1], n[2]};
        for(int k = 0; k < ele->getNumVertices(); k++){
          MVertex *v = ele->getVertex(k);
          ExtrudeParams::normals[index]->add(v->x(), v->y(), v->z(), 3, nn);
        }
      }
    }
  }
}

typedef std::set<std::pair<bool, int> > infoset;

template<class T>
static void addExtrudeNormals(std::set<T*> &entities, 
                              std::map<int, infoset> &infos,
                              std::map<int, int> &views)
{
  for(typename std::set<T*>::iterator it = entities.begin(); it != entities.end(); it++){
    T *ge = *it;
    int view = views[ge->tag()];
    infoset info = infos[ge->tag()];
    for(infoset::iterator it2 = info.begin(); it2 != info.end(); it2++){
      bool invert = it2->first;
      int index = it2->second;
      OctreePost *octree = 0;
#if defined(HAVE_POST)
      if(view >= 0){
        if(view >= 0 && view < PView::list.size())
          octree = new OctreePost(PView::list[view]);
        else
          Msg::Error("Unknown View[%d]: using normals instead", view);
      }
#endif
      if(ge->dim() == 1)
        addExtrudeNormals(((GEdge*)ge)->lines, invert, octree, index);
      else if(ge->dim() == 2){
        addExtrudeNormals(((GFace*)ge)->triangles, invert, octree, index);
        addExtrudeNormals(((GFace*)ge)->quadrangles, invert, octree, index);
      }
    }
  }
}

int Mesh2DWithBoundaryLayers(GModel *m)
{
  std::set<GFace*> sourceFaces, otherFaces;
  std::set<GEdge*> sourceEdges, otherEdges;
  std::map<int, infoset> sourceFaceInfo, sourceEdgeInfo;
  std::map<int, int> sourceFaceView, sourceEdgeView;

  // 2D boundary layers
  for(GModel::eiter it = m->firstEdge(); it != m->lastEdge(); it++){
    GEdge *ge = *it;
    if(ge->getNativeType() == GEntity::GmshModel && 
       ge->geomType() == GEntity::BoundaryLayerCurve){
      ExtrudeParams *ep = ge->meshAttributes.extrude;
      if(ep && ep->mesh.ExtrudeMesh && ep->geo.Mode == COPIED_ENTITY){
        GEdge *from = m->getEdgeByTag(std::abs(ep->geo.Source));
        if(!from){
          Msg::Error("Unknown source curve %d for boundary layer", ep->geo.Source);
          return 0;
        }
        std::pair<bool, int> tags(ep->geo.Source < 0, ep->mesh.BoundaryLayerIndex);
        sourceEdgeInfo[from->tag()].insert(tags);
        sourceEdgeView[from->tag()] = ep->mesh.ViewIndex; 
        sourceEdges.insert(from);
      }
    }
  }

  // 3D boundary layers
  for(GModel::fiter it = m->firstFace(); it != m->lastFace(); it++){
    GFace *gf = *it;
    if(gf->getNativeType() == GEntity::GmshModel && 
       gf->geomType() == GEntity::BoundaryLayerSurface){
      ExtrudeParams *ep = gf->meshAttributes.extrude;
      if(ep && ep->mesh.ExtrudeMesh && ep->geo.Mode == COPIED_ENTITY){
        GFace *from = m->getFaceByTag(std::abs(ep->geo.Source));
        if(!from){
          Msg::Error("Unknown source face %d for boundary layer", ep->geo.Source);
          return 0;
        }
        std::pair<bool, int> tags(ep->geo.Source < 0, ep->mesh.BoundaryLayerIndex);
        sourceFaceInfo[from->tag()].insert(tags);
        sourceFaceView[from->tag()] = ep->mesh.ViewIndex; 
        sourceFaces.insert(from);
        std::list<GEdge*> e = from->edges();
        sourceEdges.insert(e.begin(), e.end());
      }
    }
  }

  if(sourceEdges.empty() && sourceFaces.empty())
    return 0;

  // compute set of non-source edges and faces
  for(GModel::eiter it = m->firstEdge(); it != m->lastEdge(); it++)
    if(sourceEdges.find(*it) == sourceEdges.end()) otherEdges.insert(*it);
  for(GModel::fiter it = m->firstFace(); it != m->lastFace(); it++)
    if(sourceFaces.find(*it) == sourceFaces.end()) otherFaces.insert(*it);

  // mesh source surfaces
  std::for_each(sourceFaces.begin(), sourceFaces.end(), meshGFace());  
  
  // make sure the source surfaces for the boundary layers are
  // oriented correctly (normally we do this only after the 3D mesh is
  // done; but here it's critical since we use the normals for the
  // extrusion)
  std::for_each(sourceFaces.begin(), sourceFaces.end(), orientMeshGFace());

  // compute a normal field for all the source edges or faces
  for(int i = 0; i < 2; i++){
    if(ExtrudeParams::normals[i]) delete ExtrudeParams::normals[i];
    ExtrudeParams::normals[i] = new smooth_data();
  }
  if(sourceFaces.empty())
    addExtrudeNormals(sourceEdges, sourceEdgeInfo, sourceEdgeView);
  else
    addExtrudeNormals(sourceFaces, sourceFaceInfo, sourceFaceView);
  if(sourceEdgeView.empty() && sourceFaceView.empty())
    for(int i = 0; i < 2; i++)
      ExtrudeParams::normals[i]->normalize();

  // set the position of boundary layer points using the smooth normal
  // field 
  for(GModel::eiter it = m->firstEdge(); it != m->lastEdge(); it++){
    GEdge *ge = *it;
    if(ge->geomType() == GEntity::BoundaryLayerCurve){
      ExtrudeParams *ep = ge->meshAttributes.extrude;
      if(ep && ep->mesh.ExtrudeMesh && ep->geo.Mode == EXTRUDED_ENTITY){
        GVertex *vsrc, *vdest;
        if(ge->getBeginVertex()->geomType() == GEntity::BoundaryLayerPoint){
          vsrc = ge->getEndVertex();
          vdest = ge->getBeginVertex();
        }
        else{
          vsrc = ge->getBeginVertex();
          vdest = ge->getEndVertex();
        }
        GPoint p = vsrc->point();
        ep->Extrude(ep->mesh.NbLayer - 1, ep->mesh.NbElmLayer[ep->mesh.NbLayer - 1],
                    p.x(), p.y(), p.z());
        vdest->setPosition(p);
      }
    }
  }
  
  // remesh non-source edges (since they might have been modified by
  // the change in boundary layer points)
  std::for_each(otherFaces.begin(), otherFaces.end(), deMeshGFace());
  std::for_each(otherEdges.begin(), otherEdges.end(), meshGEdge());

  // mesh the curves bounding the boundary layers by extrusion using
  // the smooth normal field
  for(GModel::eiter it = m->firstEdge(); it != m->lastEdge(); it++){
    GEdge *ge = *it;
    if(ge->geomType() == GEntity::BoundaryLayerCurve){
      Msg::Info("Meshing curve %d", ge->tag());
      deMeshGEdge dem;
      dem(ge);
      MeshExtrudedCurve(ge);
    }
  }

  // mesh non-source surfaces
  std::for_each(otherFaces.begin(), otherFaces.end(), meshGFace());  

  // mesh the surfaces bounding the boundary layers by extrusion using
  // the smooth normal field
  for(GModel::fiter it = m->firstFace(); it != m->lastFace(); it++){
    GFace *gf = *it;
    if(gf->geomType() == GEntity::BoundaryLayerSurface){
      Msg::Info("Meshing surface %d (%s)", gf->tag(), gf->getTypeString().c_str());
      deMeshGFace dem;
      dem(gf);
      MeshExtrudedSurface(gf);
    }
  }

  return 1;
}
