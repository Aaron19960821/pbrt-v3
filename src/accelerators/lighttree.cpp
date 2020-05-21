/*************************************************************************
    > File Name: src/accelerators/lighttree.cpp
    > Author: Yuchen Wang
    > Mail: wyc8094@gmail.com 
    > Created Time: Wed May 20 15:32:36 2020
 ************************************************************************/

// accelerators/lighttree.cpp*
#include "accelerators/lighttree.h"

#include "paramset.h"
#include "stats.h"

namespace pbrt {

struct LightInfo {
  LightInfo(size_t _lightIndex, Bounds3f _bound, Float _power):
    lightIndex(_lightIndex),
    bound(_bound),
    power(_power),
    centroid(0.5f*_bound.pMin + 0.5f*_bound.pMax) {}

  size_t lightIndex;
  Bounds3f bound;
  Float power;
  Point3f centroid;
};

struct LightTreeBuildNode {
  Bounds3f bound;
  LightTreeBuildNode* children[2];
  size_t firstLightOffset;
  int splitAxis;
  int nLights;
  Float power;

  void InitLeaf(size_t _offset, int _nLights, 
      Bounds3f _bound, Float _power) {
    firstLightOffset = _offset;
    nLights = _nLights;
    bound = _bound;
    splitAxis = -1;
    children[0] = children[1] = nullptr;
    power = _power;
  }

  void InitInternal(int _splitAxis, LightTreeBuildNode* l, 
      LightTreeBuildNode* r) {
    children[0] = l;
    children[1] = r;
    splitAxis = _splitAxis;
    bound = Union(children[0]->bound, children[1]->bound);
    nLights = children[0]->nLights + children[1]->nLights;
    power = children[0]->power + children[1]->power;
  }

  bool isLeaf() const {
    return splitAxis < 0;
  }
};

struct LinearLightTreeNode {
  int splitAxis;
  Float power;
  union {
    int lightsOffset; // For leaf nodes
    int secondChildOffset; // For internal nodes
  };
  int nLight;

  bool isLeaf() const {
    return splitAxis < 0;
  }
};

LightTree::LightTree(std::vector<std::shared_ptr<Light>> lights,
    int maxLightsPerNode, 
    SplitMethod splitMethod):_lights(std::move(lights)),
                             _maxLightsPerNode(maxLightsPerNode),
                             _splitMethod(splitMethod) {
    ProfilePhase _(Prof::LightTreeConstruction);

    // Check if there are ligths;
    if (_lights.empty()) {
      Warning("LightTree::LightTree(): No light found.");
      return;
    }

    // Initializa light infos
    std::vector<LightInfo> lightsInfo;
    for (size_t i = 0; i < _lights.size(); ++i) {
      lightsInfo.push_back(LightInfo(i, 
            _lights[i]->WorldBound(), 
            _lights[i]->Power().y()));
    }

    std::cout << "all information got." << std::endl; 

    // Build light tree
    MemoryArena memory(1024*1024);
    int totalNodes = 0;
    std::vector<std::shared_ptr<Light>> orderedLights;
    orderedLights.reserve(_lights.size());

    LightTreeBuildNode* root = recursiveBuild(memory, lightsInfo, 0, lightsInfo.size(), 
        &totalNodes, orderedLights);

    std::cout << "build nodes got: " << totalNodes << std::endl; 
    _nodes = AllocAligned<LinearLightTreeNode>(totalNodes);
    int off = 0;
    flattenTree(root, off, orderedLights);
    _lights = orderedLights;

    std::cout << "flattened" << std::endl;
}

LightTree::~LightTree() {
  FreeAligned(_nodes);
}

LightTreeBuildNode* LightTree::recursiveBuild(MemoryArena& arena, std::vector<LightInfo>& lightsInfo,
    int start, int end, int* totalNodes,
    std::vector<std::shared_ptr<Light>>& orderedLights) {
  CHECK_NE(start, end);
  LightTreeBuildNode* node = arena.Alloc<LightTreeBuildNode>();
  (*totalNodes)++;
  
  // Compute Bounds
  Bounds3f bound;
  for (int i = start; i < end; ++i) {
    bound = Union(bound, lightsInfo[i].bound);
  }

  int numLights = end - start;
  if (numLights <= _maxLightsPerNode) {
    // Init a leaf when there is only one leaf.
    size_t offset = orderedLights.size();
    Float power = 0.0f;
    for (int i = start; i < end; ++i) {
      orderedLights.push_back(_lights[lightsInfo[i].lightIndex]);
      power += lightsInfo[i].power;
    }
    node->InitLeaf(offset, numLights, bound, power);
    return node;
  } else {
    Bounds3f centroidBound;
    for (int i = start; i < end; ++i) {
      centroidBound = Union(centroidBound, lightsInfo[i].centroid);
    }
    int dim = centroidBound.MaximumExtent();
    int mid;

    switch (_splitMethod) {
      case SplitMethod::Middle:
        {
          Float lmid = (centroidBound.pMin[dim] + centroidBound.pMax[dim]) / 2.0f;
          LightInfo* midLight = std::partition(&lightsInfo[start], &lightsInfo[end-1]+1, 
              [dim, lmid](const LightInfo& i) {
                return i.centroid[dim] < lmid;
              });
          mid = midLight - &lightsInfo[0];
          // If multiple bounding box overlaps each other, then use equal counts
          if (mid != start && mid != end)
            break;
        }
      case SplitMethod::EqualCounts:
        {
          mid = (start + end) / 2;
          std::nth_element(&lightsInfo[start], &lightsInfo[mid], &lightsInfo[end-1]+1,
              [dim](const LightInfo& a, const LightInfo& b) {
                return a.centroid[dim] < b.centroid[dim];
              });
          break;
        }
      case SplitMethod::SAH:
      default:
        {
        }
    }
    LightTreeBuildNode* l = recursiveBuild(arena, lightsInfo, 
        start, mid, totalNodes, orderedLights);
    LightTreeBuildNode* r = recursiveBuild(arena, lightsInfo, 
        mid, end, totalNodes, orderedLights);
    node->InitInternal(dim, l, r);
  }

  return node;
}

int LightTree::flattenTree(LightTreeBuildNode* node, 
    int& off, std::vector<std::shared_ptr<Light> >& orderedLights) {
  LinearLightTreeNode* linearNode = &_nodes[off];
  linearNode->power = node->power;
  int curOffset = off++;
  if (node->isLeaf()) {
    linearNode->lightsOffset = node->firstLightOffset;
    linearNode->splitAxis = -1;
    linearNode->nLight = node->nLights;
  } else {
    linearNode->splitAxis = node->splitAxis;
    linearNode->nLight = node->nLights;
    flattenTree(node->children[0], off, orderedLights);
    linearNode->secondChildOffset = flattenTree(node->children[1], off, orderedLights);
  }

  return curOffset;
}

void LightTree::sample(int curOffset, Float u, int& offset, Float& pdf) const {
  LinearLightTreeNode curNode = _nodes[curOffset];
  if (curNode.isLeaf()) {
    int off = static_cast<int>(curNode.nLight * u);
    int selectedIndex = curNode.lightsOffset + off;
    Light* light = _lights[selectedIndex].get();
    offset = selectedIndex;
    pdf = light->Power().y() / curNode.power;
  } else {
    Float pl = _nodes[curOffset+1].power;
    Float p = pl / curNode.power;
    if (u < p) {
      sample(curOffset+1, u/p, offset, pdf);
      pdf = p * pdf;
    } else {
      sample(curNode.secondChildOffset, (1.0f-u)/(1.0f-p), offset, pdf);
      pdf = (1.0-p) * pdf;
    }
  }
}

std::shared_ptr<Light> LightTree::getLightByIndex(int index) const {
  CHECK_NE(index, _lights.size());
  return _lights[index];
}

std::shared_ptr<LightTree> CreateLightTree(
    std::vector<std::shared_ptr<Light>> lights,
    const ParamSet& paramset) {
  LightTree::SplitMethod splitMethod;
  std::string splitMethodName = paramset.FindOneString("splitmethod", "sah");
  if (splitMethodName == "sah") {
    splitMethod = LightTree::SplitMethod::SAH;
  } else if (splitMethodName == "middle") {
    splitMethod = LightTree::SplitMethod::Middle;
  } else if (splitMethodName == "equal") {
    splitMethod = LightTree::SplitMethod::EqualCounts;
  } else {
    Warning("LightTree split method unknown. Using sah instead.");
    splitMethod = LightTree::SplitMethod::SAH;
  }

  int maxLightsPerNode = paramset.FindOneInt("maxlightspernode", 4);
  return std::make_shared<LightTree>(std::move(lights), maxLightsPerNode, splitMethod);
}

} // namespace pbrt
