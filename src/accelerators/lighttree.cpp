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
    power(_power) {}

  size_t lightIndex;
  Bounds3f bound;
  Float power;
};

struct LightTreeBuildNode {
  Bounds3f bound;
  LightTreeBuildNode* children[2];
  size_t firstLightOffset;
  int splitAxis;
  int nLights;

  void InitLeaf(size_t _offset, int _nLights, Bounds3f _bound) {
    firstLightOffset = _offset;
    nLights = _nLights;
    bound = _bound;
    splitAxis = -1;
    children[0] = children[1] = nullptr;
  }

  bool isLeaf() {
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

    // Build light tree
    MemoryArena memory(1024*1024);
    int totalNodes = 0;
    std::vector<std::shared_ptr<Light>> orderedLights;
    orderedLights.reserve(_lights.size());
}

LightTree::~LightTree() {
}

LightTreeBuildNode* recursiveBuild(MemoryArena& arena, std::vector<LightInfo>& lightsInfo,
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
  if (numLights == 1) {
    // Init a leaf when there is only one leaf.
    size_t offset = orderedLights.size();
    size_t index = lightsInfo[start].lightIndex;
    node->InitLeaf(offset, numLights, bound);
    return node;
  }

  return node;
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
