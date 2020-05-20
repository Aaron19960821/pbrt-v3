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

LightTree::LightTree(std::vector<std::shared_ptr<Light>> lights,
    int maxLightsPerNode, 
    SplitMethod splitMethod):_light(std::move(lights)),
                             _maxLightsPerNode(maxLightsPerNode),
                             _splitMethod(splitMethod) {
    ProfilePhase _(Prof::LightTreeConstruction);

    // Check if there are ligths;
    if (_light.empty()) {
      Warning("LightTree::LightTree(): No light found.");
      return;
    }
}

LightTree::~LightTree() {
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
