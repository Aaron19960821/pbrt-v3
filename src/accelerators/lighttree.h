/*************************************************************************
    > File Name: src/accelerators/lighttree.h
    > Author: Yuchen Wang
    > Mail: wyc8094@gmail.com 
    > Created Time: Wed May 20 15:14:36 2020
 ************************************************************************/

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_ACCELERATORS_LIGHTTREE_H_
#define PBRT_ACCELERATORS_LIGHTTREE_H_

#include "light.h"
#include "pbrt.h"

#include <memory>

namespace pbrt {

struct LightInfo;
struct LightTreeBuildNode;
struct LinearLightTreeNode;

class LightCone {
  public:
    LightCone();
    LightCone(const Vector3f&, Float, Float);

    Float measure() const;

    Vector3f axis() const {
      return _axis;
    }

    Float thetaO() const {
      return _thetaO;
    }

    Float thetaE() const {
      return _thetaE;
    }

    static LightCone Union(const LightCone& lc1, const LightCone& lc2);

  private:
    Vector3f _axis;
    Float _thetaO;
    Float _thetaE;
};

class LightTree {
  public:
    enum class SplitMethod {
      Middle,
      EqualCounts,
      SAH
    };

  public:
    LightTree(std::vector<std::shared_ptr<Light>> lights,
        int maxLightsPerNode = 1,
        SplitMethod splitMethod = SplitMethod::SAH);
    ~LightTree();

    void sample(int curOffset, Float u, int& lightIndex, Float& pdf) const;
    std::shared_ptr<Light> getLightByIndex(int index) const;

  private:
    LightTreeBuildNode* recursiveBuild(MemoryArena& arena, std::vector<LightInfo>&,
        int start, int end, int* totalNodes, 
        std::vector<std::shared_ptr<Light>>&);

    int flattenTree(LightTreeBuildNode* node, int& offset, 
        std::vector<std::shared_ptr<Light>>&);

  private:
    std::vector<std::shared_ptr<Light>> _lights;
    const int _maxLightsPerNode;
    const SplitMethod _splitMethod;
    LinearLightTreeNode* _nodes;
};

std::shared_ptr<LightTree> CreateLightTree(
    std::vector<std::shared_ptr<Light>> lights,
    const ParamSet& params
    );

} // namespace pbrt

#endif // PBRT_ACCELERATORS_LIGHTTREE_H_
