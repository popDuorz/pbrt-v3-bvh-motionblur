
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// accelerators/kdtreeaccel.cpp*
#include "accelerators/kdtreeaccel.h"
#include "paramset.h"
#include "interaction.h"
#include "stats.h"
#include <algorithm>


namespace pbrt {

// KdTreeAccel Local Declarations
struct KdAccelNode {

    // KdAccelNode Methods
    void InitLeaf(int *primNums, int np, std::vector<int> *primitiveIndices);

	//modified by popduorz
	void InitInterior(int axis, int ac, Float s) {
		split = s;
		flags = axis;
		aboveChild |= (ac << 2);
	}

    int nPrimitives() const { return nPrims >> 2; }
    int SplitAxis() const { return flags & 3; }
    bool IsLeaf() const { return (flags & 3) == 3; }
	int AboveChild() const { return aboveChild >> 2; }

	//add by popduorz
    int tempPrimNum() const { return temp_prim_num; }
	void setTempPrimNum(int tpn = 0) { temp_prim_num = tpn; }

    union {
        Float split;                 // Interior
        int onePrimitive;            // Leaf
        int primitiveIndicesOffset;  // Leaf
    };

	//add by popduorz
	int* temp_prims_ptr;
	int temp_prim_num;
	Float split0;
	Float split1;

  private:
    union {
        int flags;       // Both
        int nPrims;      // Leaf
        int aboveChild;  // Interior
    };
	
};

enum class EdgeType { Start, End };

struct BoundEdge {
    // BoundEdge Public Methods
    BoundEdge() {}
    BoundEdge(Float t, int primNum, bool starting) : t(t), primNum(primNum) {
        type = starting ? EdgeType::Start : EdgeType::End;
    }
    Float t;
    int primNum;
    EdgeType type;
};

inline void Interpolate(Bounds3f &boundst, const Bounds3f &bounds0, const Bounds3f &bounds1, const Float &t) {

	boundst.pMin = bounds0.pMin*(1.0f - t) + bounds1.pMin * t;
	boundst.pMax = bounds0.pMax*(1.0f - t) + bounds1.pMax * t;
}

// KdTreeAccel Method Definitions
KdTreeAccel::KdTreeAccel(const std::vector<std::shared_ptr<Primitive>> &p,
                         int isectCost, int traversalCost, Float emptyBonus,
                         int maxPrims, int maxDepth)
    : isectCost(isectCost),
      traversalCost(traversalCost),
      maxPrims(maxPrims),
      emptyBonus(emptyBonus),
      primitives(p) {
    // Build kd-tree for accelerator
    ProfilePhase _(Prof::AccelConstruction);
    nextFreeNode = nAllocedNodes = 0;
    if (maxDepth <= 0)
        maxDepth = std::round(8 + 1.3f * Log2Int(int64_t(primitives.size())));

    // Compute bounds for kd-tree construction
    std::vector<Bounds3f> primBounds;
    primBounds.reserve(primitives.size());
    for (const std::shared_ptr<Primitive> &prim : primitives) {
        Bounds3f b = prim->WorldBound();
        bounds = Union(bounds, b);
        primBounds.push_back(b);

		//add by popduorz -0.5v +0.5v???
		Bounds3f b0 = prim->WorldBound0();
		bounds0 = Union(bounds0, b0);

		Bounds3f b1 = prim->WorldBound1();
		bounds1 = Union(bounds1, b1);
    }

	//add by popduorz
	bounds = Union(bounds0, bounds);
	bounds = Union(bounds1, bounds);


    // Allocate working memory for kd-tree construction
	std::unique_ptr<BoundEdge[]> edges[3];
	for (int i = 0; i < 3; ++i)
		edges[i].reset(new BoundEdge[2 * primitives.size()]);

    std::unique_ptr<int[]> prims0(new int[primitives.size()]);
    std::unique_ptr<int[]> prims1(new int[(maxDepth + 1) * primitives.size()]);

    // Initialize _primNums_ for kd-tree construction
    std::unique_ptr<int[]> primNums(new int[primitives.size()]);
    for (size_t i = 0; i < primitives.size(); ++i) primNums[i] = i;

    // Start recursive construction of kd-tree
    buildTree(0, bounds, primBounds, primNums.get(), primitives.size(),
              maxDepth, prims0.get(), prims1.get(), false);
}


void KdAccelNode::InitLeaf(int *primNums, int np,
                           std::vector<int> *primitiveIndices) {
    flags = 3;
    nPrims |= (np << 2);
    // Store primitive ids for leaf1
    if (np == 0)
        onePrimitive = 0;
    else if (np == 1)
        onePrimitive = primNums[0];
    else {
        primitiveIndicesOffset = primitiveIndices->size();
        for (int i = 0; i < np; ++i) 
			primitiveIndices->push_back(primNums[i]);
    }
}

KdTreeAccel::~KdTreeAccel() { FreeAligned(nodes); }


inline float tri_min(float pos_0, float pos, float pos_1) {

	if (pos_0 <= pos && pos_0 <= pos_1) {
		return pos_0;
	}
	else if (pos <= pos_0&&pos <= pos_1) {
		return pos;
	}
	else {
		return pos_1;
	}
}

inline float tri_max(float pos_0, float pos, float pos_1) {

	if (pos_0 >= pos && pos_0 >= pos_1) {
		return pos_0;
	}
	else if (pos >= pos_0&&pos >= pos_1) {
		return pos;
	}
	else {
		return pos_1;
	}
}

void KdTreeAccel::buildTree(const int &nodeNum, const Bounds3f &nodeBounds,
                            const std::vector<Bounds3f> &allPrimBounds,
                            int *primNums,const int &nPrimitives,const int &depth,
                            int *prims0, int *prims1, int badRefines) {

    CHECK_EQ(nodeNum, nextFreeNode);
    // Get next free node from _nodes_ array
	if (nextFreeNode + 1 == nAllocedNodes || nextFreeNode + 1 > nAllocedNodes) {
		int nNewAllocNodes = std::max(2 * nAllocedNodes, 512);
		KdAccelNode *n = AllocAligned<KdAccelNode>(nNewAllocNodes);
		if (nAllocedNodes > 0) {
			memcpy(n, nodes, nAllocedNodes * sizeof(KdAccelNode));
			FreeAligned(nodes);
		}
		nodes = n;

		for (int i = nAllocedNodes; i < nNewAllocNodes; i++) {
			nodes[i].setTempPrimNum();
			nodes[i].temp_prims_ptr = NULL;
		}

		nAllocedNodes = nNewAllocNodes;
	}
	++nextFreeNode;


    // Initialize leaf node if termination criteria met
    if (nPrimitives <= maxPrims || depth == 0) {
        nodes[nodeNum].InitLeaf(primNums, nPrimitives, &primitiveIndices);
        return;
    }

    // Initialize interior node and continue recursion

    // Choose split axis position for interior node
    int bestAxis = -1, bestOffset = -1;
    Float bestCost = Infinity;
    Float oldCost = isectCost * Float(nPrimitives);
    Float totalSA = nodeBounds.SurfaceArea();
    Float invTotalSA = 1 / totalSA;
    Vector3f d = nodeBounds.pMax - nodeBounds.pMin;

    // Choose which axis to split along
    int axis = nodeBounds.MaximumExtent();
    int retries = 0;

	std::unique_ptr<BoundEdge[]> edges[3];
	for (int i = 0; i < 3; ++i)
		edges[i].reset(new BoundEdge[2 * nPrimitives]);

retrySplit:

    // Initialize edges for _axis_
    for (int i = 0; i < nPrimitives; ++i) {
        int pn = primNums[i];
        const Bounds3f &bounds = allPrimBounds[pn];
        edges[axis][2 * i] = BoundEdge(bounds.pMin[axis], pn, true);
        edges[axis][2 * i + 1] = BoundEdge(bounds.pMax[axis], pn, false);
    }

    // Sort _edges_ for _axis_
    std::sort(&edges[axis][0], &edges[axis][2 * nPrimitives],
              [](const BoundEdge &e0, const BoundEdge &e1) -> bool {
                  if (e0.t == e1.t)
                      return (int)e0.type < (int)e1.type;
                  else
                      return e0.t < e1.t;
              });

    // Compute cost of all splits for _axis_ to find best
    int nBelow = 0, nAbove = nPrimitives;
    for (int i = 0; i < 2 * nPrimitives; ++i) {
        if (edges[axis][i].type == EdgeType::End) --nAbove;
        Float edgeT = edges[axis][i].t;
        if (edgeT > nodeBounds.pMin[axis] && edgeT < nodeBounds.pMax[axis]) {
            // Compute cost for split at _i_th edge

            // Compute child surface areas for split at _edgeT_
            int otherAxis0 = (axis + 1) % 3, otherAxis1 = (axis + 2) % 3;
            Float belowSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                                 (edgeT - nodeBounds.pMin[axis]) *
                                     (d[otherAxis0] + d[otherAxis1]));
            Float aboveSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                                 (nodeBounds.pMax[axis] - edgeT) *
                                     (d[otherAxis0] + d[otherAxis1]));
            Float pBelow = belowSA * invTotalSA;
            Float pAbove = aboveSA * invTotalSA;
            Float eb = (nAbove == 0 || nBelow == 0) ? emptyBonus : 0;
            Float cost =
                traversalCost +
                isectCost * (1 - eb) * (pBelow * nBelow + pAbove * nAbove);

            // Update best split if this is lowest cost so far
            if (cost < bestCost) {
                bestCost = cost;
                bestAxis = axis;
                bestOffset = i;
            }
        }
        if (edges[axis][i].type == EdgeType::Start) ++nBelow;
    }
    CHECK(nBelow == nPrimitives && nAbove == 0);

    // Create leaf if no good splits were found
    if (bestAxis == -1 && retries < 2) {
        ++retries;
        axis = (axis + 1) % 3;
        goto retrySplit;
    }
    if (bestCost > oldCost) ++badRefines;
    if ((bestCost > 4 * oldCost && nPrimitives < 16) || 
		bestAxis == -1 || 
		badRefines == 3) {
        nodes[nodeNum].InitLeaf(primNums, nPrimitives, &primitiveIndices);
        return;
    }

	//add by popduorz
	float v = 0;
	int index;
	Vector3f vel0, vel1, vel2;

	//const int x = nodes[nodeNum].nPrimitives();
	for (int i = 0; i < nPrimitives; i++) {

		index =  primNums[i];
		const int* verts = primitives[index]->getShape()->getVert();
		vel0 = primitives[index]->getShape()->getMesh()->vel[verts[0]];
		vel1 = primitives[index]->getShape()->getMesh()->vel[verts[1]];
		vel2 = primitives[index]->getShape()->getMesh()->vel[verts[2]];
		
		if (bestAxis == 0) {
			v += 0.5f * vel0.x;
			v += 0.5f * vel1.x;
			v += 0.5f * vel2.x;
		}

		if (bestAxis == 1) {
			v += 0.5f * vel0.y;
			v += 0.5f * vel1.y;
			v += 0.5f * vel2.y;
		}

		if (bestAxis == 2) {
			v += 0.5f * vel0.z;
			v += 0.5f * vel1.z;
			v += 0.5f * vel2.z;
		}
	}
	v = v / 3 / nPrimitives;
	Float split_t0 = edges[bestAxis][bestOffset].t - v;
	Float split_t1 = edges[bestAxis][bestOffset].t + v;
	nodes[nodeNum].split0 = split_t0;
	nodes[nodeNum].split1 = split_t1;
	

	int n0 = 0, n1 = 0;

	//add by popudorz
	std::vector<int> temp_prims0;
	std::vector<int> temp_prims1;
	std::shared_ptr<Shape> shape;

	for (int i = 0; i < bestOffset; i++){

		index = edges[bestAxis][i].primNum;

	    shape = primitives[index]->getShape();//ignore the optimization

		if (edges[bestAxis][i].type == EdgeType::Start) {
			prims0[n0++] = edges[bestAxis][i].primNum;
		}
		
		if (edges[bestAxis][i].type == EdgeType::End) {

				const int* verts = primitives[index]->getShape()->getVert();
				Vector3f tri_v = (shape->getMesh()->vel[verts[0]] +
								  shape->getMesh()->vel[verts[1]] +
								  shape->getMesh()->vel[verts[2]])/3;

				Float split = edges[bestAxis][bestOffset].t;
				Float dist0 = shape->WorldBound0().pMin[bestAxis] - split_t0;
				Float dist1 = shape->WorldBound0().pMax[bestAxis] - split_t0;
				Float dist2 = shape->WorldBound0().pMin[bestAxis] - split_t1;
				Float dist3 = shape->WorldBound0().pMax[bestAxis] - split_t1;

				if (tri_v[bestAxis] >= 0){

					if (dist2 <= 0 && dist3 >= 0) {

						
						//first or second
						temp_prims1.push_back(edges[bestAxis][i].primNum);
						continue;
					}

					if (dist2 > 0) {

						
						temp_prims1.push_back(edges[bestAxis][i].primNum);
					}
				}

				if (tri_v[bestAxis] < 0)
				{
					if (dist0 > 0) {

						
						temp_prims1.push_back(edges[bestAxis][i].primNum);
						continue;
					}

					if (dist0 <= 0 && dist1 >= 0) {

						
						temp_prims1.push_back(edges[bestAxis][i].primNum);

					}
				}
		}
		
	}

	for (int i = bestOffset; i < 2 * nPrimitives; i++) {

		index = edges[bestAxis][i].primNum;

		shape = primitives[index]->getShape();//ignore the optimization

		if (edges[bestAxis][i].type == EdgeType::Start) {
				const int* verts = primitives[index]->getShape()->getVert();
				Vector3f tri_v = (shape->getMesh()->vel[verts[0]] +
								  shape->getMesh()->vel[verts[1]] +
								  shape->getMesh()->vel[verts[2]] ) / 3;

				Float split = edges[bestAxis][bestOffset].t;
				Float dist0 = shape->WorldBound0().pMin[bestAxis] - split_t0;
				Float dist1 = shape->WorldBound0().pMax[bestAxis] - split_t0;
				Float dist2 = shape->WorldBound0().pMin[bestAxis] - split_t1;
				Float dist3 = shape->WorldBound0().pMax[bestAxis] - split_t1;

				if (tri_v[bestAxis] >= 0) {

					if (dist1 < 0) {

						
						//first or second
						temp_prims0.push_back(edges[bestAxis][i].primNum);
						continue;
					}

					if (dist0 <= 0 && dist1 >= 0) {
					
						temp_prims0.push_back(edges[bestAxis][i].primNum);
					}
				}

				if (tri_v[bestAxis] < 0)
				{
					if (dist3 < 0) {

						temp_prims0.push_back(edges[bestAxis][i].primNum);

						continue;
					}

					if (dist2 <= 0 && dist3 >= 0) {

						
						temp_prims0.push_back(edges[bestAxis][i].primNum);

					}
				}

		}

		if (edges[bestAxis][i].type == EdgeType::End) {

			prims1[n1++] = edges[bestAxis][i].primNum;
		}
	}

	if (nodes[nodeNum].tempPrimNum()) {

		int dim_1 = 0;
		int dim_2 = 1;
		int dim_3 = 2;
		float b1_min, b1_max, b2_min, b2_max, b3_min, b3_max;

		//包围盒的三个维度的边界
		b1_min = nodeBounds.pMin.x;
		b1_max = nodeBounds.pMax.x;
		b2_min = nodeBounds.pMin.y;
		b2_max = nodeBounds.pMax.y;
		b3_min = nodeBounds.pMin.z;
		b3_max = nodeBounds.pMax.z;

		for (unsigned int i = 0; i < nodes[nodeNum].tempPrimNum(); i++)
		{
			int tp = nodes[nodeNum].temp_prims_ptr[i];
			std::shared_ptr<Shape> shape = primitives[tp]->getShape();
			float split = edges[bestAxis][bestOffset].t;

			float tri_x_min = tri_min(shape->WorldBound0().pMin[0], shape->WorldBound().pMin[0], shape->WorldBound1().pMin[0]);
			float tri_x_max = tri_max(shape->WorldBound0().pMax[0], shape->WorldBound().pMax[0], shape->WorldBound1().pMax[0]);
			float tri_y_min = tri_min(shape->WorldBound0().pMin[1], shape->WorldBound().pMin[1], shape->WorldBound1().pMin[1]);
			float tri_y_max = tri_max(shape->WorldBound0().pMax[1], shape->WorldBound().pMax[1], shape->WorldBound1().pMax[1]);
			float tri_z_min = tri_min(shape->WorldBound0().pMin[2], shape->WorldBound().pMin[2], shape->WorldBound1().pMin[2]);
			float tri_z_max = tri_max(shape->WorldBound0().pMax[2], shape->WorldBound().pMax[2], shape->WorldBound1().pMax[2]);

			bool flag = false;

			if (b1_min >= tri_x_max || b1_max <= tri_x_min){
				flag = true;
			}
			if (b2_min >= tri_y_max || b2_max <= tri_y_min)
			{
				flag = true;
			}
			if (b3_min >= tri_z_max || b3_max <= tri_z_min)
			{
				flag = true;
			}


			if (shape->WorldBound().pMin[bestAxis] <= split && shape->WorldBound().pMax[bestAxis] >= split) {
				if (!flag) {

					prims0[n0++] = tp;
					prims1[n1++] = tp;
				}
			}

			else {

				const int* verts = primitives[index]->getShape()->getVert();
				Vector3f tri_v = (shape->getMesh()->vel[verts[0]] +
								  shape->getMesh()->vel[verts[1]] +
								  shape->getMesh()->vel[verts[2]]) / 3;

				float dist0 = shape->WorldBound0().pMin[bestAxis] - split_t0;
				float dist1 = shape->WorldBound0().pMax[bestAxis] - split_t0;
				float dist2 = shape->WorldBound1().pMin[bestAxis] - split_t1;
				float dist3 = shape->WorldBound1().pMax[bestAxis] - split_t1;
				float dist4 = shape->WorldBound().pMin[bestAxis] - split;
				float dist5 = shape->WorldBound().pMax[bestAxis] - split;

				if (tri_v[bestAxis] >= 0) {

					if (dist3 < 0) {

						if (!flag) {
							prims0[n0++] = tp;
						}
						continue;
					}

					if (dist1 < 0 && dist2 <= 0 && dist3 >= 0) {

						if (!flag) {

							prims0[n0++] = tp;
							temp_prims1.push_back(tp);
						}
						continue;
					}

					if (dist1<0 && dist2>0) {
						if (dist5 < 0) {
							if (!flag) {
								
								prims0[n0++] = tp;
								temp_prims1.push_back(tp);
							}
							continue;
						}
						if (dist4>0) {
							if (!flag) {
						
								prims1[n1++] = tp;
								temp_prims0.push_back(tp);
							}
							continue;
						}
					}

					if (dist0 <= 0 && dist1 >= 0 && dist2 > 0) {
						if (!flag) {

							prims1[n1++] = tp;
							temp_prims0.push_back(tp);
						}
						continue;
					}

					if (dist0 > 0)
					{
						if (!flag) {

							prims1[n1++] = tp;
						}
					}
				}

				if (tri_v[bestAxis] < 0) {

					if (dist2 > 0) {
						if (!flag) {
						
							prims1[n1++] = tp;
						}
						continue;
					}

					if (dist0 > 0 && dist2 <= 0 && dist3 >= 0) {
						if (!flag) {
						
							prims1[n1++] = tp;
							temp_prims0.push_back(tp);
						}
						continue;
					}

					if (dist0 > 0 && dist3 < 0) {
						if (dist4>0) {
							if (!flag) {
								
								prims1[n1++] = tp;
								temp_prims0.push_back(tp);
							}
							continue;
						}
						if (dist5 < 0) {
							if (!flag) {
								
								prims0[n0++] = tp;
								temp_prims1.push_back(tp);
							}
							continue;
						}
					}

					if (dist0 <= 0 && dist1 >= 0 && dist3 < 0) {
						if (!flag) {
					
							prims0[n0++] = tp;
							temp_prims1.push_back(tp);
						}
						continue;
					}

					if (dist1 < 0) {
						if (!flag) {
						
							prims0[n0++] = tp;
						}
					}
				}
			}
		}

	}

	nodes[nodeNum + 1].temp_prims_ptr = new int[temp_prims0.size()];
	nodes[nodeNum + 1].setTempPrimNum(temp_prims0.size());
	for (int i = 0; i < temp_prims0.size(); i++) {

		nodes[nodeNum + 1].temp_prims_ptr[i] = temp_prims0[i];
	}



	// Recursively initialize children nodes
    Float tSplit = edges[bestAxis][bestOffset].t;
    Bounds3f bounds0 = nodeBounds, bounds1 = nodeBounds;
    bounds0.pMax[bestAxis] = bounds1.pMin[bestAxis] = tSplit;

    buildTree(nodeNum + 1, bounds0, allPrimBounds, prims0, n0, depth - 1, 
			   prims0, prims1 + nPrimitives, badRefines);

	int aboveChild = nextFreeNode;
	nodes[nodeNum].InitInterior(bestAxis, aboveChild, tSplit);

	//add by popduorz
	nodes[aboveChild].temp_prims_ptr = new int[temp_prims1.size()];
	nodes[aboveChild].setTempPrimNum(temp_prims1.size());
	for (int i = 0; i < temp_prims1.size(); i++) {

		nodes[aboveChild].temp_prims_ptr[i] = temp_prims1[i];
	}

	//add by popduorz
	if (nodes[nodeNum].temp_prims_ptr) {
		delete[] nodes[nodeNum].temp_prims_ptr;
		nodes[nodeNum].temp_prims_ptr = NULL;
	}

    buildTree(aboveChild, bounds1, allPrimBounds, prims1, n1, depth - 1,
			   prims0, prims1 + nPrimitives, badRefines);
}



bool KdTreeAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
    ProfilePhase p(Prof::AccelIntersect);
    // Compute initial parametric range of ray inside kd-tree extent
    Float tMin, tMax;

	//add by popduorz
	Bounds3f boundst;
	Interpolate(boundst, bounds0, bounds1, ray.getTime2());
    if (!boundst.IntersectP(ray, &tMin, &tMax)) {
        return false;
    }

    // Prepare to traverse kd-tree for ray
    Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
    PBRT_CONSTEXPR int maxTodo = 64;
    KdToDo todo[maxTodo];
    int todoPos = 0;

    // Traverse kd-tree nodes in order for ray
    bool hit = false;



    const KdAccelNode *node = &nodes[0];
    while (node != nullptr) {

        // Bail out if we found a hit closer than the current node
        if (ray.tMax < tMin) break;
        if (!node->IsLeaf()) {
            // Process kd-tree interior node

            // Compute parametric distance along ray to split plane
            int axis = node->SplitAxis();
           
			//add by popduorz
			Float interplate_split_plane = node->split0 * (1 - ray.getTime2()) +
				node->split1 * ray.getTime2();
			Float tPlane = (interplate_split_plane - ray.o[axis]) * invDir[axis];

            // Get node children pointers for ray
			const KdAccelNode *firstChild, *secondChild;
			int belowFirst =
				(ray.o[axis] < interplate_split_plane) ||
				(ray.o[axis] == interplate_split_plane && ray.d[axis] <= 0);

			if (belowFirst) {
				firstChild = node + 1;
				secondChild = &nodes[node->AboveChild()];
			}
			else {
				firstChild = &nodes[node->AboveChild()];
				secondChild = node + 1;
			}

            // Advance to next child node, possibly enqueue other child
			if (tPlane > tMax || tPlane <= 0) {
				node = firstChild;
			}
			else if (tPlane < tMin) {
				node = secondChild;
			}
            else {
                // Enqueue _secondChild_ in todo list
                todo[todoPos].node = secondChild;
                todo[todoPos].tMin = tPlane;
                todo[todoPos].tMax = tMax;
                ++todoPos;
                node = firstChild;
                tMax = tPlane;

            }
        } 
		else {
            // Check for intersections inside leaf node
            int nPrimitives = node->nPrimitives();
            if (nPrimitives == 1) {
                const std::shared_ptr<Primitive> &p =
                    primitives[node->onePrimitive];
                // Check one primitive inside leaf node
                if (p->Intersect(ray, isect)) hit = true;
            } else {
                for (int i = 0; i < nPrimitives; ++i) {
                    int index =
                        primitiveIndices[node->primitiveIndicesOffset + i];
                    const std::shared_ptr<Primitive> &p = primitives[index];
                    // Check one primitive inside leaf node
                    if (p->Intersect(ray, isect)) hit = true;
                }
            }

            // Grab next node to process from todo list
            if (todoPos > 0) {
                --todoPos;
                node = todo[todoPos].node;
                tMin = todo[todoPos].tMin;
                tMax = todo[todoPos].tMax;
            } else
                break;
        }
    }
    return hit;
}

bool KdTreeAccel::IntersectP(const Ray &ray) const {
    ProfilePhase p(Prof::AccelIntersectP);
    // Compute initial parametric range of ray inside kd-tree extent
    Float tMin, tMax;

	//add by popduorz
	Bounds3f boundst;
	Interpolate(boundst, bounds0, bounds1, ray.getTime2());
    if (!boundst.IntersectP(ray, &tMin, &tMax)) {
        return false;
    }

    // Prepare to traverse kd-tree for ray
    Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
    PBRT_CONSTEXPR int maxTodo = 64;
    KdToDo todo[maxTodo];
    int todoPos = 0;

    const KdAccelNode *node = &nodes[0];
    while (node != nullptr) {
        if (node->IsLeaf()) {
            // Check for shadow ray intersections inside leaf node
            int nPrimitives = node->nPrimitives();
            if (nPrimitives == 1) {
                const std::shared_ptr<Primitive> &p =
                    primitives[node->onePrimitive];
                if (p->IntersectP(ray)) {
                    return true;
                }
            } else {
                for (int i = 0; i < nPrimitives; ++i) {
                    int primitiveIndex =
                        primitiveIndices[node->primitiveIndicesOffset + i];
                    const std::shared_ptr<Primitive> &prim =
                        primitives[primitiveIndex];
                    if (prim->IntersectP(ray)) {
                        return true;
                    }
                }
            }

            // Grab next node to process from todo list
            if (todoPos > 0) {
                --todoPos;
                node = todo[todoPos].node;
                tMin = todo[todoPos].tMin;
                tMax = todo[todoPos].tMax;

            } else
                break;
        } 
		else {
            // Process kd-tree interior node

            // Compute parametric distance along ray to split plane
            int axis = node->SplitAxis();

			//add by popduorz
			Float interplate_split_plane = node->split0 * (1 - ray.getTime2()) + 
										   node->split1 * ray.getTime2();
			Float tPlane = (interplate_split_plane - ray.o[axis]) * invDir[axis];

            // Get node children pointers for ray
            const KdAccelNode *firstChild, *secondChild;

			//modified by popduorz
            int belowFirst =
                (ray.o[axis] < interplate_split_plane) ||
                (ray.o[axis] == interplate_split_plane && ray.d[axis] <= 0);

			if (belowFirst) {
				firstChild = node + 1;
				secondChild = &nodes[node->AboveChild()];
			}
			else {
				firstChild = &nodes[node->AboveChild()];
				secondChild = node + 1;
			}

			// Advance to next child node, possibly enqueue other child
			if (tPlane > tMax || tPlane <= 0)
				node = firstChild;
			else if (tPlane < tMin)
				node = secondChild;
			else {
				// Enqueue _secondChild_ in todo list
				todo[todoPos].node = secondChild;
				todo[todoPos].tMin = tPlane;
				todo[todoPos].tMax = tMax;
				++todoPos;
				node = firstChild;
				tMax = tPlane;
            }
        }
    }
    return false;
}

std::shared_ptr<KdTreeAccel> CreateKdTreeAccelerator(
    const std::vector<std::shared_ptr<Primitive>> &prims, const ParamSet &ps) {
    int isectCost = ps.FindOneInt("intersectcost", 80);
    int travCost = ps.FindOneInt("traversalcost", 1);
    Float emptyBonus = ps.FindOneFloat("emptybonus", 0.5f);
    int maxPrims = ps.FindOneInt("maxprims", 1);
    int maxDepth = ps.FindOneInt("maxdepth", -1);
    return std::make_shared<KdTreeAccel>(prims, isectCost, travCost, emptyBonus,
                                         maxPrims, maxDepth);
}

//add by  popduorz
const std::shared_ptr<Shape> KdTreeAccel::getShape() const {
	return NULL;
}

}  // namespace pbrt
