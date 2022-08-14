#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <stack>

namespace hippo {

static constexpr std::uint8_t kDefaultMaxNodeCount = 8;

template<
    class Data,
    class Scalar,
    std::uint8_t Dimensions,
    std::uint8_t MaxNodeCount = kDefaultMaxNodeCount,
    std::uint8_t MinNodeCount = MaxNodeCount / 2>
class RTree {
 public:
    RTree();
    RTree(const RTree& other) = delete;
    RTree(RTree&& other) = delete;

    virtual ~RTree();

    /**
     * @brief Insert a new data into the tree.
     *
     * The min must be less than max.
     */
    void Insert(const std::array<Scalar, Dimensions>& min, const std::array<Scalar, Dimensions>& max, const Data& data);

    /**
     * @brief Remove the data from the tree.
     * @return true if the data was removed, false if it was not found.
     */
    bool Remove(const std::array<Scalar, Dimensions>& min, const std::array<Scalar, Dimensions>& max, const Data& data);

    /**
     * @brief Find all by the check.
     * @tparam Overlap If the sub rect is overlapped, the rect must be overlapped.
     * @param callback Return true to continue, false to stop.
     * @return The count of found data.
     */
    template<class T, class Overlap>
    std::size_t Search(const T& check, const std::function<bool(const Data&)>& callback) const;

    /**
     * @brief Remove all entries from tree.
     */
    void RemoveAll();

    /**
     * @brief Count the data elements in this container. This is slow since no internal counter is maintained.
     */
    [[nodiscard]] std::size_t Count() const;

 private:
    static constexpr std::array<float, 21> kUnitSphereVolumes = {
        0.000000F, 2.000000F, 3.141593F,  // Dimension 0,1,2
        4.188790F, 4.934802F, 5.263789F,  // Dimension 3,4,5
        5.167713F, 4.724766F, 4.058712F,  // Dimension 6,7,8
        3.298509F, 2.550164F, 1.884104F,  // Dimension 9,10,11
        1.335263F, 0.910629F, 0.599265F,  // Dimension 12,13,14
        0.381443F, 0.235331F, 0.140981F,  // Dimension 15,16,17
        0.082146F, 0.046622F, 0.025807F,  // Dimension 18,19,20
    };

    static constexpr float kUnitSphereVolume = kUnitSphereVolumes[Dimensions];
    static constexpr std::uint8_t kNotToken = std::numeric_limits<std::uint8_t>::max();

    struct Rect {
        std::array<Scalar, Dimensions> min;
        std::array<Scalar, Dimensions> max;
    };

    struct Node;

    struct Branch {
        Rect rect;
        Data data;
        Node* child = nullptr;
    };

    // Node for each branch level.
    struct Node {
        std::uint8_t count = 0;
        std::size_t level = 0;
        std::array<Branch, MaxNodeCount> branches;
    };

    // Variables for finding a split partition.
    struct PartitionInfo {
        std::array<std::uint8_t, MaxNodeCount + 1> partitions;
        std::array<std::uint8_t, 2> counts;
        std::array<Rect, 2> covers;
        std::array<Scalar, 2> areas;
        std::array<Branch, MaxNodeCount + 1> branches;
    };

    Node* root_;

    static bool IsInternalNode(Node* node) {
        return node->level > 0;
    }

    static bool IsLeafNode(Node* node) {
        return node->level == 0;
    }

    Node* AllocNode();
    void FreeNode(Node* node);
    bool InsertRect(const Branch& branch, Node** root, std::size_t level);
    bool InsertRectRec(const Branch& branch, Node* node, Node** newNode, std::size_t level);
    Rect NodeCover(Node* node);
    bool AddBranch(const Branch& branch, Node* node, Node** newNode);
    void DisconnectBranch(Node* node, std::uint8_t index);
    std::uint8_t PickBranch(const Rect& rect, Node* node);
    Rect CombineRect(const Rect& rectA, const Rect& rectB);
    void SplitNode(Node* node, const Branch& branch, Node** newNode);
    Scalar CalcRectVolume(const Rect& rect);
    void GetBranches(Node* node, const Branch& branch, PartitionInfo* info);
    void ChoosePartition(PartitionInfo* info);
    void LoadNodes(Node* nodeA, Node* nodeB, PartitionInfo* info);
    void InitPartitionInfo(PartitionInfo* info);
    void PickSeeds(PartitionInfo* info);
    void Classify(std::uint8_t index, std::uint8_t group, PartitionInfo* info);
    bool RemoveRect(const Rect& rect, const Data& data, Node** root);
    bool RemoveRectRec(const Rect& rect, const Data& data, Node* node, std::stack<Node*>& reInsertNodes);
    bool Overlap(const Rect& rectA, const Rect& rectB) const;
    template<class T, class Overlap>
    bool Search(Node* node, const T& check, std::size_t& foundCount, const std::function<bool(const Data&)>& callback)
        const;
    void RemoveAllRec(Node* node);
    void Reset();
    void CountRec(Node* node, std::size_t& count) const;

    static_assert(std::numeric_limits<Scalar>::is_iec559, "'Scalar' accepts floating-point types only");
    static_assert(Dimensions > 0, "Dimensions must be greater than 0");
    static_assert(Dimensions < kUnitSphereVolumes.size(), "Dimensions must be less than or equal to 20");
    static_assert(MaxNodeCount > MinNodeCount, "'MaxNodeCount' must be greater than 'MinNodeCount'");
    static_assert(MinNodeCount, "'MinNodeCount' must be greater than 0");
};

#ifdef NDEBUG
#    define RTREE_ASSERT_RECT(min, max) ((void) 0)
#else
#    define RTREE_ASSERT_RECT(min, max)                 \
        for (std::uint8_t i = 0; i < Dimensions; ++i) { \
            assert((min)[i] < (max)[i]);                \
        }                                               \
        ((void) 0)
#endif

#define RTREE_TEMPLATE \
    template<class Data, class Scalar, std::uint8_t Dimensions, std::uint8_t MaxNodeCount, std::uint8_t MinNodeCount>
#define RTREE_TYPE RTree<Data, Scalar, Dimensions, MaxNodeCount, MinNodeCount>

RTREE_TEMPLATE
RTREE_TYPE::RTree() {
    root_ = AllocNode();
}

RTREE_TEMPLATE
RTREE_TYPE::~RTree() {
    Reset();
}

RTREE_TEMPLATE
void RTREE_TYPE::Insert(
    const std::array<Scalar, Dimensions>& min,
    const std::array<Scalar, Dimensions>& max,
    const Data& data
) {
    RTREE_ASSERT_RECT(min, max);

    Branch branch;
    branch.rect.min = min;
    branch.rect.max = max;
    branch.data = data;
    branch.child = nullptr;

    InsertRect(branch, &root_, 0);
}

RTREE_TEMPLATE
bool RTREE_TYPE::Remove(
    const std::array<Scalar, Dimensions>& min,
    const std::array<Scalar, Dimensions>& max,
    const Data& data
) {
    RTREE_ASSERT_RECT(min, max);

    Rect rect;
    rect.min = min;
    rect.max = max;

    return RemoveRect(rect, data, &root_);
}

RTREE_TEMPLATE
template<class T, class Overlap>
std::size_t RTREE_TYPE::Search(const T& check, const std::function<bool(const Data&)>& callback) const {
    std::size_t foundCount = 0;
    Search<T, Overlap>(root_, check, foundCount, callback);
    return foundCount;
}

RTREE_TEMPLATE
std::size_t RTREE_TYPE::Count() const {
    std::size_t count = 0;
    CountRec(root_, count);
    return count;
}

RTREE_TEMPLATE
void RTREE_TYPE::CountRec(Node* node, std::size_t& count) const {
    if (IsInternalNode(node)) {  // not a leaf node
        for (std::uint8_t index = 0; index < node->count; ++index) {
            CountRec(node->branches[index].child, count);
        }
    } else {  // A leaf node
        count += node->count;
    }
}

RTREE_TEMPLATE
void RTREE_TYPE::RemoveAll() {
    // Return if it has no node
    if (root_->count == 0) {
        return;
    }

    // Delete all existing nodes
    Reset();

    root_ = AllocNode();
}

RTREE_TEMPLATE
void RTREE_TYPE::Reset() {
    RemoveAllRec(root_);
}

RTREE_TEMPLATE
void RTREE_TYPE::RemoveAllRec(Node* node) {
    assert(node);
    assert(node->level >= 0);

    if (IsInternalNode(node)) {  // This is an internal node in the tree
        for (std::uint8_t index = 0; index < node->count; ++index) {
            RemoveAllRec(node->branches[index].child);
        }
    }
    FreeNode(node);
}

RTREE_TEMPLATE
typename RTREE_TYPE::Node* RTREE_TYPE::AllocNode() {
    // TODO: custom allocator
    return new Node;
}

RTREE_TEMPLATE
void RTREE_TYPE::FreeNode(Node* node) {
    // TODO: custom allocator
    delete node;
}

// Inserts a new data rectangle into the index structure.
// Recursively descends tree, propagates splits back up.
// Returns 0 if node was not split and old node updated.
// If node was split, returns 1 and sets the pointer pointed to by
// new_node to point to the new node. Old node updated to become one of two.
// The level argument specifies the number of steps up from the leaf
// level to insert; e.g. a data rectangle goes in at level = 0.
RTREE_TEMPLATE
bool RTREE_TYPE::InsertRectRec(const Branch& branch, Node* node, Node** newNode, std::size_t level) {
    assert(node);
    assert(newNode);
    assert(level <= node->level);

    // recurse until we reach the correct level for the new record. data records
    // will always be called with a_level == 0 (leaf)
    if (node->level > level) {
        // Still above level for insertion, go down tree recursively
        Node* otherNode;

        // find the optimal branch for this record
        std::uint8_t index = PickBranch(branch.rect, node);

        // recursively insert this record into the picked branch
        bool childWasSplit = InsertRectRec(branch, node->branches[index].child, &otherNode, level);

        if (!childWasSplit) {
            // Child was not split. Merge the bounding box of the new record with the
            // existing bounding box
            node->branches[index].rect = CombineRect(branch.rect, node->branches[index].rect);
            return false;
        } else {
            // Child was split. The old branches are now re-partitioned to two nodes,
            // so we have to re-calculate the bounding boxes of each node
            node->branches[index].rect = NodeCover(node->branches[index].child);
            Branch newBranch;
            newBranch.child = otherNode;
            newBranch.rect = NodeCover(otherNode);

            // The old node is already a child of a_node. Now add the newly-created
            // node to a_node as well. a_node might be split because of that.
            return AddBranch(newBranch, node, newNode);
        }
    } else if (node->level == level) {
        // We have reached level for insertion. Add rect, split if necessary
        return AddBranch(branch, node, newNode);
    } else {
        // Should never occur
        assert(false);
        return false;
    }
}

// Insert a data rectangle into an index structure.
// InsertRect provides for splitting the root;
// returns 1 if root was split, 0 if it was not.
// The level argument specifies the number of steps up from the leaf
// level to insert; e.g. a data rectangle goes in at level = 0.
RTREE_TEMPLATE
bool RTREE_TYPE::InsertRect(const Branch& branch, Node** root, std::size_t level) {
    assert(root);
    assert(level <= (*root)->level);
    RTREE_ASSERT_RECT(branch.rect.min, branch.rect.max);

    Node* newNode;

    if (InsertRectRec(branch, *root, &newNode, level)) {
        // Root split

        // Grow tree taller and new root
        Node* newRoot = AllocNode();
        newRoot->level = (*root)->level + 1;

        Branch newBranch;

        // add old root node as a child of the new root
        newBranch.rect = NodeCover(*root);
        newBranch.child = *root;
        // It's certain the max node count is greater than 1.
        // No new node created here. It's safe to pass nullptr.
        AddBranch(newBranch, newRoot, nullptr);

        // add the split node as a child of the new root
        newBranch.rect = NodeCover(newNode);
        newBranch.child = newNode;
        // It's certain the max node count is greater than 1.
        // No new node created here. It's safe to pass nullptr.
        AddBranch(newBranch, newRoot, nullptr);

        // set the new root as the root node
        *root = newRoot;

        return true;
    }

    return false;
}

// Find the smallest rectangle that includes all rectangles in branches of a node.
RTREE_TEMPLATE
typename RTREE_TYPE::Rect RTREE_TYPE::NodeCover(Node* node) {
    assert(node);
    assert(node->count > 0);

    Rect rect = node->branches[0].rect;
    for (std::uint8_t index = 1; index < node->count; ++index) {
        rect = CombineRect(rect, node->branches[index].rect);
    }

    return rect;
}

// Add a branch to a node. Split the node if necessary.
// Returns 0 if node not split. Old node updated.
// Returns 1 if node split, sets *new_node to address of new node.
// Old node updated, becomes one of two.
RTREE_TEMPLATE
bool RTREE_TYPE::AddBranch(const Branch& branch, Node* node, Node** newNode) {
    assert(node);

    if (node->count < MaxNodeCount) {
        // Split won't be necessary
        node->branches[node->count] = branch;
        ++node->count;
        return false;
    } else {
        assert(newNode);

        SplitNode(node, branch, newNode);
        return true;
    }
}

// Disconnect a dependent node.
// Caller must return (or stop using iteration index) after this as count has changed
RTREE_TEMPLATE
void RTREE_TYPE::DisconnectBranch(Node* node, std::uint8_t index) {
    assert(node && (index >= 0) && (index < MaxNodeCount));
    assert(node->count > 0);

    // Remove element by swapping with the last element to prevent gaps in array
    node->branches[index] = node->branches[node->count - 1];

    // Clear the data of the last element
    node->branches[node->count - 1].data = Data();

    --node->count;
}

// Pick a branch. Pick the one that will need the smallest increase
// in area to accommodate the new rectangle. This will result in the
// least total area for the covering rectangles in the current node.
// In case of a tie, pick the one which was smaller before, to get
// the best resolution when searching.
RTREE_TEMPLATE
std::uint8_t RTREE_TYPE::PickBranch(const Rect& rect, Node* node) {
    assert(node);
    assert(node->count > 0);

    std::uint8_t best;
    Scalar bestArea;
    Scalar bestIncrease;

    {
        const Rect& currentRect = node->branches[0].rect;
        Scalar area = CalcRectVolume(currentRect);
        Rect combinedRect = CombineRect(rect, currentRect);
        Scalar increase = CalcRectVolume(combinedRect) - area;

        best = 0;
        bestArea = area;
        bestIncrease = increase;
    }

    for (std::uint8_t index = 1; index < node->count; ++index) {
        const Rect& currentRect = node->branches[index].rect;
        Scalar area = CalcRectVolume(currentRect);
        Rect combinedRect = CombineRect(rect, currentRect);
        Scalar increase = CalcRectVolume(combinedRect) - area;

        if ((increase < bestIncrease) || ((increase == bestIncrease) && (area < bestArea))) {
            best = index;
            bestArea = area;
            bestIncrease = increase;
        }
    }
    return best;
}

// Combine two rectangles into larger one containing both.
RTREE_TEMPLATE
typename RTREE_TYPE::Rect RTREE_TYPE::CombineRect(const Rect& rectA, const Rect& rectB) {
    Rect newRect;

    for (std::uint8_t index = 0; index < Dimensions; ++index) {
        newRect.min[index] = std::min(rectA.min[index], rectB.min[index]);
        newRect.max[index] = std::max(rectA.max[index], rectB.max[index]);
    }

    return newRect;
}

// Split a node.
// Divides the nodes branches and the extra one between two nodes.
// Old node is one of the new ones, and one really new one is created.
// Tries more than one method for choosing a partition, uses best result.
RTREE_TEMPLATE
void RTREE_TYPE::SplitNode(Node* node, const Branch& branch, Node** newNode) {
    assert(node);

    // Could just use local here, but member or external is faster since it is reused
    PartitionInfo info;

    // Load all the branches into a buffer, initialize old node
    GetBranches(node, branch, &info);

    // Find partition
    ChoosePartition(&info);

    // Create a new node to hold (about) half of the branches
    *newNode = AllocNode();
    (*newNode)->level = node->level;

    // Reset data if it's a leaf
    if (IsLeafNode(node)) {
        for (std::uint8_t index = 0; index < node->count; ++index) {
            node->branches[index].data = Data();
        }
    }

    // Put branches from buffer into 2 nodes according to the chosen partition
    node->count = 0;
    LoadNodes(node, *newNode, &info);

    assert((node->count + (*newNode)->count) == MaxNodeCount + 1);
}

// The exact volume of the bounding sphere for the given Rect
RTREE_TEMPLATE
Scalar RTREE_TYPE::CalcRectVolume(const Rect& rect) {
    Scalar sumOfSquares = 0;
    for (std::uint8_t index = 0; index < Dimensions; ++index) {
        static constexpr Scalar kHalfScale = 0.5;
        Scalar halfExtent = (rect.max[index] - rect.min[index]) * kHalfScale;
        sumOfSquares += halfExtent * halfExtent;
    }

    Scalar radius = std::sqrt(sumOfSquares);

    // Pow maybe slow, so test for common dims like 2,3 and just use x*x, x*x*x.
    if (Dimensions == 2) {
        return (radius * radius * kUnitSphereVolume);
    } else if (Dimensions == 3) {
        return (radius * radius * radius * kUnitSphereVolume);
    } else {
        return std::pow(radius, Dimensions) * kUnitSphereVolume;
    }
}

// Load branch buffer with branches from full node plus the extra branch.
RTREE_TEMPLATE
void RTREE_TYPE::GetBranches(Node* node, const Branch& branch, PartitionInfo* info) {
    assert(node);
    assert(node->count == MaxNodeCount);

    // Load the branch buffer
    for (std::uint8_t index = 0; index < MaxNodeCount; ++index) {
        info->branches[index] = node->branches[index];
    }
    info->branches[MaxNodeCount] = branch;
}

// Method #0 for choosing a partition:
// As the seeds for the two groups, pick the two rects that would waste the
// most area if covered by a single rectangle, i.e. evidently the worst pair
// to have in the same group.
// Of the remaining, one at a time is chosen to be put in one of the two groups.
// The one chosen is the one with the greatest difference in area expansion
// depending on which group - the rect most strongly attracted to one group
// and repelled from the other.
// If one group gets too full (more would force other group to violate min
// fill requirement) then other group gets the rest.
// These last are the ones that can go in either group most easily.
RTREE_TEMPLATE
void RTREE_TYPE::ChoosePartition(PartitionInfo* info) {
    assert(info);

    InitPartitionInfo(info);
    PickSeeds(info);

    while (((info->counts[0] + info->counts[1]) < MaxNodeCount + 1) &&
           (info->counts[0] < (MaxNodeCount + 1 - MinNodeCount)) &&
           (info->counts[1] < (MaxNodeCount + 1 - MinNodeCount))) {
        bool first = true;
        Scalar biggestDiff = 0;
        std::uint8_t chosen = 0;
        std::uint8_t betterGroup = 0;

        for (std::uint8_t index = 0; index < MaxNodeCount + 1; ++index) {
            if (info->partitions[index] != kNotToken) {
                continue;
            }

            const Rect& curRect = info->branches[index].rect;
            Rect rect0 = CombineRect(curRect, info->covers[0]);
            Rect rect1 = CombineRect(curRect, info->covers[1]);
            Scalar growth0 = CalcRectVolume(rect0) - info->areas[0];
            Scalar growth1 = CalcRectVolume(rect1) - info->areas[1];
            Scalar diff = growth1 - growth0;

            std::uint8_t group;
            if (diff >= 0) {
                group = 0;
            } else {
                group = 1;
                diff = -diff;
            }

            if (first || diff > biggestDiff) {
                first = false;
                biggestDiff = diff;
                chosen = index;
                betterGroup = group;
            } else if ((diff == biggestDiff) && (info->counts[group] < info->counts[betterGroup])) {
                chosen = index;
                betterGroup = group;
            }
        }
        assert(!first);

        Classify(chosen, betterGroup, info);
    }

    // If one group too full, put remaining rects in the other
    if ((info->counts[0] + info->counts[1]) < MaxNodeCount + 1) {
        std::uint8_t group = info->counts[0] >= MaxNodeCount + 1 - MinNodeCount ? 1 : 0;
        for (std::uint8_t index = 0; index < MaxNodeCount + 1; ++index) {
            if (info->partitions[index] == kNotToken) {
                Classify(index, group, info);
            }
        }
    }

    assert((info->counts[0] + info->counts[1]) == MaxNodeCount + 1);
    assert((info->counts[0] >= MinNodeCount) && (info->counts[1] >= MinNodeCount));
}

// Copy branches from the buffer into two nodes according to the partition.
RTREE_TEMPLATE
void RTREE_TYPE::LoadNodes(Node* nodeA, Node* nodeB, PartitionInfo* info) {
    assert(nodeA);
    assert(nodeB);
    assert(info);

    for (std::uint8_t index = 0; index < MaxNodeCount + 1; ++index) {
        assert(info->partitions[index] == 0 || info->partitions[index] == 1);

        Node* targetNode = info->partitions[index] == 0 ? nodeA : nodeB;

        // It is assured that AddBranch here will not cause a node split.
        [[maybe_unused]] bool nodeWasSplit = AddBranch(info->branches[index], targetNode, nullptr);
        assert(!nodeWasSplit);
    }
}

// Initialize a PartitionInfo structure.
RTREE_TEMPLATE
void RTREE_TYPE::InitPartitionInfo(PartitionInfo* info) {
    assert(info);

    info->counts[0] = info->counts[1] = 0;
    info->areas[0] = info->areas[1] = 0;
    for (std::uint8_t index = 0; index < MaxNodeCount + 1; ++index) {
        info->partitions[index] = kNotToken;
    }
}

RTREE_TEMPLATE
void RTREE_TYPE::PickSeeds(PartitionInfo* info) {
    std::uint8_t seed0 = 0;
    std::uint8_t seed1 = 0;
    std::array<Scalar, MaxNodeCount + 1> areas;

    for (std::uint8_t index = 0; index < MaxNodeCount + 1; ++index) {
        areas[index] = CalcRectVolume(info->branches[index].rect);
    }

    bool first = true;
    Scalar worst = 0;
    for (std::uint8_t indexA = 0; indexA < MaxNodeCount; ++indexA) {
        for (std::uint8_t indexB = indexA + 1; indexB < MaxNodeCount + 1; ++indexB) {
            Rect oneRect = CombineRect(info->branches[indexA].rect, info->branches[indexB].rect);
            Scalar waste = CalcRectVolume(oneRect) - areas[indexA] - areas[indexB];
            if (first || waste > worst) {
                first = false;
                worst = waste;
                seed0 = indexA;
                seed1 = indexB;
            }
        }
    }
    assert(!first);

    Classify(seed0, 0, info);
    Classify(seed1, 1, info);
}

// Put a branch in one of the groups.
RTREE_TEMPLATE
void RTREE_TYPE::Classify(std::uint8_t index, std::uint8_t group, PartitionInfo* info) {
    assert(info);
    assert(info->partitions[index] == kNotToken);

    info->partitions[index] = group;

    // Calculate combined rect
    if (info->counts[group] == 0) {
        info->covers[group] = info->branches[index].rect;
    } else {
        info->covers[group] = CombineRect(info->branches[index].rect, info->covers[group]);
    }

    // Calculate volume of combined rect
    info->areas[group] = CalcRectVolume(info->covers[group]);

    ++info->counts[group];
}

// Delete a data rectangle from an index structure.
// Pass in a pointer to a Rect, the tid of the record, ptr to root node.
// Returns false if record not found, true if success.
// RemoveRect provides for eliminating the root.
RTREE_TEMPLATE
bool RTREE_TYPE::RemoveRect(const Rect& rect, const Data& data, Node** root) {
    assert(root);
    assert(*root);

    std::stack<Node*> reInsertNodes;

    if (RemoveRectRec(rect, data, *root, reInsertNodes)) {
        // Found and deleted a data item
        // Reinsert any branches from eliminated nodes
        while (!reInsertNodes.empty()) {
            Node* reInsertNode = reInsertNodes.top();
            reInsertNodes.pop();

            for (std::uint8_t index = 0; index < reInsertNode->count; ++index) {
                InsertRect(reInsertNode->branches[index], root, reInsertNode->level);
            }

            FreeNode(reInsertNode);
        }

        // Check for redundant root (not leaf, 1 child) and eliminate TODO replace
        // if with while? In case there is a whole branch of redundant roots...
        if ((*root)->count == 1 && IsInternalNode(*root)) {
            Node* tempNode = (*root)->branches[0].child;

            assert(tempNode);
            FreeNode(*root);
            *root = tempNode;
        }
        return true;
    } else {
        return false;
    }
}

// Delete a rectangle from non-root part of an index structure.
// Called by RemoveRect.  Descends tree recursively,
// merges branches on the way back up.
// Returns false if record not found, true if success.
RTREE_TEMPLATE
bool RTREE_TYPE::RemoveRectRec(const Rect& rect, const Data& data, Node* node, std::stack<Node*>& reInsertNodes) {
    assert(node);
    assert(node->level >= 0);

    if (IsInternalNode(node)) {  // not a leaf node
        for (std::uint8_t index = 0; index < node->count; ++index) {
            if (Overlap(rect, node->branches[index].rect)) {
                if (RemoveRectRec(rect, data, node->branches[index].child, reInsertNodes)) {
                    if (node->branches[index].child->count >= MinNodeCount) {
                        // child removed, just resize parent rect
                        node->branches[index].rect = NodeCover(node->branches[index].child);
                    } else {
                        // child removed, not enough entries in node, eliminate node
                        reInsertNodes.push(node->branches[index].child);
                        DisconnectBranch(node, index);  // Must return after this call as count has changed
                    }
                    return true;
                }
            }
        }
        return false;
    } else {  // A leaf node
        for (std::uint8_t index = 0; index < node->count; ++index) {
            if (node->branches[index].data == data) {
                DisconnectBranch(node, index);  // Must return after this call as count has changed
                return true;
            }
        }
        return false;
    }
}

// Decide whether two rectangles overlap.
RTREE_TEMPLATE
bool RTREE_TYPE::Overlap(const Rect& rectA, const Rect& rectB) const {
    for (std::uint8_t index = 0; index < Dimensions; ++index) {
        if (rectA.min[index] >= rectB.max[index] || rectB.min[index] >= rectA.max[index]) {
            return false;
        }
    }
    return true;
}

RTREE_TEMPLATE
template<class T, class Overlap>
bool RTREE_TYPE::Search(
    Node* node,
    const T& check,
    std::size_t& foundCount,
    const std::function<bool(const Data&)>& callback
) const {
    assert(node);
    assert(node->level >= 0);

    if (IsInternalNode(node)) {  // This is an internal node in the tree
        for (std::uint8_t index = 0; index < node->count; ++index) {
            if (Overlap()(check, node->branches[index].rect.min, node->branches[index].rect.max)) {
                if (!Search<T, Overlap>(node->branches[index].child, check, foundCount, callback)) {
                    // The callback indicated to stop searching
                    return false;
                }
            }
        }
    } else {  // This is a leaf node
        for (std::uint8_t index = 0; index < node->count; ++index) {
            if (Overlap()(check, node->branches[index].rect.min, node->branches[index].rect.max)) {
                const Data& data = node->branches[index].data;
                ++foundCount;

                if (callback && !callback(data)) {
                    return false;  // Don't continue searching
                }
            }
        }
    }

    return true;  // Continue searching
}

#undef RTREE_TYPE
#undef RTREE_TEMPLATE

}  // namespace hippo
