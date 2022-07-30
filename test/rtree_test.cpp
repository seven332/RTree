#include "rtree/rtree.h"

#include "gtest/gtest.h"

namespace hippo {

TEST(RTreeTest, DestoryDataAfterRemove) {
    RTree<std::shared_ptr<int>, float, 2> tree;
    auto data = std::make_shared<int>(1);
    std::array<float, 2> min = {1, 1};
    std::array<float, 2> max = {2, 2};

    tree.Insert(min, max, data);
    EXPECT_EQ(data.use_count(), 2);

    tree.Remove(min, max, data);
    EXPECT_EQ(data.use_count(), 1);

    tree.Insert(min, max, data);
    EXPECT_EQ(data.use_count(), 2);

    tree.RemoveAll();
    EXPECT_EQ(data.use_count(), 1);
}

}  // namespace hippo
