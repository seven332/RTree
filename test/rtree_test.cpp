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

TEST(RTreeTest, DoNotCopyTwiceAfterInsert) {
    RTree<std::shared_ptr<int>, float, 2> tree;

    std::array<float, 2> min = {1, 1};
    std::array<float, 2> max = {2, 2};

    static constexpr int kDataCount = 100;
    std::vector<std::shared_ptr<int>> dataVector;
    dataVector.reserve(kDataCount);
    for (int i = 0; i < kDataCount; ++i) {
        dataVector.push_back(std::make_shared<int>(i));
    }

    for (const auto& data : dataVector) {
        tree.Insert(min, max, data);
    }

    for (const auto& data : dataVector) {
        EXPECT_EQ(data.use_count(), 2);
    }
}

TEST(RTreeTest, RemoveAllOnEmptyRTree) {
    RTree<std::shared_ptr<int>, float, 2> tree;
    tree.RemoveAll();
    EXPECT_EQ(tree.Count(), 0);
}

}  // namespace hippo
