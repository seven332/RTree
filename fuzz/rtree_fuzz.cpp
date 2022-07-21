#include <array>
#include <utility>
#include <vector>

#include "rtree/rtree.h"

struct Rect {
    float l, r, b, t;
};

class RTreeTest {
 public:
    [[nodiscard]] std::size_t Count() const {
        return data_.size();
    }

    void Insert(const Rect& rect, int data) {
        data_.emplace_back(rect, data);

        std::array<float, 2> min {rect.l, rect.t};
        std::array<float, 2> max {rect.r, rect.b};
        tree_.Insert(min, max, data);

        assert(data_.size() == tree_.Count());
    }

    void Remove(int index) {
        Rect rect = data_[index].first;
        int data = data_[index].second;
        data_.erase(data_.begin() + index);

        std::array<float, 2> min {rect.l, rect.t};
        std::array<float, 2> max {rect.r, rect.b};
        tree_.Remove(min, max, data);

        assert(data_.size() == tree_.Count());
    }

    static bool Overlap(const Rect& rectA, const Rect& rectB) {
        float l = std::max(rectA.l, rectB.l);
        float r = std::min(rectA.r, rectB.r);
        float t = std::max(rectA.t, rectB.t);
        float b = std::min(rectA.b, rectB.b);
        return l < r && t < b;
    }

    void Search(const Rect& rect) {
        std::vector<int> dResult;
        for (auto& it : data_) {
            if (Overlap(it.first, rect)) {
                dResult.push_back(it.second);
            }
        }
        std::sort(dResult.begin(), dResult.end());

        std::vector<int> tResult;
        std::array<float, 2> min {rect.l, rect.t};
        std::array<float, 2> max {rect.r, rect.b};
        tree_.Search(min, max, [&tResult](const int& data) {
            tResult.push_back(data);
            return true;
        });
        std::sort(tResult.begin(), tResult.end());

        assert(dResult == tResult);
    }

    void RemoveAll() {
        data_.clear();
        tree_.RemoveAll();

        assert(data_.empty());
        assert(tree_.Count() == 0);
    }

 private:
    std::vector<std::pair<Rect, int>> data_;
    hippo::RTree<int, float, 2> tree_;
};

static bool NextChar(const std::uint8_t** data, std::size_t* size, char* value) {
    if (*size == 0) {
        return false;
    }
    *value = *reinterpret_cast<const char*>(*data);
    (*data)++;
    (*size)--;
    return true;
}

static bool NextInt(const std::uint8_t** data, std::size_t* size, int* value) {
    if (*size < sizeof(int)) {
        return false;
    }
    *value = *reinterpret_cast<const int*>(*data);
    (*data) += sizeof(int);
    (*size) -= sizeof(int);
    return true;
}

static bool NextFloat(const std::uint8_t** data, std::size_t* size, float* value) {
    while (true) {
        if (*size < sizeof(float)) {
            return false;
        }
        *value = *reinterpret_cast<const float*>(*data);
        if (std::isnan(*value) || std::isinf(*value)) {
            *data += 1;
            *size -= 1;
        } else {
            *value = *value;
            *data += sizeof(float);
            *size -= sizeof(float);
            return true;
        }
    }
}

static bool NextRect(const std::uint8_t** data, std::size_t* size, Rect* value) {
    if (!NextFloat(data, size, &value->l) || !NextFloat(data, size, &value->r) || !NextFloat(data, size, &value->b) ||
        !NextFloat(data, size, &value->t)) {
        return false;
    }
    if (value->l == value->r || value->t == value->b) {
        return false;
    }

    auto l = std::min(value->l, value->r);
    auto r = std::max(value->l, value->r);
    auto t = std::min(value->t, value->b);
    auto b = std::max(value->t, value->b);

    value->l = l;
    value->r = r;
    value->t = t;
    value->b = b;

    return true;
}

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size) {
    RTreeTest test;

    // Get the rect to search
    Rect searchRect {};
    if (!NextRect(&data, &size, &searchRect)) {
        return 0;
    }

    // Modify the tree
    int nextIntData = 0;
    while (true) {
        bool shouldInsert;
        if (test.Count() == 0) {
            shouldInsert = true;
        } else {
            char flag;
            if (!NextChar(&data, &size, &flag)) {
                break;
            }
            shouldInsert = (flag % 3) < 2;
        }

        if (shouldInsert) {
            int intData = nextIntData++;
            assert(intData >= 0);
            Rect rect {};
            if (!NextRect(&data, &size, &rect)) {
                break;
            }

            test.Insert(rect, intData);
        } else {
            int index;
            if (!NextInt(&data, &size, &index)) {
                break;
            }
            index = std::abs(index);
            index = index % static_cast<int>(test.Count());
            index = std::abs(index);

            test.Remove(index);
        }
    }

    test.Search(searchRect);

    test.RemoveAll();

    return 0;
}
