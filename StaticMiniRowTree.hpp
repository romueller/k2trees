#ifndef K2TREES_MINIBASICROWTREE_HPP
#define K2TREES_MINIBASICROWTREE_HPP

#include <queue>

#include "RowTree.hpp"
#include "Utility.hpp"

template<typename E>
class MiniRowTree : public virtual RowTree<E> {

public:
    typedef E elem_type;

    typedef typename RowTree<elem_type>::list_type list_type;


    MiniRowTree() {
        // nothing to do
    }

    MiniRowTree(const MiniRowTree& other) {

        null_ = other.null_;

        length_ = other.length_;
        positions_ = new size_type[length_];
        values_ = new elem_type[length_];
        for (auto i = 0; i < length_; i++) {

            positions_[i] = other.positions_[i];
            values_[i] = other.values_[i];

        }

    }

    MiniRowTree& operator=(const MiniRowTree& other) {

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        delete positions_;
        delete values_;

        null_ = other.null_;

        length_ = other.length_;
        positions_ = new size_type[length_];
        values_ = new elem_type[length_];
        for (auto i = 0; i < length_; i++) {

            positions_[i] = other.positions_[i];
            values_[i] = other.values_[i];

        }

        return *this;

    }

    MiniRowTree(const std::vector<elem_type>& v, const elem_type null = elem_type()) {

        null_ = null;

        length_ = 0;
        for (size_type i = 0; i < v.size(); i++) {
            length_ += (v[i] != null);
        }

        positions_ = new size_type[length_];
        values_ = new elem_type[length_];

        for (size_type i = 0, pos = 0; i < v.size(); i++) {

            if (v[i] != null) {
                positions_[pos] = i;
                values_[pos++] = v[i];
            }

        }

    }

    MiniRowTree(const list_type& list, const elem_type null = elem_type()) {

        null_ = null;

        length_ = list.size();

        positions_ = new size_type[length_];
        values_ = new elem_type[length_];

        for (size_type i = 0; i < length_; i++) {

            positions_[i] = list[i].first;
            values_[i] = list[i].second;

        }

    }

    ~MiniRowTree() {

        delete positions_;
        delete values_;

    }


    size_type getLength() override {
        return length_;
    }

    elem_type getNull() override {
        return null_;
    }


    bool isNotNull(size_type i) override {
        return std::find(positions_, positions_ + length_, i) != positions_ + length_;
    }

    elem_type getElement(size_type i) override {

        auto iter = std::find(positions_, positions_ + length_, i);

        return (iter != positions_ + length_) ? values_[iter - positions_] : null_;

    }

    size_type getFirst() override {

        size_type min = -1;
        for (auto i = 0; i < length_; i++) {
            min = std::min(min, positions_[i]);
        }

        return min;

    }

    std::vector<elem_type> getElementsInRange(size_type l, size_type r) override {

        std::vector<elem_type> elems;
        for (size_type i = 0; i < length_; i++) {
            if (l <= positions_[i] && positions_[i] <= r) {
                elems.push_back(values_[i]);
            }
        }

        return elems;

    }

    std::vector<size_type> getPositionsInRange(size_type l, size_type r) override {

        std::vector<size_type> positions;
        for (size_type i = 0; i < length_; i++) {
            if (l <= positions_[i] && positions_[i] <= r) {
                positions.push_back(i);
            }
        }

        return positions;

    }

    list_type getValuedPositionsInRange(size_type l, size_type r) override {

        list_type positions;
        for (size_type i = 0; i < length_; i++) {
            if (l <= positions_[i] && positions_[i] <= r) {
                positions.push_back(std::make_pair(i, values_[i]));
            }
        }

        return positions;

    }

    std::vector<elem_type> getAllElements() override {
        return std::vector<elem_type>(values_, values_ + length_);
    }

    std::vector<size_type> getAllPositions() override {
        return std::vector<size_type>(positions_, positions_ + length_);
    }

    list_type getAllValuedPositions() override {

        std::vector<std::pair<size_type, elem_type>> positions;
        for (size_type i = 0; i < length_; i++) {
            positions.push_back(std::make_pair(positions_[i], values_[i]));
        }

        return positions;

    }

    bool containsElement(size_type l, size_type r) override {

        bool flag = false;
        for (auto i = 0; i < length_ && !flag; i++) {
            flag = l <= positions_[i] && positions_[i] <= r;
        }

        return flag;

    }

    size_type countElements() override {
        return length_;
    }


    MiniRowTree* clone() const override {
        return new MiniRowTree<elem_type>(*this);
    }



    void print(bool all = false) override {

        std::cout << "### Parameters ###" << std::endl;
        std::cout << "length  = " << length_ << std::endl;
        std::cout << "null = " << null_ << std::endl;

        if (all) {

            std::cout << "### Positions ###" << std::endl;
            for (auto i = 0; i < length_; i++) std::cout << positions_[i];
            std::cout << std::endl << std::endl;

            std::cout << "### Values ###" << std::endl;
            for (auto i = 0; i < length_; i++) std::cout << values_[i];
            std::cout << std::endl << std::endl;

        }

    }


    void setNull(size_type i) override {

        auto iter = std::find(positions_, positions_ + length_, i);

        if (iter != positions_ + length_) {

            size_type* tmpPos = new size_type[length_ - 1];
            elem_type* tmpVal = new elem_type[length_ - 1];

            size_type pos = 0;
            size_type* p;
            elem_type* v;
            for (p = positions_, v = values_; p != positions_ + length_; p++, v++) {

                if (p != iter) {

                    tmpPos[pos] = *p;
                    tmpVal[pos] = *v;
                    pos++;

                }

            }

            delete positions_;
            delete values_;
            positions_ = tmpPos;
            values_ = tmpVal;
            length_--;

        }

    }


private:
    size_type* positions_;
    elem_type* values_;

    size_type length_;
    elem_type null_;

};

template<>
class MiniRowTree<bool> : public virtual RowTree<bool> {

public:
    typedef bool elem_type;

    typedef RelationList list_type;


    MiniRowTree() {
        // nothing to do
    }

    MiniRowTree(const MiniRowTree& other) {

        length_ = other.length_;

        positions_ = new size_type[length_];
        for (size_type i = 0; i < length_; i++) {
            positions_[i] = other.positions_[i];
        }

    }

    MiniRowTree& operator=(const MiniRowTree& other) {

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        delete positions_;

        length_ = other.length_;

        positions_ = new size_type[length_];
        for (size_type i = 0; i < length_; i++) {
            positions_[i] = other.positions_[i];
        }

        return *this;

    }

    MiniRowTree(const bit_vector_type& v) {

        length_ = 0;
        for (size_type i = 0; i < v.size(); i++) {
            length_ += v[i];
        }

        positions_ = new size_type[length_];

        for (size_type i = 0, pos = 0; i < v.size(); i++) {

            if (v[i]) {
                positions_[pos++] = i;
            }

        }

    }

    MiniRowTree(const list_type& list) {

        length_ = list.size();
        positions_ = new size_type[length_];

        for (size_type i = 0; i < list.size(); i++) {
            positions_[i] = list[i];
        }

    }

    MiniRowTree(const std::vector<std::pair<size_type, size_type>>::iterator& first, const std::vector<std::pair<size_type, size_type>>::iterator& last) {

        length_ = last - first;
        positions_ = new size_type[length_];

        size_type i = 0;
        for (auto iter = first; iter != last; iter++, i++) {
            positions_[i] = iter->second;
        }

    }

    ~MiniRowTree() {
        delete positions_;
    }

    size_type getLength() override {
        return length_;
    }

    elem_type getNull() override {
        return false;
    }


    bool isNotNull(size_type i) override {
        return std::find(positions_, positions_ + length_, i) != positions_ + length_;
    }

    elem_type getElement(size_type i) override {
        return std::find(positions_, positions_ + length_, i) != positions_ + length_;
    }

    size_type getFirst() override {

        size_type min = -1;
        for (auto i = 0; i < length_; i++) {
            min = std::min(min, positions_[i]);
        }

        return min;

    }

    std::vector<elem_type> getElementsInRange(size_type l, size_type r) override {

        std::vector<elem_type> elems;
        for (size_type i = 0; i < length_; i++) {
            if (l <= positions_[i] && positions_[i] <= r) {
                elems.push_back(true);
            }
        }

        return elems;

    }

    std::vector<size_type> getPositionsInRange(size_type l, size_type r) override {

        std::vector<size_type> positions;
        for (size_type i = 0; i < length_; i++) {
            if (l <= positions_[i] && positions_[i] <= r) {
                positions.push_back(i);
            }
        }

        return positions;

    }

    std::vector<std::pair<size_type, elem_type>> getValuedPositionsInRange(size_type l, size_type r) override {

        std::vector<std::pair<size_type, elem_type>> positions;
        for (size_type i = 0; i < length_; i++) {
            if (l <= positions_[i] && positions_[i] <= r) {
                positions.push_back(std::make_pair(i, true));
            }
        }

        return positions;

    }

    std::vector<elem_type> getAllElements() override {
        return std::vector<elem_type>(length_, true);
    }

    std::vector<size_type> getAllPositions() override {
        return std::vector<size_type>(positions_, positions_ + length_);
    }

    std::vector<std::pair<size_type, elem_type>> getAllValuedPositions() override {

        std::vector<std::pair<size_type, elem_type>> positions;
        for (size_type i = 0; i < length_; i++) {
            positions.push_back(std::make_pair(positions_[i], true));
        }

        return positions;

    }

    bool containsElement(size_type l, size_type r) override {

        bool flag = false;
        for (auto i = 0; i < length_ && !flag; i++) {
            flag = l <= positions_[i] && positions_[i] <= r;
        }

        return flag;

    }

    size_type countElements() override {
        return length_;
    }


    MiniRowTree* clone() const override {
        return new MiniRowTree<elem_type>(*this);
    }



    void print(bool all = false) override {

        std::cout << "### Parameters ###" << std::endl;
        std::cout << "length  = " << length_ << std::endl;
        std::cout << "null = " << false << std::endl;

        if (all) {

            std::cout << "### Values ###" << std::endl;
            for (auto i = 0; i < length_; i++) std::cout << positions_[i];
            std::cout << std::endl << std::endl;

        }

    }


    void setNull(size_type i) override {

        auto iter = std::find(positions_, positions_ + length_, i);

        if (iter != positions_ + length_) {

            size_type* tmp = new size_type[length_ - 1];

            size_type pos = 0;
            for (size_type* p = positions_; p != positions_ + length_; p++) {

                if (p != iter) {
                    tmp[pos++] = *p;
                }

            }

            delete positions_;
            positions_ = tmp;
            length_--;

        }

    }


private:
    size_type length_;
    size_type* positions_;

};

#endif //K2TREES_MINIBASICROWTREE_HPP
