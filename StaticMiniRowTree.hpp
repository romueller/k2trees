/*
 * Copyright (C) 2017 Robert Mueller
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contact: Robert Mueller <romueller@techfak.uni-bielefeld.de>
 * Faculty of Technology, Bielefeld University,
 * PO box 100131, DE-33501 Bielefeld, Germany
 */

#ifndef K2TREES_MINIBASICROWTREE_HPP
#define K2TREES_MINIBASICROWTREE_HPP

#include <queue>

#include "RowTree.hpp"
#include "Utility.hpp"

/**
 * Naive implementation of a universe subset with a RowTree interface for very small subsets.
 *
 * Simply contains a list of the elements.
 */
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

    /**
     * Vector-based constructor
     */
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

    /**
     * List-of-pairs-based constructor
     */
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
                positions.push_back(positions_[i]);
            }
        }

        return positions;

    }

    list_type getValuedPositionsInRange(size_type l, size_type r) override {

        list_type positions;
        for (size_type i = 0; i < length_; i++) {
            if (l <= positions_[i] && positions_[i] <= r) {
                positions.push_back(std::make_pair(positions_[i], values_[i]));
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

    size_type getFirst() override {

        size_type min = -1;
        for (auto i = 0; i < length_; i++) {
            min = std::min(min, positions_[i]);
        }

        return min;

    }



private:
    size_type* positions_; // positions of all elements
    elem_type* values_; // values of all elements

    size_type length_; // number of elements
    elem_type null_; // null element

};


/**
 * Bool specialisation of MiniRowTree.
 *
 * Has the same characteristics as the general implementation above,
 * but makes use of some simplifications since the only non-null value is 1 / true.
 */
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

    /**
     * Vector-based constructor
     */
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

    /**
     * List-of-pairs-based constructor
     */
    MiniRowTree(const list_type& list) {

        length_ = list.size();
        positions_ = new size_type[length_];

        for (size_type i = 0; i < list.size(); i++) {
            positions_[i] = list[i];
        }

    }

    /**
     * List-of-pairs-based constructor similar to the one above, but only a part of the universe is used.
     */
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
                positions.push_back(positions_[i]);
            }
        }

        return positions;

    }

    std::vector<std::pair<size_type, elem_type>> getValuedPositionsInRange(size_type l, size_type r) override {

        std::vector<std::pair<size_type, elem_type>> positions;
        for (size_type i = 0; i < length_; i++) {
            if (l <= positions_[i] && positions_[i] <= r) {
                positions.push_back(std::make_pair(positions_[i], true));
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

    size_type getFirst() override {

        size_type min = -1;
        for (auto i = 0; i < length_; i++) {
            min = std::min(min, positions_[i]);
        }

        return min;

    }



private:
    size_type length_; // number of elements
    size_type* positions_; // positions of all elements

};

#endif //K2TREES_MINIBASICROWTREE_HPP
