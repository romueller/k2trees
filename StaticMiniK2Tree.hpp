#ifndef K2TREES_STATICMINIK2TREE_HPP
#define K2TREES_STATICMINIK2TREE_HPP

#include "K2Tree.hpp"
#include "Utility.hpp"

template<typename E>
class MiniK2Tree : public virtual K2Tree<E> {

public:
    typedef E elem_type;

    typedef typename K2Tree<elem_type>::matrix_type matrix_type;
    typedef typename K2Tree<elem_type>::list_type list_type;
    typedef typename K2Tree<elem_type>::positions_type positions_type;
    typedef typename K2Tree<elem_type>::pairs_type pairs_type;


    MiniK2Tree() {
        // nothing to do
    }

    MiniK2Tree(const MiniK2Tree& other) {

        null_ = other.null_;

        length_ = other.length_;
        positions_ = new std::pair<size_type, size_type>[length_];
        values_ = new elem_type[length_];
        for (size_type k = 0; k < length_; k++) {

            positions_[k] = other.positions_[k];
            values_[k] = other.values_[k];

        }

    }

    MiniK2Tree& operator=(const MiniK2Tree& other) {

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        delete[] positions_;
        delete[] values_;

        null_ = other.null_;

        length_ = other.length_;
        positions_ = new std::pair<size_type, size_type>[length_];
        values_ = new elem_type[length_];
        for (size_type k = 0; k < length_; k++) {

            positions_[k] = other.positions_[k];
            values_[k] = other.values_[k];

        }

        return *this;

    }

    // assumes that all rows of mat are equally long
    MiniK2Tree(const matrix_type& mat, const elem_type null = elem_type()) {

        null_ = null;

        length_ = 0;
        for (size_type i = 0; i < mat.size(); i++) {
            for (size_type j = 0; j < mat[i].size(); j++) {
                length_ += (mat[i][j] != null);
            }
        }

        positions_ = new std::pair<size_type, size_type>[length_];
        values_ = new elem_type[length_];

        size_type pos = 0;
        for (size_type i = 0; i < mat.size(); i++) {
            for (size_type j = 0; j < mat[i].size(); j++) {

                if (mat[i][j] != null) {

                    positions_[pos] = std::make_pair(i, j);
                    values_[pos++] = mat[i][j];

                }

            }
        }

    }

    MiniK2Tree(const std::vector<list_type>& lists, const elem_type null = elem_type()) {

        null_ = null;

        length_ = 0;
        for (size_type i = 0; i < lists.size(); i++) {
            length_ += lists[i].size();
        }

        positions_ = new std::pair<size_type, size_type>[length_];
        values_ = new elem_type[length_];

        size_type pos = 0;
        for (size_type i = 0; i < lists.size(); i++) {
            for (size_type j = 0; j < lists[i].size(); j++) {

                positions_[pos] = std::make_pair(i, lists[i][j].first);
                values_[pos++] = lists[i][j].second;

            }
        }

    }

    MiniK2Tree(pairs_type& pairs, const elem_type null = elem_type()) {

        null_ = null;

        length_ = pairs.size();
        positions_ = new std::pair<size_type, size_type>[length_];
        values_ = new elem_type[length_];

        size_type pos = 0;
        for (auto& p : pairs) {

            positions_[pos] = std::make_pair(p.row, p.col);
            values_[pos++] = p.val;

        }

    }

    MiniK2Tree(const typename pairs_type::iterator& first, const typename pairs_type::iterator& last, const size_type x, const size_type y, const elem_type null = elem_type()) {

        null_ = null;

        length_ = last - first;
        positions_ = new std::pair<size_type, size_type>[length_];
        values_ = new elem_type[length_];

        size_type pos = 0;
        for (auto iter = first; iter != last; iter++) {

            positions_[pos] = std::make_pair(iter->row - x, iter->col - y);
            values_[pos++] = iter->val;
            pos++;

        }

    }

    ~MiniK2Tree() {

        delete[] positions_;
        delete[] values_;

    }

    size_type getNumRows() override {

        size_type max = 0;
        for (size_type k = 0; k < length_; k++) {
            max = std::max(max, positions_[k].first + 1);
        }

        return max;

    }

    size_type getNumCols() override {

        size_type max = 0;
        for (size_type k = 0; k < length_; k++) {
            max = std::max(max, positions_[k].second + 1);
        }

        return max;

    }

    elem_type getNull() override {
        return null_;
    }


    bool isNotNull(size_type i, size_type j) override {
        return std::find_if(positions_, positions_ + length_,
                            [i, j](const std::pair<size_type, size_type>& val) {
                                return val.first == i && val.second == j;
                            }) != positions_ + length_;
    }

    elem_type getElement(size_type i, size_type j) override {

        auto iter = std::find_if(positions_, positions_ + length_,
                                 [i, j](const std::pair<size_type, size_type>& val) {
                                     return val.first == i && val.second == j;
                                 });

        return (iter != positions_ + length_) ? values_[iter - positions_] : null_;

    }

    std::vector<elem_type> getSuccessorElements(size_type i) override {

        std::vector<elem_type> succs;
        for (size_type k = 0; k < length_; k++) {
            if (positions_[k].first == i) {
                succs.push_back(values_[k]);
            }
        }

        return succs;

    }

    std::vector<size_type> getSuccessorPositions(size_type i) override {

        std::vector<size_type> succs;
        for (size_type k = 0; k < length_; k++) {
            if (positions_[k].first == i) {
                succs.push_back(positions_[k].second);
            }
        }

        return succs;

    }

    pairs_type getSuccessorValuedPositions(size_type i) override {

        pairs_type succs;
        for (size_type k = 0; k < length_; k++) {
            if (positions_[k].first == i) {
                succs.push_back(ValuedPosition<elem_type>(positions_[k], values_[k]));
            }
        }

        return succs;

    }

    std::vector<elem_type> getPredecessorElements(size_type j) override {

        std::vector<elem_type> preds;
        for (size_type k = 0; k < length_; k++) {
            if (positions_[k].second == j) {
                preds.push_back(values_[k]);
            }
        }

        return preds;

    }

    std::vector<size_type> getPredecessorPositions(size_type j) override {

        std::vector<size_type> preds;
        for (size_type k = 0; k < length_; k++) {
            if (positions_[k].second == j) {
                preds.push_back(positions_[k].first);
            }
        }

        return preds;

    }

    pairs_type getPredecessorValuedPositions(size_type j) override {

        pairs_type preds;
        for (size_type k = 0; k < length_; k++) {
            if (positions_[k].second == j) {
                preds.push_back(ValuedPosition<elem_type>(positions_[k], values_[k]));
            }
        }

        return preds;

    }

    std::vector<elem_type> getElementsInRange(size_type i1, size_type i2, size_type j1, size_type j2) override {

        std::vector<elem_type> elements;
        for (size_type k = 0; k < length_; k++) {
            if (i1 <= positions_[k].first && positions_[k].first <= i2 && j1 <= positions_[k].second && positions_[k].second <= j2) {
                elements.push_back(values_[k]);
            }
        }

        return elements;

    }

    positions_type getPositionsInRange(size_type i1, size_type i2, size_type j1, size_type j2) override {

        positions_type pairs;
        for (size_type k = 0; k < length_; k++) {
            if (i1 <= positions_[k].first && positions_[k].first <= i2 && j1 <= positions_[k].second && positions_[k].second <= j2) {
                pairs.push_back(positions_[k]);
            }
        }

        return pairs;

    }

    pairs_type getValuedPositionsInRange(size_type i1, size_type i2, size_type j1, size_type j2) override {

        pairs_type pairs;
        for (size_type k = 0; k < length_; k++) {
            if (i1 <= positions_[k].first && positions_[k].first <= i2 && j1 <= positions_[k].second && positions_[k].second <= j2) {
                pairs.push_back(ValuedPosition<elem_type>(positions_[k], values_[k]));
            }
        }

        return pairs;

    }

    std::vector<elem_type> getAllElements() override {
        return std::vector<elem_type>(values_, values_ + length_);
    }

    positions_type getAllPositions() override {
        return positions_type(positions_, positions_ + length_);
    }

    pairs_type getAllValuedPositions() override {

        pairs_type pairs;
        for (size_type k = 0; k < length_; k++) {
            pairs.push_back(ValuedPosition<elem_type>(positions_[k], values_[k]));
        }

        return pairs;

    }

    bool containsElement(size_type i1, size_type i2, size_type j1, size_type j2) override {

        bool flag = false;
        for (size_type k = 0; k < length_ && !flag; k++) {
            flag = i1 <= positions_[k].first && positions_[k].first <= i2 && j1 <= positions_[k].second && positions_[k].second <= j2;
        }

        return flag;

    }

    size_type countElements() override {
        return length_;
    }


    MiniK2Tree* clone() const override {
        return new MiniK2Tree<elem_type>(*this);
    }


    void print(bool all = false) override {

        std::cout << "### Parameters ###" << std::endl;
        std::cout << "numRows  = " << getNumRows() << std::endl;
        std::cout << "numCols  = " << getNumCols() << std::endl;
        std::cout << "null = " << null_ << std::endl;

        if (all) {

            std::cout << "### Positions & Values ###" << std::endl;
            for (size_type k = 0; k < length_; k++) {
                std::cout << "(" << positions_[k].first << ", " << positions_[k].second << ", " << values_[k] << ") ";
            }
            std::cout << std::endl << std::endl;

        }

    }


    // method aliases using "relation nomenclature"

    bool areRelated(size_type i, size_type j) override {
        return isNotNull(i, j);
    }

    std::vector<size_type> getSuccessors(size_type i) override {
        return getSuccessorPositions(i);
    }

    size_type getFirstSuccessor(size_type i) override {

        size_type min = getNumCols();
        for (size_type k = 0; k < length_; k++) {
            if (positions_[k].first == i) {
                min = std::min(min, positions_[k].second);
            }
        }

        return min;

    }

    std::vector<size_type> getPredecessors(size_type j) override {
        return getPredecessorPositions(j);
    }

    positions_type getRange(size_type i1, size_type i2, size_type j1, size_type j2) override {
        return getPositionsInRange(i1, i2, j1, j2);
    }

    bool containsLink(size_type i1, size_type i2, size_type j1, size_type j2) override {
        return containsElement(i1, i2, j1, j2);
    }

    size_type countLinks() override {
        return countElements();
    }


    // note: can "invalidate" the data structure (containsLink() probably won't work correctly afterwards)
    void setNull(size_type i, size_type j) {

        auto iter = std::find_if(positions_, positions_ + length_,
                                 [i, j](const std::pair<size_type, size_type>& val) {
                                     return val.first == i && val.second == j;
                                 });

        if (iter != positions_ + length_) {

            std::pair<size_type, size_type>* tmpPos = new std::pair<size_type, size_type>[length_ - 1];
            elem_type* tmpVal = new elem_type[length_ - 1];

            size_type pos = 0;
            std::pair<size_type, size_type>* p;
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
    std::pair<size_type, size_type>* positions_;
    elem_type* values_;

    size_type length_;
    elem_type null_;

};

template<>
class MiniK2Tree<bool> : public virtual K2Tree<bool> {

public:
    typedef bool elem_type;

    typedef K2Tree<elem_type>::matrix_type matrix_type;
    typedef RelationList list_type;
    typedef K2Tree<elem_type>::positions_type positions_type;
    typedef K2Tree<elem_type>::pairs_type pairs_type;


    MiniK2Tree() {
        // nothing to do
    }

    MiniK2Tree(const MiniK2Tree& other) {

        length_ = other.length_;
        positions_ = new std::pair<size_type, size_type>[length_];
        for (size_type k = 0; k < length_; k++) {
            positions_[k] = other.positions_[k];
        }

    }

    MiniK2Tree& operator=(const MiniK2Tree& other) {

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        delete[] positions_;

        length_ = other.length_;
        positions_ = new std::pair<size_type, size_type>[length_];
        for (size_type k = 0; k < length_; k++) {
            positions_[k] = other.positions_[k];
        }

        return *this;

    }

    // assumes that all rows of mat are equally long
    MiniK2Tree(const matrix_type& mat) {

        length_ = 0;
        for (size_type i = 0; i < mat.size(); i++) {
            for (size_type j = 0; j < mat[i].size(); j++) {
                length_ += mat[i][j];
            }
        }

        positions_ = new std::pair<size_type, size_type>[length_];

        size_type pos = 0;
        for (size_type i = 0; i < mat.size(); i++) {
            for (size_type j = 0; j < mat[i].size(); j++) {

                if (mat[i][j]) {
                    positions_[pos++] = std::make_pair(i, j);
                }

            }
        }

    }

    MiniK2Tree(const std::vector<list_type>& lists) {

        length_ = 0;
        for (size_type i = 0; i < lists.size(); i++) {
            length_ += lists[i].size();
        }

        positions_ = new std::pair<size_type, size_type>[length_];

        size_type pos = 0;
        for (size_type i = 0; i < lists.size(); i++) {
            for (size_type j = 0; j < lists[i].size(); j++) {
                positions_[pos++] = std::make_pair(i, lists[i][j]);
            }
        }

    }

    MiniK2Tree(positions_type& pairs) {

        length_ = pairs.size();
        positions_ = new std::pair<size_type, size_type>[length_];

        size_type pos = 0;
        for (auto& p : pairs) {
            positions_[pos++] = p;
        }

    }

    MiniK2Tree(const positions_type::iterator& first, const positions_type::iterator& last) {

        length_ = last - first;
        positions_ = new std::pair<size_type, size_type>[length_];

        size_type pos = 0;
        for (auto iter = first; iter != last; iter++) {
            positions_[pos++] = *iter;
        }

    }

    MiniK2Tree(const positions_type::iterator& first, const positions_type::iterator& last, const size_type x, const size_type y) {

        length_ = last - first;
        positions_ = new std::pair<size_type, size_type>[length_];

        size_type pos = 0;
        for (auto iter = first; iter != last; iter++) {

            positions_[pos] = *iter;
            positions_[pos].first -= x;
            positions_[pos].second -= y;
            pos++;

        }

    }

    ~MiniK2Tree() {
        delete[] positions_;
    }

    size_type getNumRows() override {

        size_type max = 0;
        for (size_type i = 0; i < length_; i++) {
            max = std::max(max, positions_[i].first + 1);
        }

        return max;


    }

    size_type getNumCols() override {

        size_type max = 0;
        for (size_type i = 0; i < length_; i++) {
            max = std::max(max, positions_[i].second + 1);
        }

        return max;

    }

    elem_type getNull() override {
        return false;
    }


    bool isNotNull(size_type i, size_type j) override {
        return std::find_if(positions_, positions_ + length_,
                            [i, j](const std::pair<size_type, size_type>& val) {
                                return val.first == i && val.second == j;
                            }) != positions_ + length_;
    }

    elem_type getElement(size_type i, size_type j) override {
        return std::find_if(positions_, positions_ + length_,
                            [i, j](const std::pair<size_type, size_type>& val) {
                                return val.first == i && val.second == j;
                            }) != positions_ + length_;
    }

    std::vector<elem_type> getSuccessorElements(size_type i) override {

        std::vector<elem_type> succs;
        for (size_type k = 0; k < length_; k++) {
            if (positions_[k].first == i) {
                succs.push_back(true);
            }
        }

        return succs;

    }

    std::vector<size_type> getSuccessorPositions(size_type i) override {

        std::vector<size_type> succs;
        for (size_type k = 0; k < length_; k++) {
            if (positions_[k].first == i) {
                succs.push_back(positions_[k].second);
            }
        }

        return succs;

    }

    pairs_type getSuccessorValuedPositions(size_type i) override {

        pairs_type succs;
        for (size_type k = 0; k < length_; k++) {
            if (positions_[k].first == i) {
                succs.push_back(ValuedPosition<elem_type>(positions_[k], true));
            }
        }

        return succs;

    }

    std::vector<elem_type> getPredecessorElements(size_type j) override {

        std::vector<elem_type> preds;
        for (size_type k = 0; k < length_; k++) {
            if (positions_[k].second == j) {
                preds.push_back(true);
            }
        }

        return preds;

    }

    std::vector<size_type> getPredecessorPositions(size_type j) override {

        std::vector<size_type> preds;
        for (size_type k = 0; k < length_; k++) {
            if (positions_[k].second == j) {
                preds.push_back(positions_[k].first);
            }
        }

        return preds;

    }

    pairs_type getPredecessorValuedPositions(size_type j) override {

        pairs_type preds;
        for (size_type k = 0; k < length_; k++) {
            if (positions_[k].second == j) {
                preds.push_back(ValuedPosition<elem_type>(positions_[k], true));
            }
        }

        return preds;

    }

    std::vector<elem_type> getElementsInRange(size_type i1, size_type i2, size_type j1, size_type j2) override {

        std::vector<elem_type> elements;
        for (size_type k = 0; k < length_; k++) {
            if (i1 <= positions_[k].first && positions_[k].first <= i2 && j1 <= positions_[k].second && positions_[k].second <= j2) {
                elements.push_back(true);
            }
        }

        return elements;

    }

    positions_type getPositionsInRange(size_type i1, size_type i2, size_type j1, size_type j2) override {

        positions_type pairs;
        for (size_type k = 0; k < length_; k++) {
            if (i1 <= positions_[k].first && positions_[k].first <= i2 && j1 <= positions_[k].second && positions_[k].second <= j2) {
                pairs.push_back(positions_[k]);
            }
        }

        return pairs;

    }

    pairs_type getValuedPositionsInRange(size_type i1, size_type i2, size_type j1, size_type j2) override {

        pairs_type pairs;
        for (size_type k = 0; k < length_; k++) {
            if (i1 <= positions_[k].first && positions_[k].first <= i2 && j1 <= positions_[k].second && positions_[k].second <= j2) {
                pairs.push_back(ValuedPosition<elem_type>(positions_[k], true));
            }
        }

        return pairs;

    }

    std::vector<elem_type> getAllElements() override {
        return std::vector<elem_type>(length_, true);
    }

    positions_type getAllPositions() override {
        return positions_type(positions_, positions_ + length_);
    }

    pairs_type getAllValuedPositions() override {

        pairs_type pairs;
        for (size_type k = 0; k < length_; k++) {
            pairs.push_back(ValuedPosition<elem_type>(positions_[k], true));
        }

        return pairs;

    }

    bool containsElement(size_type i1, size_type i2, size_type j1, size_type j2) override {

        bool flag = false;
        for (size_type k = 0; k < length_ && !flag; k++) {
            flag = i1 <= positions_[k].first && positions_[k].first <= i2 && j1 <= positions_[k].second && positions_[k].second <= j2;
        }

        return flag;

    }

    size_type countElements() override {
        return length_;
    }


    MiniK2Tree* clone() const override {
        return new MiniK2Tree<elem_type>(*this);
    }


    void print(bool all = false) override {

        std::cout << "### Parameters ###" << std::endl;
        std::cout << "numRows  = " << getNumRows() << std::endl;
        std::cout << "numCols  = " << getNumCols() << std::endl;
        std::cout << "null = " << false << std::endl;

        if (all) {

            std::cout << "### Positions & Values ###" << std::endl;
            for (size_type k = 0; k < length_; k++) {
                std::cout << "(" << positions_[k].first << ", " << positions_[k].second << ", " << true << ") ";
            }
            std::cout << std::endl << std::endl;

        }

    }


    // method aliases using "relation nomenclature"

    bool areRelated(size_type i, size_type j) override {
        return isNotNull(i, j);
    }

    std::vector<size_type> getSuccessors(size_type i) override {
        return getSuccessorPositions(i);
    }

    size_type getFirstSuccessor(size_type i) override {

        size_type min = getNumCols();
        for (size_type k = 0; k < length_; k++) {
            if (positions_[k].first == i) {
                min = std::min(min, positions_[k].second);
            }
        }

        return min;

    }

    std::vector<size_type> getPredecessors(size_type j) override {
        return getPredecessorPositions(j);
    }

    positions_type getRange(size_type i1, size_type i2, size_type j1, size_type j2) override {
        return getPositionsInRange(i1, i2, j1, j2);
    }

    bool containsLink(size_type i1, size_type i2, size_type j1, size_type j2) override {
        return containsElement(i1, i2, j1, j2);
    }

    size_type countLinks() override {
        return countElements();
    }


    // note: can "invalidate" the data structure (containsLink() probably won't work correctly afterwards)
    void setNull(size_type i, size_type j) {

        auto iter = std::find_if(positions_, positions_ + length_,
                                 [i, j](const std::pair<size_type, size_type>& val) {
                                     return val.first == i && val.second == j;
                                 });

        if (iter != positions_ + length_) {

            std::pair<size_type, size_type>* tmpPos = new std::pair<size_type, size_type>[length_ - 1];

            size_type pos = 0;
            std::pair<size_type, size_type>* p;
            for (p = positions_; p != positions_ + length_; p++) {

                if (p != iter) {

                    tmpPos[pos] = *p;
                    pos++;

                }

            }

            delete positions_;
            positions_ = tmpPos;
            length_--;

        }

    }


private:
    std::pair<size_type, size_type>* positions_;
    size_type length_;

};

#endif //K2TREES_STATICMINIK2TREE_HPP
