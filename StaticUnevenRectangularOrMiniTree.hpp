#ifndef K2TREES_STATICUNEVENRECTANGULARORMINITREE_HPP
#define K2TREES_STATICUNEVENRECTANGULARORMINITREE_HPP


#include <queue>

#include "K2Tree.hpp"
#include "StaticBasicRectangularTree.hpp"
#include "StaticRowTree.hpp"
#include "Utility.hpp"
#include "StaticMiniK2Tree.hpp"

template<typename E>
class UnevenKrKcOrMiniTree : public virtual K2Tree<E> {

public:
    typedef E elem_type;

    typedef typename K2Tree<elem_type>::matrix_type matrix_type;
    typedef typename K2Tree<elem_type>::list_type list_type;
    typedef typename K2Tree<elem_type>::positions_type positions_type;
    typedef typename K2Tree<elem_type>::pairs_type pairs_type;


    UnevenKrKcOrMiniTree() {

        partitions_ = 0;
        numPartitions_ = 0;

    }

    UnevenKrKcOrMiniTree(const UnevenKrKcOrMiniTree& other) {

        hr_ = other.hr_;
        hc_ = other.hc_;
        kr_ = other.kr_;
        kc_ = other.kc_;
        numRows_ = other.numRows_;
        numCols_ = other.numCols_;
        null_ = other.null_;

        partitions_ = new K2Tree<elem_type>*[other.numPartitions_];
        for (size_type k = 0; k < other.numPartitions_; k++) {
            partitions_[k] = other.partitions_[k]->clone();
        }
        partitionSize_ = other.partitionSize_;
        numPartitions_ = other.numPartitions_;

    }

    UnevenKrKcOrMiniTree& operator=(const UnevenKrKcOrMiniTree& other) {

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        hr_ = other.hr_;
        hc_ = other.hc_;
        kr_ = other.kr_;
        kc_ = other.kc_;
        numRows_ = other.numRows_;
        numCols_ = other.numCols_;
        null_ = other.null_;

        for (size_type k = 0; k < numPartitions_; k++) {
            delete partitions_[k];
        }
        delete[] partitions_;
        partitions_ = new K2Tree<elem_type>*[other.numPartitions_];
        for (size_type k = 0; k < other.numPartitions_; k++) {
            partitions_[k] = other.partitions_[k]->clone();
        }
        partitionSize_ = other.partitionSize_;
        numPartitions_ = other.numPartitions_;

        return *this;

    }
#if 0
    // assumes that all rows of mat are equally long
    UnevenKrKcOrMiniTree(const matrix_type& mat, const size_type kr, const size_type kc, const elem_type null = elem_type()) {

        null_ = null;

        kr_ = kr;
        kc_ = kc;
        hr_ = std::max((size_type)1, logK(mat.size(), kr_));
        hc_ = std::max((size_type)1, logK(mat[0].size(), kc_));
        numRows_ = size_type(pow(kr_, hr_));
        numCols_ = size_type(pow(kc_, hc_));


        if (hc_ > hr_) {

            partitionSize_ = size_type(pow(kc_, hr_));
            numPartitions_ = numCols_ / partitionSize_;
            partitions_ = new K2Tree<elem_type>*[numPartitions_];

            for (size_type i = 0; i < numPartitions_; i++) {//TODO-mb count elements in submatrix, use correct constructor for MiniK2Tree
                partitions_[i] = new KrKcTree<elem_type>(mat, 0, i * partitionSize_, numRows_, partitionSize_, kr_, kc_, null);
            }


        } else {

            partitionSize_ = size_type(pow(kr_, hc_));
            numPartitions_ = numRows_ / partitionSize_;
            partitions_ = new K2Tree<elem_type>*[numPartitions_];

            for (size_type j = 0; j < numPartitions_; j++) {//TODO-mb count elements in submatrix, use correct constructor for MiniK2Tree
                partitions_[j] = new KrKcTree<elem_type>(mat, j * partitionSize_, 0, partitionSize_, numCols_, kr_, kc_, null);
            }

        }

#if 1
        for (size_type k = 0; k < numPartitions_; k++) {

            if (partitions_[k]->getNumRows() == 0) {

                delete partitions_[k];
                partitions_[k] = 0;

            }

        }
#endif

    }

    UnevenKrKcOrMiniTree(const std::vector<list_type>& lists, const size_type kr, const size_type kc, const int mode, const elem_type null = elem_type()) {

        null_ = null;

        size_type maxCol = 0;
        for (auto& row : lists) {
            for (auto& elem : row) {
                maxCol = std::max(maxCol, elem.first);
            }
        }
        maxCol++; // for number of columns

        kr_ = kr;
        kc_ = kc;
        hr_ = std::max((size_type)1, logK(lists.size(), kr_));
        hc_ = std::max((size_type)1, logK(maxCol, kc_));
        numRows_ = size_type(pow(kr_, hr_));
        numCols_ = size_type(pow(kc_, hc_));


        if (hc_ > hr_) {

            partitionSize_ = size_type(pow(kc_, hr_));
            numPartitions_ = numCols_ / partitionSize_;
            partitions_ = new K2Tree<elem_type>*[numPartitions_];

            for (size_type i = 0; i < numPartitions_; i++) {//TODO-mb count elements in submatrix, use correct constructor for MiniK2Tree
                partitions_[i] = new KrKcTree<elem_type>(lists, 0, i * partitionSize_, numRows_, partitionSize_, kr_, kc_, mode, null);
            }


        } else {

            partitionSize_ = size_type(pow(kr_, hc_));
            numPartitions_ = numRows_ / partitionSize_;
            partitions_ = new K2Tree<elem_type>*[numPartitions_];

            for (size_type j = 0; j < numPartitions_; j++) {//TODO-mb count elements in submatrix, use correct constructor for MiniK2Tree
                partitions_[j] = new KrKcTree<elem_type>(lists, j * partitionSize_, 0, partitionSize_, numCols_, kr_, kc_, mode, null);
            }

        }

#if 1
        for (size_type k = 0; k < numPartitions_; k++) {

            if (partitions_[k]->getNumRows() == 0) {

                delete partitions_[k];
                partitions_[k] = 0;

            }

        }
#endif

    }
#endif
    UnevenKrKcOrMiniTree(pairs_type& pairs, const size_type kr, const size_type kc, const size_type mb, const elem_type null = elem_type()) {

        null_ = null;

        size_type maxRow = 0;
        size_type maxCol = 0;
        for (auto p : pairs) {
            maxRow = std::max(maxRow, p.row);
            maxCol = std::max(maxCol, p.col);
        }
        maxRow++; // for number of rows
        maxCol++; // for number of columns

        kr_ = kr;
        kc_ = kc;
        hr_ = std::max((size_type)1, logK(maxRow, kr_));
        hc_ = std::max((size_type)1, logK(maxCol, kc_));
        numRows_ = size_type(pow(kr_, hr_));
        numCols_ = size_type(pow(kc_, hc_));

        if (hc_ > hr_) {

            partitionSize_ = size_type(pow(kc_, hr_));
            numPartitions_ = numCols_ / partitionSize_;
            partitions_ = new KrKcTree<elem_type>*[numPartitions_];

            Subproblem sp(0, numRows_ - 1, 0, numCols_ - 1, 0, pairs.size());
            std::vector<std::pair<size_type, size_type>> intervals(numPartitions_);
            countingSort(pairs, intervals, sp, numRows_, partitionSize_, numPartitions_);

            for (size_type i = 0; i < numPartitions_; i++) {
                if (intervals[i].second - intervals[i].first > mb) {
                    partitions_[i] = new KrKcTree<elem_type>(pairs, 0, i * partitionSize_, numRows_, partitionSize_, intervals[i].first, intervals[i].second, kr_, kc_, null);
                } else {
                    partitions_[i] = new MiniK2Tree<elem_type>(pairs.begin() + intervals[i].first, pairs.begin() + intervals[i].second, null);
                }
            }


        } else {

            partitionSize_ = size_type(pow(kr_, hc_));
            numPartitions_ = numRows_ / partitionSize_;
            partitions_ = new KrKcTree<elem_type>*[numPartitions_];

            Subproblem sp(0, numRows_ - 1, 0, numCols_ - 1, 0, pairs.size());
            std::vector<std::pair<size_type, size_type>> intervals(numPartitions_);
            countingSort(pairs, intervals, sp, partitionSize_, numCols_, numPartitions_);

            for (size_type j = 0; j < numPartitions_; j++) {
                if (intervals[j].second - intervals[j].first) {
                    partitions_[j] = new KrKcTree<elem_type>(pairs, j * partitionSize_, 0, partitionSize_, numCols_, intervals[j].first, intervals[j].second, kr_, kc_, null);
                } else {
                    partitions_[j] = new MiniK2Tree<elem_type>(pairs.begin() + intervals[j].first, pairs.begin() + intervals[j].second, null);
                }
            }

        }

#if 1
        for (size_type k = 0; k < numPartitions_; k++) {

            if (partitions_[k]->getNumRows() == 0) {

                delete partitions_[k];
                partitions_[k] = 0;

            }

        }
#endif

    }

    ~UnevenKrKcOrMiniTree() {

        for (size_type k = 0; k < numPartitions_; k++) {
            delete partitions_[k];
        }
        delete[] partitions_;

    }


    size_type getHr() {
        return hr_;
    }

    size_type getHc() {
        return hc_;
    }

    size_type getKr() {
        return kr_;
    }

    size_type getKc() {
        return kc_;
    }

    size_type getNumRows() override {
        return numRows_;
    }

    size_type getNumCols() override {
        return numCols_;
    }

    elem_type getNull() override {
        return null_;
    }


    bool isNotNull(size_type i, size_type j) override {

        auto pis = determineIndices(i, j);
        auto p = partitions_[pis.partition];

        return (p != 0) && p->isNotNull(pis.row, pis.col);

    }

    elem_type getElement(size_type i, size_type j) override {

        auto pis = determineIndices(i, j);
        auto p = partitions_[pis.partition];

        return (p != 0) ? p->getElement(pis.row, pis.col) : null_;

    }

    std::vector<elem_type> getSuccessorElements(size_type i) override {

        std::vector<elem_type> succs;

        if (hc_ > hr_) {

            for (size_type k = 0; k < numPartitions_; k++) {

                auto p = partitions_[k];
                if (p != 0) {

                    auto tmp = p->getSuccessorElements(i);
                    succs.reserve(succs.size() + tmp.size());
                    std::move(tmp.begin(), tmp.end(), std::back_inserter(succs));

                }

            }

        } else {

            auto pis = determineIndices(i, 0);
            auto p = partitions_[pis.partition];
            if (p != 0) {
                succs = p->getSuccessorElements(pis.row);
            }

        }

        return succs;

    }

    std::vector<size_type> getSuccessorPositions(size_type i) override {

        std::vector<size_type> succs;

        if (hc_ > hr_) {

            size_type offset = 0;
            for (size_type k = 0; k < numPartitions_; k++, offset += partitionSize_) {

                auto p = partitions_[k];
                if (p != 0) {

                    auto tmp = p->getSuccessorPositions(i);
                    for (size_type l = 0; l < tmp.size(); l++) {
                        tmp[l] += offset;
                    }

                    succs.reserve(succs.size() + tmp.size());
                    std::move(tmp.begin(), tmp.end(), std::back_inserter(succs));

                }

            }

        } else {

            auto pis = determineIndices(i, 0);
            auto p = partitions_[pis.partition];
            if (p != 0) {
                succs = p->getSuccessorPositions(pis.row);
            }

        }

        return succs;

    }

    pairs_type getSuccessorValuedPositions(size_type i) override {

        pairs_type succs;

        if (hc_ > hr_) {

            size_type offset = 0;
            for (size_type k = 0; k < numPartitions_; k++, offset += partitionSize_) {

                auto p = partitions_[k];
                if (p != 0) {

                    auto tmp = p->getSuccessorValuedPositions(i);
                    for (size_type l = 0; l < tmp.size(); l++) {
                        tmp[l].col += offset;
                    }

                    succs.reserve(succs.size() + tmp.size());
                    std::move(tmp.begin(), tmp.end(), std::back_inserter(succs));

                }

            }

        } else {

            auto pis = determineIndices(i, 0);
            auto p = partitions_[pis.partition];
            if (p != 0) {
                succs = p->getSuccessorValuedPositions(pis.row);
            }

            if (pis.partition > 0) {

                size_type offset = pis.partition * partitionSize_;
                for (auto& succ : succs) {
                    succ.row += offset;
                }

            }

        }

        return succs;

    }

    std::vector<elem_type> getPredecessorElements(size_type j) override {

        std::vector<elem_type> preds;

        if (hc_ < hr_) {

            for (size_type k = 0; k < numPartitions_; k++) {

                auto p = partitions_[k];
                if (p != 0) {

                    auto tmp = p->getPredecessorElements(j);
                    preds.reserve(preds.size() + tmp.size());
                    std::move(tmp.begin(), tmp.end(), std::back_inserter(preds));

                }

            }

        } else {

            auto pis = determineIndices(0, j);
            auto p = partitions_[pis.partition];
            if (p != 0) {
                preds = p->getPredecessorElements(pis.col);
            }

        }

        return preds;

    }

    std::vector<size_type> getPredecessorPositions(size_type j) override {

        std::vector<size_type> preds;

        if (hc_ < hr_) {

            size_type offset = 0;
            for (size_type k = 0; k < numPartitions_; k++, offset += partitionSize_) {

                auto p = partitions_[k];
                if (p != 0) {

                    auto tmp = p->getPredecessorPositions(j);
                    for (size_type l = 0; l < tmp.size(); l++) {
                        tmp[l] += offset;
                    }

                    preds.reserve(preds.size() + tmp.size());
                    std::move(tmp.begin(), tmp.end(), std::back_inserter(preds));

                }

            }

        } else {

            auto pis = determineIndices(0, j);
            auto p = partitions_[pis.partition];
            if (p != 0) {
                preds = p->getPredecessorPositions(pis.col);
            }

        }

        return preds;

    }

    pairs_type getPredecessorValuedPositions(size_type j) override {

        pairs_type preds;

        if (hc_ < hr_) {

            size_type offset = 0;
            for (size_type k = 0; k < numPartitions_; k++, offset += partitionSize_) {

                auto p = partitions_[k];
                if (p != 0) {

                    auto tmp = p->getPredecessorValuedPositions(j);
                    for (size_type l = 0; l < tmp.size(); l++) {
                        tmp[l].row += offset;
                    }

                    preds.reserve(preds.size() + tmp.size());
                    std::move(tmp.begin(), tmp.end(), std::back_inserter(preds));

                }

            }

        } else {

            auto pis = determineIndices(0, j);
            auto p = partitions_[pis.partition];
            if (p != 0) {
                preds = p->getPredecessorValuedPositions(pis.col);
            }

            if (pis.partition > 0) {

                size_type offset = pis.partition * partitionSize_;
                for (auto& pred : preds) {
                    pred.col += offset;
                }

            }

        }

        return preds;

    }

    std::vector<elem_type> getElementsInRange(size_type i1, size_type i2, size_type j1, size_type j2) override {

        std::vector<elem_type> elements;

        auto upperLeft = determineIndices(i1, j1);
        auto lowerRight = determineIndices(i2, j2);

        // range falls completely within one partition
        if (upperLeft.partition == lowerRight.partition) {

            auto p = partitions_[upperLeft.partition];
            if (p != 0) {
                elements = p->getElementsInRange(upperLeft.row, lowerRight.row, upperLeft.col, lowerRight.col);
            }

            return elements;

        }

        // range spans multiple partitions
        std::vector<std::vector<elem_type>> tmp;
        tmp.reserve(lowerRight.partition - upperLeft.partition + 1);

        if (hc_ > hr_) {

            // first partition (partially spanned)
            auto p = partitions_[upperLeft.partition];
            if (p != 0) {
                tmp.push_back(p->getElementsInRange(upperLeft.row, lowerRight.row, upperLeft.col, partitionSize_ - 1));
            }

            // intermediate partition (fully spanned, if any)
            for (size_type k = upperLeft.partition + 1; k < lowerRight.partition; k++) {

                p = partitions_[k];
                if (p != 0) {
                    tmp.push_back(p->getElementsInRange(upperLeft.row, lowerRight.row, 0, partitionSize_ - 1));
                }

            }

            // last partition (partially spanned)
            p = partitions_[lowerRight.partition];
            if (p != 0) {
                tmp.push_back(p->getElementsInRange(upperLeft.row, lowerRight.row, 0, lowerRight.col));
            }

            size_type toReserve = elements.size();
            for (auto& v : tmp) {
                toReserve += v.size();
            }
            elements.reserve(toReserve);

            for (auto& v : tmp) {

                std::move(v.begin(), v.end(), std::back_inserter(elements));
                v.clear();
                v.shrink_to_fit();

            }

        } else {

            // first partition (partially spanned)
            auto p = partitions_[upperLeft.partition];
            if (p != 0) {
                tmp.push_back(p->getElementsInRange(upperLeft.row, partitionSize_ - 1, upperLeft.col, lowerRight.col));
            }

            // intermediate partition (fully spanned, if any)
            for (size_type k = upperLeft.partition + 1; k < lowerRight.partition; k++) {

                p = partitions_[k];
                if (p != 0) {
                    tmp.push_back(p->getElementsInRange(0, partitionSize_ - 1, upperLeft.col, lowerRight.col));
                }

            }

            // last partition (partially spanned)
            p = partitions_[lowerRight.partition];
            if (p != 0) {
                tmp.push_back(p->getElementsInRange(0, lowerRight.row, upperLeft.col, lowerRight.col));
            }

        }

        // flatten results
        size_type toReserve = elements.size();
        for (auto& v : tmp) {
            toReserve += v.size();
        }
        elements.reserve(toReserve);

        for (auto& v : tmp) {

            std::move(v.begin(), v.end(), std::back_inserter(elements));
            v.clear();
            v.shrink_to_fit();

        }

        return elements;

    }

    positions_type getPositionsInRange(size_type i1, size_type i2, size_type j1, size_type j2) override {

        positions_type elements;

        auto upperLeft = determineIndices(i1, j1);
        auto lowerRight = determineIndices(i2, j2);

        // range falls completely within one partition
        if (upperLeft.partition == lowerRight.partition) {

            auto p = partitions_[upperLeft.partition];
            if (p != 0) {

                elements = p->getPositionsInRange(upperLeft.row, lowerRight.row, upperLeft.col, lowerRight.col);

                size_type offset = partitionSize_ * upperLeft.partition;

                if (hc_ > hr_) {

                    for (auto& e : elements) {
                        e.second += offset;
                    }

                } else {

                    for (auto& e : elements) {
                        e.first += offset;
                    }

                }

            }

            return elements;

        }

        // range spans multiple partitions
        std::vector<positions_type> tmp;
        tmp.reserve(lowerRight.partition - upperLeft.partition + 1);
        size_type offset = partitionSize_ * upperLeft.partition;

        if (hc_ > hr_) {

            // first partition (partially spanned)
            auto p = partitions_[upperLeft.partition];
            if (p != 0) {

                tmp.push_back(p->getPositionsInRange(upperLeft.row, lowerRight.row, upperLeft.col, partitionSize_ - 1));
                for (auto& e : tmp.back()) {
                    e.second += offset;
                }

            }

            // intermediate partition (fully spanned, if any)
            offset += partitionSize_;
            for (size_type k = upperLeft.partition + 1; k < lowerRight.partition; k++, offset += partitionSize_) {

                p = partitions_[k];
                if (p != 0) {

                    tmp.push_back(p->getPositionsInRange(upperLeft.row, lowerRight.row, 0, partitionSize_ - 1));
                    for (auto& e : tmp.back()) {
                        e.second += offset;
                    }

                }

            }

            // last partition (partially spanned)
            p = partitions_[lowerRight.partition];
            if (p != 0) {

                tmp.push_back(p->getPositionsInRange(upperLeft.row, lowerRight.row, 0, lowerRight.col));
                for (auto& e : tmp.back()) {
                    e.second += offset;
                }

            }

        } else {

            // first partition (partially spanned)
            auto p = partitions_[upperLeft.partition];
            if (p != 0) {

                tmp.push_back(p->getPositionsInRange(upperLeft.row, partitionSize_ - 1, upperLeft.col, lowerRight.col));
                for (auto& e : tmp.back()) {
                    e.first += offset;
                }

            }

            // intermediate partition (fully spanned, if any)
            offset += partitionSize_;
            for (size_type k = upperLeft.partition + 1; k < lowerRight.partition; k++, offset += partitionSize_) {

                p = partitions_[k];
                if (p != 0) {

                    tmp.push_back(p->getPositionsInRange(0, partitionSize_ - 1, upperLeft.col, lowerRight.col));
                    for (auto& e : tmp.back()) {
                        e.first += offset;
                    }

                }

            }

            // last partition (partially spanned)
            p = partitions_[lowerRight.partition];
            if (p != 0) {

                tmp.push_back(p->getPositionsInRange(0, lowerRight.row, upperLeft.col, lowerRight.col));
                for (auto& e : tmp.back()) {
                    e.first += offset;
                }

            }

        }

        // flatten results
        size_type toReserve = elements.size();
        for (auto& v : tmp) {
            toReserve += v.size();
        }
        elements.reserve(toReserve);

        for (auto& v : tmp) {

            std::move(v.begin(), v.end(), std::back_inserter(elements));
            v.clear();
            v.shrink_to_fit();

        }

        return elements;

    }

    pairs_type getValuedPositionsInRange(size_type i1, size_type i2, size_type j1, size_type j2) override {

        pairs_type elements;

        auto upperLeft = determineIndices(i1, j1);
        auto lowerRight = determineIndices(i2, j2);

        // range falls completely within one partition
        if (upperLeft.partition == lowerRight.partition) {

            auto p = partitions_[upperLeft.partition];
            if (p != 0) {

                elements = p->getValuedPositionsInRange(upperLeft.row, lowerRight.row, upperLeft.col, lowerRight.col);

                size_type offset = partitionSize_ * upperLeft.partition;

                if (hc_ > hr_) {

                    for (auto& e : elements) {
                        e.col += offset;
                    }

                } else {

                    for (auto& e : elements) {
                        e.row += offset;
                    }

                }

            }

            return elements;

        }

        // range spans multiple partitions
        std::vector<pairs_type> tmp;
        tmp.reserve(lowerRight.partition - upperLeft.partition + 1);
        size_type offset = partitionSize_ * upperLeft.partition;

        if (hc_ > hr_) {

            // first partition (partially spanned)
            auto p = partitions_[upperLeft.partition];
            if (p != 0) {

                tmp.push_back(p->getValuedPositionsInRange(upperLeft.row, lowerRight.row, upperLeft.col, partitionSize_ - 1));
                for (auto& e : tmp.back()) {
                    e.col += offset;
                }

            }

            // intermediate partition (fully spanned, if any)
            offset += partitionSize_;
            for (size_type k = upperLeft.partition + 1; k < lowerRight.partition; k++, offset += partitionSize_) {

                p = partitions_[k];
                if (p != 0) {

                    tmp.push_back(p->getValuedPositionsInRange(upperLeft.row, lowerRight.row, 0, partitionSize_ - 1));
                    for (auto& e : tmp.back()) {
                        e.col += offset;
                    }

                }

            }

            // last partition (partially spanned)
            p = partitions_[lowerRight.partition];
            if (p != 0) {

                tmp.push_back(p->getValuedPositionsInRange(upperLeft.row, lowerRight.row, 0, lowerRight.col));
                for (auto& e : tmp.back()) {
                    e.col += offset;
                }

            }

        } else {

            // first partition (partially spanned)
            auto p = partitions_[upperLeft.partition];
            if (p != 0) {

                tmp.push_back(p->getValuedPositionsInRange(upperLeft.row, partitionSize_ - 1, upperLeft.col, lowerRight.col));
                for (auto& e : tmp.back()) {
                    e.row += offset;
                }

            }

            // intermediate partition (fully spanned, if any)
            offset += partitionSize_;
            for (size_type k = upperLeft.partition + 1; k < lowerRight.partition; k++, offset += partitionSize_) {

                p = partitions_[k];
                if (p != 0) {

                    tmp.push_back(p->getValuedPositionsInRange(0, partitionSize_ - 1, upperLeft.col, lowerRight.col));
                    for (auto& e : tmp.back()) {
                        e.row += offset;
                    }

                }

            }

            // last partition (partially spanned)
            p = partitions_[lowerRight.partition];
            if (p != 0) {

                tmp.push_back(p->getValuedPositionsInRange(0, lowerRight.row, upperLeft.col, lowerRight.col));
                for (auto& e : tmp.back()) {
                    e.row += offset;
                }

            }

        }

        // flatten results
        size_type toReserve = elements.size();
        for (auto& v : tmp) {
            toReserve += v.size();
        }
        elements.reserve(toReserve);

        for (auto& v : tmp) {

            std::move(v.begin(), v.end(), std::back_inserter(elements));
            v.clear();
            v.shrink_to_fit();

        }

        return elements;

    }

    std::vector<elem_type> getAllElements() override {

        std::vector<elem_type> elements;

        for (size_type k = 0; k < numPartitions_; k++) {

            auto p = partitions_[k];
            if (p != 0) {

                auto tmp = p->getAllElements();
                elements.reserve(elements.size() + tmp.size());
                std::move(tmp.begin(), tmp.end(), std::back_inserter(elements));

            }

        }

        return elements;

    }

    positions_type getAllPositions() override {

        positions_type elements;
        size_type offset = 0;

        for (size_type k = 0; k < numPartitions_; k++, offset += partitionSize_) {

            auto p = partitions_[k];
            if (p != 0) {

                auto tmp = p->getAllPositions();
                if (hc_ > hr_) {

                    for (size_type l = 0; l < tmp.size(); l++) {
                        tmp[l].second += offset;
                    }

                } else {

                    for (size_type l = 0; l < tmp.size(); l++) {
                        tmp[l].first += offset;
                    }

                }

                elements.reserve(elements.size() + tmp.size());
                std::move(tmp.begin(), tmp.end(), std::back_inserter(elements));

            }

        }

        return elements;

    }

    pairs_type getAllValuedPositions() override {

        pairs_type elements;
        size_type offset = 0;

        for (size_type k = 0; k < numPartitions_; k++, offset += partitionSize_) {

            auto p = partitions_[k];
            if (p != 0) {

                auto tmp = p->getAllValuedPositions();
                if (hc_ > hr_) {

                    for (size_type l = 0; l < tmp.size(); l++) {
                        tmp[l].col += offset;
                    }

                } else {

                    for (size_type l = 0; l < tmp.size(); l++) {
                        tmp[l].row += offset;
                    }

                }

                elements.reserve(elements.size() + tmp.size());
                std::move(tmp.begin(), tmp.end(), std::back_inserter(elements));

            }

        }

        return elements;

    }

    bool containsElement(size_type i1, size_type i2, size_type j1, size_type j2) override {

        auto upperLeft = determineIndices(i1, j1);
        auto lowerRight = determineIndices(i2, j2);
//        auto lowerRight = determineIndices(std::min(i2, numRows_ - 1), std::min(j2, numCols_ - 1));

        // range falls completely within one partition
        if (upperLeft.partition == lowerRight.partition) {

            auto p = partitions_[upperLeft.partition];
            return (p != 0) && p->containsElement(upperLeft.row, lowerRight.row, upperLeft.col, lowerRight.col);

        }

        // range spans multiple partitions
        bool found;

        if (hc_ > hr_) {

            // first partition (partially spanned)
            auto p = partitions_[upperLeft.partition];
            found = (p != 0) && p->containsElement(upperLeft.row, lowerRight.row, upperLeft.col, partitionSize_ - 1);

            // intermediate partition (fully spanned, if any)
            for (size_type k = upperLeft.partition + 1; (k < lowerRight.partition) && !found; k++) {

                p = partitions_[k];
                found = (p != 0) && p->containsElement(upperLeft.row, lowerRight.row, 0, partitionSize_ - 1);

            }

            // last partition (partially spanned)
            p = partitions_[lowerRight.partition];
            found = found || ((p != 0) && p->containsElement(upperLeft.row, lowerRight.row, 0, lowerRight.col));

        } else {

            // first partition (partially spanned)
            auto p = partitions_[upperLeft.partition];
            found =  (p != 0) && p->containsElement(upperLeft.row, partitionSize_ - 1, upperLeft.col, lowerRight.col);

            // intermediate partition (fully spanned, if any)
            for (size_type k = upperLeft.partition + 1; (k < lowerRight.partition) && !found; k++) {

                p = partitions_[k];
                found = (p != 0) && p->containsElement(0, partitionSize_ - 1, upperLeft.col, lowerRight.col);

            }

            // last partition (partially spanned)
            p = partitions_[lowerRight.partition];
            found = found ||  ((p != 0) && p->containsElement(0, lowerRight.row, upperLeft.col, lowerRight.col));

        }

        return found;

    }

    size_type countElements() override {

        size_type cnt = 0;
        for (size_type k = 0; k < numPartitions_; k++) {

            auto p = partitions_[k];
            cnt += (p != 0) ? p->countElements() : 0;

        }

        return cnt;

    }


    UnevenKrKcOrMiniTree* clone() const override {
        return new UnevenKrKcOrMiniTree<elem_type>(*this);
    }


    void print(bool all = false) override {

        std::cout << "### Parameters ###" << std::endl;
        std::cout << "hr  = " << hr_ << std::endl;
        std::cout << "hc  = " << hc_ << std::endl;
        std::cout << "kr  = " << kr_ << std::endl;
        std::cout << "kc  = " << kc_ << std::endl;
        std::cout << "numRows = " << numRows_ << std::endl;
        std::cout << "numCols = " << numCols_ << std::endl;
        std::cout << "partitionSize = " << partitionSize_ << std::endl;
        std::cout << "numPartitions = " << numPartitions_ << std::endl;
        std::cout << "null = " << null_ << std::endl;

        if (all) {

            for (size_type k = 0; k < numPartitions_; k++) {

                std::cout << "===== Partition " << k << " =====" << std::endl;
                auto p = partitions_[k];
                if (p != 0) {
                    p->print(true);
                } else {
                    std::cout << "((ALL NULL))" << std::endl;
                }

            }

        }

    }


    // method aliases using "relation nomenclature"

    bool areRelated(size_type i, size_type j) override {
        return isNotNull(i, j);
    }

    std::vector<size_type> getSuccessors(size_type i) override {
        return getSuccessorPositions(i);
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
    void setNull(size_type i, size_type j) override {

        auto pis = determineIndices(i, j);
        auto p = partitions_[pis.partition];

        if (p != 0) {
            p->setNull(pis.row, pis.col);
        }

    }


private:
    size_type hr_;
    size_type hc_;
    size_type kr_;
    size_type kc_;
    size_type numRows_;
    size_type numCols_;

    K2Tree<elem_type>** partitions_;
    size_type partitionSize_;
    size_type numPartitions_;

    elem_type null_;


    /* helper methods for mapping (overall) indices to positions in the partitions */

    size_type determinePartition(size_type i, size_type j) {
        return (hc_ > hr_) ? (j / partitionSize_) : (i / partitionSize_);
    }

    PartitionIndices determineIndices(size_type i, size_type j) {
        return (hc_ > hr_) ? PartitionIndices(j / partitionSize_, i, j % partitionSize_) : PartitionIndices(i / partitionSize_, i % partitionSize_, j);
    }

    /* helper methods for inplace construction from single list of pairs */

    size_type computeKey(const typename pairs_type::value_type& pair, const Subproblem& sp, size_type widthRow, size_type widthCol) {
        return ((pair.row - sp.firstRow) / widthRow) + (pair.col - sp.firstCol) / widthCol;
    }

    void countingSort(pairs_type& pairs, std::vector<std::pair<size_type, size_type>>& intervals, const Subproblem& sp, size_type widthRow, size_type widthCol, size_type sup) {

        std::vector<size_type> counts(sup);

        // determine key frequencies
        for (auto i = sp.left; i < sp.right; i++) {
            counts[computeKey(pairs[i], sp, widthRow, widthCol)]++;
        }

        // determine starting index for each key
        size_type total = 0;
        size_type tmp;

        for (auto key = 0; key < sup; key++) {

            tmp = counts[key];
            counts[key] = total;
            total += tmp;

            intervals[key].first = counts[key];
            intervals[key].second = total;

        }

        // reorder pairs of current subproblem
        pairs_type tmpPairs(sp.right - sp.left + 1);
        for (auto i = sp.left; i < sp.right; i++) {

            tmpPairs[counts[computeKey(pairs[i], sp, widthRow, widthCol)]] = pairs[i];
            counts[computeKey(pairs[i], sp, widthRow, widthCol)]++;

        }

        for (auto i = sp.left; i < sp.right; i++) {
            pairs[i] = tmpPairs[i - sp.left];
        }

    }

};


template<>
class UnevenKrKcOrMiniTree<bool> : public virtual K2Tree<bool> {// can handle non-quadratic relation matrices and width / heights that are no power of k

public:
    typedef bool elem_type;

    typedef typename K2Tree<elem_type>::matrix_type matrix_type;
    typedef typename K2Tree<elem_type>::list_type list_type;
    typedef typename K2Tree<elem_type>::positions_type positions_type;
    typedef typename K2Tree<elem_type>::pairs_type pairs_type;


    UnevenKrKcOrMiniTree() {

        partitions_ = 0;
        numPartitions_ = 0;

    }

    UnevenKrKcOrMiniTree(const UnevenKrKcOrMiniTree& other) {

        hr_ = other.hr_;
        hc_ = other.hc_;
        kr_ = other.kr_;
        kc_ = other.kc_;
        numRows_ = other.numRows_;
        numCols_ = other.numCols_;
        null_ = other.null_;

        partitions_ = new K2Tree<elem_type>*[other.numPartitions_];
        for (size_type k = 0; k < other.numPartitions_; k++) {
            partitions_[k] = (other.partitions_[k] != 0) ? other.partitions_[k]->clone() : 0;
        }
        partitionSize_ = other.partitionSize_;
        numPartitions_ = other.numPartitions_;

    }

    UnevenKrKcOrMiniTree& operator=(const UnevenKrKcOrMiniTree& other) {

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        hr_ = other.hr_;
        hc_ = other.hc_;
        kr_ = other.kr_;
        kc_ = other.kc_;
        numRows_ = other.numRows_;
        numCols_ = other.numCols_;
        null_ = other.null_;

        for (size_type k = 0; k < numPartitions_; k++) {
            delete partitions_[k];
        }
        delete[] partitions_;
        partitions_ = new K2Tree<elem_type>*[other.numPartitions_];
        for (size_type k = 0; k < other.numPartitions_; k++) {
            partitions_[k] = (other.partitions_[k] != 0) ? other.partitions_[k]->clone() : 0;
        }
        partitionSize_ = other.partitionSize_;
        numPartitions_ = other.numPartitions_;

        return *this;

    }
#if 0
    // assumes that all rows of mat are equally long
    UnevenKrKcOrMiniTree(const matrix_type& mat, const size_type kr, const size_type kc, const size_type mb) {

        null_ = false;

        kr_ = kr;
        kc_ = kc;
        hr_ = std::max((size_type)1, logK(mat.size(), kr_));
        hc_ = std::max((size_type)1, logK(mat[0].size(), kc_));
        numRows_ = size_type(pow(kr_, hr_));
        numCols_ = size_type(pow(kc_, hc_));


        if (hc_ > hr_) {

            partitionSize_ = size_type(pow(kc_, hr_));
            numPartitions_ = numCols_ / partitionSize_;
            partitions_ = new K2Tree<elem_type>*[numPartitions_];

            for (size_type i = 0; i < numPartitions_; i++) {//TODO-mb count elements in submatrix, use correct constructor for MiniK2Tree
                partitions_[i] = new KrKcTree<elem_type>(mat, 0, i * partitionSize_, numRows_, partitionSize_, kr_, kc_);
            }


        } else {

            partitionSize_ = size_type(pow(kr_, hc_));
            numPartitions_ = numRows_ / partitionSize_;
            partitions_ = new K2Tree<elem_type>*[numPartitions_];

            for (size_type j = 0; j < numPartitions_; j++) {//TODO-mb count elements in submatrix, use correct constructor for MiniK2Tree
                partitions_[j] = new KrKcTree<elem_type>(mat, j * partitionSize_, 0, partitionSize_, numCols_, kr_, kc_);
            }

        }

#if 1
        for (size_type k = 0; k < numPartitions_; k++) {

            if (partitions_[k]->getNumRows() == 0) {

                delete partitions_[k];
                partitions_[k] = 0;

            }

        }
#endif

    }

    UnevenKrKcOrMiniTree(const RelationLists& lists, const size_type kr, const size_type kc, const int mode) {

        null_ = false;

        size_type maxCol = 0;
        for (auto& row : lists) {
            for (auto& elem : row) {
                maxCol = std::max(maxCol, elem);
            }
        }
        maxCol++; // for number of columns

        kr_ = kr;
        kc_ = kc;
        hr_ = std::max((size_type)1, logK(lists.size(), kr_));
        hc_ = std::max((size_type)1, logK(maxCol, kc_));
        numRows_ = size_type(pow(kr_, hr_));
        numCols_ = size_type(pow(kc_, hc_));


        if (hc_ > hr_) {

            partitionSize_ = size_type(pow(kc_, hr_));
            numPartitions_ = numCols_ / partitionSize_;
            partitions_ = new K2Tree<elem_type>*[numPartitions_];

            for (size_type i = 0; i < numPartitions_; i++) {//TODO-mb count elements in submatrix, use correct constructor for MiniK2Tree
                partitions_[i] = new KrKcTree<elem_type>(lists, 0, i * partitionSize_, numRows_, partitionSize_, kr_, kc_, mode);
            }


        } else {

            partitionSize_ = size_type(pow(kr_, hc_));
            numPartitions_ = numRows_ / partitionSize_;
            partitions_ = new K2Tree<elem_type>*[numPartitions_];

            for (size_type j = 0; j < numPartitions_; j++) {//TODO-mb count elements in submatrix, use correct constructor for MiniK2Tree
                partitions_[j] = new KrKcTree<elem_type>(lists, j * partitionSize_, 0, partitionSize_, numCols_, kr_, kc_, mode);
            }

        }

#if 1
        for (size_type k = 0; k < numPartitions_; k++) {

            if (partitions_[k]->getNumRows() == 0) {

                delete partitions_[k];
                partitions_[k] = 0;

            }

        }
#endif

    }
#endif
    UnevenKrKcOrMiniTree(positions_type& pairs, const size_type kr, const size_type kc, const size_type mb) {

        null_ = false;

        size_type maxRow = 0;
        size_type maxCol = 0;
        for (auto p : pairs) {
            maxRow = std::max(maxRow, p.first);
            maxCol = std::max(maxCol, p.second);
        }
        maxRow++; // for number of rows
        maxCol++; // for number of columns

        kr_ = kr;
        kc_ = kc;
        hr_ = std::max((size_type)1, logK(maxRow, kr_));
        hc_ = std::max((size_type)1, logK(maxCol, kc_));
        numRows_ = size_type(pow(kr_, hr_));
        numCols_ = size_type(pow(kc_, hc_));

        if (hc_ > hr_) {

            partitionSize_ = size_type(pow(kc_, hr_));
            numPartitions_ = numCols_ / partitionSize_;
            partitions_ = new K2Tree<elem_type>*[numPartitions_];

            Subproblem sp(0, numRows_ - 1, 0, numCols_ - 1, 0, pairs.size());
            std::vector<std::pair<size_type, size_type>> intervals(numPartitions_);
            countingSort(pairs, intervals, sp, numRows_, partitionSize_, numPartitions_);

            for (size_type i = 0; i < numPartitions_; i++) {
                if (intervals[i].second - intervals[i].first > mb) {
                    partitions_[i] = new KrKcTree<elem_type>(pairs, 0, i * partitionSize_, numRows_, partitionSize_, intervals[i].first, intervals[i].second, kr_, kc_);
                } else {
                    partitions_[i] = new MiniK2Tree<elem_type>(pairs.begin() + intervals[i].first, pairs.begin() + intervals[i].second, 0, i * partitionSize_);
                }
            }


        } else {

            partitionSize_ = size_type(pow(kr_, hc_));
            numPartitions_ = numRows_ / partitionSize_;
            partitions_ = new K2Tree<elem_type>*[numPartitions_];

            Subproblem sp(0, numRows_ - 1, 0, numCols_ - 1, 0, pairs.size());
            std::vector<std::pair<size_type, size_type>> intervals(numPartitions_);
            countingSort(pairs, intervals, sp, partitionSize_, numCols_, numPartitions_);

            for (size_type j = 0; j < numPartitions_; j++) {
                if (intervals[j].second - intervals[j].first > mb) {
                    partitions_[j] = new KrKcTree<elem_type>(pairs, j * partitionSize_, 0, partitionSize_, numCols_, intervals[j].first, intervals[j].second, kr_, kc_);
                } else {
                    partitions_[j] = new MiniK2Tree<elem_type>(pairs.begin() + intervals[j].first, pairs.begin() + intervals[j].second, j * partitionSize_, 0);
                }
            }

        }

#if 1
        for (size_type k = 0; k < numPartitions_; k++) {

            if (partitions_[k]->getNumRows() == 0) {

                delete partitions_[k];
                partitions_[k] = 0;

            }

        }
#endif

    }

    ~UnevenKrKcOrMiniTree() {

        for (size_type k = 0; k < numPartitions_; k++) {
            delete partitions_[k];
        }
        delete[] partitions_;

    }


    size_type getHr() {
        return hr_;
    }

    size_type getHc() {
        return hc_;
    }

    size_type getKr() {
        return kr_;
    }

    size_type getKc() {
        return kc_;
    }

    size_type getNumRows() override {
        return numRows_;
    }

    size_type getNumCols() override {
        return numCols_;
    }

    elem_type getNull() override {
        return null_;
    }


    bool areRelated(size_type i, size_type j) override {

        auto pis = determineIndices(i, j);
        auto p = partitions_[pis.partition];

        return (p != 0) && p->areRelated(pis.row, pis.col);

    }

    std::vector<size_type> getSuccessors(size_type i) override {

        std::vector<size_type> succs;

        if (hc_ > hr_) {

            size_type offset = 0;
            for (size_type k = 0; k < numPartitions_; k++, offset += partitionSize_) {

                auto p = partitions_[k];
                if (p != 0) {

                    auto tmp = p->getSuccessors(i);
                    for (size_type l = 0; l < tmp.size(); l++) {
                        tmp[l] += offset;
                    }

                    succs.reserve(succs.size() + tmp.size());
                    std::move(tmp.begin(), tmp.end(), std::back_inserter(succs));

                }

            }

        } else {

            auto pis = determineIndices(i, 0);
            auto p = partitions_[pis.partition];
            if (p != 0) {
                succs = p->getSuccessors(pis.row);
            }

        }

        return succs;

    }

    std::vector<size_type> getPredecessors(size_type j) override {

        std::vector<size_type> preds;

        if (hc_ < hr_) {

            size_type offset = 0;
            for (size_type k = 0; k < numPartitions_; k++, offset += partitionSize_) {

                auto p = partitions_[k];
                if (p != 0) {

                    auto tmp = p->getPredecessors(j);
                    for (size_type l = 0; l < tmp.size(); l++) {
                        tmp[l] += offset;
                    }

                    preds.reserve(preds.size() + tmp.size());
                    std::move(tmp.begin(), tmp.end(), std::back_inserter(preds));

                }

            }

        } else {

            auto pis = determineIndices(0, j);
            auto p = partitions_[pis.partition];
            if (p != 0) {
                preds = p->getPredecessors(pis.col);
            }

        }

        return preds;

    }

    positions_type getRange(size_type i1, size_type i2, size_type j1, size_type j2) override {

        positions_type elements;

        auto upperLeft = determineIndices(i1, j1);
        auto lowerRight = determineIndices(i2, j2);

        // range falls completely within one partition
        if (upperLeft.partition == lowerRight.partition) {

            auto p = partitions_[upperLeft.partition];
            if (p != 0) {

                elements = p->getRange(upperLeft.row, lowerRight.row, upperLeft.col, lowerRight.col);

                size_type offset = partitionSize_ * upperLeft.partition;

                if (hc_ > hr_) {

                    for (auto& e : elements) {
                        e.second += offset;
                    }

                } else {

                    for (auto& e : elements) {
                        e.first += offset;
                    }

                }

            }

            return elements;

        }

        // range spans multiple partitions
        std::vector<positions_type> tmp;
        tmp.reserve(lowerRight.partition - upperLeft.partition + 1);
        size_type offset = partitionSize_ * upperLeft.partition;

        if (hc_ > hr_) {

            // first partition (partially spanned)
            auto p = partitions_[upperLeft.partition];
            if (p != 0) {

                tmp.push_back(p->getRange(upperLeft.row, lowerRight.row, upperLeft.col, partitionSize_ - 1));
                for (auto& e : tmp.back()) {
                    e.second += offset;
                }

            }

            // intermediate partition (fully spanned, if any)
            offset += partitionSize_;
            for (size_type k = upperLeft.partition + 1; k < lowerRight.partition; k++, offset += partitionSize_) {

                p = partitions_[k];
                if (p != 0) {

                    tmp.push_back(p->getRange(upperLeft.row, lowerRight.row, 0, partitionSize_ - 1));
                    for (auto& e : tmp.back()) {
                        e.second += offset;
                    }

                }

            }

            // last partition (partially spanned)
            p = partitions_[lowerRight.partition];
            if (p != 0) {

                tmp.push_back(p->getRange(upperLeft.row, lowerRight.row, 0, lowerRight.col));
                for (auto& e : tmp.back()) {
                    e.second += offset;
                }

            }

        } else {

            // first partition (partially spanned)
            auto p = partitions_[upperLeft.partition];
            if (p != 0) {

                tmp.push_back(p->getRange(upperLeft.row, partitionSize_ - 1, upperLeft.col, lowerRight.col));
                for (auto& e : tmp.back()) {
                    e.first += offset;
                }

            }

            // intermediate partition (fully spanned, if any)
            offset += partitionSize_;
            for (size_type k = upperLeft.partition + 1; k < lowerRight.partition; k++, offset += partitionSize_) {

                p = partitions_[k];
                if (p != 0) {

                    tmp.push_back(p->getRange(0, partitionSize_ - 1, upperLeft.col, lowerRight.col));
                    for (auto& e : tmp.back()) {
                        e.first += offset;
                    }

                }

            }

            // last partition (partially spanned)
            p = partitions_[lowerRight.partition];
            if (p != 0) {

                tmp.push_back(p->getRange(0, lowerRight.row, upperLeft.col, lowerRight.col));
                for (auto& e : tmp.back()) {
                    e.first += offset;
                }

            }

        }

        // flatten results
        size_type toReserve = elements.size();
        for (auto& v : tmp) {
            toReserve += v.size();
        }
        elements.reserve(toReserve);

        for (auto& v : tmp) {

            std::move(v.begin(), v.end(), std::back_inserter(elements));
            v.clear();
            v.shrink_to_fit();

        }

        return elements;

    }

    bool containsLink(size_type i1, size_type i2, size_type j1, size_type j2) override {

        auto upperLeft = determineIndices(i1, j1);
        auto lowerRight = determineIndices(i2, j2);
//        auto lowerRight = determineIndices(std::min(i2, numRows_ - 1), std::min(j2, numCols_ - 1));

        // range falls completely within one partition
        if (upperLeft.partition == lowerRight.partition) {

            auto p = partitions_[upperLeft.partition];
            return (p != 0) && p->containsLink(upperLeft.row, lowerRight.row, upperLeft.col, lowerRight.col);

        }

        // range spans multiple partitions
        bool found;

        if (hc_ > hr_) {

            // first partition (partially spanned)
            auto p = partitions_[upperLeft.partition];
            found = (p != 0) && p->containsLink(upperLeft.row, lowerRight.row, upperLeft.col, partitionSize_ - 1);

            // intermediate partition (fully spanned, if any)
            for (size_type k = upperLeft.partition + 1; (k < lowerRight.partition) && !found; k++) {

                p = partitions_[k];
                found = (p != 0) && p->containsLink(upperLeft.row, lowerRight.row, 0, partitionSize_ - 1);

            }

            // last partition (partially spanned)
            p = partitions_[lowerRight.partition];
            found = found || ((p != 0) && p->containsLink(upperLeft.row, lowerRight.row, 0, lowerRight.col));

        } else {

            // first partition (partially spanned)
            auto p = partitions_[upperLeft.partition];
            found =  (p != 0) && p->containsLink(upperLeft.row, partitionSize_ - 1, upperLeft.col, lowerRight.col);

            // intermediate partition (fully spanned, if any)
            for (size_type k = upperLeft.partition + 1; (k < lowerRight.partition) && !found; k++) {

                p = partitions_[k];
                found = (p != 0) && p->containsLink(0, partitionSize_ - 1, upperLeft.col, lowerRight.col);

            }

            // last partition (partially spanned)
            p = partitions_[lowerRight.partition];
            found = found ||  ((p != 0) && p->containsLink(0, lowerRight.row, upperLeft.col, lowerRight.col));

        }

        return found;

    }

    size_type countLinks() override {

        size_type cnt = 0;
        for (size_type k = 0; k < numPartitions_; k++) {

            auto p = partitions_[k];
            cnt += (p != 0) ? p->countLinks() : 0;

        }

        return cnt;

    }


    // general methods for completeness' sake (are redundant / useless for bool)

    bool isNotNull(size_type i, size_type j) override {
        return areRelated(i, j);
    }

    elem_type getElement(size_type i, size_type j) override {
        return areRelated(i, j);
    }

    std::vector<elem_type> getSuccessorElements(size_type i) override {
        return std::vector<elem_type>(getSuccessors(i).size(), true);
    }

    std::vector<size_type> getSuccessorPositions(size_type i) override {
        return getSuccessors(i);
    }

    pairs_type getSuccessorValuedPositions(size_type i) override {

        auto pos = getSuccessors(i);

        pairs_type succs;
        for (auto j : pos) {
            succs.push_back(ValuedPosition<elem_type>(i, j, true));
        }

        return succs;

    }

    std::vector<elem_type> getPredecessorElements(size_type j) override {
        return std::vector<elem_type>(getPredecessors(j).size(), true);
    }

    std::vector<size_type> getPredecessorPositions(size_type j) override {
        return getPredecessors(j);
    }

    pairs_type getPredecessorValuedPositions(size_type j) override {

        auto pos = getPredecessors(j);

        pairs_type preds;
        for (auto i : pos) {
            preds.push_back(ValuedPosition<elem_type>(i, j, true));
        }

        return preds;

    }

    std::vector<elem_type> getElementsInRange(size_type i1, size_type i2, size_type j1, size_type j2) override {
        return std::vector<elem_type>(getRange(i1, i2, j1, j2).size(), true);
    }

    positions_type getPositionsInRange(size_type i1, size_type i2, size_type j1, size_type j2) override {
        return getRange(i1, i2, j1, j2);
    }

    pairs_type getValuedPositionsInRange(size_type i1, size_type i2, size_type j1, size_type j2) override {

        auto pos = getRange(i1, i2, j1, j2);

        pairs_type pairs;
        for (auto& p : pos) {
            pairs.push_back(ValuedPosition<elem_type>(p.first, p.second, true));
        }

        return pairs;

    }

    std::vector<elem_type> getAllElements() override {
        return std::vector<elem_type>(countLinks(), true);
    }

    positions_type getAllPositions() override {
        return getRange(0, numRows_ - 1, 0, numCols_ - 1);
    }

    pairs_type getAllValuedPositions() override {

        auto pos = getAllPositions();

        pairs_type pairs;
        for (auto& p : pos) {
            pairs.push_back(ValuedPosition<elem_type>(p.first, p.second, true));
        }

        return pairs;

    }

    bool containsElement(size_type i1, size_type i2, size_type j1, size_type j2) override {
        return containsLink(i1, i2, j1, j2);
    }

    size_type countElements() override {
        return countLinks();
    }


    UnevenKrKcOrMiniTree* clone() const override {
        return new UnevenKrKcOrMiniTree<elem_type>(*this);
    }


    void print(bool all = false) override {

        std::cout << "### Parameters ###" << std::endl;
        std::cout << "hr  = " << hr_ << std::endl;
        std::cout << "hc  = " << hc_ << std::endl;
        std::cout << "kr  = " << kr_ << std::endl;
        std::cout << "kc  = " << kc_ << std::endl;
        std::cout << "numRows = " << numRows_ << std::endl;
        std::cout << "numCols = " << numCols_ << std::endl;
        std::cout << "partitionSize = " << partitionSize_ << std::endl;
        std::cout << "numPartitions = " << numPartitions_ << std::endl;
        std::cout << "null = " << null_ << std::endl;

        if (all) {

            for (size_type k = 0; k < numPartitions_; k++) {

                std::cout << "===== Partition " << k << " =====" << std::endl;
                auto p = partitions_[k];
                if (p != 0) {
                    p->print(true);
                } else {
                    std::cout << "((ALL NULL))" << std::endl;
                }

            }

        }

    }


    // note: can "invalidate" the data structure (containsLink() probably won't work correctly afterwards)
    void setNull(size_type i, size_type j) override {

        auto pis = determineIndices(i, j);
        auto p = partitions_[pis.partition];

        if (p != 0) {
            p->setNull(pis.row, pis.col);
        }

    }


private:
    size_type hr_;
    size_type hc_;
    size_type kr_;
    size_type kc_;
    size_type numRows_;
    size_type numCols_;

    K2Tree<elem_type>** partitions_;
    size_type partitionSize_;
    size_type numPartitions_;

    elem_type null_;


    /* helper methods for mapping (overall) indices to positions in the partitions */

    size_type determinePartition(size_type i, size_type j) {
        return (hc_ > hr_) ? (j / partitionSize_) : (i / partitionSize_);
    }

    PartitionIndices determineIndices(size_type i, size_type j) {
        return (hc_ > hr_) ? PartitionIndices(j / partitionSize_, i, j % partitionSize_) : PartitionIndices(i / partitionSize_, i % partitionSize_, j);
    }

    /* helper methods for inplace construction from single list of pairs */

    size_type computeKey(const positions_type::value_type& pair, const Subproblem& sp, size_type widthRow, size_type widthCol) {
        return ((pair.first - sp.firstRow) / widthRow) + (pair.second - sp.firstCol) / widthCol;
    }

    void countingSort(positions_type& pairs, std::vector<std::pair<size_type, size_type>>& intervals, const Subproblem& sp, size_type widthRow, size_type widthCol, size_type sup) {

        std::vector<size_type> counts(sup);

        // determine key frequencies
        for (auto i = sp.left; i < sp.right; i++) {
            counts[computeKey(pairs[i], sp, widthRow, widthCol)]++;
        }

        // determine starting index for each key
        size_type total = 0;
        size_type tmp;

        for (auto key = 0; key < sup; key++) {

            tmp = counts[key];
            counts[key] = total;
            total += tmp;

            intervals[key].first = counts[key];
            intervals[key].second = total;

        }

        // reorder pairs of current subproblem
        positions_type tmpPairs(sp.right - sp.left + 1);
        for (auto i = sp.left; i < sp.right; i++) {

            tmpPairs[counts[computeKey(pairs[i], sp, widthRow, widthCol)]] = pairs[i];
            counts[computeKey(pairs[i], sp, widthRow, widthCol)]++;

        }

        for (auto i = sp.left; i < sp.right; i++) {
            pairs[i] = tmpPairs[i - sp.left];
        }

    }

};

#endif //K2TREES_STATICUNEVENRECTANGULARORMINITREE_HPP
