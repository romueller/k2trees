#ifndef K2TREES_STATICHYBRIDTREE_HPP
#define K2TREES_STATICHYBRIDTREE_HPP

#include <queue>

#include "K2Tree.hpp"
#include "Utility.hpp"

template<typename E>
class HybridK2Tree : public virtual K2Tree<E> {

public:
    typedef E elem_type;

    typedef typename K2Tree<elem_type>::matrix_type matrix_type;
    typedef typename K2Tree<elem_type>::list_type list_type;
    typedef typename K2Tree<elem_type>::positions_type positions_type;
    typedef typename K2Tree<elem_type>::pairs_type pairs_type;


    HybridK2Tree() {
        // nothing to do
    }

    HybridK2Tree(const HybridK2Tree& other) {

        h_ = other.h_;
        upperH_ = other.upperH_;
        upperOnes_ = other.upperOnes_;
        upperLength_ = other.upperLength_;
        upperK_ = other.upperK_;
        lowerK_ = other.lowerK_;
        nPrime_ = other.nPrime_;
        null_ = other.null_;

        T_ = other.T_;
        L_ = other.L_;
        R_ = rank_type(&T_);

    }

    HybridK2Tree& operator=(const HybridK2Tree& other) {

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        h_ = other.h_;
        upperH_ = other.upperH_;
        upperOnes_ = other.upperOnes_;
        upperLength_ = other.upperLength_;
        upperK_ = other.upperK_;
        lowerK_ = other.lowerK_;
        nPrime_ = other.nPrime_;
        null_ = other.null_;

        T_ = other.T_;
        L_ = other.L_;
        R_ = rank_type(&T_);

        return *this;

    }

    // assumes that all rows of mat are equally long
    HybridK2Tree(const matrix_type& mat, const size_type upperK, const size_type upperH, const size_type lowerK, const elem_type null = elem_type()) {

        null_ = null;

        upperK_ = upperK;
        lowerK_ = lowerK;

        size_type maxDim = std::max(mat.size(), mat[0].size());

        nPrime_ = size_type(ceil((1.0 * maxDim) / lowerK_));
        upperH_ = std::min(upperH, std::max((size_type)1, logK(nPrime_, upperK_)));

        nPrime_ = size_type(pow(upperK_, upperH_));
        h_ = upperH_;
        do {

            nPrime_ *= lowerK_;
            h_++;

        } while (nPrime_ < maxDim);

        std::vector<std::vector<bool>> levels(h_ - 1);
        buildFromMatrix(mat, levels, mat.size(), mat[0].size(), nPrime_, 1, 0, 0);

        // helper variables (describe top section of conceptual k^2-tree) for navigation on T and L
        upperOnes_ = 0;
        upperLength_ = 0;

        if (upperH_ > 0) {

            for (size_t l = 0; l < upperH_ - 1; l++) {
                for (auto i = 0; i < levels[l].size(); i++) {
                    upperOnes_ += levels[l][i];
                }
            }
            upperLength_ = (upperOnes_ + 1) * upperK_ * upperK_;

        }

        size_type total = 0;
        for (size_type l = 0; l < h_ - 1; l++) {
            total += levels[l].size();
        }
        T_ = bit_vector_type(total);

        bit_vector_type::iterator outIter = T_.begin();
        for (size_type l = 0; l < h_ - 1; l++) {

            outIter = std::move(levels[l].begin(), levels[l].end(), outIter);
            levels[l].clear();
            levels[l].shrink_to_fit();

        }

        R_ = rank_type(&T_);

    }

    HybridK2Tree(const std::vector<list_type>& lists, const size_type upperK, const size_type upperH, const size_type lowerK, const int mode, const elem_type null = elem_type()) {

        null_ = null;

        upperK_ = upperK;
        lowerK_ = lowerK;

        size_type maxDim = 0;
        for (auto& row : lists) {
            for (auto& elem : row) {
                maxDim = std::max(maxDim, elem.first);
            }
        }
        maxDim++; // for number of columns
        maxDim = std::max(maxDim, lists.size());

        nPrime_ = size_type(ceil((1.0 * maxDim) / lowerK_));
        upperH_ = std::min(upperH, std::max((size_type)1, logK(nPrime_, upperK_)));

        nPrime_ = size_type(pow(upperK_, upperH_));
        h_ = upperH_;
        do {

            nPrime_ *= lowerK_;
            h_++;

        } while (nPrime_ < maxDim);

        switch (mode) {

            case 0: {

                std::vector<std::vector<bool>> levels(h_ - 1);
                std::vector<typename list_type::const_iterator> cursors;
                for (auto iter = lists.begin(); iter != lists.end(); iter++) {
                    cursors.push_back(iter->begin());
                }

                buildFromLists(lists, cursors, levels, nPrime_, 1, 0, 0);

                upperOnes_ = 0;
                upperLength_ = 0;

                if (upperH_ > 0) {

                    for (size_t l = 0; l < upperH_ - 1; l++) {
                        for (auto i = 0; i < levels[l].size(); i++) {
                            upperOnes_ += levels[l][i];
                        }
                    }
                    upperLength_ = (upperOnes_ + 1) * upperK_ * upperK_;

                }

                size_type total = 0;
                for (size_type l = 0; l < h_ - 1; l++) {
                    total += levels[l].size();
                }
                T_ = bit_vector_type(total);

                bit_vector_type::iterator outIter = T_.begin();
                for (size_type l = 0; l < h_ - 1; l++) {

                    outIter = std::move(levels[l].begin(), levels[l].end(), outIter);
                    levels[l].clear();
                    levels[l].shrink_to_fit();

                }

                R_ = rank_type(&T_);

                break;

            }

            case 1: {

                buildFromListsViaTree(lists);

                R_ = rank_type(&T_);

                break;

            }

            case 2: {

                buildFromListsDynamicBitmaps(lists);

                break;

            }

            default: {

                buildFromListsDynamicBitmaps(lists);
                break;

            }

        }

    }

    HybridK2Tree(pairs_type& pairs, const size_type upperK, const size_type upperH, const size_type lowerK, const elem_type null = elem_type()) {

        null_ = null;

        upperK_ = upperK;
        lowerK_ = lowerK;

        size_type maxDim = 0;
        for (auto& p : pairs) {
            maxDim = std::max({maxDim, p.row, p.col});
        }
        maxDim++; // for number of rows resp. columns

        nPrime_ = size_type(ceil((1.0 * maxDim) / lowerK_));
        upperH_ = std::min(upperH, std::max((size_type)1, logK(nPrime_, upperK_)));

        nPrime_ = size_type(pow(upperK_, upperH_));
        h_ = upperH_;
        do {

            nPrime_ *= lowerK_;
            h_++;

        } while (nPrime_ < maxDim);

        if (pairs.size() != 0) {
            buildFromListsInplace(pairs);
        }

        R_ = rank_type(&T_);

    }


    size_type getUpperK() {
        return upperK_;
    }

    size_type getLowerK() {
        return lowerK_;
    }

    size_type getUpperH() {
        return upperH_;
    }

    size_type getUpperOnes() {
        return upperOnes_;
    }

    size_type getUpperLength() {
        return upperLength_;
    }

    size_type getH() {
        return upperK_;
    }

    size_type getNumRows() override {
        return nPrime_;
    }

    size_type getNumCols() override {
        return nPrime_;
    }

    elem_type getNull() override {
        return null_;
    }


    bool isNotNull(size_type i, size_type j) override {
        return checkInit(i, j);
    }

    elem_type getElement(size_type i, size_type j) override {
        return getInit(i, j);
    }

    std::vector<elem_type> getSuccessorElements(size_type i) override {

        std::vector<elem_type> succs;
//        successorsElemInit(succs, i);
        allSuccessorElementsIterative(succs, i);

        return succs;

    }

    std::vector<size_type> getSuccessorPositions(size_type i) override {

        std::vector<size_type> succs;
//        successorsPosInit(succs, i);
        allSuccessorPositionsIterative(succs, i);

        return succs;

    }

    pairs_type getSuccessorValuedPositions(size_type i) override {

        pairs_type succs;
//        successorsValPosInit(succs, i);
        allSuccessorValuedPositionsIterative(succs, i);

        return succs;

    }

    std::vector<elem_type> getPredecessorElements(size_type j) override {

        std::vector<elem_type> preds;
        predecessorsElemInit(preds, j);

        return preds;

    }

    std::vector<size_type> getPredecessorPositions(size_type j) override {

        std::vector<size_type> preds;
        predecessorsPosInit(preds, j);

        return preds;

    }

    pairs_type getPredecessorValuedPositions(size_type j) override {

        pairs_type preds;
        predecessorsValPosInit(preds, j);

        return preds;

    }

    std::vector<elem_type> getElementsInRange(size_type i1, size_type i2, size_type j1, size_type j2) override {

        std::vector<elem_type> elements;
        rangeElemInit(elements, i1, i2, j1, j2);
//        rangeElemInit(elements, i1, std::min(i2, nPrime_ - 1), j1, std::min(j2, nPrime_ - 1));

        return elements;

    }

    positions_type getPositionsInRange(size_type i1, size_type i2, size_type j1, size_type j2) override {

        positions_type pairs;
        rangePosInit(pairs, i1, i2, j1, j2);
//        rangePosInit(pairs, i1, std::min(i2, nPrime_ - 1), j1, std::min(j2, nPrime_ - 1));

        return pairs;

    }

    pairs_type getValuedPositionsInRange(size_type i1, size_type i2, size_type j1, size_type j2) override {

        pairs_type pairs;
        rangeValPosInit(pairs, i1, i2, j1, j2);
//        rangeValPosInit(pairs, i1, std::min(i2, nPrime_ - 1), j1, std::min(j2, nPrime_ - 1));

        return pairs;

    }

    std::vector<elem_type> getAllElements() override {
        return getElementsInRange(0, nPrime_ - 1, 0, nPrime_ - 1);
    }

    positions_type getAllPositions() override {
        return getPositionsInRange(0, nPrime_ - 1, 0, nPrime_ - 1);
    }

    pairs_type getAllValuedPositions() override {
        return getValuedPositionsInRange(0, nPrime_ - 1, 0, nPrime_ - 1);
    }

    bool containsElement(size_type i1, size_type i2, size_type j1, size_type j2) override {
        return elemInRangeInit(i1, i2, j1, j2);
//        return elemInRangeInit(i1, std::min(i2, nPrime_ - 1), j1, std::min(j2, nPrime_ - 1));
    }

    size_type countElements() override {

        size_type cnt = 0;
        for (size_type i = 0; i < L_.size(); i++) {
            cnt += (L_[i] != null_);
        }

        return cnt;

    }


    HybridK2Tree* clone() const override {
        return new HybridK2Tree<elem_type>(*this);
    }


    void print(bool all = false) override {

        std::cout << "### Parameters ###" << std::endl;
        std::cout << "h  = " << h_ << std::endl;
        std::cout << "upperH  = " << upperH_ << std::endl;
        std::cout << "upperOnes  = " << upperOnes_ << std::endl;
        std::cout << "upperLength  = " << upperLength_ << std::endl;
        std::cout << "upperK  = " << upperK_ << std::endl;
        std::cout << "lowerK  = " << lowerK_ << std::endl;
        std::cout << "n' = " << nPrime_ << std::endl;
        std::cout << "null = " << null_ << std::endl;

        if (all) {

            std::cout << "### T ###" << std::endl;
            for (size_type i = 0; i < T_.size(); i++) std::cout << T_[i];
            std::cout << std::endl << std::endl;

            std::cout << "### L ###" << std::endl;
            for (size_type i = 0; i < L_.size(); i++) std::cout << L_[i];
            std::cout << std::endl << std::endl;

            std::cout << "### R ###" << std::endl;
            printRanks(R_);
            std::cout << std::endl;

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
        setInit(i, j);
    }

private:
    bit_vector_type T_;
    std::vector<elem_type> L_;
    rank_type R_;

    size_type upperK_;
    size_type lowerK_;
    size_type upperH_;
    size_type upperOnes_;
    size_type upperLength_;
    size_type h_;
    size_type nPrime_;

    elem_type null_;


    /* helper method for construction from relation matrix */

    bool buildFromMatrix(const matrix_type& mat, std::vector<std::vector<bool>>& levels, size_type numRows, size_type numCols, size_type n, size_type l, size_type p, size_type q) {// 3.3.1 / Algorithm 1

        auto k = (l <= upperH_) ? upperK_ : lowerK_;

        if (l == h_) {

            std::vector<elem_type> C;

            for (size_type i = 0; i < k; i++) {
                for (size_type j = 0; j < k; j++) {
                    C.push_back((((p + i) < numRows) && ((q + j) < numCols)) ? mat[p + i][q + j] : null_);
                }
            }

            if (isAll(C, null_)) {
                return false;
            } else {

                L_.insert(L_.end(), C.begin(), C.end());
                return true;

            }

        } else {

            std::vector<bool> C;

            for (size_type i = 0; i < k; i++) {
                for (size_type j = 0; j < k; j++) {
                    C.push_back(buildFromMatrix(mat, levels, numRows, numCols, n / k, l + 1, p + i * (n / k), q + j * (n / k)));
                }
            }

            if (isAllZero(C)) {
                return false;
            } else {

                levels[l - 1].insert(levels[l - 1].end(), C.begin(), C.end());
                return true;

            }

        }

    }

    /* helper method for construction from relation lists */

    bool buildFromLists(const std::vector<list_type>& lists, std::vector<typename list_type::const_iterator>& cursors, std::vector<std::vector<bool>>& levels, size_type n, size_type l, size_type p, size_type q) {// 3.3.2

        auto k = (l <= upperH_) ? upperK_ : lowerK_;

        if (l == h_) {

            std::vector<elem_type> C;

            for (size_type i = 0; i < k; i++) {
                for (size_type j = 0; j < k; j++) {

                    C.push_back(((p + i) < lists.size()) && (cursors[p + i] != lists[p + i].end()) && ((q + j) == cursors[p + i]->first) ? cursors[p + i]->second : null_);
                    if (C.back()) cursors[p + i]++;

                }
            }

            if (isAll(C, null_)) {
                return false;
            } else {

                L_.insert(L_.end(), C.begin(), C.end());
                return true;

            }

        } else {

            std::vector<bool> C;

            for (size_type i = 0; i < k; i++) {
                for (size_type j = 0; j < k; j++) {
                    C.push_back(buildFromLists(lists, cursors, levels, n / k, l + 1, p + i * (n / k), q + j * (n / k)));
                }
            }

            if (isAllZero(C)) {
                return false;
            } else {

                levels[l - 1].insert(levels[l - 1].end(), C.begin(), C.end());
                return true;

            }

        }

    }

    /* helper methods for construction from relation lists via temporary tree */

    void buildFromListsViaTree(const std::vector<list_type>& lists) {// 3.3.3, so far without special bit vectors without initialisation

        Node<elem_type>* root = new Node<elem_type>(null_);

        for (size_type i = 0; i < lists.size(); i++) {
            for (size_type j = 0; j < lists[i].size(); j++) {
                insert(root, nPrime_, i, lists[i][j].first, lists[i][j].second, 0);
            }
        }

        if (!root->isLeaf()) {

            std::vector<bool> T;

            // traverse over tree and generate T and L while doing it
            std::queue<std::pair<Node<elem_type>*, size_type>> queue;
            std::pair<Node<elem_type>*, size_type> node;
            Node<elem_type>* child;
            size_type k;
            queue.push(std::make_pair(root, 0));

            upperOnes_ = 0;
            upperLength_ = 0;

            while (!queue.empty()) {

                node = queue.front();
                queue.pop();

                k = (node.second < upperH_) ? upperK_ : lowerK_;

                for (size_type i = 0; i < k * k; i++) {

                    child = node.first->getChild(i);

                    upperOnes_ += ((upperH_ > 0) && (node.second < (upperH_ - 1))) * (child != 0);

                    if (child != 0 && child->isLeaf()) {
                        L_.push_back(child->getLabel());
                    } else {

                        T.push_back(child != 0);

                        if (T.back()) {
                            queue.push(std::make_pair(child, node.second + 1));
                        }

                    }

                }

                upperLength_ = (upperH_ > 0) ? (upperOnes_ + 1) * upperK_ * upperK_ : 0;

            }

            T_ = bit_vector_type(T.size());
            std::move(T.begin(), T.end(), T_.begin());

        }

        delete root;

    }

    void insert(Node<elem_type>* node, size_type n, size_type p, size_type q, elem_type val, size_type l) {

        auto k = (l < upperH_) ? upperK_ : lowerK_;

        if (n == k) {

            if (node->isLeaf()) {
                node->turnInternal(k * k, true);
            }

            node->addChild(p * k + q, val);

        } else {

            if (node->isLeaf()) {
                node->turnInternal(k * k, false);
            }

            size_type z = (p / (n / k)) * k + q / (n / k);

            insert(node->hasChild(z) ? node->getChild(z) : node->addChild(z, null_), n / k, p % (n / k), q % (n / k), val, l + 1);

        }

    }

    /* helper methods for construction from relation lists via dynamic bitmap representations */

    void buildFromListsDynamicBitmaps(const std::vector<list_type>& lists) {// 3.3.4, currently no succinct dynamic bitmaps

        if (h_ == 1) {

            L_ = std::vector<elem_type>(lowerK_ * lowerK_, null_);

            for (size_type i = 0; i < lists.size(); i++) {
                for (size_type j = 0; j < lists[i].size(); j++) {
                    L_[i * lowerK_ + lists[i][j].first] = lists[i][j].second;
                }
            }

            if (isAll(L_, null_)) {
                L_ = std::vector<elem_type>(0);
            }

        } else {

            std::vector<bool> T;
            NaiveDynamicRank R;

            for (size_type i = 0; i < lists.size(); i++) {
                for (size_type j = 0; j < lists[i].size(); j++) {
                    insertInit(T, R, i, lists[i][j].first, lists[i][j].second);
                }
            }

            T_ = bit_vector_type(T.size());
            std::move(T.begin(), T.end(), T_.begin());

        }

        R_ = rank_type(&T_);

    }

    void insertInit(std::vector<bool>& T, NaiveDynamicRank& R, size_type p, size_type q, elem_type val) {

        auto k = (upperH_ > 0) ? upperK_ : lowerK_;

        if (T.empty()) {

            T = std::vector<bool>(k * k);
            R = NaiveDynamicRank(T);

            upperOnes_ = 0;
            upperLength_ = (upperH_ > 0) ? k * k : 0;

        }

        insert(T, R, nPrime_ / k, p % (nPrime_ / k), q % (nPrime_ / k), val, (p / (nPrime_ / k)) * k + q / (nPrime_ / k), 1);

    }

    void insert(std::vector<bool>& T, NaiveDynamicRank& R, size_type n, size_type p, size_type q, elem_type val, size_type z, size_type l) {

        auto k = (l < upperH_) ? upperK_ : lowerK_;

        if (!T[z]) {

            T[z] = 1;
            R.increaseFrom(z + 1);

            if (l < upperH_) {

                upperOnes_++;
                upperLength_ += k;

            }

            size_type y = (l >= upperH_) * upperLength_ + (R.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k * k;

            if ((l + 1) == h_) {

                L_.insert(L_.begin() + y - T.size(), k * k, null_);
                L_[y + (p / (n / k)) * k + q / (n / k) - T.size()] = val;

            } else {

                T.insert(T.begin() + y, k * k, 0);
                R.insert(y + 1, k * k);

                insert(T, R, n / k, p % (n / k), q % (n / k), val, y, l + 1);

            }

        } else {

            size_type y = (l >= upperH_) * upperLength_ + (R.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k * k + (p / (n / k)) * k + q / (n / k);

            if ((l + 1) == h_) {
                L_[y - T.size()] = val;
            } else {
                insert(T, R, n / k, p % (n / k), q % (n / k), val, y, l + 1);
            }

        }

    }

    /* helper methods for inplace construction from single list of pairs */

    size_type computeKey(const typename pairs_type::value_type& pair, const Subproblem& sp, size_type width, size_type k) {
        return ((pair.row - sp.firstRow) / width) * k + (pair.col - sp.firstCol) / width;
    }

    void countingSort(pairs_type& pairs, std::vector<std::pair<size_type, size_type>>& intervals, const Subproblem& sp, size_type width, size_type sup, size_type k) {

        std::vector<size_type> counts(sup);

        // determine key frequencies
        for (size_type i = sp.left; i < sp.right; i++) {
            counts[computeKey(pairs[i], sp, width, k)]++;
        }

        // determine starting index for each key
        size_type total = 0;
        size_type tmp;

        for (size_type key = 0; key < sup; key++) {

            tmp = counts[key];
            counts[key] = total;
            total += tmp;

            intervals[key].first = counts[key];
            intervals[key].second = total;

        }

        // reorder pairs of current subproblem
        pairs_type tmpPairs(sp.right - sp.left + 1);
        for (size_type i = sp.left; i < sp.right; i++) {

            tmpPairs[counts[computeKey(pairs[i], sp, width, k)]] = pairs[i];
            counts[computeKey(pairs[i], sp, width, k)]++;

        }

        for (size_type i = sp.left; i < sp.right; i++) {
            pairs[i] = tmpPairs[i - sp.left];
        }

    }

    void buildFromListsInplace(pairs_type& pairs) {// 3.3.5

        std::queue<std::pair<Subproblem, size_type>> queue;
        std::pair<Subproblem, size_type> sp;
        size_type k, S;
        std::vector<std::pair<size_type, size_type>> intervals(std::max(upperK_, lowerK_) * std::max(upperK_, lowerK_));
        std::vector<bool> T;
        std::vector<elem_type> appToL;

        upperOnes_ = 0;
        upperLength_ = 0;

        queue.push(std::make_pair(Subproblem(0, nPrime_ - 1, 0, nPrime_ - 1, 0, pairs.size()), 0));

        while (!queue.empty()) {

            sp = queue.front();
            queue.pop();

            k = (sp.second < upperH_) ? upperK_ : lowerK_;
            S = sp.first.lastRow - sp.first.firstRow + 1;

            if (S > k) {

                countingSort(pairs, intervals, sp.first, S / k, k * k, k);

                for (size_type i = 0; i < k * k; i++) {

                    if (intervals[i].first < intervals[i].second) {

                        T.push_back(true);
                        queue.push(std::make_pair(
                                Subproblem(
                                sp.first.firstRow + (i / k) * (S / k),
                                sp.first.firstRow + (i / k + 1) * (S / k) - 1,
                                sp.first.firstCol + (i % k) * (S / k),
                                sp.first.firstCol + (i % k + 1) * (S / k) - 1,
                                sp.first.left + intervals[i].first,
                                sp.first.left + intervals[i].second
                                ),
                                sp.second + 1
                        ));

                        upperOnes_+= ((upperH_ > 0) && (sp.second < upperH_ - 1));

                    } else {
                        T.push_back(false);
                    }

                }

            } else {

                appToL = std::vector<elem_type>(k * k);

                for (size_type i = sp.first.left; i < sp.first.right; i++) {
                    appToL[(pairs[i].row - sp.first.firstRow) * k + (pairs[i].col - sp.first.firstCol)] = pairs[i].val;
                }

                L_.insert(L_.end(), appToL.begin(), appToL.end());

            }

        }

        upperLength_ = (upperH_ > 0) ? (upperOnes_ + 1) * upperK_ * upperK_ : 0;

        T_ = bit_vector_type(T.size());
        std::move(T.begin(), T.end(), T_.begin());

    }


    /* isNotNull() */

    bool checkInit(size_type p, size_type q) {

        size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

        return (L_.empty()) ? false : check(nPrime_ / k, p % (nPrime_ / k), q % (nPrime_ / k), (p / (nPrime_ / k)) * k + q / (nPrime_ / k), 1);

    }

    bool check(size_type n, size_type p, size_type q, size_type z, size_type l) {

        if (z >= T_.size()) {
            return (L_[z - T_.size()] != null_);
        } else {

            auto k = (l < upperH_) ? upperK_ : lowerK_;

            return T_[z] ? check(n / k, p % (n / k), q % (n / k), (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k * k + (p / (n / k)) * k + q / (n / k), l + 1) : false;

        }

    }

    /* getElement() */

    elem_type getInit(size_type p, size_type q) {

        size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

        return (L_.empty()) ? null_ : get(nPrime_ / k, p % (nPrime_ / k), q % (nPrime_ / k), (p / (nPrime_ / k)) * k + q / (nPrime_ / k), 1);

    }

    elem_type get(size_type n, size_type p, size_type q, size_type z, size_type l) {

        if (z >= T_.size()) {
            return L_[z - T_.size()];
        } else {

            auto k = (l < upperH_) ? upperK_ : lowerK_;

            return T_[z] ? get(n / k, p % (n / k), q % (n / k), (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k * k + (p / (n / k)) * k + q / (n / k), l + 1) : null_;

        }

    }


    /* getSuccessorElements() */

    void allSuccessorElementsIterative(std::vector<elem_type>& succs, size_type p) {

        if (L_.empty()) return;

        std::queue<SubrowInfo> queue, nextLevelQueue;
        size_type lenT = T_.size();

        size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

        if (lenT == 0) {

            size_type offset = p * nPrime_;
            for (size_type i = 0; i < nPrime_; i++) {
                if (L_[offset + i] != null_) {
                    succs.push_back(L_[offset + i]);
                }
            }

        } else {

            // successorsPosInit
            size_type n = nPrime_/ k;
            size_type l = 1;
            size_type relP = p;
            for (size_type j = 0, dq = 0, z = k * (relP / n); j < k; j++, dq += n, z++) {
                queue.push(SubrowInfo(dq, z));
            }

            // successorsPos
            relP %= n;
            k = (l < upperH_) ? upperK_ : lowerK_;
            n /= k;
            for (; n > 1; l++, relP %= n, k = (l < upperH_) ? upperK_ : lowerK_, n /= k) {

                size_type a = (l >= upperH_) * upperLength_;
                size_type b = (l >= upperH_) * (upperOnes_ + 1);

                while (!queue.empty()) {

                    auto& cur = queue.front();

                    if (T_[cur.z]) {

                        size_type y = a + (R_.rank(cur.z + 1) - b) * k * k + k * (relP / n);

                        for (size_type j = 0, newDq = cur.dq; j < k; j++, newDq += n, y++) {
                            nextLevelQueue.push(SubrowInfo(newDq, y));
                        }

                    }

                    queue.pop();

                }

                queue.swap(nextLevelQueue);

            }

            size_type a = (l >= upperH_) * upperLength_;
            size_type b = (l >= upperH_) * (upperOnes_ + 1);

            while (!queue.empty()) {

                auto& cur = queue.front();

                if (T_[cur.z]) {

                    size_type y = a + (R_.rank(cur.z + 1) - b) * k * k + k * (relP / n) - lenT;

                    for (size_type j = 0, newDq = cur.dq; j < k; j++, newDq += n, y++) {
                        if (L_[y] != null_) {
                            succs.push_back(L_[y]);
                        }
                    }

                }

                queue.pop();

            }

        }

    }

    void successorsElemInit(std::vector<elem_type>& succs, size_type p) {

        if (!L_.empty()) {

            size_type k = (upperH_ > 0) ? upperK_ : lowerK_;
            size_type y = k * (p / (nPrime_ / k));

            for (size_type j = 0; j < k; j++) {
                successorsElem(succs, nPrime_ / k, p % (nPrime_ / k), (nPrime_ / k) * j, y + j, 1);
            }

        }

    }

    void successorsElem(std::vector<elem_type>& succs, size_type n, size_type p, size_type q, size_type z, size_type l) {

        if (z >= T_.size()) {

            if (L_[z - T_.size()] != null_) {
                succs.push_back(L_[z - T_.size()]);
            }

        } else {

            if (T_[z]) {

                auto k = (l < upperH_) ? upperK_ : lowerK_;
                size_type y = (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k * k + k * (p / (n / k));

                for (size_type j = 0; j < k; j++) {
                    successorsElem(succs, n / k, p % (n / k), q + (n / k) * j, y + j, l + 1);
                }

            }

        }


    }

    /* getSuccessorPositions() */

    void allSuccessorPositionsIterative(std::vector<size_type>& succs, size_type p) {

        if (L_.empty()) return;

        std::queue<SubrowInfo> queue, nextLevelQueue;
        size_type lenT = T_.size();

        size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

        if (lenT == 0) {

            size_type offset = p * nPrime_;
            for (size_type i = 0; i < nPrime_; i++) {
                if (L_[offset + i] != null_) {
                    succs.push_back(i);
                }
            }

        } else {

            // successorsPosInit
            size_type n = nPrime_/ k;
            size_type l = 1;
            size_type relP = p;
            for (size_type j = 0, dq = 0, z = k * (relP / n); j < k; j++, dq += n, z++) {
                queue.push(SubrowInfo(dq, z));
            }

            // successorsPos
            relP %= n;
            k = (l < upperH_) ? upperK_ : lowerK_;
            n /= k;
            for (; n > 1; l++, relP %= n, k = (l < upperH_) ? upperK_ : lowerK_, n /= k) {

                size_type a = (l >= upperH_) * upperLength_;
                size_type b = (l >= upperH_) * (upperOnes_ + 1);

                while (!queue.empty()) {

                    auto& cur = queue.front();

                    if (T_[cur.z]) {

                        size_type y = a + (R_.rank(cur.z + 1) - b) * k * k + k * (relP / n);

                        for (size_type j = 0, newDq = cur.dq; j < k; j++, newDq += n, y++) {
                            nextLevelQueue.push(SubrowInfo(newDq, y));
                        }

                    }

                    queue.pop();

                }

                queue.swap(nextLevelQueue);

            }

            size_type a = (l >= upperH_) * upperLength_;
            size_type b = (l >= upperH_) * (upperOnes_ + 1);

            while (!queue.empty()) {

                auto& cur = queue.front();

                if (T_[cur.z]) {

                    size_type y = a + (R_.rank(cur.z + 1) - b) * k * k + k * (relP / n) - lenT;

                    for (size_type j = 0, newDq = cur.dq; j < k; j++, newDq += n, y++) {
                        if (L_[y] != null_) {
                            succs.push_back(newDq);
                        }
                    }

                }

                queue.pop();

            }

        }

    }

    void successorsPosInit(std::vector<size_type>& succs, size_type p) {

        if (!L_.empty()) {

            size_type k = (upperH_ > 0) ? upperK_ : lowerK_;
            size_type y = k * (p / (nPrime_ / k));

            for (size_type j = 0; j < k; j++) {
                successorsPos(succs, nPrime_ / k, p % (nPrime_ / k), (nPrime_ / k) * j, y + j, 1);
            }

        }

    }

    void successorsPos(std::vector<size_type>& succs, size_type n, size_type p, size_type q, size_type z, size_type l) {

        if (z >= T_.size()) {

            if (L_[z - T_.size()] != null_) {
                succs.push_back(q);
            }

        } else {

            if (T_[z]) {

                auto k = (l < upperH_) ? upperK_ : lowerK_;
                size_type y = (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k * k + k * (p / (n / k));

                for (size_type j = 0; j < k; j++) {
                    successorsPos(succs, n / k, p % (n / k), q + (n / k) * j, y + j, l + 1);
                }

            }

        }


    }

    /* getSuccessorValuedPositions() */

    void allSuccessorValuedPositionsIterative(pairs_type& succs, size_type p) {

        if (L_.empty()) return;

        std::queue<SubrowInfo> queue, nextLevelQueue;
        size_type lenT = T_.size();

        size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

        if (lenT == 0) {

            size_type offset = p * nPrime_;
            for (size_type i = 0; i < nPrime_; i++) {
                if (L_[offset + i] != null_) {
                    succs.push_back(ValuedPosition<elem_type>(p, i, L_[offset + i]));
                }
            }

        } else {

            // successorsPosInit
            size_type n = nPrime_/ k;
            size_type l = 1;
            size_type relP = p;
            for (size_type j = 0, dq = 0, z = k * (relP / n); j < k; j++, dq += n, z++) {
                queue.push(SubrowInfo(dq, z));
            }

            // successorsPos
            relP %= n;
            k = (l < upperH_) ? upperK_ : lowerK_;
            n /= k;
            for (; n > 1; l++, relP %= n, k = (l < upperH_) ? upperK_ : lowerK_, n /= k) {

                size_type a = (l >= upperH_) * upperLength_;
                size_type b = (l >= upperH_) * (upperOnes_ + 1);

                while (!queue.empty()) {

                    auto& cur = queue.front();

                    if (T_[cur.z]) {

                        size_type y = a + (R_.rank(cur.z + 1) - b) * k * k + k * (relP / n);

                        for (size_type j = 0, newDq = cur.dq; j < k; j++, newDq += n, y++) {
                            nextLevelQueue.push(SubrowInfo(newDq, y));
                        }

                    }

                    queue.pop();

                }

                queue.swap(nextLevelQueue);

            }

            size_type a = (l >= upperH_) * upperLength_;
            size_type b = (l >= upperH_) * (upperOnes_ + 1);

            while (!queue.empty()) {

                auto& cur = queue.front();

                if (T_[cur.z]) {

                    size_type y = a + (R_.rank(cur.z + 1) - b) * k * k + k * (relP / n) - lenT;

                    for (size_type j = 0, newDq = cur.dq; j < k; j++, newDq += n, y++) {
                        if (L_[y] != null_) {
                            succs.push_back(ValuedPosition<elem_type>(p, newDq, L_[y]));
                        }
                    }

                }

                queue.pop();

            }

        }

    }

    void successorsValPosInit(pairs_type& succs, size_type p) {

        if (!L_.empty()) {

            size_type k = (upperH_ > 0) ? upperK_ : lowerK_;
            size_type y = k * (p / (nPrime_ / k));

            for (size_type j = 0; j < k; j++) {
                successorsValPos(succs, nPrime_ / k, p % (nPrime_ / k), (nPrime_ / k) * j, y + j, 1);
            }

            for (auto& s : succs) {
                s.row = p;
            }

        }

    }

    void successorsValPos(pairs_type& succs, size_type n, size_type p, size_type q, size_type z, size_type l) {

        if (z >= T_.size()) {

            if (L_[z - T_.size()] != null_) {
                succs.push_back(ValuedPosition<elem_type>(0, q, L_[z - T_.size()]));
            }

        } else {

            if (T_[z]) {

                auto k = (l < upperH_) ? upperK_ : lowerK_;
                size_type y = (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k * k + k * (p / (n / k));

                for (size_type j = 0; j < k; j++) {
                    successorsValPos(succs, n / k, p % (n / k), q + (n / k) * j, y + j, l + 1);
                }

            }

        }


    }


    /* getPredecessorElements() */

    void predecessorsElemInit(std::vector<elem_type>& preds, size_type q) {

        if (!L_.empty()) {

            size_type k = (upperH_ > 0) ? upperK_ : lowerK_;
            size_type y = q / (nPrime_ / k);

            for (size_type i = 0; i < k; i++) {
                predecessorsElem(preds, nPrime_ / k, q % (nPrime_ / k), (nPrime_ / k) * i, y + i * k, 1);
            }

        }

    }

    void predecessorsElem(std::vector<elem_type>& preds, size_type n, size_type q, size_type p, size_type z, size_type l) {

        if (z >= T_.size()) {

            if (L_[z - T_.size()] != null_) {
                preds.push_back(L_[z - T_.size()]);
            }

        } else {

            if (T_[z]) {

                auto k = (l < upperH_) ? upperK_ : lowerK_;
                size_type y = (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k * k + q / (n / k);

                for (size_type i = 0; i < k; i++) {
                    predecessorsElem(preds, n / k, q % (n / k), p + (n / k) * i, y + i * k, l + 1);
                }

            }

        }

    }

    /* getPredecessorPositions() */

    void predecessorsPosInit(std::vector<size_type>& preds, size_type q) {

        if (!L_.empty()) {

            size_type k = (upperH_ > 0) ? upperK_ : lowerK_;
            size_type y = q / (nPrime_ / k);

            for (size_type i = 0; i < k; i++) {
                predecessorsPos(preds, nPrime_ / k, q % (nPrime_ / k), (nPrime_ / k) * i, y + i * k, 1);
            }

        }

    }

    void predecessorsPos(std::vector<size_type>& preds, size_type n, size_type q, size_type p, size_type z, size_type l) {

        if (z >= T_.size()) {

            if (L_[z - T_.size()] != null_) {
                preds.push_back(p);
            }

        } else {

            if (T_[z]) {

                auto k = (l < upperH_) ? upperK_ : lowerK_;
                size_type y = (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k * k + q / (n / k);

                for (size_type i = 0; i < k; i++) {
                    predecessorsPos(preds, n / k, q % (n / k), p + (n / k) * i, y + i * k, l + 1);
                }

            }

        }

    }

    /* getPredecessorValuedPositions() */

    void predecessorsValPosInit(pairs_type& preds, size_type q) {

        if (!L_.empty()) {

            size_type k = (upperH_ > 0) ? upperK_ : lowerK_;
            size_type y = q / (nPrime_ / k);

            for (size_type i = 0; i < k; i++) {
                predecessorsValPos(preds, nPrime_ / k, q % (nPrime_ / k), (nPrime_ / k) * i, y + i * k, 1);
            }

            for (auto& p : preds) {
                p.col = q;
            }

        }

    }

    void predecessorsValPos(pairs_type& preds, size_type n, size_type q, size_type p, size_type z, size_type l) {

        if (z >= T_.size()) {

            if (L_[z - T_.size()] != null_) {
                preds.push_back(ValuedPosition<elem_type>(p, 0, L_[z - T_.size()]));
            }

        } else {

            if (T_[z]) {

                auto k = (l < upperH_) ? upperK_ : lowerK_;
                size_type y = (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k * k + q / (n / k);

                for (size_type i = 0; i < k; i++) {
                    predecessorsValPos(preds, n / k, q % (n / k), p + (n / k) * i, y + i * k, l + 1);
                }

            }

        }

    }


    /* getElementsInRange() */

    void rangeElemInit(std::vector<elem_type>& elements, size_type p1, size_type p2, size_type q1, size_type q2) {

        if (!L_.empty()) {

            size_type p1Prime, p2Prime;

            size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

            for (size_type i = p1 / (nPrime_ / k); i <= p2 / (nPrime_ / k); i++) {

                p1Prime = (i == p1 / (nPrime_ / k)) * (p1 % (nPrime_ / k));
                p2Prime = (i == p2 / (nPrime_ / k)) ? p2 % (nPrime_ / k) : (nPrime_ / k) - 1;

                for (size_type j = q1 / (nPrime_ / k); j <= q2 / (nPrime_ / k); j++) {
                    rangeElem(
                            elements,
                            nPrime_ / k,
                            p1Prime,
                            p2Prime,
                            (j == q1 / (nPrime_ / k)) * (q1 % (nPrime_ / k)),
                            (j == q2 / (nPrime_ / k)) ? q2 % (nPrime_ / k) : (nPrime_ / k) - 1,
                            (nPrime_ / k) * i,
                            (nPrime_ / k) * j,
                            k * i + j,
                            1
                    );
                }

            }

        }

    }

    void rangeElem(std::vector<elem_type>& elements, size_type n, size_type p1, size_type p2, size_type q1, size_type q2, size_type dp, size_type dq, size_type z, size_type l) {

        if (z >= T_.size()) {

            if (L_[z - T_.size()] != null_) {
                elements.push_back(L_[z - T_.size()]);
            }

        } else {

            if (T_[z]) {

                auto k = (l < upperH_) ? upperK_ : lowerK_;
                size_type y = (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k * k;
                size_type p1Prime, p2Prime;

                for (size_type i = p1 / (n / k); i <= p2 / (n / k); i++) {

                    p1Prime = (i == p1 / (n / k)) * (p1 % (n / k));
                    p2Prime = (i == p2 / (n / k)) ? p2 % (n / k) : n / k - 1;

                    for (size_type j = q1 / (n / k); j <= q2 / (n / k); j++) {
                        rangeElem(
                                elements,
                                n / k,
                                p1Prime,
                                p2Prime,
                                (j == q1 / (n / k)) * (q1 % (n / k)),
                                (j == q2 / (n / k)) ? q2 % (n / k) : n / k - 1,
                                dp + (n / k) * i,
                                dq + (n / k) * j,
                                y + k * i + j,
                                l + 1
                        );
                    }

                }

            }

        }

    }

    /* getPositionsInRange() */

    void rangePosInit(positions_type& pairs, size_type p1, size_type p2, size_type q1, size_type q2) {

        if (!L_.empty()) {

            size_type p1Prime, p2Prime;

            size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

            for (size_type i = p1 / (nPrime_ / k); i <= p2 / (nPrime_ / k); i++) {

                p1Prime = (i == p1 / (nPrime_ / k)) * (p1 % (nPrime_ / k));
                p2Prime = (i == p2 / (nPrime_ / k)) ? p2 % (nPrime_ / k) : (nPrime_ / k) - 1;

                for (size_type j = q1 / (nPrime_ / k); j <= q2 / (nPrime_ / k); j++) {
                    rangePos(
                            pairs,
                            nPrime_ / k,
                            p1Prime,
                            p2Prime,
                            (j == q1 / (nPrime_ / k)) * (q1 % (nPrime_ / k)),
                            (j == q2 / (nPrime_ / k)) ? q2 % (nPrime_ / k) : (nPrime_ / k) - 1,
                            (nPrime_ / k) * i,
                            (nPrime_ / k) * j,
                            k * i + j,
                            1
                    );
                }

            }

        }

    }

    void rangePos(positions_type& pairs, size_type n, size_type p1, size_type p2, size_type q1, size_type q2, size_type dp, size_type dq, size_type z, size_type l) {

        if (z >= T_.size()) {

            if (L_[z - T_.size()] != null_) {
                pairs.push_back(std::make_pair(dp, dq));
            }

        } else {

            if (T_[z]) {

                auto k = (l < upperH_) ? upperK_ : lowerK_;
                size_type y = (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k * k;
                size_type p1Prime, p2Prime;

                for (size_type i = p1 / (n / k); i <= p2 / (n / k); i++) {

                    p1Prime = (i == p1 / (n / k)) * (p1 % (n / k));
                    p2Prime = (i == p2 / (n / k)) ? p2 % (n / k) : n / k - 1;

                    for (size_type j = q1 / (n / k); j <= q2 / (n / k); j++) {
                        rangePos(
                                pairs,
                                n / k,
                                p1Prime,
                                p2Prime,
                                (j == q1 / (n / k)) * (q1 % (n / k)),
                                (j == q2 / (n / k)) ? q2 % (n / k) : n / k - 1,
                                dp + (n / k) * i,
                                dq + (n / k) * j,
                                y + k * i + j,
                                l + 1
                        );
                    }

                }

            }

        }

    }

    /* getValuedPositionsRange() */

    void rangeValPosInit(pairs_type& pairs, size_type p1, size_type p2, size_type q1, size_type q2) {

        if (!L_.empty()) {

            size_type p1Prime, p2Prime;

            size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

            for (size_type i = p1 / (nPrime_ / k); i <= p2 / (nPrime_ / k); i++) {

                p1Prime = (i == p1 / (nPrime_ / k)) * (p1 % (nPrime_ / k));
                p2Prime = (i == p2 / (nPrime_ / k)) ? p2 % (nPrime_ / k) : (nPrime_ / k) - 1;

                for (size_type j = q1 / (nPrime_ / k); j <= q2 / (nPrime_ / k); j++) {
                    rangeValPos(
                            pairs,
                            nPrime_ / k,
                            p1Prime,
                            p2Prime,
                            (j == q1 / (nPrime_ / k)) * (q1 % (nPrime_ / k)),
                            (j == q2 / (nPrime_ / k)) ? q2 % (nPrime_ / k) : (nPrime_ / k) - 1,
                            (nPrime_ / k) * i,
                            (nPrime_ / k) * j,
                            k * i + j,
                            1
                    );
                }

            }

        }

    }

    void rangeValPos(pairs_type& pairs, size_type n, size_type p1, size_type p2, size_type q1, size_type q2, size_type dp, size_type dq, size_type z, size_type l) {

        if (z >= T_.size()) {

            if (L_[z - T_.size()] != null_) {
                pairs.push_back(ValuedPosition<elem_type>(dp, dq, L_[z - T_.size()]));
            }

        } else {

            if (T_[z]) {

                auto k = (l < upperH_) ? upperK_ : lowerK_;
                size_type y = (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k * k;
                size_type p1Prime, p2Prime;

                for (size_type i = p1 / (n / k); i <= p2 / (n / k); i++) {

                    p1Prime = (i == p1 / (n / k)) * (p1 % (n / k));
                    p2Prime = (i == p2 / (n / k)) ? p2 % (n / k) : n / k - 1;

                    for (size_type j = q1 / (n / k); j <= q2 / (n / k); j++) {
                        rangeValPos(
                                pairs,
                                n / k,
                                p1Prime,
                                p2Prime,
                                (j == q1 / (n / k)) * (q1 % (n / k)),
                                (j == q2 / (n / k)) ? q2 % (n / k) : n / k - 1,
                                dp + (n / k) * i,
                                dq + (n / k) * j,
                                y + k * i + j,
                                l + 1
                        );
                    }

                }

            }

        }

    }


    /* containsElement() */

    bool elemInRangeInit(size_type p1, size_type p2, size_type q1, size_type q2) {

        if (!L_.empty()) {

            // dividing by k_ (as stated in the paper) in not correct,
            // because it does not use the size of the currently considered submatrix but of its submatrices
            if ((p1 == 0) && (q1 == 0) && (p2 == (nPrime_ /*/ k_*/ - 1)) && (q2 == (nPrime_ /*/ k_*/ - 1))) {
                return 1;
            }

            size_type p1Prime, p2Prime;

            size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

            for (size_type i = p1 / (nPrime_ / k); i <= p2 / (nPrime_ / k); i++) {

                p1Prime = (i == p1 / (nPrime_ / k)) * (p1 % (nPrime_ / k));
                p2Prime = (i == p2 / (nPrime_ / k)) ? (p2 % (nPrime_ / k)) : nPrime_ / k - 1;

                for (size_type j = q1 / (nPrime_ / k); j <= q2 / (nPrime_ / k); j++) {

                    if (elemInRange(nPrime_ / k, p1Prime, p2Prime, (j == q1 / (nPrime_ / k)) * (q1 % (nPrime_ / k)), (j == q2 / (nPrime_ / k)) ? q2 % (nPrime_ / k) : nPrime_ / k - 1, k * i + j, 1)) {
                        return true;
                    }

                }

            }

        }

        return false;

    }

    bool elemInRange(size_type n, size_type p1, size_type p2, size_type q1, size_type q2, size_type z, size_type l) {

        if (z >= T_.size()) {

            return L_[z - T_.size()];

        } else {

            if (T_[z]) {

                // dividing by k_ (as stated in the paper) in not correct,
                // because it does not use the size of the currently considered submatrix but of its submatrices
                if ((p1 == 0) && (q1 == 0) && (p2 == (n /*/ k_*/ - 1)) && (q2 == (n /*/ k_*/ - 1))) {
                    return true;
                }

                auto k = (l < upperH_) ? upperK_ : lowerK_;
                size_type y = (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k * k;
                size_type p1Prime, p2Prime;

                for (size_type i = p1 / (n / k); i <= p2 / (n / k); i++) {

                    p1Prime = (i == p1 / (n / k)) * (p1 % (n / k));
                    p2Prime = (i == p2 / (n / k)) ? (p2 % (n / k)) : n / k - 1;

                    for (size_type j = q1 / (n / k); j <= q2 / (n / k); j++) {

                        if (elemInRange(n / k, p1Prime, p2Prime, (j == q1 / (n / k)) * (q1 % (n / k)), (j == q2 / (n / k)) ? q2 % (n / k) : n / k - 1, y + k * i + j, l + 1)) {
                            return true;
                        }

                    }

                }

            }

            return false;

        }

    }


    /* setNull() */

    void setInit(size_type p, size_type q) {

        size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

        if (!L_.empty()) {
            set(nPrime_ / k, p % (nPrime_ / k), q % (nPrime_ / k), (p / (nPrime_ / k)) * k + q / (nPrime_ / k), 1);
        }

    }

    void set(size_type n, size_type p, size_type q, size_type z, size_type l) {

        if (z >= T_.size()) {
            L_[z - T_.size()] = null_;
        } else {

            auto k = (l < upperH_) ? upperK_ : lowerK_;

            if (T_[z]) {
                set(n / k, p % (n / k), q % (n / k), (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k * k + (p / (n / k)) * k + q / (n / k), l + 1);
            }

        }

    }

};

template<>
class HybridK2Tree<bool> : public virtual K2Tree<bool> {

public:
    typedef bool elem_type;

    typedef K2Tree<elem_type>::matrix_type matrix_type;
    typedef K2Tree<elem_type>::list_type list_type;
    typedef K2Tree<elem_type>::positions_type positions_type;
    typedef K2Tree<elem_type>::pairs_type pairs_type;


    HybridK2Tree() {
        // nothing to do
    }

    HybridK2Tree(const HybridK2Tree& other) {

        h_ = other.h_;
        upperH_ = other.upperH_;
        upperOnes_ = other.upperOnes_;
        upperLength_ = other.upperLength_;
        upperK_ = other.upperK_;
        lowerK_ = other.lowerK_;
        nPrime_ = other.nPrime_;
        null_ = other.null_;

        T_ = other.T_;
        L_ = other.L_;
        R_ = rank_type(&T_);

    }

    HybridK2Tree& operator=(const HybridK2Tree& other) {

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        h_ = other.h_;
        upperH_ = other.upperH_;
        upperOnes_ = other.upperOnes_;
        upperLength_ = other.upperLength_;
        upperK_ = other.upperK_;
        lowerK_ = other.lowerK_;
        nPrime_ = other.nPrime_;
        null_ = other.null_;

        T_ = other.T_;
        L_ = other.L_;
        R_ = rank_type(&T_);

        return *this;

    }

    // assumes that all rows of mat are equally long
    HybridK2Tree(const matrix_type& mat, const size_type upperK, const size_type upperH, const size_type lowerK) {

        null_ = false;

        upperK_ = upperK;
        lowerK_ = lowerK;

        size_type maxDim = std::max(mat.size(), mat[0].size());

        nPrime_ = size_type(ceil((1.0 * maxDim) / lowerK_));
        upperH_ = std::min(upperH, std::max((size_type)1, logK(nPrime_, upperK_)));

        nPrime_ = size_type(pow(upperK_, upperH_));
        h_ = upperH_;
        do {

            nPrime_ *= lowerK_;
            h_++;

        } while (nPrime_ < maxDim);

        std::vector<std::vector<bool>> levels(h_);
        buildFromMatrix(mat, levels, mat.size(), mat[0].size(), nPrime_, 1, 0, 0);

        // helper variables (describe top section of conceptual k^2-tree) for navigation on T and L
        upperOnes_ = 0;
        upperLength_ = 0;

        if (upperH_ > 0) {

            for (size_t l = 0; l < upperH_ - 1; l++) {
                for (auto i = 0; i < levels[l].size(); i++) {
                    upperOnes_ += levels[l][i];
                }
            }
            upperLength_ = (upperOnes_ + 1) * upperK_ * upperK_;

        }

        size_type total = 0;
        for (size_type l = 0; l < h_ - 1; l++) {
            total += levels[l].size();
        }
        T_ = bit_vector_type(total);

        bit_vector_type::iterator outIter = T_.begin();
        for (size_type l = 0; l < h_ - 1; l++) {

            outIter = std::move(levels[l].begin(), levels[l].end(), outIter);
            levels[l].clear();
            levels[l].shrink_to_fit();

        }

        L_ = bit_vector_type(levels[h_ - 1].size());
        std::move(levels[h_ - 1].begin(), levels[h_ - 1].end(), L_.begin());
        levels[h_ - 1].clear();
        levels[h_ - 1].shrink_to_fit();

        R_ = rank_type(&T_);

    }

    HybridK2Tree(const RelationLists& lists, const size_type upperK, const size_type upperH, const size_type lowerK, const int mode) {

        null_ = false;

        upperK_ = upperK;
        lowerK_ = lowerK;

        size_type maxDim = 0;
        for (auto& row : lists) {
            for (auto elem : row) {
                maxDim = std::max(maxDim, elem);
            }
        }
        maxDim++; // for number of columns
        maxDim = std::max(maxDim, lists.size());

        nPrime_ = size_type(ceil((1.0 * maxDim) / lowerK_));
        upperH_ = std::min(upperH, std::max((size_type)1, logK(nPrime_, upperK_)));

        nPrime_ = size_type(pow(upperK_, upperH_));
        h_ = upperH_;
        do {

            nPrime_ *= lowerK_;
            h_++;

        } while (nPrime_ < maxDim);

        switch (mode) {

            case 0: {

                std::vector<std::vector<bool>> levels(h_);
                std::vector<RelationList::const_iterator> cursors;
                for (auto iter = lists.begin(); iter != lists.end(); iter++) {
                    cursors.push_back(iter->begin());
                }

                buildFromLists(lists, cursors, levels, nPrime_, 1, 0, 0);

                upperOnes_ = 0;
                upperLength_ = 0;

                if (upperH_ > 0) {

                    for (size_t l = 0; l < upperH_ - 1; l++) {
                        for (auto i = 0; i < levels[l].size(); i++) {
                            upperOnes_ += levels[l][i];
                        }
                    }
                    upperLength_ = (upperOnes_ + 1) * upperK_ * upperK_;

                }

                size_type total = 0;
                for (size_type l = 0; l < h_ - 1; l++) {
                    total += levels[l].size();
                }
                T_ = bit_vector_type(total);

                bit_vector_type::iterator outIter = T_.begin();
                for (size_type l = 0; l < h_ - 1; l++) {

                    outIter = std::move(levels[l].begin(), levels[l].end(), outIter);
                    levels[l].clear();
                    levels[l].shrink_to_fit();

                }

                L_ = bit_vector_type(levels[h_ - 1].size());
                std::move(levels[h_ - 1].begin(), levels[h_ - 1].end(), L_.begin());
                levels[h_ - 1].clear();
                levels[h_ - 1].shrink_to_fit();

                R_ = rank_type(&T_);

                break;

            }

            case 1: {

                buildFromListsViaTree(lists);

                R_ = rank_type(&T_);

                break;

            }

            case 2: {

                buildFromListsDynamicBitmaps(lists);

                break;

            }

            default: {

                buildFromListsDynamicBitmaps(lists);
                break;

            }

        }

    }

    HybridK2Tree(positions_type& pairs, const size_type upperK, const size_type upperH, const size_type lowerK) {

        null_ = false;

        upperK_ = upperK;
        lowerK_ = lowerK;

        size_type maxDim = 0;
        for (auto& p : pairs) {
            maxDim = std::max({maxDim, p.first, p.second});
        }
        maxDim++; // for number of rows resp. columns

        nPrime_ = size_type(ceil((1.0 * maxDim) / lowerK_));
        upperH_ = std::min(upperH, std::max((size_type)1, logK(nPrime_, upperK_)));

        nPrime_ = size_type(pow(upperK_, upperH_));
        h_ = upperH_;
        do {

            nPrime_ *= lowerK_;
            h_++;

        } while (nPrime_ < maxDim);

        if (pairs.size() != 0) {
            buildFromListsInplace(pairs);
        }

        R_ = rank_type(&T_);

    }


    size_type getUpperK() {
        return upperK_;
    }

    size_type getLowerK() {
        return lowerK_;
    }

    size_type getUpperH() {
        return upperH_;
    }

    size_type getUpperOnes() {
        return upperOnes_;
    }

    size_type getUpperLength() {
        return upperLength_;
    }

    size_type getH() {
        return upperK_;
    }

    size_type getNumRows() override {
        return nPrime_;
    }

    size_type getNumCols() override {
        return nPrime_;
    }

    elem_type getNull() override {
        return null_;
    }


    bool areRelated(size_type i, size_type j) override {
        return checkLinkInit(i, j);
    }

    std::vector<size_type> getSuccessors(size_type i) override {

        std::vector<size_type> succs;
//        successorsInit(succs, i);
        allSuccessorPositionsIterative(succs, i);

        return succs;

    }

    std::vector<size_type> getPredecessors(size_type j) override {

        std::vector<size_type> preds;
        predecessorsInit(preds, j);

        return preds;

    }

    positions_type getRange(size_type i1, size_type i2, size_type j1, size_type j2) override {

        positions_type pairs;
        rangeInit(pairs, i1, i2, j1, j2);
//        rangeInit(pairs, i1, std::min(i2, nPrime_ - 1), j1, std::min(j2, nPrime_ - 1));

        return pairs;

    }

    bool containsLink(size_type i1, size_type i2, size_type j1, size_type j2) override {
        return linkInRangeInit(i1, i2, j1, j2);
//        return linkInRangeInit(i1, std::min(i2, nPrime_ - 1), j1, std::min(j2, nPrime_ - 1));
    }

    size_type countLinks() override {

        size_type res = 0;
        for (size_type i = 0; i < L_.size(); i++) {
            res += L_[i];
        }

        return res;

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
        return getRange(0, nPrime_ - 1, 0, nPrime_ - 1);
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
        return linkInRangeInit(i1, i2, j1, j2);
//        return linkInRangeInit(i1, std::min(i2, nPrime_ - 1), j1, std::min(j2, nPrime_ - 1));
    }

    size_type countElements() override {
        return countLinks();
    }


    HybridK2Tree* clone() const override {
        return new HybridK2Tree<elem_type>(*this);
    }


    void print(bool all = false) override {

        std::cout << "### Parameters ###" << std::endl;
        std::cout << "h  = " << h_ << std::endl;
        std::cout << "upperH  = " << upperH_ << std::endl;
        std::cout << "upperOnes  = " << upperOnes_ << std::endl;
        std::cout << "upperLength  = " << upperLength_ << std::endl;
        std::cout << "upperK  = " << upperK_ << std::endl;
        std::cout << "lowerK  = " << lowerK_ << std::endl;
        std::cout << "n' = " << nPrime_ << std::endl;
        std::cout << "null = " << null_ << std::endl;

        if (all) {

            std::cout << "### T ###" << std::endl;
            for (size_type i = 0; i < T_.size(); i++) std::cout << T_[i];
            std::cout << std::endl << std::endl;

            std::cout << "### L ###" << std::endl;
            for (size_type i = 0; i < L_.size(); i++) std::cout << L_[i];
            std::cout << std::endl << std::endl;

            std::cout << "### R ###" << std::endl;
            printRanks(R_);
            std::cout << std::endl;

        }

    }


    // note: can "invalidate" the data structure (containsLink() probably won't work correctly afterwards)
    void setNull(size_type i, size_type j) override {
        setInit(i, j);
    }


private:
    bit_vector_type T_;
    bit_vector_type L_;
    rank_type R_;

    size_type upperK_;
    size_type lowerK_;
    size_type upperH_;
    size_type upperOnes_;
    size_type upperLength_;
    size_type h_;
    size_type nPrime_;

    elem_type null_;


    /* helper method for construction from relation matrix */

    bool buildFromMatrix(const RelationMatrix& mat, std::vector<std::vector<bool>>& levels, size_type numRows, size_type numCols, size_type n, size_type l, size_type p, size_type q) {// 3.3.1 / Algorithm 1

        auto k = (l <= upperH_) ? upperK_ : lowerK_;

        std::vector<bool> C;

        if (l == h_) {

            for (size_type i = 0; i < k; i++) {
                for (size_type j = 0; j < k; j++) {
                    C.push_back((((p + i) < numRows) && ((q + j) < numCols)) ? mat[p + i][q + j] : false);
                }
            }

        } else {

            for (size_type i = 0; i < k; i++) {
                for (size_type j = 0; j < k; j++) {
                    C.push_back(buildFromMatrix(mat, levels, numRows, numCols, n / k, l + 1, p + i * (n / k), q + j * (n / k)));
                }
            }

        }

        if (isAllZero(C)) {
            return false;
        } else {

            levels[l - 1].insert(levels[l - 1].end(), C.begin(), C.end());
            return true;

        }

    }

    /* helper method for construction from relation lists */

    bool buildFromLists(const RelationLists& lists, std::vector<RelationList::const_iterator>& cursors, std::vector<std::vector<bool>>& levels, size_type n, size_type l, size_type p, size_type q) {// 3.3.2

        auto k = (l <= upperH_) ? upperK_ : lowerK_;

        std::vector<bool> C;

        if (l == h_) {

            for (size_type i = 0; i < k; i++) {
                for (size_type j = 0; j < k; j++) {

                    C.push_back(((p + i) < lists.size()) && (cursors[p + i] != lists[p + i].end()) && ((q + j) == *(cursors[p + i])));
                    if (C.back()) cursors[p + i]++;

                }
            }

        } else {

            for (size_type i = 0; i < k; i++) {
                for (size_type j = 0; j < k; j++) {
                    C.push_back(buildFromLists(lists, cursors, levels, n / k, l + 1, p + i * (n / k), q + j * (n / k)));
                }
            }

        }

        if (isAllZero(C)) {
            return false;
        } else {

            levels[l - 1].insert(levels[l - 1].end(), C.begin(), C.end());
            return true;

        }

    }

    /* helper methods for construction from relation lists via temporary tree */

    void buildFromListsViaTree(const RelationLists& lists) {// 3.3.3, so far without special bit vectors without initialisation

        Node<bool>* root = new Node<bool>(false);

        for (size_type i = 0; i < lists.size(); i++) {
            for (size_type j = 0; j < lists[i].size(); j++) {
                insert(root, nPrime_, i, lists[i][j], 0);
            }
        }

        if (!root->isLeaf()) {

            std::vector<bool> T, L;

            // traverse over tree and generate T and L while doing it
            std::queue<std::pair<Node<bool>*, size_type>> queue;
            std::pair<Node<bool>*, size_type> node;
            Node<bool>* child;
            size_type k;
            queue.push(std::make_pair(root, 0));

            upperOnes_ = 0;
            upperLength_ = 0;

            while (!queue.empty()) {

                node = queue.front();
                queue.pop();

                k = (node.second < upperH_) ? upperK_ : lowerK_;

                for (size_type i = 0; i < k * k; i++) {

                    child = node.first->getChild(i);

                    upperOnes_ += ((upperH_ > 0) && (node.second < (upperH_ - 1))) * (child != 0);

                    if (child != 0 && child->isLeaf()) {
                        L.push_back(child->getLabel());
                    } else {

                        T.push_back(child != 0);

                        if (T.back()) {
                            queue.push(std::make_pair(child, node.second + 1));
                        }

                    }

                }

                upperLength_ = (upperH_ > 0) ? (upperOnes_ + 1) * upperK_ * upperK_ : 0;

            }

            L_ = bit_vector_type(L.size());
            std::move(L.begin(), L.end(), L_.begin());
            L.clear();
            L.shrink_to_fit();

            T_ = bit_vector_type(T.size());
            std::move(T.begin(), T.end(), T_.begin());

        }

        delete root;

    }

    void insert(Node<bool>* node, size_type n, size_type p, size_type q, size_type l) {

        auto k = (l < upperH_) ? upperK_ : lowerK_;

        if (n == k) {

            if (node->isLeaf()) {
                node->turnInternal(k * k, true);
            }

            node->addChild(p * k + q, true);

        } else {

            if (node->isLeaf()) {
                node->turnInternal(k * k, false);
            }

            size_type z = (p / (n / k)) * k + q / (n / k);

            insert(node->hasChild(z) ? node->getChild(z) : node->addChild(z, true), n / k, p % (n / k), q % (n / k), l + 1);

        }

    }

    /* helper methods for construction from relation lists via dynamic bitmap representations */

    void buildFromListsDynamicBitmaps(const RelationLists& lists) {// 3.3.4, currently no succinct dynamic bitmaps

        if (h_ == 1) {

            L_ = bit_vector_type(lowerK_ * lowerK_, 0);

            for (size_type i = 0; i < lists.size(); i++) {
                for (size_type j = 0; j < lists[i].size(); j++) {
                    L_[i * lowerK_ + lists[i][j]] = 1;
                }
            }

            if (isAllZero(L_)) {
                L_ = bit_vector_type(0);
            }

        } else {

            std::vector<bool> T, L;
            NaiveDynamicRank R;

            for (size_type i = 0; i < lists.size(); i++) {
                for (size_type j = 0; j < lists[i].size(); j++) {
                    insertInit(T, L, R, i, lists[i][j]);
                }
            }

            L_ = bit_vector_type(L.size());
            std::move(L.begin(), L.end(), L_.begin());
            L.clear();
            L.shrink_to_fit();

            T_ = bit_vector_type(T.size());
            std::move(T.begin(), T.end(), T_.begin());

        }

        R_ = rank_type(&T_);

    }

    void insertInit(std::vector<bool>& T, std::vector<bool>& L, NaiveDynamicRank& R, size_type p, size_type q) {

        auto k = (upperH_ > 0) ? upperK_ : lowerK_;

        if (T.empty()) {

            T = std::vector<bool>(k * k);
            R = NaiveDynamicRank(T);

            upperOnes_ = 0;
            upperLength_ = (upperH_ > 0) ? k * k : 0;

        }

        insert(T, L, R, nPrime_ / k, p % (nPrime_ / k), q % (nPrime_ / k), (p / (nPrime_ / k)) * k + q / (nPrime_ / k), 1);

    }

    void insert(std::vector<bool>& T, std::vector<bool>& L, NaiveDynamicRank& R, size_type n, size_type p, size_type q, size_type z, size_type l) {

        auto k = (l < upperH_) ? upperK_ : lowerK_;

        if (!T[z]) {

            T[z] = 1;
            R.increaseFrom(z + 1);

            if (l < upperH_) {

                upperOnes_++;
                upperLength_ += k;

            }

            size_type y = (l >= upperH_) * upperLength_ + (R.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k * k;

            if ((l + 1) == h_) {

                L.insert(L.begin() + y - T.size(), k * k, 0);
                L[y + (p / (n / k)) * k + q / (n / k) - T.size()] = 1;

            } else {

                T.insert(T.begin() + y, k * k, 0);
                R.insert(y + 1, k * k);

                insert(T, L, R, n / k, p % (n / k), q % (n / k), y, l + 1);

            }

        } else {

            size_type y = (l >= upperH_) * upperLength_ + (R.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k * k + (p / (n / k)) * k + q / (n / k);

            if ((l + 1) == h_) {
                L[y - T.size()] = 1;
            } else {
                insert(T, L, R, n / k, p % (n / k), q % (n / k), y, l + 1);
            }

        }

    }

    /* helper methods for inplace construction from single list of pairs */

    size_type computeKey(const positions_type::value_type& pair, const Subproblem& sp, size_type width, size_type k) {
        return ((pair.first - sp.firstRow) / width) * k + (pair.second - sp.firstCol) / width;
    }

    void countingSort(positions_type& pairs, std::vector<std::pair<size_type, size_type>>& intervals, const Subproblem& sp, size_type width, size_type sup, size_type k) {

        std::vector<size_type> counts(sup);

        // determine key frequencies
        for (size_type i = sp.left; i < sp.right; i++) {
            counts[computeKey(pairs[i], sp, width, k)]++;
        }

        // determine starting index for each key
        size_type total = 0;
        size_type tmp;

        for (size_type key = 0; key < sup; key++) {

            tmp = counts[key];
            counts[key] = total;
            total += tmp;

            intervals[key].first = counts[key];
            intervals[key].second = total;

        }

        // reorder pairs of current subproblem
        positions_type tmpPairs(sp.right - sp.left + 1);
        for (size_type i = sp.left; i < sp.right; i++) {

            tmpPairs[counts[computeKey(pairs[i], sp, width, k)]] = pairs[i];
            counts[computeKey(pairs[i], sp, width, k)]++;

        }

        for (size_type i = sp.left; i < sp.right; i++) {
            pairs[i] = tmpPairs[i - sp.left];
        }

    }

    void buildFromListsInplace(positions_type& pairs) {// 3.3.5

        std::queue<std::pair<Subproblem, size_type>> queue;
        std::pair<Subproblem, size_type> sp;
        size_type k, S;
        std::vector<std::pair<size_type, size_type>> intervals(std::max(upperK_, lowerK_) * std::max(upperK_, lowerK_));
        std::vector<bool> T, L;
        bit_vector_type appToL;

        upperOnes_ = 0;
        upperLength_ = 0;

        queue.push(std::make_pair(Subproblem(0, nPrime_ - 1, 0, nPrime_ - 1, 0, pairs.size()), 0));

        while (!queue.empty()) {

            sp = queue.front();
            queue.pop();

            k = (sp.second < upperH_) ? upperK_ : lowerK_;
            S = sp.first.lastRow - sp.first.firstRow + 1;

            if (S > k) {

                countingSort(pairs, intervals, sp.first, S / k, k * k, k);

                for (size_type i = 0; i < k * k; i++) {

                    if (intervals[i].first < intervals[i].second) {

                        T.push_back(true);
                        queue.push(std::make_pair(
                                Subproblem(
                                sp.first.firstRow + (i / k) * (S / k),
                                sp.first.firstRow + (i / k + 1) * (S / k) - 1,
                                sp.first.firstCol + (i % k) * (S / k),
                                sp.first.firstCol + (i % k + 1) * (S / k) - 1,
                                sp.first.left + intervals[i].first,
                                sp.first.left + intervals[i].second
                                ),
                                sp.second + 1
                        ));

                        upperOnes_+= ((upperH_ > 0) && (sp.second < upperH_ - 1));

                    } else {
                        T.push_back(false);
                    }

                }

            } else {

                appToL = bit_vector_type(k * k);

                for (size_type i = sp.first.left; i < sp.first.right; i++) {
                    appToL[(pairs[i].first - sp.first.firstRow) * k + (pairs[i].second - sp.first.firstCol)] = true;
                }

                L.insert(L.end(), appToL.begin(), appToL.end());

            }

        }

        upperLength_ = (upperH_ > 0) ? (upperOnes_ + 1) * upperK_ * upperK_ : 0;

        L_ = bit_vector_type(L.size());
        std::move(L.begin(), L.end(), L_.begin());
        L.clear();
        L.shrink_to_fit();

        T_ = bit_vector_type(T.size());
        std::move(T.begin(), T.end(), T_.begin());

    }


    /* areRelated() */

    bool checkLinkInit(size_type p, size_type q) {

        size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

        return (L_.empty()) ? false : checkLink(nPrime_ / k, p % (nPrime_ / k), q % (nPrime_ / k), (p / (nPrime_ / k)) * k + q / (nPrime_ / k), 1);

    }

    bool checkLink(size_type n, size_type p, size_type q, size_type z, size_type l) {

        if (z >= T_.size()) {
            return L_[z - T_.size()];
        } else {

            auto k = (l < upperH_) ? upperK_ : lowerK_;

            return T_[z] ? checkLink(n / k, p % (n / k), q % (n / k), (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k * k + (p / (n / k)) * k + q / (n / k), l + 1) : false;

        }

    }

    /* getSuccessors() */

    void allSuccessorPositionsIterative(std::vector<size_type>& succs, size_type p) {

        if (L_.empty()) return;

        std::queue<SubrowInfo> queue, nextLevelQueue;
        size_type lenT = T_.size();

        size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

        if (lenT == 0) {

            size_type offset = p * nPrime_;
            for (size_type i = 0; i < nPrime_; i++) {
                if (L_[offset + i]) {
                    succs.push_back(i);
                }
            }

        } else {

            // successorsPosInit
            size_type n = nPrime_/ k;
            size_type l = 1;
            size_type relP = p;
            for (size_type j = 0, dq = 0, z = k * (relP / n); j < k; j++, dq += n, z++) {
                queue.push(SubrowInfo(dq, z));
            }

            // successorsPos
            relP %= n;
            k = (l < upperH_) ? upperK_ : lowerK_;
            n /= k;
            for (; n > 1; l++, relP %= n, k = (l < upperH_) ? upperK_ : lowerK_, n /= k) {

                size_type a = (l >= upperH_) * upperLength_;
                size_type b = (l >= upperH_) * (upperOnes_ + 1);

                while (!queue.empty()) {

                    auto& cur = queue.front();

                    if (T_[cur.z]) {

                        size_type y = a + (R_.rank(cur.z + 1) - b) * k * k + k * (relP / n);

                        for (size_type j = 0, newDq = cur.dq; j < k; j++, newDq += n, y++) {
                            nextLevelQueue.push(SubrowInfo(newDq, y));
                        }

                    }

                    queue.pop();

                }

                queue.swap(nextLevelQueue);

            }

            size_type a = (l >= upperH_) * upperLength_;
            size_type b = (l >= upperH_) * (upperOnes_ + 1);

            while (!queue.empty()) {

                auto& cur = queue.front();

                if (T_[cur.z]) {

                    size_type y = a + (R_.rank(cur.z + 1) - b) * k * k + k * (relP / n) - lenT;

                    for (size_type j = 0, newDq = cur.dq; j < k; j++, newDq += n, y++) {
                        if (L_[y]) {
                            succs.push_back(newDq);
                        }
                    }

                }

                queue.pop();

            }

        }

    }

    void successorsInit(std::vector<size_type>& succs, size_type p) {

        if (!L_.empty()) {

            size_type k = (upperH_ > 0) ? upperK_ : lowerK_;
            size_type y = k * (p / (nPrime_ / k));

            for (size_type j = 0; j < k; j++) {
                successors(succs, nPrime_ / k, p % (nPrime_ / k), (nPrime_ / k) * j, y + j, 1);
            }

        }

    }

    void successors(std::vector<size_type>& succs, size_type n, size_type p, size_type q, size_type z, size_type l) {

        if (z >= T_.size()) {

            if (L_[z - T_.size()]) {
                succs.push_back(q);
            }

        } else {

            if (T_[z]) {

                auto k = (l < upperH_) ? upperK_ : lowerK_;
                size_type y = (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k * k + k * (p / (n / k));

                for (size_type j = 0; j < k; j++) {
                    successors(succs, n / k, p % (n / k), q + (n / k) * j, y + j, l + 1);
                }

            }

        }


    }

    /* getPredecessors() */

    void predecessorsInit(std::vector<size_type>& preds, size_type q) {

        if (!L_.empty()) {

            size_type k = (upperH_ > 0) ? upperK_ : lowerK_;
            size_type y = q / (nPrime_ / k);

            for (size_type i = 0; i < k; i++) {
                predecessors(preds, nPrime_ / k, q % (nPrime_ / k), (nPrime_ / k) * i, y + i * k, 1);
            }

        }

    }

    void predecessors(std::vector<size_type>& preds, size_type n, size_type q, size_type p, size_type z, size_type l) {

        if (z >= T_.size()) {

            if (L_[z - T_.size()]) {
                preds.push_back(p);
            }

        } else {

            if (T_[z]) {

                auto k = (l < upperH_) ? upperK_ : lowerK_;
                size_type y = (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k * k + q / (n / k);

                for (size_type i = 0; i < k; i++) {
                    predecessors(preds, n / k, q % (n / k), p + (n / k) * i, y + i * k, l + 1);
                }

            }

        }

    }

    /* getRange() */

    void rangeInit(positions_type& pairs, size_type p1, size_type p2, size_type q1, size_type q2) {

        if (!L_.empty()) {

            size_type p1Prime, p2Prime;

            size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

            for (size_type i = p1 / (nPrime_ / k); i <= p2 / (nPrime_ / k); i++) {

                p1Prime = (i == p1 / (nPrime_ / k)) * (p1 % (nPrime_ / k));
                p2Prime = (i == p2 / (nPrime_ / k)) ? p2 % (nPrime_ / k) : (nPrime_ / k) - 1;

                for (size_type j = q1 / (nPrime_ / k); j <= q2 / (nPrime_ / k); j++) {
                    range(
                            pairs,
                            nPrime_ / k,
                            p1Prime,
                            p2Prime,
                            (j == q1 / (nPrime_ / k)) * (q1 % (nPrime_ / k)),
                            (j == q2 / (nPrime_ / k)) ? q2 % (nPrime_ / k) : (nPrime_ / k) - 1,
                            (nPrime_ / k) * i,
                            (nPrime_ / k) * j,
                            k * i + j,
                            1
                    );
                }

            }

        }

    }

    void range(positions_type& pairs, size_type n, size_type p1, size_type p2, size_type q1, size_type q2, size_type dp, size_type dq, size_type z, size_type l) {

        if (z >= T_.size()) {

            if (L_[z - T_.size()]) {
                pairs.push_back(std::make_pair(dp, dq));
            }

        } else {

            if (T_[z]) {

                auto k = (l < upperH_) ? upperK_ : lowerK_;
                size_type y = (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k * k;
                size_type p1Prime, p2Prime;

                for (size_type i = p1 / (n / k); i <= p2 / (n / k); i++) {

                    p1Prime = (i == p1 / (n / k)) * (p1 % (n / k));
                    p2Prime = (i == p2 / (n / k)) ? p2 % (n / k) : n / k - 1;

                    for (size_type j = q1 / (n / k); j <= q2 / (n / k); j++) {
                        range(
                                pairs,
                                n / k,
                                p1Prime,
                                p2Prime,
                                (j == q1 / (n / k)) * (q1 % (n / k)),
                                (j == q2 / (n / k)) ? q2 % (n / k) : n / k - 1,
                                dp + (n / k) * i,
                                dq + (n / k) * j,
                                y + k * i + j,
                                l + 1
                        );
                    }

                }

            }

        }

    }

    /* linkInRange() */

    bool linkInRangeInit(size_type p1, size_type p2, size_type q1, size_type q2) {

        if (!L_.empty()) {

            // dividing by k_ (as stated in the paper) in not correct,
            // because it does not use the size of the currently considered submatrix but of its submatrices
            if ((p1 == 0) && (q1 == 0) && (p2 == (nPrime_ /*/ k_*/ - 1)) && (q2 == (nPrime_ /*/ k_*/ - 1))) {
                return 1;
            }

            size_type p1Prime, p2Prime;

            size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

            for (size_type i = p1 / (nPrime_ / k); i <= p2 / (nPrime_ / k); i++) {

                p1Prime = (i == p1 / (nPrime_ / k)) * (p1 % (nPrime_ / k));
                p2Prime = (i == p2 / (nPrime_ / k)) ? (p2 % (nPrime_ / k)) : nPrime_ / k - 1;

                for (size_type j = q1 / (nPrime_ / k); j <= q2 / (nPrime_ / k); j++) {

                    if (linkInRange(nPrime_ / k, p1Prime, p2Prime, (j == q1 / (nPrime_ / k)) * (q1 % (nPrime_ / k)), (j == q2 / (nPrime_ / k)) ? q2 % (nPrime_ / k) : nPrime_ / k - 1, k * i + j, 1)) {
                        return true;
                    }

                }

            }

        }

        return false;

    }

    bool linkInRange(size_type n, size_type p1, size_type p2, size_type q1, size_type q2, size_type z, size_type l) {

        if (z >= T_.size()) {

            return L_[z - T_.size()];

        } else {

            if (T_[z]) {

                // dividing by k_ (as stated in the paper) in not correct,
                // because it does not use the size of the currently considered submatrix but of its submatrices
                if ((p1 == 0) && (q1 == 0) && (p2 == (n /*/ k_*/ - 1)) && (q2 == (n /*/ k_*/ - 1))) {
                    return true;
                }

                auto k = (l < upperH_) ? upperK_ : lowerK_;
                size_type y = (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k * k;
                size_type p1Prime, p2Prime;

                for (size_type i = p1 / (n / k); i <= p2 / (n / k); i++) {

                    p1Prime = (i == p1 / (n / k)) * (p1 % (n / k));
                    p2Prime = (i == p2 / (n / k)) ? (p2 % (n / k)) : n / k - 1;

                    for (size_type j = q1 / (n / k); j <= q2 / (n / k); j++) {

                        if (linkInRange(n / k, p1Prime, p2Prime, (j == q1 / (n / k)) * (q1 % (n / k)), (j == q2 / (n / k)) ? q2 % (n / k) : n / k - 1, y + k * i + j, l + 1)) {
                            return true;
                        }

                    }

                }

            }

            return false;

        }

    }


    /* setNull() */

    void setInit(size_type p, size_type q) {

        size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

        if (!L_.empty()) {
            set(nPrime_ / k, p % (nPrime_ / k), q % (nPrime_ / k), (p / (nPrime_ / k)) * k + q / (nPrime_ / k), 1);
        }

    }

    void set(size_type n, size_type p, size_type q, size_type z, size_type l) {

        if (z >= T_.size()) {
            L_[z - T_.size()] = null_;
        } else {

            auto k = (l < upperH_) ? upperK_ : lowerK_;

            if (T_[z]) {
                set(n / k, p % (n / k), q % (n / k), (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k * k + (p / (n / k)) * k + q / (n / k), l + 1);
            }

        }

    }

};

#endif //K2TREES_STATICHYBRIDTREE_HPP
