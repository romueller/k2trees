#ifndef K2TREES_STATICHYBRIDROWTREE_HPP
#define K2TREES_STATICHYBRIDROWTREE_HPP

#include <queue>

#include "RowTree.hpp"
#include "Utility.hpp"

template<typename E>
class HybridRowTree : public virtual RowTree<E> {

public:
    typedef E elem_type;

    typedef typename RowTree<elem_type>::list_type list_type; // position ("column number") + value


    HybridRowTree() {
        // nothing to do
    }

    HybridRowTree(const HybridRowTree& other) {

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

    HybridRowTree& operator=(const HybridRowTree& other) {

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

    HybridRowTree(const std::vector<elem_type>& v, const size_type upperK, const size_type upperH, const size_type lowerK, const elem_type null = elem_type()) {

        // at least one lowerK-level
        // at least zero upperK-levels but never more than upperH upperK-levels

        null_ = null;

        upperK_ = upperK;
        lowerK_ = lowerK;

        nPrime_ = size_type(ceil((1.0 * v.size()) / lowerK_));
        upperH_ = std::min(upperH, std::max((size_type)1, logK(nPrime_, upperK_)));

        nPrime_ = size_type(pow(upperK_, upperH_));
        h_ = upperH_;
        do {

            nPrime_ *= lowerK_;
            h_++;

        } while (nPrime_ < v.size());

        std::vector<std::vector<bool>> levels(h_ - 1);
        buildFromVector(v, levels, nPrime_, 1, 0);

        // helper variables (describe top section of conceptual k^2-tree) for navigation on T and L
        upperOnes_ = 0;
        upperLength_ = 0;

        if (upperH_ > 0) {

            for (size_t l = 0; l < upperH_ - 1; l++) {
                for (auto i = 0; i < levels[l].size(); i++) {
                    upperOnes_ += levels[l][i];
                }
            }
            upperLength_ = (upperOnes_ + 1) * upperK_;

        }


        size_type total = 0;
        for (auto l = 0; l < h_ - 1; l++) {
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

    HybridRowTree(const list_type& list, const size_type upperK, const size_type upperH, const size_type lowerK, const int mode, const elem_type null = elem_type()) {

        // at least one lowerK-level
        // at least zero upperK-levels but never more than upperH upperK-levels

        null_ = null;

        upperK_ = upperK;
        lowerK_ = lowerK;

        size_type maxCol = 0;
        for (auto& e : list) {
            maxCol = std::max(maxCol, e.first);
        }
        maxCol++; // for number of columns

        nPrime_ = size_type(ceil((1.0 * maxCol) / lowerK_));
        upperH_ = std::min(upperH, std::max((size_type)1, logK(nPrime_, upperK_)));

        nPrime_ = size_type(pow(upperK_, upperH_));
        h_ = upperH_;
        do {

            nPrime_ *= lowerK_;
            h_++;

        } while (nPrime_ < maxCol);

        switch (mode) {

            case 0: {

                std::vector<std::vector<bool>> levels(h_ - 1);
                auto cursor = list.begin();

                buildFromList(list, cursor, levels, nPrime_, 1, 0);

                upperOnes_ = 0;
                upperLength_ = 0;

                if (upperH_ > 0) {

                    for (size_t l = 0; l < upperH_ - 1; l++) {
                        for (auto i = 0; i < levels[l].size(); i++) {
                            upperOnes_ += levels[l][i];
                        }
                    }
                    upperLength_ = (upperOnes_ + 1) * upperK_;

                }

                size_type total = 0;
                for (auto l = 0; l < h_ - 1; l++) {
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

                buildFromListsViaTree(list);

                R_ = rank_type(&T_);

                break;

            }

            case 2: {

                buildFromListsDynamicBitmaps(list);

                break;

            }

            default: {

                buildFromListsDynamicBitmaps(list);
                break;

            }

        }

    }

    HybridRowTree(list_type& pairs, const size_type upperK, const size_type upperH, const size_type lowerK, const elem_type null = elem_type()) {

        // at least one lowerK-level
        // at least zero upperK-levels but never more than upperH upperK-levels

        null_ = null;

        upperK_ = upperK;
        lowerK_ = lowerK;

        size_type maxCol = 0;
        for (auto& e : pairs) {
            maxCol = std::max(maxCol, e.first);
        }
        maxCol++; // for number of columns

        nPrime_ = size_type(ceil((1.0 * maxCol) / lowerK_));
        upperH_ = std::min(upperH, std::max((size_type)1, logK(nPrime_, upperK_)));

        nPrime_ = size_type(pow(upperK_, upperH_));
        h_ = upperH_;
        do {

            nPrime_ *= lowerK_;
            h_++;

        } while (nPrime_ < maxCol);

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

    size_type getLength() override {
        return nPrime_;
    }

    elem_type getNull() override {
        return null_;
    }


    bool isNotNull(size_type i) override {
        return checkInit(i);
    }

    elem_type getElement(size_type i) override {
        return getInit(i);
    }

    std::vector<elem_type> getElementsInRange(size_type l, size_type r) override {

        std::vector<elem_type> elems;
        rangeElemInit(elems, l, r);
//        rangeElemInit(elems, l, std::min(r, nPrime_ - 1));

        return elems;

    }

    std::vector<size_type> getPositionsInRange(size_type l, size_type r) override {

        std::vector<size_type> positions;
        rangePosInit(positions, l, r);
//        rangePosInit(positions, l, std::min(r, nPrime_ - 1));

        return positions;

    }

    list_type getValuedPositionsInRange(size_type l, size_type r) override {

        list_type positions;
        rangeValPosInit(positions, l, r);
//        rangeValPosInit(positions, l, std::min(r, nPrime_ - 1));

        return positions;

    }

    std::vector<elem_type> getAllElements() override {
//        return getElementsInRange(0, nPrime_ - 1);
        std::vector<elem_type> elements;
        fullRangeElemIterative(elements);
        return elements;

    }

    std::vector<size_type> getAllPositions() override {
//        return getPositionsInRange(0, nPrime_ - 1);
        std::vector<size_type> positions;
        fullRangePosIterative(positions);
        return positions;

    }

    list_type getAllValuedPositions() override {
//        return getValuedPositionsInRange(0, nPrime_ - 1);
        list_type positions;
        fullRangeValPosIterative(positions);
        return positions;

    }

    bool containsElement(size_type l, size_type r) override {
        return elemInRangeInit(l, r);
//        return elemInRangeInit(l, std::min(r, nPrime_ - 1));
    }

    size_type countElements() override {

        size_type res = 0;
        for (auto i = 0; i < L_.size(); i++) {
            res += (L_[i] != null_);
        }

        return res;

    }


    HybridRowTree* clone() const override {
        return new HybridRowTree<elem_type>(*this);
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
            for (auto i = 0; i < T_.size(); i++) std::cout << T_[i];
            std::cout << std::endl << std::endl;

            std::cout << "### L ###" << std::endl;
            for (auto i = 0; i < L_.size(); i++) std::cout << L_[i];
            std::cout << std::endl << std::endl;

            std::cout << "### R ###" << std::endl;
            printRanks(R_);
            std::cout << std::endl;

        }

    }


    // note: can "invalidate" the data structure (containsLink() probably won't work correctly afterwards)
    void setNull(size_type i) override {
        setInit(i);
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

    /* helper method for construction from vector */

    bool buildFromVector(const std::vector<elem_type>& v, std::vector<std::vector<bool>>& levels, size_type numCols, size_type l, size_type q) {

        auto k = (l <= upperH_) ? upperK_ : lowerK_;

        if (l == h_) {

            std::vector<elem_type> C;

            for (auto j = 0; j < k; j++) {
                C.push_back(((q + j) < v.size()) ? v[q + j] : null_);
            }

            if (isAll(C, null_)) {
                return false;
            } else {

                L_.insert(L_.end(), C.begin(), C.end());
                return true;

            }

        } else {

            std::vector<bool> C;

            for (auto j = 0; j < k; j++) {
                C.push_back(buildFromVector(v, levels, numCols / k, l + 1, q + j * (numCols / k)));
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

    bool buildFromList(const list_type& list, typename list_type::const_iterator& cursor, std::vector<std::vector<bool>>& levels, size_type n, size_type l, size_type q) {// 3.3.2

        auto k = (l <= upperH_) ? upperK_ : lowerK_;

        if (l == h_) {

            std::vector<elem_type> C;

            for (auto j = 0; j < k; j++) {

                C.push_back((cursor != list.end()) && ((q + j) == cursor->first) ? cursor->second : null_);
                if (C.back()) cursor++;

            }

            if (isAll(C, null_)) {
                return false;
            } else {

                L_.insert(L_.end(), C.begin(), C.end());
                return true;

            }

        } else {

            std::vector<bool> C;

            for (auto j = 0; j < k; j++) {
                C.push_back(buildFromList(list, cursor, levels, n / k, l + 1, q + j * (n / k)));
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

    void buildFromListsViaTree(const list_type& list) {// 3.3.3, so far without special bit vectors without initialisation

        Node<elem_type>* root = new Node<elem_type>(null_);

        for (auto j = 0; j < list.size(); j++) {
            insert(root, nPrime_, list[j].first, list[j].second, 0);
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

                for (size_type i = 0; i < k; i++) {

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

                upperLength_ = (upperH_ > 0) ? (upperOnes_ + 1) * upperK_ : 0;

            }

            T_ = bit_vector_type(T.size());
            std::move(T.begin(), T.end(), T_.begin());

        }

        delete root;

    }

    void insert(Node<elem_type>* node, size_type n, size_type q, elem_type val, size_type l) {

        auto k = (l < upperH_) ? upperK_ : lowerK_;

        if (n == k) {

            if (node->isLeaf()) {
                node->turnInternal(k, true);
            }

            node->addChild(q, val);

        } else {

            if (node->isLeaf()) {
                node->turnInternal(k, false);
            }

            size_type z = q / (n / k);

            insert(node->hasChild(z) ? node->getChild(z) : node->addChild(z, null_), n / k, q % (n / k), val, l + 1);

        }

    }

    /* helper methods for construction from relation lists via dynamic bitmap representations */

    void buildFromListsDynamicBitmaps(const list_type& list) {// 3.3.4, currently no succinct dynamic bitmaps

        if (h_ == 1) {

            L_ = std::vector<elem_type>(lowerK_, null_);

            for (auto j = 0; j < list.size(); j++) {
                L_[list[j].first] = list[j].second;
            }

            if (isAll(L_, null_)) {
                L_ = std::vector<elem_type>(0);
            }

        } else {

            std::vector<bool> T;
            NaiveDynamicRank R;

            for (auto j = 0; j < list.size(); j++) {
                insertInit(T, R, list[j].first, list[j].second);
            }

            T_ = bit_vector_type(T.size());
            std::move(T.begin(), T.end(), T_.begin());

        }

        R_ = rank_type(&T_);

    }

    void insertInit(std::vector<bool>& T, NaiveDynamicRank& R, size_type q, elem_type val) {

        auto k = (upperH_ > 0) ? upperK_ : lowerK_;

        if (T.empty()) {

            T = std::vector<bool>(k);
            R = NaiveDynamicRank(T);

            upperOnes_ = 0;
            upperLength_ = (upperH_ > 0) ? k : 0;

        }

        insert(T, R, nPrime_ / k, q % (nPrime_ / k), val, q / (nPrime_ / k), 1);

    }

    void insert(std::vector<bool>& T, NaiveDynamicRank& R, size_type n, size_type q, elem_type val, size_type z, size_type l) {

        auto k = (l < upperH_) ? upperK_ : lowerK_;

        if (!T[z]) {

            T[z] = 1;
            R.increaseFrom(z + 1);

            if (l < upperH_) {

                upperOnes_++;
                upperLength_ += k;

            }

            size_type y = (l >= upperH_) * upperLength_ + (R.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k;

            if ((l + 1) == h_) {

                L_.insert(L_.begin() + y - T.size(), k, null_);
                L_[y + q / (n / k) - T.size()] = val;

            } else {

                T.insert(T.begin() + y, k, null_);
                R.insert(y + 1, k);

                insert(T, R, n / k, q % (n / k), val, y + q / (n / k), l + 1);

            }

        } else {

            size_type y = (l >= upperH_) * upperLength_ + (R.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k + q / (n / k);

            if ((l + 1) == h_) {
                L_[y - T.size()] = val;
            } else {
                insert(T, R, n / k, q % (n / k), val, y, l + 1);
            }

        }

    }

    /* helper methods for inplace construction from single list of pairs */

    size_type computeKey(const typename list_type::value_type& pair, const Subproblem& sp, size_type width) {
        return (pair.first - sp.firstCol) / width;
    }

    void countingSort(list_type& pairs, std::vector<std::pair<size_type, size_type>>& intervals, const Subproblem& sp, size_type width, size_type sup) {

        std::vector<size_type> counts(sup);

        // determine key frequencies
        for (auto i = sp.left; i < sp.right; i++) {
            counts[computeKey(pairs[i], sp, width)]++;
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
        list_type tmpPairs(sp.right - sp.left + 1);
        for (auto i = sp.left; i < sp.right; i++) {

            tmpPairs[counts[computeKey(pairs[i], sp, width)]] = pairs[i];
            counts[computeKey(pairs[i], sp, width)]++;

        }

        for (auto i = sp.left; i < sp.right; i++) {
            pairs[i] = tmpPairs[i - sp.left];
        }

    }

    void buildFromListsInplace(list_type& pairs) {// 3.3.5

        std::queue<std::pair<Subproblem, size_type>> queue;
        std::pair<Subproblem, size_type> sp;
        size_type k, S;
        std::vector<std::pair<size_type, size_type>> intervals(std::max(upperK_, lowerK_) * std::max(upperK_, lowerK_));
        std::vector<bool> T;
        std::vector<elem_type> appToL;

        upperOnes_ = 0;
        upperLength_ = 0;

        queue.push(std::make_pair(Subproblem(0, 0, 0, nPrime_ - 1, 0, pairs.size()), 0));

        while (!queue.empty()) {

            sp = queue.front();
            queue.pop();

            k = (sp.second < upperH_) ? upperK_ : lowerK_;
            S = sp.first.lastCol - sp.first.firstCol + 1;

            if (S > k) {

                countingSort(pairs, intervals, sp.first, S / k, k);

                for (auto i = 0; i < k; i++) {

                    if (intervals[i].first < intervals[i].second) {

                        T.push_back(true);
                        queue.push(std::make_pair(
                                Subproblem(
                                        0,
                                        0,
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

                appToL = std::vector<elem_type>(k);

                for (auto i = sp.first.left; i < sp.first.right; i++) {
                    appToL[pairs[i].first - sp.first.firstCol] = pairs[i].second;
                }

                L_.insert(L_.end(), appToL.begin(), appToL.end());

            }

        }

        upperLength_ = (upperH_ > 0) ? (upperOnes_ + 1) * upperK_ : 0;

        T_ = bit_vector_type(T.size());
        std::move(T.begin(), T.end(), T_.begin());

    }


    /* isNotNull() */

    bool checkInit(size_type q) {

        size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

        return (L_.empty()) ? false : check(nPrime_ / k, q % (nPrime_ / k), q / (nPrime_ / k), 1);

    }

    bool check(size_type n, size_type q, size_type z, size_type l) {

        if (z >= T_.size()) {
            return (L_[z - T_.size()] != null_);
        } else {

            auto k = (l < upperH_) ? upperK_ : lowerK_;

            return T_[z] ? check(n / k, q % (n / k), (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k + q / (n / k), l + 1) : false;

        }

    }

    /* getElement() */

    elem_type getInit(size_type q) {

        size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

        return (L_.empty()) ? null_ : get(nPrime_ / k, q % (nPrime_ / k), q / (nPrime_ / k), 1);

    }

    elem_type get(size_type n, size_type q, size_type z, size_type l) {

        if (z >= T_.size()) {
            return L_[z - T_.size()];
        } else {

            auto k = (l < upperH_) ? upperK_ : lowerK_;

            return T_[z] ? get(n / k, q % (n / k), (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k + q / (n / k), l + 1) : null_;

        }

    }

    /* getElementsInRange() */

    void fullRangeElemIterative(std::vector<elem_type>& elems) {

        if (L_.empty()) return;

        std::queue<SubrowInfo> queue, nextLevelQueue;
        size_type lenT = T_.size();
        size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

        if (lenT == 0) {

            for (size_type i = 0; i < nPrime_; i++) {
                if (L_[i] != null_) {
                    elems.push_back(L_[i]);
                }
            }

        } else {

            // rangeInit
            size_type n = nPrime_/ k;
            size_type l = 1;

            for (size_type z = 0, dq = 0; z < k; z++, dq += n) {
                queue.push(SubrowInfo(dq, z));
            }

            // range
            k = (l < upperH_) ? upperK_ : lowerK_;
            n /= k;
            for (; n > 1; l++, k = (l < upperH_) ? upperK_ : lowerK_, n /= k) {

                size_type a = (l >= upperH_) * upperLength_;
                size_type b = (l >= upperH_) * (upperOnes_ + 1);

                while (!queue.empty()) {

                    auto& cur = queue.front();

                    if (T_[cur.z]) {

                        auto y = a + (R_.rank(cur.z + 1) - b) * k;

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

                    auto y = a + (R_.rank(cur.z + 1) - b) * k - lenT;

                    for (size_type j = 0, newDq = cur.dq; j < k; j++, newDq += n, y++) {
                        if (L_[y] != null_) {
                            elems.push_back(L_[y]);
                        }
                    }

                }

                queue.pop();

            }

        }

    }

    void rangeElemInit(std::vector<elem_type>& elems, size_type l, size_type r) {

        if (!L_.empty()) {

            size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

            for (auto j = l / (nPrime_ / k); j <= r / (nPrime_ / k); j++) {
                rangeElem(
                        elems,
                        nPrime_ / k,
                        (j == l / (nPrime_ / k)) * (l % (nPrime_ / k)),
                        (j == r / (nPrime_ / k)) ? r % (nPrime_ / k) : (nPrime_ / k) - 1,
                        (nPrime_ / k) * j,
                        j,
                        1
                );
            }

        }

    }

    void rangeElem(std::vector<elem_type>& elems, size_type n, size_type l, size_type r, size_type dq, size_type z, size_type level) {

        if (z >= T_.size()) {

            if (L_[z - T_.size()] != null_) {
                elems.push_back(L_[z - T_.size()]);
            }

        } else {

            if (T_[z]) {

                auto k = (level < upperH_) ? upperK_ : lowerK_;
                auto y = (level >= upperH_) * upperLength_ + (R_.rank(z + 1) - (level >= upperH_) * (upperOnes_ + 1)) * k;

                for (auto j = l / (n / k); j <= r / (n / k); j++) {
                    rangeElem(
                            elems,
                            n / k,
                            (j == l / (n / k)) * (l % (n / k)),
                            (j == r / (n / k)) ? r % (n / k) : n / k - 1,
                            dq + (n / k) * j,
                            y + j,
                            level + 1
                    );
                }

            }

        }

    }

    /* getPositionsInRange() */

    void fullRangePosIterative(std::vector<size_type>& elems) {

        if (L_.empty()) return;

        std::queue<SubrowInfo> queue, nextLevelQueue;
        size_type lenT = T_.size();
        size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

        if (lenT == 0) {

            for (size_type i = 0; i < nPrime_; i++) {
                if (L_[i] != null_) {
                    elems.push_back(i);
                }
            }

        } else {

            // rangeInit
            size_type n = nPrime_/ k;
            size_type l = 1;

            for (size_type z = 0, dq = 0; z < k; z++, dq += n) {
                queue.push(SubrowInfo(dq, z));
            }

            // range
            k = (l < upperH_) ? upperK_ : lowerK_;
            n /= k;
            for (; n > 1; l++, k = (l < upperH_) ? upperK_ : lowerK_, n /= k) {

                size_type a = (l >= upperH_) * upperLength_;
                size_type b = (l >= upperH_) * (upperOnes_ + 1);

                while (!queue.empty()) {

                    auto& cur = queue.front();

                    if (T_[cur.z]) {

                        auto y = a + (R_.rank(cur.z + 1) - b) * k;

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

                    auto y = a + (R_.rank(cur.z + 1) - b) * k - lenT;

                    for (size_type j = 0, newDq = cur.dq; j < k; j++, newDq += n, y++) {
                        if (L_[y] != null_) {
                            elems.push_back(newDq);
                        }
                    }

                }

                queue.pop();

            }

        }

    }

    void rangePosInit(std::vector<size_type>& elems, size_type l, size_type r) {

        if (!L_.empty()) {

            size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

            for (auto j = l / (nPrime_ / k); j <= r / (nPrime_ / k); j++) {
                rangePos(
                        elems,
                        nPrime_ / k,
                        (j == l / (nPrime_ / k)) * (l % (nPrime_ / k)),
                        (j == r / (nPrime_ / k)) ? r % (nPrime_ / k) : (nPrime_ / k) - 1,
                        (nPrime_ / k) * j,
                        j,
                        1
                );
            }

        }

    }

    void rangePos(std::vector<size_type>& elems, size_type n, size_type l, size_type r, size_type dq, size_type z, size_type level) {

        if (z >= T_.size()) {

            if (L_[z - T_.size()] != null_) {
                elems.push_back(dq);
            }

        } else {

            if (T_[z]) {

                auto k = (level < upperH_) ? upperK_ : lowerK_;
                auto y = (level >= upperH_) * upperLength_ + (R_.rank(z + 1) - (level >= upperH_) * (upperOnes_ + 1)) * k;

                for (auto j = l / (n / k); j <= r / (n / k); j++) {
                    rangePos(
                            elems,
                            n / k,
                            (j == l / (n / k)) * (l % (n / k)),
                            (j == r / (n / k)) ? r % (n / k) : n / k - 1,
                            dq + (n / k) * j,
                            y + j,
                            level + 1
                    );
                }

            }

        }

    }

    /* getValuedPositionsInRange() */

    void fullRangeValPosIterative(list_type& elems) {

        if (L_.empty()) return;

        std::queue<SubrowInfo> queue, nextLevelQueue;
        size_type lenT = T_.size();
        size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

        if (lenT == 0) {

            for (size_type i = 0; i < nPrime_; i++) {
                if (L_[i] != null_) {
                    elems.push_back(std::make_pair(i, L_[i]));
                }
            }

        } else {

            // rangeInit
            size_type n = nPrime_/ k;
            size_type l = 1;

            for (size_type z = 0, dq = 0; z < k; z++, dq += n) {
                queue.push(SubrowInfo(dq, z));
            }

            // range
            k = (l < upperH_) ? upperK_ : lowerK_;
            n /= k;
            for (; n > 1; l++, k = (l < upperH_) ? upperK_ : lowerK_, n /= k) {

                size_type a = (l >= upperH_) * upperLength_;
                size_type b = (l >= upperH_) * (upperOnes_ + 1);

                while (!queue.empty()) {

                    auto& cur = queue.front();

                    if (T_[cur.z]) {

                        auto y = a + (R_.rank(cur.z + 1) - b) * k;

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

                    auto y = a + (R_.rank(cur.z + 1) - b) * k - lenT;

                    for (size_type j = 0, newDq = cur.dq; j < k; j++, newDq += n, y++) {
                        if (L_[y] != null_) {
                            elems.push_back(std::make_pair(newDq, L_[y]));
                        }
                    }

                }

                queue.pop();

            }

        }

    }

    void rangeValPosInit(list_type& elems, size_type l, size_type r) {

        if (!L_.empty()) {

            size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

            for (auto j = l / (nPrime_ / k); j <= r / (nPrime_ / k); j++) {
                rangeValPos(
                        elems,
                        nPrime_ / k,
                        (j == l / (nPrime_ / k)) * (l % (nPrime_ / k)),
                        (j == r / (nPrime_ / k)) ? r % (nPrime_ / k) : (nPrime_ / k) - 1,
                        (nPrime_ / k) * j,
                        j,
                        1
                );
            }

        }

    }

    void rangeValPos(list_type& elems, size_type n, size_type l, size_type r, size_type dq, size_type z, size_type level) {

        if (z >= T_.size()) {

            if (L_[z - T_.size()] != null_) {
                elems.push_back(std::make_pair(dq, L_[z - T_.size()]));
            }

        } else {

            if (T_[z]) {

                auto k = (level < upperH_) ? upperK_ : lowerK_;
                auto y = (level >= upperH_) * upperLength_ + (R_.rank(z + 1) - (level >= upperH_) * (upperOnes_ + 1)) * k;

                for (auto j = l / (n / k); j <= r / (n / k); j++) {
                    rangeValPos(
                            elems,
                            n / k,
                            (j == l / (n / k)) * (l % (n / k)),
                            (j == r / (n / k)) ? r % (n / k) : n / k - 1,
                            dq + (n / k) * j,
                            y + j,
                            level + 1
                    );
                }

            }

        }

    }

    /* containsElement() */

    bool elemInRangeInit(size_type l, size_type r) {

        if (!L_.empty()) {

            if ((l == 0) && (r == (nPrime_ - 1))) {
                return 1;
            }

            size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

            for (auto j = l / (nPrime_ / k); j <= r / (nPrime_ / k); j++) {

                if (elemInRange(nPrime_ / k, (j == l / (nPrime_ / k)) * (l % (nPrime_ / k)), (j == r / (nPrime_ / k)) ? r % (nPrime_ / k) : nPrime_ / k - 1, j, 1)) {
                    return true;
                }

            }

        }

        return false;

    }

    bool elemInRange(size_type n, size_type l, size_type r, size_type z, size_type level) {

        if (z >= T_.size()) {

            return L_[z - T_.size()] != null_;

        } else {

            if (T_[z]) {

                if ((l == 0) && (r == (n - 1))) {
                    return true;
                }

                auto k = (level < upperH_) ? upperK_ : lowerK_;
                auto y = (level >= upperH_) * upperLength_ + (R_.rank(z + 1) - (level >= upperH_) * (upperOnes_ + 1)) * k;

                for (auto j = l / (n / k); j <= r / (n / k); j++) {

                    if (elemInRange(n / k, (j == l / (n / k)) * (l % (n / k)), (j == r / (n / k)) ? r % (n / k) : n / k - 1, y + j, level + 1)) {
                        return true;
                    }

                }

            }

            return false;

        }

    }


    /* setNull() */

    void setInit(size_type q) {

        size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

        if (!L_.empty()) {
            set(nPrime_ / k, q % (nPrime_ / k), q / (nPrime_ / k), 1);
        }

    }

    void set(size_type n, size_type q, size_type z, size_type l) {

        if (z >= T_.size()) {
            L_[z - T_.size()] = null_;
        } else {

            auto k = (l < upperH_) ? upperK_ : lowerK_;

            if (T_[z]) {
                set(n / k, q % (n / k), (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k + q / (n / k), l + 1);
            }

        }

    }

};

template<>
class HybridRowTree<bool> : public virtual RowTree<bool> {

public:
    typedef bool elem_type;

    typedef RelationList list_type;


    HybridRowTree() {
        // nothing to do
    }

    HybridRowTree(const HybridRowTree& other) {

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

    HybridRowTree& operator=(const HybridRowTree& other) {

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

    HybridRowTree(const bit_vector_type& v, const size_type upperK, const size_type upperH, const size_type lowerK) {

        // at least one lowerK-level
        // at least zero upperK-levels but never more than upperH upperK-levels

        null_ = false;

        upperK_ = upperK;
        lowerK_ = lowerK;

        nPrime_ = size_type(ceil((1.0 * v.size()) / lowerK_));
        upperH_ = std::min(upperH, std::max((size_type)1, logK(nPrime_, upperK_)));

        nPrime_ = size_type(pow(upperK_, upperH_));
        h_ = upperH_;
        do {

            nPrime_ *= lowerK_;
            h_++;

        } while (nPrime_ < v.size());

        std::vector<std::vector<bool>> levels(h_);
        buildFromVector(v, levels, nPrime_, 1, 0);

        // helper variables (describe top section of conceptual k^2-tree) for navigation on T and L
        upperOnes_ = 0;
        upperLength_ = 0;

        if (upperH_ > 0) {

            for (size_t l = 0; l < upperH_ - 1; l++) {
                for (auto i = 0; i < levels[l].size(); i++) {
                    upperOnes_ += levels[l][i];
                }
            }
            upperLength_ = (upperOnes_ + 1) * upperK_;

        }

        size_type total = 0;
        for (auto l = 0; l < h_ - 1; l++) {
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

    HybridRowTree(const list_type& list, const size_type upperK, const size_type upperH, const size_type lowerK, const int mode) {

        // at least one lowerK-level
        // at least zero upperK-levels but never more than upperH upperK-levels

        null_ = false;

        upperK_ = upperK;
        lowerK_ = lowerK;

        size_type maxCol = 0;
        for (auto& e : list) {
            maxCol = std::max(maxCol, e);
        }
        maxCol++; // for number of columns

        nPrime_ = size_type(ceil((1.0 * maxCol) / lowerK_));
        upperH_ = std::min(upperH, std::max((size_type)1, logK(nPrime_, upperK_)));

        nPrime_ = size_type(pow(upperK_, upperH_));
        h_ = upperH_;
        do {

            nPrime_ *= lowerK_;
            h_++;

        } while (nPrime_ < maxCol);

        switch (mode) {

            case 0: {

                std::vector<std::vector<bool>> levels(h_);
                list_type::const_iterator cursor = list.begin();

                buildFromList(list, cursor, levels, nPrime_, 1, 0);

                upperOnes_ = 0;
                upperLength_ = 0;

                if (upperH_ > 0) {

                    for (size_t l = 0; l < upperH_ - 1; l++) {
                        for (auto i = 0; i < levels[l].size(); i++) {
                            upperOnes_ += levels[l][i];
                        }
                    }
                    upperLength_ = (upperOnes_ + 1) * upperK_;

                }

                size_type total = 0;
                for (auto l = 0; l < h_ - 1; l++) {
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

                buildFromListsViaTree(list);

                R_ = rank_type(&T_);

                break;

            }

            case 2: {

                buildFromListsDynamicBitmaps(list);

                break;

            }

            default: {

                buildFromListsDynamicBitmaps(list);
                break;

            }

        }

    }

    HybridRowTree(list_type& pairs, const size_type upperK, const size_type upperH, const size_type lowerK) {

        // at least one lowerK-level
        // at least zero upperK-levels but never more than upperH upperK-levels

        null_ = false;

        upperK_ = upperK;
        lowerK_ = lowerK;

        size_type maxIndex = 0;
        for (auto p : pairs) {
            maxIndex = std::max({maxIndex, p});
        }
        maxIndex++; // for number of columns

        nPrime_ = size_type(ceil((1.0 * maxIndex) / lowerK_));
        upperH_ = std::min(upperH, std::max((size_type)1, logK(nPrime_, upperK_)));

        nPrime_ = size_type(pow(upperK_, upperH_));
        h_ = upperH_;
        do {

            nPrime_ *= lowerK_;
            h_++;

        } while (nPrime_ < maxIndex);


        if (pairs.size() != 0) {
            buildFromListsInplace(pairs);
        }

        R_ = rank_type(&T_);

    }

    HybridRowTree(const std::vector<std::pair<size_type, size_type>>::iterator& first, const std::vector<std::pair<size_type, size_type>>::iterator& last, const size_type upperK, const size_type upperH, const size_type lowerK) {

        // at least one lowerK-level
        // at least zero upperK-levels but never more than upperH upperK-levels

        null_ = false;

        upperK_ = upperK;
        lowerK_ = lowerK;

        size_type maxIndex = 0;
        for (auto iter = first; iter != last; iter++) {
            maxIndex = std::max({maxIndex, iter->second});
        }
        maxIndex++; // for number of columns

        nPrime_ = size_type(ceil((1.0 * maxIndex) / lowerK_));
        upperH_ = std::min(upperH, std::max((size_type)1, logK(nPrime_, upperK_)));

        nPrime_ = size_type(pow(upperK_, upperH_));
        h_ = upperH_;
        do {

            nPrime_ *= lowerK_;
            h_++;

        } while (nPrime_ < maxIndex);


        if (last - first != 0) {
            buildFromListsInplace(first, last);
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

    size_type getLength() override {
        return nPrime_;
    }

    elem_type getNull() override {
        return null_;
    }


    bool isNotNull(size_type i) override {
        return checkInit(i);
    }

    elem_type getElement(size_type i) override {
        return isNotNull(i);
    }

    std::vector<elem_type> getElementsInRange(size_type l, size_type r) override {
        return std::vector<elem_type>(getPositionsInRange(l, r).size(), 1);
    }

    std::vector<size_type> getPositionsInRange(size_type l, size_type r) override {

        std::vector<size_type> positions;
        rangeInit(positions, l, r);
//        rangePosInit(positions, l, std::min(r, nPrime_ - 1));

        return positions;

    }

    std::vector<std::pair<size_type, elem_type>> getValuedPositionsInRange(size_type l, size_type r) override {

        std::vector<size_type> positions;
        rangeInit(positions, l, r);
//        rangePosInit(positions, l, std::min(r, nPrime_ - 1));

        std::vector<std::pair<size_type, elem_type>> res;
        for (auto p : positions) {
            res.push_back(std::make_pair(p, 1));
        }

        return res;

    }

    std::vector<elem_type> getAllElements() override {
        return std::vector<elem_type>(countElements(), 1);
    }

    std::vector<size_type> getAllPositions() override {
//        return getPositionsInRange(0, nPrime_ - 1);
        std::vector<size_type> positions;
        fullRangeIterative(positions);
        return positions;

    }

    std::vector<std::pair<size_type, elem_type>> getAllValuedPositions() override {
//        return getValuedPositionsInRange(0, nPrime_ - 1);
        std::vector<std::pair<size_type, elem_type>> positions;
        fullRangeValPosIterative(positions);
        return positions;

    }

    bool containsElement(size_type l, size_type r) override {
        return elemInRangeInit(l, r);
//        return elemInRangeInit(l, std::min(r, nPrime_ - 1));
    }

    size_type countElements() override {

        size_type res = 0;
        for (auto i = 0; i < L_.size(); i++) {
            res += L_[i];
        }

        return res;

    }


    HybridRowTree* clone() const override {
        return new HybridRowTree<elem_type>(*this);
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
            for (auto i = 0; i < T_.size(); i++) std::cout << T_[i];
            std::cout << std::endl << std::endl;

            std::cout << "### L ###" << std::endl;
            for (auto i = 0; i < L_.size(); i++) std::cout << L_[i];
            std::cout << std::endl << std::endl;

            std::cout << "### R ###" << std::endl;
            printRanks(R_);
            std::cout << std::endl;

        }

    }


    // note: can "invalidate" the data structure (containsLink() probably won't work correctly afterwards)
    void setNull(size_type i) override {
        setInit(i);
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


    /* helper method for construction from vector */

    bool buildFromVector(const bit_vector_type& v, std::vector<std::vector<bool>>& levels, size_type numCols, size_type l, size_type q) {

        std::vector<bool> C;
        auto k = (l <= upperH_) ? upperK_ : lowerK_;

        if (l == h_) {

            for (auto j = 0; j < k; j++) {
                C.push_back(((q + j) < v.size()) ? v[q + j] : false);
            }

        } else {

            for (auto j = 0; j < k; j++) {
                C.push_back(buildFromVector(v, levels, numCols / k, l + 1, q + j * (numCols / k)));
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

    bool buildFromList(const list_type& list, list_type::const_iterator& cursor, std::vector<std::vector<bool>>& levels, size_type n, size_type l, size_type q) {// 3.3.2

        std::vector<bool> C;
        auto k = (l <= upperH_) ? upperK_ : lowerK_;

        if (l == h_) {

            for (auto j = 0; j < k; j++) {

                C.push_back((cursor != list.end()) && ((q + j) == *cursor));
                if (C.back()) cursor++;

            }

        } else {

            for (auto j = 0; j < k; j++) {
                C.push_back(buildFromList(list, cursor, levels, n / k, l + 1, q + j * (n / k)));
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

    void buildFromListsViaTree(const list_type& list) {// 3.3.3, so far without special bit vectors without initialisation

        Node<bool>* root = new Node<bool>(false);

        for (auto j = 0; j < list.size(); j++) {
            insert(root, nPrime_, list[j], 0);
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

                for (auto i = 0; i < k; i++) {

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

                upperLength_ = (upperH_ > 0) ? (upperOnes_ + 1) * upperK_ : 0;

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

    void insert(Node<bool>* node, size_type n, size_type q, size_type l) {

        auto k = (l < upperH_) ? upperK_ : lowerK_;

        if (n == k) {

            if (node->isLeaf()) {
                node->turnInternal(k, true);
            }

            node->addChild(q, true);

        } else {

            if (node->isLeaf()) {
                node->turnInternal(k, false);
            }

            size_type z = q / (n / k);

            insert(node->hasChild(z) ? node->getChild(z) : node->addChild(z, false), n / k, q % (n / k), l + 1);

        }

    }

    /* helper methods for construction from relation lists via dynamic bitmap representations */

    void buildFromListsDynamicBitmaps(const list_type& list) {// 3.3.4, currently no succinct dynamic bitmaps

        if (h_ == 1) {

            L_ = bit_vector_type(lowerK_, 0);

            for (auto j = 0; j < list.size(); j++) {
                L_[list[j]] = 1;
            }

            if (isAllZero(L_)) {
                L_ = bit_vector_type(0);
            }

        } else {

            std::vector<bool> T, L;
            NaiveDynamicRank R;

            for (auto j = 0; j < list.size(); j++) {
                insertInit(T, L, R, list[j]);
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

    void insertInit(std::vector<bool>& T, std::vector<bool>& L, NaiveDynamicRank& R, size_type q) {

        auto k = (upperH_ > 0) ? upperK_ : lowerK_;

        if (T.empty()) {

            T = std::vector<bool>(k);
            R = NaiveDynamicRank(T);

            upperOnes_ = 0;
            upperLength_ = (upperH_ > 0) ? k : 0;

        }

        insert(T, L, R, nPrime_ / k, q % (nPrime_ / k), q / (nPrime_ / k), 1);

    }

    void insert(std::vector<bool>& T, std::vector<bool>& L, NaiveDynamicRank& R, size_type n, size_type q, size_type z, size_type l) {

        auto k = (l < upperH_) ? upperK_ : lowerK_;

        if (!T[z]) {

            T[z] = 1;
            R.increaseFrom(z + 1);

            if (l < upperH_) {

                upperOnes_++;
                upperLength_ += k;

            }

            size_type y = (l >= upperH_) * upperLength_ + (R.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k;

            if ((l + 1) == h_) {

                L.insert(L.begin() + y - T.size(), k, 0);
                L[y + q / (n / k) - T.size()] = 1;

            } else {

                T.insert(T.begin() + y, k, 0);
                R.insert(y + 1, k);

                insert(T, L, R, n / k, q % (n / k), y + q / (n / k), l + 1);

            }

        } else {

            size_type y = (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k + q / (n / k);

            if ((l + 1) == h_) {
                L[y - T.size()] = 1;
            } else {
                insert(T, L, R, n / k, q % (n / k), y, l + 1);
            }

        }

    }

    /* helper methods for inplace construction from single list of pairs */

    size_type computeKey(const list_type::value_type& pair, const Subproblem& sp, size_type width) {
        return (pair - sp.firstCol) / width;
    }

    void countingSort(list_type& pairs, std::vector<std::pair<size_type, size_type>>& intervals, const Subproblem& sp, size_type width, size_type sup) {

        std::vector<size_type> counts(sup);

        // determine key frequencies
        for (auto i = sp.left; i < sp.right; i++) {
            counts[computeKey(pairs[i], sp, width)]++;
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
        list_type tmpPairs(sp.right - sp.left + 1);
        for (auto i = sp.left; i < sp.right; i++) {

            tmpPairs[counts[computeKey(pairs[i], sp, width)]] = pairs[i];
            counts[computeKey(pairs[i], sp, width)]++;

        }

        for (auto i = sp.left; i < sp.right; i++) {
            pairs[i] = tmpPairs[i - sp.left];
        }

    }

    void buildFromListsInplace(list_type& pairs) {// 3.3.5

        std::queue<std::pair<Subproblem, size_type>> queue;
        std::pair<Subproblem, size_type> sp;
        size_type k, S;
        std::vector<std::pair<size_type, size_type>> intervals(std::max(upperK_, lowerK_) * std::max(upperK_, lowerK_));
        std::vector<bool> T, L;
        bit_vector_type appToL;

        upperOnes_ = 0;
        upperLength_ = 0;

        queue.push(std::make_pair(Subproblem(0, 0, 0, nPrime_ - 1, 0, pairs.size()), 0));

        while (!queue.empty()) {

            sp = queue.front();
            queue.pop();

            k = (sp.second < upperH_) ? upperK_ : lowerK_;
            S = sp.first.lastCol - sp.first.firstCol + 1;

            if (S > k) {

                countingSort(pairs, intervals, sp.first, S / k, k);

                for (auto i = 0; i < k; i++) {

                    if (intervals[i].first < intervals[i].second) {

                        T.push_back(true);
                        queue.push(std::make_pair(
                                Subproblem(
                                        0,
                                        0,
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

                appToL = bit_vector_type(k);

                for (auto i = sp.first.left; i < sp.first.right; i++) {
                    appToL[pairs[i] - sp.first.firstCol] = true;
                }

                L.insert(L.end(), appToL.begin(), appToL.end());

            }

        }

        upperLength_ = (upperH_ > 0) ? (upperOnes_ + 1) * upperK_ : 0;

        L_ = bit_vector_type(L.size());
        std::move(L.begin(), L.end(), L_.begin());
        L.clear();
        L.shrink_to_fit();

        T_ = bit_vector_type(T.size());
        std::move(T.begin(), T.end(), T_.begin());

    }

    /* helper methods for inplace construction from single list of pairs given as iterators */

    size_type computeKey(const std::vector<std::pair<size_type, size_type>>::iterator::value_type& pair, const Subproblem& sp, size_type width) {
        return (pair.second - sp.firstCol) / width;
    }

    void countingSort(const std::vector<std::pair<size_type, size_type>>::iterator& first, std::vector<std::pair<size_type, size_type>>& intervals, const Subproblem& sp, size_type width, size_type sup) {

        std::vector<size_type> counts(sup);

        // determine key frequencies
        for (auto i = sp.left; i < sp.right; i++) {
            counts[computeKey(first[i], sp, width)]++;
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
        std::vector<std::pair<size_type, size_type>> tmpPairs(sp.right - sp.left + 1);
        for (auto i = sp.left; i < sp.right; i++) {

            tmpPairs[counts[computeKey(first[i], sp, width)]] = first[i];
            counts[computeKey(first[i], sp, width)]++;

        }

        for (auto i = sp.left; i < sp.right; i++) {
            first[i] = tmpPairs[i - sp.left];
        }

    }

    void buildFromListsInplace(const std::vector<std::pair<size_type, size_type>>::iterator& first, const std::vector<std::pair<size_type, size_type>>::iterator& last) {// 3.3.5

        std::queue<std::pair<Subproblem, size_type>> queue;
        std::pair<Subproblem, size_type> sp;
        size_type k, S;
        std::vector<std::pair<size_type, size_type>> intervals(std::max(upperK_, lowerK_) * std::max(upperK_, lowerK_));
        std::vector<bool> T, L;
        bit_vector_type appToL;

        upperOnes_ = 0;
        upperLength_ = 0;

        queue.push(std::make_pair(Subproblem(0, 0, 0, nPrime_ - 1, 0, last - first), 0));

        while (!queue.empty()) {

            sp = queue.front();
            queue.pop();

            k = (sp.second < upperH_) ? upperK_ : lowerK_;
            S = sp.first.lastCol - sp.first.firstCol + 1;

            if (S > k) {

                countingSort(first, intervals, sp.first, S / k, k);

                for (auto i = 0; i < k; i++) {

                    if (intervals[i].first < intervals[i].second) {

                        T.push_back(true);
                        queue.push(std::make_pair(
                                Subproblem(
                                        0,
                                        0,
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

                appToL = bit_vector_type(k);

                for (auto i = sp.first.left; i < sp.first.right; i++) {
                    appToL[first[i].second - sp.first.firstCol] = true;
                }

                L.insert(L.end(), appToL.begin(), appToL.end());

            }

        }

        upperLength_ = (upperH_ > 0) ? (upperOnes_ + 1) * upperK_ : 0;

        L_ = bit_vector_type(L.size());
        std::move(L.begin(), L.end(), L_.begin());
        L.clear();
        L.shrink_to_fit();

        T_ = bit_vector_type(T.size());
        std::move(T.begin(), T.end(), T_.begin());

    }



    /* isNotNull() */

    bool checkInit(size_type q) {

        size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

        return (L_.empty()) ? false : check(nPrime_ / k, q % (nPrime_ / k), q / (nPrime_ / k), 1);

    }

    bool check(size_type n, size_type q, size_type z, size_type l) {

        if (z >= T_.size()) {
            return L_[z - T_.size()];
        } else {

            auto k = (l < upperH_) ? upperK_ : lowerK_;

            return T_[z] ? check(n / k, q % (n / k), (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k + q / (n / k), l + 1) : false;

        }

    }

    /* getRange() */

    void fullRangeIterative(std::vector<size_type>& elems) {

        if (L_.empty()) return;

        std::queue<SubrowInfo> queue, nextLevelQueue;
        size_type lenT = T_.size();
        size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

        if (lenT == 0) {

            for (size_type i = 0; i < nPrime_; i++) {
                if (L_[i]) {
                    elems.push_back(i);
                }
            }

        } else {

            // rangeInit
            size_type n = nPrime_/ k;
            size_type l = 1;

            for (size_type z = 0, dq = 0; z < k; z++, dq += n) {
                queue.push(SubrowInfo(dq, z));
            }

            // range
            k = (l < upperH_) ? upperK_ : lowerK_;
            n /= k;
            for (; n > 1; l++, k = (l < upperH_) ? upperK_ : lowerK_, n /= k) {

                size_type a = (l >= upperH_) * upperLength_;
                size_type b = (l >= upperH_) * (upperOnes_ + 1);

                while (!queue.empty()) {

                    auto& cur = queue.front();

                    if (T_[cur.z]) {

                        auto y = a + (R_.rank(cur.z + 1) - b) * k;

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

                    auto y = a + (R_.rank(cur.z + 1) - b) * k - lenT;

                    for (size_type j = 0, newDq = cur.dq; j < k; j++, newDq += n, y++) {
                        if (L_[y]) {
                            elems.push_back(newDq);
                        }
                    }

                }

                queue.pop();

            }

        }

    }

    void fullRangeValPosIterative(std::vector<std::pair<size_type, elem_type>>& elems) {

        if (L_.empty()) return;

        std::queue<SubrowInfo> queue, nextLevelQueue;
        size_type lenT = T_.size();
        size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

        if (lenT == 0) {

            for (size_type i = 0; i < nPrime_; i++) {
                if (L_[i]) {
                    elems.push_back(std::make_pair(i, 1));
                }
            }

        } else {

            // rangeInit
            size_type n = nPrime_/ k;
            size_type l = 1;

            for (size_type z = 0, dq = 0; z < k; z++, dq += n) {
                queue.push(SubrowInfo(dq, z));
            }

            // range
            k = (l < upperH_) ? upperK_ : lowerK_;
            n /= k;
            for (; n > 1; l++, k = (l < upperH_) ? upperK_ : lowerK_, n /= k) {

                size_type a = (l >= upperH_) * upperLength_;
                size_type b = (l >= upperH_) * (upperOnes_ + 1);

                while (!queue.empty()) {

                    auto& cur = queue.front();

                    if (T_[cur.z]) {

                        auto y = a + (R_.rank(cur.z + 1) - b) * k;

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

                    auto y = a + (R_.rank(cur.z + 1) - b) * k - lenT;

                    for (size_type j = 0, newDq = cur.dq; j < k; j++, newDq += n, y++) {
                        if (L_[y]) {
                            elems.push_back(std::make_pair(newDq, 1));
                        }
                    }

                }

                queue.pop();

            }

        }

    }

    void rangeInit(std::vector<size_type>& elems, size_type l, size_type r) {

        if (!L_.empty()) {

            size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

            for (auto j = l / (nPrime_ / k); j <= r / (nPrime_ / k); j++) {
                range(
                        elems,
                        nPrime_ / k,
                        (j == l / (nPrime_ / k)) * (l % (nPrime_ / k)),
                        (j == r / (nPrime_ / k)) ? r % (nPrime_ / k) : (nPrime_ / k) - 1,
                        (nPrime_ / k) * j,
                        j,
                        1
                );
            }

        }

    }

    void range(std::vector<size_type>& elems, size_type n, size_type l, size_type r, size_type dq, size_type z, size_type level) {

        if (z >= T_.size()) {

            if (L_[z - T_.size()]) {
                elems.push_back(dq);
            }

        } else {

            if (T_[z]) {

                auto k = (level < upperH_) ? upperK_ : lowerK_;
                auto y = (level >= upperH_) * upperLength_ + (R_.rank(z + 1) - (level >= upperH_) * (upperOnes_ + 1)) * k;

                for (auto j = l / (n / k); j <= r / (n / k); j++) {
                    range(
                            elems,
                            n / k,
                            (j == l / (n / k)) * (l % (n / k)),
                            (j == r / (n / k)) ? r % (n / k) : n / k - 1,
                            dq + (n / k) * j,
                            y + j,
                            level + 1
                    );
                }

            }

        }

    }

    /* linkInRange() */

    bool elemInRangeInit(size_type l, size_type r) {

        if (!L_.empty()) {

            if ((l == 0) && (r == (nPrime_ - 1))) {
                return 1;
            }

            size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

            for (auto j = l / (nPrime_ / k); j <= r / (nPrime_ / k); j++) {

                if (elemInRange(nPrime_ / k, (j == l / (nPrime_ / k)) * (l % (nPrime_ / k)), (j == r / (nPrime_ / k)) ? r % (nPrime_ / k) : nPrime_ / k - 1, j, 1)) {
                    return true;
                }

            }

        }

        return false;

    }

    bool elemInRange(size_type n, size_type l, size_type r, size_type z, size_type level) {

        if (z >= T_.size()) {

            return L_[z - T_.size()];

        } else {

            if (T_[z]) {

                if ((l == 0) && (r == (n - 1))) {
                    return true;
                }

                auto k = (level < upperH_) ? upperK_ : lowerK_;
                auto y = (level >= upperH_) * upperLength_ + (R_.rank(z + 1) - (level >= upperH_) * (upperOnes_ + 1)) * k;

                for (auto j = l / (n / k); j <= r / (n / k); j++) {

                    if (elemInRange(n / k, (j == l / (n / k)) * (l % (n / k)), (j == r / (n / k)) ? r % (n / k) : n / k - 1, y + j, level + 1)) {
                        return true;
                    }

                }

            }

            return false;

        }

    }


    /* setNull() */

    void setInit(size_type q) {

        size_type k = (upperH_ > 0) ? upperK_ : lowerK_;

        if (!L_.empty()) {
            set(nPrime_ / k, q % (nPrime_ / k), q / (nPrime_ / k), 1);
        }

    }

    void set(size_type n, size_type q, size_type z, size_type l) {

        if (z >= T_.size()) {
            L_[z - T_.size()] = null_;
        } else {

            auto k = (l < upperH_) ? upperK_ : lowerK_;

            if (T_[z]) {
                set(n / k, q % (n / k), (l >= upperH_) * upperLength_ + (R_.rank(z + 1) - (l >= upperH_) * (upperOnes_ + 1)) * k + q / (n / k), l + 1);
            }

        }

    }

};

#endif //K2TREES_STATICHYBRIDROWTREE_HPP
