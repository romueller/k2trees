#ifndef K2TREES_STATICROWTREE_HPP
#define K2TREES_STATICROWTREE_HPP

#include <queue>

#include "RowTree.hpp"
#include "Utility.hpp"

template<typename E>
class BasicRowTree : public virtual RowTree<E> {

public:
    typedef E elem_type;

    typedef typename RowTree<elem_type>::list_type list_type;


    BasicRowTree() {
        // nothing to do
    }

    BasicRowTree(const BasicRowTree& other) {

        h_ = other.h_;
        k_ = other.k_;
        nPrime_ = other.nPrime_;
        null_ = other.null_;

        T_ = other.T_;
        L_ = other.L_;
        R_ = rank_type(&T_);

    }

    BasicRowTree& operator=(const BasicRowTree& other) {

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        h_ = other.h_;
        k_ = other.k_;
        nPrime_ = other.nPrime_;
        null_ = other.null_;

        T_ = other.T_;
        L_ = other.L_;
        R_ = rank_type(&T_);

        return *this;

    }

    BasicRowTree(const std::vector<elem_type>& v, const size_type k, const elem_type null = elem_type()) {

        null_ = null;

        k_ = k;
        h_ = std::max((size_type)1, logK(v.size(), k_));
        nPrime_ = size_type(pow(k_, h_));

        std::vector<std::vector<bool>> levels(h_ - 1);
        buildFromVector(v, levels, nPrime_, 1, 0);

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

    BasicRowTree(const list_type& list, const size_type k, const int mode, const elem_type null = elem_type()) {

        null_ = null;

        size_type maxCol = 0;
        for (auto& e : list) {
            maxCol = std::max(maxCol, e.first);
        }

        k_ = k;
        h_ = std::max((size_type)1, logK(maxCol + 1, k_));
        nPrime_ = size_type(pow(k_, h_));

        switch (mode) {

            case 0: {

                std::vector<std::vector<bool>> levels(h_ - 1);
                auto cursor = list.begin();

                buildFromList(list, cursor, levels, nPrime_, 1, 0);

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

    BasicRowTree(list_type& pairs, const size_type k, const elem_type null = elem_type()) {

        null_ = null;

        size_type maxCol = 0;
        for (auto& e : pairs) {
            maxCol = std::max(maxCol, e.first);
        }

        k_ = k;
        h_ = std::max((size_type)1, logK(maxCol + 1, k_));
        nPrime_ = size_type(pow(k_, h_));

        if (pairs.size() != 0) {
            buildFromListsInplace(pairs);
        }

        R_ = rank_type(&T_);

    }


    size_type getH() {
        return h_;
    }

    size_type getK() {
        return k_;
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


    BasicRowTree* clone() const override {
        return new BasicRowTree<elem_type>(*this);
    }


    void print(bool all = false) override {

        std::cout << "### Parameters ###" << std::endl;
        std::cout << "h  = " << h_ << std::endl;
        std::cout << "k  = " << k_ << std::endl;
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

    size_type k_;
    size_type h_;
    size_type nPrime_;

    elem_type null_;


    /* helper method for construction from vector */

    bool buildFromVector(const std::vector<elem_type>& v, std::vector<std::vector<bool>>& levels, size_type numCols, size_type l, size_type q) {

        if (l == h_) {

            std::vector<elem_type> C;

            for (auto j = 0; j < k_; j++) {
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

            for (auto j = 0; j < k_; j++) {
                C.push_back(buildFromVector(v, levels, numCols / k_, l + 1, q + j * (numCols / k_)));
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

        if (l == h_) {

            std::vector<elem_type> C;

            for (auto j = 0; j < k_; j++) {

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

            for (auto j = 0; j < k_; j++) {
                C.push_back(buildFromList(list, cursor, levels, n / k_, l + 1, q + j * (n / k_)));
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
            insert(root, nPrime_, list[j].first, list[j].second);
        }

        if (!root->isLeaf()) {

            std::vector<bool> T;

            // traverse over tree and generate T and L while doing it
            std::queue<Node<elem_type>*> queue;
            Node<elem_type>* node;
            Node<elem_type>* child;
            queue.push(root);

            while (!queue.empty()) {

                node = queue.front();
                queue.pop();

                for (size_type i = 0; i < k_; i++) {

                    child = node->getChild(i);

                    if (child != 0 && child->isLeaf()) {
                        L_.push_back(child->getLabel());
                    } else {

                        T.push_back(child != 0);

                        if (T.back()) {
                            queue.push(child);
                        }

                    }

                }

            }

            T_ = bit_vector_type(T.size());
            std::move(T.begin(), T.end(), T_.begin());

        }

        delete root;

    }

    void insert(Node<elem_type>* node, size_type n, size_type q, elem_type val) {

        if (n == k_) {

            if (node->isLeaf()) {
                node->turnInternal(k_, true);
            }

            node->addChild(q, val);

        } else {

            if (node->isLeaf()) {
                node->turnInternal(k_, false);
            }

            size_type z = q / (n / k_);

            insert(node->hasChild(z) ? node->getChild(z) : node->addChild(z, null_), n / k_, q % (n / k_), val);

        }

    }

    /* helper methods for construction from relation lists via dynamic bitmap representations */

    void buildFromListsDynamicBitmaps(const list_type& list) {// 3.3.4, currently no succinct dynamic bitmaps

        if (h_ == 1) {

            L_ = std::vector<elem_type>(k_, null_);

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

        if (T.empty()) {

            T = std::vector<bool>(k_);
            R = NaiveDynamicRank(T);

        }

        insert(T, R, nPrime_ / k_, q % (nPrime_ / k_), val, q / (nPrime_ / k_), 1);

    }

    void insert(std::vector<bool>& T, NaiveDynamicRank& R, size_type n, size_type q, elem_type val, size_type z, size_type l) {

        if (!T[z]) {

            T[z] = 1;
            R.increaseFrom(z + 1);

            size_type y = R.rank(z + 1) * k_ + q / (n / k_);

            if ((l + 1) == h_) {

                L_.insert(L_.begin() + R.rank(z + 1) * k_ - T.size(), k_, null_);
                L_[y - T.size()] = val;

            } else {

                T.insert(T.begin() + R.rank(z + 1) * k_, k_, null_);
                R.insert(R.rank(z + 1) * k_ + 1, k_);

                insert(T, R, n / k_, q % (n / k_), val, y, l + 1);

            }

        } else {

            size_type y = R.rank(z) * k_ + q / (n / k_);

            if ((l + 1) == h_) {
                L_[y - T.size()] = val;
            } else {
                insert(T, R, n / k_, q % (n / k_), val, y, l + 1);
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

        std::queue<Subproblem> queue;
        Subproblem sp;
        size_type S;
        std::vector<std::pair<size_type, size_type>> intervals(k_);
        std::vector<bool> T;
        std::vector<elem_type> appToL;

        queue.push(Subproblem(0, 0, 0, nPrime_ - 1, 0, pairs.size()));

        while (!queue.empty()) {

            sp = queue.front();
            queue.pop();

            S = sp.lastCol - sp.firstCol + 1;

            if (S > k_) {

                countingSort(pairs, intervals, sp, S / k_, k_);

                for (auto i = 0; i < k_; i++) {

                    if (intervals[i].first < intervals[i].second) {

                        T.push_back(true);
                        queue.push(Subproblem(
                                0,
                                0,
                                sp.firstCol + (i % k_) * (S / k_),
                                sp.firstCol + (i % k_ + 1) * (S / k_) - 1,
                                sp.left + intervals[i].first,
                                sp.left + intervals[i].second
                        ));

                    } else {
                        T.push_back(false);
                    }

                }

            } else {

                appToL = std::vector<elem_type>(k_);

                for (auto i = sp.left; i < sp.right; i++) {
                    appToL[pairs[i].first - sp.firstCol] = pairs[i].second;
                }

                L_.insert(L_.end(), appToL.begin(), appToL.end());

            }

        }

        T_ = bit_vector_type(T.size());
        std::move(T.begin(), T.end(), T_.begin());

    }


    /* isNotNull() */

    bool checkInit(size_type q) {
        return (L_.empty()) ? false : check(nPrime_ / k_, q % (nPrime_ / k_), q / (nPrime_ / k_));
    }

    bool check(size_type n, size_type q, size_type z) {

        if (z >= T_.size()) {
            return (L_[z - T_.size()] != null_);
        } else {
            return T_[z] ? check(n / k_, q % (n / k_), R_.rank(z + 1) * k_ + q / (n / k_)) : false;
        }

    }

    /* getElement() */

    elem_type getInit(size_type q) {
        return (L_.empty()) ? null_ : get(nPrime_ / k_, q % (nPrime_ / k_), q / (nPrime_ / k_));
    }

    elem_type get(size_type n, size_type q, size_type z) {

        if (z >= T_.size()) {
            return L_[z - T_.size()];
        } else {
            return T_[z] ? get(n / k_, q % (n / k_), R_.rank(z + 1) * k_ + q / (n / k_)) : null_;
        }

    }

    /* getElementsInRange() */

    void fullRangeElemIterative(std::vector<elem_type>& elems) {

        if (L_.empty()) return;

        std::queue<SubrowInfo> queue, nextLevelQueue;
        size_type lenT = T_.size();

        if (lenT == 0) {

            for (size_type i = 0; i < nPrime_; i++) {
                if (L_[i] != null_) {
                    elems.push_back(L_[i]);
                }
            }

        } else {

            // rangeInit
            size_type n = nPrime_/ k_;

            for (size_type z = 0, dq = 0; z < k_; z++, dq += n) {
                queue.push(SubrowInfo(dq, z));
            }

            // range
            n /= k_;
            for (; n > 1; n /= k_) {

                while (!queue.empty()) {

                    auto& cur = queue.front();

                    if (T_[cur.z]) {

                        auto y = R_.rank(cur.z + 1) * k_;

                        for (size_type j = 0, newDq = cur.dq; j < k_; j++, newDq += n, y++) {
                            nextLevelQueue.push(SubrowInfo(newDq, y));
                        }

                    }

                    queue.pop();

                }

                queue.swap(nextLevelQueue);

            }

            while (!queue.empty()) {

                auto& cur = queue.front();

                if (T_[cur.z]) {

                    auto y = R_.rank(cur.z + 1) * k_ - lenT;

                    for (size_type j = 0, newDq = cur.dq; j < k_; j++, newDq += n, y++) {
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

            for (auto j = l / (nPrime_ / k_); j <= r / (nPrime_ / k_); j++) {
                rangeElem(
                        elems,
                        nPrime_ / k_,
                        (j == l / (nPrime_ / k_)) * (l % (nPrime_ / k_)),
                        (j == r / (nPrime_ / k_)) ? r % (nPrime_ / k_) : (nPrime_ / k_) - 1,
                        (nPrime_ / k_) * j,
                        j
                );
            }

        }

    }

    void rangeElem(std::vector<elem_type>& elems, size_type n, size_type l, size_type r, size_type dq, size_type z) {

        if (z >= T_.size()) {

            if (L_[z - T_.size()] != null_) {
                elems.push_back(L_[z - T_.size()]);
            }

        } else {

            if (T_[z]) {

                auto y = R_.rank(z + 1) * k_;

                for (auto j = l / (n / k_); j <= r / (n / k_); j++) {
                    rangeElem(
                            elems,
                            n / k_,
                            (j == l / (n / k_)) * (l % (n / k_)),
                            (j == r / (n / k_)) ? r % (n / k_) : n / k_ - 1,
                            dq + (n / k_) * j,
                            y + j
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

        if (lenT == 0) {

            for (size_type i = 0; i < nPrime_; i++) {
                if (L_[i] != null_) {
                    elems.push_back(i);
                }
            }

        } else {

            // rangeInit
            size_type n = nPrime_/ k_;

            for (size_type z = 0, dq = 0; z < k_; z++, dq += n) {
                queue.push(SubrowInfo(dq, z));
            }

            // range
            n /= k_;
            for (; n > 1; n /= k_) {

                while (!queue.empty()) {

                    auto& cur = queue.front();

                    if (T_[cur.z]) {

                        auto y = R_.rank(cur.z + 1) * k_;

                        for (size_type j = 0, newDq = cur.dq; j < k_; j++, newDq += n, y++) {
                            nextLevelQueue.push(SubrowInfo(newDq, y));
                        }

                    }

                    queue.pop();

                }

                queue.swap(nextLevelQueue);

            }

            while (!queue.empty()) {

                auto& cur = queue.front();

                if (T_[cur.z]) {

                    auto y = R_.rank(cur.z + 1) * k_ - lenT;

                    for (size_type j = 0, newDq = cur.dq; j < k_; j++, newDq += n, y++) {
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

            for (auto j = l / (nPrime_ / k_); j <= r / (nPrime_ / k_); j++) {
                rangePos(
                        elems,
                        nPrime_ / k_,
                        (j == l / (nPrime_ / k_)) * (l % (nPrime_ / k_)),
                        (j == r / (nPrime_ / k_)) ? r % (nPrime_ / k_) : (nPrime_ / k_) - 1,
                        (nPrime_ / k_) * j,
                        j
                );
            }

        }

    }

    void rangePos(std::vector<size_type>& elems, size_type n, size_type l, size_type r, size_type dq, size_type z) {

        if (z >= T_.size()) {

            if (L_[z - T_.size()] != null_) {
                elems.push_back(dq);
            }

        } else {

            if (T_[z]) {

                auto y = R_.rank(z + 1) * k_;

                for (auto j = l / (n / k_); j <= r / (n / k_); j++) {
                    rangePos(
                            elems,
                            n / k_,
                            (j == l / (n / k_)) * (l % (n / k_)),
                            (j == r / (n / k_)) ? r % (n / k_) : n / k_ - 1,
                            dq + (n / k_) * j,
                            y + j
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

        if (lenT == 0) {

            for (size_type i = 0; i < nPrime_; i++) {
                if (L_[i] != null_) {
                    elems.push_back(std::make_pair(i, L_[i]));
                }
            }

        } else {

            // rangeInit
            size_type n = nPrime_/ k_;

            for (size_type z = 0, dq = 0; z < k_; z++, dq += n) {
                queue.push(SubrowInfo(dq, z));
            }

            // range
            n /= k_;
            for (; n > 1; n /= k_) {

                while (!queue.empty()) {

                    auto& cur = queue.front();

                    if (T_[cur.z]) {

                        auto y = R_.rank(cur.z + 1) * k_;

                        for (size_type j = 0, newDq = cur.dq; j < k_; j++, newDq += n, y++) {
                            nextLevelQueue.push(SubrowInfo(newDq, y));
                        }

                    }

                    queue.pop();

                }

                queue.swap(nextLevelQueue);

            }

            while (!queue.empty()) {

                auto& cur = queue.front();

                if (T_[cur.z]) {

                    auto y = R_.rank(cur.z + 1) * k_ - lenT;

                    for (size_type j = 0, newDq = cur.dq; j < k_; j++, newDq += n, y++) {
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

            for (auto j = l / (nPrime_ / k_); j <= r / (nPrime_ / k_); j++) {
                rangeValPos(
                        elems,
                        nPrime_ / k_,
                        (j == l / (nPrime_ / k_)) * (l % (nPrime_ / k_)),
                        (j == r / (nPrime_ / k_)) ? r % (nPrime_ / k_) : (nPrime_ / k_) - 1,
                        (nPrime_ / k_) * j,
                        j
                );
            }

        }

    }

    void rangeValPos(list_type& elems, size_type n, size_type l, size_type r, size_type dq, size_type z) {

        if (z >= T_.size()) {

            if (L_[z - T_.size()] != null_) {
                elems.push_back(std::make_pair(dq, L_[z - T_.size()]));
            }

        } else {

            if (T_[z]) {

                auto y = R_.rank(z + 1) * k_;

                for (auto j = l / (n / k_); j <= r / (n / k_); j++) {
                    rangeValPos(
                            elems,
                            n / k_,
                            (j == l / (n / k_)) * (l % (n / k_)),
                            (j == r / (n / k_)) ? r % (n / k_) : n / k_ - 1,
                            dq + (n / k_) * j,
                            y + j
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

            for (auto j = l / (nPrime_ / k_); j <= r / (nPrime_ / k_); j++) {

                if (elemInRange(nPrime_ / k_, (j == l / (nPrime_ / k_)) * (l % (nPrime_ / k_)), (j == r / (nPrime_ / k_)) ? r % (nPrime_ / k_) : nPrime_ / k_ - 1, j)) {
                    return true;
                }

            }

        }

        return false;

    }

    bool elemInRange(size_type n, size_type l, size_type r, size_type z) {

        if (z >= T_.size()) {

            return L_[z - T_.size()] != null_;

        } else {

            if (T_[z]) {

                if ((l == 0) && (r == (n - 1))) {
                    return true;
                }

                auto y = R_.rank(z + 1) * k_;
                size_type p1Prime, p2Prime;

                for (auto j = l / (n / k_); j <= r / (n / k_); j++) {

                    if (elemInRange(n / k_, (j == l / (n / k_)) * (l % (n / k_)), (j == r / (n / k_)) ? r % (n / k_) : n / k_ - 1, y + j)) {
                        return true;
                    }

                }

            }

            return false;

        }

    }


    /* getElement() */

    void setInit(size_type q) {

        if (!L_.empty()) {
            set(nPrime_ / k_, q % (nPrime_ / k_), q / (nPrime_ / k_));
        }

    }

    void set(size_type n, size_type q, size_type z) {

        if (z >= T_.size()) {
            L_[z - T_.size()] = null_;
        } else {
            if (T_[z]) {
                set(n / k_, q % (n / k_), R_.rank(z + 1) * k_ + q / (n / k_));
            }
        }

    }

};

template<>
class BasicRowTree<bool> : public virtual RowTree<bool> {

public:
    typedef bool elem_type;

    typedef RelationList list_type;


    BasicRowTree() {
        // nothing to do
    }

    BasicRowTree(const BasicRowTree& other) {

        h_ = other.h_;
        k_ = other.k_;
        nPrime_ = other.nPrime_;
        null_ = other.null_;

        T_ = other.T_;
        L_ = other.L_;
        R_ = rank_type(&T_);

    }

    BasicRowTree& operator=(const BasicRowTree& other) {

        // check for self-assignment
        if (&other == this) {
            return *this;
        }

        h_ = other.h_;
        k_ = other.k_;
        nPrime_ = other.nPrime_;
        null_ = other.null_;

        T_ = other.T_;
        L_ = other.L_;
        R_ = rank_type(&T_);

        return *this;

    }

    BasicRowTree(const bit_vector_type& v, const size_type k) {

        null_ = false;

        k_ = k;
        h_ = std::max((size_type)1, logK(v.size(), k_));
        nPrime_ = size_type(pow(k_, h_));

        std::vector<std::vector<bool>> levels(h_);
        buildFromVector(v, levels, nPrime_, 1, 0);

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

    BasicRowTree(const list_type& list, const size_type k, const int mode) {

        null_ = false;

        size_type maxCol = 0;
        for (auto& e : list) {
            maxCol = std::max(maxCol, e);
        }

        k_ = k;
        h_ = std::max((size_type)1, logK(maxCol + 1, k_));
        nPrime_ = size_type(pow(k_, h_));

        switch (mode) {

            case 0: {

                std::vector<std::vector<bool>> levels(h_);
                list_type::const_iterator cursor = list.begin();

                buildFromList(list, cursor, levels, nPrime_, 1, 0);

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

    BasicRowTree(list_type& pairs, const size_type k) {

        null_ = false;

        size_type maxIndex = 0;
        for (auto p : pairs) {
            maxIndex = std::max({maxIndex, p});
        }

        k_ = k;
        h_ = std::max((size_type)1, logK(maxIndex + 1, k_));
        nPrime_ = size_type(pow(k_, h_));

        if (pairs.size() != 0) {
            buildFromListsInplace(pairs);
        }

        R_ = rank_type(&T_);

    }

    BasicRowTree(const std::vector<std::pair<size_type, size_type>>::iterator& first, const std::vector<std::pair<size_type, size_type>>::iterator& last, const size_type k) {

        null_ = false;

        size_type maxCol = 0;
        for (auto iter = first; iter != last; iter++) {
            maxCol = std::max(maxCol, iter->second);
        }

        k_ = k;
        h_ = std::max((size_type)1, logK(maxCol + 1, k_));
        nPrime_ = size_type(pow(k_, h_));

        if (last - first != 0) {
            buildFromListsInplace(first, last);
        }

        R_ = rank_type(&T_);

    }


    size_type getH() {
        return h_;
    }

    size_type getK() {
        return k_;
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


    BasicRowTree* clone() const override {
        return new BasicRowTree<elem_type>(*this);
    }


    void print(bool all = false) override {

        std::cout << "### Parameters ###" << std::endl;
        std::cout << "h  = " << h_ << std::endl;
        std::cout << "k  = " << k_ << std::endl;
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

    size_type k_;
    size_type h_;
    size_type nPrime_;

    elem_type null_;

    /* helper method for construction from vector */

    bool buildFromVector(const bit_vector_type& v, std::vector<std::vector<bool>>& levels, size_type numCols, size_type l, size_type q) {

        std::vector<bool> C;

        if (l == h_) {

            for (auto j = 0; j < k_; j++) {
                C.push_back(((q + j) < v.size()) ? v[q + j] : false);
            }

        } else {

            for (auto j = 0; j < k_; j++) {
                C.push_back(buildFromVector(v, levels, numCols / k_, l + 1, q + j * (numCols / k_)));
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

        if (l == h_) {

            for (auto j = 0; j < k_; j++) {

                C.push_back((cursor != list.end()) && ((q + j) == *cursor));
                if (C.back()) cursor++;

            }

        } else {

            for (auto j = 0; j < k_; j++) {
                C.push_back(buildFromList(list, cursor, levels, n / k_, l + 1, q + j * (n / k_)));
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
            insert(root, nPrime_, list[j]);
        }

        if (!root->isLeaf()) {

            std::vector<bool> T, L;

            // traverse over tree and generate T and L while doing it
            std::queue<Node<bool>*> queue;
            Node<bool>* node;
            Node<bool>* child;
            queue.push(root);

            while (!queue.empty()) {

                node = queue.front();
                queue.pop();

                for (auto i = 0; i < k_; i++) {

                    child = node->getChild(i);

                    if (child != 0 && child->isLeaf()) {
                        L.push_back(child->getLabel());
                    } else {

                        T.push_back(child != 0);

                        if (T.back()) {
                            queue.push(child);
                        }

                    }

                }

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

    void insert(Node<bool>* node, size_type n, size_type q) {

        if (n == k_) {

            if (node->isLeaf()) {
                node->turnInternal(k_, true);
            }

            node->addChild(q, true);

        } else {

            if (node->isLeaf()) {
                node->turnInternal(k_, false);
            }

            size_type z = q / (n / k_);

            insert(node->hasChild(z) ? node->getChild(z) : node->addChild(z, false), n / k_, q % (n / k_));

        }

    }

    /* helper methods for construction from relation lists via dynamic bitmap representations */

    void buildFromListsDynamicBitmaps(const list_type& list) {// 3.3.4, currently no succinct dynamic bitmaps

        if (h_ == 1) {

            L_ = bit_vector_type(k_, 0);

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

        if (T.empty()) {

            T = std::vector<bool>(k_);
            R = NaiveDynamicRank(T);

        }

        insert(T, L, R, nPrime_ / k_, q % (nPrime_ / k_), q / (nPrime_ / k_), 1);

    }

    void insert(std::vector<bool>& T, std::vector<bool>& L, NaiveDynamicRank& R, size_type n, size_type q, size_type z, size_type l) {

        if (!T[z]) {

            T[z] = 1;
            R.increaseFrom(z + 1);

            size_type y = R.rank(z + 1) * k_ + q / (n / k_);

            if ((l + 1) == h_) {

                L.insert(L.begin() + R.rank(z + 1) * k_ - T.size(), k_, 0);
                L[y - T.size()] = 1;

            } else {

                T.insert(T.begin() + R.rank(z + 1) * k_, k_, 0);
                R.insert(R.rank(z + 1) * k_ + 1, k_);

                insert(T, L, R, n / k_, q % (n / k_), y, l + 1);

            }

        } else {

            size_type y = R.rank(z + 1) * k_ + q / (n / k_);

            if ((l + 1) == h_) {
                L[y - T.size()] = 1;
            } else {
                insert(T, L, R, n / k_, q % (n / k_), y, l + 1);
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

        std::queue<Subproblem> queue;
        Subproblem sp;
        size_type S;
        std::vector<std::pair<size_type, size_type>> intervals(k_ * k_);
        std::vector<bool> T, L;
        bit_vector_type appToL;

        queue.push(Subproblem(0, 0, 0, nPrime_ - 1, 0, pairs.size()));

        while (!queue.empty()) {

            sp = queue.front();
            queue.pop();

            S = sp.lastCol - sp.firstCol + 1;

            if (S > k_) {

                countingSort(pairs, intervals, sp, S / k_, k_);

                for (auto i = 0; i < k_; i++) {

                    if (intervals[i].first < intervals[i].second) {

                        T.push_back(true);
                        queue.push(Subproblem(
                                0,
                                0,
                                sp.firstCol + (i % k_) * (S / k_),
                                sp.firstCol + (i % k_ + 1) * (S / k_) - 1,
                                sp.left + intervals[i].first,
                                sp.left + intervals[i].second
                        ));

                    } else {
                        T.push_back(false);
                    }

                }

            } else {

                appToL = bit_vector_type(k_);

                for (auto i = sp.left; i < sp.right; i++) {
                    appToL[pairs[i] - sp.firstCol] = true;
                }

                L.insert(L.end(), appToL.begin(), appToL.end());

            }

        }

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

        std::queue<Subproblem> queue;
        Subproblem sp;
        size_type S;
        std::vector<std::pair<size_type, size_type>> intervals(k_);
        std::vector<bool> T, L;
        std::vector<elem_type> appToL;

        queue.push(Subproblem(0, 0, 0, nPrime_ - 1, 0, last - first));

        while (!queue.empty()) {

            sp = queue.front();
            queue.pop();

            S = sp.lastCol - sp.firstCol + 1;

            if (S > k_) {

                countingSort(first, intervals, sp, S / k_, k_);

                for (auto i = 0; i < k_; i++) {

                    if (intervals[i].first < intervals[i].second) {

                        T.push_back(true);
                        queue.push(Subproblem(
                                0,
                                0,
                                sp.firstCol + (i % k_) * (S / k_),
                                sp.firstCol + (i % k_ + 1) * (S / k_) - 1,
                                sp.left + intervals[i].first,
                                sp.left + intervals[i].second
                        ));

                    } else {
                        T.push_back(false);
                    }

                }

            } else {

                appToL = std::vector<elem_type>(k_);

                for (auto i = sp.left; i < sp.right; i++) {
                    appToL[first[i].second - sp.firstCol] = true;
                }

                L.insert(L.end(), appToL.begin(), appToL.end());

            }

        }

        L_ = bit_vector_type(L.size());
        std::move(L.begin(), L.end(), L_.begin());
        L.clear();
        L.shrink_to_fit();

        T_ = bit_vector_type(T.size());
        std::move(T.begin(), T.end(), T_.begin());

    }


    /* isNotNull() */

    bool checkInit(size_type q) {
        return (L_.empty()) ? false : check(nPrime_ / k_, q % (nPrime_ / k_), q / (nPrime_ / k_));
    }

    bool check(size_type n, size_type q, size_type z) {

        if (z >= T_.size()) {
            return L_[z - T_.size()];
        } else {
            return T_[z] ? check(n / k_, q % (n / k_), R_.rank(z + 1) * k_ + q / (n / k_)) : false;
        }

    }

    /* getRange() */

    void fullRangeIterative(std::vector<size_type>& elems) {

        if (L_.empty()) return;

        std::queue<SubrowInfo> queue, nextLevelQueue;
        size_type lenT = T_.size();

        if (lenT == 0) {

            for (size_type i = 0; i < nPrime_; i++) {
                if (L_[i]) {
                    elems.push_back(i);
                }
            }

        } else {

            // rangeInit
            size_type n = nPrime_/ k_;

            for (size_type z = 0, dq = 0; z < k_; z++, dq += n) {
                queue.push(SubrowInfo(dq, z));
            }

            // range
            n /= k_;
            for (; n > 1; n /= k_) {

                while (!queue.empty()) {

                    auto& cur = queue.front();

                    if (T_[cur.z]) {

                        auto y = R_.rank(cur.z + 1) * k_;

                        for (size_type j = 0, newDq = cur.dq; j < k_; j++, newDq += n, y++) {
                            nextLevelQueue.push(SubrowInfo(newDq, y));
                        }

                    }

                    queue.pop();

                }

                queue.swap(nextLevelQueue);

            }

            while (!queue.empty()) {

                auto& cur = queue.front();

                if (T_[cur.z]) {

                    auto y = R_.rank(cur.z + 1) * k_ - lenT;

                    for (size_type j = 0, newDq = cur.dq; j < k_; j++, newDq += n, y++) {
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

        if (lenT == 0) {

            for (size_type i = 0; i < nPrime_; i++) {
                if (L_[i]) {
                    elems.push_back(std::make_pair(i, 1));
                }
            }

        } else {

            // rangeInit
            size_type n = nPrime_/ k_;

            for (size_type z = 0, dq = 0; z < k_; z++, dq += n) {
                queue.push(SubrowInfo(dq, z));
            }

            // range
            n /= k_;
            for (; n > 1; n /= k_) {

                while (!queue.empty()) {

                    auto& cur = queue.front();

                    if (T_[cur.z]) {

                        auto y = R_.rank(cur.z + 1) * k_;

                        for (size_type j = 0, newDq = cur.dq; j < k_; j++, newDq += n, y++) {
                            nextLevelQueue.push(SubrowInfo(newDq, y));
                        }

                    }

                    queue.pop();

                }

                queue.swap(nextLevelQueue);

            }

            while (!queue.empty()) {

                auto& cur = queue.front();

                if (T_[cur.z]) {

                    auto y = R_.rank(cur.z + 1) * k_ - lenT;

                    for (size_type j = 0, newDq = cur.dq; j < k_; j++, newDq += n, y++) {
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

            for (auto j = l / (nPrime_ / k_); j <= r / (nPrime_ / k_); j++) {
                range(
                        elems,
                        nPrime_ / k_,
                        (j == l / (nPrime_ / k_)) * (l % (nPrime_ / k_)),
                        (j == r / (nPrime_ / k_)) ? r % (nPrime_ / k_) : (nPrime_ / k_) - 1,
                        (nPrime_ / k_) * j,
                        j
                );
            }

        }

    }

    void range(std::vector<size_type>& elems, size_type n, size_type l, size_type r, size_type dq, size_type z) {

        if (z >= T_.size()) {

            if (L_[z - T_.size()]) {
                elems.push_back(dq);
            }

        } else {

            if (T_[z]) {

                auto y = R_.rank(z + 1) * k_;

                for (auto j = l / (n / k_); j <= r / (n / k_); j++) {
                    range(
                            elems,
                            n / k_,
                            (j == l / (n / k_)) * (l % (n / k_)),
                            (j == r / (n / k_)) ? r % (n / k_) : n / k_ - 1,
                            dq + (n / k_) * j,
                            y + j
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

            for (auto j = l / (nPrime_ / k_); j <= r / (nPrime_ / k_); j++) {

                if (elemInRange(nPrime_ / k_, (j == l / (nPrime_ / k_)) * (l % (nPrime_ / k_)), (j == r / (nPrime_ / k_)) ? r % (nPrime_ / k_) : nPrime_ / k_ - 1, j)) {
                    return true;
                }

            }

        }

        return false;

    }

    bool elemInRange(size_type n, size_type l, size_type r, size_type z) {

        if (z >= T_.size()) {

            return L_[z - T_.size()];

        } else {

            if (T_[z]) {

                if ((l == 0) && (r == (n - 1))) {
                    return true;
                }

                auto y = R_.rank(z + 1) * k_;
                size_type p1Prime, p2Prime;

                for (auto j = l / (n / k_); j <= r / (n / k_); j++) {

                    if (elemInRange(n / k_, (j == l / (n / k_)) * (l % (n / k_)), (j == r / (n / k_)) ? r % (n / k_) : n / k_ - 1, y + j)) {
                        return true;
                    }

                }

            }

            return false;

        }

    }


    /* getElement() */

    void setInit(size_type q) {

        if (!L_.empty()) {
            set(nPrime_ / k_, q % (nPrime_ / k_), q / (nPrime_ / k_));
        }

    }

    void set(size_type n, size_type q, size_type z) {

        if (z >= T_.size()) {
            L_[z - T_.size()] = null_;
        } else {
            if (T_[z]) {
                set(n / k_, q % (n / k_), R_.rank(z + 1) * k_ + q / (n / k_));
            }
        }

    }

};

#endif //K2TREES_STATICROWTREE_HPP
