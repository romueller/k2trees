#ifndef K2TREES_UTILITY_HPP
#define K2TREES_UTILITY_HPP

#include <algorithm>
#include <cmath>
#include <vector>

#include <sdsl/rank_support_v.hpp>

typedef unsigned long size_type;

typedef sdsl::bit_vector bit_vector_type;
typedef sdsl::rank_support_v<> rank_type;

template<typename T>
struct ValuedPosition {

    ValuedPosition() {
        // nothing to do
    }

    ValuedPosition(size_type r, size_type c, T v) {

        row = r;
        col = c;
        val = v;

    }

    ValuedPosition(const std::pair<size_type, size_type>& pos, T v) {

        row = pos.first;
        col = pos.second;
        val = v;

    }

    bool operator==(const ValuedPosition& other) {
        return (row == other.row) && (col == other.col) && (val == other.val);
    }

    size_type row;
    size_type col;
    T val;

};

struct PartitionIndices {

    size_type partition;
    size_type row;
    size_type col;

    PartitionIndices() {
        // nothing to do
    }

    PartitionIndices(size_type p, size_type r, size_type c) {

        partition = p;
        row = r;
        col = c;

    }

};


// parameters handed over in iterative versions of getting all positions in a row
struct SubrowInfo {

    size_type dq;
    size_type z;

    SubrowInfo(size_type dqq, size_type zz) {

        dq = dqq;
        z = zz;

    }

};

// parameters handed over in iterative versions of getting first positions in a row
struct ExtendedSubrowInfo {

    size_type nr;
    size_type nc;
    size_type p;
    size_type dq;
    size_type z;
    size_type j;

    ExtendedSubrowInfo(size_type nrr, size_type ncc, size_type pp, size_type dqq, size_type zz, size_type jj) {

        nr = nrr;
        nc = ncc;
        p = pp;
        dq = dqq;
        z = zz;
        j = jj;

    }

};


size_type logK(const size_type n, const size_type b);

template<typename T>
bool isAll(const std::vector<T>& v, const T val) {
    return std::all_of(v.begin(), v.end(), [&val](const T& elem) {return val == elem;});
}

bool isAllZero(const std::vector<bool>& v);
bool isAllZero(const bit_vector_type& v);

void printRanks(const rank_type& r);



/* Data structures for representing a relation R = A x B & conversion methods between them */

// Rectangular binary matrix (mat[i][j] == true iff (i,j) in R)
typedef std::vector<std::vector<bool>> RelationMatrix;

// List of lists, one row for each element from A; lists[i] empty iff there is no (i,_) in R, each row is sorted (ascending)
typedef std::vector<size_type> RelationList;
typedef std::vector<RelationList> RelationLists;

// List of (i,j)-pairs; contains (i,j) iff (i,j) in R
typedef std::vector<std::pair<size_type, size_type>> RelationPairs;

// ===== From matrix ... =====

template<typename T>
std::vector<std::vector<std::pair<size_type, T>>> matrixToLists(std::vector<std::vector<T>>& mat, T null = T()) {

    std::vector<std::vector<std::pair<size_type, T>>> lists;

    for (size_type i = 0; i < mat.size(); i++) {

        lists.push_back(std::vector<std::pair<size_type, T>>());

        for (size_type j = 0; j < mat[i].size(); j++) {
            if (mat[i][j] != null) {
                lists[i].push_back(std::make_pair(j, mat[i][j]));
            }
        }
    }

    return lists;

}

template<typename T>
std::vector<ValuedPosition<T>> matrixToPairs(std::vector<std::vector<T>>& mat, T null = T()) {

    std::vector<ValuedPosition<T>> pairs;

    for (size_type i = 0; i < mat.size(); i++) {
        for (size_type j = 0; j < mat[i].size(); j++) {

            if (mat[i][j] != null) {
                pairs.push_back(ValuedPosition<T>(i, j, mat[i][j]));
            }

        }
    }

    return pairs;

}

RelationLists boolMatrixToLists(RelationMatrix& mat);

RelationPairs boolMatrixToPairs(RelationMatrix& mat);

// ===== From lists ... =====

// assumes that all pairs in the lists fit into a matrix of size numRows x numCols
template<typename T>
std::vector<std::vector<T>> listsToMatrix(std::vector<std::vector<std::pair<size_type, T>>>& lists, size_type numRows, size_type numCols, T null = T()) {

    std::vector<std::vector<T>> mat(numRows, std::vector<T>(numCols, null));

    for (size_type i = 0; i < lists.size(); i++) {
        for (size_type j = 0; j < lists[i].size(); j++) {
            mat[i][lists[i][j].first] = lists[i][j].second;
        }
    }

    return mat;

}

template<typename T>
std::vector<std::vector<T>> listsToMatrix(std::vector<std::vector<std::pair<size_type, T>>>& lists, T null = T()) {

    size_type maxCol = 0;

    for (auto& l : lists) {
        for (size_type j = 0; j < l.size(); j++) {
            maxCol = std::max(maxCol, l[j].first);
        }
    }

    return listsToMatrix(lists, lists.size(), maxCol + 1);

}

template<typename T>
std::vector<ValuedPosition<T>> listsToPairs(std::vector<std::vector<std::pair<size_type, T>>>& lists, T null = T()) {

    std::vector<ValuedPosition<T>> pairs;

    for (size_type i = 0; i < lists.size(); i++) {
        for (size_type j = 0; j < lists[i].size(); j++) {
            pairs.push_back(ValuedPosition<T>(i, lists[i][j].first, lists[i][j].second));
        }
    }

    return pairs;

}

// assumes that all pairs in the lists fit into a matrix of size numRows x numCols
RelationMatrix boolListsToMatrix(RelationLists& lists, size_type numRows, size_type numCols);

RelationMatrix boolListsToMatrix(RelationLists& lists);

RelationPairs boolListsToPairs(RelationLists& lists);

// ===== From pairs ... =====

template<typename T>
std::vector<std::vector<T>> pairsToMatrix(std::vector<ValuedPosition<T>>& pairs, size_type numRows, size_type numCols, T null = T()) {

    std::vector<std::vector<T>> mat(numRows, std::vector<T>(numCols, null));

    for (auto& p : pairs) {
        mat[p.row][p.col] = p.val;
    }

    return mat;

}

template<typename T>
std::vector<std::vector<T>> pairsToMatrix(std::vector<ValuedPosition<T>>& pairs, T null = T()) {

    size_type maxRow = 0;
    size_type maxCol = 0;

    for (auto& p : pairs) {

        maxRow = std::max(maxRow, p.row);
        maxCol = std::max(maxCol, p.col);

    }

    return pairsToMatrix(pairs, maxRow + 1, maxCol + 1, null);

}

template<typename T>
std::vector<std::vector<std::pair<size_type, T>>> pairsToList(std::vector<ValuedPosition<T>>& pairs, size_type numRows, T null = T()) {

    std::vector<std::vector<std::pair<size_type, T>>> lists(numRows);

    for (auto& p : pairs) {
        lists[p.row].push_back(std::make_pair(p.col, p.val));
    }

    return lists;

}

template<typename T>
std::vector<std::vector<std::pair<size_type, T>>> pairsToList(std::vector<ValuedPosition<T>>& pairs, T null = T()) {

    size_type maxRow = 0;

    for (auto& p : pairs) {
        maxRow = std::max(maxRow, p.row);
    }

    return pairsToList(pairs, maxRow + 1);

}

// assumes that all pairs in the list fit into a matrix of size numRows x numCols
RelationMatrix boolPairsToMatrix(RelationPairs& pairs, size_type numRows, size_type numCols);

RelationMatrix boolPairsToMatrix(RelationPairs& pairs);

// assumes that all pairs in the list fit into a matrix with numRows rows
RelationLists boolPairsToList(RelationPairs& pairs, size_type numRows);

RelationLists boolPairsToList(RelationPairs& pairs);


template<typename F, typename S>
struct sortPairs {
    bool operator()(const std::pair<F, S>& a, const std::pair<F, S>& b) {
        return (a.first < b.first) || ((a.first == b.first) && (a.second < b.second));
    }
};

template<typename T>
struct sortValuedPositions {
    bool operator()(const ValuedPosition<T>& a, const ValuedPosition<T>& b) {
        return (a.row < b.row) || ((a.row == b.row) && (a.col < b.col)) || ((a.row == b.row) && (a.col == b.col) && (a.val < b.val));
    }
};

template<typename T>
std::vector<size_type> getSuccessorsMat(std::vector<std::vector<T>>& mat, size_type i, T null = T()) {

    std::vector<size_type> res;
    for (size_type j = 0; j < mat[i].size(); j++) {
        if (mat[i][j] != null) {
            res.push_back(j);
        }
    }

    return res;

}

template<typename T>
std::vector<size_type> getPredecessorsMat(std::vector<std::vector<T>>& mat, size_type j, T null = T()) {

    std::vector<size_type> res;
    for (size_type i = 0; i < mat.size(); i++) {
        if (mat[i][j] != null) {
            res.push_back(i);
        }
    }

    return res;

}

template<typename T>
std::vector<std::pair<size_type, size_type>> getRangeMat(std::vector<std::vector<T>>& mat, size_type i1, size_type i2, size_type j1, size_type j2, T null = T()) {

    std::vector<std::pair<size_type, size_type>> res;

    for (size_type i = i1; i <= i2; i++) {
        for (size_type j = j1; j <= j2; j++) {
            if (mat[i][j] != null) res.push_back(std::make_pair(i, j));
        }
    }

    return res;

}

template<typename T>
bool containsLinkMat(std::vector<std::vector<T>>& mat, size_type i1, size_type i2, size_type j1, size_type j2, T null = T()) {

    for (size_type i = i1; i <= i2; i++) {
        for (size_type j = j1; j <= j2; j++) {
            if (mat[i][j] != null) return true;
        }
    }

    return false;

}

template<typename T>
size_type countLinksMat(std::vector<std::vector<T>>& mat, T null = T()) {

    size_type res = 0;
    for (auto i = 0; i < mat.size(); i++) {
        for (auto j = 0; j < mat[i].size(); j++) {
            res += (mat[i][j] != null);
        }
    }

    return res;

}

template<typename T>
std::vector<T> getSuccessorElementsMat(std::vector<std::vector<T>>& mat, size_type i, T null = T()) {

    std::vector<T> succs;
    for (size_type j = 0; j < mat[i].size(); j++) {
        if (mat[i][j] != null) {
            succs.push_back(mat[i][j]);
        }
    }

    return succs;

}

template<typename T>
std::vector<size_type> getSuccessorPositionsMat(std::vector<std::vector<T>>& mat, size_type i, T null = T()) {

    std::vector<size_type> succs;
    for (size_type j = 0; j < mat[i].size(); j++) {
        if (mat[i][j] != null) {
            succs.push_back(j);
        }
    }

    return succs;

}

template<typename T>
std::vector<ValuedPosition<T>> getSuccessorValuedPositionsMat(std::vector<std::vector<T>>& mat, size_type i, T null = T()) {

    std::vector<ValuedPosition<T>> succs;
    for (size_type j = 0; j < mat[i].size(); j++) {
        if (mat[i][j] != null) {
            succs.push_back(ValuedPosition<T>(i, j, mat[i][j]));
        }
    }

    return succs;

}

template<typename T>
std::vector<T> getPredecessorElementsMat(std::vector<std::vector<T>>& mat, size_type j, T null = T()) {

    std::vector<T> preds;
    for (size_type i = 0; i < mat.size(); i++) {
        if (mat[i][j] != null) {
            preds.push_back(mat[i][j]);
        }
    }

    return preds;

}

template<typename T>
std::vector<size_type> getPredecessorPositionsMat(std::vector<std::vector<T>>& mat, size_type j, T null = T()) {

    std::vector<size_type> preds;
    for (size_type i = 0; i < mat.size(); i++) {
        if (mat[i][j] != null) {
            preds.push_back(i);
        }
    }

    return preds;

}

template<typename T>
std::vector<ValuedPosition<T>> getPredecessorValuedPositionsMat(std::vector<std::vector<T>>& mat, size_type j, T null = T()) {

    std::vector<ValuedPosition<T>> preds;
    for (size_type i = 0; i < mat.size(); i++) {
        if (mat[i][j] != null) {
            preds.push_back(ValuedPosition<T>(i, j, mat[i][j]));
        }
    }

    return preds;

}

template<typename T>
std::vector<T> getElementsInRangeMat(std::vector<std::vector<T>>& mat, size_type i1, size_type i2, size_type j1, size_type j2, T null = T()) {

    std::vector<T> elems;
    for (size_type i = i1; i <= i2; i++) {
        for (size_type j = j1; j <= j2; j++) {
            if (mat[i][j] != null) {
                elems.push_back(mat[i][j]);
            }
        }
    }

    return elems;

}

template<typename T>
std::vector<std::pair<size_type, size_type>> getPositionsInRangeMat(std::vector<std::vector<T>>& mat, size_type i1, size_type i2, size_type j1, size_type j2, T null = T()) {

    std::vector<std::pair<size_type, size_type>> elems;
    for (size_type i = i1; i <= i2; i++) {
        for (size_type j = j1; j <= j2; j++) {
            if (mat[i][j] != null) {
                elems.push_back(std::make_pair(i, j));
            }
        }
    }

    return elems;

}

template<typename T>
std::vector<ValuedPosition<T>> getValuedPositionsInRangeMat(std::vector<std::vector<T>>& mat, size_type i1, size_type i2, size_type j1, size_type j2, T null = T()) {

    std::vector<ValuedPosition<T>> elems;
    for (size_type i = i1; i <= i2; i++) {
        for (size_type j = j1; j <= j2; j++) {
            if (mat[i][j] != null) {
                elems.push_back(ValuedPosition<T>(i, j, mat[i][j]));
            }
        }
    }

    return elems;

}

template<typename T>
std::vector<T> getAllElementsMat(std::vector<std::vector<T>>& mat, T null = T()) {
    return getElementsInRangeMat(mat, 0, mat.size() - 1, 0, mat[0].size() - 1, null);
}

template<typename T>
std::vector<std::pair<size_type, size_type>> getAllPositionsMat(std::vector<std::vector<T>>& mat, T null = T()) {
    return getPositionsInRangeMat(mat, 0, mat.size() - 1, 0, mat[0].size() - 1, null);
}

template<typename T>
std::vector<ValuedPosition<T>> getAllValuedPositionsMat(std::vector<std::vector<T>>& mat, T null = T()) {
    return getValuedPositionsInRangeMat(mat, 0, mat.size() - 1, 0, mat[0].size() - 1, null);
}

template<typename T>
bool containsElementMat(std::vector<std::vector<T>>& mat, size_type i1, size_type i2, size_type j1, size_type j2, T null = T()) {

    for (size_type i = i1; i <= i2; i++) {
        for (size_type j = j1; j <= j2; j++) {
            if (mat[i][j] != null) return true;
        }
    }

    return false;

}

template<typename T>
size_type countElemsMat(std::vector<std::vector<T>>& mat, T null = T()) {

    size_type res = 0;
    for (auto i = 0; i < mat.size(); i++) {
        for (auto j = 0; j < mat[i].size(); j++) {
            res += (mat[i][j] != null);
        }
    }

    return res;

}


/* Representation of nodes in intermediate tree representation of k^2-trees */

template<typename T>
class Node {
public:
    Node(T lab) {
        lab_ = lab;
    }

    ~Node() {

        for (auto c : children_) {
            delete c;
        }

    }

    bool isLeaf() {
        return children_.size() == 0;
    }

    T getLabel() {
        return lab_;
    }

    bool hasChild(size_type i) {
        return (i < children_.size()) && (children_[i] != 0);
    }

    Node* getChild(size_type i) {
        return (i < children_.size()) ? children_[i] : 0;
    }

    Node* addChild(size_type i, T lab) {

        if (children_[i] == 0) {
            children_[i] = new Node(lab);
        } else {
            children_[i]->lab_ = lab;
        }
        return children_[i];

    }

    void turnInternal(size_type arity, bool f) {

        children_ = std::vector<Node*>(arity);

        if (f) {
            for (auto i = 0; i < arity; i++) {
                children_[i] = new Node(T());
            }
        }

    }

    size_type getArity() {
        return children_.size();
    }

private:
    T lab_;
    std::vector<Node*> children_;

};


/* Representation of subproblems / submatrices during k^2-tree construction */

struct Subproblem {

    size_type firstRow; // x1
    size_type lastRow; // x2
    size_type firstCol; // y1
    size_type lastCol; // y2

    size_type left; // a
    size_type right; // b

    Subproblem();

    Subproblem(size_type fr, size_type lr, size_type fc, size_type lc, size_type l, size_type r);

};


/* Dynamic, but naive rank data structure for intermediate steps. */

class NaiveDynamicRank {

public:
    NaiveDynamicRank();

    NaiveDynamicRank(const std::vector<bool>& arr);

    size_type rank(size_type pos);

    size_type rankSafe(size_type pos);

    void increaseFrom(size_type pos, size_type inc = 1);

    void decreaseFrom(size_type pos, size_type dec = 1);

    void insert(size_type pos, size_type num);


private:
    std::vector<size_type> ranks_;

};

#endif //K2TREES_UTILITY_HPP
