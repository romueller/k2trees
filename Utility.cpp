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

#include "Utility.hpp"

size_type logK(const size_type n, const size_type k) {

    size_type pow = 0;
    for (size_type val = 1; val < n; val *= k, pow++) { }
    return pow;

}

bool isAllZero(const std::vector<bool>& v) {
    return std::find(v.begin(), v.end(), true) == v.end();
}

bool isAllZero(const bit_vector_type& v) {
    return std::find(v.begin(), v.end(), true) == v.end();
}

void printRanks(const rank_type& r) {

    for (size_type i = 0; i < r.size(); i++) std::cout << r.rank(i + 1) << " " << std::flush;
    std::cout << std::endl;

}





// ===== From matrix ... =====

RelationLists boolMatrixToLists(RelationMatrix& mat) {

    RelationLists lists;

    for (size_type i = 0; i < mat.size(); i++) {

        lists.push_back(RelationList());

        for (size_type j = 0; j < mat[i].size(); j++) {
            if (mat[i][j]) {
                lists[i].push_back(j);
            }
        }
    }

    return lists;

}

RelationPairs boolMatrixToPairs(RelationMatrix& mat) {

    RelationPairs pairs;

    for (size_type i = 0; i < mat.size(); i++) {

        for (size_type j = 0; j < mat[i].size(); j++) {

            if (mat[i][j]) {
                pairs.push_back(std::make_pair(i, j));
            }

        }

    }

    return pairs;

}

// ===== From lists ... =====

RelationMatrix boolListsToMatrix(RelationLists& lists, size_type numRows, size_type numCols) {

    RelationMatrix mat(numRows, std::vector<bool>(numCols));

    for (size_type i = 0; i < lists.size(); i++) {

        for (size_type j = 0; j < lists[i].size(); j++) {
            mat[i][lists[i][j]] = true;
        }

    }

    return mat;

}

RelationMatrix boolListsToMatrix(RelationLists& lists) {

    size_type maxCol = 0;

    for (auto l : lists) {

        for (size_type j = 0; j < l.size(); j++) {
            maxCol = std::max(maxCol, l[j]);
        }
    }

    return boolListsToMatrix(lists, lists.size(), maxCol + 1);

}

RelationPairs boolListsToPairs(RelationLists& lists) {

    RelationPairs pairs;

    for (size_type i = 0; i < lists.size(); i++) {

        for (size_type j = 0; j < lists[i].size(); j++) {
            pairs.push_back(std::make_pair(i, lists[i][j]));
        }
    }

    return pairs;

}

// ===== From pairs ... =====

RelationMatrix boolPairsToMatrix(RelationPairs& pairs, size_type numRows, size_type numCols) {

    RelationMatrix mat(numRows, std::vector<bool>(numCols));

    for (auto p : pairs) {
        mat[p.first][p.second] = true;
    }

    return mat;

}

RelationMatrix boolPairsToMatrix(RelationPairs& pairs) {

    size_type maxRow = 0;
    size_type maxCol = 0;

    for (auto p : pairs) {

        maxRow = std::max(maxRow, p.first);
        maxCol = std::max(maxCol, p.second);

    }

    return boolPairsToMatrix(pairs, maxRow + 1, maxCol + 1);

}

RelationLists boolPairsToList(RelationPairs& pairs, size_type numRows) {

    RelationLists lists(numRows);

    for (auto p : pairs) {
        lists[p.first].push_back(p.second);
    }

    return lists;

}

RelationLists boolPairsToList(RelationPairs& pairs) {

    size_type maxRow = 0;

    for (auto p : pairs) {
        maxRow = std::max(maxRow, p.first);
    }

    return boolPairsToList(pairs, maxRow + 1);

}




std::vector<size_type> getSuccessorsMat(RelationMatrix& mat, size_type i) {

    std::vector<size_type> res;
    for (size_type j = 0; j < mat[i].size(); j++) {
        if (mat[i][j]) {
            res.push_back(j);
        }
    }

    return res;

}

std::vector<size_type> getPredecessorsMat(RelationMatrix& mat, size_type j) {

    std::vector<size_type> res;
    for (size_type i = 0; i < mat.size(); i++) {
        if (mat[i][j]) {
            res.push_back(i);
        }
    }

    return res;

}

std::vector<std::pair<size_type, size_type>> getRangeMat(std::vector<std::vector<bool>>& mat, size_type i1, size_type i2, size_type j1, size_type j2) {

    std::vector<std::pair<size_type, size_type>> res;

    for (size_type i = i1; i <= i2; i++) {
        for (size_type j = j1; j <= j2; j++) {
            if (mat[i][j]) res.push_back(std::make_pair(i, j));
        }
    }

    return res;

}

bool containsLinkMat(std::vector<std::vector<bool>>& mat, size_type i1, size_type i2, size_type j1, size_type j2) {

    for (size_type i = i1; i <= i2; i++) {
        for (size_type j = j1; j <= j2; j++) {
            if (mat[i][j]) return true;
        }
    }

    return false;

}

size_type countLinksMat(std::vector<std::vector<bool>>& mat) {

    size_type res = 0;
    for (auto i = 0; i < mat.size(); i++) {
        for (auto j = 0; j < mat[i].size(); j++) {
            res += mat[i][j];
        }
    }

    return res;

}





Subproblem::Subproblem() {
    firstRow = lastRow = firstCol = lastCol = left = right = 0;
}

Subproblem::Subproblem(size_type fr, size_type lr, size_type fc, size_type lc, size_type l, size_type r) {

    firstRow = fr;
    lastRow = lr;
    firstCol = fc;
    lastCol = lc;

    left = l;
    right = r;

}





NaiveDynamicRank::NaiveDynamicRank() {
    // nothing to do
}

NaiveDynamicRank::NaiveDynamicRank(const std::vector<bool>& arr) {

    ranks_ = std::vector<size_type>(arr.size() + 1);
    ranks_[0] = 0;
    for (auto i = 0; i < arr.size(); i++) {
        ranks_[i + 1] = ranks_[i] + arr[i];
    }

}

size_type NaiveDynamicRank::rank(size_type pos) {
    return ranks_[pos];
}

size_type NaiveDynamicRank::rankSafe(size_type pos) {
    return ranks_[std::min(pos, ranks_.size() - 1)];
}

void NaiveDynamicRank::increaseFrom(size_type pos, size_type inc) {

    for (auto i = pos; i < ranks_.size(); i++) {
        ranks_[i] += inc;
    }

}

void NaiveDynamicRank::decreaseFrom(size_type pos, size_type dec) {

    for (auto i = pos; i < ranks_.size(); i++) {
        ranks_[i] = (ranks_[i] > dec) * (ranks_[i] - dec);
    }

}

void NaiveDynamicRank::insert(size_type pos, size_type num) {
    ranks_.insert(ranks_.begin() + pos, num, (pos > 0) * ranks_[pos - 1]);
}
