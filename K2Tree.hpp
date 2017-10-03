#ifndef K2TREES_K2TREE_HPP
#define K2TREES_K2TREE_HPP

#include "Utility.hpp"

template<typename E>
class K2Tree {

public:
    typedef E elem_type;

    typedef std::vector<std::vector<elem_type>> matrix_type;
    typedef std::vector<std::pair<size_type, elem_type>> list_type;
    typedef std::vector<std::pair<size_type, size_type>> positions_type;
    typedef std::vector<ValuedPosition<elem_type>> pairs_type;

    virtual ~K2Tree() { }

    virtual size_type getNumRows() = 0;

    virtual size_type getNumCols() = 0;

    virtual elem_type getNull() = 0;


    virtual bool isNotNull(size_type i, size_type j) = 0;

    virtual elem_type getElement(size_type i, size_type j) = 0;

    virtual std::vector<elem_type> getSuccessorElements(size_type i) = 0;

    virtual std::vector<size_type> getSuccessorPositions(size_type i) = 0;

    virtual pairs_type getSuccessorValuedPositions(size_type i) = 0;

    virtual std::vector<elem_type> getPredecessorElements(size_type j) = 0;

    virtual std::vector<size_type> getPredecessorPositions(size_type j) = 0;

    virtual pairs_type getPredecessorValuedPositions(size_type j) = 0;

    virtual std::vector<elem_type> getElementsInRange(size_type i1, size_type i2, size_type j1, size_type j2) = 0;

    virtual positions_type getPositionsInRange(size_type i1, size_type i2, size_type j1, size_type j2) = 0;

    virtual pairs_type getValuedPositionsInRange(size_type i1, size_type i2, size_type j1, size_type j2) = 0;

    virtual std::vector<elem_type> getAllElements() = 0;

    virtual positions_type getAllPositions() = 0;

    virtual pairs_type getAllValuedPositions() = 0;

    virtual bool containsElement(size_type i1, size_type i2, size_type j1, size_type j2) = 0;

    virtual size_type countElements() = 0;


    virtual K2Tree* clone() const = 0;

    virtual void print(bool all = false) = 0;

    virtual bool compare(matrix_type& mat, elem_type null, bool silent) {

        bool overallEqual = true;
        bool equal;
        unsigned long cnt;

#if 1
        equal = true;
        for (size_type i = 0; i < mat.size(); i++) {
            for (size_type j = 0; j < mat[i].size(); j++) {
                equal = equal && ((mat[i][j] != null) == isNotNull(i, j));
            }
        }
        if (!silent) std::cout << "isNotNull: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_type i = 0; i < mat.size(); i++) {
            for (size_type j = 0; j < mat[i].size(); j++) {
                equal = equal && (mat[i][j] == getElement(i, j));
            }
        }
        if (!silent) std::cout << "getElement: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_type i = 0; i < mat.size(); i++) {

            auto matSuccs = getSuccessorElementsMat(mat, i, null);
            auto treeSuccs = getSuccessorElements(i);

            equal = equal && (matSuccs.size() == treeSuccs.size());
            if (equal) {
                for (size_type k = 0; k < matSuccs.size(); k++) {
                    equal = equal && (matSuccs[k] == treeSuccs[k]);
                }
            }

        }
        if (!silent) std::cout << "getSuccessorElements: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_type i = 0; i < mat.size(); i++) {

            auto matSuccs = getSuccessorPositionsMat(mat, i, null);
            auto treeSuccs = getSuccessorPositions(i);

            equal = equal && (matSuccs.size() == treeSuccs.size());
            if (equal) {
                for (size_type k = 0; k < matSuccs.size(); k++) {
                    equal = equal && (matSuccs[k] == treeSuccs[k]);
                }
            }

        }
        if (!silent) std::cout << "getSuccessorPositions: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_type i = 0; i < mat.size(); i++) {

            auto matSuccs = getSuccessorValuedPositionsMat(mat, i, null);
            auto treeSuccs = getSuccessorValuedPositions(i);

            equal = equal && (matSuccs.size() == treeSuccs.size());
            if (equal) {
                for (size_type k = 0; k < matSuccs.size(); k++) {
                    equal = equal && (matSuccs[k] == treeSuccs[k]);
                }
            }

        }
        if (!silent) std::cout << "getSuccessorValuedPositions: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_t j = 0; j < mat[0].size(); j++) {

            auto matPreds = getPredecessorElementsMat(mat, j, null);
            auto treePreds = getPredecessorElements(j);

            equal = equal && (matPreds.size() == treePreds.size());
            if (equal) {
                for (size_type k = 0; k < matPreds.size(); k++) {
                    equal = equal && (matPreds[k] == treePreds[k]);
                }
            }

        }
        if (!silent) std::cout << "getPredecessorElements: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_t j = 0; j < mat[0].size(); j++) {

            auto matPreds = getPredecessorPositionsMat(mat, j, null);
            auto treePreds = getPredecessorPositions(j);

            equal = equal && (matPreds.size() == treePreds.size());
            if (equal) {
                for (size_type k = 0; k < matPreds.size(); k++) {
                    equal = equal && (matPreds[k] == treePreds[k]);
                }
            }

        }
        if (!silent) std::cout << "getPredecessorPositions: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_t j = 0; j < mat[0].size(); j++) {

            auto matPreds = getPredecessorValuedPositionsMat(mat, j, null);
            auto treePreds = getPredecessorValuedPositions(j);

            equal = equal && (matPreds.size() == treePreds.size());
            if (equal) {
                for (size_type k = 0; k < matPreds.size(); k++) {
                    equal = equal && (matPreds[k] == treePreds[k]);
                }
            }

        }
        if (!silent) std::cout << "getPredecessorValuedPositions: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        cnt = 0;
        if (!silent) std::cout << "getElementsInRange: " << std::flush;
        for (size_t i1 = 0; i1 < mat.size(); i1++) {
            for (size_t i2 = i1; i2 < mat.size(); i2++) {
                for (size_t j1 = 0; j1 < mat[i1].size(); j1++) {
                    for (size_t j2 = j1; j2 < mat[i1].size(); j2++) {

                        cnt++;
                        if (!silent && (cnt > 0) && (cnt % 10000 == 0)) {
                            std::cout << "." << std::flush;
                        }

                        auto matPairs = getElementsInRangeMat(mat, i1, i2, j1, j2, null); std::sort(matPairs.begin(), matPairs.end());
                        auto treePairs = getElementsInRange(i1, i2, j1, j2); std::sort(treePairs.begin(), treePairs.end());

                        equal = equal && (matPairs.size() == treePairs.size());
                        if (equal) {
                            for (size_type k = 0; k < matPairs.size(); k++) {
                                equal = equal && (matPairs[k] == treePairs[k]);
                            }
                        }

                    }
                }
            }
        }
        if (!silent) std::cout << "\rgetElementsInRange: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif


#if 1
        equal = true;
        cnt = 0;
        if (!silent) std::cout << "getPositionsInRange: " << std::flush;
        for (size_t i1 = 0; i1 < mat.size(); i1++) {
            for (size_t i2 = i1; i2 < mat.size(); i2++) {
                for (size_t j1 = 0; j1 < mat[i1].size(); j1++) {
                    for (size_t j2 = j1; j2 < mat[i1].size(); j2++) {

                        cnt++;
                        if (!silent && (cnt > 0) && (cnt % 10000 == 0)) {
                            std::cout << "." << std::flush;
                        }

                        auto matPairs = getPositionsInRangeMat(mat, i1, i2, j1, j2, null); std::sort(matPairs.begin(), matPairs.end(), sortPairs<size_type, size_type>());
                        auto treePairs = getPositionsInRange(i1, i2, j1, j2); std::sort(treePairs.begin(), treePairs.end(), sortPairs<size_type, size_type>());

                        equal = equal && (matPairs.size() == treePairs.size());
                        if (equal) {
                            for (size_type k = 0; k < matPairs.size(); k++) {
                                equal = equal && (matPairs[k] == treePairs[k]);
                            }
                        }

                    }
                }
            }
        }
        if (!silent) std::cout << "\rgetPositionsInRange: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif


#if 1
        equal = true;
        cnt = 0;
        if (!silent) std::cout << "getValuedPositionsInRange: " << std::flush;
        for (size_t i1 = 0; i1 < mat.size(); i1++) {
            for (size_t i2 = i1; i2 < mat.size(); i2++) {
                for (size_t j1 = 0; j1 < mat[i1].size(); j1++) {
                    for (size_t j2 = j1; j2 < mat[i1].size(); j2++) {

                        cnt++;
                        if (!silent && (cnt > 0) && (cnt % 10000 == 0)) {
                            std::cout << "." << std::flush;
                        }

                        auto matPairs = getValuedPositionsInRangeMat(mat, i1, i2, j1, j2, null); std::sort(matPairs.begin(), matPairs.end(), sortValuedPositions<elem_type>());
                        auto treePairs = getValuedPositionsInRange(i1, i2, j1, j2); std::sort(treePairs.begin(), treePairs.end(), sortValuedPositions<elem_type>());

                        equal = equal && (matPairs.size() == treePairs.size());
                        if (equal) {
                            for (size_type k = 0; k < matPairs.size(); k++) {
                                equal = equal && (matPairs[k] == treePairs[k]);
                            }
                        }

                    }
                }
            }
        }
        if (!silent) std::cout << "\rgetValuedPositionsInRange: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        {
            auto matElems = getAllElementsMat(mat, null); std::sort(matElems.begin(), matElems.end());
            auto treeElems = getAllElements(); std::sort(treeElems.begin(), treeElems.end());

            equal = equal && (matElems.size() == treeElems.size());
            if (equal) {
                for (size_type k = 0; k < matElems.size(); k++) {
                    equal = equal && (matElems[k] == treeElems[k]);
                }
            }
        }
        if (!silent) std::cout << "getAllElements: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        {
            auto matElems = getAllPositionsMat(mat, null); std::sort(matElems.begin(), matElems.end(), sortPairs<size_type, size_type>());
            auto treeElems = getAllPositions(); std::sort(treeElems.begin(), treeElems.end(), sortPairs<size_type, size_type>());

            equal = equal && (matElems.size() == treeElems.size());
            if (equal) {
                for (size_type k = 0; k < matElems.size(); k++) {
                    equal = equal && (matElems[k] == treeElems[k]);
                }
            }
        }
        if (!silent) std::cout << "getAllPositions: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        {
            auto matElems = getAllValuedPositionsMat(mat, null); std::sort(matElems.begin(), matElems.end(), sortValuedPositions<elem_type>());
            auto treeElems = getAllValuedPositions(); std::sort(treeElems.begin(), treeElems.end(), sortValuedPositions<elem_type>());

            equal = equal && (matElems.size() == treeElems.size());
            if (equal) {
                for (size_type k = 0; k < matElems.size(); k++) {
                    equal = equal && (matElems[k] == treeElems[k]);
                }
            }
        }
        if (!silent) std::cout << "getAllValuedPositions: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_t i1 = 0; i1 < mat.size(); i1++) {
            for (size_t i2 = i1; i2 < mat.size(); i2++) {
                for (size_t j1 = 0; j1 < mat[i1].size(); j1++) {
                    for (size_t j2 = j1; j2 < mat[i1].size(); j2++) {

                        bool matFlag = containsElementMat(mat, i1, i2, j1, j2, null);
                        bool treeFlag = containsElement(i1, i2, j1, j2);

                        equal = equal && (matFlag == treeFlag);

                    }
                }
            }
        }
        if (!silent) std::cout << "containsElement: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        {
            size_type matCnt = countElemsMat(mat, null);
            size_type treeCnt = countElements();

            equal = (matCnt == treeCnt);
        }
        if (!silent) std::cout << "countElements: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_type i = 0; i < mat.size(); i++) {
            for (size_type j = 0; j < mat[i].size(); j++) {
                equal = equal && (mat[i][j] == areRelated(i, j));
            }
        }
        if (!silent) std::cout << "areRelated: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_t i = 0; i < mat.size(); i++) {

            auto matSuccs = getSuccessorsMat(mat, i, null);
            auto treeSuccs = getSuccessors(i);

            equal = equal && (matSuccs.size() == treeSuccs.size());
            if (equal) {
                for (size_type k = 0; k < matSuccs.size(); k++) {
                    equal = equal && (matSuccs[k] == treeSuccs[k]);
                }
            }

        }
        if (!silent) std::cout << "getSuccessors: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_t j = 0; j < mat[0].size(); j++) {

            auto matPreds = getPredecessorsMat(mat, j, null);
            auto treePreds = getPredecessors(j);

            equal = equal && (matPreds.size() == treePreds.size());
            if (equal) {
                for (size_type k = 0; k < matPreds.size(); k++) {
                    equal = equal && (matPreds[k] == treePreds[k]);
                }
            }

        }
        if (!silent) std::cout << "getPredecessors: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        cnt = 0;
        if (!silent) std::cout << "getRange: " << std::flush;
        for (size_t i1 = 0; i1 < mat.size(); i1++) {
            for (size_t i2 = i1; i2 < mat.size(); i2++) {
                for (size_t j1 = 0; j1 < mat[i1].size(); j1++) {
                    for (size_t j2 = j1; j2 < mat[i1].size(); j2++) {

                        cnt++;
                        if (!silent && (cnt > 0) && (cnt % 10000 == 0)) {
                            std::cout << "." << std::flush;
                        }

                        auto matPairs = getRangeMat(mat, i1, i2, j1, j2, null); std::sort(matPairs.begin(), matPairs.end(), sortPairs<size_type, size_type>());
                        auto treePairs = getRange(i1, i2, j1, j2); std::sort(treePairs.begin(), treePairs.end(), sortPairs<size_type, size_type>());

                        equal = equal && (matPairs.size() == treePairs.size());
                        if (equal) {
                            for (size_type k = 0; k < matPairs.size(); k++) {
                                equal = equal && (matPairs[k] == treePairs[k]);
                            }
                        }

                    }
                }
            }
        }
        if (!silent) std::cout << "\rgetRange: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_t i1 = 0; i1 < mat.size(); i1++) {
            for (size_t i2 = i1; i2 < mat.size(); i2++) {
                for (size_t j1 = 0; j1 < mat[i1].size(); j1++) {
                    for (size_t j2 = j1; j2 < mat[i1].size(); j2++) {

                        bool matFlag = containsLinkMat(mat, i1, i2, j1, j2, null);
                        bool treeFlag = containsLink(i1, i2, j1, j2);

                        equal = equal && (matFlag == treeFlag);

                    }
                }
            }
        }
        if (!silent) std::cout << "containsLink: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        {
            size_type matCnt = countLinksMat(mat, null);
            size_type treeCnt = countLinks();

            equal = (matCnt == treeCnt);
        }
        if (!silent) std::cout << "countLinks: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

        return overallEqual;

    }

    virtual bool compare(K2Tree& other, bool silent) {

        bool overallEqual = true;
        bool equal;
        unsigned long cnt;

#if 1
        equal = true;
        for (size_type i = 0; i < getNumRows(); i++) {
            for (size_type j = 0; j < getNumCols(); j++) {
                equal = equal && (isNotNull(i, j) == other.isNotNull(i, j));
            }
        }
        if (!silent) std::cout << "isNotNull: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_type i = 0; i < getNumRows(); i++) {
            for (size_type j = 0; j < getNumCols(); j++) {
                equal = equal && (getElement(i, j) == other.getElement(i, j));
            }
        }
        if (!silent) std::cout << "getElement: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_type i = 0; i < getNumRows(); i++) {

            auto otherSuccs = other.getSuccessorElements(i);
            auto treeSuccs = getSuccessorElements(i);

            equal = equal && (otherSuccs.size() == treeSuccs.size());
            if (equal) {
                for (size_type k = 0; k < treeSuccs.size(); k++) {
                    equal = equal && (otherSuccs[k] == treeSuccs[k]);
                }
            }

        }
        if (!silent) std::cout << "getSuccessorElements: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_type i = 0; i < getNumRows(); i++) {

            auto otherSuccs = other.getSuccessorPositions(i);
            auto treeSuccs = getSuccessorPositions(i);

            equal = equal && (otherSuccs.size() == treeSuccs.size());
            if (equal) {
                for (size_type k = 0; k < treeSuccs.size(); k++) {
                    equal = equal && (otherSuccs[k] == treeSuccs[k]);
                }
            }

        }
        if (!silent) std::cout << "getSuccessorPositions: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_type i = 0; i < getNumRows(); i++) {

            auto otherSuccs = other.getSuccessorValuedPositions(i);
            auto treeSuccs = getSuccessorValuedPositions(i);

            equal = equal && (otherSuccs.size() == treeSuccs.size());
            if (equal) {
                for (size_type k = 0; k < treeSuccs.size(); k++) {
                    equal = equal && (otherSuccs[k] == treeSuccs[k]);
                }
            }

        }
        if (!silent) std::cout << "getSuccessorValuedPositions: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_t j = 0; j < getNumCols(); j++) {

            auto otherPreds = other.getPredecessorElements(j);
            auto treePreds = getPredecessorElements(j);

            equal = equal && (otherPreds.size() == treePreds.size());
            if (equal) {
                for (size_type k = 0; k < treePreds.size(); k++) {
                    equal = equal && (otherPreds[k] == treePreds[k]);
                }
            }

        }
        if (!silent) std::cout << "getPredecessorElements: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_t j = 0; j < getNumCols(); j++) {

            auto otherPreds = other.getPredecessorPositions(j);
            auto treePreds = getPredecessorPositions(j);

            equal = equal && (otherPreds.size() == treePreds.size());
            if (equal) {
                for (size_type k = 0; k < treePreds.size(); k++) {
                    equal = equal && (otherPreds[k] == treePreds[k]);
                }
            }

        }
        if (!silent) std::cout << "getPredecessorPositions: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_t j = 0; j < getNumCols(); j++) {

            auto otherPreds = other.getPredecessorValuedPositions(j);
            auto treePreds = getPredecessorValuedPositions(j);

            equal = equal && (otherPreds.size() == treePreds.size());
            if (equal) {
                for (size_type k = 0; k < treePreds.size(); k++) {
                    equal = equal && (otherPreds[k] == treePreds[k]);
                }
            }

        }
        if (!silent) std::cout << "getPredecessorValuedPositions: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        cnt = 0;
        if (!silent) std::cout << "getElementsInRange: " << std::flush;
        for (size_t i1 = 0; i1 < getNumRows(); i1++) {
            for (size_t i2 = i1; i2 < getNumRows(); i2++) {
                for (size_t j1 = 0; j1 < getNumCols(); j1++) {
                    for (size_t j2 = j1; j2 < getNumCols(); j2++) {

                        cnt++;
                        if (!silent && (cnt > 0) && (cnt % 10000 == 0)) {
                            std::cout << "." << std::flush;
                        }

                        auto otherPairs = other.getElementsInRange(i1, i2, j1, j2); std::sort(otherPairs.begin(), otherPairs.end());
                        auto treePairs = getElementsInRange(i1, i2, j1, j2); std::sort(treePairs.begin(), treePairs.end());

                        equal = equal && (otherPairs.size() == treePairs.size());
                        if (equal) {
                            for (size_type k = 0; k < treePairs.size(); k++) {
                                equal = equal && (otherPairs[k] == treePairs[k]);
                            }
                        }

                    }
                }
            }
        }
        if (!silent) std::cout << "\rgetElementsInRange: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif


#if 1
        equal = true;
        cnt = 0;
        if (!silent) std::cout << "getPositionsInRange: " << std::flush;
        for (size_t i1 = 0; i1 < getNumRows(); i1++) {
            for (size_t i2 = i1; i2 < getNumRows(); i2++) {
                for (size_t j1 = 0; j1 < getNumCols(); j1++) {
                    for (size_t j2 = j1; j2 < getNumCols(); j2++) {

                        cnt++;
                        if (!silent && (cnt > 0) && (cnt % 10000 == 0)) {
                            std::cout << "." << std::flush;
                        }

                        auto otherPairs = other.getPositionsInRange(i1, i2, j1, j2); std::sort(otherPairs.begin(), otherPairs.end(), sortPairs<size_type, size_type>());
                        auto treePairs = getPositionsInRange(i1, i2, j1, j2); std::sort(treePairs.begin(), treePairs.end(), sortPairs<size_type, size_type>());

                        equal = equal && (otherPairs.size() == treePairs.size());
                        if (equal) {
                            for (size_type k = 0; k < treePairs.size(); k++) {
                                equal = equal && (otherPairs[k] == treePairs[k]);
                            }
                        }

                    }
                }
            }
        }
        if (!silent) std::cout << "\rgetPositionsInRange: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif


#if 1
        equal = true;
        cnt = 0;
        if (!silent) std::cout << "getValuedPositionsInRange: " << std::flush;
        for (size_t i1 = 0; i1 < getNumRows(); i1++) {
            for (size_t i2 = i1; i2 < getNumRows(); i2++) {
                for (size_t j1 = 0; j1 < getNumCols(); j1++) {
                    for (size_t j2 = j1; j2 < getNumCols(); j2++) {

                        cnt++;
                        if (!silent && (cnt > 0) && (cnt % 10000 == 0)) {
                            std::cout << "." << std::flush;
                        }

                        auto otherPairs = other.getValuedPositionsInRange(i1, i2, j1, j2); std::sort(otherPairs.begin(), otherPairs.end(), sortValuedPositions<elem_type>());
                        auto treePairs = getValuedPositionsInRange(i1, i2, j1, j2); std::sort(treePairs.begin(), treePairs.end(), sortValuedPositions<elem_type>());

                        equal = equal && (otherPairs.size() == treePairs.size());
                        if (equal) {
                            for (size_type k = 0; k < treePairs.size(); k++) {
                                equal = equal && (otherPairs[k] == treePairs[k]);
                            }
                        }

                    }
                }
            }
        }
        if (!silent) std::cout << "\rgetValuedPositionsInRange: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        {
            auto otherElems = other.getAllElements(); std::sort(otherElems.begin(), otherElems.end());
            auto treeElems = getAllElements(); std::sort(treeElems.begin(), treeElems.end());

            equal = equal && (otherElems.size() == treeElems.size());
            if (equal) {
                for (size_type k = 0; k < treeElems.size(); k++) {
                    equal = equal && (otherElems[k] == treeElems[k]);
                }
            }
        }
        if (!silent) std::cout << "getAllElements: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        {
            auto otherElems = other.getAllPositions(); std::sort(otherElems.begin(), otherElems.end(), sortPairs<size_type, size_type>());
            auto treeElems = getAllPositions(); std::sort(treeElems.begin(), treeElems.end(), sortPairs<size_type, size_type>());

            equal = equal && (otherElems.size() == treeElems.size());
            if (equal) {
                for (size_type k = 0; k < treeElems.size(); k++) {
                    equal = equal && (otherElems[k] == treeElems[k]);
                }
            }
        }
        if (!silent) std::cout << "getAllPositions: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        {
            auto otherElems = other.getAllValuedPositions(); std::sort(otherElems.begin(), otherElems.end(), sortValuedPositions<elem_type>());
            auto treeElems = getAllValuedPositions(); std::sort(treeElems.begin(), treeElems.end(), sortValuedPositions<elem_type>());

            equal = equal && (otherElems.size() == treeElems.size());
            if (equal) {
                for (size_type k = 0; k < treeElems.size(); k++) {
                    equal = equal && (otherElems[k] == treeElems[k]);
                }
            }
        }
        if (!silent) std::cout << "getAllValuedPositions: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_t i1 = 0; i1 < getNumRows(); i1++) {
            for (size_t i2 = i1; i2 < getNumRows(); i2++) {
                for (size_t j1 = 0; j1 < getNumCols(); j1++) {
                    for (size_t j2 = j1; j2 < getNumCols(); j2++) {

                        bool otherFlag = other.containsElement(i1, i2, j1, j2);
                        bool treeFlag = containsElement(i1, i2, j1, j2);

                        equal = equal && (otherFlag == treeFlag);

                    }
                }
            }
        }
        if (!silent) std::cout << "containsElement: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        {
            size_type otherCnt = other.countElements();
            size_type treeCnt = countElements();

            equal = (otherCnt == treeCnt);
        }
        if (!silent) std::cout << "countElements: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_type i = 0; i < getNumRows(); i++) {
            for (size_type j = 0; j < getNumCols(); j++) {
                equal = equal && (areRelated(i, j) == other.areRelated(i, j));
            }
        }
        if (!silent) std::cout << "areRelated: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_t i = 0; i < getNumRows(); i++) {

            auto otherSuccs = other.getSuccessors(i);
            auto treeSuccs = getSuccessors(i);

            equal = equal && (otherSuccs.size() == treeSuccs.size());
            if (equal) {
                for (size_type k = 0; k < treeSuccs.size(); k++) {
                    equal = equal && (otherSuccs[k] == treeSuccs[k]);
                }
            }

        }
        if (!silent) std::cout << "getSuccessors: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_t j = 0; j < getNumCols(); j++) {

            auto otherPreds = other.getPredecessors(j);
            auto treePreds = getPredecessors(j);

            equal = equal && (otherPreds.size() == treePreds.size());
            if (equal) {
                for (size_type k = 0; k < treePreds.size(); k++) {
                    equal = equal && (otherPreds[k] == treePreds[k]);
                }
            }

        }
        if (!silent) std::cout << "getPredecessors: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        cnt = 0;
        if (!silent) std::cout << "getRange: " << std::flush;
        for (size_t i1 = 0; i1 < getNumRows(); i1++) {
            for (size_t i2 = i1; i2 < getNumRows(); i2++) {
                for (size_t j1 = 0; j1 < getNumCols(); j1++) {
                    for (size_t j2 = j1; j2 < getNumCols(); j2++) {

                        cnt++;
                        if (!silent && (cnt > 0) && (cnt % 10000 == 0)) {
                            std::cout << "." << std::flush;
                        }

                        auto otherPairs = other.getRange(i1, i2, j1, j2); std::sort(otherPairs.begin(), otherPairs.end(), sortPairs<size_type, size_type>());
                        auto treePairs = getRange(i1, i2, j1, j2); std::sort(treePairs.begin(), treePairs.end(), sortPairs<size_type, size_type>());

                        equal = equal && (otherPairs.size() == treePairs.size());
                        if (equal) {
                            for (size_type k = 0; k < treePairs.size(); k++) {
                                equal = equal && (otherPairs[k] == treePairs[k]);
                            }
                        }

                    }
                }
            }
        }
        if (!silent) std::cout << "\rgetRange: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_t i1 = 0; i1 < getNumRows(); i1++) {
            for (size_t i2 = i1; i2 < getNumRows(); i2++) {
                for (size_t j1 = 0; j1 < getNumCols(); j1++) {
                    for (size_t j2 = j1; j2 < getNumCols(); j2++) {

                        bool otherFlag = other.containsLink(i1, i2, j1, j2);
                        bool treeFlag = containsLink(i1, i2, j1, j2);

                        equal = equal && (otherFlag == treeFlag);

                    }
                }
            }
        }
        if (!silent) std::cout << "containsLink: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        {
            size_type otherCnt = other.countLinks();
            size_type treeCnt = countLinks();

            equal = (otherCnt == treeCnt);
        }
        if (!silent) std::cout << "countLinks: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

        return overallEqual;

    }

    virtual void setNull(size_type i, size_type j) = 0;


    // method aliases using "relation nomenclature"

    virtual bool areRelated(size_type i, size_type j) = 0;

    virtual std::vector<size_type> getSuccessors(size_type i) = 0;

    virtual size_type getFirstSuccessor(size_type i) = 0;

    virtual std::vector<size_type> getPredecessors(size_type j) = 0;

    virtual positions_type getRange(size_type i1, size_type i2, size_type j1, size_type j2) = 0;

    virtual bool containsLink(size_type i1, size_type i2, size_type j1, size_type j2) = 0;

    virtual size_type countLinks() = 0;

};

#endif //K2TREES_K2TREE_HPP
