#ifndef K2TREES_ROWTREE_HPP
#define K2TREES_ROWTREE_HPP

#include "Utility.hpp"

template<typename E>
class RowTree {

public:
    typedef E elem_type;

    typedef std::vector<std::pair<size_type, elem_type>> list_type; // position ("column number") + value

    virtual ~RowTree() { }

    virtual size_type getLength() = 0;

    virtual elem_type getNull() = 0;


    virtual bool isNotNull(size_type i) = 0;

    virtual elem_type getElement(size_type i) = 0;

    virtual std::vector<elem_type> getElementsInRange(size_type l, size_type r) = 0;

    virtual std::vector<size_type> getPositionsInRange(size_type l, size_type r) = 0;

    virtual list_type getValuedPositionsInRange(size_type l, size_type r) = 0;

    virtual std::vector<elem_type> getAllElements() = 0;

    virtual std::vector<size_type> getAllPositions() = 0;

    virtual list_type getAllValuedPositions() = 0;

    virtual bool containsElement(size_type l, size_type r) = 0;

    virtual size_type countElements() = 0;


    virtual RowTree* clone() const = 0;

    virtual void print(bool all = false) = 0;

    virtual bool compare(std::vector<elem_type>& v, elem_type null, bool silent) {

        bool overallEqual = true;
        bool equal;

#if 1
        equal = true;
        for (size_type i = 0; i < v.size(); i++) {
            equal = equal && ((v[i] != null) == isNotNull(i));
        }
        if (!silent) std::cout << "isNotNull: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_type i = 0; i < v.size(); i++) {
            equal = equal && (v[i] == getElement(i));
        }
        if (!silent) std::cout << "getElement: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_t l = 0; l < v.size(); l++) {
            for (size_t r = l; r < v.size(); r++) {

                auto elems = getElementsInRange(l, r);
                std::sort(elems.begin(), elems.end());
                auto vecElems = std::vector<elem_type>();
                for (auto i = l; i <= r; i++) {
                    if (v[i] != null) vecElems.push_back(v[i]);
                }
                std::sort(vecElems.begin(), vecElems.end());

                equal = equal && (elems.size() == vecElems.size());
                if (equal) {
                    for (size_type i = 0; i < elems.size(); i++) {
                        equal = equal && (elems[i] == vecElems[i]);
                    }
                }

            }
        }
        if (!silent) std::cout << "getElementsInRange: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_t l = 0; l < v.size(); l++) {
            for (size_t r = l; r < v.size(); r++) {

                auto positions = getPositionsInRange(l, r);
                std::sort(positions.begin(), positions.end());
                auto vecPositions = std::vector<size_type>();
                for (auto i = l; i <= r; i++) {
                    if (v[i] != null) vecPositions.push_back(i);
                }
                std::sort(vecPositions.begin(), vecPositions.end());

                equal = equal && (positions.size() == vecPositions.size());
                if (equal) {
                    for (size_type i = 0; i < positions.size(); i++) {
                        equal = equal && (positions[i] == vecPositions[i]);
                    }
                }

            }
        }
        if (!silent) std::cout << "getPositionsInRange: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_t l = 0; l < v.size(); l++) {
            for (size_t r = l; r < v.size(); r++) {

                auto positions = getValuedPositionsInRange(l, r);
                std::sort(positions.begin(), positions.end(), sortPairs<size_type, elem_type>());
                auto vecPositions = list_type();
                for (auto i = l; i <= r; i++) {
                    if (v[i] != null) vecPositions.push_back(std::make_pair(i, v[i]));
                }
                std::sort(vecPositions.begin(), vecPositions.end(), sortPairs<size_type, elem_type>());

                equal = equal && (positions.size() == vecPositions.size());
                if (equal) {
                    for (size_type i = 0; i < positions.size(); i++) {
                        equal = equal && (positions[i] == vecPositions[i]);
                    }
                }

            }
        }
        if (!silent) std::cout << "getValuedPositionsInRange: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        {
            auto elems = getAllElements();
            std::sort(elems.begin(), elems.end());
            auto vecElems = std::vector<elem_type>();
            for (size_type i = 0; i < v.size(); i++) {
                if (v[i] != null) vecElems.push_back(v[i]);
            }
            std::sort(vecElems.begin(), vecElems.end());

            equal = equal && (elems.size() == vecElems.size());
            if (equal) {
                for (size_type i = 0; i < elems.size(); i++) {
                    equal = equal && (elems[i] == vecElems[i]);
                }
            }

        }
        if (!silent) std::cout << "getAllElements: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        {
            auto positions = getAllPositions();
            std::sort(positions.begin(), positions.end());
            auto vecPositions = std::vector<size_type>();
            for (size_type i = 0; i < v.size(); i++) {
                if (v[i] != null) vecPositions.push_back(i);
            }
            std::sort(vecPositions.begin(), vecPositions.end());

            equal = equal && (positions.size() == vecPositions.size());
            if (equal) {
                for (size_type i = 0; i < positions.size(); i++) {
                    equal = equal && (positions[i] == vecPositions[i]);
                }
            }

        }
        if (!silent) std::cout << "getAllPositions: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        {
            auto positions = getAllValuedPositions();
            std::sort(positions.begin(), positions.end(), sortPairs<size_type, elem_type>());
            auto vecPositions = list_type();
            for (size_type i = 0; i < v.size(); i++) {
                if (v[i] != null) vecPositions.push_back(std::make_pair(i, v[i]));
            }
            std::sort(vecPositions.begin(), vecPositions.end(), sortPairs<size_type, elem_type>());

            equal = equal && (positions.size() == vecPositions.size());
            if (equal) {
                for (size_type i = 0; i < positions.size(); i++) {
                    equal = equal && (positions[i] == vecPositions[i]);
                }
            }

        }
        if (!silent) std::cout << "getAllValuedPositions: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_t l = 0; l < v.size(); l++) {
            for (size_t r = l; r < v.size(); r++) {

                bool flag = containsElement(l, r);
                bool vecFlag = false;
                for (auto i = l; i <= r && !vecFlag; i++) {
                    if (v[i] != null) vecFlag = true;
                }

                equal = equal && (flag == vecFlag);

            }
        }
        if (!silent) std::cout << "containsElement: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        {
            size_type cnt = countElements();
            size_type vecCnt = 0;
            for (auto i = 0; i < v.size(); i++) {
                vecCnt += (v[i] != null);
            }
            equal = (cnt == vecCnt);
        }
        if (!silent) std::cout << "countElements: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

        return overallEqual;

    }

    virtual bool compare(RowTree<elem_type>& other, bool silent) {

        bool overallEqual = true;
        bool equal;

#if 1
        equal = true;
        for (size_type i = 0; i < getLength(); i++) {
            equal = equal && (isNotNull(i) == other.isNotNull(i));
        }
        if (!silent) std::cout << "isNotNull: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_type i = 0; i < getLength(); i++) {
            equal = equal && (getElement(i) == other.getElement(i));
        }
        if (!silent) std::cout << "getElement: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_t l = 0; l < getLength(); l++) {
            for (size_t r = l; r < getLength(); r++) {

                auto elems = getElementsInRange(l, r); std::sort(elems.begin(), elems.end());
                auto otherElems = other.getElementsInRange(l, r); std::sort(otherElems.begin(), otherElems.end());

                equal = equal && (elems.size() == otherElems.size());
                if (equal) {
                    for (size_type i = 0; i < elems.size(); i++) {
                        equal = equal && (elems[i] == otherElems[i]);
                    }
                }

            }
        }
        if (!silent) std::cout << "getElementsInRange: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_t l = 0; l < getLength(); l++) {
            for (size_t r = l; r < getLength(); r++) {

                auto positions = getPositionsInRange(l, r); std::sort(positions.begin(), positions.end());
                auto otherPositions = other.getPositionsInRange(l, r); std::sort(otherPositions.begin(), otherPositions.end());

                equal = equal && (positions.size() == otherPositions.size());
                if (equal) {
                    for (size_type i = 0; i < positions.size(); i++) {
                        equal = equal && (positions[i] == otherPositions[i]);
                    }
                }

            }
        }
        if (!silent) std::cout << "getPositionsInRange: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_t l = 0; l < getLength(); l++) {
            for (size_t r = l; r < getLength(); r++) {

                auto positions = getValuedPositionsInRange(l, r); std::sort(positions.begin(), positions.end(), sortPairs<size_type, elem_type>());
                auto otherPositions = other.getValuedPositionsInRange(l, r); std::sort(otherPositions.begin(), otherPositions.end(), sortPairs<size_type, elem_type>());

                equal = equal && (positions.size() == otherPositions.size());
                if (equal) {
                    for (size_type i = 0; i < positions.size(); i++) {
                        equal = equal && (positions[i] == otherPositions[i]);
                    }
                }

            }
        }
        if (!silent) std::cout << "getValuedPositionsInRange: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        {
            auto elems = getAllElements(); std::sort(elems.begin(), elems.end());
            auto otherElems = other.getAllElements(); std::sort(otherElems.begin(), otherElems.end());

            equal = equal && (elems.size() == otherElems.size());
            if (equal) {
                for (size_type i = 0; i < elems.size(); i++) {
                    equal = equal && (elems[i] == otherElems[i]);
                }
            }

        }
        if (!silent) std::cout << "getAllElements: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        {
            auto positions = getAllPositions(); std::sort(positions.begin(), positions.end());
            auto otherPositions = other.getAllPositions(); std::sort(otherPositions.begin(), otherPositions.end());

            equal = equal && (positions.size() == otherPositions.size());
            if (equal) {
                for (size_type i = 0; i < positions.size(); i++) {
                    equal = equal && (positions[i] == otherPositions[i]);
                }
            }

        }
        if (!silent) std::cout << "getAllPositions: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        {
            auto positions = getAllValuedPositions(); std::sort(positions.begin(), positions.end(), sortPairs<size_type, elem_type>());
            auto otherPositions = other.getAllValuedPositions(); std::sort(otherPositions.begin(), otherPositions.end(), sortPairs<size_type, elem_type>());

            equal = equal && (positions.size() == otherPositions.size());
            if (equal) {
                for (size_type i = 0; i < positions.size(); i++) {
                    equal = equal && (positions[i] == otherPositions[i]);
                }
            }

        }
        if (!silent) std::cout << "getAllValuedPositions: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        for (size_t l = 0; l < getLength(); l++) {
            for (size_t r = l; r < getLength(); r++) {

                bool flag = containsElement(l, r);
                bool otherFlag = other.containsElement(l, r);

                equal = equal && (flag == otherFlag);

            }
        }
        if (!silent) std::cout << "containsElement: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

#if 1
        equal = true;
        {
            size_type cnt = countElements();
            size_type otherCnt = other.countElements();
            equal = (cnt == otherCnt);
        }
        if (!silent) std::cout << "countElements: " << (equal ? "OK" : "NOT OK") << std::endl;
        overallEqual = overallEqual && equal;
#endif

        return overallEqual;

    }

    virtual void setNull(size_type i) = 0;

};

#endif //K2TREES_ROWTREE_HPP
