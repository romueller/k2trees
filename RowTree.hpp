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

#ifndef K2TREES_ROWTREE_HPP
#define K2TREES_ROWTREE_HPP

#include "Utility.hpp"

/**
 * Representation of a (weighted / valued) subset of a universe.
 * One-dimensional adaptation of K2Tree.
 *
 * A RowTree of length n describes a subset S of the universe [0 : n - 1].
 *
 * The weights / values of an entry in the relation are of type E
 * and one value is designated the null element ("element is not in the set").
 *
 * The data structure is static (with the exception of the setNull() method).
 */
template<typename E>
class RowTree {

public:
    // weights / values
    typedef E elem_type;

    // input type (collection of valued universe members)
    typedef std::vector<std::pair<size_type, elem_type>> list_type;

    virtual ~RowTree() { }

    // returns the size of the universe (length of the "row", n)
    virtual size_type getLength() = 0;

    // returns the null element
    virtual elem_type getNull() = 0;


    // checks whether i is in S
    virtual bool isNotNull(size_type i) = 0;

    // returns the value of i, if the element is in S, null otherwise
    virtual elem_type getElement(size_type i) = 0;

    // returns the smallest (left-most) element in S, or a value >= n if S is empty
    virtual size_type getFirst() = 0;

    // returns the values of all elements i in S with l <= i <= r
    virtual std::vector<elem_type> getElementsInRange(size_type l, size_type r) = 0;

    // returns the positions of all elements i in S with l <= i <= r
    virtual std::vector<size_type> getPositionsInRange(size_type l, size_type r) = 0;

    // returns the positions and values of all elements i in S with l <= i <= r
    virtual list_type getValuedPositionsInRange(size_type l, size_type r) = 0;

    // returns the values of all elements in S
    virtual std::vector<elem_type> getAllElements() = 0;

    // returns the positions of all elements in S
    virtual std::vector<size_type> getAllPositions() = 0;

    // returns the positions and values of all elements in S
    virtual list_type getAllValuedPositions() = 0;

    // checks whether S contains an element i with l <= i <= r
    virtual bool containsElement(size_type l, size_type r) = 0;

    // counts the number of elements in S
    virtual size_type countElements() = 0;


    // creates a deep copy
    virtual RowTree* clone() const = 0;

    // prints the parameters (and contents) of the RowTree
    virtual void print(bool all = false) = 0;

    // compares the RowTree with a given vector representation
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

    // compares the RowTree with another RowTree
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

    // sets the value of element i to null, i.e. removes it from the set
    virtual void setNull(size_type i) = 0;

};

#endif //K2TREES_ROWTREE_HPP
