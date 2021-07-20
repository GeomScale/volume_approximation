// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2021 Vaibhav Thakkar

// Contributed and/or modified by Vaibhav Thakkar, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef POSET_H
#define POSET_H

#include <iostream>
#include <vector>


class Poset {
public:
    typedef std::pair<unsigned int, unsigned int> RT; // relation type
    typedef std::vector<RT> RV;  // relations vector

private:
    unsigned int n;  // elements will be from 0 to n-1
    RV order_relations;     // pairs of form a <= b

public:
    Poset(unsigned int _n, RV& _order_relations) : 
        n(_n), order_relations(verify(_order_relations, n))
    {
    }


    // verify if the relations are valid
    static RV verify(const RV& relations, unsigned int n)
    { 
        for (int i = 0; i < relations.size(); i++) {
            if (relations[i].first < 0 || relations[i].first >= n)
                throw "invalid elements in order relations";

            if (relations[i].second < 0 || relations[i].second >= n)
                throw "invalid elements in order relations";
        }

        // TODO: Check if corresponding DAG is actually acyclic
        
        return relations;
    }


    void print()
    {
        std::cout << "Number of elements: " << num_elem() << std::endl;

        unsigned int count = num_relations();
        for(int i=0; i<count; ++i) {
            std::cout <<"A" << order_relations[i].first << " <= A" << order_relations[i].second << std::endl;
        }
    }


    unsigned int num_elem() const
    {
        return n;
    }


    unsigned int num_relations() const
    {
        return order_relations.size();
    }


    RT get_relation(unsigned int idx) const
    {
        return order_relations[idx];
    }

    template <typename NT>
    bool is_in(const std::vector<NT>& pt_coeffs, NT tol=NT(0)) const
    {
        for (int i = 0; i < order_relations.size(); i++) {
            unsigned int a = order_relations[i].first;
            unsigned int b = order_relations[i].second;
            
            if (! (pt_coeffs[a] - pt_coeffs[b] <= tol) ) {
                return false;
            }
        }

        return true;
    }
};

#endif