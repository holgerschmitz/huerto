/*
 * visitor.hpp
 *
 *  Created on: 30 Jan 2020
 *      Author: Holger Schmitz
 */

#ifndef HUERTO_UTIL_VISITOR_HPP_
#define HUERTO_UTIL_VISITOR_HPP_

// Copied from Andrew Durward's answer
// https://stackoverflow.com/a/11802080

// Visitor template declaration
template<typename... Types>
class Visitor;

// specialization for single type    
template<typename T>
class Visitor<T> {
public:
    virtual void visit(T & visitable) = 0;
};

// specialization for multiple types
template<typename T, typename... Types>
class Visitor<T, Types...> : public Visitor<Types...> {
public:
    // promote the function(s) from the base class
    using Visitor<Types...>::visit;

    virtual void visit(T & visitable) = 0;
};

template<typename... Types>
class Visitable {
public:
    virtual void accept(Visitor<Types...>& visitor) = 0;
};

template<typename Derived, typename... Types>
class VisitableImpl : public Visitable<Types...> {
public:
    virtual void accept(Visitor<Types...>& visitor) {
        visitor.visit(static_cast<Derived&>(*this));
    }
};

#endif // HUERTO_UTIL_VISITOR_HPP_