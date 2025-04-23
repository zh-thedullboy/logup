#pragma once

#include "../goldilocks/src/goldilocks_base_field.hpp"
#include <cstdint>
#include <vector>
#include <array>


#define FIELD_EXTENSION 2

/*
This is a field extension 2 of the goldilocks:
Prime: 0xFFFFFFFF00000001
Irreducible polynomial: x^2 - 7
*/

class Goldilocks2
{
public:
    // typedef Goldilocks::Element Element[FIELD_EXTENSION];
    using Element = std::array<Goldilocks::Element, 2>;

private:
    static const Element ZERO;
    static const Element ONE;
    static const Element NEGONE;

public:
    static constexpr uint64_t p = GOLDILOCKS_PRIME;

    static inline const Element &zero() {return ZERO;}
    static inline void zero(Element &result)
    {
        result[0] = Goldilocks::zero();
        result[1] = Goldilocks::zero();
    }

    static inline const Element &one() {return ONE;}
    static inline void one(Element &result)
    {
        result[0] = Goldilocks::one();
        result[1] = Goldilocks::zero();
    }
    
    static inline const Element &negone() {return NEGONE;}
    static inline void negone(Element &result)
    {
        result[0] = Goldilocks::negone();
        result[1] = Goldilocks::zero();
    }

    static void copy(Element &dst, const Element &src)
    {
        for (uint64_t i = 0; i < FIELD_EXTENSION; i++)
        {
            Goldilocks::copy(dst[i], src[i]);
        }
    };
    static void copy(Element *dst, const Element *src)
    {
        for (uint64_t i = 0; i < FIELD_EXTENSION; i++)
        {
            Goldilocks::copy((*dst)[i], (*src)[i]);
        }
    };
    static inline Element fromU64(uint64_t in1[FIELD_EXTENSION]){
        Element result;
        Goldilocks2::fromU64(result, in1);
        return result;
    }
    static inline void fromU64(Element &result, uint64_t in1[FIELD_EXTENSION])
    {
        for (uint64_t i = 0; i < FIELD_EXTENSION; i++)
        {
            result[i] = Goldilocks::fromU64(in1[i]);
        }
    }
    static inline Element fromS32(int32_t in1[FIELD_EXTENSION]){
        Element result;
        Goldilocks2::fromS32(result, in1);
        return result;
    }
    static inline void fromS32(Element &result, int32_t in1[FIELD_EXTENSION])
    {
        for (uint64_t i = 0; i < FIELD_EXTENSION; i++)
        {
            result[i] = Goldilocks::fromS32(in1[i]);
        }
    }
    static inline std::array<uint64_t, 2> toU64(const Element &in1){
        std::array<uint64_t, 2> res;
        Goldilocks2::toU64(res.data(), in1);
        return res;
    }
    static inline void toU64(uint64_t *result, const Element &in1)
    {
        for (uint64_t i = 0; i < FIELD_EXTENSION; i++)
        {
            result[i] = Goldilocks::toU64(in1[i]);
        }
    }

    static inline std::string toString(const Element &in1, int radix = 10)
    {
        std::string res;
        toString(res, in1, radix);
        return res;
    }
    static inline void toString(std::string &result, const Element &in1, int radix = 10)
    {
        for (uint64_t i = 0; i < FIELD_EXTENSION; i++)
        {
            result += Goldilocks::toString(in1[i]);
            (i != FIELD_EXTENSION - 1) ? result += " , " : "";
        }
    }
    static inline void toString(std::string (&result)[FIELD_EXTENSION], const Element &in1, int radix = 10)
    {
        for (uint64_t i = 0; i < FIELD_EXTENSION; i++)
        {
            result[i] = Goldilocks::toString(in1[i]);
        }
    }
    static std::string toString(const Element *in1, const uint64_t size, int radix)
    {
        std::string result = "";
        for (uint64_t i = 0; i < size; i++)
        {
            result += Goldilocks2::toString(in1[i], 10);
            result += "\n";
        }
        return result;
    }
    static inline void fromString(Element &result, const std::string (&in1)[FIELD_EXTENSION], int radix = 10)
    {
        for (uint64_t i = 0; i < FIELD_EXTENSION; i++)
        {
            result[i] = Goldilocks::fromString(in1[i]);
        }
    }

    // ======== ADD ========
    static inline void add(Element &result, const Element &a, const uint64_t &b)
    {
        result[0] = a[0] + Goldilocks::fromU64(b);
        result[1] = a[1];
    }
    static inline void add(Element &result, const Element &a, const Goldilocks::Element b)
    {
        result[0] = a[0] + b;
        result[1] = a[1];
    }
    static inline void add(Element &result, const Goldilocks::Element a, const Element &b)
    {
        add(result, b, a);
    }
    static inline void add(Element &result, const Element &a, const Element &b)
    {
        for (uint64_t i = 0; i < FIELD_EXTENSION; i++)
        {
            result[i] = a[i] + b[i];
        }
    }

    // ======== SUB ========
    static inline void sub(Element &result,const Element &a,const uint64_t &b)
    {
        result[0] = a[0] - Goldilocks::fromU64(b);
        result[1] = a[1];
    }
    static inline void sub(Element &result,const Goldilocks::Element a,const Element &b)
    {
        result[0] = a - b[0];
        result[1] = Goldilocks::neg(b[1]);
    }
    static inline void sub(Element &result,const Element &a,const Goldilocks::Element b)
    {
        result[0] = a[0] - b;
        result[1] = a[1];
    }
    static inline void sub(Element &result,const Element &a,const Element &b)
    {
        for (uint64_t i = 0; i < FIELD_EXTENSION; i++)
        {
            result[i] = a[i] - b[i];
        }
    }

    // ======== NEG ========
    static inline void neg(Element &result, Element &a)
    {
        sub(result, (Element &)zero(), a);
    }

    // ======== MUL ========
    // static inline void mul(Element *result, Element *a, Element *b)
    // {
    //     mul(*result, *a, *b);
    // }
    static inline void mul(Element &result,const Element &a,const Element &b)
    {
        Goldilocks::Element A = (a[0] + a[1]) * (b[0] + b[1]);
        Goldilocks::Element B = a[0] * b[0];
        Goldilocks::Element C = a[1] * b[1];
        Goldilocks::Element D = A - B - C;

        result[0] = B + C * Goldilocks::fromU64(7ull);
        result[1] = D;
    };
    static inline void mul(Element &result,const Element &a,const Goldilocks::Element &b)
    {
        result[0] = a[0] * b;
        result[1] = a[1] * b;
    }
    static inline void mul(Element &result,const Goldilocks::Element a,const Element &b)
    {
        mul(result, b, a);
    }
    static inline void mul(Element &result,const Element &a,const uint64_t b)
    {
        result[0] = a[0] * Goldilocks::fromU64(b);
        result[1] = a[1] * Goldilocks::fromU64(b);
    }

    // ======== DIV ========
    static inline void div(Element &result,const Element &a,const Goldilocks::Element b)
    {
        Goldilocks::Element b_inv = Goldilocks::inv(b);
        mul(result, a, b_inv);
    }
    static inline void div(Element &result,const Element &a,const Element &b)
    {
        Goldilocks2::Element b_inv; 
        Goldilocks2::inv(b_inv, b);
        mul(result, a, b_inv);
    }

    // ======== MULSCALAR ========
    static inline void mulScalar(Element &result,const Element &a,const std::string &b)
    {
        result[0] = a[0] * Goldilocks::fromString(b);
        result[1] = a[1] * Goldilocks::fromString(b);
    }

    // ======== SQUARE ========
    static inline void square(Element &result, Element &a)
    {
        mul(result, a, a);
    }

    // ======== INV ========
    static inline void inv(Element *result,const Element *a)
    {
        inv(*result, *a);
    }
    static inline void inv(Element &result,const Element &a)
    {
        Goldilocks::Element x = a[0];
        Goldilocks::Element y = a[1];
        Goldilocks::Element denom = x * x - y * y * Goldilocks::fromU64(7ull);
        Goldilocks::Element dinv = Goldilocks::inv(denom);
        result[0] = x * dinv;
        result[1] = Goldilocks::neg(y) * dinv;
    }
};


// /*
//     Operator Overloading
// */
// inline Goldilocks2::Element operator+(const Goldilocks2::Element &in1, const Goldilocks2::Element &in2) { return Goldilocks2::add(in1, in2); }
// inline Goldilocks2::Element operator*(const Goldilocks2::Element &in1, const Goldilocks2::Element &in2) { return Goldilocks2::mul(in1, in2); }
// inline Goldilocks2::Element operator-(const Goldilocks2::Element &in1, const Goldilocks2::Element &in2) { return Goldilocks2::sub(in1, in2); }
// inline Goldilocks2::Element operator/(const Goldilocks2::Element &in1, const Goldilocks2::Element &in2) { return Goldilocks2::div(in1, in2); }
// inline bool operator==(const Goldilocks2::Element &in1, const Goldilocks2::Element &in2) { return in1[0] == in2[0] && in1[1] == in2[1]; }
// inline Goldilocks2::Element operator-(const Goldilocks2::Element &in1) { return Goldilocks2::neg(in1); }
// inline Goldilocks2::Element operator+(const Goldilocks2::Element &in1) { return in1; }