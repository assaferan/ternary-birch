// implementation file for header Genus.h

template<typename R, size_t Rank>
Z Genus<R,rank>::get_mass(const QuadForm<R>& q,
                  const std::vector<PrimeSymbol<R>>& symbols)
    {		 
        Z mass = 2 * this->disc;
        Z a = q.h() * q.h() - 4 * q.a() * q.b();
        Z b = -q.a() * this->disc;

        for (const PrimeSymbol<R>& symb : symbols)
        {
            mass *= (symb.p + Math<Z>::hilbert_symbol(a, b, symb.p));
            mass /= 2;
            mass /= symb.p;
        }

        return mass;
    }