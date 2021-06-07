template<typename R>
NumberFieldElement<R> NumberFieldElement<R>::operator-() const
{
  NumberFieldElement<R> ;
}

template<typename R>
NumberFieldElement<R> NumberFieldElement<R>::operator+(const NumberFieldElement<R> & ) const;

template<typename R>
NumberFieldElement<R> NumberFieldElement<R>::operator-(const NumberFieldElement<R> & ) const;

template<typename R>
NumberFieldElement<R> NumberFieldElement<R>::operator*(const NumberFieldElement<R> & ) const;
template<typename R>
NumberFieldElement<R> NumberFieldElement<R>::operator/(const NumberFieldElement<R> & ) const;

template<typename R>
NumberFieldElement<R> NumberFieldElement<R>::operator*(const R & ) const;
template<typename R>
NumberFieldElement<R> NumberFieldElement<R>::operator/(const R & ) const;
template<typename R>
NumberFieldElement<R> & NumberFieldElement<R>::operator+=(const NumberFieldElement<R> & );
template<typename R>
NumberFieldElement<R> & NumberFieldElement<R>::operator-=(const NumberFieldElement<R> & );
template<typename R>
NumberFieldElement<R> & NumberFieldElement<R>::operator*=(const NumberFieldElement<R> & );
template<typename R>
NumberFieldElement<R> & NumberFieldElement<R>::operator/=(const NumberFieldElement<R> & );
template<typename R>
NumberFieldElement<R>& NumberFieldElement<R>::operator*=(const R & );
template<typename R>
NumberFieldElement<R>& NumberFieldElement<R>::operator/=(const R & );
