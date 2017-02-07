function[withNaN] = zero2NaN( A )

withNaN = A;
withNaN( withNaN == 0) = NaN;