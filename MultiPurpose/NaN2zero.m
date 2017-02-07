function[with0] = NaN2zero( A )

with0 = A;
with0( isnan(with0) ) = 0;