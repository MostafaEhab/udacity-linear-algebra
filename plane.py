'''
This module is a collection of functions for working with planes.
'''
from decimal import Decimal, getcontext
from vector import Vector

# Decimal precision
getcontext().prec = 30

class Plane(object):
    '''
    This class is a collection of functions for working with planes.
    '''

    # Error messages
    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

    def __init__(self, normal_vector=None, constant_term=None):
        self.dimension = 3

        if not normal_vector:
            all_zeros = ['0'] * self.dimension
            normal_vector = Vector(all_zeros)
        self.normal_vector = normal_vector

        if not constant_term:
            constant_term = Decimal('0')
        self.constant_term = Decimal(constant_term)

        self.set_basepoint()


    def set_basepoint(self):
        try:
            n = self.normal_vector
            c = self.constant_term
            basepoint_coords = ['0']*self.dimension

            initial_index = Plane.first_nonzero_index(n)
            initial_coefficient = n[initial_index]

            basepoint_coords[initial_index] = c/initial_coefficient
            self.basepoint = Vector(basepoint_coords)

        except Exception as e:
            if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                self.basepoint = None
            else:
                raise e


    def __str__(self):

        num_decimal_places = 3

        def write_coefficient(coefficient, is_initial_term=False):
            coefficient = round(coefficient, num_decimal_places)
            if coefficient % 1 == 0:
                coefficient = int(coefficient)

            output = ''

            if coefficient < 0:
                output += '-'
            if coefficient > 0 and not is_initial_term:
                output += '+'

            if not is_initial_term:
                output += ' '

            if abs(coefficient) != 1:
                output += '{}'.format(abs(coefficient))

            return output

        n = self.normal_vector

        try:
            initial_index = Plane.first_nonzero_index(n)
            terms = [write_coefficient(n[i], is_initial_term=(i==initial_index)) + 'x_{}'.format(i+1)
                     for i in range(self.dimension) if round(n[i], num_decimal_places) != 0]
            output = ' '.join(terms)

        except Exception as e:
            if str(e) == self.NO_NONZERO_ELTS_FOUND_MSG:
                output = '0'
            else:
                raise e

        constant = round(self.constant_term, num_decimal_places)
        if constant % 1 == 0:
            constant = int(constant)
        output += ' = {}'.format(constant)

        return output

    def __eq__(self, p):
        '''
        Checks if two planes are coincident.

        :param Plane p: second plane
        :return: whether the two planes are coincident
        :rtype: Boolean
        '''
        # If one of the normal vectors zero vector
        if self.normal_vector.is_zero():
            if not p.normal_vector.is_zero():
                return False
            else:
                diff = self.constant_term - p.constant_term
                return MyDecimal(diff).is_near_zero()
        elif p.normal_vector.is_zero():
            return False

        # If two lines not parallel
        if not self.is_parallel(p):
            return False

        # Vector between points on the two lines
        basepoint_diff = self.basepoint.substract(p.basepoint)

        # This vector should be orthogonal to the normal vector of both lines
        # Only need to check one, as we know the lines are parallel
        return basepoint_diff.is_orthogonal(self.normal_vector)

    def __iter__(self):
        self.current = 0
        return self

    def __next__(self):
        if self.current >= len(self.normal_vector):
            raise StopIteration
        else:
            current_value = self.normal_vector[self.current]
            self.current += 1
            return current_value

    def __len__(self):
        return len(self.normal_vector)

    def __getitem__(self, i):
        return self.normal_vector[i]

    @staticmethod
    def first_nonzero_index(iterable):
        for k, item in enumerate(iterable):
            if not MyDecimal(item).is_near_zero():
                return k
        raise Exception(Plane.NO_NONZERO_ELTS_FOUND_MSG)

    def is_parallel(self, p):
        '''
        Checks if two planes are parallel.

        :param Plane p: second plane
        :return: whether the two lines are parallel
        :rtype: Boolean
        '''
        return self.normal_vector.is_parallel(p.normal_vector)

    def is_orthogonal(self, p):
        '''
        Checks if two planes are orthogonal.

        :param Plane p: second plane
        :return: whether the two lines are orthogonal
        :rtype: Boolean
        '''
        return self.normal_vector.is_orthogonal(p.normal_vector)

class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps