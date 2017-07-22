'''
This module is a collection of functions for working with linear systems.
'''
from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from plane import Plane

# Decimal precision
getcontext().prec = 30

class LinearSystem(object):
    '''
    This class is a collection of functions for working with linear systems.
    '''

    # Error messages
    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    def __len__(self):
        return len(self.planes)

    def __getitem__(self, i):
        return self.planes[i]

    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i+1,p) for i,p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret

    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        plane1 = self.planes[row_to_add]
        normal_vector1 = plane1.normal_vector
        constant_term1 = plane1.constant_term

        plane2 = self.planes[row_to_be_added_to]
        normal_vector2 = plane2.normal_vector
        constant_term2 = plane2.constant_term

        new_normal_vector = plane1.normal_vector.scalar_multiply(coefficient).add(normal_vector2)
        new_constant_term = constant_term1 * coefficient + constant_term2

        self.planes[row_to_be_added_to] = Plane(normal_vector=new_normal_vector, constant_term=new_constant_term)

    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)
        num_variables = self.dimension

        indices = [-1] * num_equations

        for i,p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices

    def swap_rows(self, row1, row2):
        self.planes[row2], self.planes[row1] = self.planes[row1], self.planes[row2]

    def multiply_coefficient_and_row(self, coefficient, row):
        plane = self.planes[row]
        normal_vector = plane.normal_vector
        constant_term = plane.constant_term

        new_normal_vector = plane.normal_vector.scalar_multiply(coefficient)
        new_constant_term = constant_term * coefficient

        self.planes[row] = Plane(normal_vector=new_normal_vector, constant_term=new_constant_term)

    def compute_triangular_form(self):
        system = deepcopy(self)

        n_equations = len(system)
        n_variables = self.dimension
        col = 0

        for row in range(n_equations):
            while col < n_variables:
                c = MyDecimal(system[row].normal_vector[col])
                if c.is_near_zero():
                    swapped = system.swap_with_row_below_for_nonzero_coefficient(row, col)
                    if not swapped:
                        col += 1
                        continue
            
                system.clear_coeffs_below(row, col)
                col += 1
                break

        return system

    def swap_with_row_below_for_nonzero_coefficient(self, row, col):
        num_equations = len(self)

        for k in range(row + 1, num_equations):
            coeff = MyDecimal(self[k].normal_vector[col])
            if not coeff.is_near_zero():
                self.swap_rows(row, k)
                return True

        return False

    def compute_rref(self):
        triangular_form = self.compute_triangular_form()

        num_equations = len(triangular_form)
        pivot_indices = triangular_form.indices_of_first_nonzero_terms_in_each_row()

        for row in range(num_equations)[::-1]:
            pivot_var = pivot_indices[row]
            if pivot_var < 0:
                continue
            triangular_form.scale_row_to_make_coefficient_equal_one(row, pivot_var)
            triangular_form.clear_coefficients_above(row, pivot_var)

        return triangular_form

    def clear_coeffs_below(self, row, col):
        num_equations = len(self)
        beta = MyDecimal(self[row].normal_vector[col])

        for k in range(row + 1, num_equations):
            n = self[k].normal_vector
            gamma = n[col]
            alpha = -gamma / beta
            self.add_multiple_times_row_to_row(alpha, row, k)

    def clear_coefficients_above(self, row, col):
        for k in range(row)[::-1]:
            n = self[k].normal_vector
            alpha = -(n[col])
            self.add_multiple_times_row_to_row(alpha, row, k)

    def scale_row_to_make_coefficient_equal_one(self, row, col):
        n = self[row].normal_vector
        beta = Decimal('1.0') / n[col]
        self.multiply_coefficient_and_row(beta, row)

    def raise_excepion_if_contradictory_equation(self):
        for plane in self.planes:
            try:
                plane.first_nonzero_index(plane.normal_vector)

            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    constant_term = MyDecimal(plane.constant_term)
                    if not constant_term.is_near_zero():
                        raise Exception(self.NO_SOLUTIONS_MSG)
                else:
                    raise e

    def raise_excepion_if_too_few_pivots(self):
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        num_pivots = sum([1 if index >= 0 else 0 for index in pivot_indices])
        num_variables = self.dimension

        if num_pivots < num_variables:
            raise Exception(self.INF_SOLUTIONS_MSG)

    def compute_solution(self):
        try:
            return self.do_gaussian_elimination_and_parametrization()

        except Exception as e:
            if (str(e) == self.NO_SOLUTIONS_MSG or
                str(e) == self.INF_SOLUTIONS_MSG):
                return str(e)
            else:
                raise e

    def do_gaussian_elimination_and_parametrization(self):
        rref = self.compute_rref()
        rref.raise_excepion_if_contradictory_equation()

        direction_vectors = rref.extract_direction_vectors_for_parametrization()
        basepoint = rref.extract_basepoint_for_parametrization()

        return Parametrization(basepoint, direction_vectors)

    def extract_direction_vectors_for_parametrization(self):
        num_variables = self.dimension
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        free_variable_indices = set(range(num_variables)) - set(pivot_indices)

        direction_vectors = []

        for free_var in free_variable_indices:
            vector_coords = [0] * num_variables
            vector_coords[free_var] = 1
            for index, plane in enumerate(self.planes):
                pivot_var = pivot_indices[index]
                if pivot_var < 0:
                    break
                vector_coords[pivot_var] = -plane.normal_vector[free_var]

            direction_vectors.append(Vector(vector_coords))

        return direction_vectors

    def extract_basepoint_for_parametrization(self):
        num_variables = self.dimension
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()

        basepoint_coords = [0] * num_variables

        for index, plane in enumerate(self.planes):
            pivot_var = pivot_indices[index]
            if pivot_var < 0:
                break
            basepoint_coords[pivot_var] = plane.constant_term

        return Vector(basepoint_coords)

class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps

class Parametrization(object):

    BASEPT_AND_DIR_VECTORS_MUST_BE_IN_SAME_DIM = ('The basepoint and \
            direction vectors should all live in the same dimension')

    def __init__(self, basepoint, direction_vectors):

        self.basepoint = basepoint
        self.direction_vectors = direction_vectors
        self.dimension = self.basepoint.dimension

        try:
            for v in direction_vectors:
                assert v.dimension == self.dimension

        except AssertionError:
            raise Exception(self.BASEPT_AND_DIR_VECTORS_MUST_BE_IN_SAME_DIM)

    def __str__(self):

        output = ''
        for coord in range(self.dimension):
            output += 'x_{} = {} '.format(coord + 1,
                                          round(self.basepoint[coord], 3))
            for free_var, vector in enumerate(self.direction_vectors):
                output += '+ {} t_{}'.format(round(vector[coord], 3),
                                             free_var + 1)
            output += '\n'
        return output