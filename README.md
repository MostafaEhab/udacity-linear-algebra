# Udacity Linear Algebra Coursework
This repository contains my work for Udacity's Linear Algebra course.

## Course errors
I came across a few mistakes as I went through the lesson videos.

1. In the starter code for `line.py` (Intersections lesson), `n = self.normal_vector` should be changed to `n = self.normal_vector.coords` in the definition of `set_basepoint` and `__str__`. Otherwise it'll throw `TypeError: 'Vector' object is not iterable`
2. In the starter code for `plane.py` (Intersections lesson), `n = self.normal_vector` should be changed to `n = self.normal_vector.coords` in the definition of `set_basepoint` and `__str__`. Otherwise it'll throw `TypeError: 'Vector' object is not iterable`
3. In the starter code for `linsys.py` (Intersections lesson), `indices[i] = p.first_nonzero_index(p.normal_vector)` should be changed to `indices[i] = p.first_nonzero_index(p.normal_vector.coords)` in the definition of `indices_of_first_nonzero_terms_in_each_row`. Otherwise it'll throw `TypeError: 'Vector' object is not iterable`
4. In video 7 of the Intersections lesson, the z-coordinate for the first equation in the last pair should be -7.212, not -7.217.