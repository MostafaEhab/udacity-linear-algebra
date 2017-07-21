# Udacity Linear Algebra Coursework
This repository contains my work for Udacity's Linear Algebra course.

## Course errors
I came across a few mistakes as I went through the lesson videos.

1. In the starter code for `line.py` (Intersections lesson), `n = self.normal_vector` should be changed to `n = self.normal_vector.coords` in the definition of `set_basepoint` and `__str__`. Otherwise it'll throw `TypeError: 'Vector' object is not iterable`