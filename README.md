Generally useful python code (`vecutil`) for N-dimensional vector math (polylines, rays, spheres, ellipses)
along with some geometric random sampling functions and miscellaneous routines (`sampling`).

## Dependencies

* Python >= 2.6 and < 3
* numpy
* (optional, for `sampling.friendly_Image*()` functions) PIL (Python Image Library)

---

## Troubleshooting

### If the *extremely optional* Image loading code fails

Following http://jaredforsyth.com/blog/2010/apr/28/accessinit-hash-collision-3-both-1-and-1/
add the following lines at the top of `helpers.py`:

    import sys
    import PIL.Image
    sys.modules['Image'] = PIL.Image
