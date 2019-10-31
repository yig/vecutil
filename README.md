Generally useful python code (`vecutil`) for N-dimensional vector math (polylines, rays, spheres, ellipses),
along with some geometric random sampling functions (`sampling`),
line strip distances (`edge_distances`) and smoothed distance field tracing (`smooth_edge_distances`),
and miscellaneous routines (`helpers`).

## Dependencies

* Python >= 2.6 (only `vecutil.py` itself has been tested with Python 3.x)
* numpy
* (optional, for `helpers.friendly_Image*()` functions) PIL (Python Image Library)

---

## Troubleshooting

### If the *extremely optional* Image loading code fails

Following http://jaredforsyth.com/blog/2010/apr/28/accessinit-hash-collision-3-both-1-and-1/
add the following lines at the top of `helpers.py`:

    import sys
    import PIL.Image
    sys.modules['Image'] = PIL.Image
