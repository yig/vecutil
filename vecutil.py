from math import *
from numpy import *

## From http://stackoverflow.com/questions/1237379/how-do-i-set-sys-excepthook-to-invoke-pdb-globally-in-python
import pdb, sys, traceback
def info(type, value, tb):
	traceback.print_exception(type, value, tb)
	pdb.pm()
sys.excepthook = info

## To match pydb.debugger() function:
debugger = pdb.set_trace


# a small set of helper functions to
# call common array creation functions
# these are useful to ensure that
# all arrays are created as double-precision
# floats, no matter what data are provided
# as argument. For example array([1,3,4]) normally returns
# an array with data of type int, but arrayf([1,3,4])
# always creates an array of floats 

kFloatType = float64

def arrayf( arg ):
	return array( arg, kFloatType )
def asarrayf( arg ):
	return asarray( arg, kFloatType )
def zerosf( arg ):
	return zeros( arg, kFloatType )
def onesf( arg ):
	return ones( arg, kFloatType )
def identityf( arg ):
	return identity( arg, kFloatType )
def emptyf( arg ):
	return empty( arg, kFloatType )

def mag2( vec ):
	vec = asarray(vec)
	return dot(vec,vec)
def mag( vec ):
	return sqrt(mag2(vec))
def dir( vec ):
	assert mag(vec) != 0
	vec = asarray(vec)
	return vec * 1./mag(vec)

kEps = 1e-7

def expect( value ):
    if not value:
        print 'Expect failed!'
        import traceback
        traceback.print_stack()
        
        #import pydb
        #pydb.debugger()
    
    #if not value:
    #    import sys
    #    print >> sys.stderr, 'Expect failed!'

def rotate( vec, axis, radians ):
    '''
    Returns 3-vector 'vec' rotated about the given 'axis' by 'radians'.
    
    tested
    '''
    
    vec = asarrayf( vec )
    axis = asarrayf( axis )
    
    ## This code comes from adt/cvec3t.h
    
    c1 = asarray( vec )
    c2 = asarray( axis )
    s = radians
    
    c2l  = mag(c2)
    assert c2l != 0.
    unitc2 = c2/c2l
    c1par  = dot(c1,unitc2)*unitc2
    c1perp = c1-c1par
    return c1par + cos(s)*c1perp + sin(s)*cross(unitc2,c1perp)

def orthogonal( vec ):
    '''
    Returns an arbitrary (deterministic) unit vector orthogonal to 'vec'.
    
    tested
    '''
    
    vec = asarrayf( vec )
    
    ## This code comes from scene2problem.cpp
    
    assert( mag2( vec ) > kEps )
    
    ## Find the smallest component of 'vec'.
    mini = 0
    if( fabs( vec[1] ) < fabs( vec[ mini ] ) ): mini = 1
    if( fabs( vec[2] ) < fabs( vec[ mini ] ) ): mini = 2
    
    ## Add one to it.
    orth = vec.copy()
    orth[ mini ] += 1.
    
    ## Return the normalized cross product of 'vec' and 'orth',
    ## since it must be perpendicular to 'vec'.
    orth = dir( cross( vec, orth ) )
    assert( mag2( orth ) > kEps )
    return orth

def lerp( a, b, t ):
    '''
    Returns the vector 't' along the line from 'a' to 'b'.
    In effect, when 't' is 0, 'a' is returned, and when 't' is 1, 'b' is returned.
    
    tested
    '''
    
    a = asarrayf( a )
    b = asarrayf( b )
    
    return a + t * ( b - a )

def lerps( a, b, ts ):
    '''
    Vectorized lerp.
    Returns [ lerp( a, b, t ) for t in ts ].
    
    tested
    '''
    
    a = asarrayf( a )
    b = asarrayf( b )
    ts = asarrayf( ts )
    
    return a[newaxis,:] + ts[:,newaxis] * ( b - a )[newaxis,:]

def get_axis_angle_rotation( start, end ):
    '''
    Returns the tuple (axis, radians) about which to rotate the vector 'start' so that it becomes
    parallel to, and pointing in the same direction as, the vector 'end'.
    
    tested
    '''
    
    start = asarrayf( start )
    end = asarrayf( end )
    
    ## This code comes from scene2problem.cpp
    
    startn = dir( start )
    endn = dir( end )
    
    ## Now twist (rotate around the eye vector) until 'cur_problems_eye_refl_line_dir'
    ## matches 'problem_eye_refl_line_dir'.
    axis = cross( startn, endn )
    radians = dot( startn, endn )
    
    ## If the rotate to and from vectors have 0 cross product, they're either
    ## the same or opposite.
    if( mag2( axis ) < kEps ):
        ## In this case, the axis is underspecified.  Any orthogonal vector will do.
        axis = orthogonal( startn );
        
        ## The dot product is either 1 or -1
        assert( fabs( radians ) - 1. < kEps )
    
    ## YOTAM: Even though the input vectors to the cross product are
    ##        normalized, without a call to dir(), the magnitude of the
    ##        cross product was off by ~1e-5 on my machine.
    axis = dir( axis )
    
    ## Finally, take the arc cosine of the dot product to get the radians.
    radians = acos( min( 1., max( -1., radians ) ) )
    
    return axis, radians

def slerp( start, end, t ):
    '''
    Returns the vector 't' along the rotation parameterized between [0,1] that rotates
    vector 'start' to vector 'end'.  If 'start' and 'end' differ in length, the length
    is interpolated linearly.
    
    tested
    '''
    
    start = asarrayf( start )
    end = asarrayf( end )
    
    ## This code comes from scene2problem.cpp
    
    axis, radians = get_axis_angle_rotation( start, end )
    
    result = rotate( dir( start ), axis, t * radians ) * lerp( mag( start ), mag( end ), t )
    return result

def reflect( n, r ):
    '''
    Returns the vector 'n' reflected about 'r'.
    In other words, negate the portion of n that is perpendicular to r
    and keep the part of n that is parallel to r.
    
    tested
    '''
    
    n = asarrayf( n )
    r = asarrayf( r )
    
    ## This code comes from scene2problem.cpp
    
    r = dir( r )
    
    ## To reflect a vector n about a unit vector r, we want to negate the portion of n
    ## that is perpendicular to r and keep the part of n that is parallel to r.
    ## The part that is parallel to r is: dot( n, r ) * r
    ## The part that is perpendicular to r is: n - dot( n, r ) * r
    ## Adding the parallel part and negative the perpendicular part, we get:
    ##  dot( n, r ) * r - ( n - dot( n, r ) * r )
    ## <=>
    ##  2 * dot( n, r ) * r - n
    return ( 2 * dot( n, r ) ) * r - n


def get_rotation_matrix( axis, radians ):
    '''
    Given the axis, angle rotation defined by 'axis' and 'radians', returns a
    3x3 numpy.matrix M suitable for multiplying a 3-vector v: M*v.
    
    NOTE: the angle is in radians, not degrees as glRotate*() takes.
    
    EXAMPLE:
        ## A CCW rotation of pi/2 radians about the z axis in the xy plane. Maps <1,0,0> to <0,1,0>.
        dot( asarray( get_rotation_matrix( (0,0,1), pi/2 ) ), (1,0,0) )
        
        ## The matrix transforming a plane with arbitrary normal to the xy plane:
        asarray( get_rotation_matrix( *get_axis_angle_rotation( plane_normal, (0,0,1) ) ) )
    
    tested
    '''
    
    axis = asarrayf( axis )
    
    ## This code comes from hmatrix.h / adt/mat3t.h
    
    u = dir( axis )
    
    u1 = u[0]
    u2 = u[1]
    u3 = u[2]
    
    U = matrix([
        [u1*u1, u1*u2, u1*u3],
        [u2*u1, u2*u2, u2*u3],
        [u3*u1, u3*u2, u3*u3]
        ])
    
    S = matrix([
        [0, -u3, u2],
        [u3,  0, -u1],
        [-u2,  u1, 0]
        ])
    
    result = U + (asmatrix(identity(3)) - U) * cos(radians)  + S * sin(radians)
    return result

def deg2rad( degrees ):
    '''
    Given degrees, returns equivalent radians.
    
    tested
    '''
    return degrees * (pi/180.)

def rad2deg( radians ):
    '''
    Given radians, returns equivalent degrees.
    
    tested
    '''
    return radians * (180./pi)

def sphere_contains_point( sphere, point ):
    '''
    Given a sphere ( 2/3d center, radius ) and a 2/3d point,
    returns whether or not the point is inside the sphere.
    
    tested
    '''
    
    return sphere_contains_points( sphere, (point,) )[0]

def sphere_contains_points( sphere, points ):
    '''
    Given a sphere ( 2/3d center, radius ) and a sequence of 2/3d points,
    returns a sequence of booleans, whether or not each point is inside the sphere.
    
    tested
    '''
    
    sphere_pt = asarrayf( sphere[0] ).ravel()
    radius = sphere[1]
    assert len( sphere_pt ) in (2,3)
    
    if len( points ) == 0: return []
    points = asarrayf( points )
    assert points.shape[1] in (2,3)
    
    assert points.shape[1] == len( sphere_pt )
    
    return (( points - sphere_pt )**2).sum(1) < radius**2

def sphere_normal_at_point( sphere, point ):
    '''
    Given a sphere ( 2/3d center, radius ) and a 2/3d point,
    returns the normal for the point.
    Note that the point need not be on the surface of the sphere.
    
    untested
    '''
    
    return sphere_normal_at_points( sphere, (point,) )[0]

def sphere_normal_at_points( sphere, points ):
    '''
    Given a sphere ( 2/3d center, radius ) and a sequence of 2/3d points,
    returns a sequence of normals for each point.
    Note that the points need not be on the surface of the sphere.
    
    untested
    '''
    
    sphere_pt = asarrayf( sphere[0] ).ravel()
    radius = sphere[1]
    assert len( sphere_pt ) in (2,3)
    
    if len( points ) == 0: return []
    points = asarrayf( points )
    assert points.shape[1] in (2,3)
    
    assert points.shape[1] == len( sphere_pt )
    
    result = points - sphere_pt
    result *= 1 / sqrt( ( result**2 ).sum(1) )
    return result

def ray_sphere_intersection_ts( sphere, ray, ray_is_line = False ):
    '''
    Given a ray ( 2/3d point, 2/3d direction ) and sphere ( 2/3d center, radius ),
    returns a list of 't' values such that ( ray[0] + t * ray[1] ) lies on the sphere.
    If 'ray_is_line' is True, then returned 't' values may be negative.
    
    tested
    '''
    
    ray_pt = asarrayf( ray[0] ).ravel()
    ## Don't normalize ray_dir, because the result is dependent on the scale of ray_dir.
    ray_dir = asarrayf( ray[1] ).ravel()
    assert len( ray_pt ) in (2,3)
    assert len( ray_dir ) in (2,3)
    assert len( ray_pt ) == len( ray_dir )
    
    sphere_pt = asarrayf( sphere[0] ).ravel()
    radius = sphere[1]
    assert len( sphere_pt ) in (2,3)
    
    assert len( sphere_pt ) == len( ray_pt )
    
    return __ray_sphere_intersection_ts( sphere_pt, radius, ray_pt, ray_dir, ray_is_line )
def __ray_sphere_intersection_ts( sphere_pt, radius, ray_pt, ray_dir, ray_is_line ):
    '''
    Parameters already unpacked and type-checked version of ray_sphere_intersection_ts(),
    used internally by this module.
    '''
    ## We want 't' such that | ( ray_pt + t * ray_dir ) - sphere_pt |^2 = radius^2
    ## Let d = ray_dir, c = sphere_pt - ray_pt, and r = radius:
    d = ray_dir
    c = sphere_pt - ray_pt
    r = radius
    ##  | td - c |^2 = r^2
    ## <=>
    ##  (t d.x - c.x )^2 + (t d.y - c.y)^2 + (t d.z - c.z)^2 = r^2
    ## <=>
    ##  t^2 dot( d,d ) - 2t dot( d,c ) + dot( c,c ) = r^2
    ## <=> quadratic equation x = ( -b +- sqrt( b^2 - 4ac ) )/( 2a )
    qa = dot( d,d )
    qb = -2 * dot( d,c )
    qc = dot( c,c ) - r**2
    qrad = qb**2 - 4*qa*qc
    if qrad < 0.: return []
    qrad = sqrt( max( 0., qrad ) )
    ts = [ ( -qb + qrad )/( 2*qa ), ( -qb - qrad )/( 2*qa ) ]
    if not ray_is_line: ts = [ t for t in ts if t >= 0. ]
    #assert len( ts ) > 0
    return ts

def closest_point_to_sphere_on_ray( sphere, ray, ray_is_line = False ):
    '''
    Given a ray ( 2/3d point, 2/3d direction ) and sphere ( 2/3d center, radius ),
    returns a tuple containing the closest point to the sphere along ray ---
    which is the intersection point if they intersect --- and whether or not
    the ray intersected the sphere.
    If 'ray_is_line' is True, then the closest point may lie along the negative ray direction.
    
    tested
    '''
    
    ray_pt = asarrayf( ray[0] ).ravel()
    ray_dir = dir( asarrayf( ray[1] ).ravel() )
    assert len( ray_pt ) in (2,3)
    assert len( ray_dir ) in (2,3)
    assert len( ray_pt ) == len( ray_dir )
    
    sphere_pt = asarrayf( sphere[0] ).ravel()
    radius = sphere[1]
    assert len( sphere_pt ) in (2,3)
    
    assert len( sphere_pt ) == len( ray_pt )
    
    ## Does the ray intersect the sphere?
    ### Adapted from distance_to_line_squared()
    sphere_to_ray = ray_pt - sphere_pt
    ## subtract the component of the ray's direction from the vector towards the ray point
    sphere_to_ray -= dot( sphere_to_ray, ray_dir ) * ray_dir
    
    ## If the magnitude of sphere_to_ray is greater than radius, then the ray doesn't intersect
    ## the sphere; return the point.
    if mag2( sphere_to_ray ) > radius*radius: return sphere_to_ray + sphere_pt, False
    
    #from pydb import debugger
    #debugger()
    
    ts = __ray_sphere_intersection_ts( sphere_pt, radius, ray_pt, ray_dir, ray_is_line )
    assert len( ts ) > 0
    ts = [ ( fabs(t), t ) for t in ts ]
    t = min( ts )[1]
    
    #print 't:', t, 'd:', d, 'ray_pt:', ray_pt
    result = t*ray_dir + ray_pt
    return result, True

def closest_point_on_sphere_to_ray( sphere, ray, ray_is_line = False ):
    '''
    Similar to closest_point_to_sphere_on_ray().
    Differs only in the case where the ray does not intersect the sphere, in which case
    the closest point on the sphere, instead of on the ray, is returned.
    
    tested
    '''
    
    result = closest_point_to_sphere_on_ray( sphere, ray, ray_is_line )
    if result[1]: return result
    else:
        ## The closest point on the sphere is the intersection of the closest point on the ray
        ## with the surface of the sphere.
        ## To get this point, simply return the point radius along the line between the
        ## sphere center and the point.
        sphere_center = asarray( sphere[0] ).ravel()
        sphere_to_pt = asarray( result[0] ).ravel() - sphere_center
        return sphere_center + sphere[1] * dir( sphere_to_pt ), False

def closest_point_on_sphere_to_point( sphere, pt ):
    '''
    Given a point 'pt' and a sphere ( 2/3d center, radius ),
    returns the closest point on the surface of the sphere to 'pt'.
    
    tested
    '''
    
    sphere_center = asarray( sphere[0] ).ravel()
    sphere_to_pt = asarray( pt ).ravel() - sphere_center
    return sphere_center + sphere[1] * dir( sphere_to_pt )

def axes_inverse( axes ):
    '''
    Given the axes of an ellipse ( list of orthogonal but not necessarily unit length 2/3d vectors ),
    returns the inverse transform matrix which transforms points into a space where
    the ellipse is a unit circle/sphere.
    '''
    
    axes = asarray( axes )
    
    ## axes.T is the sphere to ellipse rotation and scale matrix.
    ## We want the inverse.
    
    ## We could do this:
    return linalg.inv( axes.T )
    ## But I think it's more efficient to do this:
    ## UPDATE: It's probably more efficient, but only works if the axes are orthonormal.
    rot_scale_mat = axes / ( axes**2 ).sum(1).reshape( ( axes.shape[0], 1 ) )
    #print 'inverse?:', mag( rot_scale_mat.ravel() - linalg.inv( axes.T ).ravel() )
    #return rot_scale_mat

def ray_ellipse_intersection_ts( ellipse, ray, ray_is_line = False ):
    '''
    Given a ray ( 2/3d point, 2/3d direction ) and ellipse ( 2/3d center, list of 2/3d axes ),
    returns a list of 't' values such that ( ray[0] + t * ray[1] ) lies on the sphere.
    If 'ray_is_line' is True, then returned 't' values may be negative.
    
    untested
    '''
    
    center = asarrayf( ellipse[0] ).ravel()
    axes = asarrayf( ellipse[1] )
    assert len( center ) in (2,3)
    assert axes.shape[0] == axes.shape[1]
    assert axes.shape[1] in (2,3)
    assert len( center ) == axes.shape[1]
    
    ## Copy ray[0], since we will use -=.
    ray_pt = arrayf( ray[0] ).ravel()
    ## Don't normalize ray_dir, because the result is dependent on the scale of ray_dir.
    ray_dir = asarrayf( ray[1] ).ravel()
    assert len( ray_pt ) in (2,3)
    assert len( ray_dir ) in (2,3)
    assert len( ray_pt ) == len( ray_dir )
    assert len( ray_pt ) == len( center )
    
    ## Translate so that the ellipsoid's center is at the origin.
    ray_pt -= center
    
    ## self.axes.T is the sphere to ellipse rotation and scale matrix.
    ## We want the inverse.
    rot_scale_mat = axes_inverse( axes )
    
    ## Transform the ray.
    ray_pt = dot( rot_scale_mat, ray_pt )
    ray_dir = dot( rot_scale_mat, ray_dir )
    
    ## Intersect with a unit sphere.
    ts = ray_sphere_intersection_ts(
        ( zeros(ray_pt.shape), 1 ),
        ( ray_pt, ray_dir ),
        ray_is_line
        )
    return ts

def closest_point_to_ellipse_on_ray( ellipse, ray, ray_is_line = False ):
    '''
    Given a ray ( 2/3d point, 2/3d direction ) and ellipse ( 2/3d center, list of 2/3d axes ),
    returns a tuple containing the closest point to the ellipsoid along ray ---
    which is the intersection point if they intersect --- and whether or not
    the ray intersected the ellipsoid.
    If 'ray_is_line' is True, then the closest point may lie along the negative ray direction.
    
    tested
    '''
    
    center = asarrayf( ellipse[0] ).ravel()
    axes = asarrayf( ellipse[1] )
    assert len( center ) in (2,3)
    assert axes.shape[0] == axes.shape[1]
    assert axes.shape[1] in (2,3)
    assert len( center ) == axes.shape[1]
    
    ## Copy ray[0], since we will use -=.
    ray_pt = arrayf( ray[0] ).ravel()
    ray_dir = dir( asarrayf( ray[1] ) ).ravel()
    assert len( ray_pt ) in (2,3)
    assert len( ray_dir ) in (2,3)
    assert len( ray_pt ) == len( ray_dir )
    assert len( ray_pt ) == len( center )
    
    ## Translate so that the ellipsoid's center is at the origin.
    ray_pt -= center
    
    ## self.axes.T is the sphere to ellipse rotation and scale matrix.
    ## We want the inverse.
    rot_scale_mat = axes_inverse( axes )
    
    ## Transform the ray.
    ray_pt = dot( rot_scale_mat, ray_pt )
    ray_dir = dot( rot_scale_mat, ray_dir )
    
    ## Intersect with a unit sphere.
    closest_pt, intersected = closest_point_to_sphere_on_ray(
        ( zeros(ray_pt.shape), 1 ),
        ( ray_pt, ray_dir ),
        ray_is_line
        )
    
    ## Transform the closest point back to world (ellipse) space.
    closest_pt = dot( axes.T, closest_pt ) + center
    
    return closest_pt, intersected

def closest_point_on_ellipse_to_ray( ellipse, ray, ray_is_line = False ):
    '''
    Similar to closest_point_to_ellipse_on_ray().
    Differs only in the case where the ray does not intersect the ellipse, in which case
    the closest point on the ellipse, instead of on the ray, is returned.
    
    untested
    '''
    
    center = asarrayf( ellipse[0] ).ravel()
    axes = asarrayf( ellipse[1] )
    assert len( center ) in (2,3)
    assert axes.shape[0] == axes.shape[1]
    assert axes.shape[1] in (2,3)
    assert len( center ) == axes.shape[1]
    
    ## Copy ray[0], since we will use -=.
    ray_pt = arrayf( ray[0] ).ravel()
    ray_dir = dir( asarrayf( ray[1] ) ).ravel()
    assert len( ray_pt ) in (2,3)
    assert len( ray_dir ) in (2,3)
    assert len( ray_pt ) == len( ray_dir )
    assert len( ray_pt ) == len( center )
    
    ## Translate so that the ellipsoid's center is at the origin.
    ray_pt -= center
    
    ## self.axes.T is the sphere to ellipse rotation and scale matrix.
    ## We want the inverse.
    rot_scale_mat = axes_inverse( axes )
    
    ## Transform the ray.
    ray_pt = dot( rot_scale_mat, ray_pt )
    ray_dir = dot( rot_scale_mat, ray_dir )
    
    ## Intersect with a unit sphere.
    closest_pt, intersected = closest_point_on_sphere_to_ray(
        ( zeros(ray_pt.shape), 1 ),
        ( ray_pt, ray_dir ),
        ray_is_line
        )
    
    ## Transform the closest point back to world (ellipse) space.
    closest_pt = dot( axes.T, closest_pt ) + center
    
    return closest_pt, intersected

def closest_point_on_ellipse_to_point( ellipse, pt ):
    '''
    Given a point 'pt' and an ellipse ( 2/3d center, list of 2/3d axes ),
    returns the closest point on the surface of the ellipse to 'pt'.
    
    lightly tested
    '''
    
    center = asarrayf( ellipse[0] ).ravel()
    axes = asarrayf( ellipse[1] )
    assert len( center ) in (2,3)
    assert axes.shape[0] == axes.shape[1]
    assert axes.shape[1] in (2,3)
    assert len( center ) == axes.shape[1]
    
    ## Copy pt, since we will use -=.
    pt = arrayf( pt ).ravel()
    assert len( pt ) in (2,3)
    assert len( pt ) == len( center )
    
    ## Translate so that the ellipsoid's center is at the origin.
    pt -= center
    
    ## self.axes.T is the sphere to ellipse rotation and scale matrix.
    ## We want the inverse.
    rot_scale_mat = axes_inverse( axes )
    
    ## Transform the point.
    pt = dot( rot_scale_mat, pt )
    
    ## Find the closest point to a unit sphere.
    closest_pt = closest_point_on_sphere_to_point(
        ( zeros(pt.shape), 1 ),
        pt
        )
    
    ## Transform the closest point back to world (ellipse) space.
    closest_pt = dot( axes.T, closest_pt ) + center
    
    return closest_pt

def ellipse_contains_point( ellipse, point ):
    '''
    Given an ellipse ( 2/3d center, list of 2/3 axes ) and a 2/3d point,
    returns whether or not the point is inside the ellipse.
    
    tested
    '''
    
    return ellipse_contains_points( ellipse, (point,) )[0]

def ellipse_contains_points( ellipse, points ):
    '''
    Given an ellipse ( 2/3d center, list of 2/3 axes ) and a sequence of 2/3d points,
    returns a sequence of booleans, whether or not each point is inside the ellipse.
    
    tested
    '''
    
    center = asarrayf( ellipse[0] ).ravel()
    axes = asarrayf( ellipse[1] )
    ## 2d or 3d
    assert len( center ) in (2,3)
    ## Same number of axes as dimensions.
    assert axes.shape[0] == axes.shape[1]
    ## 2d or 3d
    assert axes.shape[1] in (2,3)
    ## The center and the axes have the same number of dimensions.
    assert len( center ) == axes.shape[1]
    
    if len( points ) == 0: return []
    ## Copy points, since we will use -=.
    points = arrayf( points )
    ## 2d or 3d
    assert points.shape[1] in (2,3)
    ## The point and the ellipse must have the same number of dimensions.
    assert points.shape[1] == len( center )
    
    ## Translate so that the ellipsoid's center is at the origin.
    points -= center
    
    ## axes.T is the sphere to ellipse rotation and scale matrix.
    ## We want the inverse.
    rot_scale_mat = axes_inverse( axes )
    
    ## Transform the points.
    points = dot( rot_scale_mat, points.T ).T
    
    ## Test against a unit sphere.
    return ( points**2 ).sum(1) < 1


def ellipse_normal_at_point( ellipse, point ):
    '''
    Given an ellipse ( 2/3d center, list of 2/3d axes ) and a 2/3d point,
    returns the normal for the point.
    Note that the point need not be on the surface of the ellipse.
    
    untested
    '''
    
    return ellipse_normal_at_points( ellipse, (point,) )[0]

def ellipse_normal_at_points( ellipse, points ):
    '''
    Given an ellipse ( 2/3d center, list of 2/3d axes ) and a sequence of 2/3d points,
    returns a sequence of normals for each point.
    Note that the points need not be on the surface of the ellipse.
    
    untested
    '''
    
    center = asarrayf( ellipse[0] ).ravel()
    axes = asarrayf( ellipse[1] )
    ## 2d or 3d
    assert len( center ) in (2,3)
    ## Same number of axes as dimensions.
    assert axes.shape[0] == axes.shape[1]
    ## 2d or 3d
    assert axes.shape[1] in (2,3)
    ## The center and the axes have the same number of dimensions.
    assert len( center ) == axes.shape[1]
    
    if len( points ) == 0: return []
    ## Copy points, since we will use -=.
    points = arrayf( points )
    ## 2d or 3d
    assert points.shape[1] in (2,3)
    ## The point and the ellipse must have the same number of dimensions.
    assert points.shape[1] == len( center )
    
    ## Translate so that the ellipsoid's center is at the origin.
    points -= center
    
    ## self.axes.T is the sphere to ellipse rotation and scale matrix.
    ## We want the inverse.
    rot_scale_mat = axes_inverse( axes )
    
    ## Transform the points into object-space.
    points = dot( rot_scale_mat, points.T ).T
    
    ## Get normals for the unit sphere.
    sphere_normals = points / sqrt( ( points**2 ).sum(1) )
    assert sphere_normals.shape == points.shape
    
    ## Normals must be transformed using the transpose of the object-to-world matrix.
    result = dot( axes, sphere_normals.T ).T
    assert result.shape == points.shape
    return result

def ellipse_points_for_thetas( ellipse, thetas ):
    '''
    Given an ellipse ( 2d center, list of 2d axes ) and a sequence of theta values,
    returns a sequence of 2d points:
    ellipse[0] + cos( thetas ) * ellipse[1][0] + sin( thetas ) * ellipse[1][1]
    
    The inverse routine is ellipse_thetas_for_points().
    
    tested
    '''
    
    d = 2
    N = len(thetas)
    
    center = asarray( ellipse[0] ).reshape((1,d))
    axis0 = asarray( ellipse[1][0] ).reshape((1,d))
    axis1 = asarray( ellipse[1][1] ).reshape((1,d))
    
    return center + cos( thetas ).reshape((N,1)) * axis0 + sin( thetas ).reshape((N,1)) * axis1

def ellipse_thetas_for_points( ellipse, points ):
    '''
    Given an ellipse ( 2d center, list of 2d axes ) and a sequence of 2d points on the ellipse,
    returns a sequence of theta values 'thetas' such that
    ellipse_points_for_thetas( ellipse, thetas ) returns 'points'.
    
    tested
    '''
    
    center = asarrayf( ellipse[0] ).ravel()
    axes = asarrayf( ellipse[1] )
    ## 2d
    assert len( center ) in (2,)
    ## Same number of axes as dimensions.
    assert axes.shape[0] == axes.shape[1]
    ## 2d
    assert axes.shape[1] in (2,)
    ## The center and the axes have the same number of dimensions.
    assert len( center ) == axes.shape[1]
    
    if len( points ) == 0: return []
    ## Copy points, since we will use -=.
    points = arrayf( points )
    ## 2d
    assert points.shape[1] in (2,)
    ## The point and the ellipse must have the same number of dimensions.
    assert points.shape[1] == len( center )
    
    ## Translate so that the ellipsoid's center is at the origin.
    points -= center
    
    ## axes.T is the sphere to ellipse rotation and scale matrix.
    ## We want the inverse.
    rot_scale_mat = axes_inverse( axes )
    
    ## Transform the points.
    points = dot( rot_scale_mat, points.T ).T
    
    ## Returns atan2 thetas for the unit circle.
    return [ atan2( y,x ) for x,y in points ]

def ellipse_volume( axes ):
    '''
    Given a sequence of ellipse axes ( orthogonal 2/3d vectors ), returns its area/volume.
    
    tested
    '''
    
    axes = asarrayf( axes )
    assert axes.shape[0] == axes.shape[1]
    abc = sqrt( ( axes**2 ).sum(1) ).prod()
    
    if axes.shape[0] == 2:
        return pi * abc
    elif axes.shape[0] == 3:
        return 4./3. * pi * abc

def sphere_volume( sphere ):
    '''
    Given a sphere ( 2/3d center, radius ), returns its area/volume.
    
    untested
    '''
    
    radius = sphere[1]
    
    dim = len( sphere[0] )
    if dim == 2:
        return pi * radius**3
    elif dim == 3:
        return 4./3. * pi * radius**3

def distance_to_line_squared( pt, line ):
    '''
    Given a point 'pt' and a tuple of two points on a line 'line',
    returns the distance from the point to the line, squared.
    
    tested
    '''
    
    ## This code is adapted from navi/.../yigYigCode/Geometry.cpp
    
    pt = asarrayf( pt )
    linePt = asarrayf( line[0] )
    lineVec = dir( asarrayf( line[1] ) - linePt )
    
    ## find the vector towards the point
    toPt = pt - linePt
    
    ## subtract the component of the line's direction from the vector towards the point
    toPt -= dot( toPt, lineVec ) * lineVec
    
    ## the result is the vector from the closest point on the line towards the pt
    return mag2( toPt )

def distance_to_ray_squared( pt, ray ):
    '''
    Given a point 'pt'
    and a tuple of a point and a direction vector 'ray',
    returns the distance from the point to the ray, squared.
    
    tested
    '''
    
    return mag2( asarray(pt) - closest_point_on_ray_to_point( ray, pt ) )

def closest_point_on_ray_to_point( ray, pt ):
    '''
    Given a tuple of a point and a direction vector 'ray'
    and a point 'pt',
    returns the closest point on the ray to the point.
    
    tested
    '''
    
    ray_pt, ray_dir = ray
    ray_pt = asarray( ray_pt )
    ray_dir = asarray( ray_dir )
    pt = asarray( pt )
    
    ## Find the closest point on the ray as a line.
    ## NOTE: This could be more efficient if we had a way to tell
    ##       closest_point_on_line_to_point() to only perform 'max( t, 0 )'.
    closest_pt = closest_point_on_line_to_point( ( ray_pt, ray_pt + ray_dir ), pt )
    ## If that closest point is *behind* the ray,
    ## the ray point itself should be the closest point.
    if dot( closest_pt - ray_pt, ray_dir ) < 0:
        closest_pt = array( ray_pt )
    
    return closest_pt

def closest_point_on_line_to_point( line, pt, line_is_segment = False ):
    '''
    Given a tuple of two points on a line 'line',
    and a point 'pt',
    returns the closest point on the line to the point.
    If 'line_is_segment' is True, then the line is interpreted as a line segment.
    
    tested
    '''
    
    assert len( line ) == 2
    line = [ asarrayf( linept ).ravel() for linept in line ]
    pt = asarrayf( pt )
    ## All points on the line must have the same dimension.
    assert len( set( [ len(linept) for linept in line ] ) ) == 1
    ## The point and the line endpoints must have the same dimension.
    assert len( pt ) == len( line[0] )
    
    ## Get the normalized line direction.
    linevec = line[1] - line[0]
    linelen2 = mag2( linevec )
    
    ## Project the point onto the line.
    t = dot( pt - line[0], linevec )
    
    ## Clip t to [0,1] if 'line_is_segment'
    if line_is_segment: t = min( linelen2, max( 0., t ) )
    
    result = line[0] + (t/linelen2) * linevec
    return result

def closest_distsqr_and_point_and_edge_index_on_line_strip_to_point( line_strip, pt ):
    '''
    Given a line strip (sequence of n-dimensional points) 'linestrip' and n-dimensional point 'pt',
    returns the tuple (
        distance to the closest point on 'line_strip' to 'pt', squared,
        closest point on 'line_strip' to 'pt',
        index of point in 'line_strip' where the edge (index, index+1) contains the closest point )
        ).
    
    used
    '''
    
    assert len( line_strip ) > 0
    pt = asarrayf( pt )
    
    ## Degenerate case, one point in the line strip.
    if len( line_strip ) == 1: return ( mag2( pt - line_strip[0] ), line_strip[0], 0 )
    
    def gen_distsqr_pt_edge( i ):
        closest_pt = closest_point_on_line_to_point( line_strip[i:i+2], pt, line_is_segment = True )
        return ( mag2( pt - closest_pt ), closest_pt, i )
    
    return min( [ gen_distsqr_pt_edge( i ) for i in xrange( len( line_strip )-1 ) ], key = lambda el: ( el[0], el[2] ) )

def closest_distsqr_and_point_and_edge_index_on_line_loop_to_point( line_loop, pt ):
    '''
    Same as closest_distsqr_and_point_and_edge_index_on_line_strip_to_point(), but
    takes a line loop (closed path) instead of a line strip (open path).
    
    NOTE: The index of the closing edge is len( line_loop )-1.
    '''
    ## Append the first point to the end to turn the line loop into a line strip.
    return closest_distsqr_and_point_and_edge_index_on_line_strip_to_point(
        list( line_loop ) + [line_loop[0]], pt
        )

def closest_normal_with_cosangle_to_vector_on_plane( plane_normal, with_vector, I, closest_to, closest_not_farthest = True ):
    '''
    NOTE: Adapted from LaplacianEditingHelper.cpp in the my private git-devel meshopt branch.
    
    Returns a normal (unit vector) such that the cosine of the angle with the vector 'with_vector'
    is 'I', but constrained to lie on the plane whose normal is 'plane_normal'.
    There may be two solutions, in which case we choose the one closer to 'closest_to'.
    
    Q: How to compute this?
    A: Find the intersection of the cone and the circle on the plane.
       If there is no intersection, take the closest point to the cone, which will
       be the normalized projection of the light direction on the 'plane_normal' plane.
       Let normal = orthogonal( plane_normal ),
       Let right = dir( normal x plane_normal ),
       Let l = dir( with_vector ) . normal,
       Let m = dir( with_vector ) . right.
       Then (l,m) is the projection of the light direction onto the normal,right axes,
       and we want to find x,y such that x^2 + y^2 = 1 (x,y are on a unit circle)
       and (x,y) . (l,m) = x*l + y*m = I.
       There may be no solution, in which case the best we can do is the point on the circle
       closest to the cone, which is the unit vector in the direction of (l,m).
       (This assumes I >= 0.)
       There are ordinarily two solutions intersecting the circle, and Maxima gave them to me.
       We choose the solution closest to the 'closest_to' vector.
    
    (%i3) solve( [ x^2 + y^2 = 1, x*l + y*m = I], [x,y] );
    
                            2    2    2                       2    2    2
                  m sqrt(- I  + m  + l ) - l I      l sqrt(- I  + m  + l ) + m I
    (%o3) [[x = - ----------------------------, y = ----------------------------], 
                             2    2                            2    2
                            m  + l                            m  + l
    
                           2    2    2                         2    2    2
                 m sqrt(- I  + m  + l ) + l I        l sqrt(- I  + m  + l ) - m I
            [x = ----------------------------, y = - ----------------------------]]
                            2    2                              2    2
                           m  + l                              m  + l
    
    tested
    '''
    
    plane_normal = asarrayf( plane_normal )
    with_vector = asarrayf( with_vector )
    closest_to = asarrayf( closest_to )
    
    ## I must be positive.  Actually, we can handle negative by negating both I and with_vector.
    expect( I >= 0. )
    if( I < 0 ):
        with_vector = -with_vector
        I = -I
    
    ## I expect I to be less than 1, but if it's greater than one we will choose the "closest"
    ## normal, which will simply be the normal in the same direction as the projection of with_vector.
    expect( I <= 1. )
    
    with_vector = dir( with_vector )
    normal = orthogonal( plane_normal )
    right = dir( cross( plane_normal, normal ) )
    l = dot( with_vector, normal )
    m = dot( with_vector, right )
    
    denom = l*l + m*m
    ## Map the no solution case onto the 1 solution case.
    I = min( I, sqrt( denom ) )
    radical = sqrt( max( denom - I*I, 0. ) )
    
    x1 = -( m*radical - l*I )/denom
    y1 = ( l*radical + m*I )/denom
    
    x2 = ( m*radical + l*I )/denom
    y2 = -( l*radical - m*I )/denom
    
    ## We want the solution closest to the parameter 'closest_to'.
    xc = dot( closest_to, normal )
    yc = dot( closest_to, right )
    if( ( (x1-xc)*(x1-xc) + (y1-yc)*(y1-yc) < (x2-xc)*(x2-xc) + (y2-yc)*(y2-yc) ) == closest_not_farthest ):
        return x1 * normal + y1 * right
    else:
        return x2 * normal + y2 * right

def closest_point_along_ray_with_cosangle_to_vector( ray, I, with_vector, ray_is_line = False ):
    '''
    Given a ( point, vector ) tuple 'ray', returns a point along 'ray' that, when
    viewed as a vector (from the origin), has the given cosine angle 'I' with
    the vector 'with_vector'.
    If 'ray_is_line' is True, then the closest point may lie along the negative ray direction.
    
    tested
    '''
    
    ## We want 't' such that dot( ( ray_point + t * ray_dir ) / | ray_point + t * ray_dir |, with_vector ) = I
    ## Let d = ray_dir, c = ray_point, and v = with_vector.
    d = asarrayf( ray[1] )
    c = asarrayf( ray[0] )
    v = dir( asarrayf( with_vector ) )
    ##  dot( c + td / |c + td|, v ) = I
    ## <=>
    ##  dot( c + td, v ) / |c + td| = I
    ## <=>
    ##  ( dot( c + td, v ) / |c + td| )^2 = I^2
    ## <=>
    ##  dot( c + td, v )^2 / dot( c + td, c + td ) = I^2
    ## <=>
    ##  ( dot( c, v ) + t dot( d, v ) )^2 = I^2 dot( c + td, c + td )
    ## <=>
    ##  ( dot( c, v ) + t dot( d, v ) )^2 = I^2 ( dot( c, c + td ) + dot( td, c + td ) )
    ## <=>
    ##  ( dot( c, v ) + t dot( d, v ) )^2 = I^2 ( dot( c, c ) + dot( c, td ) + dot( td, c ) + dot( td, td ) )
    ## <=>
    ##  ( dot( c, v ) + t dot( d, v ) )^2 = I^2 ( dot( c, c ) + 2 dot( c, td ) + dot( td, td ) )
    ## <=>
    ##  dot( c, v )^2 + 2 t dot( c, v ) dot( d, v ) + t^2 dot( d, v )^2 = I^2 ( dot( c, c ) + 2 dot( c, td ) + dot( td, td ) )
    ## <=>
    ##  dot( c, v )^2 + 2 t dot( c, v ) dot( d, v ) + t^2 dot( d, v )^2 = I^2 ( dot( c, c ) + 2 t dot( c, d ) + t^2 dot( d, d ) )
    ## Let cv = dot( c, v ), dv = dot( d, v ), cc = dot( c, c ), cd = dot( c, d ), dd = dot( d, d )
    cv = dot( c, v )
    dv = dot( d, v )
    cc = dot( c, c )
    cd = dot( c, d )
    dd = dot( d, d )
    ## <=>
    ##  cv^2 + 2 t cv dv + t^2 dv^2 = I^2 ( cc + 2 t cd + t^2 dd )
    ## <=>
    ##  t^2 ( dv^2 - I^2 dd ) + t 2 ( cv dv - I^2 cd ) + ( cv^2 - I^2 cc ) = 0
    ## <=> quadratic equation x = ( -b +- sqrt( b^2 - 4ac ) )/( 2a )
    qa = dv**2 - I**2 * dd
    qb = 2 * ( cv * dv - I**2 * cd )
    qc = cv**2 - I**2 * cc
    qrad = qb**2 - 4*qa*qc
    ## This should be ~zero, or else we may have a degenerate have used the "closest
    ## point outside the sphere".
    assert qrad >= -1e-5
    qrad = sqrt( max( 0., qrad ) )
    ts = ( -qb + qrad )/( 2*qa ), ( -qb - qrad )/( 2*qa )
    if not ray_is_line: ts = [ t for t in ts if t >= 0. ]
    assert len( ts ) > 0
    ts = [ ( fabs(t), t ) for t in ts ]
    t = min( ts )[1]
    
    #print 't:', t, 'd:', d, 'c:', c
    result = c + t * d
    return result
    
    ## If we were interested in the dot product and not the cosine angle (that is, un-normalized),
    ## then we could do the much simpler:
    ## NOTE: This code is untested.
    '''
    ## We want 't' such that dot( ( ray_point + t * ray_dir ), with_vector ) = I
    ## Let d = ray_dir, c = ray_point, and v = with_vector.
    d = asarrayf( ray[1] )
    c = asarrayf( ray[0] )
    v = asarrayf( with_vector )
    ##  dot( c + td, v ) = I
    ## <=>
    ##  dot( c, v ) + t dot( d, v ) = I
    ## <=>
    ##  t dot( d, v ) = I - dot( c, v )
    ## <=>
    ##  t = ( I - dot( c, v ) ) / dot( d, v )
    t = ( I - dot( c, v ) ) / dot( d, v )
    
    return c + t * d
    '''

def project_point_to_plane( point, plane ):
    '''
    Given a point 'point' and a tuple 'plane' consisting of ( point on plane, normal to plane ),
    returns the closest point in the plane to 'point'.
    
    tested
    '''
    
    point = asarrayf( point )
    plane_p = asarrayf( plane[0] )
    plane_n = asarrayf( plane[1] )
    
    assert len( point ) == len( plane_p )
    assert len( plane_p ) == len( plane_n )
    
    ## Make sure the plane normal is normalized.
    plane_n = dir( plane_n )
    
    plane_to_point = point - plane_p
    result = plane_p + ( plane_to_point - dot( plane_to_point, plane_n ) * plane_n )
    return result

def reflect_point_across_plane( point, plane ):
    '''
    Given a point 'point' and a tuple 'plane' consisting of ( point on plane, normal to plane ),
    returns 'point' reflected across 'plane'.
    
    tested
    '''
    
    point = asarrayf( point ).ravel()
    plane_point = project_point_to_plane( point, plane )
    ## Reflection of 'point' across 'plane' is adding to 'point' twice the vector from
    ## 'point' to the closest point on the plane, 'plane_point'.
    ##   point + 2 * ( plane_point - point )
    ## <=>
    ##   2 * plane_point - point
    return 2 * plane_point - point

def ray_plane_intersection_t( ray, plane ):
    '''
    Given a ( point, vector ) tuple 'ray'
    and a tuple 'plane' consisting of ( point on plane, normal to plane ),
    returns the t value such that ray[0] + t*ray[1] intersects the plane
    or None if they do not intersect.
    
    tested
    '''
    
    ## A ray o + t*d intersects a plane passing through point p with normal N when
    ##   dot( o + t*d - p, N ) = 0 = dot( o - p, N ) + t*dot( d, N )
    ## <=>
    ##   t*dot( d, N ) = dot( p - o, N )
    ## <=>
    ##   t = dot( p - o, N ) / dot( d, N )
    
    o = ray[0]
    d = ray[1]
    p = plane[0]
    N = plane[1]
    
    assert len( o ) == len( d )
    assert len( p ) == len( N )
    assert len( o ) == len( p )
    
    denom = dot( d, N )
    if abs( denom ) < kEps:
        #print "Ray and plane don't intersect!"
        return None
    
    t = dot( asarray(p) - o, N ) / float(denom)
    return t

def ray_plane_intersection( ray, plane, ray_is_line = False ):
    '''
    Given a ( point, vector ) tuple 'ray'
    and a tuple 'plane' consisting of ( point on plane, normal to plane ),
    returns the closest point along ray that intersects the plane
    or None if they do not intersect.
    If 'ray_is_line' is True, then the closest point may lie along the negative ray direction.
    
    tested
    '''
    
    t = ray_plane_intersection_t( ray, plane )
    if t is None or t < 0. and not ray_is_line: return None
    
    return asarray( ray[0] ) + float(t) * asarray( ray[1] )

def print_matlab( arr ):
    '''
    Given an array or matrix, prints it out in such a way that it can be copied and pasted into
    Matlab.
    
    tested
    '''
    
    arr = asarray( arr )
    N = 1
    if len( arr.shape ) > 1: N = arr.shape[1]
    
    print '[ ',
    count = 0
    for v in arr.ravel():
        if count == N:
            print '; ',
            count = 0
        print '%s ' % (v,),
        count += 1
    print ']'

def lsq_ellipse( pts ):
    '''
    Given a list of points 'pts', returns a least squares fit of an ellipse to the points as a
    tuple ( ellipse_center, [ ellipse_axis1, ellipse_axis2 ] ).
    The long axis is ellipse_axis1.
    
    tested
    '''
    
    try:
        from fitellipse import fitellipse
        #z, a, b, alpha = fitellipse( asarray( pts ), 'nonlinear', constraint = 'trace' )
        #print_matlab( pts )
        #z, a, b, alpha = fitellipse( asarray( pts ), 'linear', constraint = 'trace' )
        #z, a, b, alpha = fitellipse( asarray( pts ), 'linear', constraint = 'bookstein' )
        z, a, b, alpha = fitellipse( asarray( pts ), 'linear' )
        #z, a, b, alpha = fitellipse( asarray( pts ) )
        #print z, a, b, alpha
        
        center = asarray( z ).reshape((2,))
        axis1 = a * cos( alpha ), a * sin( alpha )
        axis2 = b * -sin( alpha ), b * cos( alpha )
        
        if b > a: axis1, axis2 = axis2, axis1
        
        return center, asarray( [ asarray( axis1 ), asarray( axis2 ) ] )
    
    except linalg.linalg.LinAlgError:
        return None
    except RuntimeError, e:
        return None
        
        #print e
        
        pt0 = asarray( pts[0] )
        pt1 = asarray( pts[-1] )
        r = .5 * mag( pt1 - pt0 )
        return ( .5 * (pt0 + pt1), asarray( [ r * array((1,0)), r * array((0,1)) ] ) )

def ellipse_bbox( ellipse ):
    '''
    Returns the bounding box of an ellipse.
    Each ellipse is a tuple ( 2/3d center, list of 2/3d axes ),
    and the bounding box result is the tuple ( minimum coordinates, maximum coordinates ).
    
    used
    '''
    
    ellipse = asarray( ellipse[0] ), asarray( ellipse[1] )
    ## center +- ( take the absolute value of the axes, and then take the max coordinate among all axes )
    return ellipse[0] + abs( ellipse[1] ).max(0), ellipse[0] - abs( ellipse[1] ).max(0)

def bbox_union( bbox1, bbox2 ):
    '''
    Returns the union of two bounding boxes 'bbox1', 'bbox2', where each bounding box
    parameter is the tuple ( maximum coordinates, minimum coordinates ),
    where the coordinates are, for example, xyz.
    
    used
    '''
    
    return (
        asarray( [ bbox1[0], bbox2[0] ] ).max(0),
        asarray( [ bbox1[1], bbox2[1] ] ).min(0)
        )

def bbox_intersection( bbox1, bbox2 ):
    '''
    Returns the intersection of two bounding boxes 'bbox1', 'bbox2', where each bounding box
    parameter is the tuple ( maximum coordinates, minimum coordinates ),
    where the coordinates are, for example, xyz.
    
    used
    '''
    
    return (
        asarray( [ bbox1[0], bbox2[0] ] ).min(0),
        asarray( [ bbox1[1], bbox2[1] ] ).max(0)
        )

def ellipses_intersect( ellipse1, ellipse2 ):
    '''
    Returns whether or not 'ellipse1' and 'ellipse2' intersect.
    Each ellipse is a tuple ( 2/3d center, list of 2/3d axes ).
    Currently only 2d ellipses are implemented, so center must be 2d and there must be two 2d axes.
    
    tested
    '''
    
    ellipse1 = ( asarray( ellipse1[0] ).ravel(), asarray( ellipse1[1] ) )
    ellipse2 = ( asarray( ellipse2[0] ).ravel(), asarray( ellipse2[1] ) )
    
    assert len( ellipse1[0] ) == 2
    assert ellipse1[1].shape[1] == 2
    assert len( ellipse2[0] ) == 2
    assert ellipse2[1].shape[1] == 2
    
    ellipse1_bbox = ellipse_bbox( ellipse1 )
    ellipse2_bbox = ellipse_bbox( ellipse2 )
    both_bbox = bbox_intersection( ellipse1_bbox, ellipse2_bbox )
    
    maxx = both_bbox[0][0]
    maxy = both_bbox[0][1]
    minx = both_bbox[1][0]
    miny = both_bbox[1][1]
    
    ## Sample a uniform grid inside the bounding box
    N = 1000
    ## We want to sample space uniformly.
    ## This means N samples should cover the area (maxx-minx) * (maxy-miny) uniformly.
    ## Let NX be the number of samples in x.
    ## Then we want NX * (maxy-miny)/(maxx-minx) samples in y
    ## and N = NX * NX * (maxy-miny)/(maxx-minx).
    ## So
    NX = max( 3, sqrt( N / (maxy-miny) * (maxx-minx) ) )
    NY = max( 3, NX * (maxy-miny)/(maxx-minx) )
    #print 'NX, NY:', NX, NY
    xs = linspace( minx, maxx, int(round( NX )) )
    ys = linspace( miny, maxy, int(round( NY )) )
    ## This line performs the cartesian product of 'xs' and 'ys', giving us a list of 2d points.
    pts = asarray( map( lambda a: a.ravel(), meshgrid( xs, ys ) ) ).T
    ellipse1_contains = ellipse_contains_points( ellipse1, pts )
    ellipse2_contains = ellipse_contains_points( ellipse2, pts )
    return logical_and( ellipse1_contains, ellipse2_contains ).any()

def ellipses_overlap_area( ellipse1, ellipse2, accuracy = 1e-3 ):
    '''
    Returns the area of overlap between 'ellipse1' and 'ellipse2'.
    Each ellipse is a tuple ( 2/3d center, list of 2/3d axes ).
    Currently only 2d ellipses are implemented, so center must be 2d and there must be two 2d axes.
    
    Optional parameters:
        accuracy, which determines the floating point accuracy to achieve (default 1e-3)
    
    used
    '''
    
    assert accuracy > 0.
    oo_accuracy = 1. / accuracy
    
    ellipse1 = ( asarray( ellipse1[0] ).ravel(), asarray( ellipse1[1] ) )
    ellipse2 = ( asarray( ellipse2[0] ).ravel(), asarray( ellipse2[1] ) )
    
    assert len( ellipse1[0] ) == 2
    assert ellipse1[1].shape[1] == 2
    assert len( ellipse2[0] ) == 2
    assert ellipse2[1].shape[1] == 2
    
    ellipse1_bbox = ellipse_bbox( ellipse1 )
    ellipse2_bbox = ellipse_bbox( ellipse2 )
    both_bbox = bbox_intersection( ellipse1_bbox, ellipse2_bbox )
    
    maxx = both_bbox[0][0]
    maxy = both_bbox[0][1]
    minx = both_bbox[1][0]
    miny = both_bbox[1][1]
    
    ## Sample a uniform grid inside the bounding box
    N = max( 1, int( round( oo_accuracy ) ) )
    ## We want to sample space uniformly.
    ## This means N samples should cover the area (maxx-minx) * (maxy-miny) uniformly.
    ## Let NX be the number of samples in x.
    ## Then we want NX * (maxy-miny)/(maxx-minx) samples in y
    ## and N = NX * NX * (maxy-miny)/(maxx-minx).
    ## So
    NX = max( 3, sqrt( N / (maxy-miny) * (maxx-minx) ) )
    NY = max( 3, NX * (maxy-miny)/(maxx-minx) )
    ## Round to the nearest integer.
    NX = int( round( NX ) )
    NY = int( round( NY ) )
    #print 'NX, NY:', NX, NY
    xs = linspace( minx, maxx, NX )
    ys = linspace( miny, maxy, NY )
    ## This line performs the cartesian product of 'xs' and 'ys', giving us a list of 2d points.
    pts = asarray( map( lambda a: a.ravel(), meshgrid( xs, ys ) ) ).T
    ellipse1_contains = ellipse_contains_points( ellipse1, pts )
    ellipse2_contains = ellipse_contains_points( ellipse2, pts )
    both_contain_samples = sum( logical_and( ellipse1_contains, ellipse2_contains ) )
    ## Each sample represents a little square with area ( (maxx-minx)/NX ) * ( (maxy-miny)/NY )
    return both_contain_samples * ( (maxx-minx)/NX ) * ( (maxy-miny)/NY )

def center_of_mass_of_ellipses_intersection( ellipse1, ellipse2, accuracy = 1e-3, timeout = 1e20 ):
    '''
    Returns the center of mass of the intersection of 'ellipse1' and 'ellipse2' as a 2/3d point.
    Each ellipse is a tuple ( 2/3d center, list of 2/3d axes ).
    NOTE: This routine will not terminate if the ellipses don't overlap.
    
    Optional parameters:
        accuracy, which determines the floating point accuracy to achieve (default 1e-3)
        timeout, which determines the number of seconds after which to terminate via exception with two args, number of hits and current center-of-mass
    
    used
    '''
    
    assert accuracy > 0.
    oo_accuracy = 1. / accuracy
    
    if timeout != 1e20: raise NotImplementedError
    
    ellipse1 = ( asarray( ellipse1[0] ).ravel(), asarray( ellipse1[1] ) )
    ellipse2 = ( asarray( ellipse2[0] ).ravel(), asarray( ellipse2[1] ) )
    
    ellipse1_bbox = ellipse_bbox( ellipse1 )
    ellipse2_bbox = ellipse_bbox( ellipse2 )
    both_bbox = bbox_intersection( ellipse1_bbox, ellipse2_bbox )
    
    #from pydb import debugger
    #debugger()
    
    ## Random sample many times
    
    ## The number of points to test at once
    M = max( 1, int( round( oo_accuracy ) ) )
    
    tries = 0
    hits = 0
    com = zeros( ellipse1[0].shape )
    while hits < oo_accuracy and (tries == 0 or float(hits)/tries > accuracy):
        
        ## Test M points at once using ellipse_contains_points()
        #pt = [ random.uniform( minval, maxval ) for minval, maxval in zip( both_bbox[1], both_bbox[0] ) ]
        pts = random.random_sample( ( M, len( ellipse1[0] ) ) )
        for i, (minval, maxval) in enumerate( zip( both_bbox[1], both_bbox[0] ) ):
            pts[:,i] *= (maxval - minval)
            pts[:,i] += minval
        
        new_hits = logical_and(
            ellipse_contains_points( ( ellipse1[0], ellipse1[1] ), pts ),
            ellipse_contains_points( ( ellipse2[0], ellipse2[1] ), pts )
            )
        
        hits += sum( new_hits )
        com += pts[ new_hits ].sum(0)
        tries += M
        
        #print 'hits:', hits
    
    if hits == 0: raise RuntimeError, "Ellipses don't intersect"
    com *= 1./hits
    return com

def closest_point_on_tilted_circle2d_to_line( tilted_circle, line ):
    '''
    Given a tilted 2d circle in space (3d center, radius, 3d normal)
    and a 3d line defined by a pair of 3d points,
    returns the closest 3d point on the circle to the line.
    
    untested
    '''
    
    center, radius, normal = tilted_circle
    
    assert len( center ) == 3
    assert float( radius ) == radius
    assert len( normal ) == 3
    
    normal = asarray( normal ).ravel()
    
    axis0 = radius * orthogonal( normal )
    axis1 = radius * dir( cross( normal, axis0 ) )
    
    tilted_ellipse = ( center, ( axis0, axis1 ) )
    return closest_point_on_tilted_ellipse2d_to_line( tilted_ellipse, line )

def closest_point_on_tilted_ellipse2d_to_line( tilted_ellipse, line ):
    '''
    Given a tilted 2d ellipse in space (3d center, list of two 3d axes)
    and a 3d line defined by a pair of 3d points,
    returns the closest 3d point on the ellipse to the line.
    
    untested
    '''
    
    tilted_ellipse = ( asarray( tilted_ellipse[0] ).ravel(), asarray( tilted_ellipse[1] ) )
    cos_theta, sin_theta = closest_point_coeffs_on_tilted_ellipse2d_to_line( tilted_ellipse, line )
    return tilted_ellipse[0] + cos_theta * tilted_ellipse[1][0] + sin_theta * tilted_ellipse[1][1]

def closest_point_coeffs_on_tilted_ellipse2d_to_line( tilted_ellipse, line ):
    '''
    Given a tilted 2d ellipse in space (3d center, list of two 3d axes)
    and a 3d line defined by a pair of 3d points,
    returns the ellipse axis coefficients of the closest point on the ellipse to the line,
    where the ellipse axis coefficients are a pair of values (a,b) such that
        tilted_ellipse[0] + a * tilted_ellipse[1][0] + b * tilted_ellipse[1][1]
    is the closest point.
    
    NOTE: If the ellipse is parameterized in radians from 0 to 2*pi, where
            0 is the ellipse center + first axis,
            pi/2 is the ellipse center + second axis,
            pi is the ellipse center - first axis, and
            3*pi/2 is the ellipse center 0 second axis,
          then the parameter in radians is atan2( b, a ).
    
    tested
    '''
    
    tilted_ellipse = ( asarray( tilted_ellipse[0] ).ravel(), asarray( tilted_ellipse[1] ) )
    line = ( asarray( line[0] ).ravel(), asarray( line[1] ).ravel() )
    
    assert len( tilted_ellipse[0] ) == 3
    assert tilted_ellipse[1].shape == (2,3)
    assert len( line[0] ) == 3
    assert len( line[0] ) == len( line[1] )
    
    ## 1 Transform (by rotation) into the space where the line is a point in the
    ##   xy plane (its direction is parallel to the z axis).
    ## 2 Find the closest point on the transformed ellipse projected onto the xy plane.
    
    ## 1
    ray = ( line[0], line[1] - line[0] )
    ray2point_rot_mat = array( get_rotation_matrix( *get_axis_angle_rotation( ray[1], (0,0,-1) ) ) )
    ray_point = dot( ray2point_rot_mat, ray[0] )
    
    tilted_ellipse_rot = (
        dot( ray2point_rot_mat, tilted_ellipse[0] ),
        dot( ray2point_rot_mat, tilted_ellipse[1].T ).T
        )
    tilted_ellipse2d = (
        tilted_ellipse_rot[0][:2],
        tilted_ellipse_rot[1][:,:2]
        )
    
    ## 2
    closest_pt_ellipse2d = closest_point_on_ellipse_to_point( tilted_ellipse2d, ray_point[:2] )
    ## NOTE: The desired 3d point is *not* obtained by inverse rotation of the 2d point,
    ##         dot( ray2point_rot_mat.T, append( tilted_ellipse2d_pt, 0. ) ),
    ##       because we have stripped away the z component.
    ##       Instead, we have to find the coefficients for each axis, which is the same in both
    ##       2d and 3d.  The coefficients can then be used to find the 3d point.
    ## NOTE: We have to invert the ellipse axes --- we can't just take the transpose ---
    ##       because they aren't orthonormal.
    cos_theta, sin_theta = linalg.solve(
        asarray( tilted_ellipse2d[1] ).T,
        closest_pt_ellipse2d - tilted_ellipse2d[0]
        )
    return cos_theta, sin_theta

def line_loop_is_CCW( line_loop ):
    '''
    Given a 2d line loop (closed sequence of 2d points) 'line_loop',
    returns true if the line loop is oriented in a counter clock-wise (CCW) fashion.
    
    used
    '''
    
    line_loop = asarray( line_loop )
    N = line_loop.shape[0]
    
    ## A 2d line loop is an N-by-2 array, where N is at least three.
    assert len( line_loop.shape ) == 2
    assert N >= 3
    assert line_loop.shape[1] == 2
    
    ## Based on comp.graphics.algorithms FAQ 2.07:
    ##     http://www.faqs.org/faqs/graphics/algorithms-faq/
    ##     http://maven.smith.edu/~orourke/Code/polyorient.C
    
    ## First find the bottom vertex (and right-most, in case of a tie).
    ## Find the smallest y coordinate.
    miny = min( line_loop[:,1] )
    ## Find all points with y coordinate equal to the smallest y coordinate.
    candidates = [
        ( pt, index )
        for index, pt in enumerate( line_loop )
        if pt[1] == miny
        ]
    ## Get the index of the point with the largest x coordinate, among the candidates.
    br = max( candidates, key = lambda el: el[0][0] )[1]
    
    a = line_loop[ br - 1 ]
    b = line_loop[ br ]
    c = line_loop[ ( br + 1 ) % N ]
    
    two_area = (
        a[0] * b[1] - a[1] * b[0] +
        a[1] * c[0] - a[0] * c[1] +
        b[0] * c[1] - c[0] * b[1]
        )
    return two_area >= 0.
