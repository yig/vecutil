## This file is from: ~/Work/TAU/hasm-code/experiments/src/helpers.py

from math import *
from numpy import *

## For predictable debugging.
kDebugRandom = False
if kDebugRandom:
    import random as modrandom
    random = modrandom.Random()
    random.seed(31337)
    ## I was debugging something and this seed value would reproduce it.
    #random.seed(1)
else:
    import random

def normalized_vector3( v ):
    v = array( v, dtype = float ).ravel()
    assert v.shape == (3,)
    
    kEps = 1e-5
    vMag2 = dot( v,v )
    assert vMag2 > kEps
    v *= 1./sqrt(vMag2)
    
    return v

def normalize_vectors( vecs, eps = 1e-8 ):
    '''
    Returns the array-like sequence of vectors, normalized.  If the vector has magnitude less
    than optional parameter 'eps', it will be returned unchanged.
    '''
    
    vecs = asarray( vecs, dtype = float )
    assert len( vecs.shape ) == 2
    
    vec_norms = sqrt( ( vecs**2 ).sum( axis = 1 ) )
    vec_norms[ vec_norms < eps ] = 1.
    return vecs / vec_norms[:,newaxis]

def vector_is_normalized( vec, eps = 1e-5 ):
    '''
    Returns True if the vector parameter 'vec'
    is already normalized.
    Returns False otherwise.
    '''
    
    return vectors_are_normalized( [vec], eps )

def vectors_are_normalized( vecs, eps = 1e-5 ):
    '''
    Returns True if all vectors in the sequence of vectors
    'vecs' are already normalized.
    Returns False otherwise.
    '''
    
    assert len( vecs ) >= 1
    vecs = asarray( vecs )
    
    return max( abs( 1. - ( vecs**2 ).sum( axis = 1 ) ) ) < eps

def reflect_vectors_across_vector( vectors, r ):
    '''
    Given an array 'vectors' of vectors, and 
    a vector 'r' dimension is the same as the last dimension of 'vectors.shape[-1]',
    returns the vectors in 'vectors' reflected across 'r'.
    
    tested:
    >>> reflect_vector_across_normals( [( 0, 0, 1 )], ( 0, 0, 1 ) )
    [( 0, 0, 1 )]
    >>> reflect_vector_across_normals( [( 0, 0, 1 )], ( 0, 0, 2 ) )
    [( 0, 0, 1 )]
    >>> reflect_vector_across_normals( [( 0, 0, 2 )], ( 0, 0, 1 ) )
    [( 0, 0, 2 )]
    >>> reflect_vector_across_normals( [( 0, 0, 2 )], ( 0, 0, 2 ) )
    [( 0, 0, 2 )]
    >>> reflect_vector_across_normals( [( 1, 0, 0 )], ( 0, 0, 1 ) )
    [( -1, 0, 0 )]
    >>> reflect_vector_across_normals( [( 1, 0, 1 )], ( 0, 0, 2 ) )
    [( -1, 0, 1 )]
    >>> reflect_vector_across_normals( [( .5, 1, 0 )], ( 0, 0, 1 ) )
    [( -.5, -1, 0 )]
    >>> reflect_vector_across_normals( [( .5, 1, 1 )], ( 0, 0, 2 ) )
    [( -.5, -1, 1 )]
    >>> reflect_vector_across_normals( [( 1, 2, 0 )], ( 0, 0, 1 ) )
    [( -1, -2, 0 )]
    >>> reflect_vector_across_normals( [( 1, 2, 1 )], ( 0, 0, 2 ) )
    [( -1, -2, 1 )]
    
    tested more:
    >>> assert abs( reflect_vectors_across_vector( [( .5, 1, 1 )], ( 0, 0, 2 ) ) - asarray([( -.5, -1, 1 )]) ).max() < 1e-5
    >>> assert abs( reflect_vectors_across_vector(
        [( 0, 0, 1 ), ( 0, 0, 2 ), ( 1, 0, 0 ), ( 1, 0, 1 ), ( .5, 1, 0 ), ( .5, 1, 1 ), ( 1, 2, 0 ), ( 1, 2, 1 )],
        [(0,0,1)]
        )
        -
        asarray(
        [[( 0, 0, 1 ), ( 0, 0, 2 ), ( -1, 0, 0 ), ( -1, 0, 1 ), ( -.5, -1, 0 ), ( -.5, -1, 1 ), ( -1, -2, 0 ), ( -1, -2, 1 )]]
        )
        ).max() < 1e-5
    >>> assert abs( reflect_vectors_across_vector(
        [( 0, 0, 1 ), ( 0, 0, 2 ), ( 1, 0, 0 ), ( 1, 0, 1 ), ( .5, 1, 0 ), ( .5, 1, 1 ), ( 1, 2, 0 ), ( 1, 2, 1 )],
        [(0,0,2)]
        )
        -
        asarray(
        [[( 0, 0, 1 ), ( 0, 0, 2 ), ( -1, 0, 0 ), ( -1, 0, 1 ), ( -.5, -1, 0 ), ( -.5, -1, 1 ), ( -1, -2, 0 ), ( -1, -2, 1 )]]
        )
        ).max() < 1e-5
    >>> assert abs( reflect_vectors_across_vector(
        [[( 0, 0, 1 ), ( 0, 0, 2 ), ( 1, 0, 0 ), ( 1, 0, 1 )],
        [( .5, 1, 0 ), ( .5, 1, 1 ), ( 1, 2, 0 ), ( 1, 2, 1 )]],
        [(0,0,1)]
        )
        -
        asarray(
        [[[( 0, 0, 1 ), ( 0, 0, 2 ), ( -1, 0, 0 ), ( -1, 0, 1 )],
        [( -.5, -1, 0 ), ( -.5, -1, 1 ), ( -1, -2, 0 ), ( -1, -2, 1 )]]]
        )
        ).max() < 1e-5
    '''
    
    ## We want to move each vector in 'vectors' twice the distance
    ## to their projection on 'ur == normalized( r ) == r / sqrt( dot( r, r ) )':
    ##   vectors + 2 * ( dot( vectors, ur )*ur - vectors )
    ## =
    ##   vectors + 2 * dot( vectors, ur )*ur - 2 * vectors
    ## =
    ##   2 * dot( vectors, ur )*ur - vectors
    ## =
    ##   dot( vectors, ur )*2*ur - vectors
    ## =
    ##   dot( vectors, r )*2*r/dot(r,r) - vectors
    
    ## 'r' must be a non-zero vector.
    ## UPDATE: Enforce this to be a float, because the (r/dot(r,r)) may truncate.
    r = asarray( r, dtype = float ).squeeze()
    assert len( r.shape ) == 1
    assert abs(r).max() > 1e-5
    
    ## 'vectors' must be an N-by-len(r) array.
    ## NOTE: The fact that 'vectors' can't be simply another len(r)-array
    ##       makes sure the programmer can't screw up the order of the
    ##       parameters to this function.
    vectors = asarray( vectors )
    assert len( vectors.shape ) > 1
    assert vectors.shape[-1] == r.shape[0]
    
    result = dot( vectors, 2*r )[...,newaxis] * ( r/dot(r,r) ) - vectors
    
    #print result
    assert result.shape == vectors.shape
    
    return result

def geometric_median_normal( normals ):
    '''
    Given a sequence of normalized k-dimensional vectors `normals`,
    returns the one whose sum of angles with all other normals is smallest.
    
    tested
    >>> geometric_median_normal( [ (0,0,1), (sqrt(2)/2,0,sqrt(2)/2), (-sqrt(2)/2,0,-sqrt(2)/2) ] )
    (1,0,0)
    '''
    
    normals = asarray( normals )
    ## The first dimension is the number of normals.
    ## The second dimension is k.
    assert len( normals.shape ) == 2
    
    ## Assume input is normalized.
    # normals = normalize_normals( normals )
    total_deviation_from_normalized = abs( ( asarray( normals )**2 ).sum( axis = 1 ) - 1. ).sum()
    if total_deviation_from_normalized > 1e-5:
        print( "WARNING: geometric_median_normal() called with non-normalized normals. Sum of magnitudes - 1:", total_deviation_from_normalized )
    
    assert len( normals ) > 0
    
    ## Compute all pairs of dot products.
    ## Take arc cosine to get angles.
    all_way_angles = arccos( dot( normals, normals.T ).clip(-1,1) )
    assert all_way_angles.min() >= 0.
    ## Sum them up.
    ## The normal corresponding to the smallest value is the one we want.
    all_way_angles = all_way_angles.sum( axis = 0 )
    best_normal_index = argmin( all_way_angles )
    return normals[ best_normal_index ]

def sample_value( r = None ):
    '''
    Returns a uniformly sampled value in [0,1].
    
    Takes optional parameter 'r', an instance of random.Random.
    '''
    
    return sample_value_in_range( 0, 1, r )

def sample_value_in_range( range_start, range_end, r = None ):
    '''
    Returns a uniformly sampled value in [range_start, range_end].
    
    Takes optional parameter 'r', an instance of random.Random.
    '''
    
    if r is None: r = random
    return r.uniform( range_start, range_end )

def sample_N_values_more_than_distance_from_each_other( N, distance, value_generator = sample_value ):
    '''
    Given an integer 'N', and
    a scalar 'distance',
    returns a sequence of 'N' values, where the distance between any two values
    is greater than 'distance'.
    
    Takes an optional parameter 'value_generator', a function that takes no arguments
    and returns a scalar.  The default uniformly samples [0,1].
    
    tested
    This should be slightly greater than .1:
    min( [ ( lambda a: abs( a[1:] - a[:-1] ).min() )( asarray( sorted( sample_N_values_more_than_distance_from_each_other( 5, .1 ) ) ) ) for i in xrange(1000) ] )
    '''
    
    assert N == int(N)
    assert N > 0
    assert distance == float(distance)
    
    while True:
        values = asarray( [ value_generator() for i in xrange(N) ] )
        if all_values_are_more_than_distance_from_each_other( values, distance ): return values

def all_values_are_more_than_distance_from_each_other( values, distance ):
    '''
    Given a list of N scalars 'values'
    and a scalar 'distance',
    returns False if the distance between any two values is less than 'distance', and
    returns True otherwise.
    '''
    
    ## This function is ill-defined if there are less than 2 values.
    assert len( values ) >= 2
    values = asarray( values ).squeeze()
    assert len( values.shape ) == 1
    
    assert distance == float(distance)
    
    for i in xrange(values.shape[0]):
        for j in xrange(i+1,values.shape[0]):
            ## If the distance between the two value is
            ## less than or equal to distance, return False.
            if abs( values[i] - values[j] ) <= distance: return False
    
    return True

def sample_value_more_than_distance_from_v( distance, v, value_generator = sample_value ):
    '''
    Given a scalar 'distance' and a scalar 'v',
    returns a random value more than 'distance' from 'v'.
    
    Takes an optional parameter 'value_generator', a function that takes no arguments
    and returns a scalar.  The default uniformly samples [0,1].
    
    tested
    e.g. should be in [0,.25) and (.75,1].
    ## This should be slightly less than 1
    asarray( [ sample_value_more_than_distance_from_v( .25, .5 ) for i in xrange(1000) ] ).max()
    ## This should be slightly more than 0
    asarray( [ sample_value_more_than_distance_from_v( .25, .5 ) for i in xrange(1000) ] ).min()
    ## This should be slightly less than .25
    asarray([ v for v in asarray( [ sample_value_more_than_distance_from_v( .25, .5 ) for i in xrange(1000) ] ) if v < .5 ]).max()
    ## This should be slightly more than .75
    asarray([ v for v in asarray( [ sample_value_more_than_distance_from_v( .25, .5 ) for i in xrange(1000) ] ) if v > .5 ]).min()
    '''
    
    assert distance == float(distance)
    
    while True:
        n = value_generator()
        ## If the distance between the two value is more than distance, we're done.
        if abs( v - n ) > distance: return n

def sample_point_inside_unit_circle( r = None ):
    if r is None: r = random
    x, y = 1,1
    while x**2 + y**2 > 1.:
        x,y = r.uniform( -1, 1 ), r.uniform( -1, 1 )
    return x,y

def sample_point_inside_unit_circle_far_from_point( x0,y0, too_close, r = None ):
    if r is None: r = random
    x, y = 1,1
    while ( x*x + y*y > 1. ) or ( (x-x0)**2 + (y-y0)**2 < too_close**2 ):
        x,y = r.uniform( -1, 1 ), r.uniform( -1, 1 )
    return x,y

def sample_3D_point( r = None ):
    '''
    Returns a uniformly sampled, 3D point in [0,1]^3.
    
    Takes optional parameter 'r', an instance of random.Random.
    '''
    
    if r is None: r = random
    
    p = zeros(3)
    p[:] = r.uniform( 0, 1 ), r.uniform( 0, 1 ), r.uniform( 0, 1 )
    return p

def sample_2D_integer_point_in_mask( mask, r = None ):
    '''
    Given a 2D boolean array 'mask',
    returns a uniformly sampled, 2D point 'p' whose coordinates are
    integers such that 'mask[ p ]' evaluates to True.
    
    Takes optional parameter 'r', an instance of random.Random.
    '''
    
    if r is None: r = random
    
    assert len( mask.shape ) == 2
    assert mask.any()
    
    p = zeros( 2, dtype = int )
    while True:
        p[:] = r.randint( 0, mask.shape[0]-1 ), r.randint( 0, mask.shape[1]-1 )
        if mask[ p[0], p[1] ]: return p

def sample_normal( r = None ):
    '''
    Returns a uniformly sampled, 3D unit normal vector.
    
    Takes optional parameter 'r', an instance of random.Random.
    
    tested
    e.g. should be 1:
    (sample_normal()**2).sum()
    '''
    
    if r is None: r = random
    
    kEps = 1e-5
    n = zeros(3)
    nMag2 = 0.
    while nMag2 < kEps:
        n[:] = r.uniform( -1, 1 ), r.uniform( -1, 1 ), r.uniform( -1, 1 )
        nMag2 = dot( n, n )
    
    n *= 1./sqrt(nMag2)
    return n

def sample_normal_within_cos_phi_of_v( cos_phi, v, normal_generator = sample_normal ):
    '''
    Given the cosine of an angle phi 'cos_phi' and a 3D vector 'v',
    returns a 3D unit normal vector where the angle between the returned vector
    and 'v' is less than phi.
    
    Takes an optional parameter 'normal_generator', a function that takes no arguments
    and returns a unit 3D vector.  The default uniformly samples 3D unit vectors.
    
    tested
    e.g. should be greater than sqrt(.5):
    asarray( [ sample_normal_within_cos_phi_of_v( cos(pi/4), (0,0,1) ) for i in xrange(1000) ] )[:,2].min()
    -asarray( [ sample_normal_within_cos_phi_of_v( cos(pi/4), (-1,0,0) ) for i in xrange(1000) ] )[:,0].max()
    '''
    
    assert cos_phi == float(cos_phi)
    cos_phi = clip( cos_phi, -1, 1 )
    
    v = normalized_vector3( v )
    
    while True:
        n = normal_generator()
        ## If the angle between two vectors is less than phi,
        ## then the dot product between them is greater than cos( phi ).
        if dot( v,n ) > cos_phi: return n

def sample_normal_more_than_cos_phi_from_v( cos_phi, v, normal_generator = sample_normal ):
    '''
    Given the cosine of an angle phi 'cos_phi' and a 3D vector 'v',
    returns a 3D unit normal vector where the angle between the returned vector
    and 'v' is greater than phi.
    
    Takes an optional parameter 'normal_generator', a function that takes no arguments
    and returns a unit 3D vector.  The default uniformly samples 3D unit vectors.
    
    tested
    e.g. should be less than sqrt(.5):
    asarray( [ sample_normal_more_than_cos_phi_from_v( cos(pi/4), (0,0,1) ) for i in xrange(1000) ] )[:,2].max()
    -asarray( [ sample_normal_more_than_cos_phi_from_v( cos(pi/4), (-1,0,0) ) for i in xrange(1000) ] )[:,0].min()
    '''
    
    assert cos_phi == float(cos_phi)
    cos_phi = clip( cos_phi, -1, 1 )
    
    v = normalized_vector3( v )
    
    while True:
        n = normal_generator()
        ## If the angle between two vectors is less than phi,
        ## then the dot product between them is greater than cos( phi ).
        if dot( v,n ) < cos_phi: return n

def sample_front_facing_friendly_normal( xy_radius, normal_generator = sample_normal ):
    '''
    Returns a 3D unit normal vector whose z value is non-negative and who x,y values
    lie within the circle at the origin with radius 'xy_radius'.
    
    Takes an optional parameter 'normal_generator', a function that takes no arguments
    and returns a unit 3D vector.  The default uniformly samples 3D unit vectors.
    
    tested
    e.g. should be greater than 0:
    asarray([ sample_front_facing_friendly_normal( 1. ) for i in xrange(1000) ])[:,2].min()
    e.g. should be greater than sqrt(.5):
    asarray([ sample_front_facing_friendly_normal( sqrt(.5) ) for i in xrange(1000) ])[:,2].min()
    '''
    
    assert xy_radius == float(xy_radius)
    
    while True:
        n = normal_generator()
        if dot( n[:2], n[:2] ) < xy_radius**2 and n[2] >= 0.:
            return n

def all_normals_are_more_than_cos_theta_from_each_other( normals, cos_theta ):
    '''
    Given a list of N-dimensional unit normals 'normals'
    and the cosine of an angle 'cos_theta',
    returns False if the angle between any two normals is less than theta, and
    returns True otherwise.
    
    NOTE: Assumes that 'normals' is normalized!
    '''
    
    ## We can't squeeze() 'normals', because it might have only 1 element!
    ## UPDATE: Actually, this function is ill-defined if there are less than 2 normals.
    assert len( normals ) >= 2
    normals = asarray( normals ).squeeze()
    
    ## 3D normals
    assert len( normals.shape ) == 2
    #assert normals.shape[1] == 3
    
    assert cos_theta == float(cos_theta)
    assert cos_theta > -1.
    
    ## Are the vectors normalized?
    assert vectors_are_normalized( normals )
    ## Normalize the vectors.
    #import helpers
    #normals = helpers.normalize_vectors( normals )
    
    for i in xrange(normals.shape[0]):
        for j in xrange(i+1,normals.shape[0]):
            ## If the angle between two vectors is greater than theta,
            ## then the dot product between them is less than cos( theta ).
            if dot( normals[i], normals[j] ) >= cos_theta: return False
    
    return True

def all_normals_are_less_than_cos_theta_from_each_other( normals, cos_theta ):
    '''
    Given a list of N-dimensional unit normals 'normals'
    and the cosine of an angle 'cos_theta',
    returns True if the angle between every two normals is less than theta, and
    returns False otherwise.
    
    NOTE: Assumes that 'normals' is normalized!
    '''
    
    ## We can't squeeze() 'normals', because it might have only 1 element!
    ## UPDATE: Actually, this function is ill-defined if there are less than 2 normals.
    assert len( normals ) >= 2
    normals = asarray( normals )
    
    ## 3D normals
    assert len( normals.shape ) == 2
    assert normals.shape[1] == 3
    
    assert cos_theta == float(cos_theta)
    assert cos_theta > -1.
    
    ## Are the vectors normalized?
    assert vectors_are_normalized( normals )
    ## Normalize the vectors.
    #import helpers
    #normals = helpers.normalize_vectors( normals )
    
    for i in xrange(normals.shape[0]):
        for j in xrange(i+1,normals.shape[0]):
            ## If the angle between two vectors is greater than theta,
            ## then the dot product between them is less than cos( theta ).
            if dot( normals[i], normals[j] ) <= cos_theta: return False
    
    return True

def all_normals_are_less_than_cos_theta_from_n( normals, cos_theta, n ):
    '''
    Given a list of N-dimensional unit normals 'normals'
    the cosine of an angle 'cos_theta',
    and an additional normal 'n',
    returns True if the angle between every normal in 'normals' and 'n' is
    less than theta, and
    returns False otherwise.
    
    NOTE: Assumes that 'normals' and 'n' is normalized!
    
    tested
    '''
    
    ## We can't squeeze() 'normals', because it might have only 1 element!
    ## This function is ill-defined if there is less than 1 normal in 'normals'.
    assert len( normals ) >= 1
    normals = asarray( normals )
    n = asarray( n ).squeeze()
    
    ## 3D normals
    assert len( normals.shape ) == 2
    assert normals.shape[1] == 3
    assert n.shape == (3,)
    
    assert cos_theta == float(cos_theta)
    assert cos_theta > -1.
    
    ## Are the vectors normalized?
    assert vectors_are_normalized( normals )
    assert vector_is_normalized( n )
    ## Normalize the vectors.
    #import helpers
    #normals = helpers.normalize_vectors( normals )
    
    for i in xrange(normals.shape[0]):
        ## If the angle between two vectors is greater than theta,
        ## then the dot product between them is less than cos( theta ).
        if dot( normals[i], n ) <= cos_theta: return False
    
    return True

def normals_min_and_max_cos_angle_from_each_other( normals ):
    '''
    Given a list of N-dimensional unit normals 'normals',
    returns two values, the minimum and maximum cosine of the angle between
    any two normals in 'normals'.
    
    NOTE: Assumes that 'normals' is normalized!
    
    tested:
    >>> normals_min_and_max_cos_angle_from_each_other( [(1,0,0),(1,0,0),(0,0,1.)] )
    (0.0, 1.0)
    >>> normals_min_and_max_cos_angle_from_each_other( [(1,0,0),(1,0,0),(.5,0,0)] )
    AssertionError
    >>> normals_min_and_max_cos_angle_from_each_other( [(1,0,0),(1,0,0),(.5,0,sqrt(1-.5**2))] )
    (0.5, 1.0)
    >>> normals_min_and_max_cos_angle_from_each_other( [(1,0,0),(1,0,0),(-.5,0,sqrt(1-.5**2))] )
    (-0.5, 1.0)
    '''
    
    ## We can't squeeze() 'normals', because it might have only 1 element!
    ## UPDATE: Actually, this function is ill-defined if there are less than 2 normals.
    assert len( normals ) >= 2
    normals = asarray( normals )
    
    ## 3D normals
    assert len( normals.shape ) == 2
    assert normals.shape[1] == 3
    
    ## Are the vectors normalized?
    assert vectors_are_normalized( normals )
    ## Normalize the vectors.
    #import helpers
    #normals = helpers.normalize_vectors( normals )
    
    ## The cosine of an angle will always be in [-1,1],
    ## so these initial values are out-of-range and will
    ## be set by the first comparison.
    kTooLarge = 31337
    min_cos_angle = kTooLarge
    max_cos_angle = -kTooLarge
    
    for i in xrange(normals.shape[0]):
        for j in xrange(i+1,normals.shape[0]):
            cos_angle = dot( normals[i], normals[j] )
            min_cos_angle = min( cos_angle, min_cos_angle )
            max_cos_angle = max( cos_angle, max_cos_angle )
    
    assert min_cos_angle <= kTooLarge
    assert max_cos_angle >= -kTooLarge
    
    return min_cos_angle, max_cos_angle

def normals_min_and_max_cos_angle_from_n( normals, n ):
    '''
    Given a list of N-dimensional unit normals 'normals' and
    an additional normal 'n',
    returns two values, the minimum and maximum cosine of the angle between
    any two normals in 'normals'.
    
    NOTE: Assumes that 'normals' and 'n' is normalized!
    
    tested:
    >>> normals_min_and_max_cos_angle_from_n( [(1,0,0),(1,0,0),(0,0,1.)], (0,0,1.) )
    (0.0, 1.0)
    >>> normals_min_and_max_cos_angle_from_n( [(1,0,0),(1,0,0),(0,0,1.)], (0,0,.5) )
    AssertionError
    >>> normals_min_and_max_cos_angle_from_n( [(1,0,0),(1,0,0),(0,0,1.)], (1.,0,0) )
    (0.0, 1.0)
    >>> normals_min_and_max_cos_angle_from_n( [(1,0,0),(1,0,0),(0,0,1.)], (0,1.,0) )
    (0.0, 0.0)
    >>> normals_min_and_max_cos_angle_from_n( [(1,0,0),(1,0,0),(0,0,1.)], (.5,0,sqrt(1-.5**2)) )
    (0.5, 0.8660254037844386)
    '''
    
    ## We can't squeeze() 'normals', because it might have only 1 element!
    ## This function is ill-defined if there is less than 1 normal in 'normals'.
    assert len( normals ) >= 1
    normals = asarray( normals )
    n = asarray( n ).squeeze()
    
    ## 3D normals
    assert len( normals.shape ) == 2
    assert normals.shape[1] == 3
    assert n.shape == (3,)
    
    ## Are the vectors normalized?
    assert vectors_are_normalized( normals )
    assert vector_is_normalized( n )
    ## Normalize the vectors.
    #import helpers
    #normals = helpers.normalize_vectors( normals )
    
    ## The cosine of an angle will always be in [-1,1],
    ## so these initial values are out-of-range and will
    ## be set by the first comparison.
    kTooLarge = 31337
    min_cos_angle = kTooLarge
    max_cos_angle = -kTooLarge
    
    for i in xrange(normals.shape[0]):
        cos_angle = dot( normals[i], n )
        
        min_cos_angle = min( cos_angle, min_cos_angle )
        max_cos_angle = max( cos_angle, max_cos_angle )
    
    assert min_cos_angle <= kTooLarge
    assert max_cos_angle >= -kTooLarge
    
    return min_cos_angle, max_cos_angle

def sample_N_normals_more_than_cos_theta_from_each_other( N, cos_theta, normal_generator = sample_normal ):
    '''
    Given an integer 'N', and
    the cosine of an angle theta 'cos_theta',
    returns a sequence of 'N' 3D unit normals, where the angle between any two vectors
    is greater than theta.
    
    Takes an optional parameter 'normal_generator', a function that takes no arguments
    and returns a unit 3D vector.  The default uniformly samples 3D unit vectors.
    
    tested
    e.g. should be less than 0:
    asarray( [ dot( *sample_N_normals_more_than_cos_theta_from_each_other( 2, 0, sample_normal ) ) for i in xrange(1000) ] ).max()
    e.g. should be less than -.5:
    asarray( [ dot( *sample_N_normals_more_than_cos_theta_from_each_other( 2, -.5, sample_normal ) ) for i in xrange(1000) ] ).max()
    e.g. should be less than R:
    asarray( [ dot( *sample_N_normals_more_than_cos_theta_from_each_other( 2, R, sample_normal ) ) for i in xrange(1000) ] ).max()
    '''
    
    assert N == int(N)
    assert N > 0
    assert cos_theta == float(cos_theta)
    assert cos_theta > -1.
    
    while True:
        normals = asarray( [ normal_generator() for i in xrange(N) ] )
        if all_normals_are_more_than_cos_theta_from_each_other( normals, cos_theta ): return normals

def sample_N_points_more_than_d_from_each_other( N, d, point_generator = sample_3D_point ):
    '''
    Given an integer 'N', and
    a distance 'd',
    returns a sequence of 'N' k-dimensional points, where the distance between any two
    points is greater than 'd'.
    
    Takes an optional parameter 'point_generator', a function that takes no arguments
    and returns a point.  The default uniformly samples 3D points in [0,1]^3.
    '''
    
    assert N == int(N)
    assert N > 0
    assert d == float(d)
    assert d >= 0.
    
    while True:
        points = asarray( [ point_generator() for i in xrange(N) ] )
        if all_points_are_more_than_d_from_each_other( points, d ): return points

def sample_point_more_than_d_from_p( d, p, point_generator = sample_3D_point ):
    '''
    Given a distance 'd'
    and another point 'p',
    return a point whose distance to 'p' is greater than 'd'.
    
    Takes an optional parameter 'point_generator', a function that takes no arguments
    and returns a point.  The default uniformly samples 3D points in [0,1]^3.
    '''
    
    assert d == float(d)
    assert d >= 0.
    
    p = asarray( p ).squeeze()
    assert len( p.shape ) == 1
    
    while True:
        point = point_generator()
        if point_to_point_distance_squared( point, p ) > d*d: return point

def point_to_point_distance_squared( p1, p2 ):
    '''
    Given two N-dimensional points, returns the squared distance between them.
    '''
    
    p1 = asarray( p1 ).squeeze()
    p2 = asarray( p2 ).squeeze()
    assert p1.shape == p2.shape
    assert len( p1.shape ) == 1
    
    d = p1 - p2
    distsqr = dot( d, d )
    return distsqr

def point_to_point_distance( p1, p2 ):
    '''
    Given two N-dimensional points, returns the distance between them.
    '''
    
    return sqrt( point_to_point_distance_squared( p1, p2 ) )

def all_points_are_less_than_d_from_each_other( points, d ):
    '''
    Given a list of N-dimensional points 'points'
    and a scalar 'd',
    returns False if the distance between any two points is greater
    than or equal to 'd'.
    Returns True otherwise.
    '''
    
    ## We can't squeeze() 'points', because it might have only 1 element!
    ## UPDATE: Actually, this function is ill-defined if there are less than 2 points.
    assert len( points ) >= 2
    points = asarray( points ).squeeze()
    
    assert len( points.shape ) == 2
    
    assert d == float(d)
    assert d >= 0.
    
    for i in xrange(points.shape[0]):
        for j in xrange(i+1,points.shape[0]):
            if point_to_point_distance_squared( points[i], points[j] ) >= d*d: return False
    
    return True

def all_points_are_more_than_d_from_each_other( points, d ):
    '''
    Given a list of N-dimensional points 'points'
    and a scalar 'd',
    returns False if the distance between any two points is less
    than or equal to 'd'.
    Returns True otherwise.
    '''
    
    ## We can't squeeze() 'points', because it might have only 1 element!
    ## UPDATE: Actually, this function is ill-defined if there are less than 2 points.
    assert len( points ) >= 2
    points = asarray( points ).squeeze()
    
    assert len( points.shape ) == 2
    
    assert d == float(d)
    assert d >= 0.
    
    for i in xrange(points.shape[0]):
        for j in xrange(i+1,points.shape[0]):
            if point_to_point_distance_squared( points[i], points[j] ) <= d*d: return False
    
    return True

def all_points_are_less_than_d_from_p( points, d, p ):
    '''
    Given a list of N-dimensional points 'points',
    a scalar 'd',
    and another N-dimensional point 'p',
    returns False if the distance between any point
    in 'points' and 'p' is greater than or equal to 'd'.
    Returns True otherwise.
    '''
    
    assert len( points ) >= 1
    points = asarray( points )
    
    assert len( points.shape ) == 2
    
    p = asarray( p ).squeeze()
    assert len( p.shape ) == 1
    assert p.shape[0] == points.shape[1]
    
    assert d == float(d)
    assert d >= 0.
    
    for point in points:
        if point_to_point_distance_squared( point, p ) >= d*d: return False
    
    return True

def all_points_are_more_than_d_from_p( points, d, p ):
    '''
    Given a list of N-dimensional points 'points',
    a scalar 'd',
    and another N-dimensional point 'p',
    returns False if the distance between any point
    in 'points' and 'p' is less than or equal to 'd'.
    Returns True otherwise.
    '''
    
    assert len( points ) >= 1
    points = asarray( points )
    
    assert len( points.shape ) == 2
    
    p = asarray( p ).squeeze()
    assert len( p.shape ) == 1
    assert p.shape[0] == points.shape[1]
    
    assert d == float(d)
    assert d >= 0.
    
    for point in points:
        if point_to_point_distance_squared( point, p ) <= d*d: return False
    
    return True

def points_min_and_max_distance_from_each_other( points ):
    '''
    Given a list of N-dimensional points 'points'
    and a scalar 'd',
    returns two values, the minimum and maximum distance between
    any two points in 'points'.
    '''
    
    ## We can't squeeze() 'points', because it might have only 1 element!
    ## UPDATE: Actually, this function is ill-defined if there are less than 2 points.
    assert len( points ) >= 2
    points = asarray( points ).squeeze()
    
    assert len( points.shape ) == 2
    
    min_dist_sqr = point_to_point_distance_squared( points[0], points[1] )
    max_dist_sqr = min_dist_sqr
    
    for i in xrange(points.shape[0]):
        for j in xrange(i+1,points.shape[0]):
            d2 = point_to_point_distance_squared( points[i], points[j] )
            min_dist_sqr = min( d2, min_dist_sqr )
            max_dist_sqr = max( d2, max_dist_sqr )
    
    return sqrt( min_dist_sqr ), sqrt( max_dist_sqr )

def points_min_and_max_distance_from_p( points, p ):
    '''
    Given a list of N-dimensional points 'points',
    a scalar 'd',
    and another N-dimensional point 'p',
    returns two values, the minimum and maximum distance between
    any points in 'points' and 'p'.
    '''
    
    assert len( points ) >= 1
    points = asarray( points )
    
    assert len( points.shape ) == 2
    
    p = asarray( p ).squeeze()
    assert len( p.shape ) == 1
    assert p.shape[0] == points.shape[1]
    
    min_dist_sqr = point_to_point_distance_squared( points[0], p )
    max_dist_sqr = min_dist_sqr
    
    for point in points:
        d2 = point_to_point_distance_squared( point, p )
        min_dist_sqr = min( d2, min_dist_sqr )
        max_dist_sqr = max( d2, max_dist_sqr )
    
    return sqrt( min_dist_sqr ), sqrt( max_dist_sqr )


def random_interleave( *lists ):
    '''
    Given a list of lists [ l0, l1, l2, ..., ln[0] ]
    returns a flattened list [ l0[0], l1[0], l2[0], ..., ln[0], l0[1], l1[1], ... ].
    If the lengths of an input list is less than the maximum input list length,
    it will randomly skip its designated location such that it appears
    uniformly distributed in the output.
    
    NOTE: Each list is not itself shuffled.
    '''
    
    ## Handle a base case.
    if len( lists ) == 0: return []
    
    lists = [ list( l ) for l in lists ]
    
    total_length = sum([ len( l ) for l in lists ])
    
    result = []
    while True:
        remaining = max([ len( l ) for l in lists ])
        if remaining == 0: break
        
        for l in lists:
            if len( l ) == 0:
                continue
            elif len( l ) == remaining or random.uniform( 0, 1 ) < float( len( l ) ) / remaining:
                result.append( l[0] )
                del l[0]
    
    assert list(set([ len( l ) for l in lists ])) == [0]
    assert len( result ) == total_length
    return result

def all_rows_cols_inside_image_path( rows_cols, image_path ):
    '''
    Given a path to an image 'image_path',
    returns True if every row, column tuple in 'rows_cols'
    is within the image bounds.
    Returns False otherwise.
    '''
    
    import Image
    img = Image.open( image_path )
    return all( [ row >= 0 and row < img.size[1] and col >= 0 and col < img.size[0] for row, col in rows_cols ] )

def cluster( values, close_enough_together, disjoint = False ):
    '''
    Given a list of 'values' and
    a function taking a list of values and returning True if the values are
    close enough together and False otherwise,
    returns a set of clusters in the form of sets of indices into 'values'
    that are close enough together as determined by 'close_enough_together',
    where every index into 'values' appears at least once.
    
    If optional parameter 'disjoint' is True, then the sets in the
    result will be disjoint and every index into 'values' will appear
    exactly once.
    The clusters are culled in a greedy manner, not by solving the
    knapsack problem.
    One effect of the greedy algorithm is that the largest cluster
    is always returned.
    
    tested:
    cluster( [1,1,1,4,5,3,2,2,9,-1], lambda vals: len(set(vals)) <= 1 )
    [[0,1,2],[3],[4],[5],[6,7],[8],[9]]
    '''
    
    if len( values ) == 0: return []
    
    clusters = set([ frozenset([i]) for i, a in enumerate( values ) ])
    
    while True:
        new_clusters = set()
        for cluster in clusters:
            for i in xrange(len( values )):
                candidate = cluster.union(frozenset([i]))
                if candidate in clusters: continue
                if close_enough_together( [ values[i] for i in candidate ] ):
                    new_clusters.add( candidate )
        clusters = clusters.union( new_clusters )
        #print 'len( new_clusters ):', len( new_clusters )
        
        ## break if no progress.
        if len( new_clusters ) == 0: break
    
    if disjoint:
        ## sort clusters by size in descending order.
        clusters = list( clusters )
        clusters.sort( key = lambda x: len( x ) )
        clusters.reverse()
        
        seen = set()
        ## walk clusters from largest on down, taking the cluster only if we haven't yet 'seen' it.
        result = []
        for cluster in clusters:
            if len( seen.intersection( cluster ) ) == 0:
                seen.update( cluster )
                result.append( cluster )
        
        clusters = set( result )
        ## If sets are disjoint, the total length of all sets should be the
        ## length of the input values.
        assert sum( [ len( c ) for c in clusters ] ) == len( values )
        ## Every element in the input should be in exactly cluster.
        assert all( [ sum( [ i in c for c in clusters ] ) == 1 for i in xrange(len( values )) ] )
    
    ## Every element in the input should be in at least one cluster.
    assert all( [ sum( [ i in c for c in clusters ] ) >= 1 for i in xrange(len( values )) ] )
    
    return clusters

def n_choose_k( n, k ):
    '''
    Given positive integers 'n' and 'k',
    returns n-choose-k, the number of size-k
    subsets of a set of size n.
    
    This code comes from:
    http://stackoverflow.com/questions/2096573/counting-combinations-and-permutations-efficiently
    '''
    
    from itertools import izip
    return reduce(lambda x, y: x * y[0] / y[1], izip(xrange(n - k + 1, n+1), xrange(1, k+1)), 1)
