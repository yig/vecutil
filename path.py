from numpy import *

def get_point_on_path( t, pts ):
    '''
    Returns point 't' along the given sequence of k-dimensional points 'pts',
    as measured by path lengths. Interpolates the given points linearly.
    
    NOTE: The value 't' is clamped to lie within the range [0,1]
    
    tested
    pts = [ (0,0), (2,0), (2,1), (2,1), (0,1) ]
    get_point_on_path( 0, pts )
    get_point_on_path( -1, pts )
    '''
    
    eps = 1e-7
    
    t = max( 0., min( 1., t ) )
    
    ## Numpy array
    pts = asfarray( pts )
    ## Get all segment lengths
    lengths = sqrt( ( ( pts[1:] - pts[:-1] )**2 ).sum( axis = 1 ) )
    
    ## Get the cumulative length along the path
    cumlengths = cumsum( lengths )
    ## Find where t along the path lands
    tlen = t * cumlengths[-1]
    lands = searchsorted( cumlengths, tlen )
    
    ## If we fall off the right edge, return the last point.
    if lands == len( cumlengths ):
        return pts[-1]
    
    ## If the t falls inside an epsilon-sized line segment, return an endpoint
    ## rather than dividing by the length.
    if lengths[ lands ] < eps:
        return pts[ lands ]
    
    ## Find the percent along the next segment.
    t_in_segment = float( tlen - ( cumlengths[ lands-1 ] if lands > 0 else 0. ) ) / lengths[ lands ]
    
    ## Linearly interpolate t_in_segment in the last segment
    result = pts[ lands ] + t_in_segment * ( pts[ lands + 1 ] - pts[ lands ] )
    
    return result

def compute_largest_minimum_distance( points1, points2 ):
    '''
    Given two sequences of N-dimensional points,
    computes the Hausdorff distance.
    The Hausdorff distance is the largest minimum distance between any two points.
    Returns the distance d, the index into points1, and the index into points2,
    such that
        d == distance( points1[ index into points1 ], points2[ index into points2 ] )
    is the largest minimum distance between any two points.
    
    tested
    '''
    
    assert len( points1 ) == len( points2 )
    points1 = asarray( points1 )
    points2 = asarray( points2 )
    assert points1.shape == points2.shape
    
    ## allDistSqrs[i][j] is the distance squared from points1[i] to points2[j].
    allDistSqrs = ( (points2[newaxis,...] - points1[:,newaxis,:])**2 ).sum(-1)
    ## Hausdorff distance is the longest shortest distance from either to either.
    dist2 = max( allDistSqrs.min(0).max(), allDistSqrs.min(1).max() )
    indices = where( allDistSqrs == dist2 )
    assert len( indices ) > 0
    
    points1_index = indices[0][0]
    points2_index = indices[1][0]
    dist = sqrt( dist2 )
    
    return dist, points1_index, points2_index

def all_distances( points1, points2 ):
    '''
    Given a sequence of n K-dimensional points 'points1' and
    a sequence of m K-dimensional points 'points2',
    returns an n-by-m matrix where the i,j entry is the distance-squared between
    the i-th point of 'points1' and the j-th point of 'points2'.
    
    To find all distances below some threshold, use:
        D2 = all_distances( points1, points2 )
        entries = where( D2 < threshold**2 )
        for i,j in zip( *entries ):
            ## Do something with D2[i,j]
    
    tested
    '''
    
    points1 = asarray( points1 )
    points2 = asarray( points2 )
    
    allDistSqrs = ( (points2[newaxis,...] - points1[:,newaxis,:])**2 ).sum(-1)
    assert allDistSqrs.shape == ( len( points1 ), len( points2 ) )
    return allDistSqrs

def test_all_distances():
    points1 = [ ( 0,0 ), ( 1,0 ), ( 2,0 ) ]
    points2 = [ ( 0,0 ), ( 0,1 ) ]
    
    threshold = .01
    
    D2 = all_distances( points1, points2 )
    entries = where( D2 < threshold**2 )
    for i,j in zip( *entries ):
        print i,j, D2[i,j]

def main():
    test_all_distances()

if __name__ == '__main__': main()
