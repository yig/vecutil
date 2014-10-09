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
