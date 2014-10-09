from numpy import *

def closest_distsqr_and_edge_index_and_t_on_edges_to_point( edges, pt ):
    '''
    Given a sequence of edges (pairs of points) 'edges' and n-dimensional point 'pt',
    returns the tuple (
        squared distance to the closest edge on 'line_strip',
        index of edge in 'edges' which contains the closest point,
        t along the edge such that the closest point is (1-t)*edges[index][0] + t*(edges[index][1])
        ).
    
    tested (see test_line_strip_and_loop_distances())
    '''
    
    assert len( edges ) > 0
    
    ## This function takes a sequence of points, so we have to pack 'pt' into a length-1 list.
    ## It also takes edges as pairs of points, not as a line loop, so we have to duplicate points.
    result = min_distanceSqr_edge_t_to_edges( [ pt ], edges )
    ## Unpack the result, since we passed a length-1 list containing 'pt'.
    ## Also, since min_distanceSqr_edge_t_to_edges supports arbitrary outer dimensions,
    ## squeeze() everything.
    return result[0][0].squeeze(), result[1][0].squeeze(), result[2][0].squeeze()

def closest_distsqr_and_edge_index_and_t_on_line_strip_to_point( line_strip, pt ):
    '''
    Given a line strip (sequence of n-dimensional points) 'linestrip' and n-dimensional point 'pt',
    returns the tuple (
        squared distance to the closest point on 'line_strip',
        index of point in 'line_strip' where the edge (index, index+1) contains the closest point ),
        t along the line strip such that the closest point is (1-t)*line_strip[index] + t*(line_strip[index+1])
        ).
    
    tested (see test_line_strip_and_loop_distances())
    '''
    
    assert len( line_strip ) >= 2
    
    return closest_distsqr_and_edge_index_and_t_on_edges_to_point(
        list( zip( line_strip[:-1], line_strip[1:] ) ),
        pt
        )

def closest_distsqr_and_edge_index_and_t_on_line_loop_to_point( line_loop, pt ):
    '''
    Same as closest_distsqr_and_point_and_edge_index_on_line_strip_to_point(), but
    takes a line loop (closed path) instead of a line strip (open path).
    
    NOTE: The index of the closing edge is len( line_loop )-1.
    
    tested (see test_line_strip_and_loop_distances())
    '''
    
    assert len( line_loop ) >= 3
    
    ## Append the first point to the end to turn the line loop into a line strip.
    return closest_distsqr_and_edge_index_and_t_on_line_strip_to_point(
        list( line_loop ) + [line_loop[0]], pt
        )

def test_line_strip_and_loop_distances():
    pt = (0,0)
    line_strip = [ ( 0, -.1 ), ( 1, -.1 ), ( 1, .9 ), ( 0, .9 ) ]
    print 'pt:', pt
    print 'line_strip:', line_strip
    print 'closest_distsqr_and_edge_index_and_t_on_line_strip_to_point():'
    print closest_distsqr_and_edge_index_and_t_on_line_strip_to_point( line_strip, pt )
    
    ## Interpreted as a line loop, 'line_strip' should pass directly through 'pt'.
    print 'closest_distsqr_and_edge_index_and_t_on_line_loop_to_point():'
    print closest_distsqr_and_edge_index_and_t_on_line_loop_to_point( line_strip, pt )

def distancesSqr_and_t_to_edges( pts, edges ):
    '''
    Input parameter 'pts' has dimensions #pts x (x,y,...).
    Input parameter 'edges' has dimensions ... x #edges x 2 endpoints x N coordinates (x,y,...).
    Returns an array of distances squared with dimensions #edges x #pts.
    '''
    
    pts = asfarray( pts ).T
    ## pts has dimensions N (x,y,...) x #pts
    edges = asfarray( edges )
    ## edges has dimensions ... x #edges x 2 endpoints x N coordinates (x,y,...)
    
    N = pts.shape[0]
    
    assert len( pts.shape ) == 2 and pts.shape[0] == N
    assert edges.shape[-2] == 2 and edges.shape[-1] == N
    #print 'pts.shape:', pts.shape
    #print 'edges.shape:', edges.shape
    
    
    ## get distance squared to each edge:
    ##   let p = black_pixel_pos, a = endpoint0, b = endpoint1, d = ( b-a ) / dot( b-a,b-a )
    ##   dot( p-a, d ) < 0 => dot( p-a, p-a )
    ##   dot( p-a, d ) > 1 => dot( p-b, p-b )
    ##   else              => dot( dot( p-a, d ) * (b-a) - p, same )
    p_a = pts[newaxis,...] - edges[...,:,0,:,newaxis]
    p_b = pts[newaxis,...] - edges[...,:,1,:,newaxis]
    ## p_a and p_b have dimensions ... x #edges x N coordinates (x,y,...) x #pts
    b_a = edges[...,:,1,:] - edges[...,:,0,:]
    ## b_a has dimensions ... x #edges x N coordinates (x,y,...)
    d = b_a / ( b_a**2 ).sum( -1 )[...,newaxis]
    ## d has same dimensions as b_a
    assert b_a.shape == d.shape
    cond = ( p_a * d[...,newaxis] ).sum( -2 )
    ## cond has dimensions ... x #edges x #pts
    assert cond.shape[-2:] == (edges.shape[-3], pts.shape[-1])
    
    ## clip cond so that values are between 0 and 1
    cond.clip( 0, 1, cond )
    
    #distancesSqr = empty( cond.shape, Real )
    ## distancesSqr has dimensions ... x #edges x #pts
    #assert distancesSqr.shape[-2:] == (edges.shape[-3], pts.shape[-1])
    
    # distancesSqr = p_a - cond[:,newaxis,:] * b_a[...,newaxis]
    # distancesSqr = ( distancesSqr**2 ).sum( 1 )
    # <=>
    distancesSqr = ( ( p_a - cond[:,newaxis,:] * b_a[...,newaxis] )**2 ).sum( -2 )
    
    #print 'distancesSqr:', distancesSqr
    #print 'distances:', sqrt( distancesSqr.min(0) )
    
    return distancesSqr, cond

def min_distanceSqr_edge_t_to_edges( pts, edges ):
    '''
    Input parameter 'pts' has dimensions #pts x 2 (x,y).
    Input parameter 'edges' has dimensions ... x #edges x 2 endpoints x 2 coordinates (x,y).
    Returns the tuple three things, each of which has length #pts and stores:
        (
        squared distance to the closest edge on 'line_strip',
        index of edge in 'edges' which contains the closest point,
        t along the edge such that the closest point is (1-t)*edges[index][0] + t*(edges[index][1])
        ).
    '''
    
    distancesSqr, cond = distancesSqr_and_t_to_edges( pts, edges )
    edge_index = distancesSqr.argmin( -2 )
    return distancesSqr[ edge_index ], edge_index, cond[ edge_index ]

def test_timing():
    import timeit, time
    
    N = 10
    dim = 3
    npts = 300
    nedges = 100
    
    #random.uniform( 0, 1 )
    
    setup2d_basic = \
'''
from numpy import asarray, array
from edge_distances import distancesSqr_and_t_to_edges
pts = [(.5, 0.), (.5,1.), (0.,.25), (1.,.25), (0.,0.), (1.,0.), (.3,.3)]
pts = asarray( pts )
edges = [[(.25, .25), (.75, .25)], [(.25, .25), (.25, .75)], [(.75, .25), (.75, .75)], [(.25, .75), (.75, .75)]]
'''
    
    setupNd = \
'''
from numpy import asarray, array
from edge_distances import distancesSqr_and_t_to_edges
import random
r = random.Random( 7 )
u = r.uniform
dim = %s
npts = %s
nedges = %s
pts = [ [ u( 0,1 ) for d in range(dim) ] for i in xrange( npts ) ]
pts = asarray( pts )
edges = [ [ [ u( 0,1 ) for d in range(dim) ], [ u( 0,1 ) for d in range(dim) ] ] for i in xrange( nedges ) ]
''' % (dim, npts, nedges)
    
    print 'repititions:', N
    print 'dimension:', dim
    print 'num points:', npts
    print 'num edges:', nedges
    
    #timeit.Timer( 'print abs( distancesSqr_to_edges( pts, edges ) - distancesSqr_to_edges_Nd( pts, edges ) ).sum(),', setup2d, time.clock ).timeit(100)
    #print 
    
    ## I have to interleave calls to the two functions I'm comparing, since whichever executes
    ## second has an advantage.  I verified this by calling identical functions.
    
    
    print 'distancesSqr_to_edges:',
    print timeit.Timer( 'distancesSqr_and_t_to_edges( pts, edges )', setupNd ).repeat(4,N)

def test1():
    pt = [ 0,0 ]
    edges = [
        ## distance should be 0 at t = 0
        [ (0,0), (1,0) ],
        [ (0,0), (0,1) ],
        ## distance should be 0 at t = 1
        [ (1,0), (0,0) ],
        [ (0,1), (0,0) ],
        ## distance should be 1 at t = 0
        [ (0,1), (1,1) ],
        ## distance should be 1 at t = 1
        [ (1,1), (0,1) ],
        ## distance should be .3 at t = .75
        [ (-.75,.3), (.25,.3) ],
        ## distance should be 0 at t = .5
        [ (-1,0), (1,0) ],
        ## distance should be 0 at t = .25
        [ (-.25,0), (.75,0) ],
        ]
    distSqrs, ts = distancesSqr_and_t_to_edges( [ pt ], edges )
    print 'distSqrs:'
    print distSqrs
    print 'ts:'
    print ts
    minDistSqr, min_edge_index, min_t = min_distanceSqr_edge_t_to_edges( [ pt ], edges )
    print 'min_distanceSqr_edge_t_to_edges():', minDistSqr, min_edge_index, min_t

def test2():
    pts = [(.5, 0.), (.5,1.), (0.,.25), (1.,.25), (0.,0.), (1.,0.), (.3,.3)]
    pts = asarray( pts )
    edges = [[(.25, .25), (.75, .25)], [(.25, .25), (.25, .75)], [(.75, .25), (.75, .75)], [(.25, .75), (.75, .75)]]
    distSqrs, ts = distancesSqr_and_t_to_edges( pts, edges )
    print 'distSqrs:'
    print distSqrs
    print 'ts:'
    print ts
    print 'min_distanceSqr_edge_t_to_edges():', min_distanceSqr_edge_t_to_edges( pts, edges )

def main():
    test1()
    #test_timing()

if __name__ == '__main__':
    main()
