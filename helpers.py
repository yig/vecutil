## This file is composed of some functions from ~/Work/ERATO/tracing/code/helpers.py
## UPDATE: With the addition of more functions.

from numpy import *

def color_spiral( i, N ):
    '''
    Takes parameters 'i' and 'N' where 'i' is a number in the range [0,'N'], where N > 0.
    Returns a 3-tuple ( red, green, blue ) of well-stratified colors, each value in the range [0,1].
    '''
    #from math import sin, sqrt
    #col = ( sin( 1 + i ), sin( 1 + 7 * i ), sin( 1 + 11 * i ) )
    #from random import random
    #col = ( random(), random(), random() )
    from colorsys import hsv_to_rgb
    t = float( i ) / N
    col = hsv_to_rgb(
        ( 2 * sqrt(t) * t ) % 1,
        sqrt( 1. - .8 * t ),
        1
        )
    
    col = tuple( map( abs, list( col ) ) )
    
    return col

def isiterable( a ):
    '''
    Returns True is 'a' is iterable, False otherwise.
    '''
    ## http://bytes.com/forum/thread514838.html
    try:
        it = iter(a)
        return True
    except TypeError:
        return False

def unique( in_fname ):
    import os
    
    orig_fname = in_fname.rstrip('/')
    
    count = 1
    fname = orig_fname
    
    while os.path.exists( fname ):
        fname = os.path.splitext( orig_fname )[0] + ' ' + str(count) + os.path.splitext( orig_fname )[1]
        count += 1
    
    return fname

kImageFormatsLossy = ['jpg','jpeg']
kImageFormatsLossless = ['png','tif','tiff','bmp','gif','psd','pdf','tga','pnm','ppm']
def seems_image( path ):
    import os
    return os.path.splitext( path )[1] in kImageFormatsLossy + kImageFormatsLossless

def all_image_paths_in_dir( dirpath ):
    return all_paths_in_dir( dirpath, filter = seems_image )

def all_paths_in_dir( dirpath, filter = None ):
    import os
    
    if filter is None:
        filter = lambda x: True
    
    image_paths = []
    ## http://stackoverflow.com/questions/120656/directory-listing-in-python
    for root, dirs, files in os.walk( dirpath ):
        image_paths.extend( [ os.path.join( root, file ) for file in files if filter( file ) ] )
    
    return image_paths

def friendly_Image_open_asarray( filename ):
    import Image
    img = Image.open( filename )
    ## Merge down to a greyscale (Luminance) image.
    ## UPDATE: Actually, don't do this.  We want to work on a color image!
    #img = img.convert( 'L' )
    ## Convert back and forth to numpy array with numpy.asarray( img ) and Image.fromarray( arr )
    ## Source: http://stackoverflow.com/questions/384759/pil-and-numpy
    ## One-liner: arr = asarray( img )
    ## One-liner: Image.fromarray( arr ).save( 'foo.png' )
    ## One-liner: Image.fromarray( arr ).show()
    ## One-liner (back and forth): Image.fromarray( asarray( img, dtype = uint8 ) ).save( 'foo.png' )
    ## One-liner (back and forth): Image.fromarray( asarray( img, dtype = uint8 ) ).show()
    #arr = asarray( img ) / 255.
    arr = asarray( img, dtype = uint8 )
    ## Ignore the alpha channel if there is one.
    assert len( arr.shape ) == 2 or ( len( arr.shape ) == 3 and arr.shape[2] in (3,4) )
    if len( arr.shape ) == 3: arr = arr[:,:,:3]
    
    return arr

def friendly_Image_from_float_array( arr ):
    ## Use Image.fromarray() if your values are in [0,255].
    assert arr.dtype in ( float32, float64, float )
    
    import Image
    return Image.fromarray( asarray( ( arr * 255 ).clip( 0, 255 ), dtype = uint8 ) )

def normal_map_to_color_Image( arr ):
    '''
    Given a normal map 'arr' as an N-by-M-by-3 array where every 3-dimensional
    vector has unit length,
    returns a color Image obtained by shifting the numbers to lie within the range [0,255].
    '''
    
    import Image
    assert len( arr.shape ) == 3 and arr.shape[2] == 3
    colornormals = .5*arr + .5
    return Image.fromarray( asarray( (255*colornormals).clip( 0, 255 ).round(), dtype = uint8 ) )

def normalize_image( img ):
    '''
    Given a 2D array, linearly maps the smallest value to 0 and the largest value to 1.
    '''
    
    ## Convert to a floating point array because of the division.
    if img.dtype.kind != 'f':
        img = asarray( img, dtype = float )
    
    assert len( img.shape ) == 2
    
    min_val = img.min()
    max_val = img.max()
    result = ( img - min_val ) / ( max_val - min_val )
    return result

def float_img2char_img( img ):
    '''
    Given a 2D array representing an image with values from 0 to 1, returns a uint8 image ranging from 0 to 255.
    '''
    
    ## Specifying dtype=uint8 is crucial, otherwise the png came out as garbage.
    result = asarray( ( 255.*img ).clip( 0, 255 ).round(), dtype = uint8 )
    return result

def normalize_to_char_img( img ):
    '''
    Given a 2D array representing an image, returns a uint8 array representing
    the image after linearly mapping the values to lie within the
    range [0,255].
    '''
    return float_img2char_img( normalize_image( img ) )

def gen_uuid():
    '''
    Returns a uuid.UUID object according to my tastes.
    '''
    
    ## WARNING: uuid.uuid4(), which allegedly generates a random uuid, generates the same IDs in sequence when using multiprocessing processes (which fork(), I think).
    ## See http://stackoverflow.com/questions/2759644/python-multiprocessing-doesnt-play-nicely-with-uuid-uuid4
    ## Also see my python bug report: http://bugs.python.org/issue8621
    import uuid, os
    try:
        return uuid.UUID( bytes = os.urandom( 16 ), version = 4 )
    except NotImplementedError:
        return uuid.uuid1()

def histogram( vals ):
    '''
    Given an iterable sequence of values 'vals', return a dictionary mapping a value to the number
    of times it appears in 'vals'.
    '''
    
    hist = {}
    for val in vals:
        hist.setdefault( val, 0 )
        hist[ val ] += 1
    return hist

def dictlist( pairs ):
    '''
    Given a sequence of ( key, value ) pairs 'pairs',
    return a dictionary mapping each key to a list of all its values.
    '''
    
    d = {}
    for key, val in pairs:
        d.setdefault( key, [] )
        d[ key ].append( val )
    return d

class Struct( object ): pass
