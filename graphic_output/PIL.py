#
# The Python Imaging Library
# $Id$
#
# drawing interface operations
#
# History:
# 1996-04-13 fl   Created (experimental)
# 1996-08-07 fl   Filled polygons, ellipses.
# 1996-08-13 fl   Added text support
# 1998-06-28 fl   Handle I and F images
# 1998-12-29 fl   Added arc; use arc primitive to draw ellipses
# 1999-01-10 fl   Added shape stuff (experimental)
# 1999-02-06 fl   Added bitmap support
# 1999-02-11 fl   Changed all primitives to take options
# 1999-02-20 fl   Fixed backwards compatibility
# 2000-10-12 fl   Copy on write, when necessary
# 2001-02-18 fl   Use default ink for bitmap/text also in fill mode
# 2002-10-24 fl   Added support for CSS-style color strings
# 2002-12-10 fl   Added experimental support for RGBA-on-RGB drawing
# 2002-12-11 fl   Refactored low-level drawing API (work in progress)
# 2004-08-26 fl   Made Draw() a factory function, added getdraw() support
# 2004-09-04 fl   Added width support to line primitive
# 2004-09-10 fl   Added font mode handling
# 2006-06-19 fl   Added font bearing support (getmask2)
#
# Copyright (c) 1997-2006 by Secret Labs AB
# Copyright (c) 1996-2006 by Fredrik Lundh
#
# See the README file for information on usage and redistribution.
#

import numbers

from PIL import Image, ImageColor
from PIL._util import isStringType

try:
    import warnings
except ImportError:
    warnings = None


##
# A simple 2D drawing interface for PIL images.
# <p>
# Application code should use the <b>Draw</b> factory, instead of
# directly.

class ImageDraw:

    ##
    # Create a drawing instance.
    #
    # @param im The image to draw in.
    # @param mode Optional mode to use for color values.  For RGB
    #    images, this argument can be RGB or RGBA (to blend the
    #    drawing into the image).  For all other modes, this argument
    #    must be the same as the image mode.  If omitted, the mode
    #    defaults to the mode of the image.

    def __init__(self, im, mode=None):
        im.load()
        if im.readonly:
            im._copy()  # make it writeable
        blend = 0
        if mode is None:
            mode = im.mode
        if mode != im.mode:
            if mode == "RGBA" and im.mode == "RGB":
                blend = 1
            else:
                raise ValueError("mode mismatch")
        if mode == "P":
            self.palette = im.palette
        else:
            self.palette = None
        self.im = im.im
        self.draw = Image.core.draw(self.im, blend)
        self.mode = mode
        if mode in ("I", "F"):
            self.ink = self.draw.draw_ink(1, mode)
        else:
            self.ink = self.draw.draw_ink(-1, mode)
        if mode in ("1", "P", "I", "F"):
            # FIXME: fix Fill2 to properly support matte for I+F images
            self.fontmode = "1"
        else:
            self.fontmode = "L"  # aliasing is okay for other modes
        self.fill = 0
        self.font = None

    ##
    # Set the default pen color.

    def setink(self, ink):
        # compatibility
        if warnings:
            warnings.warn(
                "'setink' is deprecated; use keyword arguments instead",
                DeprecationWarning, stacklevel=2
                )
        if isStringType(ink):
            ink = ImageColor.getcolor(ink, self.mode)
        if self.palette and not isinstance(ink, numbers.Number):
            ink = self.palette.getcolor(ink)
        self.ink = self.draw.draw_ink(ink, self.mode)

    ##
    # Set the default background color.

    def setfill(self, onoff):
        # compatibility
        if warnings:
            warnings.warn(
                "'setfill' is deprecated; use keyword arguments instead",
                DeprecationWarning, stacklevel=2
                )
        self.fill = onoff

    ##
    # Set the default font.

    def setfont(self, font):
        # compatibility
        self.font = font

    ##
    # Get the current default font.

    def getfont(self):
        if not self.font:
            # FIXME: should add a font repository
            from PIL import ImageFont
            self.font = ImageFont.load_default()
        return self.font

    def _getink(self, ink, fill=None):
        if ink is None and fill is None:
            if self.fill:
                fill = self.ink
            else:
                ink = self.ink
        else:
            if ink is not None:
                if isStringType(ink):
                    ink = ImageColor.getcolor(ink, self.mode)
                if self.palette and not isinstance(ink, numbers.Number):
                    ink = self.palette.getcolor(ink)
                ink = self.draw.draw_ink(ink, self.mode)
            if fill is not None:
                if isStringType(fill):
                    fill = ImageColor.getcolor(fill, self.mode)
                if self.palette and not isinstance(fill, numbers.Number):
                    fill = self.palette.getcolor(fill)
                fill = self.draw.draw_ink(fill, self.mode)
        return ink, fill

    ##
    # Draw an arc.

    def arc(self, xy, start, end, fill=None):
        ink, fill = self._getink(fill)
        if ink is not None:
            self.draw.draw_arc(xy, start, end, ink)

    ##
    # Draw a bitmap.

    def bitmap(self, xy, bitmap, fill=None):
        bitmap.load()
        ink, fill = self._getink(fill)
        if ink is None:
            ink = fill
        if ink is not None:
            self.draw.draw_bitmap(xy, bitmap.im, ink)

    ##
    # Draw a chord.

    def chord(self, xy, start, end, fill=None, outline=None):
        ink, fill = self._getink(outline, fill)
        if fill is not None:
            self.draw.draw_chord(xy, start, end, fill, 1)
        if ink is not None:
            self.draw.draw_chord(xy, start, end, ink, 0)

    ##
    # Draw an ellipse.

    def ellipse(self, xy, fill=None, outline=None):
        ink, fill = self._getink(outline, fill)
        if fill is not None:
            self.draw.draw_ellipse(xy, fill, 1)
        if ink is not None:
            self.draw.draw_ellipse(xy, ink, 0)

    ##
    # Draw a line, or a connected sequence of line segments.

    def line(self, xy, fill=None, width=0):
        ink, fill = self._getink(fill)
        if ink is not None:
            self.draw.draw_lines(xy, ink, width)

    ##
    # (Experimental) Draw a shape.

    def shape(self, shape, fill=None, outline=None):
        # experimental
        shape.close()
        ink, fill = self._getink(outline, fill)
        if fill is not None:
            self.draw.draw_outline(shape, fill, 1)
        if ink is not None:
            self.draw.draw_outline(shape, ink, 0)

    ##
    # Draw a pieslice.

    def pieslice(self, xy, start, end, fill=None, outline=None):
        ink, fill = self._getink(outline, fill)
        if fill is not None:
            self.draw.draw_pieslice(xy, start, end, fill, 1)
        if ink is not None:
            self.draw.draw_pieslice(xy, start, end, ink, 0)

    ##
    # Draw one or more individual pixels.

    def point(self, xy, fill=None):
        ink, fill = self._getink(fill)
        if ink is not None:
            self.draw.draw_points(xy, ink)

    ##
    # Draw a polygon.

    def polygon(self, xy, fill=None, outline=None):
        ink, fill = self._getink(outline, fill)
        if fill is not None:
            self.draw.draw_polygon(xy, fill, 1)
        if ink is not None:
            self.draw.draw_polygon(xy, ink, 0)

    ##
    # Draw a rectangle.

    def rectangle(self, xy, fill=None, outline=None):
        ink, fill = self._getink(outline, fill)
        if fill is not None:
            self.draw.draw_rectangle(xy, fill, 1)
        if ink is not None:
            self.draw.draw_rectangle(xy, ink, 0)

    ##
    # Draw text.

    def text(self, xy, text, fill=None, font=None, anchor=None):
        ink, fill = self._getink(fill)
        if font is None:
            font = self.getfont()
        if ink is None:
            ink = fill
        if ink is not None:
            try:
                mask, offset = font.getmask2(text, self.fontmode)
                xy = xy[0] + offset[0], xy[1] + offset[1]
            except AttributeError:
                try:
                    mask = font.getmask(text, self.fontmode)
                except TypeError:
                    mask = font.getmask(text)
            self.draw.draw_bitmap(xy, mask, ink)

    ##
    # Get the size of a given string, in pixels.

    def textsize(self, text, font=None):
        if font is None:
            font = self.getfont()
        return font.getsize(text)


##
# A simple 2D drawing interface for PIL images.
#
# @param im The image to draw in.
# @param mode Optional mode to use for color values.  For RGB
#    images, this argument can be RGB or RGBA (to blend the
#    drawing into the image).  For all other modes, this argument
#    must be the same as the image mode.  If omitted, the mode
#    defaults to the mode of the image.

def Draw(im, mode=None):
    try:
        return im.getdraw(mode)
    except AttributeError:
        return ImageDraw(im, mode)

# experimental access to the outline API
try:
    Outline = Image.core.outline
except:
    Outline = None


##
# (Experimental) A more advanced 2D drawing interface for PIL images,
# based on the WCK interface.
#
# @param im The image to draw in.
# @param hints An optional list of hints.
# @return A (drawing context, drawing resource factory) tuple.

def getdraw(im=None, hints=None):
    # FIXME: this needs more work!
    # FIXME: come up with a better 'hints' scheme.
    handler = None
    if not hints or "nicest" in hints:
        try:
            from PIL import _imagingagg as handler
        except ImportError:
            pass
    if handler is None:
        from PIL import ImageDraw2 as handler
    if im:
        im = handler.Draw(im)
    return im, handler


##
# (experimental) Fills a bounded region with a given color.
#
# @param image Target image.
# @param xy Seed position (a 2-item coordinate tuple).
# @param value Fill color.
# @param border Optional border value.  If given, the region consists of
#     pixels with a color different from the border color.  If not given,
#     the region consists of pixels having the same color as the seed
#     pixel.

def floodfill(image, xy, value, border=None):
    "Fill bounded region."
    # based on an implementation by Eric S. Raymond
    pixel = image.load()
    x, y = xy
    try:
        background = pixel[x, y]
        if background == value:
            return  # seed point already has fill color
        pixel[x, y] = value
    except IndexError:
        return  # seed point outside image
    edge = [(x, y)]
    if border is None:
        while edge:
            newedge = []
            for (x, y) in edge:
                for (s, t) in ((x+1, y), (x-1, y), (x, y+1), (x, y-1)):
                    try:
                        p = pixel[s, t]
                    except IndexError:
                        pass
                    else:
                        if p == background:
                            pixel[s, t] = value
                            newedge.append((s, t))
            edge = newedge
    else:
        while edge:
            newedge = []
            for (x, y) in edge:
                for (s, t) in ((x+1, y), (x-1, y), (x, y+1), (x, y-1)):
                    try:
                        p = pixel[s, t]
                    except IndexError:
                        pass
                    else:
                        if p != value and p != border:
                            pixel[s, t] = value
                            newedge.append((s, t))
            edge = newedge




#
# The Python Imaging Library.
# $Id$
#
# the Image class wrapper
#
# partial release history:
# 1995-09-09 fl   Created
# 1996-03-11 fl   PIL release 0.0 (proof of concept)
# 1996-04-30 fl   PIL release 0.1b1
# 1999-07-28 fl   PIL release 1.0 final
# 2000-06-07 fl   PIL release 1.1
# 2000-10-20 fl   PIL release 1.1.1
# 2001-05-07 fl   PIL release 1.1.2
# 2002-03-15 fl   PIL release 1.1.3
# 2003-05-10 fl   PIL release 1.1.4
# 2005-03-28 fl   PIL release 1.1.5
# 2006-12-02 fl   PIL release 1.1.6
# 2009-11-15 fl   PIL release 1.1.7
#
# Copyright (c) 1997-2009 by Secret Labs AB.  All rights reserved.
# Copyright (c) 1995-2009 by Fredrik Lundh.
#
# See the README file for information on usage and redistribution.
#

from __future__ import print_function

from PIL import VERSION, PILLOW_VERSION, _plugins

import warnings


class DecompressionBombWarning(RuntimeWarning):
    pass


class _imaging_not_installed:
    # module placeholder
    def __getattr__(self, id):
        raise ImportError("The _imaging C module is not installed")


# Limit to around a quarter gigabyte for a 24 bit (3 bpp) image
MAX_IMAGE_PIXELS = int(1024 * 1024 * 1024 / 4 / 3)

try:
    # give Tk a chance to set up the environment, in case we're
    # using an _imaging module linked against libtcl/libtk (use
    # __import__ to hide this from naive packagers; we don't really
    # depend on Tk unless ImageTk is used, and that module already
    # imports Tkinter)
    __import__("FixTk")
except ImportError:
    pass

try:
    # If the _imaging C module is not present, you can still use
    # the "open" function to identify files, but you cannot load
    # them.  Note that other modules should not refer to _imaging
    # directly; import Image and use the Image.core variable instead.
    from PIL import _imaging as core
    if PILLOW_VERSION != getattr(core, 'PILLOW_VERSION', None):
        raise ImportError("The _imaging extension was built for another "
                          " version of Pillow or PIL")

except ImportError as v:
    core = _imaging_not_installed()
    # Explanations for ways that we know we might have an import error
    if str(v).startswith("Module use of python"):
        # The _imaging C module is present, but not compiled for
        # the right version (windows only).  Print a warning, if
        # possible.
        warnings.warn(
            "The _imaging extension was built for another version "
            "of Python.",
            RuntimeWarning
            )
    elif str(v).startswith("The _imaging extension"):
        warnings.warn(str(v), RuntimeWarning)
    elif "Symbol not found: _PyUnicodeUCS2_FromString" in str(v):
        warnings.warn(
            "The _imaging extension was built for Python with UCS2 support; "
            "recompile PIL or build Python --without-wide-unicode. ",
            RuntimeWarning
            )
    elif "Symbol not found: _PyUnicodeUCS4_FromString" in str(v):
        warnings.warn(
            "The _imaging extension was built for Python with UCS4 support; "
            "recompile PIL or build Python --with-wide-unicode. ",
            RuntimeWarning
            )
    # Fail here anyway. Don't let people run with a mostly broken Pillow.
    raise

try:
    import builtins
except ImportError:
    import __builtin__
    builtins = __builtin__

from PIL import ImageMode
from PIL._binary import i8
from PIL._util import isPath
from PIL._util import isStringType
from PIL._util import deferred_error

import os
import sys

# type stuff
import collections
import numbers

# works everywhere, win for pypy, not cpython
USE_CFFI_ACCESS = hasattr(sys, 'pypy_version_info')
try:
    import cffi
    HAS_CFFI = True
except:
    HAS_CFFI = False


def isImageType(t):
    """
    Checks if an object is an image object.

    .. warning::

       This function is for internal use only.

    :param t: object to check if it's an image
    :returns: True if the object is an image
    """
    return hasattr(t, "im")

#
# Debug level

DEBUG = 0

#
# Constants (also defined in _imagingmodule.c!)

NONE = 0

# transpose
FLIP_LEFT_RIGHT = 0
FLIP_TOP_BOTTOM = 1
ROTATE_90 = 2
ROTATE_180 = 3
ROTATE_270 = 4

# transforms
AFFINE = 0
EXTENT = 1
PERSPECTIVE = 2
QUAD = 3
MESH = 4

# resampling filters
NONE = 0
NEAREST = 0
ANTIALIAS = 1  # 3-lobed lanczos
LINEAR = BILINEAR = 2
CUBIC = BICUBIC = 3

# dithers
NONE = 0
NEAREST = 0
ORDERED = 1  # Not yet implemented
RASTERIZE = 2  # Not yet implemented
FLOYDSTEINBERG = 3  # default

# palettes/quantizers
WEB = 0
ADAPTIVE = 1

MEDIANCUT = 0
MAXCOVERAGE = 1
FASTOCTREE = 2

# categories
NORMAL = 0
SEQUENCE = 1
CONTAINER = 2

if hasattr(core, 'DEFAULT_STRATEGY'):
    DEFAULT_STRATEGY = core.DEFAULT_STRATEGY
    FILTERED = core.FILTERED
    HUFFMAN_ONLY = core.HUFFMAN_ONLY
    RLE = core.RLE
    FIXED = core.FIXED


# --------------------------------------------------------------------
# Registries

ID = []
OPEN = {}
MIME = {}
SAVE = {}
EXTENSION = {}

# --------------------------------------------------------------------
# Modes supported by this version

_MODEINFO = {
    # NOTE: this table will be removed in future versions.  use
    # getmode* functions or ImageMode descriptors instead.

    # official modes
    "1": ("L", "L", ("1",)),
    "L": ("L", "L", ("L",)),
    "I": ("L", "I", ("I",)),
    "F": ("L", "F", ("F",)),
    "P": ("RGB", "L", ("P",)),
    "RGB": ("RGB", "L", ("R", "G", "B")),
    "RGBX": ("RGB", "L", ("R", "G", "B", "X")),
    "RGBA": ("RGB", "L", ("R", "G", "B", "A")),
    "CMYK": ("RGB", "L", ("C", "M", "Y", "K")),
    "YCbCr": ("RGB", "L", ("Y", "Cb", "Cr")),
    "LAB": ("RGB", "L", ("L", "A", "B")),
    "HSV": ("RGB", "L", ("H", "S", "V")),

    # Experimental modes include I;16, I;16L, I;16B, RGBa, BGR;15, and
    # BGR;24.  Use these modes only if you know exactly what you're
    # doing...

}

if sys.byteorder == 'little':
    _ENDIAN = '<'
else:
    _ENDIAN = '>'

_MODE_CONV = {
    # official modes
    "1": ('|b1', None),  # broken
    "L": ('|u1', None),
    "I": (_ENDIAN + 'i4', None),
    "F": (_ENDIAN + 'f4', None),
    "P": ('|u1', None),
    "RGB": ('|u1', 3),
    "RGBX": ('|u1', 4),
    "RGBA": ('|u1', 4),
    "CMYK": ('|u1', 4),
    "YCbCr": ('|u1', 3),
    "LAB": ('|u1', 3),  # UNDONE - unsigned |u1i1i1
    # I;16 == I;16L, and I;32 == I;32L
    "I;16": ('<u2', None),
    "I;16B": ('>u2', None),
    "I;16L": ('<u2', None),
    "I;16S": ('<i2', None),
    "I;16BS": ('>i2', None),
    "I;16LS": ('<i2', None),
    "I;32": ('<u4', None),
    "I;32B": ('>u4', None),
    "I;32L": ('<u4', None),
    "I;32S": ('<i4', None),
    "I;32BS": ('>i4', None),
    "I;32LS": ('<i4', None),
}


def _conv_type_shape(im):
    shape = im.size[1], im.size[0]
    typ, extra = _MODE_CONV[im.mode]
    if extra is None:
        return shape, typ
    else:
        return shape+(extra,), typ


MODES = sorted(_MODEINFO.keys())

# raw modes that may be memory mapped.  NOTE: if you change this, you
# may have to modify the stride calculation in map.c too!
_MAPMODES = ("L", "P", "RGBX", "RGBA", "CMYK", "I;16", "I;16L", "I;16B")


def getmodebase(mode):
    """
    Gets the "base" mode for given mode.  This function returns "L" for
    images that contain grayscale data, and "RGB" for images that
    contain color data.

    :param mode: Input mode.
    :returns: "L" or "RGB".
    :exception KeyError: If the input mode was not a standard mode.
    """
    return ImageMode.getmode(mode).basemode


def getmodetype(mode):
    """
    Gets the storage type mode.  Given a mode, this function returns a
    single-layer mode suitable for storing individual bands.

    :param mode: Input mode.
    :returns: "L", "I", or "F".
    :exception KeyError: If the input mode was not a standard mode.
    """
    return ImageMode.getmode(mode).basetype


def getmodebandnames(mode):
    """
    Gets a list of individual band names.  Given a mode, this function returns
    a tuple containing the names of individual bands (use
    :py:method:`~PIL.Image.getmodetype` to get the mode used to store each
    individual band.

    :param mode: Input mode.
    :returns: A tuple containing band names.  The length of the tuple
        gives the number of bands in an image of the given mode.
    :exception KeyError: If the input mode was not a standard mode.
    """
    return ImageMode.getmode(mode).bands


def getmodebands(mode):
    """
    Gets the number of individual bands for this mode.

    :param mode: Input mode.
    :returns: The number of bands in this mode.
    :exception KeyError: If the input mode was not a standard mode.
    """
    return len(ImageMode.getmode(mode).bands)

# --------------------------------------------------------------------
# Helpers

_initialized = 0


def preinit():
    "Explicitly load standard file format drivers."

    global _initialized
    if _initialized >= 1:
        return

    try:
        from PIL import BmpImagePlugin
    except ImportError:
        pass
    try:
        from PIL import GifImagePlugin
    except ImportError:
        pass
    try:
        from PIL import JpegImagePlugin
    except ImportError:
        pass
    try:
        from PIL import PpmImagePlugin
    except ImportError:
        pass
    try:
        from PIL import PngImagePlugin
    except ImportError:
        pass
#   try:
#       import TiffImagePlugin
#   except ImportError:
#       pass

    _initialized = 1


def init():
    """
    Explicitly initializes the Python Imaging Library. This function
    loads all available file format drivers.
    """

    global _initialized
    if _initialized >= 2:
        return 0

    for plugin in _plugins:
        try:
            if DEBUG:
                print ("Importing %s" % plugin)
            __import__("PIL.%s" % plugin, globals(), locals(), [])
        except ImportError:
            if DEBUG:
                print("Image: failed to import", end=' ')
                print(plugin, ":", sys.exc_info()[1])

    if OPEN or SAVE:
        _initialized = 2
        return 1


# --------------------------------------------------------------------
# Codec factories (used by tobytes/frombytes and ImageFile.load)

def _getdecoder(mode, decoder_name, args, extra=()):

    # tweak arguments
    if args is None:
        args = ()
    elif not isinstance(args, tuple):
        args = (args,)

    try:
        # get decoder
        decoder = getattr(core, decoder_name + "_decoder")
        # print(decoder, mode, args + extra)
        return decoder(mode, *args + extra)
    except AttributeError:
        raise IOError("decoder %s not available" % decoder_name)


def _getencoder(mode, encoder_name, args, extra=()):

    # tweak arguments
    if args is None:
        args = ()
    elif not isinstance(args, tuple):
        args = (args,)

    try:
        # get encoder
        encoder = getattr(core, encoder_name + "_encoder")
        # print(encoder, mode, args + extra)
        return encoder(mode, *args + extra)
    except AttributeError:
        raise IOError("encoder %s not available" % encoder_name)


# --------------------------------------------------------------------
# Simple expression analyzer

def coerce_e(value):
    return value if isinstance(value, _E) else _E(value)


class _E:
    def __init__(self, data):
        self.data = data

    def __add__(self, other):
        return _E((self.data, "__add__", coerce_e(other).data))

    def __mul__(self, other):
        return _E((self.data, "__mul__", coerce_e(other).data))


def _getscaleoffset(expr):
    stub = ["stub"]
    data = expr(_E(stub)).data
    try:
        (a, b, c) = data  # simplified syntax
        if (a is stub and b == "__mul__" and isinstance(c, numbers.Number)):
            return c, 0.0
        if a is stub and b == "__add__" and isinstance(c, numbers.Number):
            return 1.0, c
    except TypeError:
        pass
    try:
        ((a, b, c), d, e) = data  # full syntax
        if (a is stub and b == "__mul__" and isinstance(c, numbers.Number) and
                d == "__add__" and isinstance(e, numbers.Number)):
            return c, e
    except TypeError:
        pass
    raise ValueError("illegal expression")


# --------------------------------------------------------------------
# Implementation wrapper

class Image:
    """
    This class represents an image object.  To create
    :py:class:`~PIL.Image.Image` objects, use the appropriate factory
    functions.  There's hardly ever any reason to call the Image constructor
    directly.

    * :py:func:`~PIL.Image.open`
    * :py:func:`~PIL.Image.new`
    * :py:func:`~PIL.Image.frombytes`
    """
    format = None
    format_description = None

    def __init__(self):
        # FIXME: take "new" parameters / other image?
        # FIXME: turn mode and size into delegating properties?
        self.im = None
        self.mode = ""
        self.size = (0, 0)
        self.palette = None
        self.info = {}
        self.category = NORMAL
        self.readonly = 0
        self.pyaccess = None

    def _new(self, im):
        new = Image()
        new.im = im
        new.mode = im.mode
        new.size = im.size
        new.palette = self.palette
        if im.mode == "P" and not new.palette:
            from PIL import ImagePalette
            new.palette = ImagePalette.ImagePalette()
        try:
            new.info = self.info.copy()
        except AttributeError:
            # fallback (pre-1.5.2)
            new.info = {}
            for k, v in self.info:
                new.info[k] = v
        return new

    _makeself = _new  # compatibility

    # Context Manager Support
    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        """
        Closes the file pointer, if possible.

        This operation will destroy the image core and release its memory.
        The image data will be unusable afterward.

        This function is only required to close images that have not
        had their file read and closed by the
        :py:meth:`~PIL.Image.Image.load` method.
        """
        try:
            self.fp.close()
        except Exception as msg:
            if DEBUG:
                print ("Error closing: %s" % msg)

        # Instead of simply setting to None, we're setting up a
        # deferred error that will better explain that the core image
        # object is gone.
        self.im = deferred_error(ValueError("Operation on closed image"))

    def _copy(self):
        self.load()
        self.im = self.im.copy()
        self.pyaccess = None
        self.readonly = 0

    def _dump(self, file=None, format=None):
        import tempfile
        suffix = ''
        if format:
            suffix = '.'+format
        if not file:
            f, file = tempfile.mkstemp(suffix)
            os.close(f)

        self.load()
        if not format or format == "PPM":
            self.im.save_ppm(file)
        else:
            if not file.endswith(format):
                file = file + "." + format
            self.save(file, format)
        return file

    def __eq__(self, other):
        if self.__class__.__name__ != other.__class__.__name__:
            return False
        a = (self.mode == other.mode)
        b = (self.size == other.size)
        c = (self.getpalette() == other.getpalette())
        d = (self.info == other.info)
        e = (self.category == other.category)
        f = (self.readonly == other.readonly)
        g = (self.tobytes() == other.tobytes())
        return a and b and c and d and e and f and g

    def __ne__(self, other):
        eq = (self == other)
        return not eq

    def __repr__(self):
        return "<%s.%s image mode=%s size=%dx%d at 0x%X>" % (
            self.__class__.__module__, self.__class__.__name__,
            self.mode, self.size[0], self.size[1],
            id(self)
            )

    def __getattr__(self, name):
        if name == "__array_interface__":
            # numpy array interface support
            new = {}
            shape, typestr = _conv_type_shape(self)
            new['shape'] = shape
            new['typestr'] = typestr
            new['data'] = self.tobytes()
            return new
        raise AttributeError(name)

    def __getstate__(self):
        return [
            self.info,
            self.mode,
            self.size,
            self.getpalette(),
            self.tobytes()]

    def __setstate__(self, state):
        Image.__init__(self)
        self.tile = []
        info, mode, size, palette, data = state
        self.info = info
        self.mode = mode
        self.size = size
        self.im = core.new(mode, size)
        if mode in ("L", "P"):
            self.putpalette(palette)
        self.frombytes(data)

    def tobytes(self, encoder_name="raw", *args):
        """
        Return image as a bytes object

        :param encoder_name: What encoder to use.  The default is to
                             use the standard "raw" encoder.
        :param args: Extra arguments to the encoder.
        :rtype: A bytes object.
        """

        # may pass tuple instead of argument list
        if len(args) == 1 and isinstance(args[0], tuple):
            args = args[0]

        if encoder_name == "raw" and args == ():
            args = self.mode

        self.load()

        # unpack data
        e = _getencoder(self.mode, encoder_name, args)
        e.setimage(self.im)

        bufsize = max(65536, self.size[0] * 4)  # see RawEncode.c

        data = []
        while True:
            l, s, d = e.encode(bufsize)
            data.append(d)
            if s:
                break
        if s < 0:
            raise RuntimeError("encoder error %d in tobytes" % s)

        return b"".join(data)

    # Declare tostring as alias to tobytes
    def tostring(self, *args, **kw):
        """Deprecated alias to tobytes.

        .. deprecated:: 2.0
        """
        warnings.warn(
            'tostring() is deprecated. Please call tobytes() instead.',
            DeprecationWarning,
            stacklevel=2,
        )
        return self.tobytes(*args, **kw)

    def tobitmap(self, name="image"):
        """
        Returns the image converted to an X11 bitmap.

        .. note:: This method only works for mode "1" images.

        :param name: The name prefix to use for the bitmap variables.
        :returns: A string containing an X11 bitmap.
        :raises ValueError: If the mode is not "1"
        """

        self.load()
        if self.mode != "1":
            raise ValueError("not a bitmap")
        data = self.tobytes("xbm")
        return b"".join([
            ("#define %s_width %d\n" % (name, self.size[0])).encode('ascii'),
            ("#define %s_height %d\n" % (name, self.size[1])).encode('ascii'),
            ("static char %s_bits[] = {\n" % name).encode('ascii'), data, b"};"
            ])

    def frombytes(self, data, decoder_name="raw", *args):
        """
        Loads this image with pixel data from a bytes object.

        This method is similar to the :py:func:`~PIL.Image.frombytes` function,
        but loads data into this image instead of creating a new image object.
        """

        # may pass tuple instead of argument list
        if len(args) == 1 and isinstance(args[0], tuple):
            args = args[0]

        # default format
        if decoder_name == "raw" and args == ():
            args = self.mode

        # unpack data
        d = _getdecoder(self.mode, decoder_name, args)
        d.setimage(self.im)
        s = d.decode(data)

        if s[0] >= 0:
            raise ValueError("not enough image data")
        if s[1] != 0:
            raise ValueError("cannot decode image data")

    def fromstring(self, *args, **kw):
        """Deprecated alias to frombytes.

        .. deprecated:: 2.0
        """
        warnings.warn(
            'fromstring() is deprecated. Please call frombytes() instead.',
            DeprecationWarning)
        return self.frombytes(*args, **kw)

    def load(self):
        """
        Allocates storage for the image and loads the pixel data.  In
        normal cases, you don't need to call this method, since the
        Image class automatically loads an opened image when it is
        accessed for the first time. This method will close the file
        associated with the image.

        :returns: An image access object.
        """
        if self.im and self.palette and self.palette.dirty:
            # realize palette
            self.im.putpalette(*self.palette.getdata())
            self.palette.dirty = 0
            self.palette.mode = "RGB"
            self.palette.rawmode = None
            if "transparency" in self.info:
                if isinstance(self.info["transparency"], int):
                    self.im.putpalettealpha(self.info["transparency"], 0)
                else:
                    self.im.putpalettealphas(self.info["transparency"])
                self.palette.mode = "RGBA"

        if self.im:
            if HAS_CFFI and USE_CFFI_ACCESS:
                if self.pyaccess:
                    return self.pyaccess
                from PIL import PyAccess
                self.pyaccess = PyAccess.new(self, self.readonly)
                if self.pyaccess:
                    return self.pyaccess
            return self.im.pixel_access(self.readonly)

    def verify(self):
        """
        Verifies the contents of a file. For data read from a file, this
        method attempts to determine if the file is broken, without
        actually decoding the image data.  If this method finds any
        problems, it raises suitable exceptions.  If you need to load
        the image after using this method, you must reopen the image
        file.
        """
        pass

    def convert(self, mode=None, matrix=None, dither=None,
                palette=WEB, colors=256):
        """
        Returns a converted copy of this image. For the "P" mode, this
        method translates pixels through the palette.  If mode is
        omitted, a mode is chosen so that all information in the image
        and the palette can be represented without a palette.

        The current version supports all possible conversions between
        "L", "RGB" and "CMYK." The **matrix** argument only supports "L"
        and "RGB".

        When translating a color image to black and white (mode "L"),
        the library uses the ITU-R 601-2 luma transform::

            L = R * 299/1000 + G * 587/1000 + B * 114/1000

        The default method of converting a greyscale ("L") or "RGB"
        image into a bilevel (mode "1") image uses Floyd-Steinberg
        dither to approximate the original image luminosity levels. If
        dither is NONE, all non-zero values are set to 255 (white). To
        use other thresholds, use the :py:meth:`~PIL.Image.Image.point`
        method.

        :param mode: The requested mode.
        :param matrix: An optional conversion matrix.  If given, this
           should be 4- or 16-tuple containing floating point values.
        :param dither: Dithering method, used when converting from
           mode "RGB" to "P" or from "RGB" or "L" to "1".
           Available methods are NONE or FLOYDSTEINBERG (default).
        :param palette: Palette to use when converting from mode "RGB"
           to "P".  Available palettes are WEB or ADAPTIVE.
        :param colors: Number of colors to use for the ADAPTIVE palette.
           Defaults to 256.
        :rtype: :py:class:`~PIL.Image.Image`
        :returns: An :py:class:`~PIL.Image.Image` object.
        """

        if not mode:
            # determine default mode
            if self.mode == "P":
                self.load()
                if self.palette:
                    mode = self.palette.mode
                else:
                    mode = "RGB"
            else:
                return self.copy()

        self.load()

        if matrix:
            # matrix conversion
            if mode not in ("L", "RGB"):
                raise ValueError("illegal conversion")
            im = self.im.convert_matrix(mode, matrix)
            return self._new(im)

        if mode == "P" and self.mode == "RGBA":
            return self.quantize(colors)

        trns = None
        delete_trns = False
        # transparency handling
        if "transparency" in self.info and \
                self.info['transparency'] is not None:
            if self.mode in ('L', 'RGB') and mode == 'RGBA':
                # Use transparent conversion to promote from transparent
                # color to an alpha channel.
                return self._new(self.im.convert_transparent(
                    mode, self.info['transparency']))
            elif self.mode in ('L', 'RGB', 'P') and mode in ('L', 'RGB', 'P'):
                t = self.info['transparency']
                if isinstance(t, bytes):
                    # Dragons. This can't be represented by a single color
                    warnings.warn('Palette images with Transparency  ' +
                                  ' expressed in bytes should be converted ' +
                                  'to RGBA images')
                    delete_trns = True
                else:
                    # get the new transparency color.
                    # use existing conversions
                    trns_im = Image()._new(core.new(self.mode, (1, 1)))
                    if self.mode == 'P':
                        trns_im.putpalette(self.palette)
                    trns_im.putpixel((0, 0), t)

                    if mode in ('L', 'RGB'):
                        trns_im = trns_im.convert(mode)
                    else:
                        # can't just retrieve the palette number, got to do it
                        # after quantization.
                        trns_im = trns_im.convert('RGB')
                    trns = trns_im.getpixel((0, 0))

            elif self.mode == 'P' and mode == 'RGBA':
                t = self.info['transparency']
                delete_trns = True
                
                if isinstance(t, bytes):
                    self.im.putpalettealphas(t)
                elif isinstance(t, int):
                    self.im.putpalettealpha(t,0)
                else:
                    raise ValueError("Transparency for P mode should" +
                                     " be bytes or int")


        if mode == "P" and palette == ADAPTIVE:
            im = self.im.quantize(colors)
            new = self._new(im)
            from PIL import ImagePalette
            new.palette = ImagePalette.raw("RGB", new.im.getpalette("RGB"))
            if delete_trns:
                # This could possibly happen if we requantize to fewer colors.
                # The transparency would be totally off in that case.
                del(new.info['transparency'])
            if trns is not None:
                try:
                    new.info['transparency'] = new.palette.getcolor(trns)
                except:
                    # if we can't make a transparent color, don't leave the old
                    # transparency hanging around to mess us up.
                    del(new.info['transparency'])
                    warnings.warn("Couldn't allocate palette entry " +
                                  "for transparency")
            return new

        # colorspace conversion
        if dither is None:
            dither = FLOYDSTEINBERG

        try:
            im = self.im.convert(mode, dither)
        except ValueError:
            try:
                # normalize source image and try again
                im = self.im.convert(getmodebase(self.mode))
                im = im.convert(mode, dither)
            except KeyError:
                raise ValueError("illegal conversion")

        new_im = self._new(im)
        if delete_trns:
            # crash fail if we leave a bytes transparency in an rgb/l mode.
            del(new_im.info['transparency'])
        if trns is not None:
            if new_im.mode == 'P':
                try:
                    new_im.info['transparency'] = new_im.palette.getcolor(trns)
                except:
                    del(new_im.info['transparency'])
                    warnings.warn("Couldn't allocate palette entry " +
                                  "for transparency")
            else:
                new_im.info['transparency'] = trns
        return new_im

    def quantize(self, colors=256, method=None, kmeans=0, palette=None):

        # methods:
        #    0 = median cut
        #    1 = maximum coverage
        #    2 = fast octree

        # NOTE: this functionality will be moved to the extended
        # quantizer interface in a later version of PIL.

        self.load()

        if method is None:
            # defaults:
            method = 0
            if self.mode == 'RGBA':
                method = 2

        if self.mode == 'RGBA' and method != 2:
            # Caller specified an invalid mode.
            raise ValueError('Fast Octree (method == 2) is the ' +
                             ' only valid method for quantizing RGBA images')

        if palette:
            # use palette from reference image
            palette.load()
            if palette.mode != "P":
                raise ValueError("bad mode for palette image")
            if self.mode != "RGB" and self.mode != "L":
                raise ValueError(
                    "only RGB or L mode images can be quantized to a palette"
                    )
            im = self.im.convert("P", 1, palette.im)
            return self._makeself(im)

        im = self.im.quantize(colors, method, kmeans)
        return self._new(im)

    def copy(self):
        """
        Copies this image. Use this method if you wish to paste things
        into an image, but still retain the original.

        :rtype: :py:class:`~PIL.Image.Image`
        :returns: An :py:class:`~PIL.Image.Image` object.
        """
        self.load()
        im = self.im.copy()
        return self._new(im)

    def crop(self, box=None):
        """
        Returns a rectangular region from this image. The box is a
        4-tuple defining the left, upper, right, and lower pixel
        coordinate.

        This is a lazy operation.  Changes to the source image may or
        may not be reflected in the cropped image.  To break the
        connection, call the :py:meth:`~PIL.Image.Image.load` method on
        the cropped copy.

        :param box: The crop rectangle, as a (left, upper, right, lower)-tuple.
        :rtype: :py:class:`~PIL.Image.Image`
        :returns: An :py:class:`~PIL.Image.Image` object.
        """

        self.load()
        if box is None:
            return self.copy()

        # lazy operation
        return _ImageCrop(self, box)

    def draft(self, mode, size):
        """
        NYI

        Configures the image file loader so it returns a version of the
        image that as closely as possible matches the given mode and
        size.  For example, you can use this method to convert a color
        JPEG to greyscale while loading it, or to extract a 128x192
        version from a PCD file.

        Note that this method modifies the :py:class:`~PIL.Image.Image` object
        in place.  If the image has already been loaded, this method has no
        effect.

        :param mode: The requested mode.
        :param size: The requested size.
        """
        pass

    def _expand(self, xmargin, ymargin=None):
        if ymargin is None:
            ymargin = xmargin
        self.load()
        return self._new(self.im.expand(xmargin, ymargin, 0))

    def filter(self, filter):
        """
        Filters this image using the given filter.  For a list of
        available filters, see the :py:mod:`~PIL.ImageFilter` module.

        :param filter: Filter kernel.
        :returns: An :py:class:`~PIL.Image.Image` object.  """

        self.load()

        if isinstance(filter, collections.Callable):
            filter = filter()
        if not hasattr(filter, "filter"):
            raise TypeError("filter argument should be ImageFilter.Filter " +
                            "instance or class")

        if self.im.bands == 1:
            return self._new(filter.filter(self.im))
        # fix to handle multiband images since _imaging doesn't
        ims = []
        for c in range(self.im.bands):
            ims.append(self._new(filter.filter(self.im.getband(c))))
        return merge(self.mode, ims)

    def getbands(self):
        """
        Returns a tuple containing the name of each band in this image.
        For example, **getbands** on an RGB image returns ("R", "G", "B").

        :returns: A tuple containing band names.
        :rtype: tuple
        """
        return ImageMode.getmode(self.mode).bands

    def getbbox(self):
        """
        Calculates the bounding box of the non-zero regions in the
        image.

        :returns: The bounding box is returned as a 4-tuple defining the
           left, upper, right, and lower pixel coordinate. If the image
           is completely empty, this method returns None.

        """

        self.load()
        return self.im.getbbox()

    def getcolors(self, maxcolors=256):
        """
        Returns a list of colors used in this image.

        :param maxcolors: Maximum number of colors.  If this number is
           exceeded, this method returns None.  The default limit is
           256 colors.
        :returns: An unsorted list of (count, pixel) values.
        """

        self.load()
        if self.mode in ("1", "L", "P"):
            h = self.im.histogram()
            out = []
            for i in range(256):
                if h[i]:
                    out.append((h[i], i))
            if len(out) > maxcolors:
                return None
            return out
        return self.im.getcolors(maxcolors)

    def getdata(self, band=None):
        """
        Returns the contents of this image as a sequence object
        containing pixel values.  The sequence object is flattened, so
        that values for line one follow directly after the values of
        line zero, and so on.

        Note that the sequence object returned by this method is an
        internal PIL data type, which only supports certain sequence
        operations.  To convert it to an ordinary sequence (e.g. for
        printing), use **list(im.getdata())**.

        :param band: What band to return.  The default is to return
           all bands.  To return a single band, pass in the index
           value (e.g. 0 to get the "R" band from an "RGB" image).
        :returns: A sequence-like object.
        """

        self.load()
        if band is not None:
            return self.im.getband(band)
        return self.im  # could be abused

    def getextrema(self):
        """
        Gets the the minimum and maximum pixel values for each band in
        the image.

        :returns: For a single-band image, a 2-tuple containing the
           minimum and maximum pixel value.  For a multi-band image,
           a tuple containing one 2-tuple for each band.
        """

        self.load()
        if self.im.bands > 1:
            extrema = []
            for i in range(self.im.bands):
                extrema.append(self.im.getband(i).getextrema())
            return tuple(extrema)
        return self.im.getextrema()

    def getim(self):
        """
        Returns a capsule that points to the internal image memory.

        :returns: A capsule object.
        """

        self.load()
        return self.im.ptr

    def getpalette(self):
        """
        Returns the image palette as a list.

        :returns: A list of color values [r, g, b, ...], or None if the
           image has no palette.
        """

        self.load()
        try:
            if bytes is str:
                return [i8(c) for c in self.im.getpalette()]
            else:
                return list(self.im.getpalette())
        except ValueError:
            return None  # no palette

    def getpixel(self, xy):
        """
        Returns the pixel value at a given position.

        :param xy: The coordinate, given as (x, y).
        :returns: The pixel value.  If the image is a multi-layer image,
           this method returns a tuple.
        """

        self.load()
        if self.pyaccess:
            return self.pyaccess.getpixel(xy)
        return self.im.getpixel(xy)

    def getprojection(self):
        """
        Get projection to x and y axes

        :returns: Two sequences, indicating where there are non-zero
            pixels along the X-axis and the Y-axis, respectively.
        """

        self.load()
        x, y = self.im.getprojection()
        return [i8(c) for c in x], [i8(c) for c in y]

    def histogram(self, mask=None, extrema=None):
        """
        Returns a histogram for the image. The histogram is returned as
        a list of pixel counts, one for each pixel value in the source
        image. If the image has more than one band, the histograms for
        all bands are concatenated (for example, the histogram for an
        "RGB" image contains 768 values).

        A bilevel image (mode "1") is treated as a greyscale ("L") image
        by this method.

        If a mask is provided, the method returns a histogram for those
        parts of the image where the mask image is non-zero. The mask
        image must have the same size as the image, and be either a
        bi-level image (mode "1") or a greyscale image ("L").

        :param mask: An optional mask.
        :returns: A list containing pixel counts.
        """
        self.load()
        if mask:
            mask.load()
            return self.im.histogram((0, 0), mask.im)
        if self.mode in ("I", "F"):
            if extrema is None:
                extrema = self.getextrema()
            return self.im.histogram(extrema)
        return self.im.histogram()

    def offset(self, xoffset, yoffset=None):
        """
        .. deprecated:: 2.0

        .. note:: New code should use :py:func:`PIL.ImageChops.offset`.

        Returns a copy of the image where the data has been offset by the given
        distances. Data wraps around the edges. If **yoffset** is omitted, it
        is assumed to be equal to **xoffset**.

        :param xoffset: The horizontal distance.
        :param yoffset: The vertical distance.  If omitted, both
           distances are set to the same value.
        :returns: An :py:class:`~PIL.Image.Image` object.
        """
        if warnings:
            warnings.warn(
                "'offset' is deprecated; use 'ImageChops.offset' instead",
                DeprecationWarning, stacklevel=2
                )
        from PIL import ImageChops
        return ImageChops.offset(self, xoffset, yoffset)

    def paste(self, im, box=None, mask=None):
        """
        Pastes another image into this image. The box argument is either
        a 2-tuple giving the upper left corner, a 4-tuple defining the
        left, upper, right, and lower pixel coordinate, or None (same as
        (0, 0)).  If a 4-tuple is given, the size of the pasted image
        must match the size of the region.

        If the modes don't match, the pasted image is converted to the mode of
        this image (see the :py:meth:`~PIL.Image.Image.convert` method for
        details).

        Instead of an image, the source can be a integer or tuple
        containing pixel values.  The method then fills the region
        with the given color.  When creating RGB images, you can
        also use color strings as supported by the ImageColor module.

        If a mask is given, this method updates only the regions
        indicated by the mask.  You can use either "1", "L" or "RGBA"
        images (in the latter case, the alpha band is used as mask).
        Where the mask is 255, the given image is copied as is.  Where
        the mask is 0, the current value is preserved.  Intermediate
        values can be used for transparency effects.

        Note that if you paste an "RGBA" image, the alpha band is
        ignored.  You can work around this by using the same image as
        both source image and mask.

        :param im: Source image or pixel value (integer or tuple).
        :param box: An optional 4-tuple giving the region to paste into.
           If a 2-tuple is used instead, it's treated as the upper left
           corner.  If omitted or None, the source is pasted into the
           upper left corner.

           If an image is given as the second argument and there is no
           third, the box defaults to (0, 0), and the second argument
           is interpreted as a mask image.
        :param mask: An optional mask image.
        """

        if isImageType(box) and mask is None:
            # abbreviated paste(im, mask) syntax
            mask = box
            box = None

        if box is None:
            # cover all of self
            box = (0, 0) + self.size

        if len(box) == 2:
            # lower left corner given; get size from image or mask
            if isImageType(im):
                size = im.size
            elif isImageType(mask):
                size = mask.size
            else:
                # FIXME: use self.size here?
                raise ValueError(
                    "cannot determine region size; use 4-item box"
                    )
            box = box + (box[0]+size[0], box[1]+size[1])

        if isStringType(im):
            from PIL import ImageColor
            im = ImageColor.getcolor(im, self.mode)

        elif isImageType(im):
            im.load()
            if self.mode != im.mode:
                if self.mode != "RGB" or im.mode not in ("RGBA", "RGBa"):
                    # should use an adapter for this!
                    im = im.convert(self.mode)
            im = im.im

        self.load()
        if self.readonly:
            self._copy()

        if mask:
            mask.load()
            self.im.paste(im, box, mask.im)
        else:
            self.im.paste(im, box)

    def point(self, lut, mode=None):
        """
        Maps this image through a lookup table or function.

        :param lut: A lookup table, containing 256 (or 65336 if
           self.mode=="I" and mode == "L") values per band in the
           image.  A function can be used instead, it should take a
           single argument. The function is called once for each
           possible pixel value, and the resulting table is applied to
           all bands of the image.
        :param mode: Output mode (default is same as input).  In the
           current version, this can only be used if the source image
           has mode "L" or "P", and the output has mode "1" or the
           source image mode is "I" and the output mode is "L".
        :returns: An :py:class:`~PIL.Image.Image` object.
        """

        self.load()

        if isinstance(lut, ImagePointHandler):
            return lut.point(self)

        if callable(lut):
            # if it isn't a list, it should be a function
            if self.mode in ("I", "I;16", "F"):
                # check if the function can be used with point_transform
                # UNDONE wiredfool -- I think this prevents us from ever doing
                # a gamma function point transform on > 8bit images.
                scale, offset = _getscaleoffset(lut)
                return self._new(self.im.point_transform(scale, offset))
            # for other modes, convert the function to a table
            lut = [lut(i) for i in range(256)] * self.im.bands

        if self.mode == "F":
            # FIXME: _imaging returns a confusing error message for this case
            raise ValueError("point operation not supported for this mode")

        return self._new(self.im.point(lut, mode))

    def putalpha(self, alpha):
        """
        Adds or replaces the alpha layer in this image.  If the image
        does not have an alpha layer, it's converted to "LA" or "RGBA".
        The new layer must be either "L" or "1".

        :param alpha: The new alpha layer.  This can either be an "L" or "1"
           image having the same size as this image, or an integer or
           other color value.
        """

        self.load()
        if self.readonly:
            self._copy()

        if self.mode not in ("LA", "RGBA"):
            # attempt to promote self to a matching alpha mode
            try:
                mode = getmodebase(self.mode) + "A"
                try:
                    self.im.setmode(mode)
                    self.pyaccess = None
                except (AttributeError, ValueError):
                    # do things the hard way
                    im = self.im.convert(mode)
                    if im.mode not in ("LA", "RGBA"):
                        raise ValueError  # sanity check
                    self.im = im
                    self.pyaccess = None
                self.mode = self.im.mode
            except (KeyError, ValueError):
                raise ValueError("illegal image mode")

        if self.mode == "LA":
            band = 1
        else:
            band = 3

        if isImageType(alpha):
            # alpha layer
            if alpha.mode not in ("1", "L"):
                raise ValueError("illegal image mode")
            alpha.load()
            if alpha.mode == "1":
                alpha = alpha.convert("L")
        else:
            # constant alpha
            try:
                self.im.fillband(band, alpha)
            except (AttributeError, ValueError):
                # do things the hard way
                alpha = new("L", self.size, alpha)
            else:
                return

        self.im.putband(alpha.im, band)

    def putdata(self, data, scale=1.0, offset=0.0):
        """
        Copies pixel data to this image.  This method copies data from a
        sequence object into the image, starting at the upper left
        corner (0, 0), and continuing until either the image or the
        sequence ends.  The scale and offset values are used to adjust
        the sequence values: **pixel = value*scale + offset**.

        :param data: A sequence object.
        :param scale: An optional scale value.  The default is 1.0.
        :param offset: An optional offset value.  The default is 0.0.
        """

        self.load()
        if self.readonly:
            self._copy()

        self.im.putdata(data, scale, offset)

    def putpalette(self, data, rawmode="RGB"):
        """
        Attaches a palette to this image.  The image must be a "P" or
        "L" image, and the palette sequence must contain 768 integer
        values, where each group of three values represent the red,
        green, and blue values for the corresponding pixel
        index. Instead of an integer sequence, you can use an 8-bit
        string.

        :param data: A palette sequence (either a list or a string).
        """
        from PIL import ImagePalette

        if self.mode not in ("L", "P"):
            raise ValueError("illegal image mode")
        self.load()
        if isinstance(data, ImagePalette.ImagePalette):
            palette = ImagePalette.raw(data.rawmode, data.palette)
        else:
            if not isinstance(data, bytes):
                if bytes is str:
                    data = "".join(chr(x) for x in data)
                else:
                    data = bytes(data)
            palette = ImagePalette.raw(rawmode, data)
        self.mode = "P"
        self.palette = palette
        self.palette.mode = "RGB"
        self.load()  # install new palette

    def putpixel(self, xy, value):
        """
        Modifies the pixel at the given position. The color is given as
        a single numerical value for single-band images, and a tuple for
        multi-band images.

        Note that this method is relatively slow.  For more extensive changes,
        use :py:meth:`~PIL.Image.Image.paste` or the :py:mod:`~PIL.ImageDraw`
        module instead.

        See:

        * :py:meth:`~PIL.Image.Image.paste`
        * :py:meth:`~PIL.Image.Image.putdata`
        * :py:mod:`~PIL.ImageDraw`

        :param xy: The pixel coordinate, given as (x, y).
        :param value: The pixel value.
        """

        self.load()
        if self.readonly:
            self._copy()
            self.pyaccess = None
            self.load()

        if self.pyaccess:
            return self.pyaccess.putpixel(xy, value)
        return self.im.putpixel(xy, value)

    def resize(self, size, resample=NEAREST):
        """
        Returns a resized copy of this image.

        :param size: The requested size in pixels, as a 2-tuple:
           (width, height).
        :param resample: An optional resampling filter.  This can be
           one of :py:attr:`PIL.Image.NEAREST` (use nearest neighbour),
           :py:attr:`PIL.Image.BILINEAR` (linear interpolation in a 2x2
           environment), :py:attr:`PIL.Image.BICUBIC` (cubic spline
           interpolation in a 4x4 environment), or
           :py:attr:`PIL.Image.ANTIALIAS` (a high-quality downsampling filter).
           If omitted, or if the image has mode "1" or "P", it is
           set :py:attr:`PIL.Image.NEAREST`.
        :returns: An :py:class:`~PIL.Image.Image` object.
        """

        if resample not in (NEAREST, BILINEAR, BICUBIC, ANTIALIAS):
            raise ValueError("unknown resampling filter")

        self.load()

        size=tuple(size)
        if self.size == size:
            return self._new(self.im)

        if self.mode in ("1", "P"):
            resample = NEAREST

        if self.mode == 'RGBA':
            return self.convert('RGBa').resize(size, resample).convert('RGBA')

        if resample == ANTIALIAS:
            # requires stretch support (imToolkit & PIL 1.1.3)
            try:
                im = self.im.stretch(size, resample)
            except AttributeError:
                raise ValueError("unsupported resampling filter")
        else:
            im = self.im.resize(size, resample)

        return self._new(im)

    def rotate(self, angle, resample=NEAREST, expand=0):
        """
        Returns a rotated copy of this image.  This method returns a
        copy of this image, rotated the given number of degrees counter
        clockwise around its centre.

        :param angle: In degrees counter clockwise.
        :param resample: An optional resampling filter.  This can be
           one of :py:attr:`PIL.Image.NEAREST` (use nearest neighbour),
           :py:attr:`PIL.Image.BILINEAR` (linear interpolation in a 2x2
           environment), or :py:attr:`PIL.Image.BICUBIC`
           (cubic spline interpolation in a 4x4 environment).
           If omitted, or if the image has mode "1" or "P", it is
           set :py:attr:`PIL.Image.NEAREST`.
        :param expand: Optional expansion flag.  If true, expands the output
           image to make it large enough to hold the entire rotated image.
           If false or omitted, make the output image the same size as the
           input image.
        :returns: An :py:class:`~PIL.Image.Image` object.
        """

        if expand:
            import math
            angle = -angle * math.pi / 180
            matrix = [
                math.cos(angle), math.sin(angle), 0.0,
                -math.sin(angle), math.cos(angle), 0.0
                ]

            def transform(x, y, matrix=matrix):
                (a, b, c, d, e, f) = matrix
                return a*x + b*y + c, d*x + e*y + f

            # calculate output size
            w, h = self.size
            xx = []
            yy = []
            for x, y in ((0, 0), (w, 0), (w, h), (0, h)):
                x, y = transform(x, y)
                xx.append(x)
                yy.append(y)
            w = int(math.ceil(max(xx)) - math.floor(min(xx)))
            h = int(math.ceil(max(yy)) - math.floor(min(yy)))

            # adjust center
            x, y = transform(w / 2.0, h / 2.0)
            matrix[2] = self.size[0] / 2.0 - x
            matrix[5] = self.size[1] / 2.0 - y

            return self.transform((w, h), AFFINE, matrix, resample)

        if resample not in (NEAREST, BILINEAR, BICUBIC):
            raise ValueError("unknown resampling filter")

        self.load()

        if self.mode in ("1", "P"):
            resample = NEAREST

        return self._new(self.im.rotate(angle, resample))

    def save(self, fp, format=None, **params):
        """
        Saves this image under the given filename.  If no format is
        specified, the format to use is determined from the filename
        extension, if possible.

        Keyword options can be used to provide additional instructions
        to the writer. If a writer doesn't recognise an option, it is
        silently ignored. The available options are described later in
        this handbook.

        You can use a file object instead of a filename. In this case,
        you must always specify the format. The file object must
        implement the **seek**, **tell**, and **write**
        methods, and be opened in binary mode.

        :param file: File name or file object.
        :param format: Optional format override.  If omitted, the
           format to use is determined from the filename extension.
           If a file object was used instead of a filename, this
           parameter should always be used.
        :param options: Extra parameters to the image writer.
        :returns: None
        :exception KeyError: If the output format could not be determined
           from the file name.  Use the format option to solve this.
        :exception IOError: If the file could not be written.  The file
           may have been created, and may contain partial data.
        """

        if isPath(fp):
            filename = fp
        else:
            if hasattr(fp, "name") and isPath(fp.name):
                filename = fp.name
            else:
                filename = ""

        # may mutate self!
        self.load()

        self.encoderinfo = params
        self.encoderconfig = ()

        preinit()

        ext = os.path.splitext(filename)[1].lower()

        if not format:
            try:
                format = EXTENSION[ext]
            except KeyError:
                init()
                try:
                    format = EXTENSION[ext]
                except KeyError:
                    raise KeyError(ext)  # unknown extension

        try:
            save_handler = SAVE[format.upper()]
        except KeyError:
            init()
            save_handler = SAVE[format.upper()]  # unknown format

        if isPath(fp):
            fp = builtins.open(fp, "wb")
            close = 1
        else:
            close = 0

        try:
            save_handler(self, fp, filename)
        finally:
            # do what we can to clean up
            if close:
                fp.close()

    def seek(self, frame):
        """
        Seeks to the given frame in this sequence file. If you seek
        beyond the end of the sequence, the method raises an
        **EOFError** exception. When a sequence file is opened, the
        library automatically seeks to frame 0.

        Note that in the current version of the library, most sequence
        formats only allows you to seek to the next frame.

        See :py:meth:`~PIL.Image.Image.tell`.

        :param frame: Frame number, starting at 0.
        :exception EOFError: If the call attempts to seek beyond the end
            of the sequence.
        """

        # overridden by file handlers
        if frame != 0:
            raise EOFError

    def show(self, title=None, command=None):
        """
        Displays this image. This method is mainly intended for
        debugging purposes.

        On Unix platforms, this method saves the image to a temporary
        PPM file, and calls the **xv** utility.

        On Windows, it saves the image to a temporary BMP file, and uses
        the standard BMP display utility to show it (usually Paint).

        :param title: Optional title to use for the image window,
           where possible.
        :param command: command used to show the image
        """

        _show(self, title=title, command=command)

    def split(self):
        """
        Split this image into individual bands. This method returns a
        tuple of individual image bands from an image. For example,
        splitting an "RGB" image creates three new images each
        containing a copy of one of the original bands (red, green,
        blue).

        :returns: A tuple containing bands.
        """

        self.load()
        if self.im.bands == 1:
            ims = [self.copy()]
        else:
            ims = []
            for i in range(self.im.bands):
                ims.append(self._new(self.im.getband(i)))
        return tuple(ims)

    def tell(self):
        """
        Returns the current frame number. See :py:meth:`~PIL.Image.Image.seek`.

        :returns: Frame number, starting with 0.
        """
        return 0

    def thumbnail(self, size, resample=ANTIALIAS):
        """
        Make this image into a thumbnail.  This method modifies the
        image to contain a thumbnail version of itself, no larger than
        the given size.  This method calculates an appropriate thumbnail
        size to preserve the aspect of the image, calls the
        :py:meth:`~PIL.Image.Image.draft` method to configure the file reader
        (where applicable), and finally resizes the image.

        Note that the bilinear and bicubic filters in the current
        version of PIL are not well-suited for thumbnail generation.
        You should use :py:attr:`PIL.Image.ANTIALIAS` unless speed is much more
        important than quality.

        Also note that this function modifies the :py:class:`~PIL.Image.Image`
        object in place.  If you need to use the full resolution image as well,
        apply this method to a :py:meth:`~PIL.Image.Image.copy` of the original
        image.

        :param size: Requested size.
        :param resample: Optional resampling filter.  This can be one
           of :py:attr:`PIL.Image.NEAREST`, :py:attr:`PIL.Image.BILINEAR`,
           :py:attr:`PIL.Image.BICUBIC`, or :py:attr:`PIL.Image.ANTIALIAS`
           (best quality).  If omitted, it defaults to
           :py:attr:`PIL.Image.ANTIALIAS`. (was :py:attr:`PIL.Image.NEAREST`
           prior to version 2.5.0)
        :returns: None
        """

        # preserve aspect ratio
        x, y = self.size
        if x > size[0]:
            y = int(max(y * size[0] / x, 1))
            x = int(size[0])
        if y > size[1]:
            x = int(max(x * size[1] / y, 1))
            y = int(size[1])
        size = x, y

        if size == self.size:
            return

        self.draft(None, size)

        self.load()

        try:
            im = self.resize(size, resample)
        except ValueError:
            if resample != ANTIALIAS:
                raise
            im = self.resize(size, NEAREST)  # fallback

        self.im = im.im
        self.mode = im.mode
        self.size = size

        self.readonly = 0
        self.pyaccess = None

    # FIXME: the different tranform methods need further explanation
    # instead of bloating the method docs, add a separate chapter.
    def transform(self, size, method, data=None, resample=NEAREST, fill=1):
        """
        Transforms this image.  This method creates a new image with the
        given size, and the same mode as the original, and copies data
        to the new image using the given transform.

        :param size: The output size.
        :param method: The transformation method.  This is one of
          :py:attr:`PIL.Image.EXTENT` (cut out a rectangular subregion),
          :py:attr:`PIL.Image.AFFINE` (affine transform),
          :py:attr:`PIL.Image.PERSPECTIVE` (perspective transform),
          :py:attr:`PIL.Image.QUAD` (map a quadrilateral to a rectangle), or
          :py:attr:`PIL.Image.MESH` (map a number of source quadrilaterals
          in one operation).
        :param data: Extra data to the transformation method.
        :param resample: Optional resampling filter.  It can be one of
           :py:attr:`PIL.Image.NEAREST` (use nearest neighbour),
           :py:attr:`PIL.Image.BILINEAR` (linear interpolation in a 2x2
           environment), or :py:attr:`PIL.Image.BICUBIC` (cubic spline
           interpolation in a 4x4 environment). If omitted, or if the image
           has mode "1" or "P", it is set to :py:attr:`PIL.Image.NEAREST`.
        :returns: An :py:class:`~PIL.Image.Image` object.
        """

        if self.mode == 'RGBA':
            return self.convert('RGBa').transform(
                size, method, data, resample, fill).convert('RGBA')

        if isinstance(method, ImageTransformHandler):
            return method.transform(size, self, resample=resample, fill=fill)
        if hasattr(method, "getdata"):
            # compatibility w. old-style transform objects
            method, data = method.getdata()
        if data is None:
            raise ValueError("missing method data")

        im = new(self.mode, size, None)
        if method == MESH:
            # list of quads
            for box, quad in data:
                im.__transformer(box, self, QUAD, quad, resample, fill)
        else:
            im.__transformer((0, 0)+size, self, method, data, resample, fill)

        return im

    def __transformer(self, box, image, method, data,
                      resample=NEAREST, fill=1):

        # FIXME: this should be turned into a lazy operation (?)

        w = box[2]-box[0]
        h = box[3]-box[1]

        if method == AFFINE:
            # change argument order to match implementation
            data = (data[2], data[0], data[1],
                    data[5], data[3], data[4])
        elif method == EXTENT:
            # convert extent to an affine transform
            x0, y0, x1, y1 = data
            xs = float(x1 - x0) / w
            ys = float(y1 - y0) / h
            method = AFFINE
            data = (x0 + xs/2, xs, 0, y0 + ys/2, 0, ys)
        elif method == PERSPECTIVE:
            # change argument order to match implementation
            data = (data[2], data[0], data[1],
                    data[5], data[3], data[4],
                    data[6], data[7])
        elif method == QUAD:
            # quadrilateral warp.  data specifies the four corners
            # given as NW, SW, SE, and NE.
            nw = data[0:2]
            sw = data[2:4]
            se = data[4:6]
            ne = data[6:8]
            x0, y0 = nw
            As = 1.0 / w
            At = 1.0 / h
            data = (x0, (ne[0]-x0)*As, (sw[0]-x0)*At,
                    (se[0]-sw[0]-ne[0]+x0)*As*At,
                    y0, (ne[1]-y0)*As, (sw[1]-y0)*At,
                    (se[1]-sw[1]-ne[1]+y0)*As*At)
        else:
            raise ValueError("unknown transformation method")

        if resample not in (NEAREST, BILINEAR, BICUBIC):
            raise ValueError("unknown resampling filter")

        image.load()

        self.load()

        if image.mode in ("1", "P"):
            resample = NEAREST

        self.im.transform2(box, image.im, method, data, resample, fill)

    def transpose(self, method):
        """
        Transpose image (flip or rotate in 90 degree steps)

        :param method: One of :py:attr:`PIL.Image.FLIP_LEFT_RIGHT`,
          :py:attr:`PIL.Image.FLIP_TOP_BOTTOM`, :py:attr:`PIL.Image.ROTATE_90`,
          :py:attr:`PIL.Image.ROTATE_180`, or :py:attr:`PIL.Image.ROTATE_270`.
        :returns: Returns a flipped or rotated copy of this image.
        """

        self.load()
        im = self.im.transpose(method)
        return self._new(im)

    def effect_spread(self, distance):
        """
        Randomly spread pixels in an image.

        :param distance: Distance to spread pixels.
        """
        self.load()
        im = self.im.effect_spread(distance)
        return self._new(im)


# --------------------------------------------------------------------
# Lazy operations

class _ImageCrop(Image):

    def __init__(self, im, box):

        Image.__init__(self)

        x0, y0, x1, y1 = box
        if x1 < x0:
            x1 = x0
        if y1 < y0:
            y1 = y0

        self.mode = im.mode
        self.size = x1-x0, y1-y0

        self.__crop = x0, y0, x1, y1

        self.im = im.im

    def load(self):

        # lazy evaluation!
        if self.__crop:
            self.im = self.im.crop(self.__crop)
            self.__crop = None

        if self.im:
            return self.im.pixel_access(self.readonly)

        # FIXME: future versions should optimize crop/paste
        # sequences!


# --------------------------------------------------------------------
# Abstract handlers.

class ImagePointHandler:
    # used as a mixin by point transforms (for use with im.point)
    pass


class ImageTransformHandler:
    # used as a mixin by geometry transforms (for use with im.transform)
    pass


# --------------------------------------------------------------------
# Factories

#
# Debugging

def _wedge():
    "Create greyscale wedge (for debugging only)"

    return Image()._new(core.wedge("L"))


def new(mode, size, color=0):
    """
    Creates a new image with the given mode and size.

    :param mode: The mode to use for the new image.
    :param size: A 2-tuple, containing (width, height) in pixels.
    :param color: What color to use for the image.  Default is black.
       If given, this should be a single integer or floating point value
       for single-band modes, and a tuple for multi-band modes (one value
       per band).  When creating RGB images, you can also use color
       strings as supported by the ImageColor module.  If the color is
       None, the image is not initialised.
    :returns: An :py:class:`~PIL.Image.Image` object.
    """

    if color is None:
        # don't initialize
        return Image()._new(core.new(mode, size))

    if isStringType(color):
        # css3-style specifier

        from PIL import ImageColor
        color = ImageColor.getcolor(color, mode)

    return Image()._new(core.fill(mode, size, color))


def frombytes(mode, size, data, decoder_name="raw", *args):
    """
    Creates a copy of an image memory from pixel data in a buffer.

    In its simplest form, this function takes three arguments
    (mode, size, and unpacked pixel data).

    You can also use any pixel decoder supported by PIL.  For more
    information on available decoders, see the section
    **Writing Your Own File Decoder**.

    Note that this function decodes pixel data only, not entire images.
    If you have an entire image in a string, wrap it in a
    :py:class:`~io.BytesIO` object, and use :py:func:`~PIL.Image.open` to load
    it.

    :param mode: The image mode.
    :param size: The image size.
    :param data: A byte buffer containing raw data for the given mode.
    :param decoder_name: What decoder to use.
    :param args: Additional parameters for the given decoder.
    :returns: An :py:class:`~PIL.Image.Image` object.
    """

    # may pass tuple instead of argument list
    if len(args) == 1 and isinstance(args[0], tuple):
        args = args[0]

    if decoder_name == "raw" and args == ():
        args = mode

    im = new(mode, size)
    im.frombytes(data, decoder_name, args)
    return im


def fromstring(*args, **kw):
    """Deprecated alias to frombytes.

    .. deprecated:: 2.0
    """
    warnings.warn(
        'fromstring() is deprecated. Please call frombytes() instead.',
        DeprecationWarning,
        stacklevel=2
    )
    return frombytes(*args, **kw)


def frombuffer(mode, size, data, decoder_name="raw", *args):
    """
    Creates an image memory referencing pixel data in a byte buffer.

    This function is similar to :py:func:`~PIL.Image.frombytes`, but uses data
    in the byte buffer, where possible.  This means that changes to the
    original buffer object are reflected in this image).  Not all modes can
    share memory; supported modes include "L", "RGBX", "RGBA", and "CMYK".

    Note that this function decodes pixel data only, not entire images.
    If you have an entire image file in a string, wrap it in a
    **BytesIO** object, and use :py:func:`~PIL.Image.open` to load it.

    In the current version, the default parameters used for the "raw" decoder
    differs from that used for :py:func:`~PIL.Image.fromstring`.  This is a
    bug, and will probably be fixed in a future release.  The current release
    issues a warning if you do this; to disable the warning, you should provide
    the full set of parameters.  See below for details.

    :param mode: The image mode.
    :param size: The image size.
    :param data: A bytes or other buffer object containing raw
        data for the given mode.
    :param decoder_name: What decoder to use.
    :param args: Additional parameters for the given decoder.  For the
        default encoder ("raw"), it's recommended that you provide the
        full set of parameters::

            frombuffer(mode, size, data, "raw", mode, 0, 1)

    :returns: An :py:class:`~PIL.Image.Image` object.

    .. versionadded:: 1.1.4
    """

    # may pass tuple instead of argument list
    if len(args) == 1 and isinstance(args[0], tuple):
        args = args[0]

    if decoder_name == "raw":
        if args == ():
            if warnings:
                warnings.warn(
                    "the frombuffer defaults may change in a future release; "
                    "for portability, change the call to read:\n"
                    "  frombuffer(mode, size, data, 'raw', mode, 0, 1)",
                    RuntimeWarning, stacklevel=2
                )
            args = mode, 0, -1  # may change to (mode, 0, 1) post-1.1.6
        if args[0] in _MAPMODES:
            im = new(mode, (1, 1))
            im = im._new(
                core.map_buffer(data, size, decoder_name, None, 0, args)
                )
            im.readonly = 1
            return im

    return frombytes(mode, size, data, decoder_name, args)


def fromarray(obj, mode=None):
    """
    Creates an image memory from an object exporting the array interface
    (using the buffer protocol).

    If obj is not contiguous, then the tobytes method is called
    and :py:func:`~PIL.Image.frombuffer` is used.

    :param obj: Object with array interface
    :param mode: Mode to use (will be determined from type if None)
    :returns: An image memory.

    .. versionadded:: 1.1.6
    """
    arr = obj.__array_interface__
    shape = arr['shape']
    ndim = len(shape)
    try:
        strides = arr['strides']
    except KeyError:
        strides = None
    if mode is None:
        try:
            typekey = (1, 1) + shape[2:], arr['typestr']
            mode, rawmode = _fromarray_typemap[typekey]
        except KeyError:
            # print typekey
            raise TypeError("Cannot handle this data type")
    else:
        rawmode = mode
    if mode in ["1", "L", "I", "P", "F"]:
        ndmax = 2
    elif mode == "RGB":
        ndmax = 3
    else:
        ndmax = 4
    if ndim > ndmax:
        raise ValueError("Too many dimensions: %d > %d." % (ndim, ndmax))

    size = shape[1], shape[0]
    if strides is not None:
        if hasattr(obj, 'tobytes'):
            obj = obj.tobytes()
        else:
            obj = obj.tostring()

    return frombuffer(mode, size, obj, "raw", rawmode, 0, 1)

_fromarray_typemap = {
    # (shape, typestr) => mode, rawmode
    # first two members of shape are set to one
    # ((1, 1), "|b1"): ("1", "1"), # broken
    ((1, 1), "|u1"): ("L", "L"),
    ((1, 1), "|i1"): ("I", "I;8"),
    ((1, 1), "<i2"): ("I", "I;16"),
    ((1, 1), ">i2"): ("I", "I;16B"),
    ((1, 1), "<i4"): ("I", "I;32"),
    ((1, 1), ">i4"): ("I", "I;32B"),
    ((1, 1), "<f4"): ("F", "F;32F"),
    ((1, 1), ">f4"): ("F", "F;32BF"),
    ((1, 1), "<f8"): ("F", "F;64F"),
    ((1, 1), ">f8"): ("F", "F;64BF"),
    ((1, 1, 3), "|u1"): ("RGB", "RGB"),
    ((1, 1, 4), "|u1"): ("RGBA", "RGBA"),
    }

# shortcuts
_fromarray_typemap[((1, 1), _ENDIAN + "i4")] = ("I", "I")
_fromarray_typemap[((1, 1), _ENDIAN + "f4")] = ("F", "F")


def _decompression_bomb_check(size):
    if MAX_IMAGE_PIXELS is None:
        return

    pixels = size[0] * size[1]

    if pixels > MAX_IMAGE_PIXELS:
        warnings.warn(
            "Image size (%d pixels) exceeds limit of %d pixels, "
            "could be decompression bomb DOS attack." %
            (pixels, MAX_IMAGE_PIXELS),
            DecompressionBombWarning)


def open(fp, mode="r"):
    """
    Opens and identifies the given image file.

    This is a lazy operation; this function identifies the file, but
    the file remains open and the actual image data is not read from
    the file until you try to process the data (or call the
    :py:meth:`~PIL.Image.Image.load` method).  See
    :py:func:`~PIL.Image.new`.

    :param file: A filename (string) or a file object.  The file object
       must implement :py:meth:`~file.read`, :py:meth:`~file.seek`, and
       :py:meth:`~file.tell` methods, and be opened in binary mode.
    :param mode: The mode.  If given, this argument must be "r".
    :returns: An :py:class:`~PIL.Image.Image` object.
    :exception IOError: If the file cannot be found, or the image cannot be
       opened and identified.
    """

    if mode != "r":
        raise ValueError("bad mode %r" % mode)

    if isPath(fp):
        filename = fp
        fp = builtins.open(fp, "rb")
    else:
        filename = ""

    prefix = fp.read(16)

    preinit()

    for i in ID:
        try:
            factory, accept = OPEN[i]
            if not accept or accept(prefix):
                fp.seek(0)
                im = factory(fp, filename)
                _decompression_bomb_check(im.size)
                return im
        except (SyntaxError, IndexError, TypeError):
            # import traceback
            # traceback.print_exc()
            pass

    if init():

        for i in ID:
            try:
                factory, accept = OPEN[i]
                if not accept or accept(prefix):
                    fp.seek(0)
                    im = factory(fp, filename)
                    _decompression_bomb_check(im.size)
                    return im
            except (SyntaxError, IndexError, TypeError):
                # import traceback
                # traceback.print_exc()
                pass

    raise IOError("cannot identify image file %r"
                  % (filename if filename else fp))


#
# Image processing.

def alpha_composite(im1, im2):
    """
    Alpha composite im2 over im1.

    :param im1: The first image.
    :param im2: The second image.  Must have the same mode and size as
       the first image.
    :returns: An :py:class:`~PIL.Image.Image` object.
    """

    im1.load()
    im2.load()
    return im1._new(core.alpha_composite(im1.im, im2.im))


def blend(im1, im2, alpha):
    """
    Creates a new image by interpolating between two input images, using
    a constant alpha.::

        out = image1 * (1.0 - alpha) + image2 * alpha

    :param im1: The first image.
    :param im2: The second image.  Must have the same mode and size as
       the first image.
    :param alpha: The interpolation alpha factor.  If alpha is 0.0, a
       copy of the first image is returned. If alpha is 1.0, a copy of
       the second image is returned. There are no restrictions on the
       alpha value. If necessary, the result is clipped to fit into
       the allowed output range.
    :returns: An :py:class:`~PIL.Image.Image` object.
    """

    im1.load()
    im2.load()
    return im1._new(core.blend(im1.im, im2.im, alpha))


def composite(image1, image2, mask):
    """
    Create composite image by blending images using a transparency mask.

    :param image1: The first image.
    :param image2: The second image.  Must have the same mode and
       size as the first image.
    :param mask: A mask image.  This image can can have mode
       "1", "L", or "RGBA", and must have the same size as the
       other two images.
    """

    image = image2.copy()
    image.paste(image1, None, mask)
    return image


def eval(image, *args):
    """
    Applies the function (which should take one argument) to each pixel
    in the given image. If the image has more than one band, the same
    function is applied to each band. Note that the function is
    evaluated once for each possible pixel value, so you cannot use
    random components or other generators.

    :param image: The input image.
    :param function: A function object, taking one integer argument.
    :returns: An :py:class:`~PIL.Image.Image` object.
    """

    return image.point(args[0])


def merge(mode, bands):
    """
    Merge a set of single band images into a new multiband image.

    :param mode: The mode to use for the output image.
    :param bands: A sequence containing one single-band image for
        each band in the output image.  All bands must have the
        same size.
    :returns: An :py:class:`~PIL.Image.Image` object.
    """

    if getmodebands(mode) != len(bands) or "*" in mode:
        raise ValueError("wrong number of bands")
    for im in bands[1:]:
        if im.mode != getmodetype(mode):
            raise ValueError("mode mismatch")
        if im.size != bands[0].size:
            raise ValueError("size mismatch")
    im = core.new(mode, bands[0].size)
    for i in range(getmodebands(mode)):
        bands[i].load()
        im.putband(bands[i].im, i)
    return bands[0]._new(im)


# --------------------------------------------------------------------
# Plugin registry

def register_open(id, factory, accept=None):
    """
    Register an image file plugin.  This function should not be used
    in application code.

    :param id: An image format identifier.
    :param factory: An image file factory method.
    :param accept: An optional function that can be used to quickly
       reject images having another format.
    """
    id = id.upper()
    ID.append(id)
    OPEN[id] = factory, accept


def register_mime(id, mimetype):
    """
    Registers an image MIME type.  This function should not be used
    in application code.

    :param id: An image format identifier.
    :param mimetype: The image MIME type for this format.
    """
    MIME[id.upper()] = mimetype


def register_save(id, driver):
    """
    Registers an image save function.  This function should not be
    used in application code.

    :param id: An image format identifier.
    :param driver: A function to save images in this format.
    """
    SAVE[id.upper()] = driver


def register_extension(id, extension):
    """
    Registers an image extension.  This function should not be
    used in application code.

    :param id: An image format identifier.
    :param extension: An extension used for this format.
    """
    EXTENSION[extension.lower()] = id.upper()


# --------------------------------------------------------------------
# Simple display support.  User code may override this.

def _show(image, **options):
    # override me, as necessary
    _showxv(image, **options)


def _showxv(image, title=None, **options):
    from PIL import ImageShow
    ImageShow.show(image, title, **options)


# --------------------------------------------------------------------
# Effects

def effect_mandelbrot(size, extent, quality):
    """
    Generate a Mandelbrot set covering the given extent.

    :param size: The requested size in pixels, as a 2-tuple:
       (width, height).
    :param extent: The extent to cover, as a 4-tuple:
       (x0, y0, x1, y2).
    :param quality: Quality.
    """
    return Image()._new(core.effect_mandelbrot(size, extent, quality))


def effect_noise(size, sigma):
    """
    Generate Gaussian noise centered around 128.

    :param size: The requested size in pixels, as a 2-tuple:
       (width, height).
    :param sigma: Standard deviation of noise.
    """
    return Image()._new(core.effect_noise(size, sigma))

# End of file


#
# The Python Imaging Library.
# $Id$
#
# PIL raster font management
#
# History:
# 1996-08-07 fl   created (experimental)
# 1997-08-25 fl   minor adjustments to handle fonts from pilfont 0.3
# 1999-02-06 fl   rewrote most font management stuff in C
# 1999-03-17 fl   take pth files into account in load_path (from Richard Jones)
# 2001-02-17 fl   added freetype support
# 2001-05-09 fl   added TransposedFont wrapper class
# 2002-03-04 fl   make sure we have a "L" or "1" font
# 2002-12-04 fl   skip non-directory entries in the system path
# 2003-04-29 fl   add embedded default font
# 2003-09-27 fl   added support for truetype charmap encodings
#
# Todo:
# Adapt to PILFONT2 format (16-bit fonts, compressed, single file)
#
# Copyright (c) 1997-2003 by Secret Labs AB
# Copyright (c) 1996-2003 by Fredrik Lundh
#
# See the README file for information on usage and redistribution.
#

from __future__ import print_function

from PIL import Image
from PIL._util import isDirectory, isPath
import os
import sys

try:
    import warnings
except ImportError:
    warnings = None


class _imagingft_not_installed:
    # module placeholder
    def __getattr__(self, id):
        raise ImportError("The _imagingft C module is not installed")

try:
    from PIL import _imagingft as core
except ImportError:
    core = _imagingft_not_installed()

# FIXME: add support for pilfont2 format (see FontFile.py)

# --------------------------------------------------------------------
# Font metrics format:
#       "PILfont" LF
#       fontdescriptor LF
#       (optional) key=value... LF
#       "DATA" LF
#       binary data: 256*10*2 bytes (dx, dy, dstbox, srcbox)
#
# To place a character, cut out srcbox and paste at dstbox,
# relative to the character position.  Then move the character
# position according to dx, dy.
# --------------------------------------------------------------------


class ImageFont:
    "PIL font wrapper"

    def _load_pilfont(self, filename):

        file = open(filename, "rb")

        for ext in (".png", ".gif", ".pbm"):
            try:
                fullname = os.path.splitext(filename)[0] + ext
                image = Image.open(fullname)
            except:
                pass
            else:
                if image and image.mode in ("1", "L"):
                    break
        else:
            raise IOError("cannot find glyph data file")

        self.file = fullname

        return self._load_pilfont_data(file, image)

    def _load_pilfont_data(self, file, image):

        # read PILfont header
        if file.readline() != b"PILfont\n":
            raise SyntaxError("Not a PILfont file")
        file.readline().split(b";")
        self.info = []  # FIXME: should be a dictionary
        while True:
            s = file.readline()
            if not s or s == b"DATA\n":
                break
            self.info.append(s)

        # read PILfont metrics
        data = file.read(256*20)

        # check image
        if image.mode not in ("1", "L"):
            raise TypeError("invalid font image mode")

        image.load()

        self.font = Image.core.font(image.im, data)

        # delegate critical operations to internal type
        self.getsize = self.font.getsize
        self.getmask = self.font.getmask


##
# Wrapper for FreeType fonts.  Application code should use the
# <b>truetype</b> factory function to create font objects.

class FreeTypeFont:
    "FreeType font wrapper (requires _imagingft service)"

    def __init__(self, font=None, size=10, index=0, encoding="", file=None):
        # FIXME: use service provider instead
        if file:
            if warnings:
                warnings.warn(
                    'file parameter deprecated, '
                    'please use font parameter instead.',
                    DeprecationWarning)
            font = file

        if isPath(font):
            self.font = core.getfont(font, size, index, encoding)
        else:
            self.font_bytes = font.read()
            self.font = core.getfont(
                "", size, index, encoding, self.font_bytes)

    def getname(self):
        return self.font.family, self.font.style

    def getmetrics(self):
        return self.font.ascent, self.font.descent

    def getsize(self, text):
        size, offset = self.font.getsize(text)
        return (size[0] + offset[0], size[1] + offset[1])

    def getoffset(self, text):
        return self.font.getsize(text)[1]

    def getmask(self, text, mode=""):
        return self.getmask2(text, mode)[0]

    def getmask2(self, text, mode="", fill=Image.core.fill):
        size, offset = self.font.getsize(text)
        im = fill("L", size, 0)
        self.font.render(text, im.id, mode == "1")
        return im, offset

##
# Wrapper that creates a transposed font from any existing font
# object.
#
# @param font A font object.
# @param orientation An optional orientation.  If given, this should
#     be one of Image.FLIP_LEFT_RIGHT, Image.FLIP_TOP_BOTTOM,
#     Image.ROTATE_90, Image.ROTATE_180, or Image.ROTATE_270.


class TransposedFont:
    "Wrapper for writing rotated or mirrored text"

    def __init__(self, font, orientation=None):
        self.font = font
        self.orientation = orientation  # any 'transpose' argument, or None

    def getsize(self, text):
        w, h = self.font.getsize(text)
        if self.orientation in (Image.ROTATE_90, Image.ROTATE_270):
            return h, w
        return w, h

    def getmask(self, text, mode=""):
        im = self.font.getmask(text, mode)
        if self.orientation is not None:
            return im.transpose(self.orientation)
        return im


def load(filename):
    """
    Load a font file.  This function loads a font object from the given
    bitmap font file, and returns the corresponding font object.

    :param filename: Name of font file.
    :return: A font object.
    :exception IOError: If the file could not be read.
    """
    f = ImageFont()
    f._load_pilfont(filename)
    return f


def truetype(font=None, size=10, index=0, encoding="", filename=None):
    """
    Load a TrueType or OpenType font file, and create a font object.
    This function loads a font object from the given file, and creates
    a font object for a font of the given size.

    This function requires the _imagingft service.

    :param filename: A truetype font file. Under Windows, if the file
                     is not found in this filename, the loader also looks in
                     Windows :file:`fonts/` directory.
    :param size: The requested size, in points.
    :param index: Which font face to load (default is first available face).
    :param encoding: Which font encoding to use (default is Unicode). Common
                     encodings are "unic" (Unicode), "symb" (Microsoft
                     Symbol), "ADOB" (Adobe Standard), "ADBE" (Adobe Expert),
                     and "armn" (Apple Roman). See the FreeType documentation
                     for more information.
    :return: A font object.
    :exception IOError: If the file could not be read.
    """

    if filename:
        if warnings:
            warnings.warn(
                'filename parameter deprecated, '
                'please use font parameter instead.',
                DeprecationWarning)
        font = filename

    try:
        return FreeTypeFont(font, size, index, encoding)
    except IOError:
        if sys.platform == "win32":
            # check the windows font repository
            # NOTE: must use uppercase WINDIR, to work around bugs in
            # 1.5.2's os.environ.get()
            windir = os.environ.get("WINDIR")
            if windir:
                filename = os.path.join(windir, "fonts", font)
                return FreeTypeFont(filename, size, index, encoding)
        raise


def load_path(filename):
    """
    Load font file. Same as :py:func:`~PIL.ImageFont.load`, but searches for a
    bitmap font along the Python path.

    :param filename: Name of font file.
    :return: A font object.
    :exception IOError: If the file could not be read.
    """
    for dir in sys.path:
        if isDirectory(dir):
            if not isinstance(filename, str):
                if bytes is str:
                    filename = filename.encode("utf-8")
                else:
                    filename = filename.decode("utf-8")
            try:
                return load(os.path.join(dir, filename))
            except IOError:
                pass
    raise IOError("cannot find font file")


def load_default():
    """Load a "better than nothing" default font.

    .. versionadded:: 1.1.4

    :return: A font object.
    """
    from io import BytesIO
    import base64
    f = ImageFont()
    f._load_pilfont_data(
        # courB08
        BytesIO(base64.decodestring(b'''
UElMZm9udAo7Ozs7OzsxMDsKREFUQQoAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAYAAAAA//8AAQAAAAAAAAABAAEA
BgAAAAH/+gADAAAAAQAAAAMABgAGAAAAAf/6AAT//QADAAAABgADAAYAAAAA//kABQABAAYAAAAL
AAgABgAAAAD/+AAFAAEACwAAABAACQAGAAAAAP/5AAUAAAAQAAAAFQAHAAYAAP////oABQAAABUA
AAAbAAYABgAAAAH/+QAE//wAGwAAAB4AAwAGAAAAAf/5AAQAAQAeAAAAIQAIAAYAAAAB//kABAAB
ACEAAAAkAAgABgAAAAD/+QAE//0AJAAAACgABAAGAAAAAP/6AAX//wAoAAAALQAFAAYAAAAB//8A
BAACAC0AAAAwAAMABgAAAAD//AAF//0AMAAAADUAAQAGAAAAAf//AAMAAAA1AAAANwABAAYAAAAB
//kABQABADcAAAA7AAgABgAAAAD/+QAFAAAAOwAAAEAABwAGAAAAAP/5AAYAAABAAAAARgAHAAYA
AAAA//kABQAAAEYAAABLAAcABgAAAAD/+QAFAAAASwAAAFAABwAGAAAAAP/5AAYAAABQAAAAVgAH
AAYAAAAA//kABQAAAFYAAABbAAcABgAAAAD/+QAFAAAAWwAAAGAABwAGAAAAAP/5AAUAAABgAAAA
ZQAHAAYAAAAA//kABQAAAGUAAABqAAcABgAAAAD/+QAFAAAAagAAAG8ABwAGAAAAAf/8AAMAAABv
AAAAcQAEAAYAAAAA//wAAwACAHEAAAB0AAYABgAAAAD/+gAE//8AdAAAAHgABQAGAAAAAP/7AAT/
/gB4AAAAfAADAAYAAAAB//oABf//AHwAAACAAAUABgAAAAD/+gAFAAAAgAAAAIUABgAGAAAAAP/5
AAYAAQCFAAAAiwAIAAYAAP////oABgAAAIsAAACSAAYABgAA////+gAFAAAAkgAAAJgABgAGAAAA
AP/6AAUAAACYAAAAnQAGAAYAAP////oABQAAAJ0AAACjAAYABgAA////+gAFAAAAowAAAKkABgAG
AAD////6AAUAAACpAAAArwAGAAYAAAAA//oABQAAAK8AAAC0AAYABgAA////+gAGAAAAtAAAALsA
BgAGAAAAAP/6AAQAAAC7AAAAvwAGAAYAAP////oABQAAAL8AAADFAAYABgAA////+gAGAAAAxQAA
AMwABgAGAAD////6AAUAAADMAAAA0gAGAAYAAP////oABQAAANIAAADYAAYABgAA////+gAGAAAA
2AAAAN8ABgAGAAAAAP/6AAUAAADfAAAA5AAGAAYAAP////oABQAAAOQAAADqAAYABgAAAAD/+gAF
AAEA6gAAAO8ABwAGAAD////6AAYAAADvAAAA9gAGAAYAAAAA//oABQAAAPYAAAD7AAYABgAA////
+gAFAAAA+wAAAQEABgAGAAD////6AAYAAAEBAAABCAAGAAYAAP////oABgAAAQgAAAEPAAYABgAA
////+gAGAAABDwAAARYABgAGAAAAAP/6AAYAAAEWAAABHAAGAAYAAP////oABgAAARwAAAEjAAYA
BgAAAAD/+gAFAAABIwAAASgABgAGAAAAAf/5AAQAAQEoAAABKwAIAAYAAAAA//kABAABASsAAAEv
AAgABgAAAAH/+QAEAAEBLwAAATIACAAGAAAAAP/5AAX//AEyAAABNwADAAYAAAAAAAEABgACATcA
AAE9AAEABgAAAAH/+QAE//wBPQAAAUAAAwAGAAAAAP/7AAYAAAFAAAABRgAFAAYAAP////kABQAA
AUYAAAFMAAcABgAAAAD/+wAFAAABTAAAAVEABQAGAAAAAP/5AAYAAAFRAAABVwAHAAYAAAAA//sA
BQAAAVcAAAFcAAUABgAAAAD/+QAFAAABXAAAAWEABwAGAAAAAP/7AAYAAgFhAAABZwAHAAYAAP//
//kABQAAAWcAAAFtAAcABgAAAAD/+QAGAAABbQAAAXMABwAGAAAAAP/5AAQAAgFzAAABdwAJAAYA
AP////kABgAAAXcAAAF+AAcABgAAAAD/+QAGAAABfgAAAYQABwAGAAD////7AAUAAAGEAAABigAF
AAYAAP////sABQAAAYoAAAGQAAUABgAAAAD/+wAFAAABkAAAAZUABQAGAAD////7AAUAAgGVAAAB
mwAHAAYAAAAA//sABgACAZsAAAGhAAcABgAAAAD/+wAGAAABoQAAAacABQAGAAAAAP/7AAYAAAGn
AAABrQAFAAYAAAAA//kABgAAAa0AAAGzAAcABgAA////+wAGAAABswAAAboABQAGAAD////7AAUA
AAG6AAABwAAFAAYAAP////sABgAAAcAAAAHHAAUABgAAAAD/+wAGAAABxwAAAc0ABQAGAAD////7
AAYAAgHNAAAB1AAHAAYAAAAA//sABQAAAdQAAAHZAAUABgAAAAH/+QAFAAEB2QAAAd0ACAAGAAAA
Av/6AAMAAQHdAAAB3gAHAAYAAAAA//kABAABAd4AAAHiAAgABgAAAAD/+wAF//0B4gAAAecAAgAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAYAAAAB
//sAAwACAecAAAHpAAcABgAAAAD/+QAFAAEB6QAAAe4ACAAGAAAAAP/5AAYAAAHuAAAB9AAHAAYA
AAAA//oABf//AfQAAAH5AAUABgAAAAD/+QAGAAAB+QAAAf8ABwAGAAAAAv/5AAMAAgH/AAACAAAJ
AAYAAAAA//kABQABAgAAAAIFAAgABgAAAAH/+gAE//sCBQAAAggAAQAGAAAAAP/5AAYAAAIIAAAC
DgAHAAYAAAAB//kABf/+Ag4AAAISAAUABgAA////+wAGAAACEgAAAhkABQAGAAAAAP/7AAX//gIZ
AAACHgADAAYAAAAA//wABf/9Ah4AAAIjAAEABgAAAAD/+QAHAAACIwAAAioABwAGAAAAAP/6AAT/
+wIqAAACLgABAAYAAAAA//kABP/8Ai4AAAIyAAMABgAAAAD/+gAFAAACMgAAAjcABgAGAAAAAf/5
AAT//QI3AAACOgAEAAYAAAAB//kABP/9AjoAAAI9AAQABgAAAAL/+QAE//sCPQAAAj8AAgAGAAD/
///7AAYAAgI/AAACRgAHAAYAAAAA//kABgABAkYAAAJMAAgABgAAAAH//AAD//0CTAAAAk4AAQAG
AAAAAf//AAQAAgJOAAACUQADAAYAAAAB//kABP/9AlEAAAJUAAQABgAAAAH/+QAF//4CVAAAAlgA
BQAGAAD////7AAYAAAJYAAACXwAFAAYAAP////kABgAAAl8AAAJmAAcABgAA////+QAGAAACZgAA
Am0ABwAGAAD////5AAYAAAJtAAACdAAHAAYAAAAA//sABQACAnQAAAJ5AAcABgAA////9wAGAAAC
eQAAAoAACQAGAAD////3AAYAAAKAAAAChwAJAAYAAP////cABgAAAocAAAKOAAkABgAA////9wAG
AAACjgAAApUACQAGAAD////4AAYAAAKVAAACnAAIAAYAAP////cABgAAApwAAAKjAAkABgAA////
+gAGAAACowAAAqoABgAGAAAAAP/6AAUAAgKqAAACrwAIAAYAAP////cABQAAAq8AAAK1AAkABgAA
////9wAFAAACtQAAArsACQAGAAD////3AAUAAAK7AAACwQAJAAYAAP////gABQAAAsEAAALHAAgA
BgAAAAD/9wAEAAACxwAAAssACQAGAAAAAP/3AAQAAALLAAACzwAJAAYAAAAA//cABAAAAs8AAALT
AAkABgAAAAD/+AAEAAAC0wAAAtcACAAGAAD////6AAUAAALXAAAC3QAGAAYAAP////cABgAAAt0A
AALkAAkABgAAAAD/9wAFAAAC5AAAAukACQAGAAAAAP/3AAUAAALpAAAC7gAJAAYAAAAA//cABQAA
Au4AAALzAAkABgAAAAD/9wAFAAAC8wAAAvgACQAGAAAAAP/4AAUAAAL4AAAC/QAIAAYAAAAA//oA
Bf//Av0AAAMCAAUABgAA////+gAGAAADAgAAAwkABgAGAAD////3AAYAAAMJAAADEAAJAAYAAP//
//cABgAAAxAAAAMXAAkABgAA////9wAGAAADFwAAAx4ACQAGAAD////4AAYAAAAAAAoABwASAAYA
AP////cABgAAAAcACgAOABMABgAA////+gAFAAAADgAKABQAEAAGAAD////6AAYAAAAUAAoAGwAQ
AAYAAAAA//gABgAAABsACgAhABIABgAAAAD/+AAGAAAAIQAKACcAEgAGAAAAAP/4AAYAAAAnAAoA
LQASAAYAAAAA//gABgAAAC0ACgAzABIABgAAAAD/+QAGAAAAMwAKADkAEQAGAAAAAP/3AAYAAAA5
AAoAPwATAAYAAP////sABQAAAD8ACgBFAA8ABgAAAAD/+wAFAAIARQAKAEoAEQAGAAAAAP/4AAUA
AABKAAoATwASAAYAAAAA//gABQAAAE8ACgBUABIABgAAAAD/+AAFAAAAVAAKAFkAEgAGAAAAAP/5
AAUAAABZAAoAXgARAAYAAAAA//gABgAAAF4ACgBkABIABgAAAAD/+AAGAAAAZAAKAGoAEgAGAAAA
AP/4AAYAAABqAAoAcAASAAYAAAAA//kABgAAAHAACgB2ABEABgAAAAD/+AAFAAAAdgAKAHsAEgAG
AAD////4AAYAAAB7AAoAggASAAYAAAAA//gABQAAAIIACgCHABIABgAAAAD/+AAFAAAAhwAKAIwA
EgAGAAAAAP/4AAUAAACMAAoAkQASAAYAAAAA//gABQAAAJEACgCWABIABgAAAAD/+QAFAAAAlgAK
AJsAEQAGAAAAAP/6AAX//wCbAAoAoAAPAAYAAAAA//oABQABAKAACgClABEABgAA////+AAGAAAA
pQAKAKwAEgAGAAD////4AAYAAACsAAoAswASAAYAAP////gABgAAALMACgC6ABIABgAA////+QAG
AAAAugAKAMEAEQAGAAD////4AAYAAgDBAAoAyAAUAAYAAP////kABQACAMgACgDOABMABgAA////
+QAGAAIAzgAKANUAEw==
''')), Image.open(BytesIO(base64.decodestring(b'''
iVBORw0KGgoAAAANSUhEUgAAAx4AAAAUAQAAAAArMtZoAAAEwElEQVR4nABlAJr/AHVE4czCI/4u
Mc4b7vuds/xzjz5/3/7u/n9vMe7vnfH/9++vPn/xyf5zhxzjt8GHw8+2d83u8x27199/nxuQ6Od9
M43/5z2I+9n9ZtmDBwMQECDRQw/eQIQohJXxpBCNVE6QCCAAAAD//wBlAJr/AgALyj1t/wINwq0g
LeNZUworuN1cjTPIzrTX6ofHWeo3v336qPzfEwRmBnHTtf95/fglZK5N0PDgfRTslpGBvz7LFc4F
IUXBWQGjQ5MGCx34EDFPwXiY4YbYxavpnhHFrk14CDAAAAD//wBlAJr/AgKqRooH2gAgPeggvUAA
Bu2WfgPoAwzRAABAAAAAAACQgLz/3Uv4Gv+gX7BJgDeeGP6AAAD1NMDzKHD7ANWr3loYbxsAD791
NAADfcoIDyP44K/jv4Y63/Z+t98Ovt+ub4T48LAAAAD//wBlAJr/AuplMlADJAAAAGuAphWpqhMx
in0A/fRvAYBABPgBwBUgABBQ/sYAyv9g0bCHgOLoGAAAAAAAREAAwI7nr0ArYpow7aX8//9LaP/9
SjdavWA8ePHeBIKB//81/83ndznOaXx379wAAAD//wBlAJr/AqDxW+D3AABAAbUh/QMnbQag/gAY
AYDAAACgtgD/gOqAAAB5IA/8AAAk+n9w0AAA8AAAmFRJuPo27ciC0cD5oeW4E7KA/wD3ECMAn2tt
y8PgwH8AfAxFzC0JzeAMtratAsC/ffwAAAD//wBlAJr/BGKAyCAA4AAAAvgeYTAwHd1kmQF5chkG
ABoMIHcL5xVpTfQbUqzlAAAErwAQBgAAEOClA5D9il08AEh/tUzdCBsXkbgACED+woQg8Si9VeqY
lODCn7lmF6NhnAEYgAAA/NMIAAAAAAD//2JgjLZgVGBg5Pv/Tvpc8hwGBjYGJADjHDrAwPzAjv/H
/Wf3PzCwtzcwHmBgYGcwbZz8wHaCAQMDOwMDQ8MCBgYOC3W7mp+f0w+wHOYxO3OG+e376hsMZjk3
AAAAAP//YmCMY2A4wMAIN5e5gQETPD6AZisDAwMDgzSDAAPjByiHcQMDAwMDg1nOze1lByRu5/47
c4859311AYNZzg0AAAAA//9iYGDBYihOIIMuwIjGL39/fwffA8b//xv/P2BPtzzHwCBjUQAAAAD/
/yLFBrIBAAAA//9i1HhcwdhizX7u8NZNzyLbvT97bfrMf/QHI8evOwcSqGUJAAAA//9iYBB81iSw
pEE170Qrg5MIYydHqwdDQRMrAwcVrQAAAAD//2J4x7j9AAMDn8Q/BgYLBoaiAwwMjPdvMDBYM1Tv
oJodAAAAAP//Yqo/83+dxePWlxl3npsel9lvLfPcqlE9725C+acfVLMEAAAA//9i+s9gwCoaaGMR
evta/58PTEWzr21hufPjA8N+qlnBwAAAAAD//2JiWLci5v1+HmFXDqcnULE/MxgYGBj+f6CaJQAA
AAD//2Ji2FrkY3iYpYC5qDeGgeEMAwPDvwQBBoYvcTwOVLMEAAAA//9isDBgkP///0EOg9z35v//
Gc/eeW7BwPj5+QGZhANUswMAAAD//2JgqGBgYGBgqEMXlvhMPUsAAAAA//8iYDd1AAAAAP//AwDR
w7IkEbzhVQAAAABJRU5ErkJggg==
'''))))
    return f

# End of file


#
# The Python Imaging Library.
# $Id$
#
# standard image operations
#
# History:
# 2001-10-20 fl   Created
# 2001-10-23 fl   Added autocontrast operator
# 2001-12-18 fl   Added Kevin's fit operator
# 2004-03-14 fl   Fixed potential division by zero in equalize
# 2005-05-05 fl   Fixed equalize for low number of values
#
# Copyright (c) 2001-2004 by Secret Labs AB
# Copyright (c) 2001-2004 by Fredrik Lundh
#
# See the README file for information on usage and redistribution.
#

from PIL import Image
from PIL._util import isStringType
import operator
from functools import reduce


#
# helpers

def _border(border):
    if isinstance(border, tuple):
        if len(border) == 2:
            left, top = right, bottom = border
        elif len(border) == 4:
            left, top, right, bottom = border
    else:
        left = top = right = bottom = border
    return left, top, right, bottom


def _color(color, mode):
    if isStringType(color):
        from PIL import ImageColor
        color = ImageColor.getcolor(color, mode)
    return color


def _lut(image, lut):
    if image.mode == "P":
        # FIXME: apply to lookup table, not image data
        raise NotImplementedError("mode P support coming soon")
    elif image.mode in ("L", "RGB"):
        if image.mode == "RGB" and len(lut) == 256:
            lut = lut + lut + lut
        return image.point(lut)
    else:
        raise IOError("not supported for this image mode")

#
# actions


def autocontrast(image, cutoff=0, ignore=None):
    """
    Maximize (normalize) image contrast. This function calculates a
    histogram of the input image, removes **cutoff** percent of the
    lightest and darkest pixels from the histogram, and remaps the image
    so that the darkest pixel becomes black (0), and the lightest
    becomes white (255).

    :param image: The image to process.
    :param cutoff: How many percent to cut off from the histogram.
    :param ignore: The background pixel value (use None for no background).
    :return: An image.
    """
    histogram = image.histogram()
    lut = []
    for layer in range(0, len(histogram), 256):
        h = histogram[layer:layer+256]
        if ignore is not None:
            # get rid of outliers
            try:
                h[ignore] = 0
            except TypeError:
                # assume sequence
                for ix in ignore:
                    h[ix] = 0
        if cutoff:
            # cut off pixels from both ends of the histogram
            # get number of pixels
            n = 0
            for ix in range(256):
                n = n + h[ix]
            # remove cutoff% pixels from the low end
            cut = n * cutoff // 100
            for lo in range(256):
                if cut > h[lo]:
                    cut = cut - h[lo]
                    h[lo] = 0
                else:
                    h[lo] -= cut
                    cut = 0
                if cut <= 0:
                    break
            # remove cutoff% samples from the hi end
            cut = n * cutoff // 100
            for hi in range(255, -1, -1):
                if cut > h[hi]:
                    cut = cut - h[hi]
                    h[hi] = 0
                else:
                    h[hi] -= cut
                    cut = 0
                if cut <= 0:
                    break
        # find lowest/highest samples after preprocessing
        for lo in range(256):
            if h[lo]:
                break
        for hi in range(255, -1, -1):
            if h[hi]:
                break
        if hi <= lo:
            # don't bother
            lut.extend(list(range(256)))
        else:
            scale = 255.0 / (hi - lo)
            offset = -lo * scale
            for ix in range(256):
                ix = int(ix * scale + offset)
                if ix < 0:
                    ix = 0
                elif ix > 255:
                    ix = 255
                lut.append(ix)
    return _lut(image, lut)


def colorize(image, black, white):
    """
    Colorize grayscale image.  The **black** and **white**
    arguments should be RGB tuples; this function calculates a color
    wedge mapping all black pixels in the source image to the first
    color, and all white pixels to the second color.

    :param image: The image to colorize.
    :param black: The color to use for black input pixels.
    :param white: The color to use for white input pixels.
    :return: An image.
    """
    assert image.mode == "L"
    black = _color(black, "RGB")
    white = _color(white, "RGB")
    red = []
    green = []
    blue = []
    for i in range(256):
        red.append(black[0]+i*(white[0]-black[0])//255)
        green.append(black[1]+i*(white[1]-black[1])//255)
        blue.append(black[2]+i*(white[2]-black[2])//255)
    image = image.convert("RGB")
    return _lut(image, red + green + blue)


def crop(image, border=0):
    """
    Remove border from image.  The same amount of pixels are removed
    from all four sides.  This function works on all image modes.

    .. seealso:: :py:meth:`~PIL.Image.Image.crop`

    :param image: The image to crop.
    :param border: The number of pixels to remove.
    :return: An image.
    """
    left, top, right, bottom = _border(border)
    return image.crop(
        (left, top, image.size[0]-right, image.size[1]-bottom)
        )


def deform(image, deformer, resample=Image.BILINEAR):
    """
    Deform the image.

    :param image: The image to deform.
    :param deformer: A deformer object.  Any object that implements a
                    **getmesh** method can be used.
    :param resample: What resampling filter to use.
    :return: An image.
    """
    return image.transform(
        image.size, Image.MESH, deformer.getmesh(image), resample
        )


def equalize(image, mask=None):
    """
    Equalize the image histogram. This function applies a non-linear
    mapping to the input image, in order to create a uniform
    distribution of grayscale values in the output image.

    :param image: The image to equalize.
    :param mask: An optional mask.  If given, only the pixels selected by
                 the mask are included in the analysis.
    :return: An image.
    """
    if image.mode == "P":
        image = image.convert("RGB")
    h = image.histogram(mask)
    lut = []
    for b in range(0, len(h), 256):
        histo = [_f for _f in h[b:b+256] if _f]
        if len(histo) <= 1:
            lut.extend(list(range(256)))
        else:
            step = (reduce(operator.add, histo) - histo[-1]) // 255
            if not step:
                lut.extend(list(range(256)))
            else:
                n = step // 2
                for i in range(256):
                    lut.append(n // step)
                    n = n + h[i+b]
    return _lut(image, lut)


def expand(image, border=0, fill=0):
    """
    Add border to the image

    :param image: The image to expand.
    :param border: Border width, in pixels.
    :param fill: Pixel fill value (a color value).  Default is 0 (black).
    :return: An image.
    """
    "Add border to image"
    left, top, right, bottom = _border(border)
    width = left + image.size[0] + right
    height = top + image.size[1] + bottom
    out = Image.new(image.mode, (width, height), _color(fill, image.mode))
    out.paste(image, (left, top))
    return out


def fit(image, size, method=Image.NEAREST, bleed=0.0, centering=(0.5, 0.5)):
    """
    Returns a sized and cropped version of the image, cropped to the
    requested aspect ratio and size.

    This function was contributed by Kevin Cazabon.

    :param size: The requested output size in pixels, given as a
                 (width, height) tuple.
    :param method: What resampling method to use. Default is
                   :py:attr:`PIL.Image.NEAREST`.
    :param bleed: Remove a border around the outside of the image (from all
                  four edges. The value is a decimal percentage (use 0.01 for
                  one percent). The default value is 0 (no border).
    :param centering: Control the cropping position.  Use (0.5, 0.5) for
                      center cropping (e.g. if cropping the width, take 50% off
                      of the left side, and therefore 50% off the right side).
                      (0.0, 0.0) will crop from the top left corner (i.e. if
                      cropping the width, take all of the crop off of the right
                      side, and if cropping the height, take all of it off the
                      bottom).  (1.0, 0.0) will crop from the bottom left
                      corner, etc. (i.e. if cropping the width, take all of the
                      crop off the left side, and if cropping the height take
                      none from the top, and therefore all off the bottom).
    :return: An image.
    """

    # by Kevin Cazabon, Feb 17/2000
    # kevin@cazabon.com
    # http://www.cazabon.com

    # ensure inputs are valid
    if not isinstance(centering, list):
        centering = [centering[0], centering[1]]

    if centering[0] > 1.0 or centering[0] < 0.0:
        centering[0] = 0.50
    if centering[1] > 1.0 or centering[1] < 0.0:
        centering[1] = 0.50

    if bleed > 0.49999 or bleed < 0.0:
        bleed = 0.0

    # calculate the area to use for resizing and cropping, subtracting
    # the 'bleed' around the edges

    # number of pixels to trim off on Top and Bottom, Left and Right
    bleedPixels = (
        int((float(bleed) * float(image.size[0])) + 0.5),
        int((float(bleed) * float(image.size[1])) + 0.5)
        )

    liveArea = (0, 0, image.size[0], image.size[1])
    if bleed > 0.0:
        liveArea = (
            bleedPixels[0], bleedPixels[1], image.size[0] - bleedPixels[0] - 1,
            image.size[1] - bleedPixels[1] - 1
            )

    liveSize = (liveArea[2] - liveArea[0], liveArea[3] - liveArea[1])

    # calculate the aspect ratio of the liveArea
    liveAreaAspectRatio = float(liveSize[0])/float(liveSize[1])

    # calculate the aspect ratio of the output image
    aspectRatio = float(size[0]) / float(size[1])

    # figure out if the sides or top/bottom will be cropped off
    if liveAreaAspectRatio >= aspectRatio:
        # liveArea is wider than what's needed, crop the sides
        cropWidth = int((aspectRatio * float(liveSize[1])) + 0.5)
        cropHeight = liveSize[1]
    else:
        # liveArea is taller than what's needed, crop the top and bottom
        cropWidth = liveSize[0]
        cropHeight = int((float(liveSize[0])/aspectRatio) + 0.5)

    # make the crop
    leftSide = int(liveArea[0] + (float(liveSize[0]-cropWidth) * centering[0]))
    if leftSide < 0:
        leftSide = 0
    topSide = int(liveArea[1] + (float(liveSize[1]-cropHeight) * centering[1]))
    if topSide < 0:
        topSide = 0

    out = image.crop(
        (leftSide, topSide, leftSide + cropWidth, topSide + cropHeight)
        )

    # resize the image and return it
    return out.resize(size, method)


def flip(image):
    """
    Flip the image vertically (top to bottom).

    :param image: The image to flip.
    :return: An image.
    """
    return image.transpose(Image.FLIP_TOP_BOTTOM)


def grayscale(image):
    """
    Convert the image to grayscale.

    :param image: The image to convert.
    :return: An image.
    """
    return image.convert("L")


def invert(image):
    """
    Invert (negate) the image.

    :param image: The image to invert.
    :return: An image.
    """
    lut = []
    for i in range(256):
        lut.append(255-i)
    return _lut(image, lut)


def mirror(image):
    """
    Flip image horizontally (left to right).

    :param image: The image to mirror.
    :return: An image.
    """
    return image.transpose(Image.FLIP_LEFT_RIGHT)


def posterize(image, bits):
    """
    Reduce the number of bits for each color channel.

    :param image: The image to posterize.
    :param bits: The number of bits to keep for each channel (1-8).
    :return: An image.
    """
    lut = []
    mask = ~(2**(8-bits)-1)
    for i in range(256):
        lut.append(i & mask)
    return _lut(image, lut)


def solarize(image, threshold=128):
    """
    Invert all pixel values above a threshold.

    :param image: The image to solarize.
    :param threshold: All pixels above this greyscale level are inverted.
    :return: An image.
    """
    lut = []
    for i in range(256):
        if i < threshold:
            lut.append(i)
        else:
            lut.append(255-i)
    return _lut(image, lut)


# --------------------------------------------------------------------
# PIL USM components, from Kevin Cazabon.

def gaussian_blur(im, radius=None):
    """ PIL_usm.gblur(im, [radius])"""

    if radius is None:
        radius = 5.0

    im.load()

    return im.im.gaussian_blur(radius)

gblur = gaussian_blur


def unsharp_mask(im, radius=None, percent=None, threshold=None):
    """ PIL_usm.usm(im, [radius, percent, threshold])"""

    if radius is None:
        radius = 5.0
    if percent is None:
        percent = 150
    if threshold is None:
        threshold = 3

    im.load()

    return im.im.unsharp_mask(radius, percent, threshold)

usm = unsharp_mask
