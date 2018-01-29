from functools import wraps

import numpy as np

from imagepipe.raw_functions import f_3d_stack_2d_filter, f_2d_stack_2d_filter
from imagepipe.tools.helpers import safe_dir_create, PipeArgError, list_not_string

safe_dir_create('verification')  # guarantees existence, although might not be the best location to do it


def doublewrap(f):
    """
    a decorator decorator, allowing the decorator to be used as:
    @decorator(with, arguments, and=kwargs)
    or
    @decorator

    credits: http://stackoverflow.com/questions/653368/how-to-create-a-python-decorator-that-can-be-used-either-with-or-without-paramet

    """
    @wraps(f)
    def new_dec(*args, **kwargs):
        if len(args) == 1 and len(kwargs) == 0 and callable(args[0]):
            # actual decorated function
            return f(args[0])
        else:
            # decorator arguments
            return lambda realf: f(realf, *args, **kwargs)

    return new_dec


@doublewrap
def generator_wrapper(my_function, in_dims=(3,), out_dims=None):

    if out_dims is None:
            out_dims = in_dims

    @wraps(my_function)
    def inner_wrapper(*args, **kwargs):
        """
        converts a function to accepting a generator of named dicts and adds channel selection logic
        """

        iterator = args[0]
        args = args[1:]
        print my_function.__name__

        if 'in_channel' in kwargs:
            #Start in/out channel logic

            in_chan = kwargs['in_channel']  # Multiple arguments
            del kwargs['in_channel']

            if not list_not_string(in_chan):  # convert string to list of strings
                in_chan = [in_chan]

            if 'out_channel' in kwargs:
                out_chan = kwargs['out_channel']  # output explicitely provided
                del kwargs['out_channel']
                if not list_not_string(out_chan):  # convert string to list of strings
                    out_chan = [out_chan]

            else:  # implicit output, bound to in_channel only if a single input is provided
                if len(in_chan) == 1:
                    print 'Input %s will be overwritten by function %s' % (in_chan[0], my_function.__name__)
                    out_chan = in_chan
                else:
                    print my_function.__name__
                    print in_chan, in_dims
                    raise PipeArgError('Please provide out_channel argument')

            if len(in_chan) != len(in_dims):
                print my_function.__name__
                print in_chan, in_dims
                print len(in_chan), len(in_dims)
                raise PipeArgError('%s inbound channels are piped, function allows %s' %
                                   (len(in_chan), len(in_dims)))

            if len(out_chan) != len(out_dims):
                print my_function.__name__
                print out_chan, out_dims
                print len(out_chan), len(out_dims)
                raise PipeArgError('%s outbound channels are piped, function allows %s' %
                                   (len(out_chan), len(out_dims)))
            # end in/out channel logic

            for name_space in iterator:
                # start args prepare
                print "pullin function %s " % my_function.__name__ 
                print args, kwargs
                args_puck = []

                for i, chan in enumerate(in_chan):
                    if in_dims[i] and len(name_space[chan].shape) != in_dims[i]:
                        print my_function.__name__
                        print chan, len(name_space[chan].shape), in_dims[i]
                        raise PipeArgError('Mismatched inbound channel dimension for channel. %s is of dim %s, expected %s' %
                                           (chan, len(name_space[chan].shape), in_dims[i]))
                    args_puck.append(name_space[chan])


                local_args = tuple(args_puck) + args
                # end args prepare
                print "local args ready"
                print my_function.__name__
                return_puck = my_function(*local_args, **kwargs)
                print "return puck ready"

                if return_puck is None and out_chan[0] == '_':
                    print my_function.__name__, "yields"
                    yield name_space  # unlike return, yield is probably non-blocking....

                else:
                    # start output prepare
                    if not isinstance(return_puck, tuple):
                        return_puck = (return_puck, )

                    for i, chan in enumerate(out_chan):
                        if out_dims[i] and len(return_puck[i].shape) != out_dims[i]:
                            print my_function.__name__
                            print chan
                            raise PipeArgError('Mismatched outgoing channel dimension for channel. %s is of dim %s, expected %s' %
                                               (chan, len(return_puck[i].shape), out_dims[i]))
                        if chan != '_':
                            name_space[chan] = return_puck[i]
                    # end output prepare
                    print my_function.__name__, "yields"
                    yield name_space

        else:
            for name_space in iterator:

                local_args = (name_space,) + args
                name_space = my_function(*local_args, **kwargs)
                print my_function.__name__, "yields"
                yield name_space

    return inner_wrapper


def pad_skipping_iterator(secondary_namespace):
    for key, value in secondary_namespace.iteritems():
        if key != '_pad':
            yield value


def splitter(outer_generator, to, sources, mask):
    """
    Creates a secondary namespace by using mask as a pad to conserve only certain segments in sources

    :param outer_generator:
    :param to:
    :param sources:
    :param mask:
    :return:
    """
    for primary_namespace in outer_generator:
        primary_namespace[to] = {}
        unique_vals = np.unique(primary_namespace[mask])
        unique_vals = unique_vals[unique_vals > 0]

        primary_namespace[to]['_pad'] = (unique_vals, primary_namespace[mask])  # used to rebuild padded images

        for val in unique_vals:
            secondary_namespace = {}
            primary_namespace[to][val] = secondary_namespace

            for chan in sources:
                local_mask = primary_namespace[mask] == val
                if len(primary_namespace[chan].shape) == 2:
                    base_chan = f_2d_stack_2d_filter(primary_namespace[chan], local_mask)
                elif len(primary_namespace[chan].shape) == 3:
                    base_chan = f_3d_stack_2d_filter(primary_namespace[chan], local_mask)
                else:
                    raise PipeArgError('masking impossible: dims not match, base channel %s is of dim %s' %
                                       (chan, len(primary_namespace[chan].shape)))

                secondary_namespace[chan] = base_chan
        yield primary_namespace


def for_each(outer_generator, embedded_transformer, inside, **kwargs):

    for primary_namespace in outer_generator:

        secondary_generator = embedded_transformer(pad_skipping_iterator(primary_namespace[inside]), **kwargs)
        for i, _ in enumerate(secondary_generator):  # forces secondary generator to evaluate
            pass
        yield primary_namespace


def paint_from_mask(outer_generator, based_on, in_anchor, out_channel=None):

    if out_channel is None:
        out_channel = in_anchor

    for primary_namespace in outer_generator:
        secondary_namespace = primary_namespace[based_on]
        mask = secondary_namespace['_pad'][1]
        mask_values = secondary_namespace['_pad'][0]
        accumulator = np.zeros_like(mask)

        for i, unique_value in enumerate(mask_values):
            if i == 0:
                accumulator = accumulator.astype(secondary_namespace[unique_value][in_anchor].dtype)
            accumulator[mask == unique_value] = secondary_namespace[unique_value][in_anchor]

        primary_namespace[out_channel] = accumulator
        yield primary_namespace


def tile_from_mask(outer_generator, based_on, in_anchor, out_channel=None):

    if out_channel is None:
        out_channel = in_anchor

    for primary_namespace in outer_generator:

        secondary_namespace = primary_namespace[based_on]
        mask = secondary_namespace['_pad'][1]
        mask_values = secondary_namespace['_pad'][0]
        accumulator = np.zeros_like(mask)

        for i, unique_value in enumerate(mask_values):

            accumulator = accumulator.astype(secondary_namespace[unique_value][in_anchor].dtype)
            accumulator[mask == unique_value] = secondary_namespace[unique_value][in_anchor][mask == unique_value]

        primary_namespace[out_channel] = accumulator
        yield primary_namespace


