from functools import wraps

import numpy as np

from supporting_tools import PipeArgError, pad_skipping_iterator, doublewrap, list_not_string, \
    _3d_stack_2d_filter, _2d_stack_2d_filter


@doublewrap
def generator_wrapper(f, in_dims=(3,), out_dims=None):

    if out_dims is None:
            out_dims = in_dims

    @wraps(f)
    def inner_wrapper(*args, **kwargs):
        """
        converts a function to accepting a generator of named dicts and adds channel selection logic
        """

        iterator = args[0]
        args = args[1:]

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
                    print 'Input %s will be overwritten by function %s' % (in_chan[0], f.__name__)
                    out_chan = in_chan
                else:
                    print f.__name__
                    print in_chan, in_dims
                    raise PipeArgError('Please provide out_channel argument')

            if len(in_chan) != len(in_dims):
                print f.__name__
                print in_chan, in_dims
                print len(in_chan), len(in_dims)
                raise PipeArgError('%s inbound channels are piped, function allows %s' %
                                   (len(in_chan), len(in_dims)))

            if len(out_chan) != len(out_dims):
                print f.__name__
                print out_chan, out_dims
                print len(out_chan), len(out_dims)
                raise PipeArgError('%s outbound channels are piped, function allows %s' %
                                   (len(out_chan), len(out_dims)))
            # end in/out channel logic

            for name_space in iterator:
                # start args prepare
                args_puck = []

                for i, chan in enumerate(in_chan):
                    if in_dims[i] and len(name_space[chan].shape) != in_dims[i]:
                        print f.__name__
                        print chan, len(name_space[chan].shape), in_dims[i]
                        raise PipeArgError('Mismatched inbound channel dimension for channel. %s is of dim %s, expected %s' %
                                           (chan, len(name_space[chan].shape), in_dims[i]))
                    args_puck.append(name_space[chan])

                local_args = tuple(args_puck) + args
                # end args prepare
                return_puck = f(*local_args, **kwargs)

                if return_puck is None and out_chan[0] == '_':
                    yield name_space  # unlike return, yield is probably non-blocking....

                else:
                    # start output prepare
                    if not isinstance(return_puck, tuple):
                        return_puck = (return_puck, )

                    for i, chan in enumerate(out_chan):
                        if out_dims[i] and len(return_puck[i].shape) != out_dims[i]:
                            print f.__name__
                            print chan
                            raise PipeArgError('Mismatched outgoing channel dimension for channel. %s is of dim %s, expected %s' %
                                               (chan, len(return_puck[i].shape), out_dims[i]))
                        if chan != '_':
                            name_space[chan] = return_puck[i]
                    # end output prepare

                    yield name_space

        else:
            for name_space in iterator:

                local_args = (name_space,) + args
                name_space = f(*local_args, **kwargs)
                yield name_space

    return inner_wrapper


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
                    base_chan = _2d_stack_2d_filter(primary_namespace[chan], local_mask)
                elif len(primary_namespace[chan].shape) == 3:
                    base_chan = _3d_stack_2d_filter(primary_namespace[chan], local_mask)
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

# original: in dims below = 2
# maybe not set dims and covert to 2 dims


# To try: multiscale percentile edge finding.


