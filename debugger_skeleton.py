import logging
from abc import abstractmethod, abstractproperty
import dill


class Debugger(object):

    def __init__(self, logger_name='Default Debug Logger', debug_file='debug_log.log'):
        self.logger = logging.getLogger(logger_name)
        self.logger.setLevel(logging.DEBUG)

        fh = logging.FileHandler(debug_file)
        fh.setLevel(logging.DEBUG)

        # create console handler with a higher log level
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)

        # create formatter and add it to the handlers
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)

        # add the handlers to the logger
        self.logger.addHandler(fh)
        self.logger.addHandler(ch)

    def pickle(self):
        return dill.dumps(self)

    def unpickle(self, dill_object):
        self.__dict__.update(dill.loads(dill_object).__dict__)

    def render(self, show=False):
        for _property in [a for a in dir(self)
                          if not a.startswith('__') and not callable(getattr(self, a))]:

            if isinstance(_property, DebugFrame):
                if _property.is_complete():
                    _property.render()


class DebugFrame(object):

    _fields = []  # contains one tuple for each
    # check if all the data_store fields have been set

    def __init__(self):
        for name in self._fields:
            setattr(self, name, None)

    @abstractmethod
    def render(self, show=False):
        pass

    def is_complete(self):
        for name in self._fields:
            if name is None:
                return False
        else:
            return True
