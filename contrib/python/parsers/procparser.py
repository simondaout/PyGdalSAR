################################################################################
# 
# NSBAS - New Small Baseline Chain
# 
################################################################################
# Author        : Matthieu Volat (ISTerre)
################################################################################
"""Proc file parser.

This module propose a very similar interface to the ConfigParser
module in python to read, manipulate and write proc files.

Aside from the fact they handle different syntaxes, there are a few 
differences on how the ProcParser behave compared to the ConfigParser 
objects:
* as there is no section in proc files, the object exposes directly 
  its options, and there is never a section argument in functions.
* it is possible that an option is repeated with different values 
  in a proc file. In that case, the parser will set the value of the option 
  to a list of the sucessive values found in the file with that option.
* there is no substitution in the options values

Example:
>>> p = procparser.ProcParser()
>>> p.read("/foo/bar.proc")
>>> print p.get("OrbitType")
HDR
>>> print len(p.get("SarDir"))
12

class:

ProcParser -- responsible for parsing a list of proc files and managing
                  the parsed database.
"""

import sys

if sys.version_info[0] == 2:
    from io import open
if sys.version_info[0] >= 3:
    from collections.abc import MutableMapping
    from configparser import ParsingError, Error
    from collections import OrderedDict
elif sys.version_info[0] >= 2:
    from collections import MutableMapping
    from ConfigParser import ParsingError, Error
    if sys.version_info[1] >= 7:
        from collections import OrderedDict
    elif sys.version_info[1] >= 4:
        from ordereddict import OrderedDict
import itertools
import re


class NoOptionError(Error):
    def __init__(self, option):
        Error.__init__(self, "No option %r" % (option))
        self.option = option
        self.args = (option)

class ProcParser(MutableMapping):
    BOOLEAN_STATES = {'1': True,  'yes': True, 'true': True,   'on': True,
                      '0': False, 'no': False, 'false': False, 'off': False}

    LINE_EXPR = r"(?P<option>.*?)\s*(?P<vi>=)\s*(?P<value>.*)$"
    _re = re.compile(LINE_EXPR.format(delim="="), re.VERBOSE)

    def __init__(self, defaults=None):
        self._defaults = OrderedDict()
        self._options = OrderedDict()
        if defaults:
            for key, value in defaults.items():
                self._defaults[key] = value

    def __delitem__(self, key):
        if not key in self._options and not key in self._defaults:
            raise KeyError(key)
        if key in self._options:
            del self._options[key]
        if key in self._defaults:
            del self._defaults[key]

    def __getitem__(self, key):
        if key in self._options:
            return self._options[key]
        elif key in self._defaults:
            return self._defaults[key]
        else:
            raise KeyError(key)

    def __iter__(self):
        opts = self._options.copy()
        opts.update(self._defaults)
        return itertools.chain(opts.keys())

    def __len__(self):
        opts = self._options.copy()
        opts.update(self._defaults)
        return opts.len()
        
    def __setitem__(self, key, value):
        self._options[key] = value

    def _handle_error(self, exc, fpname, lineno, line):
        if not exc:
            exc = ParsingError(fpname)
            exc.append(lineno, repr(line))
        return exc

    def _read(self, fp, fpname):
        e = None
        for lineno, line in enumerate(fp, start=1):
            # skip full comment
            if line.strip().startswith("#"):
                continue
            # strip inline comment
            index = line.find("#", 1)
            if index == 0 or (index > 0 and line[index-1].isspace()):
                line = line[:index]
            # strip remaining white spaces
            line = line.strip()
            # skip empty
            if (not line): 
                continue
            # match 
            m = self._re.match(line)
            if m:
                name, val = m.group("option", "value")
                if not name in self._options:
                    self._options[name] = val
                elif isinstance(self._options[name], list):
                    self._options[name].append(val)
                else:
                    self._options[name] = [self._options[name], val]
            else:
                # Handle things like the RawConfigParser
                e = self._handle_error(e, fpname, lineno, line)
        if e:
            raise e

    def has_option(self, option):
        return (option in self._options or option in self._defaults)

    def options(self):
        opts = self._options.copy()
        opts.update(self._defaults)
        return opts.keys()

    def read(self, filenames, encoding=None):
        read_ok = []
        if isinstance(filenames, str) or isinstance(filenames, unicode):
            filenames = [filenames]
        for filename in filenames:
            try:
                with open(filename, encoding=encoding) as fp:
                    self._read(fp, filename)
            except IOError:
                continue
            read_ok.append(filename)
        return read_ok

    def read_file(self, f, source=None):
        if source is None:
            try:
                source = f.name
            except AttributeError:
                source = "<???>"
        self._read(f, source)

    def get(self, option):
        if option in self:
            return self[option]
        else:
            raise NoOptionError()

    def getint(self, option):
        return int(self.get(option))

    def getfloat(self, option):
        return float(self.get(option))

    def getboolean(self, option):
        value = self.get(option)
        if value.lower() not in self.BOOLEAN_STATES:
            raise ValueError('Not a boolean: %s' % value)
        return self.BOOLEAN_STATES[value.lower()]

    def items(self):
        d = self._defaults.copy()
        d.update(self._options)
        if "__name__" in d:
            del d["__name__"]
        return d.items()

    def remove_option(self, option):
        existed = option in self._options
        if existed:
            del self._options[option]
        return existed

    def set(self, option, value):
        self._options[option] = value

    def write(self, f, space_around_delimiters=True):
        delimiter = " = " if space_around_delimiters else "="
        for key, value in self._options.items():
            if isinstance(value, list):
                for subvalue in value:
                    f.write("{}{}{}\n".format(key, delimiter, subvalue))
            else:
                if not value:
                    value = ""
                f.write("{}{}{}\n".format(key, delimiter, value))
