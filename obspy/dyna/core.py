# -*- coding: utf-8 -*-
"""
DYNA and ITACA bindings to ObsPy core module.

:copyright:
    The ObsPy Development Team (devs@obspy.org)
:license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA

from StringIO import StringIO
from obspy.core import Stream, Trace, UTCDateTime, Stats
import numpy as np
import re


def isDYNA(filename):
    """
    Checks whether a file is a DYNA ASCII file or not.

    :type filename: str
    :param filename: Name of the DYNA file to be checked.
    :rtype: bool
    :return: ``True`` if a DYNA ASCII file.

    .. rubric:: Example

    >>> isDYNA("/path/to/IT.AQG..HNE.D.20090330.215717.X.ACC.ASC")
    True
    """
    # first eleven chars should contain 'EVENT_NAME:'
    try:
        temp = open(filename, 'rt').read(11)
    except:
        return False
    if temp != 'EVENT_NAME:':
        return False
    # first 6 chars of line 54 should contain 'USER4:'
    f = open(filename, 'rt')
    for i in xrange(55):
        line = f.readline()
        if i == 54 and line[:6] != 'USER4:':
            return False
    return True


def isITACA(filename):
    """
    Checks whether a file is a ITACA ASCII file or not.

    :type filename: str
    :param filename: Name of the ITACA file to be checked.
    :rtype: bool
    :return: ``True`` if a ITACA ASCII file.

    .. rubric:: Example

    >>> isITACA("/path/to/20090330_215717ITDPC_AQG__WEX.DAT")
    True
    """
    # first eleven chars should contain 'EVENT_NAME:'
    try:
        temp = open(filename, 'rt').read(11)
    except:
        return False
    if temp != 'EVENT_NAME:':
        return False
    # first 12 chars of line 8 should contain 'MAGNITUDE_S:'
    f = open(filename, 'rt')
    for i in xrange(8):
        line = f.readline()
        if i == 7 and line[:12] != 'MAGNITUDE_S:':
            return False
    # filename must match the following regexp:
    # ^\d{8}_\d{6}\w{5}_\w{5}\w{3}.\w{3}$
    fname_regexp = '.*\d{8}_\d{6}\w{5}_\w{5}\w{3}.\w{3}$'
    if not re.match(fname_regexp, filename):
        print("filename does not match ITACA format")
        return False
    return True


def readDYNA(filename, headonly=False, **kwargs):  # @UnusedVariable
    """
    Reads a DYNA ASCII file and returns an ObsPy Stream object.

    .. warning::
        This function should NOT be called directly, it registers via the
        ObsPy :func:`~obspy.core.stream.read` function, call this instead.

    :type filename: str
    :param filename: ASCII file to be read.
    :type headonly: bool, optional
    :param headonly: If set to True, read only the head. This is most useful
        for scanning available data in huge (temporary) data sets.
    :rtype: :class:`~obspy.core.stream.Stream`
    :return: A ObsPy Stream object.

    .. rubric:: Example

    >>> from obspy.core import read
    >>> st = read("/path/to/IT.AQG..HNE.D.20090330.215717.X.ACC.ASC")
    >>> st  # doctest: +ELLIPSIS
    <obspy.core.stream.Stream object at 0x...>
    >>> print(st)  # doctest: +ELLIPSIS
    1 Trace(s) in Stream:
    IT.AQG..E | 2009-03-30T21:56:50.967000Z - 2009-03-30T21:57:55.962000Z | 200.0 Hz, 13000 samples
    """
    headers = {}
    data = StringIO()

    # read file
    fh = open(filename, 'rt')
    for i in xrange(55):
        key, value = fh.readline().strip().split(':')
        headers[key.strip()] = value.strip()

    # create ObsPy stream object
    stream = Stream()
    header = Stats()
    header['dyna'] = {}

    header['network'] = headers['NETWORK']
    header['station'] = headers['STATION_CODE']
    header['location'] = headers['LOCATION']
    header['channel'] = headers['STREAM'][2:3]
    try:
        # use toUTCDateTime to convert from DYNA format
        header['starttime'] = toUTCDateTime(
            headers['DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS'])
    except:
        header['starttime'] = toUTCDateTime('19700101_000000')
    header['sampling_rate'] = 1 / float(headers['SAMPLING_INTERVAL_S'])
    header['delta'] = float(headers['SAMPLING_INTERVAL_S'])
    header['npts'] = int(headers['NDATA'])
    header['calib'] = 1  # not in file header

    # DYNA dict float data
    header['dyna']['EVENT_LATITUDE_DEGREE'] = strtofloat(
        headers['EVENT_LATITUDE_DEGREE'])
    header['dyna']['EVENT_LONGITUDE_DEGREE'] = strtofloat(
        headers['EVENT_LONGITUDE_DEGREE'])
    header['dyna']['EVENT_DEPTH_KM'] = strtofloat(headers['EVENT_DEPTH_KM'])
    header['dyna']['HYPOCENTER_REFERENCE'] = headers['HYPOCENTER_REFERENCE']
    header['dyna']['MAGNITUDE_W'] = strtofloat(headers['MAGNITUDE_W'])
    header['dyna']['MAGNITUDE_L'] = strtofloat(headers['MAGNITUDE_L'])
    header['dyna']['STATION_LATITUDE_DEGREE'] = strtofloat(
        headers['STATION_LATITUDE_DEGREE'])
    header['dyna']['STATION_LONGITUDE_DEGREE'] = strtofloat(
        headers['STATION_LONGITUDE_DEGREE'])
    header['dyna']['VS30_M_S'] = strtofloat(headers['VS30_M/S'])
    header['dyna']['EPICENTRAL_DISTANCE_KM'] = strtofloat(
        headers['EPICENTRAL_DISTANCE_KM'])
    header['dyna']['EARTHQUAKE_BACKAZIMUTH_DEGREE'] = strtofloat(
        headers['EARTHQUAKE_BACKAZIMUTH_DEGREE'])
    header['dyna']['DURATION_S'] = strtofloat(headers['DURATION_S'])
    header['dyna']['INSTRUMENTAL_FREQUENCY_HZ'] = strtofloat(
        headers['INSTRUMENTAL_FREQUENCY_HZ'])
    header['dyna']['INSTRUMENTAL_DAMPING'] = strtofloat(
        headers['INSTRUMENTAL_DAMPING'])
    header['dyna']['FULL_SCALE_G'] = strtofloat(headers['FULL_SCALE_G'])

    # data type is acceleration
    if headers['DATA_TYPE'] == "ACCELERATION" \
            or headers['DATA_TYPE'] == "ACCELERATION RESPONSE SPECTRUM":
        header['dyna']['PGA_CM_S_2'] = strtofloat(headers['PGA_CM/S^2'])
        header['dyna']['TIME_PGA_S'] = strtofloat(headers['TIME_PGA_S'])
    # data type is velocity
    if headers['DATA_TYPE'] == "VELOCITY" \
            or headers['DATA_TYPE'] == "PSEUDO-VELOCITY RESPONSE SPECTRUM":
        header['dyna']['PGV_CM_S'] = strtofloat(headers['PGV_CM/S'])
        header['dyna']['TIME_PGV_S'] = strtofloat(headers['TIME_PGV_S'])
    # data type is displacement
    if headers['DATA_TYPE'] == "DISPLACEMENT" \
            or headers['DATA_TYPE'] == "DISPLACEMENT RESPONSE SPECTRUM":
        header['dyna']['PGD_CM'] = strtofloat(headers['PGD_CM'])
        header['dyna']['TIME_PGD_S'] = strtofloat(headers['TIME_PGD_S'])

    header['dyna']['LOW_CUT_FREQUENCY_HZ'] = strtofloat(
        headers['LOW_CUT_FREQUENCY_HZ'])
    header['dyna']['HIGH_CUT_FREQUENCY_HZ'] = strtofloat(
        headers['HIGH_CUT_FREQUENCY_HZ'])

    # DYNA dict int data
    header['dyna']['STATION_ELEVATION_M'] = strtoint(
        headers['STATION_ELEVATION_M'])
    header['dyna']['N_BIT_DIGITAL_CONVERTER'] = strtoint(
        headers['N_BIT_DIGITAL_CONVERTER'])
    header['dyna']['FILTER_ORDER'] = strtoint(headers['FILTER_ORDER'])

    # DYNA dict string data
    header['dyna']['EVENT_NAME'] = headers['EVENT_NAME']
    header['dyna']['EVENT_ID'] = headers['EVENT_ID']
    header['dyna']['EVENT_DATE_YYYYMMDD'] = headers['EVENT_DATE_YYYYMMDD']
    header['dyna']['EVENT_TIME_HHMMSS'] = headers['EVENT_TIME_HHMMSS']
    header['dyna']['MAGNITUDE_W_REFERENCE'] = headers['MAGNITUDE_W_REFERENCE']
    header['dyna']['MAGNITUDE_L_REFERENCE'] = headers['MAGNITUDE_L_REFERENCE']
    header['dyna']['FOCAL_MECHANISM'] = headers['FOCAL_MECHANISM']
    header['dyna']['STATION_NAME'] = headers['STATION_NAME']
    header['dyna']['SITE_CLASSIFICATION_EC8'] = headers[
        'SITE_CLASSIFICATION_EC8']
    header['dyna']['MORPHOLOGIC_CLASSIFICATION'] = headers[
        'MORPHOLOGIC_CLASSIFICATION']
    header['dyna']['DATE_TIME_FIRST_SAMPLE_PRECISION'] = headers[
        'DATE_TIME_FIRST_SAMPLE_PRECISION']
    header['dyna']['STREAM'] = headers['STREAM']
    header['dyna']['UNITS'] = headers['UNITS']
    header['dyna']['INSTRUMENT'] = headers['INSTRUMENT']
    header['dyna']['INSTRUMENT_ANALOG_DIGITAL'] = headers[
        'INSTRUMENT_ANALOG/DIGITAL']
    header['dyna']['BASELINE_CORRECTION'] = headers['BASELINE_CORRECTION']
    header['dyna']['FILTER_TYPE'] = headers['FILTER_TYPE']
    header['dyna']['LATE_NORMAL_TRIGGERED'] = headers['LATE/NORMAL_TRIGGERED']
    header['dyna']['HEADER_FORMAT'] = headers['HEADER_FORMAT']
    header['dyna']['DATABASE_VERSION'] = headers['DATABASE_VERSION']
    header['dyna']['DATA_TYPE'] = headers['DATA_TYPE']
    header['dyna']['PROCESSING'] = headers['PROCESSING']
    header['dyna']['DATA_TIMESTAMP_YYYYMMDD_HHMMSS'] = headers[
        'DATA_TIMESTAMP_YYYYMMDD_HHMMSS']
    header['dyna']['USER1'] = headers['USER1']
    header['dyna']['USER2'] = headers['USER2']
    header['dyna']['USER3'] = headers['USER3']
    header['dyna']['USER4'] = headers['USER4']

    if headonly:
        # skip data
        stream.append(Trace(header=header))
    else:
        # read data
        data = np.loadtxt(fh, dtype='float32')
        stream.append(Trace(data=data, header=header))

    fh.close()
    return stream


def readITACA(filename, headonly=False, **kwargs):  # @UnusedVariable
    """
    Reads a ITACA ASCII file and returns an ObsPy Stream object.

    .. warning::
        This function should NOT be called directly, it registers via the
        ObsPy :func:`~obspy.core.stream.read` function, call this instead.

    :type filename: str
    :param filename: ASCII file to be read.
    :type headonly: bool, optional
    :param headonly: If set to True, read only the head. This is most useful
        for scanning available data in huge (temporary) data sets.
    :rtype: :class:`~obspy.core.stream.Stream`
    :return: A ObsPy Stream object.

    .. rubric:: Example

    >>> from obspy.core import read
    >>> st = read("/path/to/20090330_215717ITDPC_AQG__WEX.DAT")
    >>> st  # doctest: +ELLIPSIS
    <obspy.core.stream.Stream object at 0x...>
    >>> print(st)  # doctest: +ELLIPSIS
    1 Trace(s) in Stream:
    IT.AQG..E | 2009-03-30T21:57:17.967000Z - 2009-03-30T21:57:55.962000Z | 200.0 Hz, 13000 samples
    """
    headers = {}
    data = StringIO()

    # read file
    fh = open(filename, 'rt')
    for i in xrange(43):
        key, value = fh.readline().strip().split(':')
        headers[key.strip()] = value.strip()

    # create ObsPy stream object
    stream = Stream()
    header = Stats()
    header['itaca'] = {}

    header['network'] = filename[-18:-16]
    header['station'] = headers['STATION_CODE']
    header['location'] = ''
    if headers['COMPONENT'] == 'WE':
        header['channel'] = 'HNE'
    # EW should *NEVER* appear, but we handle it anyway
    if headers['COMPONENT'] == 'EW':
        header['channel'] = 'HNE'
    #  -just in case ;)
    if headers['COMPONENT'] == 'NS':
        header['channel'] = 'HNN'
    if headers['COMPONENT'] == 'UP':
        header['channel'] = 'HNZ'
    try:
        tfs = headers['EVENT_DATE_YYYYMMDD'] + \
            '_' + headers['TIME_FIRST_SAMPLE_S']
        # use toUTCDateTime to convert from DYNA format
        header['starttime'] = toUTCDateTime(tfs)
        if re.match('^00', headers['TIME_FIRST_SAMPLE_S']) and re.match('^23', headers['EVENT_TIME_HHMMSS']):
            header['starttime'] = header['starttime'] + 86400
        if re.match('^23', headers['TIME_FIRST_SAMPLE_S']) and re.match('^00', headers['EVENT_TIME_HHMMSS']):
            header['starttime'] = header['starttime'] - 86400
    except:
        header['starttime'] = toUTCDateTime('19700101_000000')
    header['sampling_rate'] = 1 / float(headers['SAMPLING_INTERVAL_S'])
    header['delta'] = float(headers['SAMPLING_INTERVAL_S'])
    header['npts'] = int(headers['NDATA'])
    header['calib'] = 1  # not in file header

    # ITACA dict float data
    header['itaca']['EVENT_LATITUDE_DEGREE'] = strtofloat(
        headers['EVENT_LATITUDE_DEGREE'])
    header['itaca']['EVENT_LONGITUDE_DEGREE'] = strtofloat(
        headers['EVENT_LONGITUDE_DEGREE'])
    header['itaca']['EVENT_DEPTH_KM'] = strtofloat(headers['EVENT_DEPTH_KM'])
    header['itaca']['MAGNITUDE_L'] = strtofloat(headers['MAGNITUDE_L'])
    header['itaca']['MAGNITUDE_S'] = strtofloat(headers['MAGNITUDE_S'])
    header['itaca']['MAGNITUDE_W'] = strtofloat(headers['MAGNITUDE_W'])
    header['itaca']['STATION_LATITUDE_DEGREE'] = strtofloat(
        headers['STATION_LATITUDE_DEGREE'])
    header['itaca']['STATION_LONGITUDE_DEGREE'] = strtofloat(
        headers['STATION_LONGITUDE_DEGREE'])
    header['itaca']['EPICENTRAL_DISTANCE_KM'] = strtofloat(
        headers['EPICENTRAL_DISTANCE_KM'])
    header['itaca']['EARTHQUAKE_BACKAZIMUTH_DEGREE'] = strtofloat(
        headers['EARTHQUAKE_BACKAZIMUTH_DEGREE'])
    header['itaca']['DURATION_S'] = strtofloat(headers['DURATION_S'])
    header['itaca']['INSTRUMENTAL_FREQUENCY_HZ'] = strtofloat(
        headers['INSTRUMENTAL_FREQUENCY_HZ'])
    header['itaca']['INSTRUMENTAL_DAMPING'] = strtofloat(
        headers['INSTRUMENTAL_DAMPING'])
    header['itaca']['FULL_SCALE_G'] = strtofloat(headers['FULL_SCALE_G'])

    # data type is acceleration
    if headers['DATA_TYPE'] == "UNPROCESSED ACCELERATION" \
            or headers['DATA_TYPE'] == "PROCESSED ACCELERATION" \
            or headers['DATA_TYPE'][-8:] == "SPECTRUM":
        header['itaca']['PGA_CM_S_2'] = strtofloat(headers['PGA_CM/S^2'])
        header['itaca']['TIME_PGA_S'] = strtofloat(headers['TIME_PGA_S'])
    # data type is velocity
    if headers['DATA_TYPE'] == "VELOCITY":
        header['itaca']['PGV_CM_S'] = strtofloat(headers['PGV_CM/S'])
        header['itaca']['TIME_PGV_S'] = strtofloat(headers['TIME_PGV_S'])
    # data type is displacement
    if headers['DATA_TYPE'] == "DISPLACEMENT":
        header['itaca']['PGD_CM'] = strtofloat(headers['PGD_CM'])
        header['itaca']['TIME_PGD_S'] = strtofloat(headers['TIME_PGD_S'])

    header['itaca']['LOW_CUT_FREQUENCY_HZ'] = strtofloat(
        headers['LOW_CUT_FREQUENCY_HZ'])
    header['itaca']['HIGH_CUT_FREQUENCY_HZ'] = strtofloat(
        headers['HIGH_CUT_FREQUENCY_HZ'])

    # ITACA dict int data
    header['itaca']['STATION_ELEVATION_M'] = strtoint(
        headers['STATION_ELEVATION_M'])
    header['itaca']['N_BIT_DIGITAL_CONVERTER'] = strtoint(
        headers['N_BIT_DIGITAL_CONVERTER'])
    header['itaca']['FILTER_ORDER'] = strtoint(headers['FILTER_ORDER'])

    # ITACA dict string data
    header['itaca']['EVENT_NAME'] = headers['EVENT_NAME']
    header['itaca']['EVENT_DATE_YYYYMMDD'] = headers['EVENT_DATE_YYYYMMDD']
    header['itaca']['EVENT_TIME_HHMMSS'] = headers['EVENT_TIME_HHMMSS']
    header['itaca']['FOCAL_MECHANISM'] = headers['FOCAL_MECHANISM']
    header['itaca']['STATION_NAME'] = headers['STATION_NAME']
    header['itaca']['SITE_CLASSIFICATION_EC8'] = headers[
        'SITE_CLASSIFICATION_EC8']
    header['itaca']['MORPHOLOGIC_CLASSIFICATION'] = headers[
        'MORPHOLOGIC_CLASSIFICATION']
    header['itaca']['COMPONENT'] = headers['COMPONENT']
    header['itaca']['UNITS'] = headers['UNITS']
    header['itaca']['INSTRUMENT'] = headers['INSTRUMENT']
    header['itaca']['INSTRUMENT_ANALOG_DIGITAL'] = headers[
        'INSTRUMENT_ANALOG/DIGITAL']
    header['itaca']['BASELINE_CORRECTION'] = headers['BASELINE_CORRECTION']
    header['itaca']['FILTER_TYPE'] = headers['FILTER_TYPE']
    header['itaca']['LATE_NORMAL_TRIGGERED'] = headers['LATE/NORMAL_TRIGGERED']
    header['itaca']['DATA_VERSION'] = headers['DATA_VERSION']
    header['itaca']['DATA_TYPE'] = headers['DATA_TYPE']

    if headonly:
        # skip data
        stream.append(Trace(header=header))
    else:
        # read data
        data = np.loadtxt(fh, dtype='float32')
        stream.append(Trace(data=data, header=header))

    fh.close()
    return stream


def writeDYNA(stream, filename, **kwargs):  # @UnusedVariable
    """
    Writes a DYNA ASCII file from given ObsPy Stream object.

    .. warning::
        This function should NOT be called directly, it registers via the
        the :meth:`~obspy.core.stream.Stream.write` method of an
        ObsPy :class:`~obspy.core.stream.Stream` object, call this instead.

    :type stream: :class:`~obspy.core.stream.Stream`
    :param stream: The ObsPy Stream object to write.
    :type filename: str
    :param filename: Name of the ASCII file to write.
    """
    fh = open(filename, 'wb')

    for trace in stream:
        fh.write("EVENT_NAME: %s\n" % trace.stats.dyna.EVENT_NAME)
        fh.write("EVENT_ID: %s\n" % trace.stats.dyna.EVENT_ID)
        fh.write("EVENT_DATE_YYYYMMDD: %s\n" %
                 trace.stats.dyna.EVENT_DATE_YYYYMMDD)
        fh.write("EVENT_TIME_HHMMSS: %s\n" %
                 trace.stats.dyna.EVENT_TIME_HHMMSS)
        fh.write("EVENT_LATITUDE_DEGREE: %s\n" %
                 floattostr(trace.stats.dyna.EVENT_LATITUDE_DEGREE, 4))
        fh.write("EVENT_LONGITUDE_DEGREE: %s\n" %
                 floattostr(trace.stats.dyna.EVENT_LONGITUDE_DEGREE, 4))
        fh.write("EVENT_DEPTH_KM: %s\n" %
                 floattostr(trace.stats.dyna.EVENT_DEPTH_KM, 1))
        fh.write("HYPOCENTER_REFERENCE: %s\n" %
                 trace.stats.dyna.HYPOCENTER_REFERENCE)
        fh.write("MAGNITUDE_W: %s\n" %
                 floattostr(trace.stats.dyna.MAGNITUDE_W, 1))
        fh.write("MAGNITUDE_W_REFERENCE: %s\n" %
                 trace.stats.dyna.MAGNITUDE_W_REFERENCE)
        fh.write("MAGNITUDE_L: %s\n" %
                 floattostr(trace.stats.dyna.MAGNITUDE_L, 1))
        fh.write("MAGNITUDE_L_REFERENCE: %s\n" %
                 trace.stats.dyna.MAGNITUDE_L_REFERENCE)
        fh.write("FOCAL_MECHANISM: %s\n" % trace.stats.dyna.FOCAL_MECHANISM)
        fh.write("NETWORK: %s\n" % trace.stats.network)
        fh.write("STATION_CODE: %s\n" % trace.stats.station)
        fh.write("STATION_NAME: %s\n" % trace.stats.dyna.STATION_NAME)
        fh.write("STATION_LATITUDE_DEGREE: %s\n" %
                 floattostr(trace.stats.dyna.STATION_LATITUDE_DEGREE, 6))
        fh.write("STATION_LONGITUDE_DEGREE: %s\n" %
                 floattostr(trace.stats.dyna.STATION_LONGITUDE_DEGREE, 6))
        fh.write("STATION_ELEVATION_M: %s\n" %
                 floattostr(trace.stats.dyna.STATION_ELEVATION_M, 0))
        fh.write("LOCATION: %s\n" % trace.stats.location)
        fh.write("VS30_M/S: %s\n" % floattostr(trace.stats.dyna.VS30_M_S, 0))
        fh.write("SITE_CLASSIFICATION_EC8: %s\n" %
                 trace.stats.dyna.SITE_CLASSIFICATION_EC8)
        fh.write("MORPHOLOGIC_CLASSIFICATION: %s\n" %
                 trace.stats.dyna.MORPHOLOGIC_CLASSIFICATION)
        fh.write("EPICENTRAL_DISTANCE_KM: %s\n" %
                 floattostr(trace.stats.dyna.EPICENTRAL_DISTANCE_KM, 1))
        fh.write("EARTHQUAKE_BACKAZIMUTH_DEGREE: %s\n" %
                 floattostr(trace.stats.dyna.EARTHQUAKE_BACKAZIMUTH_DEGREE, 1))
        if trace.stats.dyna.DATE_TIME_FIRST_SAMPLE_PRECISION == 'seconds' \
                or trace.stats.dyna.DATE_TIME_FIRST_SAMPLE_PRECISION == 'milliseconds':
            fh.write("DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS: %s\n" %
                     fromUTCDateTime(trace.stats.starttime))
        else:
            fh.write("DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS: \n")
        fh.write("DATE_TIME_FIRST_SAMPLE_PRECISION: %s\n" %
                 trace.stats.dyna.DATE_TIME_FIRST_SAMPLE_PRECISION)
        fh.write("SAMPLING_INTERVAL_S: %s\n" %
                 floattostr(trace.stats.delta, 6))
        fh.write("NDATA: %s\n" % floattostr(trace.stats.npts, 0))
        fh.write("DURATION_S: %s\n" %
                 floattostr(trace.stats.dyna.DURATION_S, 6))
        fh.write("STREAM: %s\n" % trace.stats.dyna.STREAM)
        fh.write("UNITS: %s\n" % trace.stats.dyna.UNITS)
        fh.write("INSTRUMENT: %s\n" % trace.stats.dyna.INSTRUMENT)
        fh.write("INSTRUMENT_ANALOG/DIGITAL: %s\n" %
                 trace.stats.dyna.INSTRUMENT_ANALOG_DIGITAL)
        fh.write("INSTRUMENTAL_FREQUENCY_HZ: %s\n" %
                 floattostr(trace.stats.dyna.INSTRUMENTAL_FREQUENCY_HZ, 3))
        fh.write("INSTRUMENTAL_DAMPING: %s\n" %
                 floattostr(trace.stats.dyna.INSTRUMENTAL_DAMPING, 6))
        fh.write("FULL_SCALE_G: %s\n" %
                 floattostr(trace.stats.dyna.FULL_SCALE_G, 1))
        fh.write("N_BIT_DIGITAL_CONVERTER: %s\n" %
                 floattostr(trace.stats.dyna.N_BIT_DIGITAL_CONVERTER, 0))
        # data type is acceleration
        if trace.stats.dyna.DATA_TYPE == "ACCELERATION" \
                or trace.stats.dyna.DATA_TYPE == "ACCELERATION RESPONSE SPECTRUM":
            fh.write("PGA_CM/S^2: %s\n" %
                     floattostr(trace.stats.dyna.PGA_CM_S_2, 6))
            fh.write("TIME_PGA_S: %s\n" %
                     floattostr(trace.stats.dyna.TIME_PGA_S, 6))
        # data type is velocity
        elif trace.stats.dyna.DATA_TYPE == "VELOCITY" \
                or trace.stats.dyna.DATA_TYPE == "PSEUDO-VELOCITY RESPONSE SPECTRUM":
            fh.write("PGV_CM/S: %s\n" %
                     floattostr(trace.stats.dyna.PGV_CM_S, 6))
            fh.write("TIME_PGV_S: %s\n" %
                     floattostr(trace.stats.dyna.TIME_PGV_S, 6))
        # data type is displacement
        elif trace.stats.dyna.DATA_TYPE == "DISPLACEMENT" \
                or trace.stats.dyna.DATA_TYPE == "DISPLACEMENT RESPONSE SPECTRUM":
            fh.write("PGD_CM: %s\n" % floattostr(trace.stats.dyna.PGD_CM, 6))
            fh.write("TIME_PGD_S: %s\n" %
                     floattostr(trace.stats.dyna.TIME_PGD_S, 6))
        fh.write("BASELINE_CORRECTION: %s\n" %
                 trace.stats.dyna.BASELINE_CORRECTION)
        fh.write("FILTER_TYPE: %s\n" % trace.stats.dyna.FILTER_TYPE)
        fh.write("FILTER_ORDER: %s\n" %
                 floattostr(trace.stats.dyna.FILTER_ORDER, 0))
        fh.write("LOW_CUT_FREQUENCY_HZ: %s\n" %
                 floattostr(trace.stats.dyna.LOW_CUT_FREQUENCY_HZ, 3))
        fh.write("HIGH_CUT_FREQUENCY_HZ: %s\n" %
                 floattostr(trace.stats.dyna.HIGH_CUT_FREQUENCY_HZ, 3))
        fh.write("LATE/NORMAL_TRIGGERED: %s\n" %
                 trace.stats.dyna.LATE_NORMAL_TRIGGERED)
        fh.write("DATABASE_VERSION: %s\n" % trace.stats.dyna.DATABASE_VERSION)
        fh.write("HEADER_FORMAT: %s\n" % trace.stats.dyna.HEADER_FORMAT)
        fh.write("DATA_TYPE: %s\n" % trace.stats.dyna.DATA_TYPE)
        fh.write("PROCESSING: %s\n" % trace.stats.dyna.PROCESSING)
        fh.write("DATA_TIMESTAMP_YYYYMMDD_HHMMSS: %s\n" %
                 trace.stats.dyna.DATA_TIMESTAMP_YYYYMMDD_HHMMSS)
        fh.write("USER1: %s\n" % trace.stats.dyna.USER1)
        fh.write("USER2: %s\n" % trace.stats.dyna.USER2)
        fh.write("USER3: %s\n" % trace.stats.dyna.USER3)
        fh.write("USER4: %s\n" % trace.stats.dyna.USER4)

        if trace.stats.dyna.DATA_TYPE == "ACCELERATION" \
                or trace.stats.dyna.DATA_TYPE == "VELOCITY":
            for d in trace.data:
                fh.write("%-.6f\n" % d)
        elif trace.stats.dyna.DATA_TYPE == "DISPLACEMENT":
            for d in trace.data:
                fh.write("%e\n" % d)
        elif trace.stats.dyna.DATA_TYPE[-8:] == "SPECTRUM":
            for j in xrange(len(trace.data)):
                for i in xrange(2):
                    if i == 0:
                        fh.write("%12.6f" % trace.data[j][i])
                    elif i == 1:
                        fh.write("%13.6f\n" % trace.data[j][i])
    fh.close()


def toUTCDateTime(value):
    """
    Converts time string used within DYNA/ITACA File into a UTCDateTime.

    :type value: str
    :param value: A Date time string.
    :return: Converted :class:`~obspy.core.UTCDateTime` object.

    .. rubric:: Example

    >>> toUTCDateTime('20090330_215650.967')
    UTCDateTime(2009, 3, 30, 21, 56, 50, 967000)
    """
    try:
        date, time = value.split('_')
    except ValueError:
        date = value

    year = int(date[0:4])
    month = int(date[4:6])
    day = int(date[6:8])

    hour = int(time[0:2])
    mins = int(time[2:4])
    secs = float(time[4:])

    return UTCDateTime(year, month, day, hour, mins) + secs


def fromUTCDateTime(dt):
    """
    Converts UTCDateTime object into a time string used within Seismic Handler.

    :type dt: :class:`~obspy.core.UTCDateTime`
    :param dt: A UTCDateTime object.
    :return: Converted date time string as defined in DYNA file header.

    .. rubric:: Example

    >>> from obspy.core import UTCDateTime
    >>> dt = UTCDateTime(2009, 3, 30, 21, 56, 50, 967000)
    >>> fromUTCDateTime(dt)
    '20090330_215650.967'
    """
    pattern = "%4d%02d%02d_%02d%02d%02d.%03d"

    return pattern % (dt.year, dt.month, dt.day, dt.hour,
                      dt.minute, dt.second, dt.microsecond / 1000)


def strtofloat(sf):
    try:
        x = float(sf)
    except:
        return None
    return x


def strtoint(sf):
    try:
        x = int(sf)
    except:
        return None
    return x


def floattostr(fs, n):
    y = format("%-.0f" % n)
    try:
        x = eval('format("%-.' + y + 'f" % fs)')
    except:
        return ''
    return x

if __name__ == '__main__':
    import doctest
    doctest.testmod(exclude_empty=True)

'''
# SAMPLE WITH DATA TYPES FROM DYNA FORMAT (MEMO) :)

####common header data
NETWORK: IT
STATION_CODE: AQG
LOCATION:
DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS: 20090330_215650.967
STREAM: N
SAMPLING_INTERVAL_S: 0.005000
NDATA: 13000

###DYNA dict string data

###DYNA dict float data

###DYNA dict int data

'''

'''
DYNA fields not in itaca header
itaca to dyna conversion example (memo):
'''
