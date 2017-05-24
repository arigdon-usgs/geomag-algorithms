"""Algorithm that creates an averaged Dst.

"""

from Algorithm import Algorithm
from AlgorithmException import AlgorithmException
import numpy
import obspy.core


# Possible correction factors.
# Defaults to 1.0 if station not found in list.
CORR = {
    'HON': 1.0,
    'SJG': 1.0,
    'HER': 1.0,
    'KAK': 1.0,
    'GUA': 1.0
}


class AverageAlgorithm(Algorithm):
    """Algorithm that creates an averaged Dst.

    Parameters
    ----------

    """

    def __init__(self):
        Algorithm.__init__(self)
        self._npts = -1
        self._stt = -1
        self._stats = None


    def check_stream(self, timeseries):
        """checks a stream to make certain the required data
            exists.

        Parameters
        ----------
        timeseries: obspy.core.Stream
            stream to be checked.
        """

        # A stream produced by EdgeFactory should always pass these checks.

        # must have a timeseries for all observatories
        if len(timeseries) != len(self.observatories):
            raise AlgorithmException(
                'Expected data for %d stations, received %d'
                    % (len(self.observatories), len(timeseries)) )

        # timeseries starttime and number of samples must match
        for ts in timeseries:
            # grab 1st set of stats to use in output.
            # Its values will be good if checks pass.
            if self._stats is None:
                self._stats = ts.stats

            if self._npts == -1:
                self._npts = ts.stats.npts
            if ts.stats.npts != self._npts:
                raise AlgorithmException(
                    'Received timeseries have different lengths')

            if self._stt == -1:
                self._stt = ts.stats.starttime
            if ts.stats.starttime != self._stt:
                raise AlgorithmException(
                    'Received timeseries have different starttimes')


    def process(self, timeseries):
        """averages a channel across multiple stations

        Parameters
        ----------

        Returns
        -------
        out_stream:
            new stream object containing the averaged values.
        """

        self.check_stream(timeseries)

        # initialize sum and count series
        dst_sum = numpy.full(self._npts,0.)
        dst_count = numpy.full(self._npts,0.)

        # loop over stations
        for obsy in self.observatories:

            # lookup latitude correction factor, default = 1.0
            latcorr = 1.0
            if CORR.has_key(obsy):
                latcorr = CORR[obsy]

#            print '============================'
#            print 'Obsy: %s latcorr: %f' % (obsy, latcorr)
#            print '============================'

            # loop thru the station's timeseries, add live
            # data to dst_sum and dst_count arrays
            ts = timeseries.select(station=obsy)[0]
            for i in range(0, self._npts):

                if not numpy.isnan(ts.data[i]):
                    dst_sum[i] += latcorr * ts.data[i]
                    dst_count[i] += 1
#                    print 'samp: %d val*latcorr: %f' % (i, ts.data[i])

        # after looping over stations, compute average
        # FOR NOW: set average to NaN if any individual values are
        # missing, i.e., have live data for all stations
        for i in range(0, self._npts):
            if dst_count[i] == len(self.observatories):
                dst_sum[i] /= dst_count[i]
            else:
                dst_sum[i] = numpy.nan

        # return averaged values as a stream with channel = MSD
        # station will be set to USGS in get_trace
        return obspy.core.Stream((
                get_trace('MSD', self._stats, dst_sum), ))


    @classmethod
    def add_arguments(cls, parser):
        """Add command line arguments to argparse parser.

        Parameters
        ----------
        parser: ArgumentParser
            command line argument parser
        """
        parser.add_argument('--observatory-scale',
               default=(None,),
               help='Scale factor for observatories specified with ' +
                    '--observatory argument',
               nargs='*',
               type=float)


    def configure(self, arguments):
        """Configure algorithm using comand line arguments.

        Parameters
        ----------
        arguments: Namespace
            parsed command line arguments
        """

        self.observatories = arguments.observatory
        self.scales = arguments.observatory_scale
        print self.scales
        if self.scales[0] is not None:
            if len(self.observatories) != len(self.scales):
                raise AlgorithmException(
                  'Mismatch between observatories and scale factors')
            else:
                for i,obs in enumerate(self.observatories):
                   CORR[obs] = self.scales[i]



def get_trace(channel, stats, data):
    """Utility to create a new trace object.

    Parameters
    ----------
    channel : str
        channel name.
    stats : obspy.core.Stats
        channel metadata to clone.
    data : numpy.array
        channel data.

    Returns
    -------
    obspy.core.Trace
        trace containing data and metadata.
    """
    stats = obspy.core.Stats(stats)

    stats.channel = channel
#    stats.channel = 'MSD'
    stats.station = 'USGS'
    stats.network = 'NT'
    stats.location = 'R0'
    return obspy.core.Trace(data, stats)
