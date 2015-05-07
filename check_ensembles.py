# -*- coding: utf-8 -*-

'''
Author: Harald von Waldow <hvwaldow@env.ethz.ch>
Date: April 2015

Utility classes to check climate simulation output for completeness
and soundness. At the moment this is applied to the ENSEMBLES dataset,
but it should be possible to also use it for other projects, such as
CORDEX or CMIP.

Of course, only a very specific subset of file characteristics can be
evaluated. A negative check result by no means indicates that the file
or dataset is error-free!
'''

import os
import sys
import netCDF4 as nc
import numpy as np
import glob
import unittest
import datetime
import pprint as pp
import re
import cPickle
import multiprocessing as mp
from collections import Counter
from scipy.stats import describe

# Pattern that extracts start- and endyear from filename
# into group(1) and group(2) of re.MatchObject respectively.
FN_YEAR_PATTERN = r".*_(\d{4})-(\d{4})_.*.nc"

# Path-prefix to reach storage-host
PREFIX = "/"


class CFConventions(object):
    '''Some specs from CF-1.6'''
    def __init__(self):
        self.valid_periods = ["millisecs", "msecs", "seconds", "secs",
                              "minutes", "mins", "hours", "hrs", "days",
                              "weeks"]
        self.valid_calendars = ['gregorian', 'standard',
                                'proleptic_gregorian ', 'noleap', '365_day',
                                'noleap', '365_day', '360_day', 'julian']


class GlobalCheck(object):
    '''
    Produce file counts for all paths and patterns
    Holds path/pattern combinations
    '''
    def __init__(self, paths, variables, output=True):
        self.combis = []
        for path in paths:
            for v in variables:
                pat = self.mk_pattern(v)
                self.combis.append((path, pat))
                if output:
                    tc = TimeCheck(path, pat)
                    tc.get_files()
                    nofiles = len(tc.files)
                    print('Dataset: {} \t {} files:'
                          .format(path+'/'+pat, nofiles))
                    print('\n'.join(tc.files))

    def mk_pattern(self, variable):
        return(".*_{}\.nc".format(variable))


class TimeCheck(object):
    '''
    Operates on a "datasets", that is a set of netCDF files, whichg
    only differ in the time span covered. Initialize with a path and
    a regular expression that define the files of the dataset.
    '''
    def __init__(self, path, pattern):
        self.path = path
        self.pattern = pattern
        self.cfco = CFConventions()
        self.get_files()
        self.sort_files()
        self._mk_ds()
        self.MAXYEARS = 50.

    def init(self):
        '''
        Manual call to init required because already quite
        some testing happens here.
        '''
        print("checking dataset: {}/{}".format(self.path, self.pattern))
        print("Has time:units? {}"
              .format(self.check_has_attr('time', 'units')))
        print("Has same time:units ? {} - {}"
              .format(self.check_has_same_attr('time', 'units'), self.units))
        print("time:units = \"days\"? {}"
              .format(self.check_timeunit_is_days()))
        print("Has time:calendar? {}"
              .format(self.check_has_attr('time', 'calendar')))
        print("Has same calendar? {} - {}"
              .format(self.check_has_same_attr('time', 'calendar'), self.calendar))
        print("Has valid calendar unit? {}"
              .format(self.check_calendar_unit_is_valid()))

    def get_files(self):
        self.files = [os.path.join(self.path, f) for f in os.listdir(self.path)
                      if re.search(self.pattern, f)]
        self.files.sort()
        if len(self.files) == 0:
            raise RuntimeError(
                "No files found for dataset with path={} and pattern={}."
                .format(self.path, self.pattern))

    def sort_files(self):
        '''Sorts files according to first date.'''
        starttimes = [nc.Dataset(f).variables['time'][0] for f in self.files]
        filetuples = zip(self.files, starttimes)
        filetuples.sort(key=lambda f: f[1])
        self.files = [f[0] for f in filetuples]

    def _mk_ds(self):
        '''Returns netCDF4 - Objects for all files'''
        self.ds = [nc.Dataset(f) for f in self.files]

    def check_has_attr(self, variable, attribute):
        res = [attribute in d.variables[variable].ncattrs() for d in self.ds]
        if not all(res):
            nocal = [f[0] for f in zip(self.files, res) if not f[1]]
            print('''ERROR:
            Files without attribute \"{}\" for variable \"{}\":
            {}'''.format(attribute, variable, nocal))
            return("ERROR")
        else:
            return("OK")

    def check_has_same_attr(self, variable, attribute):
        atts = [d.variables[variable].getncattr(attribute) for d in self.ds]
        if len(set(atts)) > 1:
            pp.pprint("ERROR: Not all files have the same {} for {}:"
                      .format(attribute, variable))
            pp.pprint(dict(zip(self.files, atts)))
            return("ERROR")
        else:
            exec("self.{} = \"{}\"".format(attribute, atts[0]))
            return("OK")

    def check_calendar_unit_is_valid(self):
        '''Check for valid (according to CF-1.6) calendar'''
        cal = self.ds[0].variables['time'].getncattr('calendar')
        if not (cal in self.cfco.valid_calendars):
            print("ERROR: Calendar {} not valid in CF-1.6".format(cal))
            return("ERROR")
        else:
            return("OK")

    def check_timeunit_is_days(self):
        units = [d.variables['time'].getncattr('units') for d in self.ds]
        if re.match("days", units[0]):
            return("OK")
        else:
            print("ERROR: time units not \"days\".Other units currently" +
                  " not implemented.")
        return("ERROR")

    def time_get_info(self):
        for d in self.ds:
            tim = d.variables['time']
            timesteps = tim[:]
            span = (timesteps[0], timesteps[-1])
            filename = os.path.basename(d.filepath())
            yield((filename, timesteps, span))
            
    def find_interfile_problems(self):
        err = False
        for idx in range(0, len(self.files) - 1):
            a_end = self.ds[idx].variables['time'][-1]
            b_start = self.ds[idx + 1].variables['time'][0]
            if not (b_start == a_end + 1):
                fna = os.path.basename(self.ds[idx].filepath())
                fnb = os.path.basename(self.ds[idx + 1].filepath())
                print("ERROR: non-continuous between {} and {}"
                      .format(fna, fnb))
                err = True
        if not err:
            print("{} OK".format(os.path.join(self.path, self.pattern)))

    def check_timesteps(self):
        for f in self.time_get_info():
            print(f[0]),
            if self.check_timestep_span(f):
                print("ERROR - Timesteps fundamentally flawed")
                continue
            try:
                shouldsteps = np.arange(f[2][0], f[2][1] + 1)
            except Exception as e:
                print("ERROR: {}".format(e.message or "unknown"))
                continue
            equal = f[1] == shouldsteps
            if (type(equal) == bool) or (not all(equal)):
                # comparison resulted in single bool or not all True
                diffmsg = self.compare_ts(f, shouldsteps)
                if diffmsg:
                    print("ERROR - " + diffmsg)
                else:
                    print("ERROR - Timesteps somehow wrong")
            else:
                print("OK")

    def check_timestep_span(self, f):
        duration = f[2][1] - f[2][0]
        if duration < 360 or duration > self.MAXYEARS * 366:
            return(True)
        else:
            return(False)

    def compare_ts(self, f, shouldsteps):
        ts = f[1]
        sts = set(ts)
        sshould = set(shouldsteps)
        # dubletten
        if len(sts) != len(ts):
            dubs = [k for k, v in Counter(ts).items() if v > 1]
            dubsdates = nc.num2date(dubs, self.units, self.calendar)
            dubsdates = [x.strftime("%Y-%m-%d") for x in dubsdates]
            return("non-unique timesteps: {} = {}".format(dubs, dubsdates))
        if sts.issubset(sshould):
            miss = list(sshould - sts)
            missdates = nc.num2date(miss, self.units, self.calendar)
            missdates = [x.strftime("%Y-%m-%d") for x in missdates]
            return("missing timesteps: {}={}".format(miss, missdates))
        elif sshould.issubset(sts):
            extra = list(sts - sshould)
            extradates = nc.num2date(extra, self.units, self.calendar)
            extradates = [x.strftime("%Y-%m-%d") for x in extradates]
            return("extra timesteps: {}={}".format(extra, extradates))
        else:
            nomatch = list(sts ^ sshould)
            if len(nomatch) > 0:
                nomatchdates = nc.num2date(nomatch, self.units, self.calendar)
                nomatchdates = [x.strftime("%Y-%m-%d") for x in nomatchdates]
                return("non-matching timesteps: {}={}"
                       .format(nomatch, nomatchdates))
            else:
                return(False)


class FixTimeAxis(object):
    ''' Replaces time-axis of a netcdf-file base on time:calendar,
    time:units and years taken from filename'''
    def __init__(self, filename):
        print("reading file: {}".format(filename))
        self.filename = filename
        res = re.search(FN_YEAR_PATTERN,
                        os.path.basename(filename))
        self.starty = int(res.group(1))
        self.endy = int(res.group(2))
        print("startyear: {}   endyear: {}".format(self.starty, self.endy))

    def new_timeaxis(self):
        print("replacing time-axis")
        ds = nc.Dataset(self.filename, "a")
        units = ds.variables['time'].getncattr('units')
        calendar = ds.variables['time'].getncattr('calendar')
        ts = ds.variables['time'][:]
        print("calendar: {}".format(calendar))
        print("units: {}".format(units))
        if not re.match("days since", units):
            sys.exit("units of days necessary; aborting")
        print("No timesteps: {}".format(len(ts)))
        print("ts[0]: {}   ts[end]: {}".format(ts[0], ts[-1]))
        start = datetime.datetime(self.starty, 1, 1, 12)
        end = datetime.datetime(self.endy, 12, 31, 12)
        diff = end - start
        alldays = [start + datetime.timedelta(days=x)
                   for x in range(0, diff.days + 1)]
        if calendar == '360_day':
            # We assume the file contains full years !!!
            no_years = len(set([x.year for x in alldays]))
            days_360 = 360 * no_years
            ts_start = nc.date2num(start, units, calendar)
            ts_new = np.arange(ts_start, ts_start + days_360)
            print("replacing with ts[0]: {}   ts[end]: {}"
                  .format(ts_new[0], ts_new[-1]))
            ds.variables['time'][:] = ts_new
            ds.sync()
            ds.close()
            print("finished")
        elif calendar == 'standard':
            ts_start = nc.date2num(start, units, calendar)
            ts_new = np.arange(ts_start, ts_start + len(alldays))
            print("replacing with ts[0]: {}   ts[end]: {}"
                  .format(ts_new[0], ts_new[-1]))
            ds.variables['time'][:] = ts_new
            ds.sync()
            ds.close()
            print("finished")
        else:
            print("Calendar {} not implemented".format(calendar))


class CheckValues(object):
    def __init__(self, path, pattern, variable):
        self.variable = variable
        self.path = path
        self.pattern = pattern
        self.get_files()
        self.sort_files()
        self._mk_ds()

    def get_files(self):
        self.files = [os.path.join(self.path, f) for f in os.listdir(self.path)
                      if re.search(self.pattern, f)]
        self.files.sort()
        if len(self.files) == 0:
            raise RuntimeError(
                "No files found for dataset with path={} and pattern={}."
                .format(self.path, self.pattern))

    def sort_files(self):
        '''Sorts files according to first date.'''
        starttimes = [nc.Dataset(f).variables['time'][0] for f in self.files]
        filetuples = zip(self.files, starttimes)
        filetuples.sort(key=lambda f: f[1])
        self.files = [f[0] for f in filetuples]

    def _mk_ds(self):
        '''Returns netCDF4 - Objects for all files'''
        self.ds = [nc.Dataset(f) for f in self.files]

    def collect_values(self):
        mami = []
        fields = (x.variables[self.variable][:] for x in self.ds)
        mami = np.array([(np.min(x), np.max(x)) for x in fields])
        print(mami)
            

   
# Functions to run checks and fixes in parallel
def run_tests(combi):
    tc = TimeCheck(combi[0], combi[1])
    print("Initial checks")
    tc.init()
    print("Checking timesteps")
    tc.check_timesteps()
    print("Checking continuity:")
    tc.find_interfile_problems()


def run_fix(f):
    tf = FixTimeAxis(f)
    tf.new_timeaxis()


if __name__ == "__main__":

# ############## MODIFY THIS ##################################################

    # paths to model/experiment specific simulation data
    paths = ["data/ENSEMBLES-RCM/A1B/CNRM_ARPEGE_new/DM",
             "data/ENSEMBLES-RCM/A1B/DMI_ECHAM5/DM",
             "data/ENSEMBLES-RCM/A1B/ETHZ/DM",
             "data/ENSEMBLES-RCM/A1B/HadRM3Q0/DM",
             "data/ENSEMBLES-RCM/A1B/ICTP/DM",
             "data/ENSEMBLES-RCM/A1B/KNMI/DM",
             "data/ENSEMBLES-RCM/A1B/MPI/DM",
             "data/ENSEMBLES-RCM/A1B/SMHI_BCM/DM",
             "data/ENSEMBLES-RCM/A1B/SMHI_ECHAM5/DM",
             "data/ENSEMBLES-RCM/A1B/SMHI_HadCM3Q3_new/DM"]
    variables = ["tasmax", "tasmin", "pr", "rsds", "hurs", "wss"]


#  path to directory that contains files in the need of fixing
#  files that need fixing
    # fixpath = "data/ENSEMBLES-RCM/A1B/SMHI_HadCM3Q3_new/DM"

    # fixfiles = [
    #     "SMHIRCA_A1B_HadCM3Q3_DM_25km_1951-1960_tasmin.nc",
    #     "SMHIRCA_A1B_HadCM3Q3_DM_25km_1961-1970_tasmin.nc",
    #     "SMHIRCA_A1B_HadCM3Q3_DM_25km_1971-1980_pr.nc",
    #     "SMHIRCA_A1B_HadCM3Q3_DM_25km_1981-1990_pr.nc",
    #     "SMHIRCA_A1B_HadCM3Q3_DM_25km_1991-2000_pr.nc",
    #     "SMHIRCA_A1B_HadCM3Q3_DM_25km_1951-1960_tasmax.nc",
    #     "SMHIRCA_A1B_HadCM3Q3_DM_25km_1961-1970_tasmax.nc",
    #     "SMHIRCA_A1B_HadCM3Q3_DM_25km_1961-1970_rsds.nc",
    #     "SMHIRCA_A1B_HadCM3Q3_DM_25km_1971-1980_rsds.nc",
    #     "SMHIRCA_A1B_HadCM3Q3_DM_25km_1951-1960_hurs.nc",
    #     "SMHIRCA_A1B_HadCM3Q3_DM_25km_1961-1970_hurs.nc",
    #     "SMHIRCA_A1B_HadCM3Q3_DM_25km_1951-1960_wss.nc",
    #     "SMHIRCA_A1B_HadCM3Q3_DM_25km_1961-1970_wss.nc"]

    fixpath = "data/ENSEMBLES-RCM/A1B/SMHI_BCM/DM"
    fixfiles = ["SMHIRCA_A1B_BCM_DM_25km_2061-2070_hurs.nc",
                "SMHIRCA_A1B_BCM_DM_25km_2071-2080_hurs.nc",
                "SMHIRCA_A1B_BCM_DM_25km_2081-2090_hurs.nc"]


# ##############################################################################
    paths = [os.path.join(PREFIX, x) for x in paths]
    fixpaths = [os.path.join(PREFIX, fixpath, f) for f in fixfiles]
    
    # p = mp.Pool(processes=8)
    # p.map(run_fix, fixpaths)
 
    gc = GlobalCheck(paths, variables, output=True)
    p = mp.Pool()
    p.map(run_tests, gc.combis)

    # valuepaths = ["/data/ENSEMBLES-RCM/A1B/SMHI_BCM/DM"]
    # valuevars = ['tasmax']
    # gc = GlobalCheck(valuepaths, valuevars, output=True)
    # cv = CheckValues(gc.combis[0][0], gc.combis[0][1], 'tasmax')
    # cv.collect_values()

