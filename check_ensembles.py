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
        #print("Checking path: {}\n         pattern: {}".format(path, pattern))
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
        print("Has same time:units ? {}"
              .format(self.check_has_same_attr('time', 'units')))
        print("time:units = \"days\"? {}"
              .format(self.check_timeunit_is_days()))
        print("Has time:calendar? {}"
              .format(self.check_has_attr('time', 'calendar')))
        print("Has same calendar? {}"
              .format(self.check_has_same_attr('time', 'calendar')))
        print("Has valid calendar unit? {}"
              .format(self.check_calendar_unit_is_valid()))
        # self.time_get_info()
        # self.problemdict = {}
        # self.find_interfile_problems()

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

    # def check_has_variable(self, variable):
    #     '''Checks whether all files have variable "variable"'''
    #     res = [variable in d.variables.keys() for d in self.ds]
    #     if not all(res):
    #         notime = [f[0] for f in zip(self.files, res) if not f[1]]
    #         raise RuntimeError(
    #             "Files without variable \"{}\" found: {}"
    #             .format(variable, notime))
    #     else:
    #         print("All files have variable \"{}\".".format(variable))
    #         return(True)

    # def ncdftime2datetime(self, t):
    #     return(datetime.datetime(*t.timetuple()[0:5]))

    def time_get_info(self):
        for d in self.ds:
            tim = d.variables['time']
            timesteps = tim[:]
            span = (timesteps[0], timesteps[-1])
            filename = os.path.basename(d.filepath())
            yield((filename, timesteps, span))
            
        # tim = [d.variables['time'] for d in self.ds]
        # timesteps = [t[:] for t in tim]
        # spans = [(t[0], t[-1]) for t in timesteps]
        # shouldsteps = [np.arange(s[0], s[1]+1) for s in spans]
        # self.time_info = zip(tim, timesteps, spans, shouldsteps) 

    # def find_interfile_problems(self):
    #     continuity = []
    #     for idx in range(0, len(self.files) - 1):
    #         continuity.append(
    #             self.time_info[idx+1][2][0] == self.time_info[idx][2][1] + 1)
    #     if all(continuity):
    #         print("No gaps in-between files")
    #     else:
    #         print("Gaps between files detected")
    #         self.gaps = continuity

    # def _time_global_check(self):
    #     globerrors = [list(set(t[1]) ^ set(t[3])) for t in self.time_info]
    #     anygloberrors = any([len(g) > 0 for g in globerrors])
    #     return(globerrors if anygloberrors else False)

    def check_timesteps(self):
        print("Checking timesteps")
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
            return("non-unique timesteps: {}".format(dubs))
        if sts.issubset(sshould):
            miss = sshould - sts
            return("missing timesteps: {}".format(list(miss)))
        elif sshould.issubset(sts):
            extra = sts - sshould
            return("extra timesteps: {}".format(list(extra)))
        else:
            nomatch = sts ^ sshould
            if len(nomatch) > 0:
                return("non-matching timesteps: {}".format(list(nomatch)))
            else:
                return(False)




                    
            
        # if self.no_dates_should > self.tim.shape[0]:
        #     print("dates missing!")
        #     miss = self.time_find_missing()
        #     print(miss)
        #     return(-1)
        # elif self.tim.shape[0] > self.no_dates_should:
        #     print("too much dates!")
        #     return(1)
        # else:
        #     print("number of dates ok!")
        #     return(0)

    # def time_find_missing(self):
    #     sshould = set(np.arange(self.tim[0], self.tim[-1]+1))
    #     sis = set(self.timesteps)
    #     return(list(sshould - sis))

    # def time_check(self):
    #     sshould = set(np.arange(self.tim[0], self.tim[-1]+1))
    #     sis = set(self.timesteps)
    #     return(list(sshould ^ sis))

    # def check_continuity(self, date1, date2, calendar):
    #     if (calendar == '360_day' and
    #         ((date1.year == date2.year and
    #           date2.month == date1.month + 1 and
    #           date1.day == 30 & date2.day == 1) or
    #          (date1.year == date2.year - 1 and
    #           date1.month == 12 and date2.month == 1 and
    #           date1.day == 30 and date2.day == 1))):
    #         deltat = datetime.timedelta(days=2)
    #     else:
    #         deltat = datetime.timedelta(days=1)
    #     return(date1 + deltat == date2)
    
    # def get_timestep_count(self):
    #     tscount = []
    #     for d in self.ds:
    #         tc.get_files()
    #         tc._mk_ds()
    #         tscount.append(d.variables['time'][:].shape[0])

    # def mk_boxplot(self):
    #     for d in self.ds:
    #         a = d.variables['time']
    #         print(a.shape)

        


    
    

if __name__ == "__main__":
    prefix = ""
    paths = ["{}/data/ENSEMBLES-RCM/A1B/CNRM_ARPEGE_new/DM".format(prefix),
             "{}/data/ENSEMBLES-RCM/A1B/DMI_ECHAM5/DM".format(prefix),
             "{}/data/ENSEMBLES-RCM/A1B/ETHZ/DM".format(prefix),
             "{}/data/ENSEMBLES-RCM/A1B/HadRM3Q0/DM".format(prefix),
             "{}/data/ENSEMBLES-RCM/A1B/ICTP/DM".format(prefix),
             "{}/data/ENSEMBLES-RCM/A1B/KNMI/DM".format(prefix),
             "{}/data/ENSEMBLES-RCM/A1B/MPI/DM".format(prefix),
             "{}/data/ENSEMBLES-RCM/A1B/SMHI_BCM/DM".format(prefix),
             "{}/data/ENSEMBLES-RCM/A1B/SMHI_ECHAM5/DM".format(prefix),
             "{}/data/ENSEMBLES-RCM/A1B/SMHI_HadCM3Q3_new/DM".format(prefix)]

    variables = ["tasmax", "tasmin", "pr", "rsds", "hurs", "wss"]

    # sizdict = {}
    def run_tests(combi):
        tc = TimeCheck(combi[0], combi[1])
        tc.check_timesteps()
        #tc.init()
        
    gc = GlobalCheck(paths, variables, output=False)
    p = mp.Pool()
    p.map(run_tests, gc.combis)


### Ad Hoc

def rewrite_timesteps(filename, start, end):
    ds = nc.Dataset(filename, "a")
    units = ds.variables['time'].getncattr('units')
    calendar = ds.variables['time'].getncattr('calendar')
    ts = ds.variables['time'][:]
    print("Calendar: {}".format(calendar))
    print("Units: {}".format(units))
    print("No timesteps: {}".format(len(ts)))
    print("ts[0]: {}   ts[end]: {}".format(ts[0], ts[-1]))
    diff = end - start
    alldays = [start + datetime.timedelta(days=x)
               for x in range(0, diff.days + 1)]
    if calendar == '360_day':
        ## We assume the file contains full years !!!
        no_years =  len(set([x.year for x in alldays]))
        days_360 = 360 * no_years
        ts_start = nc.date2num(start, units, calendar)
        ts_new = np.arange(ts_start, ts_start + days_360)
        ds.variables['time'][:] = ts_new
        ds.sync()
        ds.close()
        #jstart = nc.date2num(start, units, calendar)
        
    else:
        print("Calendar {} not implemented".format(calendar))

    

# filename = ("/net/atmos/data/ENSEMBLES-RCM/A1B/SMHI_HadCM3Q3_new/DM/" +
#             "SMHIRCA_A1B_HadCM3Q3_DM_25km_1951-1960_tasmax.nc")
# filename = ("/net/atmos/data/ENSEMBLES-RCM/A1B/SMHI_HadCM3Q3_new/DM/" +
#             "SMHIRCA_A1B_HadCM3Q3_DM_25km_1971-1980_tasmax.nc")
# filename = ("/net/atmos/data/ENSEMBLES-RCM/A1B/SMHI_HadCM3Q3_new/DM/" +
#             "SMHIRCA_A1B_HadCM3Q3_DM_25km_1981-1990_tasmax.nc")
# filename = ("/net/atmos/data/ENSEMBLES-RCM/A1B/SMHI_BCM/DM/" +
#             "SMHIRCA_A1B_BCM_DM_25km_2061-2070_hurs.nc")
filename="SMHIRCA_A1B_HadCM3Q3_DM_25km_1971-1980_pr.nc"


start = datetime.datetime(1971,1,1,12)
end = datetime.datetime(1980,12,31,12)
res = rewrite_timesteps(filename,start,end)

LAUFT !!! Jetzt echt reparieren!



# rewrite_timesteps(

# for combi in gc.combis:
#         tc = TimeCheck(combi[0], combi[1])
#         #tc.init()
#         # for i in tc.time_get_info(): 
#         #     print(i)
        
#         del tc

        # tc.get_files()
        # tc._mk_ds()
        # for d in tc.ds:
        #     sizdict[os.path.basename(d.filepath())] = [d.variables['time'].shape,
        #                                               d.variables['time'].getncattr('calendar')] 
        # tc.mk_boxplot()
    # print(sizdict)



    # print(tc.files)
    # print(tc.check_has_time())
    # print(tc.check_calendar())

    # basepath = '
    # models = ['CNRM_ARPEGE_new', 'DMI_ECHAM5', 'ETHZ', 'HadRM3Q0',
    #           'ICTP', 'KNMI', 'MPI', 'SMHI_BCM', 'SMHI_ECHAM5',
    #           'SMHI_HadCM3Q3_new']
    # time_frequency = 'DM'
    # variables = ['tasmax', 'tasmin', 'pr', 'rsds', 'hurs', 'wss']

    # model = models[7]
    # problemdict = {}
    # resfilename = "tim_check_result_{}.cpy".format(model)
    # for variable in variables:
    #     problemdict[(model, variable)] = {}
    #     try:
    #         TC = TimeCheck(basepath, model, time_frequency, variable)
    #     except RuntimeError as rte:
    #         print(str(rte))
    #         sys.exit(1)
    #     old_enddate = False
    #     for fn in TC.filenames:
    #         TC.set_file(fn)
    #         print(fn)
    #         print("From: {} To: {}".format(TC.startdate, TC.enddate))
    #         if old_enddate:
    #             if TC.check_continuity(old_enddate, TC.startdate, TC.tim.calendar):
    #                 print("file continuity: OK!")
    #             else:
    #                 print("ERROR: files not continuous!")
    #                 if 'continuity' in problemdict[(model, variable)]:
    #                     problemdict[(model, variable)]['continuity'].append(
    #                         (os.path.basename(fn), old_enddate, TC.startdate))
    #                 else:
    #                     problemdict[(model, variable)]['continuity'] = [
    #                         (os.path.basename(fn), old_enddate, TC.startdate)]
    #         else:
    #             problemdict[(model, variable)]['fromto'] = [str(TC.startdate)]
    #         old_enddate = TC.enddate
    #         globdiff = TC.time_check()
    #         if len(globdiff) == 0:
    #             print("Timesteps as expected!")
    #         else:
    #             print("Problem with timesteps here:")
    #             print(globdiff)
    #             miss = TC. time_find_missing()
    #             misslist = [TC.ncdftime2datetime(
    #                 nc.num2date(x, TC.tim.units, TC.tim.calendar))
    #                 for x in miss]
    #             print("Missing: {}".format(misslist))
    #             problemdict[(model, variable)]['missing'] = (
    #                 os.path.basename(fn), globdiff,
    #                 [str(x) for x in misslist])
    #     problemdict[(model, variable)]['fromto'].append(str(TC.enddate))
    # cPickle.dump(problemdict, open(resfilename, 'wb'), protocol=-1)

    #
    # get wrong calendar entries
    # parallel "ncdump -h  {} | grep 'calendar = \"days\"' | sed 's/\(^.*$\)/\1 {}/p' " ::: *.nc
