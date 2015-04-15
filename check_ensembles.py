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
import collections
import cPickle
import unittest
import mock


class CFConventions(object):
    '''Some specs from CF-1.6'''
    def __init__(self):
        self.valid_periods = ["millisecs", "msecs", "seconds", "secs",
                              "minutes", "mins", "hours", "hrs", "days",
                              "weeks"]
        self.valid_calendars = ['gregorian', 'standard',
                                'proleptic_gregorian ', 'noleap', '365_day',
                                'noleap', '365_day', '360_day', 'julian']


class TimeCheck(object):
    '''
    Operates on a "datasets", that is a set of netCDF files, whichg
    only differ in the time span covered. Initialize with a path and
    a regular expression that define the files of the dataset.
    '''
    def __init__(self, path, pattern):
        print("Checking path: {}\n         pattern: {}".format(path, pattern))
        self.path = path
        self.pattern = pattern

    def init(self):
        '''
        Manual call to init required because already quite
        some testing happens here.
        '''
        self.files = self.get_files(self.path, self.pattern)
        self.sort_files()
        self.ds = self._mk_ds()
        self.check_has_time()
        # #
        # self.check_has_attr('time', 'units')
        # self.check_has_same_attr('time', 'units')
        # self.check_timeunit_is_days()
        # #
        # self.check_has_attr('time', 'calendar')
        # self.check_has_same_attr('time', 'calendar')
        # self.calendar = self.check_calendar_unit_is_valid()
        # #
        # self.time_info = self.time_get_info()
        # #
        # self.problemdict = {}
        # self.find_interfile_problems()
        # print(self.time_info[0])

    def get_files(self, path, pattern):
        files = [os.path.join(path, f) for f in os.listdir(path)
                 if re.search(pattern, f)]
        if len(files) == 0:
            raise RuntimeError(
                "No files found for dataset with path={} and pattern={}."
                .format(path, pattern))
        else:
            return(files)

    def sort_files(self):
        '''Sorts files according to first date.'''
        starttimes = [nc.Dataset(f).variables['time'][0] for f in self.files]
        filetuples = zip(self.files, starttimes)
        filetuples.sort(key=lambda f: f[1])
        self.files = [f[0] for f in filetuples]

    def _mk_ds(self):
        '''Returns netCDF4 - Objects for all files'''
        return([nc.Dataset(f) for f in self.files])

    def check_has_attr(self, variable, attribute):
        res = [attribute in d.variables[variable].ncattrs() for d in self.ds]
        if not all(res):
            nocal = [f[0] for f in zip(self.files, res) if not f[1]]
            raise RuntimeError(
                "Files without attribute \"{}\" for variable \"{}\":\n{}"
                .format(attribute, variable, nocal))
        else:
            print("All files have attribute \"{}\" for variable \"{}\"."
                  .format(attribute, variable))

    def check_has_same_attr(self, variable, attribute):
        atts = [d.variables[variable].getncattr(attribute) for d in self.ds]
        if len(set(atts)) > 1:
            pp.pprint("ERROR: Not all files have the same {} for {}:"
                      .format(attribute, variable))
            pp.pprint(dict(zip(self.files, atts)))
            raise RuntimeError()
        else:
            print("All files have the same {} for {}: {}"
                  .format(attribute, variable, atts[0]))

    def check_calendar_unit_is_valid(self):
        '''Check for valid (according to CF-1.6) calendar'''
#valcals
        cal = self.ds[0].variables['time'].getncattr('calendar')
        if not (cal in valcals):
            raise RuntimeError("Calendar \"{}\" not valid in CF-1.6".format(cal))
        else:
            print("Calendar \"{}\" is valid in CF-1.6".format(cal))
            return(cal)

    def check_timeunit_is_days(self):
        units = [d.variables['time'].getncattr('units') for d in self.ds]
        if re.match("days", units[0]):
            print("Time units are days. Good")
        else:
            raise RuntimeError('''ERROR: time units not \"days\".
            other units currently not implemented.''')

    def check_has_variable(self, variable):
        '''Checks whether all files have variable "variable"'''
        res = [variable in d.variables.keys() for d in self.ds]
        if not all(res):
            notime = [f[0] for f in zip(self.files, res) if not f[1]]
            raise RuntimeError(
                "Files without variable \"{}\" found: {}"
                .format(variable, notime))
        else:
            print("All files have variable \"{}\".".format(variable))
            return(True)

    def ncdftime2datetime(self, t):
        return(datetime.datetime(*t.timetuple()[0:5]))

    def time_get_info(self):
        tim = [d.variables['time'] for d in self.ds]
        timesteps = [t[:] for t in tim]
        spans = [(t[0], t[-1]) for t in timesteps]
        shouldsteps = [np.arange(s[0], s[1]+1) for s in spans]
        return(zip(tim, timesteps, spans, shouldsteps))
        
    def find_interfile_problems(self):
        continuity = []
        for idx in range(0, len(self.files) - 1):
            continuity.append(
                self.time_info[idx+1][2][0] == self.time_info[idx][2][1] + 1)
        if all(continuity):
            print("No gaps in-between files")
        else:
            print("Gaps between files detected")
            self.gaps = continuity


    def _time_global_check(self):
        globerrors = [list(set(t[1]) ^ set(t[3])) for t in self.time_info]
        anygloberrors = any([len(g) > 0 for g in globerrors])
        return(globerrors if anygloberrors else False)

    def time_check_count(self):
        if self.no_dates_should > self.tim.shape[0]:
            print("dates missing!")
            miss = self.time_find_missing()
            print(miss)
            return(-1)
        elif self.tim.shape[0] > self.no_dates_should:
            print("too much dates!")
            return(1)
        else:
            print("number of dates ok!")
            return(0)

    def time_find_missing(self):
        sshould = set(np.arange(self.tim[0], self.tim[-1]+1))
        sis = set(self.timesteps)
        return(list(sshould - sis))

    def time_check(self):
        sshould = set(np.arange(self.tim[0], self.tim[-1]+1))
        sis = set(self.timesteps)
        return(list(sshould ^ sis))

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

    
class TestTimeCheck(unittest.TestCase):

    def setUp(self):
        self.tc = TimeCheck("/a/path/to/somewhere", ".*hus.nc")

    @mock.patch('check_ensembles.os.listdir')
    def test_get_files(self, mock_os_listdir):
        mock_os_listdir.return_value = ["a_file_2015_hus.nc",
                                        "not_matching_husx.nc"]
        self.assertEqual(self.tc.get_files(self.tc.path, self.tc.pattern),
                         ["/a/path/to/somewhere/a_file_2015_hus.nc"])
        mock_os_listdir.return_value = []
        self.assertRaises(RuntimeError, self.tc.get_files,
                          self.tc.path, self.tc.pattern)

    @mock.patch('check_ensembles.nc.Dataset')
    def test_sort_files(self, mock_ncds):
        class NcVarClass(object):
            def __init__(self, f):
                if f == 'file1':
                    self.variables = {'time': [2, 0]}
                elif f == 'file2':
                    self.variables = {'time': [3, 1]}
                elif f == 'file3':
                    self.variables = {'time': [1, 2]}
                else:
                    self.variables = False

        mock_ncds.side_effect = lambda f: NcVarClass(f)
        self.tc.files = ["file1", "file2", "file3"]
        print(self.tc.files)
        self.tc.sort_files()
        self.assertEqual(self.tc.files, ["file3", "file1", "file2"])

    def test_check_has_variable(self):
        class NcDsClass(object):
            def __init__(self, testvar):
                    self.variables = {'one': 1, testvar: 2, 'three': 3}
        self.tc.files = ["file1", "file2"]
        self.tc.ds = [NcDsClass('testvar'), NcDsClass('testvar')]
        self.assertTrue(self.tc.check_has_variable('testvar'))
        self.tc.ds = [NcDsClass('testvar'), NcDsClass('autrevar')]
        self.assertRaises(RuntimeError, self.tc.check_has_variable, 'testvar')


    # @mock.patch('check_ensembles.nc.Dataset')
    # def test_getattribute(self, mock_ncds):
    #     my_dict = {'time': 10, 'tante': 20}
    #     def mydict(name):
    #         return(my_dict[name])
    #     mock_ncds.return_value.variables.__getitem__.side_effect = (
    #         lambda x: 'hello')

    #     #mock_ncds.variables.__getitem__ = mock.Mock(side_effect=mydict)
    #     res = self.tc.getattribute()
    #     print(res)

    # def getattribute(self):
    #     f = '/net/atmos/data/ENSEMBLES-RCM/A1B/ETHZ/DM/ETHZ-CLM_SCN_HadCM3Q0_DM_25km_2041-2050_vas.nc'
    #     #return(nc.Dataset(f).variables['time'].ncattrs())
    #     return(nc.Dataset(f).variables['time'])

    @mock.patch('check_ensembles.nc.Dataset')    
    def test_check_has_attribute(self, mock_ncds):
        mock_ds = mock.Mock('check.ensembles.TimeCheck.ds')
        mock_ds.__getitem__ = lambda x: 
        
        class Myvar(object):
            def ncattrs():
                return(['att1', 'att2'])
        myvar = Myvar()
        mock_ncds.Dataset.variables.__getitem__.side_effect = myvar
        self.tc.check_has_attr('whatevervar', 'time')
        
    # def check_has_attr(self, variable, attribute):
    #     res = [attribute in d.variables[variable].ncattrs() for d in self.ds]


if __name__ == "__main__":
    unittest.main()
    # path = "/data/ENSEMBLES-RCM/A1B/ETHZ/DM"
    #path = "/data/ENSEMBLES-RCM/A1B/SMHI_BCM/DM"
    # pattern = ".*20[234567]1.*hurs\.nc"
    # tc = TimeCheck(path, pattern)
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
