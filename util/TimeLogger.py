#------------------------------------------------------------------------------
import json
import collections
import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import xml.etree.ElementTree as ET
#------------------------------------------------------------------------------

# Example: fail if there is more than a 10% variation between min and last
#   check_performance( logfile, 0.1 )
def check_performance( fname, rtol ):
    rtol = float(rtol)
    print "Tolerance: ", rtol
    data = load_json( fname )
    isOkay = True
    for entry, timings in data.iteritems() :
        t = timings["time"]
        tmin = min(t)
        tend = t[ len(t)-1 ]
        rerr = ( tend - tmin ) / tmin
        #print tmin, tend, rerr, type(rtol), type(rerr)
        if rerr > rtol :
            print "failed performance test on ", entry
            isOkay = False
            #raise UserWarning( "error threshold exceeded on " + entry )
    if isOkay :
        print "Passed performance test"
        
    return isOkay

#------------------------------------------------------------------------------
def generate_plots( fname ) :
    data = load_json( fname )
    for entry, timings in data.iteritems() :
        print "generating plots for ", entry
        plot( entry, timings )
#------------------------------------------------------------------------------
        

#------------------------------------------------------------------------------
def plot( name, data ) :
    
    dates = data["date"]
    times = data["time"]

    ddiff = max(dates)-min(dates)
    
    plt.close()
    plt.figure()    
    
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%b-%d %X'))
    
    # set x-axis scale
    if ddiff.days > 60 :
        plt.gca().xaxis.set_major_locator(mdates.MonthLocator())
    elif ddiff.days > 2 :
        plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    else :
        plt.gca().xaxis.set_major_locator(mdates.HourLocator())

    plt.plot( dates, times, 'bo-' )
    plt.gcf().autofmt_xdate()    

    plt.title( name )
    plt.xlabel( "Date" )
    plt.ylabel( "Time (s)" )
    plt.grid(True)
    plt.setp(plt.gca().get_xmajorticklabels(), size=6,rotation=30) 
#    plt.show()
    plt.savefig( name )
    
#------------------------------------------------------------------------------
    

#------------------------------------------------------------------------------
def load_json( fname ) :
    
    jf = open( fname )
    data = json.load( jf )

    # sort by ascending date
    data = collections.OrderedDict(sorted(data.items()))
    
    jf.close()
    
    print "there are ", len(data), " different timing entries\n"

    perf = {}

    for date, timings in data.iteritems() :
#        print  datetime.datetime.strptime(date,"%Y-%b-%d %X")
        for entry, t in timings.iteritems() :
            if entry not in perf:
                perf[entry] = dict( time=[], date=[] )
#            print "adding entry: ", date, " with t=", t, " to ", entry
            perf[entry]["time"].append(float(t))
#            perf[entry]["date"].append(date)
            perf[entry]["date"].append( datetime.datetime.strptime(date,"%Y-%b-%d %X") )
    return perf
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
def load_xml( fname ):
    tree = ET.parse( fname )
    return tree.getroot()
#------------------------------------------------------------------------------

