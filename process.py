import sys, glob, numpy

def concatenate(base):
    tf = open(base+'.growth.all','w')
    for f in glob.glob(base+".*.growth"):
        print f
        tf.write(open(f).read())
    tf.close()
        
def get_clean_stats(filename):
    list = []

    f = open(filename,'r')
    for line in f:
        list.append(eval(line))

    points_removed = True

    while points_removed:
        avg = numpy.average(list)
        stdev = numpy.std(list)
    
        points_removed = False
    
        for p in list:
            if abs((p-avg)/stdev) > 2.7:
                print "Point removed: ", p
                list.remove(p)
                points_removed = True
    return [avg,stdev]
