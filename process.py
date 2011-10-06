import sys, glob, numpy

def concatenate(base):
    tf = open(base+'.force.all','w')
    for f in glob.glob(base+".*.force"):
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
        print avg, stdev

        points_removed = False
    
        for p in list:
            if abs((p-avg)/stdev) > 2.7:
                print "Point removed: ", p
                list.remove(p)
                points_removed = True
    return [avg,stdev]

#Main script
fjd = open('job_data')
fjdout = open('job_data.csv','w')

for line in fjd:
    if line.startswith("R\t"):
        fjdout.write(line)
    else:
        base = line.split("\t")[-1][0:-1]
        print base
        concatenate(base)
        stats = get_clean_stats(base+'.force.all')
        stats_string = [str(stats[0]),str(stats[1])]
        
        fjdout.write(line[0:-1]+"\t"+"\t".join(stats_string)+"\n")

fjdout.close()
