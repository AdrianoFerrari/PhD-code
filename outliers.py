import numpy, sys

list = []

f = open(sys.argv[1],'r')
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
fout = open(sys.argv[1]+".avg",'w')
s = str(avg) + "\t" + str(stdev)
fout.write(s)

f.close()
fout.close()
