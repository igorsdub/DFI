import sys
mtsite=[60,62,69,74,104,105,116,154,157,194,216,217,219]
mtsite=[60,69,74,104,105,116,154,157,194,216,217,219]
#mtsite_str=['{0:>4}'.format(i) for i in mtsite]
#print mtsite_str
#f1=open('originID', 'r')
f1=open('originID_struct', 'r')
OriginIDtoNewID={}
line=1
count=1
while line:
    line=f1.readline()
    if not line:
        break
    originid=line[:-1]
    newid='A{0:>4}'.format(count)
    count=count+1
    OriginIDtoNewID[originid]=newid
#    print line[1:-1]
f1.close()
#for key in OriginIDtoNewID.keys():
#    print key
#    print OriginIDtoNewID[key]

mtsite_new=[0]
f2=open(sys.argv[1], 'r')
line=1
count=1
while line:
    line=f2.readline()
    if not line:
        break
    if line[16]!='B':
        originid=line[21:26]
######################print renumbered PDB file####################################
        print line[:21]+OriginIDtoNewID[originid]+line[26:-1]
        resid=int(line[22:26])
	if resid in mtsite and 'CA' in line:
#	    print resid, '->', OriginIDtoNewID[originid], '->', int(OriginIDtoNewID[originid].replace('A',''))
	    mtsite_new.append(int(OriginIDtoNewID[originid].replace('A','')))
f2.close()
mtsite_new.append(int(OriginIDtoNewID[originid].replace('A',''))+1)
print len(mtsite_new), mtsite_new

f3=open('loop.core','r')
content=f3.read()
f3.close()
for i in range(len(mtsite_new)-1):
    if mtsite_new[i]+1<=mtsite_new[i+1]-1:
#        print mtsite_new[i]+1, mtsite_new[i+1]-1
       newcontent=content.replace('DIRECTION', sys.argv[2])
       newcontent=newcontent.replace('START', str(mtsite_new[i]+1))
       newcontent=newcontent.replace('STOP', str(mtsite_new[i+1]-1))
       print newcontent[:-1] 
