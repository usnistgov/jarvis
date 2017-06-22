import os,glob
all_folds=[]
count=0
record=[]
for file in glob.glob("*-*/*_*"):

    if os.path.isdir(file) :
       path=str(os.getcwd())+str("/")+str(file)
       cwd=str(os.getcwd())
       os.chdir(path)
       f_name=str('JARVIS-ID')
       if os.path.isfile(f_name):
         # print os.getcwd()
          f=open(f_name,'r')
          for rline in f:
              r=int(str((str(rline)).split('\n')[0]).split('JVASP-')[1])
              if r>count:
                 count=r
              record.append(r)
          f.close()
       os.chdir(cwd)
print count
for file in glob.glob("*-*/*_*"):

    if os.path.isdir(file) :
       path=str(os.getcwd())+str("/")+str(file)
       cwd=str(os.getcwd())
       os.chdir(path)
       f_name=str('JARVIS-ID')
       if not os.path.isfile(f_name):
         # print os.getcwd()
          f=open(f_name,'w')
          count=count+1
          if count in record:
             print "problem in ", os.getcwd()
          JARVISID=str('JVASP-')+str(count)+'\n'
          f.write(JARVISID)
          f.close()
       os.chdir(cwd)
"""
       if not os.path.isfile(f_name):
          f=open(f_name,'w')
          count=count+1
          JARVISID=str('JVASP-')+str(count)+'\n'
          f.write(JARVISID)
          f.close()
          os.chdir(cwd)
       else:
          f=open(f_name,'r')
          for rline in f:
              r=int(str((str(rline)).split('\n')[0]).split('JVASP-')[1])
              if r>count:
                 print "BIG PROBLEM=",path
                 import sys
                 f.close()
                 sys.exit()
          
"""
