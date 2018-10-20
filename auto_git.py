import os

f=open('list_git_changes','r')
lines=f.read().splitlines()
f.close()
for i in lines:
  tmp=(i.split())
  if tmp[0]=='deleted:':
    cmd=str('git rm ')+str(tmp[1])
    os.system(cmd)
  if tmp[0]=='modified:':
    cmd=str('git add ')+str(tmp[1])
    os.system(cmd)
