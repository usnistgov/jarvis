import dask_searchcv as dcv
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import f1_score, accuracy_score, roc_auc_score
from sklearn.model_selection import cross_val_predict
from sklearn.feature_selection import VarianceThreshold
from sklearn.ensemble import RandomForestRegressor,GradientBoostingRegressor
from sklearn.model_selection import learning_curve
from math import sqrt
from get_des5 import get_comp_descp
from sklearn.ensemble import RandomForestRegressor
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from sklearn.preprocessing import StandardScaler

from sklearn.pipeline import Pipeline
from sklearn.svm import SVC
from sklearn import linear_model, decomposition, datasets
from sklearn.linear_model import LogisticRegression
from sklearn import preprocessing
from sklearn.svm import LinearSVC
from sklearn.decomposition import PCA
from sklearn.feature_selection import f_regression
import warnings,itertools
warnings.filterwarnings('ignore')
from pymatgen.core.composition import Composition
from pymatgen.phasediagram.entries import PDEntry
from pymatgen.phasediagram.pdmaker import PhaseDiagram
from pymatgen.phasediagram.plotter import PDPlotter
from sklearn.feature_selection import SelectKBest
from monty.json import MontyEncoder, MontyDecoder
from monty.serialization import loadfn, dumpfn
import numpy as np
import json,pickle,os,sys,time
from sklearn.svm import SVR
from sklearn import svm
from pymatgen import Composition
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error,r2_score,mean_squared_error
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.grid_search import GridSearchCV
from sklearn.ensemble.gradient_boosting import GradientBoostingRegressor
from sklearn.datasets import make_hastie_10_2
from sklearn.model_selection import cross_val_score
from numpy import zeros, mean
from sklearn import linear_model, cross_validation, metrics, ensemble
from sklearn import metrics

def plot_learning_curve(estimator=None, title='', X=None, y=None, ylim=None, cv=3,
                        n_jobs=-1, train_sizes=np.linspace(.1, 1.0, 5),fname='fig.png'):
    plt.figure()
    #fname='fig.png'
    plt.title(title)
    if ylim is not None:
        plt.ylim(*ylim)
    plt.xlabel("Training examples")
    plt.ylabel("Score")
    train_sizes, train_scores, test_scores = learning_curve(
        estimator, X, y, cv=cv, n_jobs=n_jobs, train_sizes=train_sizes,scoring='r2')
    train_scores_mean = np.mean(train_scores, axis=1)
    train_scores_std = np.std(train_scores, axis=1)
    test_scores_mean = np.mean(test_scores, axis=1)
    test_scores_std = np.std(test_scores, axis=1)
    plt.grid()

    plt.fill_between(train_sizes, train_scores_mean - train_scores_std,
                     train_scores_mean + train_scores_std, alpha=0.1,
                     color="r")
    plt.fill_between(train_sizes, test_scores_mean - test_scores_std,
                     test_scores_mean + test_scores_std, alpha=0.1, color="g")
    plt.plot(train_sizes, train_scores_mean, 'o-', color="r",
             label="Training score")
    plt.plot(train_sizes, test_scores_mean, 'o-', color="g",
             label="Cross-validation score")

    plt.legend(loc="best")
    fname=str(fname)+str('1.png')
    plt.savefig(fname)
    plt.close()


    data=str(fname)+str('learn_data')
    f=open(data,'w')
    line=str('train_sizes ')+str('train_scores_mean ')+str('train_scores_std ')+str('test_scores_mean ')+str('test_scores_std ')+'\n'
    f.write(line)
    for i,j,k,l,m in zip(train_sizes,train_scores_mean,train_scores_std,test_scores_mean,test_scores_std):
        line=str(i)+str(' ')+str(j)+str(' ')+str(k)+str(' ')+str(l)+str(' ')+str(m)+str('\n')
        f.write(line)
    f.close()

    plt.figure()
    #fname='fig.png'
    plt.title(title)
    if ylim is not None:
        plt.ylim(*ylim)
    plt.xlabel("Training examples")
    plt.ylabel("Score")
    plt.grid()
    plt.fill_between(train_sizes, test_scores_mean - test_scores_std,
                     test_scores_mean + test_scores_std, alpha=0.1, color="g")
    plt.plot(train_sizes, test_scores_mean, 'o-', color="g",
             label="Cross-validation score")
    plt.legend(loc="best")
    fname=str(fname)+str('2.png')
    plt.savefig(fname)
    plt.close()

    #return plt

def train_gb(n_jobs=-1,x=[],y=[],filename='gb_savep',cv=3,cv_val=True,learning_curve_only=False):
        info={}

        X_train, X_test, Y_train, Y_test = train_test_split(x, y, random_state=1, test_size=0.1)
        localtime = time.asctime( time.localtime(time.time()) )
        print localtime
        if cv_val==True:
          localtime = time.asctime( time.localtime(time.time()) )
          print "beofre grid",localtime
          print 'n_jobs',n_jobs
          param_grid={'est__n_estimators':[3000],'est__loss':['ls','huber'],'est__learning_rate':[0.1,0.05,0.01],'est__max_depth':[4,6],'est__min_samples_split':[2,4,6],'est__min_samples_leaf':[2,5,9,17],'est__max_features':[1,0.3,0.1]}
          est= GradientBoostingRegressor()
          pipe=Pipeline([ ("fs", VarianceThreshold(.00001)),("est", est)])
          print "pipekeys",pipe.get_params().keys()
          gs_cv=GridSearchCV(pipe,param_grid,cv=cv,scoring='r2',n_jobs=n_jobs).fit(X_train,Y_train)
          model=gs_cv.best_estimator_
          localtime = time.asctime( time.localtime(time.time()) )
          print "after grid",localtime
          print localtime
          print "cv=",cv 
          #print 'cv results',(gs_cv.grid_scores_)
          print 'cv best score',(gs_cv.best_score_)
          print 'cv best params',(gs_cv.best_params_)
          #print 'cv results',(gs_cv.cv_results_)

          try:
             print 'cv results',(gs_cv.grid_scores_)
             print 'cv results',(gs_cv.cv_results_)
          except:
             pass
        else:
          model= GradientBoostingRegressor()

        if learning_curve_only==True:

          fname=str(filename)+str('Learnonly.png')
          plot_learning_curve(estimator=model,fname=fname,title="SelectKBest",X=x, y=y,n_jobs=n_jobs)
          #plot_learning_curve(estimator=Pipeline([("fs", SelectKBest(f_classif, k=500)), ("est", est)]),fname=fname,title="SelectKBest",X=x, y=y,n_jobs=-1)
        scores = cross_val_score(model, X_train, Y_train, cv=cv,n_jobs=n_jobs,scoring='r2')
        print "Cross-validated scores,mean:", scores,np.mean(np.array(scores))
        model.fit(X_train,Y_train)
        try:
           feature_importance =100*model.steps[1][1].feature_importances_
        except:
           feature_importance =100*model.feature_importances_
           pass 
        predictions=model.predict(X_test)
        plt.scatter(Y_test, predictions)
        fnm=str(filename)+str('PRED.png')
        plt.savefig(fnm)
        plt.close()
        data=str(filename)+str('pred_data')
        f=open(data,'w')
        line=str('y_test ')+str('predictions')+'\n'
        f.write(line)
        for i,j in zip(Y_test, predictions):
         line=str(i)+str(' ')+str(j)+str('\n')
         f.write(line)
        f.close()

        accuracy = r2_score(Y_test, predictions)
        print "Cross-Predicted Accuracy r2:", accuracy
        print "rmse", np.sqrt(mean_squared_error(Y_test, predictions))
        print 'mae',(mean_absolute_error(Y_test, predictions))
        print ('Y_test, predictions start')
        for i,j in zip(Y_test, predictions):
              print i,j
        print ('Y_test, predictions end')
        info['r2_test_pred']=accuracy
        info['rmse']=np.sqrt(mean_squared_error(Y_test, predictions))
        info['mae']=(mean_absolute_error(Y_test, predictions))
        info['len_pred']=len(Y_test)
        for k,v in info.iteritems():
          print k,v

        model.fit(x,y)
        pickle.dump( model, open( filename, "wb" ) )
        fname=str(filename)+str('-GB-full.png')
        yp=model.predict(x)
        fig, ax = plt.subplots()
        ax.set_xlabel('DFT prediction')
        ax.set_ylabel('ML prediction')
        ax.scatter(y,yp, edgecolors=(0, 0, 0))
        plt.savefig(fname)
        plt.close()
        print localtime


        fname=str(filename)+str('-GB-CV.png')
        predicted = cross_val_predict(model, x, y, cv=cv,n_jobs=n_jobs)
        fig, ax = plt.subplots()
        ax.set_xlabel('DFT prediction')
        ax.set_ylabel('ML prediction')
        ax.scatter(y, predicted, edgecolors=(0, 0, 0))
        plt.savefig(fname)
        plt.close()
        fname=str('Learning-')+str(filename)
        plot_learning_curve(estimator=model,fname=fname, title='',X=x,y=y,cv=cv,n_jobs=n_jobs)
        print 'namesd steps'
        #print pipe.steps
        #feature_importance =pipe.named_steps['est'].feature_importances_
        print 'feature_importance',feature_importance
        for ii,i in enumerate(feature_importance):
             print ii,i
        pos = np.arange(len(feature_importance)) + .5
        plt.barh(pos, feature_importance, align='center')
        fname=str(filename)+str('-Imprt.png')
        plt.savefig(fname)
        plt.close()


        print 'actual','predicted for whole data'
        for i,j in zip(y,yp):
             print i,j

def rms_val(y_actual='',y_predicted=''):
    rms = sqrt(mean_squared_error(y_actual, y_predicted))
    return rms

def grid_cv(X,Y):


  parameter_candidates = [
  {'C': [1, 10, 100, 1000], 'kernel': ['linear']},
  {'C': [1, 10, 100, 1000], 'gamma': [0.001, 0.0001], 'kernel': ['rbf']},
  ]
  clf = GridSearchCV(estimator=svm.SVC(), param_grid=parameter_candidates, n_jobs=-1)
  clf.fit(X, Y)

  print('Best score for data1:', clf.best_score_)
  print('Best C:',clf.best_estimator_.C)
  print('Best Kernel:',clf.best_estimator_.kernel)
  print('Best Gamma:',clf.best_estimator_.gamma)
  if clf.best_estimator_.kernel=='linear':
     
      final_clf=svm.SVC(C=clf.best_estimator_.C, kernel='linear').fit(X,y)
  elif clf.best_estimator_.kernel=='rbf':
      final_clf=svm.SVC(C=clf.best_estimator_.C, kernel='rbf',gamma=clf.best_estimator_.gamma).fit(X,y)
  else:
      print ('PROBLEMMMM')
  return final_clf
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False


def uni_en(el='Na'):
   en='na'
   f=open('uni.json','r')
   u=json.load(f,cls=MontyDecoder )
   f.close()
   for i,j in u.iteritems():
      #print i,j
      if str(i)==str(el):
          en=j['en']
    #      if el=='N':
    #         print el,'found',en

   return en

def form_enp(strt='',en=''):
          form_en='na'
          els=strt.composition.elements
          dd=strt.composition.get_el_amt_dict()
          ref=0.0
          for k,v in dd.iteritems():
                 e1=uni_en(k)
                 if e1=='na':
                      ref='na'
                      break
                 else:
                      ref=ref+float(v)*float(e1)
          if isfloat(ref):
              form_en=float(en)-float(ref)
              form_en=round(float(form_en)/float(strt.composition.num_atoms),3)
          return form_en


def join_str(strng=''):
    new=''
    ele=strng.split('-')
    #print ele
    for i,j in enumerate(ele):
       if j!='':
         new=new+str('-')+str(j)
    return new
def convert(xml_file, xml_attribs=True,att='',divide=False):
    with open(xml_file) as f:
        info=[]
        X=[]
        Y=[]
        d = xmltodict.parse(f, xml_attribs=xml_attribs)
        #print d
        ele= str(d['JARVIS-VASP']['elements'])
        el_list=ele.split('-')
        #print 'ellist',el_list
        search=join_str(strng=ele)
        #search= str('-'.join([item for item in el_list]))
        jid=str(d['JARVIS-VASP']['JVID'])
        form=str(d['JARVIS-VASP']['formula'])
        sgp=str(d['JARVIS-VASP']['space-group'])
        func=str(d['JARVIS-VASP']['functional'])
        strt=str(d['JARVIS-VASP']['structure'])
        p=Poscar.from_string(strt)
        nat=p.structure.composition.num_atoms
        #sgp=SpacegroupAnalyzer(p.structure).get_spacegroup_symbol()
        typ=str(d['JARVIS-VASP']['calc_type'])
        beg=str(d['JARVIS-VASP']['brill_eg-ev'])
        mbjeg=str(d['JARVIS-VASP']['mbj_eg-ev'])
        kv=str(d['JARVIS-VASP']['kv-gpa'])
        gv=str(d['JARVIS-VASP']['gv-gpa'])
        diel_x=(d['JARVIS-VASP']['realx_arr'])
        diel_z=(d['JARVIS-VASP']['realz_arr'])
        ref=str(d['JARVIS-VASP']['reference'])
        energy=d['JARVIS-VASP']['energy-ev']
        enp=str(round(float(energy)/float(nat),4))
        #print ele,jid,sgp,typ,enp
        try:
          X=get_comp_descp(struct=p.structure)
          if divide==True:
          #print 'X',len(X)
            X=X/float(np.mean(np.array(X)))
        except:
           pass
        #X=prepare_x(struct=p.structure)
        #print X,type(X),(p.structure.volume)
        if func=='optB88-vDW' and typ=='bulk':
         if att=='form_en':
           fenp=form_enp(p.structure,energy)
           if isfloat(fenp):
           #return X,[fenp]
              Y=[fenp]
         elif att=='energy':
           if isfloat(energy):
              Y=[float(energy)]
         elif att=='nx':
           tmp=diel_x[0]
           if isfloat(tmp):
              nx=float(tmp)**.5
              Y=[float(nx)]
         elif att=='nz':
           tmp=diel_z[0]
           if isfloat(tmp):
              nz=float(tmp)**.5
              Y=[float(nz)]
         elif att=='bg':
           if isfloat(beg):
              Y=[float(beg)]
           #return X,[beg]
         elif att=='mbj_bg':
           if isfloat(mbjeg):
              Y=[float(mbjeg)]
           #return X,[mbjbg]
         elif att=='kv':
           if isfloat(kv):
              Y=[float(kv)]
         elif att=='gv':
           if isfloat(gv):
              Y=[float(gv)]
           #return X,[kv,gv]
         else:
            print 'Not recognized'
            import sys
            sys.exit()
        return X,Y
       
def datacollect(att='epsx',divide=False):
      vers='v25h'
      print 'att=',att,'divide',divide
      tol=0.0
      X=[]
      Y=[]
      count=0
      memb=loadfn('jv_opt.json',cls=MontyDecoder)
      if att=='epsx': 
	for i in memb:
         opgap=i['op_gap']
         epsx=(i['epsx'])
         if isfloat(epsx): 
          if  opgap>=tol and float(epsx)>0:
           #print i['file']
           epsx=sqrt(float(epsx))
           strt=i['final_str']
           x=get_comp_descp(struct=strt)
	   y=epsx
           count=count+1
           if len(x)>=1800 :
             #typ1=str(np.array(list(y)).dtype)
             #typ2=str(np.array(list(x)).dtype)
             #if typ1=='float64' and typ2=='float64':
             X.append(x)
             Y.append(y)
        print "X,Y",len(X),len(Y)
        X=np.array(X).astype(np.float)
        Y=np.array(Y).astype(np.float)
        filename=str(att)+str('_gbsave_')+str(vers)
        train_gb(x=X,y=Y,filename=filename)  
        #filename=str(att)+str('_gbsave_v12a')
        #train_gb(cv=True,x=X,y=Y,filename=filename)  
      if att=='epsz': 
	for i in memb:
         opgap=i['op_gap']
         epsz=(i['epsz'])
         if isfloat(epsz): 
          if opgap>=tol and float(epsz)>0:
           epsz=sqrt(float(epsz))
           #strt=i['data'][0]['contcar']
           strt=i['final_str']
           x=get_comp_descp(struct=strt)
	   y=epsz
           count=count+1
           if len(x)>=1800 :
             #typ1=str(np.array(list(y)).dtype)
             #typ2=str(np.array(list(x)).dtype)
             X.append(x)
             Y.append(y)
        print "X,Y",len(X),len(Y)
        X=np.array(X).astype(np.float)
        Y=np.array(Y).astype(np.float)
        filename=str(att)+str('_gbsave_')+str(vers)
        train_gb(x=X,y=Y,filename=filename)  
        
      if att=='epsy': 
	for i in memb:
         opgap=i['op_gap']
         epsy=(i['epsy'])
         if isfloat(epsy) and opgap>=tol:
          if float(epsy)>0:
           epsy=sqrt(float(epsy))
           #strt=i['data'][0]['contcar']
           strt=i['final_str']
           x=get_comp_descp(struct=strt)
	   y=epsy
           count=count+1
           if len(x)>=1800 :
             #typ1=str(np.array(list(y)).dtype)
             #typ2=str(np.array(list(x)).dtype)
             X.append(x)
             Y.append(y)
        print "X,Y",len(X),len(Y)
        X=np.array(X).astype(np.float)
        Y=np.array(Y).astype(np.float)
        filename=str(att)+str('_gbsave_')+str(vers)
        train_gb(x=X,y=Y,filename=filename)  
      if att=='mepsx': 
	for i in memb:
         mgap=i['mbj_gap']
         mepsx=(i['mepsx'])
         if isfloat(mepsx) and mgap>=tol:
          if  float(mepsx)>0:
           mepsx=sqrt(float(mepsx))
           #strt=i['data'][0]['contcar']
           strt=i['final_str']
           x=get_comp_descp(struct=strt)
	   y=mepsx
           count=count+1
           if len(x)>=1800 :
             #typ1=str(np.array(list(y)).dtype)
             #typ2=str(np.array(list(x)).dtype)
             #if typ1=='float64' and typ2=='float64':
             X.append(x)
             Y.append(y)
        print "X,Y",len(X),len(Y)
        X=np.array(X).astype(np.float)
        Y=np.array(Y).astype(np.float)
        filename=str(att)+str('_gbsave_')+str(vers)
        train_gb(x=X,y=Y,filename=filename)  
        #filename=str(att)+str('_gbsave_v12a')
        #train_gb(cv=True,x=X,y=Y,filename=filename)  
      if att=='mepsy': 
	for i in memb:
         mgap=i['mbj_gap']
         mepsy=(i['mepsy'])
         if isfloat(mepsy) and mgap>=tol:
          if  float(mepsy)>0:
           mepsy=sqrt(float(mepsy))
           #strt=i['data'][0]['contcar']
           strt=i['final_str']
           x=get_comp_descp(struct=strt)
	   y=mepsy
           count=count+1
           if len(x)>=1800 :
             #typ1=str(np.array(list(y)).dtype)
             #typ2=str(np.array(list(x)).dtype)
             #if typ1=='float64' and typ2=='float64':
             X.append(x)
             Y.append(y)
        print "X,Y",len(X),len(Y)
        X=np.array(X).astype(np.float)
        Y=np.array(Y).astype(np.float)
        filename=str(att)+str('_gbsave_')+str(vers)
        train_gb(x=X,y=Y,filename=filename)  
      if att=='mepsz': 
	for i in memb:
         mgap=i['mbj_gap']
         mepsz=(i['mepsz'])
         if isfloat(mepsz) and mgap>=tol :
          if  float(mepsz)>0:
           mepsz=sqrt(float(mepsz))
           #strt=i['data'][0]['contcar']
           strt=i['final_str']
           x=get_comp_descp(struct=strt)
	   y=mepsz
           count=count+1
           if len(x)>=1800 :
             #typ1=str(np.array(list(y)).dtype)
             #typ2=str(np.array(list(x)).dtype)
             #if typ1=='float64' and typ2=='float64':
             X.append(x)
             Y.append(y)
        print "X,Y",len(X),len(Y)
        X=np.array(X).astype(np.float)
        Y=np.array(Y).astype(np.float)
        filename=str(att)+str('_gbsave_')+str(vers)
        train_gb(x=X,y=Y,filename=filename)  
      if att=='op_gap': 
	for i in memb:
         opgap=i['op_gap']
         if opgap!='na' :
           #strt=i['data'][0]['contcar']
           strt=i['final_str']
           x=get_comp_descp(struct=strt)
	   y=opgap
           count=count+1
           if len(x)>=1800 :
             #typ1=str(np.array(list(y)).dtype)
             #typ2=str(np.array(list(x)).dtype)
             X.append(x)
             Y.append(y)
        print "X,Y",len(X),len(Y)
        X=np.array(X).astype(np.float)
        Y=np.array(Y).astype(np.float)
        filename=str(att)+str('_gbsave_')+str(vers)
        train_gb(x=X,y=Y,filename=filename)  
      if att=='mbj_gap': 
	for i in memb:
         mgap=i['mbj_gap']
         opgap=i['op_gap']
         if mgap!='na' and opgap!='na':
          tmp_mbj=round(mgap,2)
          tmp_op=round(opgap,2)
          if tmp_mbj>=tmp_op:
           #strt=i['data'][0]['contcar']
           strt=i['final_str']
           x=get_comp_descp(struct=strt)
	   y=mgap
           count=count+1
           if len(x)>=1800 :
             #typ1=str(np.array(list(y)).dtype)
             #typ2=str(np.array(list(x)).dtype)
             X.append(x)
             Y.append(y)
        print "X,Y",len(X),len(Y)
        X=np.array(X).astype(np.float)
        Y=np.array(Y).astype(np.float)
        filename=str(att)+str('_gbsave_')+str(vers)
        train_gb(x=X,y=Y,filename=filename)  

      if att=='energy': 
	for i in memb:
         energy=i['fin_en']#[0]['final_energy']
         if energy!='na' :
           #strt=i['data'][0]['contcar']
           strt=i['final_str']
           x=get_comp_descp(struct=strt)
	   y=float(energy)/float(strt.composition.num_atoms)
           count=count+1
           if len(x)>=1800 :
             #typ1=str(np.array(list(y)).dtype)
             #typ2=str(np.array(list(x)).dtype)
             X.append(x)
             Y.append(y)
        print "X,Y",len(X),len(Y)
        X=np.array(X).astype(np.float)
        Y=np.array(Y).astype(np.float)
        filename=str(att)+str('_gbsave_')+str(vers)
        train_gb(x=X,y=Y,filename=filename)  
      if att=='kv': 
	for i in memb:
         kv=i['kv']
         if kv!='na' :
           #strt=i['data'][0]['contcar']
           strt=i['final_str']
           x=get_comp_descp(struct=strt)
	   y=kv
           count=count+1
           if len(x)>=1800 :
             #typ1=str(np.array(list(y)).dtype)
             #typ2=str(np.array(list(x)).dtype)
             X.append(x)
             Y.append(y)
        print "X,Y",len(X),len(Y)
        X=np.array(X).astype(np.float)
        Y=np.array(Y).astype(np.float)
        filename=str(att)+str('_gbsave_')+str(vers)
        train_gb(x=X,y=Y,filename=filename)  
      if att=='gv': 
	for i in memb:
         gv=i['gv']
         if gv!='na' :
           #strt=i['data'][0]['contcar']
           strt=i['final_str']
           x=get_comp_descp(struct=strt)
	   y=gv
           count=count+1
           if len(x)>=1800 :
             #typ1=str(np.array(list(y)).dtype)
             #typ2=str(np.array(list(x)).dtype)
             X.append(x)
             Y.append(y)
        print "X,Y",len(X),len(Y)
        X=np.array(X).astype(np.float)
        Y=np.array(Y).astype(np.float)
        filename=str(att)+str('_gbsave_')+str(vers)
        train_gb(x=X,y=Y,filename=filename)  
      if att=='form_en': 
	for i in memb:
         en=i['fin_en']#[0]['final_energy']
         #strt=i['data'][0]['contcar']
         strt=i['final_str']
         fe=form_enp(strt=strt,en=en)
         if fe!='na' :
           #strt=i['data'][0]['contcar']
           x=get_comp_descp(struct=strt)
	   y=fe
           count=count+1
           if len(x)>=1800 :
             #typ1=str(np.array(list(y)).dtype)
             #typ2=str(np.array(list(x)).dtype)
             X.append(x)
             Y.append(y)
        print "X,Y",len(X),len(Y)
        X=np.array(X).astype(np.float)
        Y=np.array(Y).astype(np.float)
        filename=str(att)+str('_gbsave_')+str(vers)
        train_gb(x=X,y=Y,filename=filename)  
def prepare_x(struct=None):
    descriptors = ['row', 'group', 'atomic_mass', 'atomic_radius',   'X']
    #descriptors = ['row', 'group', 'atomic_mass', 'atomic_radius', 'boiling_point', 'melting_point', 'X']
    val=[]
    for i in descriptors:
        d=np.mean(get_pymatgen_descriptor(struct.composition,i))
        val.append(d)
    vat=float(struct.volume)/float(struct.composition.num_atoms)
    val.append(vat)
    return val

def train_data():
        data = loadfn('all_info_0_data.json', cls=MontyDecoder)
        X=[]
        Y=[]
        for i in data:
            #print i['final_str']
            #print i['opt_info']['refrx'][0]
            #print i['opt_info']['refry'][0]
            #print i['opt_info']['refrz'][0]
            #print i['bbandgap']
            #print i['elastic']
            try:
               val=prepare_x(i['final_str'])
               temp=[i['opt_info']['refrx'][0],float(i['bbandgap'].split()[0]),i['elastic'][0][0]]
               typ=str(np.array(list(temp)).dtype)
               print typ
               #add=1
               #for j in temp:
               #    if str(j)=='na':
               #        print 'na'
               #        add=0
               #        #print j
               #if add==1:
               if typ=='float64':

                       X.append(val)
                       Y.append([i['opt_info']['refrx'][0],float(i['bbandgap'].split()[0]),i['elastic'][0][0]])
            except:
                    pass
        linear_regression = LinearRegression()
        X=np.array(X).astype(np.float)
        Y=np.array(Y).astype(np.float)
        linear_regression.fit(X, Y)
        pickle.dump( linear_regression, open( "save.p", "wb" ) )
def predic(att='form_en',plt_fig=False):
  #train_data()
  Y=[]
  X1=[]
  x_arr=[]
  y_arr=[]
  #filname=str(att)+str('_save.p')
  filname=str(att)+str('svr_save.p')
  object_file =pickle.load( open( filname, "rb" )) 
  ll=[]
  for i in glob.glob("*.xml"):
	  ll.append(i)
  ll=['JVASP-54.xml','JVASP-32.xml']
  ll=['JVASP-54.xml','JVASP-32.xml','JVASP-3567.xml','JVASP-140.xml']      
  for j,i in enumerate(ll):
           x,y=convert(i,att=att)
           pred=object_file.predict([x])
           #print i,pred[0],x
           if x!='na' and y!='na' and x!=[] and y!=[] and x!=['na'] and y!=['na']:
             typ1=str(np.array(list(y)).dtype)
             typ2=str(np.array(list(x)).dtype)
             if typ1=='float64' and typ2=='float64':
              if np.all(np.isfinite(x)) and np.all(np.isfinite(y)):

               X=[]
               X.append(x)
               pred=object_file.predict(X)
               x_arr.append(pred[0])
               y_arr.append(y[0])
               if plt_fig==False:
                  Y.append(y)
                  X1.append(x)
               #print pred[0],y[0]
  if plt_fig:
    plt.plot(x_arr,y_arr,'.',label='SVR')#Support Vector Regression
    plt.plot(y_arr,y_arr,'-',label='DFT')
    plt.legend()
    pl_name=str(att)+str('_SVR-DFT.png')
    plt.savefig(pl_name)
    plt.close()
  else:
   print 'Predicted,Actual'
   pred=object_file.predict(X1)
   for i,j in zip(pred,Y):
       print i,j[0]
 
  #print 'Predicted,Actual'
  #pred=object_file.predict(X)
  #for i,j in zip(pred,Y):
  #     print i,j[0]
  #print 'predicted Y',pred
  #print 'actual Y',Y
def comp_mbj_bg():
  ll=[]
  for i in glob.glob("*.xml"):
	  ll.append(i)
  ll=['JVASP-54.xml','JVASP-32.xml','JVASP-535.xml']
  ll=['JVASP-1002.xml']
  #ll=['JVASP-54.xml','JVASP-32.xml','JVASP-3567.xml','JVASP-140.xml']      
  #filname=str('mbj_bgrfm_save.p')
  #filname=str('mbj_bg_save.p')
  #filname=str('mbj_bgpipe_save.p')
  #filname=str('form_enrfm_save.p')
  #filname=str('form_en_save.p')
  filname=str('form_ensvr_save.p')
  object_file =pickle.load( open( filname, "rb" )) 
  object_file1 =pickle.load( open( 'mbj_bgsvr_save.p', "rb" )) 
  print 'MBJpred,beg,mbj'
  for j,i in enumerate(ll):

    with open(i) as f:
        info=[]
        X=[]
        Y=[]
        d = xmltodict.parse(f, xml_attribs=True)
        #print d
        ele= str(d['JARVIS-VASP']['elements'])
        el_list=ele.split('-')
        #print 'ellist',el_list
        search=join_str(strng=ele)
        #search= str('-'.join([item for item in el_list]))
        jid=str(d['JARVIS-VASP']['JVID'])
        form=str(d['JARVIS-VASP']['formula'])
        sgp=str(d['JARVIS-VASP']['space-group'])
        func=str(d['JARVIS-VASP']['functional'])
        strt=str(d['JARVIS-VASP']['structure'])
        p=Poscar.from_string(strt)
        nat=p.structure.composition.num_atoms
        #sgp=SpacegroupAnalyzer(p.structure).get_spacegroup_symbol()
        typ=str(d['JARVIS-VASP']['calc_type'])
        beg=str(d['JARVIS-VASP']['brill_eg-ev'])
        mbjeg=str(d['JARVIS-VASP']['mbj_eg-ev'])
        kv=str(d['JARVIS-VASP']['kv-gpa'])
        gv=str(d['JARVIS-VASP']['gv-gpa'])
        ref=str(d['JARVIS-VASP']['reference'])
        energy=d['JARVIS-VASP']['energy-ev']
        enp=str(round(float(energy)/float(nat),4))
        #print ele,jid,sgp,typ,enp
        warnings.filterwarnings('ignore')
        #X=prepare_x(struct=p.structure)
        X=get_comp_descp(struct=p.structure)
        X=X/float(np.mean(np.array(X)))
        #print X,len(X)
        try:
          X=get_comp_descp(struct=p.structure)
          X=X/float(np.mean(np.array(X)))
          pred=object_file.predict(X)
          #pred1=object_file1.predict(X)
        
          print i,pred[0],form_enp(p.structure,energy)
          #print i,pred[0],beg,mbjeg,len(X),np.mean(np.array(X))
        except:
          pass
        #print pred[0],pred1[0],beg,mbjeg



def pd_plot(system='Al-O'):
    

  att='energy'
  att='form_en'
  filname=str(att)+str('svr_save.p')
  object_file =pickle.load( open( filname, "rb" )) 
  output=[]
  l=system.split('-')
  comb = []
  for i in range(len(l)):
      comb += itertools.combinations(l,i+1)
  comb_list = [ list(t) for t in comb ]
#  comb_list=['Zn','O','Zn-O']
  with MPRester("") as m:
    for i in    comb_list:
        dd='-'.join(i)
        print dd
        data = m.get_data(dd)
        for d in data:
            x = {}
            x['material_id'] = str(d['material_id'])
            structure = m.get_structure_by_material_id(x['material_id'])
            X=get_comp_descp(struct=structure)
            X=X/float(np.mean(np.array(X)))
            pred=object_file.predict(X)
            print structure.composition.reduced_formula,pred[0],d['formation_energy_per_atom'],str(d['material_id'])
            output.append(PDEntry(Composition(structure.composition),float(pred[0])))
  pd = PhaseDiagram(output)
  print output
  plotter = PDPlotter(pd, show_unstable=True)
  name=str(system)+str('_phasediagram.png')
  plotter.write_image(name,image_format="png")
#comp_mbj_bg()

datacollect(att='mbj_gap')
datacollect(att='op_gap')
datacollect(att='form_en')
datacollect(att='epsx')
datacollect(att='epsy')
datacollect(att='epsz')
datacollect(att='mepsx')
datacollect(att='mepsy')
datacollect(att='mepsz')
datacollect(att='kv')
datacollect(att='gv')
datacollect(att='energy')
