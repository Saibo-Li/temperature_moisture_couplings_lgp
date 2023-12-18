#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from sko.GA import GA
import numpy as np
import pandas as pd

# Organize data
df = pd.read_csv(r'samples.csv')
df_single = df[df['name']=='single']
df_double= df[df['name']=='double']
df_triple = df[df['name']=='triple']

# Training sample
frac = 0.7
df_single_training = df_single.sample(frac=frac,axis=0)
df_double_training = df_double.sample(frac=frac,axis=0)
df_triple_training = df_triple.sample(frac=frac,axis=0)
df_training = pd.concat([df_single_training,df_double_training,df_triple_training],axis=0)

# Verification sample
df_single_vali = df_single[~df_single.index.isin(df_single_training.index)]
df_double_vali = df_double[~df_double.index.isin(df_double_training.index)]
df_triple_vali = df_triple[~df_triple.index.isin(df_triple_training.index)]
df_vali = pd.concat([df_single_vali,df_double_vali,df_triple_vali],axis=0)

def optim_a(x,df_training):
    '''
    Parameters to be calibrated
   
    Take abs () as the objective function, that is, close to 0, the higher the accuracy.
    '''
    x1,x2,x3 = x
    
    LGP = a*df_training['single']*satellite_based_scale
    df_training = pd.concat([df_training,pd.DataFrame(LGP,columns=['LGP'])],axis=1)
    
    #
    df_training_single = df_training[df_training['name']=='single']
    df_training_double= df_training[df_training['name']=='double']
    df_training_triple = df_training[df_training['name']=='triple']
    
    #Calculate the sample LGP value of triple and average it
    def label_sub_triple(LGP,site,scale):
		LGP_new = LGP*scale*x3
        return abs(LGP_new-site)
    df_training_triple.loc[:,'dif'] =df_training_triple.apply(lambda x: label_double_sub_triple(x.LGP,x.site,x.satellite_based_scale),axis=1)
	
	#Calculate the sample LGP value of doublee and average it
    def label_sub_double(LGP,site,scale):
		LGP_new = LGP*scale*x2
        return abs(LGP_new-site)

    df_training_double.loc[:,'dif'] =df_training_double.apply(lambda x: label_double_sub_triple(x.LGP,x.site,x.satellite_based_scale),axis=1)
    
    #Calculate the sample LGP value of single and average it
    def label_sub_single(LGP,site,scale):
		LGP_new = LGP*scale*x1
        return abs(LGP_new-site)

    df_training_single.loc[:,'dif']=df_training_single.apply(lambda x: label_triple_sub_single(x.LGP,x.site,x.satellite_based_scale),axis=1)

    return abs(df_training_single.[['dif']].mean()+df_training_double.[['dif']].mean()+df_training_triple.[['dif']].mean())

def optim(x):
    return optim_a(x,df_training)

# size_pop=50, max_iter=200
ga = GA(func=optim, n_dim=1, size_pop=10, max_iter=5, prob_mut=0.1, lb=[-1,-1,-1], ub=[1,1,1], precision=1e-7)

start_time = datetime.datetime.now()
best_x, best_y = ga.run()
print((datetime.datetime.now() - start_time).total_seconds())
print('best_x:', best_x, '\n', 'best_y:', best_y)


