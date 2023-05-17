#!/usr/bin/env python
# coding: utf-8

# # Importing Libraries

# In[101]:


import os
import io
import csv
import sys
import random
import numpy as np
import pandas as pd
from tqdm import tqdm
import importlib
import importlib.resources as resources
from functools import reduce
from itertools import combinations
from collections import defaultdict

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from lime.lime_tabular import LimeTabularExplainer

from signaturizer import Signaturizer
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import Chem, DataStructs
from rdkit.Chem import Crippen
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors
from rdkit.Chem import MACCSkeys
import tensorflow as tf
import pickle
import joblib

import warnings
warnings.filterwarnings('ignore', '.*X does not have valid feature names*', )
warnings.filterwarnings('ignore', '.*DataFrame is highly fragmented*', )
warnings.filterwarnings('ignore', '.*Could not cast to int64*', )
warnings.filterwarnings('ignore', '.*pandas.Int64Index is deprecated*', )
warnings.filterwarnings('ignore', '.*Trying to unpickle estimator*', )
warnings.filterwarnings('ignore', '.*Could not load dynamic library*', )
warnings.filterwarnings('ignore', '.*Ignore above cudart dlerror*', )
warnings.filterwarnings('ignore', '.*tensorflow/stream_executor*', )


print('---- Please read the documentation carefully for information about required input format (OpenBabel generated Canonical SMILES) and different functionalities of the AgeXtend Package. ----')


# # Featurizing and Preprocessing Class

# In[10]:


def Feature_Signaturizer(dat):
    sig_df=pd.DataFrame()
    sig_df['smiles']=dat['smiles'].tolist()
    desc=['A','B','C','D','E']
    print('Performing Descriptor Calculation.')
    for dsc in tqdm(desc):
        for i in range(1,6):
            #print('Performing '+dsc+str(i)+' Descriptor Calculation.')
            sign = Signaturizer(dsc+str(i))
            results = sign.predict(dat['smiles'].tolist())
            #print('Performing '+str(list(pd.DataFrame(results.signature).shape)))
            df=pd.DataFrame(results.signature)
            for clm in list(df.columns):
                df=df.rename(columns={clm:dsc+str(i)+'_'+str(clm)})
            sig_df=pd.concat([sig_df,df],axis = 1)
    return handle_missing_values(sig_df)

def handle_missing_values(Idata):
    print('Processing Missing Values in the generated features.')
    data = Idata.drop(['smiles'],axis=1)
    data = data.replace([np.inf, -np.inf, "", " "], np.nan)
    data = data.replace(["", " "], np.nan)
    for i in data.columns:
        data[i] = data[i].fillna(data[i].mean())
    data['smiles'] = Idata['smiles']
    print('Done')
    return data

def ECFP_featurizer(dat):
    data = pd.DataFrame()
    data['smiles'] = dat['smiles']
    n_bits = 2048  # Number of bits in fixed-length fingerprint
    string = []
    print('Calculating ECFPs.')

    for i in data['smiles']:
        ms = Chem.MolFromSmiles(i)
        if ms is not None:
            fp_bitvect = AllChem.GetMorganFingerprintAsBitVect(ms, 2, nBits=n_bits)
            string.append(fp_bitvect)
        if ms is None:
            print('Found RDKIT error: Skipping a SMILE since it cannot be featurized.')
            string.append(np.nan)
    data['ecfp'] = string
    data.dropna(subset=['ecfp'],inplace=True)
    data.reset_index(drop=True,inplace=True)
    return hallmark_tanimoto_test(data)

def MFP_featurizer(dat):
    testset_bits = pd.DataFrame()
    testset_bits['query_smiles'] = dat['smiles']
    list_bit_vect = []
    print('Calculating MFPs')
    for i in dat['smiles']:
        ms = Chem.MolFromSmiles(i)
        if ms is not None:
            fp_bitvect = MACCSkeys.GenMACCSKeys(ms)
            list_bit_vect.append(fp_bitvect)
        if ms is None:
            print('Found RDKIT error: Skipping a SMILE since it cannot be featurized.')
            list_bit_vect.append(np.nan)
    
    testset_bits['bit_vectors'] = list_bit_vect
    testset_bits.dropna(subset=['bit_vectors'],inplace=True)
    testset_bits.reset_index(drop=True,inplace=True)    
    return binding_data(testset_bits)


# # Lipinsky Filter

# In[103]:


def lipinsky_pass(data):
    lip_df = pd.DataFrame()
    
    Num_HDonors = []
    Num_HAcceptors = []
    Mol_Weight = []
    Mol_logp = []
    
    Lipinsky_Pass = []
    
    for i in data['Query_smiles']:
        failed = []

        mol = Chem.MolFromSmiles(i)
        if mol is None:
            raise Exception('%s is not a valid SMILES string' % i)

        num_hdonors = Lipinski.NumHDonors(mol)
        num_hacceptors = Lipinski.NumHAcceptors(mol)
        mol_weight = Descriptors.MolWt(mol)
        mol_logp = Crippen.MolLogP(mol)

        if num_hdonors > 5:
            Num_HDonors.append('Over 5 H-bond donors, found %s' % num_hdonors)
            failed.append('failed num hdonors')
        else:
            Num_HDonors.append('Found %s H-bond donors' % num_hdonors)

        if num_hacceptors > 10:
            Num_HAcceptors.append('Over 10 H-bond acceptors, found %s' \
            % num_hacceptors)
            failed.append('failed num hacceptors')
        else:
            Num_HAcceptors.append('Found %s H-bond acceptors' % num_hacceptors)

        if mol_weight >= 500:
            Mol_Weight.append('Molecular weight over 500, calculated %s'\
            % mol_weight)
            failed.append('failed molecular weight')
        else:
            Mol_Weight.append('Molecular weight: %s' % mol_weight)

        if mol_logp >= 5:
            Mol_logp.append('Log partition coefficient over 5, calculated %s' \
            % mol_logp)
            failed.append('failed log partition coefficient')
        else:
            Mol_logp.append('Log partition coefficient: %s' % mol_logp)
        if failed :
            Lipinsky_Pass.append('False')
        else :
            Lipinsky_Pass.append('True')
            
    lip_df['Num_HDonors'] = Num_HDonors
    lip_df['Num_HAcceptors'] = Num_HAcceptors
    lip_df['Mol_Weight'] = Mol_Weight
    lip_df['Mol_logp'] = Mol_logp
    lip_df['Lipinsky_Pass'] = Lipinsky_Pass
    
    lip_df = pd.concat([data,lip_df],axis=1)
    print('Done with druggability information.')
    return lip_df


# # Getting Predictions Class

# In[12]:


def get_predictions(probs,thresh):
    preds = []
    for prob in probs:
        if prob >= thresh:
            preds.append(1)    
        else:
            preds.append(0)
    return preds


# # Signaturizer Anti Aging Model Class

# In[61]:


class model_sign_aging:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.feature_names_in_
        return data[F_names]
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred        
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','AIM_D_Signaturizer_model_svm_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model, test.drop(['smiles'],axis=1))
        probs = model.predict_proba(test_filtered)    
        preds = self.get_labels(probs)
        return probs,preds


# # HallMarks Of Aging Classes

# ## Altered Intercellular Communication

# In[62]:


class model_aic:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.feature_names_in_
        return data[F_names]
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred        
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','HOLY_AIC_model_svm_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model,test.drop(['smiles'],axis=1))
        probs = model.predict_proba(test_filtered)    
        preds = self.get_labels(probs)
        return probs,preds


# ## Cellular Senescence

# In[63]:


class model_cs:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.feature_names_in_
        return data[F_names]
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred        
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','HOLY_CS_model_svm_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model,test.drop(['smiles'],axis=1))
        probs = model.predict_proba(test_filtered)    
        preds = self.get_labels(probs)
        return probs,preds


# ## Deregulated Nutrient Sensing

# In[64]:


class model_dns:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.feature_names_in_
        return data[F_names]
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred        
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','HOLY_DNS_model_svm_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model,test.drop(['smiles'],axis=1))
        probs = model.predict_proba(test_filtered)    
        preds = self.get_labels(probs)
        return probs,preds


# ## Epigenetic Alterations

# In[90]:


class model_ea:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.feature_names_in_
        return data[F_names]
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred        
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','HOLY_EA_model_svm_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model,test.drop(['smiles'],axis=1))
        probs = model.predict_proba(test_filtered)    
        preds = self.get_labels(probs)
        return probs,preds


# ## Genomic Instability

# In[66]:


class model_gi:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.feature_names_in_
        return data[F_names]
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred        
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','HOLY_GI_model_RF_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model,test.drop(['smiles'],axis=1))
        probs = model.predict_proba(test_filtered)    
        preds = self.get_labels(probs)
        return probs,preds


# ## Loss of Proteostasis

# In[67]:


class model_lp:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.feature_names_in_
        return data[F_names]
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred        
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','HOLY_LP_model_svm_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model,test.drop(['smiles'],axis=1))
        probs = model.predict_proba(test_filtered)    
        preds = self.get_labels(probs)
        return probs,preds


# ## Mitochondrial Dysfunction

# In[68]:


class model_md:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.get_booster().feature_names
        return data[F_names]
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred        
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','HOLY_MD_model_XGB_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model,test.drop(['smiles'],axis=1))
        probs = model.predict_proba(test_filtered)    
        preds = self.get_labels(probs)
        return probs,preds


# ## Stem Cell Exhaustion

# In[69]:


class model_sce:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.feature_names_in_
        return data[F_names]
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred        
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','HOLY_SCE_model_et_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model,test.drop(['smiles'],axis=1))
        probs = model.predict_proba(test_filtered)    
        preds = self.get_labels(probs)
        return probs,preds


# ## Telomere Attrition

# In[70]:


class model_ta:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.feature_names_in_
        return data[F_names]
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred        
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','HOLY_TA_model_ET_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model,test.drop(['smiles'],axis=1))
        probs = model.predict_proba(test_filtered)    
        preds = self.get_labels(probs)
        return probs,preds


# # Toxicity Models

# In[71]:


class model_AMES:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.feature_names_in_
        return data[F_names]
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred        
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','HOLY_AMES_model_MLP_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model,test.drop(['smiles'],axis=1))
        probs = model.predict_proba(test_filtered)    
        preds = self.get_labels(probs)
        return probs,preds


# In[72]:


class model_BBB:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.feature_names_in_
        return data[F_names]
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred        
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','HOLY_BBB_model_MLPBoruta_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model,test.drop(['smiles'],axis=1))
        probs = model.predict_proba(test_filtered)    
        preds = self.get_labels(probs)
        return probs,preds


# In[73]:


class model_CYP1A2:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.feature_names_in_
        return data[F_names]
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred        
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','HOLY_CYP1A2_model_MLP_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model,test.drop(['smiles'],axis=1))
        probs = model.predict_proba(test_filtered)    
        preds = self.get_labels(probs)
        return probs,preds


# In[74]:


class model_CYP2C19:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.feature_names_in_
        return data[F_names]
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred        
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','HOLY_CYP2C19_model_SVM_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model,test.drop(['smiles'],axis=1))
        probs = model.predict_proba(test_filtered)    
        preds = self.get_labels(probs)
        return probs,preds


# In[75]:


class model_CYP2C9:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.get_booster().feature_names
        return data[F_names]
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred        
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','HOLY_CYP2C9_model_XGB_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model,test.drop(['smiles'],axis=1))
        probs = model.predict_proba(test_filtered)    
        preds = self.get_labels(probs)
        return probs,preds


# In[76]:


class model_CYP2D6:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.feature_names_in_
        return data[F_names]
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred        
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','HOLY_CYP2D6_model_SVM_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model,test.drop(['smiles'],axis=1))
        probs = model.predict_proba(test_filtered)    
        preds = self.get_labels(probs)
        return probs,preds


# In[77]:


class model_CYP3A4:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.feature_names_in_
        return data[F_names]
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred        
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','HOLY_CYP3A4_model_SVM_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model,test.drop(['smiles'],axis=1))
        probs = model.predict_proba(test_filtered)    
        preds = self.get_labels(probs)
        return probs,preds


# In[78]:


class model_DILI:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.feature_names_in_
        return data[F_names]
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred        
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','HOLY_DILI_model_SVM_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model,test.drop(['smiles'],axis=1))
        probs = model.predict_proba(test_filtered)    
        preds = self.get_labels(probs)
        return probs,preds


# In[79]:


class model_Hepato:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.feature_names_in_
        return data[F_names]
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred        
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','HOLY_Hepato_model_SVM_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model,test.drop(['smiles'],axis=1))
        probs = model.predict_proba(test_filtered)    
        preds = self.get_labels(probs)
        return probs,preds


# In[80]:


class model_hERG:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.feature_names_in_
        return data[F_names]
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred        
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','HOLY_hERG_model_SVM_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model, test.drop(['smiles'],axis=1))
        probs = model.predict_proba(test_filtered)    
        preds = self.get_labels(probs)
        return probs,preds


# In[81]:


class model_HLM:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.get_booster().feature_names
        return data[F_names]
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred        
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','HOLY_HLM_model_XGB_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model, test.drop(['smiles'],axis=1))
        probs = model.predict_proba(test_filtered)    
        preds = self.get_labels(probs)
        return probs,preds


# In[82]:


class model_MMP:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.get_booster().feature_names
        return data[F_names]
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred        
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','HOLY_MMP_model_XGB_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model,test.drop(['smiles'],axis=1))
        probs = model.predict_proba(test_filtered)    
        preds = self.get_labels(probs)
        return probs,preds


# In[83]:


class model_PGPInh:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.feature_names_in_
        return data[F_names]
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred        
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','HOLY_PGPInh_model_RF_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model,test.drop(['smiles'],axis=1))
        probs = model.predict_proba(test_filtered)    
        preds = self.get_labels(probs)
        return probs,preds


# In[84]:


class model_PGPSub:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.feature_names_in_
        return data[F_names]
    def get_labels(self,pred_test): #Getting discrete labels from probability values    
        test_pred = []        
        for i in range(pred_test.shape[0]):
            if(pred_test[i][0]>pred_test[i][1]):
                test_pred.append(0)
            else:
                test_pred.append(1)
        return test_pred        
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','HOLY_PGPSub_model_SVM_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model,test.drop(['smiles'],axis=1))
        probs = model.predict_proba(test_filtered)    
        preds = self.get_labels(probs)
        return probs,preds


# In[85]:


class model_MRTD:
    def __init__(self,test):        
        self.test = test
    def extract_feature(self,model,data):
        F_names = model.feature_names_in_
        return data[F_names]
    def test_model(self):
        test = self.test
        with resources.open_binary('AgeXtend','HOLY_MRTD_model_DTR_HPTuned_fitted.pkl') as fp:
                model = fp.read()
        model = joblib.load(io.BytesIO(model))
        test_filtered = self.extract_feature(model,test.drop(['smiles'],axis=1))
        preds = model.predict(test_filtered)    
        return preds


# In[42]:

def binding_data(testset):
    with resources.open_binary('AgeXtend','BINDINGDB_db_MFPs.pkl') as fp:
        db_file = fp.read()
    bindingdb_file = pickle.load(io.BytesIO(db_file))

    with resources.open_binary('AgeXtend','BINDINGDB_MFPs.pkl') as fp:
        mfps_file = fp.read()
    bindingdb_mfps = pickle.load(io.BytesIO(mfps_file))
    
    print('Finding targets from BindingDB database. This may take a while. . .')
    sims = defaultdict(dict)
    x = 0
    for i in testset.iterrows():
        for y in bindingdb_mfps.iterrows():
            name = str(x) + '_' + str(y[0])
            sims[str(x)][name] = DataStructs.TanimotoSimilarity(DataStructs.ExplicitBitVect(y[1][0]),i[1][1])
        x += 1
    final_df = pd.DataFrame()
    for i in sims:
        values = list(sims[i].values())
        max_value = np.max(list(sims[i].values()))
        ind_list = [k for k,v in sims[i].items() if v == max_value]
        indices = []
        for t in ind_list[0:3]:
            indices.append(t.split('_')[1])
        temp=pd.DataFrame()
        temp = bindingdb_file.loc[[int(x) for x in indices]]
        temp['Query_smiles'] = testset['query_smiles'][int(i)]
        temp['Tanimoto_Similarity'] = max_value
        final_df = pd.concat([temp,final_df],axis=0,ignore_index=True)
    return lipinsky_pass(final_df)


# In[96]:

def hallmark_tanimoto_test(test_df):
    with resources.open_binary('AgeXtend','HallMark_Data_ECFPs.pkl') as fp:
        hm_ecfps_file = fp.read()
    hallmarks_ecfp = pickle.load(io.BytesIO(hm_ecfps_file))

    split_df_by_prop = hallmarks_ecfp.groupby('HallMark_MoleculeProperty')
    subcat = list(split_df_by_prop.median().index)

    print('Running Hallmarks Tanimoto Test.')
    sims = defaultdict(dict)
    final_df = pd.DataFrame()

    for p in test_df.iterrows():
        temp_df = pd.DataFrame(index=[0],columns=['Query_Smile'])
        temp_df.loc[[0],['Query_Smile']] = p[1][0]

        for i in tqdm(range(len(subcat))):
            x = 0
            for y in (list(split_df_by_prop)[i])[1].iterrows():
                name = str(y[0])
                sims[subcat[i]][y[1][0]] = DataStructs.TanimotoSimilarity(y[1][2],p[1][1])
            all_sim_values = list(sims[subcat[i]].values())
            max_value = np.max(all_sim_values)
            index_smile = [k for k,v in sims[subcat[i]].items() if v == max_value][0]

            temp_df[subcat[i]+'_TS'] = max_value
            temp_df[subcat[i]+'_smile'] = index_smile
            x += 1
        final_df = pd.concat([final_df,temp_df],axis=0)
    print('Done Hallmarks Tanimoto Test Results.')
    return final_df


# # Calculating Prediction Class

# In[88]:


def AgeXtend_Predictions(Sig_data):
    
    ####################################### Output dataframes #########################################################################
    anti_aging_preds = pd.DataFrame(columns=['smiles','Anti_Aging_Status','Anti_Aging_Prob'])
    anti_aging_preds['smiles'] = Sig_data['smiles']
    
    ################################## Signaturizer Anti Aging Model #################################################################
    print('Onto Anti-Aging Predictions.')    
    m0 = model_sign_aging(Sig_data)
    probs,preds = m0.test_model()
    anti_aging_preds['Anti_Aging_Prob'] = probs[:,1]
    anti_aging_preds['Anti_Aging_Status'] = preds 
    
#     anti_aging_preds.to_csv('anti_aging_prediction.csv',index=False)
    print('Done')
    return anti_aging_preds


# In[39]:


def AgeXtend_HallMark_Predictions(Sig_data, hallmark_probs):
    
    ####################################### Output dataframes #########################################################################
    anti_aging_preds = pd.DataFrame(columns=['smiles','Anti_Aging_Status','Anti_Aging_Prob'])
    
    predictions = pd.DataFrame(columns=['smiles','Anti_Aging_Status','AIC_Prediction_Status','CS_Prediction_Status',
                                        'DNS_Prediction_Status','EA_Prediction_Status','GI_Prediction_Status','LP_Prediction_Status',
                                        'MD_Prediction_Status','SCE_Prediction_Status','TA_Prediction_Status'])
    
    probabilities = pd.DataFrame(columns=['smiles','Anti_Aging_Prob','AIC_Prob1','CS_Prob1','DNS_Prob1',
                                          'EA_Prob1','GI_Prob1','LP_Prob1','MD_Prob1','SCE_Prob1','TA_Prob1'])

    anti_aging_preds['smiles'] = Sig_data['smiles']
    
    ################################## Signaturizer Anti Aging Model #################################################################
    print('Onto Anti-Aging Predictions.')    
    m0 = model_sign_aging(Sig_data)
    probs,preds = m0.test_model()
    anti_aging_preds['Anti_Aging_Prob'] = probs[:,1]
    anti_aging_preds['Anti_Aging_Status'] = preds    
       
    
    ################################## Hallmarks of Aging Models ####################################################################
    hoa_input_mols = (anti_aging_preds[anti_aging_preds['Anti_Aging_Status']==1]).reset_index(drop=True,inplace=False)
    hoa_input_df = Sig_data[Sig_data['smiles'].isin(hoa_input_mols['smiles'])].reset_index(drop=True,inplace=False)
    
    predictions['smiles'] = hoa_input_mols['smiles']
    probabilities['smiles'] = hoa_input_mols['smiles']
    
    predictions['Anti_Aging_Status'] = hoa_input_mols['Anti_Aging_Status']
    probabilities['Anti_Aging_Prob'] = hoa_input_mols['Anti_Aging_Prob']
    
    print('Onto Hallmarks of aging predictions.')
    m1 = model_aic(hoa_input_df)
    m2 = model_cs(hoa_input_df)
    m3 = model_dns(hoa_input_df)
    m4 = model_ea(hoa_input_df)
    m5 = model_gi(hoa_input_df)   
    m6 = model_lp(hoa_input_df)
    m7 = model_md(hoa_input_df)
    m8 = model_sce(hoa_input_df)
    m9 = model_ta(hoa_input_df)
    
    probs,preds = m1.test_model()
    probabilities['AIC_Prob1'] = probs[:,1]
    predictions['AIC_Prediction_Status'] = preds    
    probs,preds = m2.test_model()
    probabilities['CS_Prob1'] = probs[:,1]
    predictions['CS_Prediction_Status'] = preds    
    probs,preds = m3.test_model()    
    probabilities['DNS_Prob1'] = probs[:,1]
    predictions['DNS_Prediction_Status'] = preds
    probs,preds = m4.test_model()    
    probabilities['EA_Prob1'] = probs[:,1]
    predictions['EA_Prediction_Status'] = preds 
    probs,preds = m5.test_model()    
    probabilities['GI_Prob1'] = probs[:,1]
    predictions['GI_Prediction_Status'] = preds      
    probs,preds = m6.test_model()
    probabilities['LP_Prob1'] = probs[:,1]
    predictions['LP_Prediction_Status'] = preds        
    probs,preds = m7.test_model() 
    probabilities['MD_Prob1'] = probs[:,1]
    predictions['MD_Prediction_Status'] = preds
    probs,preds = m8.test_model() 
    probabilities['SCE_Prob1'] = probs[:,1]
    predictions['SCE_Prediction_Status'] = preds
    probs,preds = m9.test_model() 
    probabilities['TA_Prob1'] = probs[:,1]
    predictions['TA_Prediction_Status'] = preds
    
    print('Done with AgeXtend predictions and prediction probabilities')
    
    if hallmark_probs == True:
        return anti_aging_preds,predictions,probabilities
    else:
        return anti_aging_preds,predictions,''

# In[99]:


def AgeXtend_Tox_Predictions(Sig_data, hallmark_tox_probs):
    
    ####################################### Output dataframes #########################################################################
    anti_aging_preds = pd.DataFrame(columns=['smiles','Anti_Aging_Status','Anti_Aging_Prob'])
    
    predictions = pd.DataFrame(columns=['smiles','Anti_Aging_Status','AIC_Prediction_Status','CS_Prediction_Status',
                                        'DNS_Prediction_Status','EA_Prediction_Status','GI_Prediction_Status','LP_Prediction_Status',
                                        'MD_Prediction_Status','SCE_Prediction_Status','TA_Prediction_Status'])
    
    probabilities = pd.DataFrame(columns=['smiles','Anti_Aging_Prob','AIC_Prob1','CS_Prob1','DNS_Prob1',
                                          'EA_Prob1','GI_Prob1','LP_Prob1','MD_Prob1','SCE_Prob1','TA_Prob1'])

    anti_aging_preds['smiles'] = Sig_data['smiles']
    
    ################################## Signaturizer Anti Aging Model #################################################################
    print('Onto Anti-Aging Predictions.')
    m0 = model_sign_aging(Sig_data)
    probs,preds = m0.test_model()
    anti_aging_preds['Anti_Aging_Prob'] = probs[:,1]
    anti_aging_preds['Anti_Aging_Status'] = preds    
       
    
    ################################## Hallmarks of Aging Models ####################################################################
    hoa_input_mols = (anti_aging_preds[anti_aging_preds['Anti_Aging_Status']==1]).reset_index(drop=True,inplace=False)
    hoa_input_df = Sig_data[Sig_data['smiles'].isin(hoa_input_mols['smiles'])].reset_index(drop=True,inplace=False)
    
    predictions['smiles'] = hoa_input_mols['smiles']
    probabilities['smiles'] = hoa_input_mols['smiles']
    
    predictions['Anti_Aging_Status'] = hoa_input_mols['Anti_Aging_Status']
    probabilities['Anti_Aging_Prob'] = hoa_input_mols['Anti_Aging_Prob']

    print('Onto Hallmarks of aging predictions.')
    
    m1 = model_aic(hoa_input_df)
    m2 = model_cs(hoa_input_df)
    m3 = model_dns(hoa_input_df)
    m4 = model_ea(hoa_input_df)
    m5 = model_gi(hoa_input_df)   
    m6 = model_lp(hoa_input_df)
    m7 = model_md(hoa_input_df)
    m8 = model_sce(hoa_input_df)
    m9 = model_ta(hoa_input_df)
    
    probs,preds = m1.test_model()
    probabilities['AIC_Prob1'] = probs[:,1]
    predictions['AIC_Prediction_Status'] = preds    
    probs,preds = m2.test_model()
    probabilities['CS_Prob1'] = probs[:,1]
    predictions['CS_Prediction_Status'] = preds    
    probs,preds = m3.test_model()    
    probabilities['DNS_Prob1'] = probs[:,1]
    predictions['DNS_Prediction_Status'] = preds
    probs,preds = m4.test_model()    
    probabilities['EA_Prob1'] = probs[:,1]
    predictions['EA_Prediction_Status'] = preds 
    probs,preds = m5.test_model()    
    probabilities['GI_Prob1'] = probs[:,1]
    predictions['GI_Prediction_Status'] = preds      
    probs,preds = m6.test_model()
    probabilities['LP_Prob1'] = probs[:,1]
    predictions['LP_Prediction_Status'] = preds        
    probs,preds = m7.test_model() 
    probabilities['MD_Prob1'] = probs[:,1]
    predictions['MD_Prediction_Status'] = preds
    probs,preds = m8.test_model() 
    probabilities['SCE_Prob1'] = probs[:,1]
    predictions['SCE_Prediction_Status'] = preds
    probs,preds = m9.test_model() 
    probabilities['TA_Prob1'] = probs[:,1]
    predictions['TA_Prediction_Status'] = preds
    
    print('Onto toxicity predictions')
    
    m11 = model_AMES(hoa_input_df)
    probs,preds = m11.test_model()
    probabilities['AMES_Prob1'] = probs[:,1]
    predictions['AMES_Preds'] = preds
    
    m12 = model_BBB(hoa_input_df)
    probs,preds = m12.test_model()
    probabilities['BBB_Prob1'] = probs[:,1]
    predictions['BBB_Preds'] = preds
    
    m13 = model_CYP1A2(hoa_input_df)
    probs,preds = m13.test_model()
    probabilities['CYP1A2_Prob1'] = probs[:,1]
    predictions['CYP1A2_Preds'] = preds
    
    m14 = model_CYP2C19(hoa_input_df)
    probs,preds = m14.test_model()
    probabilities['CYP2C19_Prob1'] = probs[:,1]
    predictions['CYP2C19_Preds'] = preds
    
    m15 = model_CYP2C9(hoa_input_df)
    probs,preds = m15.test_model()
    probabilities['CYP2C9_Prob1'] = probs[:,1]
    predictions['CYP2C9_Preds'] = preds
    
    m16 = model_CYP2D6(hoa_input_df)
    probs,preds = m16.test_model()
    probabilities['CYP2D6_Prob1'] = probs[:,1]
    predictions['CYP2D6_Preds'] = preds
    
    m17 = model_CYP3A4(hoa_input_df)
    probs,preds = m17.test_model()
    probabilities['CYP3A4_Prob1'] = probs[:,1]
    predictions['CYP3A4_Preds'] = preds
    
    m18 = model_DILI(hoa_input_df)
    probs,preds = m18.test_model()
    probabilities['DILI_Prob1'] = probs[:,1]
    predictions['DILI_Preds'] = preds
    
    m19 = model_Hepato(hoa_input_df)
    probs,preds = m19.test_model()
    probabilities['Hepato_Prob1'] = probs[:,1]
    predictions['Hepato_Preds'] = preds
    
    m20 = model_hERG(hoa_input_df)
    probs,preds = m20.test_model()
    probabilities['hERG_Prob1'] = probs[:,1]
    predictions['hERG_Preds'] = preds
    
    m21 = model_HLM(hoa_input_df)
    probs,preds = m21.test_model()
    probabilities['HLM_Prob1'] = probs[:,1]
    predictions['HLM_Preds'] = preds
    
    m22 = model_MMP(hoa_input_df)
    probs,preds = m22.test_model()
    probabilities['MMP_Prob1'] = probs[:,1]
    predictions['MMP_Preds'] = preds
    
    m23 = model_PGPInh(hoa_input_df)
    probs,preds = m23.test_model()
    probabilities['P-gp_Inhibitor_Prob1'] = probs[:,1]
    predictions['P-gp_Inhibitor_Preds'] = preds
    
    m24 = model_PGPSub(hoa_input_df)
    probs,preds = m24.test_model()
    probabilities['P-gp_Substrate_Prob1'] = probs[:,1]
    predictions['P-gp_Substrate_Preds'] = preds
    
    m25 = model_MRTD(hoa_input_df)
    preds = m25.test_model()
    probabilities['MRTD (mg/day)'] = preds
    predictions['MRTD (mg/day)'] = preds
    print('Done with AgeXtend-Hallmark Toxicity predictions.')
    
    if hallmark_tox_probs == True:
        return anti_aging_preds,predictions,probabilities
    else:
        return anti_aging_preds,predictions,''

# In[43]:
    

def featurize(query_smiles=[]):
    print('Reading file. Please make sure the input SMILES are OpenBabel generated Canonical SMILES.')
    if len(query_smiles)==0:
        raise ValueError("Error: Provided query list is empty")
    else:
        input_data = pd.DataFrame()
        input_data['smiles'] = query_smiles
        print('Featurizing input SMILES.')
        smile_features_df = Feature_Signaturizer(input_data)
        return smile_features_df
          
def predict(Sig_Input, probs=False, HC=False, TS=False, BDL=False):
        result_list ={}
        if HC :
                anti_aging_preds,predictions,probabilities = AgeXtend_Tox_Predictions(Sig_Input,probs)
                result_list['Anti_Aging_Prediction']=anti_aging_preds
                result_list['HallMarks_Toxicity_Status']=predictions
                result_list['HallMarks_Toxicity_Probabilities']=probabilities
        else:
                anti_aging_preds,predictions,probabilities = AgeXtend_HallMark_Predictions(Sig_Input,probs)
                result_list['Anti_Aging_Prediction']=anti_aging_preds
                result_list['HallMarks_Status']=predictions
                result_list['HallMarks_Probabilities']=probabilities
        if TS:
                hallmark_tanimoto = ECFP_featurizer(predictions)
                result_list['HallMark_like_response']=hallmark_tanimoto 
        if BDL:
                bindingdb_druggability = MFP_featurizer(predictions)
                result_list['Druggability_and_Potential_Targets']=bindingdb_druggability 
        return result_list
