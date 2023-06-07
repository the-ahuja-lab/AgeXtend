# AgeXtend

### Introduction
Finding the Anti-Aging Potential of chemical compounds

## Environment Setup

**Major dependencies**
1. [Signaturizer(v1.1.11)](https://gitlabsbnb.irbbarcelona.org/packages/signaturizer)
2. [RDKit(v2022.3.1)](https://www.rdkit.org/)
3. [LIME](https://github.com/marcotcr/lime)

**Minor dependencies**
1. os
2. [scikit-learn v0.23.2](https://scikit-learn.org/stable/whats_new/v1.0.html)
3. [xgboost v1.5.1](https://github.com/dmlc/xgboost)
4. [pandas](https://pandas.pydata.org/)
5. [numpy](https://numpy.org)
6. [tqdm](https://tqdm.github.io)
7. [joblib](https://pypi.org/project/joblib/)
8. [matplotlib](https://pypi.org/project/matplotlib/)
9. [seaborn](https://seaborn.pydata.org/)
10. [importlib](https://pypi.org/project/importlib/)
11. [importlib-resources v5.7.1](https://github.com/python/importlib_resources)


**Quick setup**

The environment setup can be accomplished dierctly by the follwing command
```
$ pip install -r requirement.txt
```
The **requirement.txt** can be found [here](https://github.com/the-ahuja-lab/AgeXtend/blob/main/env/requirement.txt)


## How to use AgeXtend?

### Installation using pip 
```
$ pip install -i https://test.pypi.org/simple/ AgeXtend
```

## Modules
AgeXtend supports 2 distinct modules:<br/>
1. Predictor
2. Browser

### Predictor

To get predictions for Anti-Aging properties:<br/>
```
>>> from AgeXtend import Predictor
```
Prepare a list of canonical SMILES (Openbabel generated) strings
```
>>> smiles = ['ClCC=C', 'C=CCOC(=O)CC(C)C'] 
```
Create AgeXtend type object of any query compunds of interest (i.e. list of canonical SMILES)
```
>>> agex_obj = Predictor.featurize(smiles)
```
Use the AgeXtend object for Anti-Aging property predictions
```
>>> output = Predictor.predict(agex_obj)
```
To get the list of result dataframes produced as output
```
>>> output.keys()
dict_keys(['Anti_Aging_Prediction', 'HallMarks_Status', 'HallMarks_Probabilities'])
```
**Note:** HallMarks_Probabilities will be empty by default, unless selected otherwise by additional arguments

To get any specific result dataframe
```
>>> output['HallMarks_Status']
```

#### Additional arguments:
**AgeXtend** also supports the following functionalities along with the Anti-Aging Predictions, that can be chosen using boolean option (True)

| Parameter Name | Description | Type | Default value | **Output(If True)** |
| -------- | -------- | -------- | -------- | -------- |
| probs | Probabilities of having Anti-Aging / Toxic properties | boolean  | False | HallMarks_Probabilities / HallMarks_Toxicity_Probabilities |
| HC | Health/Toxicity Check | boolean  | False | HallMarks_Toxicity_Status |
| TS | HallMark of Aging Tanimoto Similarity Test | boolean  | False | HallMark_like_response |
| BDL | BindingDB Target Information and Druggability (Lipinski Rule) | boolean  | False | Druggability_and_Potential_Targets |


**Example**
```
>>> output = Predictor.predict(agex_obj, TS=True)
>>> output.keys()
dict_keys(['Anti_Aging_Prediction', 'HallMarks_Status', 'HallMarks_Probabilities', 'HallMark_like_response'])
>>> 
>>> output = Predictor.predict(agex_obj, HC=True, BDL=True)
>>> output.keys()
dict_keys(['Anti_Aging_Prediction', 'HallMarks_Toxicity_Status', 'HallMarks_Toxicity_Probabilities', 'Druggability_and_Potential_Targets'])
```


### Browser

To acess the AgeXtend pre-complied prediction data of various Databases:<br/>
```
>>> from AgeXtend import Browser
```
Use openbebl format SMILE of the query compound
```
>>> Browser.search(query='OC(=O)CCCc1c[nH]c2c1cccc2', output='/path/to/output/folder/')
```

Unzip the **AgeXtendBroswerOut.zip** file to visualise/print the generated report (HTML format)
**Note:** the **report file** (AgeXtend_BrowserOut.html) must be in the same folder with the **images** folder


#### Additional arguments:
**AgeXtend** also supports the use of locally complied Predictor module outputs (Folder containing CSV format outputs)

| Parameter Name | Description | Default Database |
| -------- | -------- | -------- |
| path | Path to the AgeXtend pre-complied prediction database | gutMGene |

**Example**
```
>>> Browser.search(path='/path/to/Database/Folder/', query='OC(=O)CCCc1c[nH]c2c1cccc2', output='/path/to/output/folder/')
```
