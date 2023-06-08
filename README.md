# AgeXtend

### Introduction
 <br>
<div align="center">
<img src="Images/Overview.png"></div>
<br>

**AgeXtend** is a multimodal, bioactivity-based, and mechanism-backed Explainable AI Predictor for geroprotectors. The first work package, **Predictor**, involves a bioactivity-based classifier for predicting geroprotective compounds, coupled with an **explainability** module providing mechanistic insights into biological processes behind the predictions, a **toxicity** module to evaluate their toxicity in the biological context, and finally, a **target** module that suggests protein targets of the putative geroprotective compounds.

## Environment Setup (done using requirement.txt)

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

The following file - [requirement.txt](https://github.com/the-ahuja-lab/AgeXtend/blob/main/env/requirement.txt) will be used for the environment setup and by running the following command.
```
$ pip install -r requirement.txt
```


## How to use AgeXtend?

### Installation using pip 
```
$ pip install -i https://test.pypi.org/simple/ AgeXtend
```

## Work Packages
AgeXtend supports 2 distinct work packages:<br/>
1. Predictor
2. Browser

### Predictor

#### Prediction Module ####

Predicts the anti aging potential for your input SMILES:<br/>
```
>>> from AgeXtend import Predictor
```
Prepare a list of canonical SMILES (Openbabel generated) strings
```
>>> smiles = ['ClCC=C', 'C=CCOC(=O)CC(C)C'] 
```
Create an AgeXtend type object to featurize your query compounds (i.e. list of canonical SMILES)
```
>>> agex_obj = Predictor.featurize(smiles)
```
Use the AgeXtend object for further predictions
```
>>> output = Predictor.predict(agex_obj)
```
Get the list of resulting dataframes that are part of the rest of the modules
```
>>> output.keys()
dict_keys(['Anti_Aging_Prediction', 'Explainability_Status', 'Explainability_Probabilities'])
```
**Note:** Explainability_Probabilities will be empty by default, unless selected otherwise by additional arguments.

To get result dataframe of any specific module
```
>>> output['Explainability_Status']
```

#### Additional arguments:
**AgeXtend** also supports the following functionalities along with the Anti-Aging Predictions, that can be chosen using boolean option (True)

| Parameter Name | Description | Type | Default value | **Output(If True)** |
| -------- | -------- | -------- | -------- | -------- |
| probs | Probabilities of having Anti-Aging Potential / Toxic properties | boolean  | False | Explainability_Probabilities / Explainability_Toxicity_Probabilities |
| HC | Health/Toxicity Check | boolean  | False | Explainability_Toxicity_Status |
| TS | Sub Explainability Level Tanimoto Similarity Test | boolean  | False | Explainability_response |
| BDL | BindingDB Target Information and Druggability (Lipinski Rule) | boolean  | False | Druggability_and_Potential_Targets |


**Example**
```
>>> output = Predictor.predict(agex_obj, TS=True)
>>> output.keys()
dict_keys(['Anti_Aging_Prediction', 'Explainability_Status', 'Explainability_Probabilities', 'Explainability_response'])
>>> 
>>> output = Predictor.predict(agex_obj, HC=True, BDL=True)
>>> output.keys()
dict_keys(['Anti_Aging_Prediction', 'Explainability_Toxicity_Status', 'Explainability_Toxicity_Probabilities', 'Druggability_and_Potential_Targets'])
```


### Browser

To access the AgeXtend pre-complied predictions data of various databases:<br/>
```
>>> from AgeXtend import Browser
```
Use openbabel format SMILE of the query compound
```
>>> Browser.search(query='OC(=O)CCCc1c[nH]c2c1cccc2', output='/path/to/output/folder/')
```

Unzip the **AgeXtendBrowserOut.zip** file to visualize/print the generated report (HTML format)

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

