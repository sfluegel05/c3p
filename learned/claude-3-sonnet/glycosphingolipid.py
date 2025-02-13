"""
Classifies: CHEBI:24402 glycosphingolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split

# Load a dataset of known glycosphingolipids and non-glycosphingolipids
dataset = pd.read_csv('glycosphingolipid_dataset.csv')
X = dataset['SMILES'].values
y = dataset['is_glycosphingolipid'].values

# Calculate molecular descriptors
nms = [x[0] for x in Descriptors._descList]
calc = MoleculeDescriptors.MolecularDescriptorCalculator(nms)
X_desc = [calc.CalcDescriptors(Chem.MolFromSmiles(smile)) for smile in X]

# Split the dataset into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X_desc, y, test_size=0.2, random_state=42)

# Train a Random Forest classifier
clf = RandomForestClassifier(n_estimators=100, random_state=42)
clf.fit(X_train, y_train)

def is_glycosphingolipid(smiles):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string and a trained Random Forest classifier.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosphingolipid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    desc = calc.CalcDescriptors(mol)
    prediction = clf.predict([desc])

    if prediction[0]:
        return True, "Classified as a glycosphingolipid by the Random Forest model"
    else:
        return False, "Classified as not a glycosphingolipid by the Random Forest model"