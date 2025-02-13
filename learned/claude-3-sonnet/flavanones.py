"""
Classifies: CHEBI:28863 flavanones
"""
"""
Classifies: CHEBI:27555 flavanone
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score

# Load training data (replace with your own data)
flavanones = [Chem.MolFromSmiles(smiles) for smiles in ["CC(C)=CCc1c(O)c(CC=C(C)C)c2O[C@@H](CC(=O)c2c1O)c1ccc(O)cc1", ...]]
non_flavanones = [Chem.MolFromSmiles(smiles) for smiles in ["CCCC", ...]]

# Prepare data
X = flavanones + non_flavanones
y = [1] * len(flavanones) + [0] * len(non_flavanones)

# Split into train and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Calculate fingerprints
fingerprint = MoleculeDescriptors.MorganFingerprint

# Train the classifier
clf = RandomForestClassifier()
clf.fit([fingerprint(mol) for mol in X_train], y_train)

# Evaluate the model
y_pred = clf.predict([fingerprint(mol) for mol in X_test])
f1 = f1_score(y_test, y_pred)
print(f"F1 score: {f1:.2f}")

def is_flavanone(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string using a trained machine learning model.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavanone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    prediction = clf.predict([fingerprint(mol)])[0]
    if prediction == 1:
        return True, "Classified as a flavanone by the trained model"
    else:
        return False, "Classified as non-flavanone by the trained model"