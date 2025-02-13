"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: CHEBI:35484 beta-lactam
A lactam in which the amide bond is contained within a four-membered ring,
which includes the amide nitrogen and the carbonyl carbon.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score

# Load dataset of known beta-lactams and non-beta-lactams
# (not included here, but assume it's a list of tuples (smiles, is_beta_lactam))
dataset = [...] 

# Extract features and labels
smiles, labels = zip(*dataset)
mols = [Chem.MolFromSmiles(s) for s in smiles]
featurizer = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in MoleculeDescriptors._describeList])
features = [featurizer.CalcDescriptors(mol) for mol in mols]

# Split into train and test sets
X_train, X_test, y_train, y_test = train_test_split(features, labels, test_size=0.2, random_state=42)

# Train a random forest classifier
clf = RandomForestClassifier(n_estimators=100, random_state=42)
clf.fit(X_train, y_train)

# Evaluate on test set
y_pred = clf.predict(X_test)
f1 = f1_score(y_test, y_pred)
print(f"F1 score on test set: {f1:.3f}")

def is_beta_lactam(smiles):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    features = featurizer.CalcDescriptors(mol)
    prediction = clf.predict([features])[0]
    
    if prediction:
        return True, "Classified as beta-lactam based on structural features"
    else:
        return False, "Not classified as beta-lactam based on structural features"