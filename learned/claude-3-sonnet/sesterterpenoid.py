"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem.Fingerprints import FingerprintMols
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score

def extract_features(mol):
    """
    Extract relevant features from a molecule for sesterterpenoid classification.
    """
    features = []
    
    # Check for C25 skeleton and isoprene units
    smi = Chem.MolToSmiles(mol)
    c25_pattern = Chem.MolFromSmarts('[C&D3]([C&D3])[C&D3]([C&D3])[C&D3]([C&D3])[C&D3]([C&D3])[C&D3]([C&D3])[C&D3]')
    isoprene_pattern = Chem.MolFromSmarts('C=C(C)C')
    features.append(mol.HasSubstructMatch(c25_pattern))
    features.append(mol.HasSubstructMatch(isoprene_pattern))
    
    # Check for specific ring systems and functional groups
    sesterterpene_rings = Chem.MolFromSmarts('C1CCC2CCCCC2C1')
    ester_group = Chem.MolFromSmarts('C(=O)OC')
    epoxide_group = Chem.MolFromSmarts('C1OC1')
    features.append(mol.HasSubstructMatch(sesterterpene_rings))
    features.append(mol.HasSubstructMatch(ester_group))
    features.append(mol.HasSubstructMatch(epoxide_group))
    
    # Include molecular descriptors
    features.append(Descriptors.MolLogP(mol))
    features.append(Descriptors.TPSA(mol))
    features.append(rdMolDescriptors.CalcNumRotatableBonds(mol))
    
    # Include Morgan fingerprints
    fps = FingerprintMols.FingerprintMol(mol, maxPaths=5, fpSize=2048)
    features.extend(fps.GetNonzeroElements().values())
    
    return features

def is_sesterterpenoid(smiles):
    """
    Determine if a molecule is a sesterterpenoid based on its SMILES string.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    features = extract_features(mol)
    prediction = model.predict([features])
    
    if prediction[0]:
        return True, "Classified as a sesterterpenoid"
    else:
        return False, "Not classified as a sesterterpenoid"

# Load or generate training data
# (implementation omitted for brevity)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Train a Random Forest classifier
model = RandomForestClassifier(n_estimators=100, random_state=42)
model.fit(X_train, y_train)

# Evaluate the model
y_pred = model.predict(X_test)
f1 = f1_score(y_test, y_pred)
print(f"F1 score: {f1:.2f}")