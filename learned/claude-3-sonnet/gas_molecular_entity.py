"""
Classifies: CHEBI:138675 gas molecular entity
"""
"""
Classifies: CHEBI:36963 gas molecular entity
"""
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
import pandas as pd

# Load training data
data = pd.read_csv("gas_molecular_entities.csv")
X = data["smiles"].values
y = data["is_gas"].values

# Split data into train and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Train Random Forest classifier
clf = RandomForestClassifier(n_estimators=100, random_state=42)
clf.fit(X_train, y_train)

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is a gas molecular entity at standard temperature and pressure (STP; 0°C and 100 kPa).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gas molecular entity at STP, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate molecular descriptors
    mol_wt = Descriptors.MolWt(mol)
    n_rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    n_hydrogen_bonds = Descriptors.NumHBonds(mol)
    n_rings = Descriptors.RingCount(mol)
    topological_polar_surface_area = Descriptors.TPSA(mol)
    vapor_pressure = Descriptors.MolUnprioritized(mol)[0]
    
    # Prepare input for prediction
    X = [mol_wt, n_rotatable_bonds, n_hydrogen_bonds, n_rings, topological_polar_surface_area, vapor_pressure]
    
    # Make prediction
    prediction = clf.predict([X])[0]
    
    if prediction:
        return True, "Classified as gas molecular entity at STP based on machine learning model"
    else:
        return False, "Not classified as gas molecular entity at STP based on machine learning model"