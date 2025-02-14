"""
Classifies: CHEBI:76579 triradylglycerol
"""
"""
Classifies: CHEBI:36373 triradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
import numpy as np
from sklearn.ensemble import RandomForestClassifier

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES string.
    A triradylglycerol has one of three possible substituent groups - 
    acyl, alkyl, or alk-1-enyl - at each of the three positions (sn-1, sn-2, sn-3) 
    on the glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for three substituents
    oxy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    if len(oxy_atoms) != 6:
        return False, "Must have exactly 6 oxygen atoms (3 ester groups)"

    # Extract fingerprints
    fingerprint = MoleculeDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=2048)
    fingerprint = np.array(fingerprint)

    # Use pre-trained random forest model
    model = RandomForestClassifier()
    model.fit(X_train, y_train)  # Assuming X_train and y_train are provided
    prediction = model.predict([fingerprint])

    if prediction[0]:
        return True, "Molecule matches the structural patterns of a triradylglycerol"
    else:
        return False, "Molecule does not match the structural patterns of a triradylglycerol"

# Example usage: (Assuming X_train and y_train are available)
smiles = "O(C(=O)CCCCCCCCCCCCCCCCCCCCC)[C@@H](COC(=O)CCCCCCCCCCCCC(C)C)COC(=O)CCCCCCC"
result, reason = is_triradylglycerol(smiles)
print(result, reason)