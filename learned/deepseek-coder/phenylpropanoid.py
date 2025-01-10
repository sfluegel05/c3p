"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: CHEBI:26195 phenylpropanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    A phenylpropanoid is an organic aromatic compound with a structure based on a phenylpropane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one aromatic ring
    aromatic_rings = [ring for ring in Chem.GetSymmSSSR(mol) if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)]
    if not aromatic_rings:
        return False, "No aromatic ring found"

    # More flexible phenylpropane pattern matching
    # Allows for variations in the propyl chain (e.g., double bonds, substitutions)
    phenylpropane_patterns = [
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4][CX4][CX4]"),  # Standard phenylpropane
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4]=[CX4]"),      # Double bond in chain
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4][CX4]"),       # Shorter chain
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4][CX4][CX4][CX4]")  # Longer chain
    ]

    has_phenylpropane = any(mol.HasSubstructMatch(pattern) for pattern in phenylpropane_patterns)
    if not has_phenylpropane:
        return False, "No phenylpropane-like skeleton found"

    # Check for common phenylpropanoid features
    # These are not required but increase confidence
    common_features = [
        Chem.MolFromSmarts("[OX2H]"),  # Hydroxyl
        Chem.MolFromSmarts("[OX2][CX4]"),  # Methoxy
        Chem.MolFromSmarts("[CX3]=[OX1]"),  # Carbonyl
        Chem.MolFromSmarts("[CX3]=[CX3]")  # Double bond
    ]

    feature_count = sum(mol.HasSubstructMatch(feat) for feat in common_features)
    if feature_count == 0:
        return False, "No common phenylpropanoid features found"

    # Check molecular weight - phenylpropanoids typically >100 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, "Molecular weight too low for phenylpropanoid"

    # If all checks pass, classify as phenylpropanoid
    return True, "Contains phenylpropane-like skeleton and common features"