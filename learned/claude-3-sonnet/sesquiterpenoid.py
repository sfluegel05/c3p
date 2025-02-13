"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
"""
Classifies: CHEBI:35411 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

# Comprehensive list of SMARTS patterns for common sesquiterpene skeletons
sesquiterpene_skeletons = [
    "[C@@]12CC[C@@](C)(O)[C@]1([H])[C@]1([H])[C@]([H])(CCC2=C)C1(C)C", # ent-spathulenol
    "[C@@]12CC[C@@](C)(O)[C@]1([H])[C@]1([H])[C@]([H])(CC[C@]2(C)C)C1(C)C", # curcumol
    "O1[C@]23[C@@](CC[C@@H]2C)(C(C[C@]1(O)[C@@H](C3)C(C)C)=C)[H]", # curcumol
    "O1[C@]23[C@@](CC[C@@H]2C)(C(C[C@]1(O)[C@@H](C3)C(C)C)=C)[H]", # cucumin H
    "[C@]12CC[C@@](C)(O)[C@]1([H])[C@]1([H])[C@]([H])(CC[C@@]2(C)C)C1(C)C", # epi-curcumol
    "O[C@]12C(=C)CC[C@@H]([C@H]1C=C(CO)CC2)C(C)C", # donacinol A
    # Add more SMARTS patterns for other common sesquiterpene skeletons
]

# Additional structural features characteristic of sesquiterpenoids
sesquiterpene_features = [
    "[OH]", # Alcohol
    "[O]", # Ether
    "[OC(=O)]", # Ester
    "[C@@]", # Stereocenter
    "[r6,r5,r4]", # Rings of size 6, 5, or 4
    "[/C=C(/C)]", # Conjugated double bonds
    "[C@@H]1[C@@H]([C@@H]([C@@H]([C@@H]1)C)C)C", # Decalin system
    "[C@@]12[C@@H]([C@@H]([C@@H]1[C@@H]2)C)C", # Bicyclo[3.3.1] system
]

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sesquiterpene skeletons
    for skeleton in sesquiterpene_skeletons:
        skeleton_mol = Chem.MolFromSmarts(skeleton)
        if mol.HasSubstructMatch(skeleton_mol):
            return True, f"Matched sesquiterpene skeleton: {skeleton}"

    # Check for additional sesquiterpenoid features
    feature_count = 0
    for feature in sesquiterpene_features:
        feature_mol = Chem.MolFromSmarts(feature)
        if mol.HasSubstructMatch(feature_mol):
            feature_count += 1

    # Require at least 3 additional features
    if feature_count >= 3:
        return True, "Matched at least 3 additional sesquiterpenoid features"

    # Check for carbon count in a reasonable range
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 13 or c_count > 17:
        return False, "Carbon count outside the expected range for sesquiterpenoids"

    # Check for a reasonable number of rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3 or n_rotatable > 10:
        return False, "Rotatable bond count outside the expected range for sesquiterpenoids"

    # If none of the above conditions are met, classify as non-sesquiterpenoid
    return False, "Does not match sesquiterpenoid criteria"