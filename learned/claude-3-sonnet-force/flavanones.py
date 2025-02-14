"""
Classifies: CHEBI:28863 flavanones
"""
"""
Classifies: CHEBI:26176 flavanone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    A flavanone is a flavan with a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for flavanone scaffold pattern
    flavanone_pattern = Chem.MolFromSmarts("[C@@]1(C2=CC=CC=C2)C(=O)CC(=O)C=2C=C(C=CC2=C1)O")
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No flavanone scaffold found"
    
    # Check for aromatic rings
    aromatic_rings = mol.GetAromaticRings()
    if len(aromatic_rings) < 2:
        return False, "Missing aromatic rings for flavanone"
    
    # Check for carbonyl and ether groups
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if len(carbonyl_matches) != 2:
        return False, f"Found {len(carbonyl_matches)} carbonyl groups, need exactly 2"
    
    ether_pattern = Chem.MolFromSmarts("[OX2]C")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if len(ether_matches) != 1:
        return False, f"Found {len(ether_matches)} ether groups, need exactly 1"

    # Count rotatable bonds - flavanones typically have 1-2
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 1 or n_rotatable > 3:
        return False, "Unexpected number of rotatable bonds for flavanone"

    # Check molecular weight - flavanones typically 200-400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 400:
        return False, "Molecular weight outside typical range for flavanone"

    return True, "Contains flavanone scaffold with carbonyl, ether, and aromatic rings"