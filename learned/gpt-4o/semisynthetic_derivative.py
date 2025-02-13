"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule could be classified as a semisynthetic derivative based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a semisynthetic derivative, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Check for stereochemistry presence
    if not any(atom.HasChiralTag() for atom in mol.GetAtoms()):
        return None, "No stereochemistry found, cannot infer semisynthetic nature"

    # Assess complexity - simple measure could be based on heavy atom count and rotatable bonds
    heavy_atoms = rdMolDescriptors.CalcNumHeavyAtoms(mol)
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    if heavy_atoms < 20 or rotatable_bonds < 5:
        return False, "Molecule seems too simple to be semisynthetic"

    # Look for the presence of functional groups typically added semisynthetically (e.g., acetyl groups)
    acetyl_group = Chem.MolFromSmarts("CC(=O)O")
    if mol.HasSubstructMatch(acetyl_group):
        return True, "Contains acetyl group, indicates possible chemical modification"
    
    ether_group = Chem.MolFromSmarts("COC")
    if mol.HasSubstructMatch(ether_group):
        return True, "Contains ether group, indicates possible chemical modification"
    
    ester_group = Chem.MolFromSmarts("C(=O)OC")
    if mol.HasSubstructMatch(ester_group):
        return True, "Contains ester group, indicates possible chemical modification"

    # Check for commonly semisynthetic backbone (e.g., lactone structures from macrolides)
    lactone_group = Chem.MolFromSmarts("C1OC(=O)OCC1")
    if mol.HasSubstructMatch(lactone_group):
        return True, "Contains lactone structure, indicates possible semisynthetic derivative"
    
    # If none of the patterns are found or clear evidence is absent, return indeterminate
    return None, "Cannot reliably classify as semisynthetic derivative based on structure alone"