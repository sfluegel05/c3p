"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is classified as a glycosaminoglycan based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool, str: True if the molecule is a glycosaminoglycan, False otherwise, with a reason
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count relevant atoms
    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    sulfur_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)

    # Calculate descriptors
    mol_wt = Descriptors.ExactMolWt(mol)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)

    # Glycosaminoglycan characteristics
    # Assume that they have repetitive polysaccharide units, often with sulfate groups.
    if nitrogen_count < 4:
        return False, "Insufficient nitrogen atoms for typical amino sugar structure."
    if oxygen_count < 8:
        return False, "Insufficient oxygen atoms for typical polysaccharide structure."
    if sulfur_count == 0:
        return False, "No sulfur atoms detected; some GAGs are heavily sulfated."
    if n_rotatable < 15:
        return False, "Insufficient flexibility suggesting it's not a polymer structure."
    if mol_wt < 500:  # GAGs are typically large
        return False, "Molecular weight too low for a glycosaminoglycan."

    reason = "Molecule meets glycosaminoglycan characteristics: "
    reason += f"Nitrogen: {nitrogen_count}, Oxygen: {oxygen_count}, "
    reason += f"Sulfur: {sulfur_count}, Molecular weight: {mol_wt:.2f}, "
    reason += f"Rotatable bonds: {n_rotatable}"

    return True, reason