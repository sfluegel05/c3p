"""
Classifies: CHEBI:25903 peptide antibiotic
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    Peptide antibiotics are characterized by peptide bonds and often contain macrocyclic structures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a peptide antibiotic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for peptide bonds (look for -C(=O)-N- pattern)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No peptide bonds detected"

    # Check for macrocyclic rings
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    macrocyclic_detected = any(len(ring) >= 12 for ring in atom_rings)
    if not macrocyclic_detected:
        return False, "No macrocyclic structures detected"

    # Check for chiral centers
    chiral_centers = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
    if chiral_centers == 0:
        return False, "No chiral centers detected, which are common in peptide antibiotics"

    return True, "Contains peptide bonds, macrocyclic rings, and chiral centers"

# Example usage:
# result, reason = is_peptide_antibiotic("CC(C)C[C@@H](N)C(=O)N[C@@H]1C(C)OC(=O)C(C)N(C)C1=O")
# print(result, reason)