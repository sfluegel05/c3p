"""
Classifies: CHEBI:25903 peptide antibiotic
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    Peptide antibiotics are characterized by peptide bonds, often have macrocyclic structures,
    possess certain antimicrobial functional groups, and typically have complex stereochemistry.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a peptide antibiotic, False otherwise
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
    
    # Check for macrocyclic rings or multiple smaller rings
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    large_ring_detected = any(len(ring) >= 12 for ring in atom_rings)
    multiple_ring_system = len(atom_rings) > 1  # Allow for multiple smaller rings

    if not (large_ring_detected or multiple_ring_system):
        return False, "No significant ring systems detected"

    # Check for chiral centers
    chiral_centers = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
    if chiral_centers < 3:
        return False, "Insufficient chiral centers detected, which are common in peptide antibiotics"

    # Check for specific antimicrobial motifs (perhaps better defined with domain knowledge)
    # For now, look for guanidyl groups, often found in arginine-rich peptides
    guanidine_pattern = Chem.MolFromSmarts("NC(=N)N")
    if not mol.HasSubstructMatch(guanidine_pattern):
        return True, "Contains peptide bonds, significant ring structures, and chiral centers but no guanidyl groups detected"

    return True, "Contains peptide bonds, significant ring structures, and sufficient chiral centers"

# Example usage:
# result, reason = is_peptide_antibiotic("CC(C)C[C@@H](N)C(=O)N[C@@H]1C(C)OC(=O)C(C)N(C)C1=O")
# print(result, reason)