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

    # Check for macrocyclic rings or larger structures
    # Assess ring sizes commonly seen in macrocyclic peptide antibiotics
    ring_info = mol.GetRingInfo()
    large_ring_detected = any(ring_info.NumAtomRingsSize(r_size) > 0 for r_size in range(12, 30))  # Looking for ring sizes between 12-30
    if not large_ring_detected:
        return False, "No macrocyclic structure or typical peptide ring sizes detected"

    # Peptides, particularly cyclic ones, may still have a considerable number of chiral centers
    # Refine the threshold based on additional data/knowledge or still check for their presence
    chiral_centers = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
    # Optional: use a threshold if deemed necessary after further verification against known data

    return True, "Contains peptide bonds and structural features typical of peptide antibiotics"

# Example usage:
# result, reason = is_peptide_antibiotic("CC(C)C[C@@H](N)C(=O)N[C@@H]1C(C)OC(=O)C(C)N(C)C1=O")
# print(result, reason)