"""
Classifies: CHEBI:25903 peptide antibiotic
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    Peptide antibiotics are characterized by peptide bonds and sometimes macrocyclic structures.

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

    # Check for macrocyclic rings or low molecular weight peptides
    ring_info = mol.GetRingInfo()
    if not ring_info.IsAtomInRingOfSize(12): # Check for a macrocycle, could adjust size
        return False, "No macrocyclic structure detected"

    # Count chiral centers - peptidic antibiotics often have several chiral sites
    chiral_centers = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
    if chiral_centers < 3:  # Arbitrary threshold for complexity
        return False, "Too few chiral centers"
    
    # Further functionality can be added to detect specific macrocyclic motifs, sugars, lipids, etc.

    return True, "Contains peptide bonds and structural features typical of peptide antibiotics"