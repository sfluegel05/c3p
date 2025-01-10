"""
Classifies: CHEBI:15841 polypeptide
"""
"""
Classifies: CHEBI:16670 polypeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide based on its SMILES string.
    A polypeptide is a peptide containing ten or more amino acid residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Broad pattern to match amide bonds (peptide bonds)
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3H0,nX3H1]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    
    # Each amide bond corresponds to one peptide linkage
    # Number of residues = number of amide bonds + 1
    num_residues = len(amide_matches) + 1
    
    if num_residues < 10:
        return False, f"Found {num_residues} amino acid residues, need at least 10"
    
    # Check for typical polypeptide elements
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if n_count < num_residues or o_count < num_residues:
        return False, "Insufficient nitrogen or oxygen atoms for polypeptide"
    
    # Check for a continuous chain of residues
    # Ensure that the amide bonds are connected in a chain
    # This is a simplified check and may not catch all edge cases
    # but should work for most polypeptides
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("[NX3H0,nX3H1][CX4H1]([CX3](=[OX1]))[NX3H0,nX3H1]")):
        return False, "Residues not connected in a continuous chain"
    
    return True, f"Contains {num_residues} amino acid residues, qualifies as a polypeptide"