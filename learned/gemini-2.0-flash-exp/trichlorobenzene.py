"""
Classifies: CHEBI:27096 trichlorobenzene
"""
from rdkit import Chem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene based on its SMILES string.
    A trichlorobenzene has at least one benzene ring with exactly three chlorine atoms attached to it
    either directly or via a single bond to another atom (i.e. -O-Cl).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trichlorobenzene, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count total number of chlorine atoms
    chlorine_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17]
    if len(chlorine_atoms) != 3:
        return False, f"Found {len(chlorine_atoms)} chlorine atoms, need exactly 3"

    # Check for benzene rings
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene_pattern):
       return False, "No benzene ring found"

    #Check for chlorine attachments (direct or via single bond) to a benzene
    chlorine_attachment_pattern = Chem.MolFromSmarts("c[X;!H]-[Cl]") # any non-hydrogen atom
    matches = mol.GetSubstructMatches(chlorine_attachment_pattern)
    
    if len(matches) != 3:
      return False, "There should be exactly 3 chlorine atoms bonded to a benzene ring directly or indirectly."
      
    return True, "Contains a benzene ring with exactly 3 chlorine atoms attached directly or via single bond."