"""
Classifies: CHEBI:16385 organic sulfide
"""
from rdkit import Chem
from rdkit.Chem import AllChem


def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide based on its SMILES string.
    An organic sulfide has the structure R-S-R', where R and R' are carbon-containing groups (not hydrogen).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organic sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of sulfur atom
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'S']
    if not sulfur_atoms:
        return False, "No sulfur atom found"
    
    # Check for sulfur bonded to two carbons
    sulfide_pattern = Chem.MolFromSmarts("[CX4,CX3][S][CX4,CX3]") # Modified SMARTS to allow for sp2/sp3 carbons
    if not mol.HasSubstructMatch(sulfide_pattern):
         return False, "Sulfur is not bonded to two carbon atoms."
     
    return True, "Has the structure R-S-R' (R and R' not hydrogen)"