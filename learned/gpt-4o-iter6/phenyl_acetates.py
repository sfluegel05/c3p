"""
Classifies: CHEBI:140310 phenyl acetates
"""
from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate is defined as an acetate ester obtained by formal condensation of the carboxy group
    of acetic acid with the hydroxy group of any phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a phenyl acetate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define acetate ester pattern: C(=O)O
    acetate_pattern = Chem.MolFromSmarts("C(=O)O")
    acetate_matches = mol.GetSubstructMatches(acetate_pattern)
    if not acetate_matches:
        return False, "No acetate ester group found"

    # Define aromatic ring pattern (phenyl or equivalent): aromatic carbon
    aromatic_pattern = Chem.MolFromSmarts("a")
    aromatic_matches = mol.GetSubstructMatches(aromatic_pattern)
    if not aromatic_matches:
        return False, "No aromatic phenyl structures found"

    # Verify the presence of acetate group on phenolic/aryl position
    phenyl_acetate_pattern = Chem.MolFromSmarts("cC(=O)O")
    phenyl_acetate_matches = mol.GetSubstructMatches(phenyl_acetate_pattern)
    if not phenyl_acetate_matches:
        return False, "Acetate group not bonded to phenyl moiety"
    
    return True, "Molecule is a phenyl acetate with acetate ester linked to phenyl group"