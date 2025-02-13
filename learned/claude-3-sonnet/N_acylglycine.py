"""
Classifies: CHEBI:16180 N-acylglycine
"""
"""
Classifies: CHEBI:35632 N-acylglycine
An N-acyl-amino acid in which amino acid specified is glycine.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylglycine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycine substructure (H2N-CH2-COOH)
    glycine_pattern = Chem.MolFromSmarts("N(C(=O)O)CC")
    if not mol.HasSubstructMatch(glycine_pattern):
        return False, "No glycine substructure found"
    
    # Look for acyl group (R-C(=O)-) attached to the amine nitrogen
    acyl_pattern = Chem.MolFromSmarts("N(C(=O)[C,c])")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No acyl group attached to amine nitrogen"
    
    # Check if the acyl group and glycine are connected
    acyl_atoms = mol.GetSubstructMatches(acyl_pattern)
    glycine_atoms = mol.GetSubstructMatches(glycine_pattern)
    
    # Find common nitrogen atom
    common_n_atoms = list(set([atom.GetBeginAtomIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]) & set([atom for match in acyl_atoms for atom in match]) & set([atom for match in glycine_atoms for atom in match]))
    
    if len(common_n_atoms) != 1:
        return False, "Acyl group and glycine not connected via amine nitrogen"
    
    return True, "Contains a glycine moiety with an acyl group attached to the amine nitrogen"