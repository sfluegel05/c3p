"""
Classifies: CHEBI:16180 N-acylglycine
"""
"""
Classifies: CHEBI:48215 N-acylglycine
An N-acyl-amino acid in which amino acid specified is glycine.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine has an acyl group (R-C(=O)-) attached to the nitrogen of glycine (NH-CH2-COOH).

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
    
    # Look for glycine substructure (NH-CH2-COOH)
    glycine_pattern = Chem.MolFromSmarts("N(C(C(=O)O))")
    glycine_atoms = mol.GetSubstructMatches(glycine_pattern)
    if not glycine_atoms:
        return False, "No glycine substructure found"
    
    # Look for acyl group (R-C(=O)-) attached to nitrogen
    acyl_pattern = Chem.MolFromSmarts("C(=O)N")
    acyl_atoms = mol.GetSubstructMatches(acyl_pattern)
    if not acyl_atoms:
        return False, "No acyl group attached to nitrogen"
    
    # Check if the nitrogen is shared between glycine and acyl group
    common_n_atoms = list(set([atom.GetBeginAtomIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]) & set([atom for match in acyl_atoms for atom in match]) & set([atom for match in glycine_atoms for atom in match]))
    if not common_n_atoms:
        return False, "Acyl group and glycine not connected via nitrogen"
    
    return True, "Contains an acyl group (R-C(=O)-) attached to the nitrogen of glycine (NH-CH2-COOH)"