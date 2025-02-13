"""
Classifies: CHEBI:33709 amino acid
"""
"""
Classifies: CHEBI:33597 amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    An amino acid is a carboxylic acid containing one or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid group (-C(=O)O)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxyl_match = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_match:
        return False, "No carboxylic acid group found"
    
    # Look for amino group (-N)
    amino_pattern = Chem.MolFromSmarts("N")
    amino_match = mol.GetSubstructMatches(amino_pattern)
    if not amino_match:
        return False, "No amino group found"
    
    # Check if amino group is attached to carboxylic acid carbon
    for amino_idx in amino_match:
        for carboxyl_idx in carboxyl_match:
            carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx[0])
            amino_atom = mol.GetAtomWithIdx(amino_idx)
            if any(bond.GetBeginAtomIdx() == amino_idx and bond.GetEndAtomIdx() == carboxyl_idx[0] or
                   bond.GetBeginAtomIdx() == carboxyl_idx[0] and bond.GetEndAtomIdx() == amino_idx
                   for bond in carboxyl_atom.GetBonds()):
                break
        else:
            continue
        break
    else:
        return False, "Amino group not attached to carboxylic acid carbon"
    
    # Check for additional functional groups if present
    # e.g., hydroxyl (-OH), thiol (-SH), phosphate (-OP(O)(O)=O), etc.
    
    return True, "Contains a carboxylic acid group and an amino group attached to the carboxylic acid carbon"