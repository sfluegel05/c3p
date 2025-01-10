"""
Classifies: CHEBI:29348 fatty amide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is defined by an amide group derived from a long chain fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for amide bond pattern (C(=O)N)
    amide_pattern = Chem.MolFromSmarts("C(=O)[NX3][*]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found"

    # Ensure long carbon chain (at least 8 carbons) is attached to carbonyl
    c_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            if len(Chem.FindAtomEnvironmentOfRadiusN(mol, 8, atom.GetIdx())) >= 8:
                c_count += 1

        # Check for a sufficiently long carbon chain attached
        if c_count >= 8:
            break
    else:
        return False, "Insufficient carbon chain length"

    return True, "Contains amide bond derived from a long chain fatty acid"