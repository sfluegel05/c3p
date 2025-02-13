"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: CHEBI:37496 nucleobase analogue
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    A nucleobase analogue is a molecule that can substitute for a normal nucleobase in nucleic acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for nucleobase ring systems
    base_patterns = ['c1ncnc[nH]1', 'c1ncnc2[nH]cnc12', 'c1c[nH]c(=O)[nH]c1', 'c1c[nH]c(=O)nc1']
    base_pattern_mols = [Chem.MolFromSmarts(pat) for pat in base_patterns]
    
    has_base = any(mol.HasSubstructMatch(pat_mol) for pat_mol in base_pattern_mols)
    if not has_base:
        return False, "No nucleobase ring system found"
    
    # Look for modifications (substitutions, additions, etc.) on the ring
    modified_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [6, 7, 8, 1]]
    if not modified_atoms:
        return False, "No modifications found on nucleobase ring"
    
    # Check for common nucleobase analogue modifications
    mod_patterns = ['[OH]', '[NH2]', '[OX2H]', '[OX2H][CX4H2]', '[NX3][CX3](=[OX1])[CX3]', '[NX3][CX3](=[OX1])N', '[SX2]']
    mod_pattern_mols = [Chem.MolFromSmarts(pat) for pat in mod_patterns]
    
    has_mod = any(mol.HasSubstructMatch(pat_mol) for pat_mol in mod_pattern_mols)
    if not has_mod:
        return False, "No common nucleobase analogue modifications found"
    
    return True, "Contains a modified nucleobase ring system"