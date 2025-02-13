"""
Classifies: CHEBI:37143 organofluorine compound
"""
"""
Classifies: CHEBI:33839 organofluorine compound

An organofluorine compound is a compound containing at least one carbon-fluorine bond.
"""
from rdkit import Chem

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organofluorine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carbon-fluorine bonds
    has_c_f_bond = False
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 9) or (atom1.GetAtomicNum() == 9 and atom2.GetAtomicNum() == 6):
            has_c_f_bond = True
            break
    
    # Check for fluorinated alkyl groups
    if not has_c_f_bond:
        smarts_patterns = ['[C]C(F)(F)(F)', '[C]C(F)(F)']
        for pattern in smarts_patterns:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                has_c_f_bond = True
                break
    
    if has_c_f_bond:
        return True, "Contains at least one carbon-fluorine bond"
    else:
        return False, "Does not contain any carbon-fluorine bonds"