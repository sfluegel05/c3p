"""
Classifies: CHEBI:75874 para-terphenyl
"""
"""
Classifies: CHEBI:39158 para-terphenyl
A ring assembly based on a 1,4-diphenylbenzene skeleton and its substituted derivatives thereof.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule is a para-terphenyl based on its SMILES string.
    A para-terphenyl is a ring assembly based on a 1,4-diphenylbenzene skeleton
    and its substituted derivatives thereof.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a para-terphenyl, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for para-terphenyl core structure
    # C1=CC(=CC=C1C2=CC=CC=C2)C3=CC=CC=C3
    terphenyl_pattern = Chem.MolFromSmarts("c1ccc(cc1)c2ccc(cc2)c3ccccc3")
    if not mol.HasSubstructMatch(terphenyl_pattern):
        return False, "No para-terphenyl core structure found"

    # Check for substituents
    allowed_atoms = [6, 8, 17, 35]  # C, O, Cl, Br
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, f"Atom {atom.GetSymbol()} not allowed in para-terphenyls"

    # Check for common functional groups like methoxy, hydroxy, ester, etc.
    functional_groups = ["CO", "OC", "CO.O", "C(=O)O"]
    for fg in functional_groups:
        fg_pattern = Chem.MolFromSmarts(fg)
        if mol.HasSubstructMatch(fg_pattern):
            return True, f"Contains para-terphenyl core and allowed functional group {fg}"

    return True, "Contains para-terphenyl core with allowed substituents"