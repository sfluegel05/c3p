"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies: CHEBI:51574 organoiodine compound
An organoiodine compound is a compound containing at least one carbon-iodine bond.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organoiodine_compound(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organoiodine compound, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of iodine atoms
    iodine_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 53]
    if not iodine_atoms:
        return False, "No iodine atoms present"

    # Check for carbon-iodine bonds
    has_c_i_bond = False
    for atom in iodine_atoms:
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                has_c_i_bond = True
                break
        if has_c_i_bond:
            break

    if has_c_i_bond:
        return True, "Contains at least one carbon-iodine bond"
    else:
        return False, "No carbon-iodine bonds found"

    # Additional checks (optional)
    # mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    # if mol_wt < 100:
    #     return False, "Molecular weight too low for an organoiodine compound"

    # return True, "Contains at least one carbon-iodine bond"