"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
"""
Classifies: CHEBI:53054 quaternary ammonium ion
A derivative of ammonium, NH4(+), in which all four of the hydrogens bonded to nitrogen
have been replaced with univalent (usually organyl) groups.
"""
from rdkit import Chem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a quaternary ammonium ion based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quaternary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find nitrogen atoms with 4 substituents
    quat_n_candidates = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetDegree() == 4]
    
    for candidate in quat_n_candidates:
        # Check if nitrogen has a formal charge of +1
        if candidate.GetFormalCharge() == 1:
            # Check if the sum of formal charges of substituents is -1
            neighbor_charges = sum(mol.GetAtomWithIdx(neighbor_idx).GetFormalCharge() for neighbor_idx in candidate.GetNeighbors())
            if neighbor_charges == -1:
                return True, "Contains a positively charged quaternary nitrogen with 4 substituents"
    
    return False, "No quaternary ammonium ion found"