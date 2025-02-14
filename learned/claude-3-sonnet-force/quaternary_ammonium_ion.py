"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
"""
Classifies: CHEBI:53054 quaternary ammonium ion
A derivative of ammonium, NH4(+), in which all four of the hydrogens bonded to nitrogen
have been replaced with univalent (usually organyl) groups.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
        # Check if all substituents are univalent (usually organyl) groups
        if all(mol.GetAtomWithIdx(neighbor_idx).GetDegree() == 4 for neighbor_idx in candidate.GetNeighbors()):
            # Check if the molecule has a net positive charge
            mol_charge = rdMolDescriptors.CalcTPSACharge(mol)
            if mol_charge > 0:
                return True, "Contains a positively charged quaternary nitrogen with 4 univalent substituents"
    
    return False, "No quaternary ammonium ion found"