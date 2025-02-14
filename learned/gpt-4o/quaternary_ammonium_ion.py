"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
from rdkit import Chem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a quaternary ammonium ion based on its SMILES string.
    A quaternary ammonium ion is characterized by a nitrogen atom with a formal positive charge
    bonded to four organic groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a quaternary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for quaternary ammonium feature: [N+] with 4 single bonds
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen
            if atom.GetFormalCharge() == 1:  # Positive charge
                if len(atom.GetBonds()) == 4:  # Four single bonds
                    return True, "Contains a quaternary ammonium nitrogen: [N+](C)(C)(C)(C)"

    return False, "Does not contain a quaternary ammonium nitrogen"