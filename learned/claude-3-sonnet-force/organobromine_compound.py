"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: CHEBI:38241 organobromine compound
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    An organobromine compound is a compound containing at least one bromine atom within a certain distance from a carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organobromine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule contains bromine atoms
    has_br = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 35:  # Bromine
            has_br = True
            break

    if not has_br:
        return False, "Does not contain any bromine atoms"

    # Check if bromine atoms are within a certain distance from carbon atoms
    max_distance = 2  # Adjust this value as needed
    has_organo_br = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 35:  # Bromine
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Carbon
                    has_organo_br = True
                    break
                else:
                    conformer = mol.GetConformer()
                    dist = AllChem.GetBondLength(conformer, atom.GetIdx(), neighbor.GetIdx())
                    if dist <= max_distance:
                        has_organo_br = True
                        break
            if has_organo_br:
                break

    if has_organo_br:
        return True, "Contains at least one bromine atom within a certain distance from a carbon atom"
    else:
        return False, "Does not contain any bromine atoms within a certain distance from carbon atoms"