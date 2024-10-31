from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monohydroxyanthraquinone(smiles: str):
    """
    Determines if a molecule is a monohydroxyanthraquinone (anthraquinone with exactly one hydroxy group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monohydroxyanthraquinone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check if molecule contains anthraquinone core
    anthraquinone_pattern = Chem.MolFromSmarts('O=C1c2ccccc2C(=O)c2ccccc12')
    if not mol.HasSubstructMatch(anthraquinone_pattern):
        return False, "No anthraquinone core found"

    # Count number of OH groups directly attached to the anthraquinone core
    oh_pattern = Chem.MolFromSmarts('[OH]')
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    
    # Get the atoms in the anthraquinone core
    anthraquinone_matches = mol.GetSubstructMatches(anthraquinone_pattern)
    if not anthraquinone_matches:
        return False, "No anthraquinone core found"
        
    anthraquinone_atoms = set(anthraquinone_matches[0])
    
    # Count OH groups attached to anthraquinone core
    oh_count = 0
    for match in oh_matches:
        oh_oxygen = match[0]
        # Get carbon atom the OH is attached to
        for bond in mol.GetAtomWithIdx(oh_oxygen).GetBonds():
            if bond.GetOtherAtomIdx(oh_oxygen) in anthraquinone_atoms:
                oh_count += 1
                break

    if oh_count == 0:
        return False, "No hydroxy groups found on anthraquinone core"
    elif oh_count > 1:
        return False, f"Found {oh_count} hydroxy groups, but monohydroxyanthraquinone must have exactly one"
    else:
        return True, "Found anthraquinone core with exactly one hydroxy group"
# Pr=1.0
# Recall=1.0