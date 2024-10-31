from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_methoxynaphthalene(smiles: str):
    """
    Determines if a molecule is a methoxynaphthalene (naphthalene with one or more methoxy substituents).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methoxynaphthalene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for naphthalene core
    naphthalene_pattern = Chem.MolFromSmarts('c1ccc2ccccc2c1')
    if not mol.HasSubstructMatch(naphthalene_pattern):
        return False, "No naphthalene core found"

    # Check for methoxy groups
    methoxy_pattern = Chem.MolFromSmarts('COc')
    methoxy_matches = mol.GetSubstructMatches(methoxy_pattern)
    
    if not methoxy_matches:
        return False, "No methoxy groups found"

    # Verify methoxy groups are attached to naphthalene core
    naphthalene_atoms = set()
    for match in mol.GetSubstructMatches(naphthalene_pattern):
        naphthalene_atoms.update(match)

    methoxy_on_naphthalene = False
    methoxy_positions = []
    
    for match in methoxy_matches:
        # match[2] is the aromatic carbon the methoxy is attached to
        if match[2] in naphthalene_atoms:
            methoxy_on_naphthalene = True
            # Get position number relative to naphthalene core
            for naphthalene_match in mol.GetSubstructMatches(naphthalene_pattern):
                if match[2] in naphthalene_match:
                    pos = naphthalene_match.index(match[2]) + 1
                    methoxy_positions.append(str(pos))
                    break

    if not methoxy_on_naphthalene:
        return False, "Methoxy groups not attached to naphthalene core"

    return True, f"Methoxynaphthalene with methoxy groups at position(s): {', '.join(methoxy_positions)}"
# Pr=0.9333333333333333
# Recall=1.0