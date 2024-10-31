from rdkit import Chem
from rdkit.Chem import AllChem

def is_dimethoxyflavanone(smiles: str):
    """
    Determines if a molecule is a dimethoxyflavanone (flavanone with two methoxy substituents).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dimethoxyflavanone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for flavanone core structure
    # Matches 2-phenyl-2,3-dihydro-4H-chromen-4-one core
    flavanone_pattern = Chem.MolFromSmarts('[O;!R][#6;R1]1[#6;R1][#6;R1](=[O;!R])[#6;R1]2[#6;R1][#6;R1][#6;R1][#6;R1][#6;R1]2[#8;R1]1')
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No flavanone core structure found"

    # Count methoxy groups (-OCH3)
    methoxy_pattern = Chem.MolFromSmarts('O-[CH3]')
    methoxy_matches = mol.GetSubstructMatches(methoxy_pattern)
    
    if len(methoxy_matches) < 2:
        return False, f"Only {len(methoxy_matches)} methoxy groups found, need at least 2"

    # Get the atoms that are part of the flavanone core
    core_matches = mol.GetSubstructMatches(flavanone_pattern)
    if not core_matches:
        return False, "Could not map flavanone core"
        
    core_atoms = set(core_matches[0])

    # Find methoxy groups connected to either ring of the flavanone
    methoxy_on_rings = 0
    for match in methoxy_matches:
        o_atom = mol.GetAtomWithIdx(match[0])
        for neighbor in o_atom.GetNeighbors():
            if neighbor.GetIdx() not in match:  # Not the methyl carbon
                if neighbor.IsInRing():
                    methoxy_on_rings += 1
                break

    if methoxy_on_rings >= 2:
        return True, "Valid dimethoxyflavanone with 2 methoxy substituents"
    else:
        return False, "Less than 2 methoxy groups connected to ring system"
# Pr=None
# Recall=0.0