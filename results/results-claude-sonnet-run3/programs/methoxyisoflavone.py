from rdkit import Chem
from rdkit.Chem import AllChem

def is_methoxyisoflavone(smiles: str):
    """
    Determines if a molecule is a methoxyisoflavone (isoflavone with at least one methoxy substituent).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methoxyisoflavone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for isoflavone core structure
    # SMARTS pattern for isoflavone core:
    # A more specific pattern that matches the 3-phenylchromen-4-one structure
    isoflavone_pattern = Chem.MolFromSmarts('[#6]1=[#6]-[#6](=[O])-c2c(occ1-c1ccccc1)cccc2')
    
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "Not an isoflavone core structure"

    # Check for methoxy group
    # Using a more specific SMARTS pattern for methoxy groups
    methoxy_pattern = Chem.MolFromSmarts('cOC')
    if not mol.HasSubstructMatch(methoxy_pattern):
        return False, "No methoxy substituent found"

    # Count methoxy groups
    methoxy_matches = mol.GetSubstructMatches(methoxy_pattern)
    num_methoxy = len(methoxy_matches)

    # Verify that the methoxy groups are actually connected to the ring system
    for match in methoxy_matches:
        o_atom = mol.GetAtomWithIdx(match[1])
        if o_atom.GetDegree() == 2:  # Oxygen should have exactly 2 bonds
            c_atom = mol.GetAtomWithIdx(match[2])
            if c_atom.GetDegree() == 1 and c_atom.GetNumImplicitHs() == 3:  # Should be CH3
                return True, f"Methoxyisoflavone with {num_methoxy} methoxy group(s)"

    return False, "No valid methoxy substituent found"
# Pr=None
# Recall=0.0