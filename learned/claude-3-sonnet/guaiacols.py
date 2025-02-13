"""
Classifies: CHEBI:134251 guaiacols
"""
"""
Classifies: CHEBI:27090 guaiacols
Any phenol carrying an additional methoxy substituent at the ortho-position.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    A guaiacol is any phenol carrying an additional methoxy substituent at the ortho-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a guaiacol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for phenol rings (aromatic ring with -OH substituent)
    phenol_pattern = Chem.MolFromSmarts("c1ccccc1O")
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    if not phenol_matches:
        return False, "No phenol moiety found"
    
    # Look for ortho-methoxy substituent on phenol rings
    ortho_methoxy_pattern = Chem.MolFromSmarts("Oc1ccc(OC)cc1")
    found_guaiacol = False
    for phenol_match in phenol_matches:
        atoms = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in phenol_match]
        if any(mol.GetAtomWithIdx(idx).HasSubstructMatch(ortho_methoxy_pattern) for idx in phenol_match):
            found_guaiacol = True
            break
    
    if found_guaiacol:
        return True, "Contains a phenol ring with a methoxy substituent at the ortho position"
    else:
        return False, "No phenol ring with ortho-methoxy substituent found"