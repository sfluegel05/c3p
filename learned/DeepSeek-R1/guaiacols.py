"""
Classifies: CHEBI:134251 guaiacols
"""
"""
Classifies: guaiacols (CHEBI:49211)
Definition: Any phenol carrying an additional methoxy substituent at the ortho-position.
"""
from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    A guaiacol is a phenol (aromatic ring with hydroxyl group) that has a methoxy group (-OCH3)
    in the ortho position (adjacent to the hydroxyl-bearing carbon).

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

    # SMARTS pattern for phenol (hydroxyl attached to aromatic carbon)
    phenol_pattern = Chem.MolFromSmarts("[OH]-[c]")
    if not mol.HasSubstructMatch(phenol_pattern):
        return False, "No phenolic hydroxyl group found"

    # SMARTS pattern for ortho methoxy substituent
    # Looks for: hydroxyl attached to aromatic carbon, adjacent to another aromatic carbon with methoxy (-O-CH3)
    ortho_methoxy_pattern = Chem.MolFromSmarts("[OH]-[c]@[c]-[OX2]-[CH3]")
    
    if mol.HasSubstructMatch(ortho_methoxy_pattern):
        return True, "Phenol with methoxy group in ortho position"
    else:
        return False, "No ortho methoxy substituent found on phenolic ring"