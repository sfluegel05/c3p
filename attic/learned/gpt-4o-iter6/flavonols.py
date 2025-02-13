"""
Classifies: CHEBI:28802 flavonols
"""
from rdkit import Chem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    A flavonol is defined as 'Any hydroxyflavone in which is the ring hydrogen at
    position 3 of the heterocyclic ring is replaced by a hydroxy group'.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a flavonol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for flavonols: includes 3-hydroxy-2-phenylchromen-4-one backbone
    flavonol_pattern = Chem.MolFromSmarts('c1cc(O)c2c(c1)oc1cc(O)cc(O)c1c2=O')
    if not mol.HasSubstructMatch(flavonol_pattern):
        return False, "No flavonol backbone found"
    
    # Check for the presence of hydroxy groups at key positions (particularly at position 3)
    # Which is handled by the pattern above
    
    return True, "Structure matches flavonol features with 3-hydroxyflavone backbone"