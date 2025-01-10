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
    
    # Update SMARTS pattern for flavonols to capture 3-hydroxyflavone backbone more generally
    # Revised pattern to not restrict overly and permit flavonol central structure
    flavonol_pattern = Chem.MolFromSmarts('c1c(oc2c(c1)ccc(c2=O)O)c1ccccc1')
    hydroxy_position_pattern = Chem.MolFromSmarts('[OH]c1oc2c(ccc(c2=O)O)c(c1)')
    
    if not mol.HasSubstructMatch(flavonol_pattern):
        return False, "No flavonol backbone found"

    # Ensure the presence of a hydroxy group at the 3 position in the core
    if not mol.HasSubstructMatch(hydroxy_position_pattern):
        return False, "No hydroxy group at the 3 position in the flavonol structure"

    return True, "Structure matches flavonol features with 3-hydroxyflavone backbone"

# This approach should cover more cases as it focuses on essential structural markers for flavonols.