"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    Catechins are flavan-3-ols, having the core structure of two benzene rings
    linked by a pyran ring and a hydroxyl group at the 3-position of the pyran ring.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is a catechin, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core flavan-3-ol substructure: a chromane ring system with a hydroxyl at position 3
    # The 2nd ring (B ring) is a simple phenyl
    # No stereochemistry is enforced at this step.
    flavan_3ol_core = Chem.MolFromSmarts("C1Cc2cc(O)cc(c2O[C]1c3ccccc3)O")
    
    if not mol.HasSubstructMatch(flavan_3ol_core):
        return False, "Core flavan-3-ol structure not found."
    
    return True, "Molecule contains the core flavan-3-ol structure."