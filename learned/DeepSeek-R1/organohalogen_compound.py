"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: organohalogen compound (CHEBI: orgHalogen)
"""
from rdkit import Chem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on the presence of at least one carbon-halogen bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organohalogen compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for any carbon (aromatic or aliphatic) bonded to a halogen (F, Cl, Br, I, At)
    halogen_bond_pattern = Chem.MolFromSmarts("[C,c]~[F,Cl,Br,I,At]")
    
    if mol.HasSubstructMatch(halogen_bond_pattern):
        return True, "Contains at least one carbon-halogen bond"
    else:
        return False, "No carbon-halogen bonds found"