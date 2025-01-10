"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.

    A 4'-hydroxyflavanone has a core flavanone structure with a 4'-hydroxy group on the phenyl ring. 

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 4'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for flavanone core: benzene ring with ketone and cyclic ether linkage
    flavanone_core_pattern = Chem.MolFromSmarts("c1ccc(cc1)-C2-C(=O)-c3ccccc3O2")
    
    if not mol.HasSubstructMatch(flavanone_core_pattern):
        return False, "No core flavanone structure found"

    # Pattern for 4'-hydroxy group on the phenyl ring
    hydroxy_pattern = Chem.MolFromSmarts("c1cc(O)ccc1")
    
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No hydroxy group at the 4' position"

    return True, "Contains core flavanone structure with 4'-hydroxy group"