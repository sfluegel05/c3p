"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
from rdkit import Chem

def is_4_hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    A 4'-hydroxyflavanone is a hydroxyflavanone with a hydroxy substituent located at position 4'.

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

    # Define core flavanone structure; C3-C6 coupling, closed C-ring with carbonyl
    flavanone_core_smarts = "[#6]1=[#8]-[C@H]2C(c3ccccc3O2)=C(c4ccccc4)C1=O"
    flavanone_core_mol = Chem.MolFromSmarts(flayanone_core_smarts)
    if not mol.HasSubstructMatch(flayanone_core_mol):
        return False, "Missing core flavanone structure"
    
    # SMARTS for specific 4'-hydroxy group on the B-ring (phenyl group opposite the ketone)
    hydroxy_4_prime_position_smarts = "[#6](=[#8])[C@@H]1O[C@@H](C=C1)c2ccc(O)cc2"  # Hydroxy on B-ring's para position
    hydroxy_4_prime_mol = Chem.MolFromSmarts(hydroxy_4_prime_position_smarts)
    if mol.HasSubstructMatch(hydroxy_4_prime_mol):
        return True, "Molecule is a 4'-hydroxyflavanone"

    return False, "No 4'-hydroxy group found on flavanone B ring"