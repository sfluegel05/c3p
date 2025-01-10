"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.
    A 3'-hydroxyflavanone has a flavanone core structure with a hydroxy group
    on the phenyl ring specifically at the 3' position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the flavanone core structure
    # Core pattern with stereochemistry considerations and connectivity flexibility
    flavanone_core_pattern = Chem.MolFromSmarts("C1=COC2=C1C(=O)C(C[C@H](C2)O)C3=CC=CC=C3")
    if not mol.HasSubstructMatch(flavanone_core_pattern):
        return False, "Flavanone core structure not found"

    # Find the phenyl ring attached (ensure it's correctly matched)
    phenyl_ring_pattern = Chem.MolFromSmarts("c1ccccc1")
    phenyl_matches = mol.GetSubstructMatches(phenyl_ring_pattern)
    
    # Check for the hydroxy group at the 3' position
    hydroxy_pattern = Chem.MolFromSmarts("c1cc(O)ccc1")
    
    # Iterate through each potential phenyl group
    for match in phenyl_matches:
        phenyl_mol = Chem.PathToSubmol(mol, match)
        if phenyl_mol.HasSubstructMatch(hydroxy_pattern):
            return True, "Contains flavanone core with a hydroxy group at 3' position"

    return False, "No hydroxy group at 3' position on any phenyl ring"