"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
from rdkit import Chem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    A 4'-hydroxyflavanone is expected to possess a flavanone structure with a 4'-hydroxy group on the phenyl ring.

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

    # Identification of core flavanone structure
    # This pattern accounts for the cyclohexanone structure attached to a benzene with an ether link, 
    # having potential variance in ring connections.
    flavanone_core_pattern = Chem.MolFromSmarts("c1ccc(cc1)C2COc3c(O)cccc3C2=O")
    if not mol.HasSubstructMatch(flavanone_core_pattern):
        return False, "No flavanone core structure found"

    # Identify 4'-hydroxy group attached to the phenyl ring particularly
    hydroxy_pattern = Chem.MolFromSmarts("c1ccc(O)cc1")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "4'-hydroxy group not identified"

    # Further check for stereochemistry alignment if needed based on examples (using chiral flag)
    chiral_centers = Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=False)
    if len(chiral_centers) < 1:
        return False, "Stereochemistry not adequately identified"

    return True, "Contains flavanone structure with 4'-hydroxy group and appropriate stereochemistry"