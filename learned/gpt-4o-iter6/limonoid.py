"""
Classifies: CHEBI:39434 limonoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    Limonoids are highly oxygenated triterpenoids with a specific core structure featuring a furan ring.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a limonoid, False otherwise
        str: Reason for classification
    """

    # Attempt to parse the SMILES string into a molecular structure
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of at least one furan ring
    furan_pattern = Chem.MolFromSmarts("c1ccoc1")
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring identified, atypical for limonoid"

    # Check for a generic steroid-like core
    steroid_pattern = Chem.MolFromSmarts("C1CCC2CCC3C4CCCC5C3C2C1C=C4CCC5")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Lacks steroid-like core structure"

    # The significant presence of ether groups can help indicate high oxygenation
    ether_pattern = Chem.MolFromSmarts("[OX2]C")
    ether_count = len(mol.GetSubstructMatches(ether_pattern))
    if ether_count < 3:
        return False, "Insufficient ether linkages typical of limonoids"

    # Consider molecular weight range
    mol_wt = Descriptors.MolWt(mol)
    if mol_wt < 450 or mol_wt > 1200:
        return False, "Molecular weight falls outside typical limonoid range"

    return True, "Molecule exhibits characteristics consistent with limonoid structural class"