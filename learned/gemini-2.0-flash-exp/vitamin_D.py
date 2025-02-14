"""
Classifies: CHEBI:27300 vitamin D
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D based on its SMILES string.
    Vitamin D is a seco-steroid with a characteristic broken B-ring, triene system,
    and hydroxyl groups at specific positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a vitamin D, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for broken B-ring - Look for a 6-membered ring with a double bond connected to a single bond, and a specific connection to the other rings.
    # The key is the presence of a C=C single bond within a 6 member ring connected to two other carbons that will continue the steroidal framework
    # We use a complex SMARTS pattern to enforce this connectivity, while being as general as possible.
    broken_b_ring_pattern = Chem.MolFromSmarts("[C;R6]1([C;R](=[C;R])[C;R])~[C;R]~[C;R]~[C;R]~[C;R]1")
    if not mol.HasSubstructMatch(broken_b_ring_pattern):
        return False, "Not a seco-steroid structure (missing broken B ring)"

    # 2. Check for triene system in the conjugated system
    triene_pattern = Chem.MolFromSmarts("C=C-C=C-C") # this is still /C=C/C=C/C
    if not mol.HasSubstructMatch(triene_pattern):
        return False, "No triene system found"
    
    # 3. Check for at least one hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[CH0;R][OH1]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 1:
        return False, "No hydroxyl group found"

    # 4. Check for a side chain of at least 3 carbons
    side_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    side_chain_matches = mol.GetSubstructMatches(side_chain_pattern)
    if len(side_chain_matches) < 1 :
        return False, "No long enough side chain"

    # 5. Check molecular weight range - Vitamin D typically between 350 and 600 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350 or mol_wt > 600:
        return False, f"Molecular weight out of range: {mol_wt}"


    return True, "Meets vitamin D structural criteria"