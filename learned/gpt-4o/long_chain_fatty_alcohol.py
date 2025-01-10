"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    A long-chain fatty alcohol has a chain length ranging from C13 to C22 with a hydroxyl group attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for hydroxyl group (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"

    # Calculate the number of carbon atoms in the longest carbon chain
    c_count = 0
    ssr = Chem.GetSymmSSSR(mol)  # Get smallest set of smallest rings
    for ring in ssr:
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:  # Count carbon atoms
                c_count += 1

    # Check if carbon count is within 13 to 22
    if c_count < 13 or c_count > 22:
        return False, f"Carbon chain length is {c_count}, expected between 13 and 22"
    
    return True, "Has a hydroxyl group and the carbon chain length is within range for long-chain fatty alcohol"