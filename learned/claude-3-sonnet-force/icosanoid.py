"""
Classifies: CHEBI:23899 icosanoid
"""
"""
Classifies: CHEBI:35976 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    Icosanoids are signaling molecules arising from oxidation of C20 essential fatty acids
    like arachidonic acid, EPA, and DGLA. They typically contain a cyclopentane/cyclopentene
    ring, multiple carbon-carbon double bonds, and long carbon chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an icosanoid, False otherwise
        str: Reason for classification
    """
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # Check for 20 carbon atoms
        if sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6) != 20:
            return False, "Does not contain 20 carbon atoms"

        # Look for cyclopentane/cyclopentene ring
        ring_pattern = Chem.MolFromSmarts("[CR5]")
        ring_matches = mol.GetSubstructMatches(ring_pattern)
        if not ring_matches:
            return False, "No cyclopentane/cyclopentene ring found"

        # Check for multiple carbon-carbon double bonds
        double_bond_pattern = Chem.MolFromSmarts("[C]=[C]")
        double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
        if len(double_bond_matches) < 3:
            return False, "Less than 3 carbon-carbon double bonds"

        # Check for long carbon chains
        chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
        chain_matches = mol.GetSubstructMatches(chain_pattern)
        if len(chain_matches) < 2:
            return False, "Missing long carbon chains"

        # Check for oxygen atoms
        if sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8) < 3:
            return False, "Too few oxygen atoms for an icosanoid"

        return True, "Contains cyclopentane/cyclopentene ring, multiple double bonds, long carbon chains, and oxygen atoms"

    except Exception as e:
        return False, f"Error processing molecule: {str(e)}"