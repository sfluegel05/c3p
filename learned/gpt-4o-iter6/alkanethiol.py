"""
Classifies: CHEBI:47908 alkanethiol
"""
from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is characterized by a sulfanyl group (-SH) attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Thiol group pattern - we need to find sulfur atoms with a single hydrogen attached
    thiol_pattern = Chem.MolFromSmarts("[SHX1]")
    if not mol.HasSubstructMatch(thiol_pattern):
        return False, "No thiol (-SH) group found"
    
    # Check if thiol is attached to an alkyl (carbon and hydrogens forming a non-aromatic chain or ring)
    alkyl_pattern = Chem.MolFromSmarts("[CX4;!R]")
    thiol_matches = mol.GetSubstructMatches(thiol_pattern)
    
    for match in thiol_matches:
        sulfur_atom = mol.GetAtomWithIdx(match[0])
        neighbors = sulfur_atom.GetNeighbors()
        for atom in neighbors:
            if atom.GetSymbol() == "C" and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                # Found an SP3 carbon, typically indicative of an alkyl or basic alkane environment
                return True, "Contains sulfanyl group attached to an SP3 hybridized carbon (alkyl group)"
        
    return False, "Thiol group not attached to an appropriate alkyl (saturated carbon)"

# Example usage
print(is_alkanethiol("SCC(CC)C"))  # Example SMILES for 2-Methyl-1-butanethiol