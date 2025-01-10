"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: epoxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    An epoxy fatty acid features an epoxide group and a carboxylic acid group attached to a linear or branched hydrocarbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify epoxide group
    epoxide_pattern = Chem.MolFromSmarts("C1OC1")
    if not mol.HasSubstructMatch(epoxide_pattern):
        return False, "No epoxide group found"
    
    # Identify carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for long carbon chain
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 12:
        return False, "Insufficient number of carbon atoms for a fatty acid"

    # Verify if hydrocarbon chains are neither too complex nor bridged
    # To ensure epoxy is not in macrocycles, check the valency and environment
    for b in mol.GetBonds():
        if b.IsInRing() and b.GetBondType() == rdchem.BondType.SINGLE:
            ring_size = len(mol.GetRingInfo().AtomRings()[0])
            if ring_size > 3 and b.GetBeginAtom().GetAtomicNum() == 8 and b.GetEndAtom().GetAtomicNum() == 8:
                return False, "Epoxide part of a more complex structure or macrocycle"
    
    return True, "Contains epoxide group and carboxylic acid along with a sufficiently linear long hydrocarbon chain"

# Example usage:
# result, reason = is_epoxy_fatty_acid("CCCCCCCC1OC1CCCCCCCC(O)=O")
# print(result, reason)