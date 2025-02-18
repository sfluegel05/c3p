"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: CHEBI:32876 tertiary amine
"""
from rdkit import Chem
from rdkit.Chem import HybridizationType

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine has a nitrogen atom bonded to three hydrocarbyl groups (carbon-based substituents).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen atom
            # Skip aromatic nitrogens (e.g., in pyridine/pyrrole)
            if atom.GetIsAromatic():
                continue
            
            # Check for SP3 hybridization (excludes conjugated/planar nitrogens)
            if atom.GetHybridization() != HybridizationType.SP3:
                continue
            
            # Check valence: exactly three bonds (excluding hydrogen)
            if atom.GetDegree() != 3:
                continue
            
            # Check all bonds are single bonds
            if any(bond.GetBondType() != Chem.BondType.SINGLE for bond in atom.GetBonds()):
                continue
            
            # Check formal charge (neutral amine)
            if atom.GetFormalCharge() != 0:
                continue
            
            # Verify all three substituents are carbon-based
            substituents = [neighbor.GetAtomicNum() for neighbor in atom.GetNeighbors()]
            if all(atomic_num == 6 for atomic_num in substituents):
                return True, "Tertiary amine with three carbon-based substituents"
    
    return False, "No tertiary amine group detected"