"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: alkanesulfonate oxoanion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion has a sulfonate group (-SO3-) attached to an aliphatic carbon.
    The carbon at position 1 can be attached to hydrogens, a carbon chain, or other groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for aliphatic carbon-sulfonate pattern: C-S([O-])(=O)=O 
    # where C is not aromatic
    alkanesulfonate_pattern = Chem.MolFromSmarts("[CX4;!$(C:*)]-[SX4](=[OX1])(=[OX1])[O-]")
    if not mol.HasSubstructMatch(alkanesulfonate_pattern):
        return False, "No aliphatic carbon-sulfonate group found"

    # Get matches
    matches = mol.GetSubstructMatches(alkanesulfonate_pattern)
    
    for match in matches:
        carbon_idx = match[0]
        sulfur_idx = match[1]
        
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
        
        # Verify carbon is not aromatic
        if carbon_atom.GetIsAromatic():
            continue
            
        # Verify sulfonate group has correct charges
        if sulfur_atom.GetFormalCharge() != 0:
            continue
            
        # Count oxygen neighbors of sulfur with correct properties
        oxygen_count = 0
        double_bond_count = 0
        negative_charge_count = 0
        
        for neighbor in sulfur_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                oxygen_count += 1
                if mol.GetBondBetweenAtoms(sulfur_idx, neighbor.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
                    double_bond_count += 1
                if neighbor.GetFormalCharge() == -1:
                    negative_charge_count += 1
        
        # Must have 3 oxygens: 2 double-bonded and 1 with negative charge
        if oxygen_count != 3 or double_bond_count != 2 or negative_charge_count != 1:
            continue
            
        # Check that carbon is not part of a conjugated system
        conjugated = False
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                bond = mol.GetBondBetweenAtoms(carbon_idx, neighbor.GetIdx())
                if bond.GetBondType() != Chem.BondType.SINGLE:
                    conjugated = True
                    break
                    
        if conjugated:
            continue
            
        return True, "Contains aliphatic carbon-sulfonate group (-CH2-SO3-)"

    return False, "No valid alkanesulfonate group found"