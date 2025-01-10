"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
"""
Classifies: CHEBI:22677 aliphatic nitrile
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile is any nitrile derived from an aliphatic compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic nitrile, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for cyano group (C≡N)
    cyano_pattern = Chem.MolFromSmarts("[CX2]#[NX1]")
    if cyano_pattern is None:
        return False, "Error in SMARTS pattern"
        
    cyano_matches = mol.GetSubstructMatches(cyano_pattern)
    
    if not cyano_matches:
        return False, "No cyano (C≡N) group found"
    
    # Pattern for conjugated systems (C=C-C≡N or Ar-C≡N or C=N-C≡N)
    conjugated_pattern = Chem.MolFromSmarts("[$(C=C),$(c),$(C=N)]-[CH2]C#N")
    conjugated_matches = mol.GetSubstructMatches(conjugated_pattern) if conjugated_pattern else []
    
    # For each cyano group, check if it's a proper aliphatic nitrile
    for match in cyano_matches:
        nitrile_c = mol.GetAtomWithIdx(match[0])  # Get the carbon of C≡N
        
        # Get the carbon atom attached to the nitrile carbon
        for neighbor in nitrile_c.GetNeighbors():
            if neighbor.GetAtomicNum() != 7:  # Skip the nitrogen atom of C≡N
                # Check multiple conditions that would disqualify an aliphatic nitrile
                
                # 1. Carbon should not be aromatic
                if neighbor.GetIsAromatic():
                    continue
                    
                # 2. Carbon should not be part of a conjugated system
                if neighbor.GetIsConjugated():
                    continue
                    
                # 3. Carbon should not be sp or sp2 hybridized
                if neighbor.GetHybridization() != Chem.HybridizationType.SP3:
                    continue
                    
                # 4. Check if the carbon is part of certain functional groups that would make it non-aliphatic
                non_aliphatic_pattern = Chem.MolFromSmarts("[$(C=O),$(C=N),$(C=C)]-[CH2]C#N")
                if non_aliphatic_pattern and mol.HasSubstructMatch(non_aliphatic_pattern):
                    continue
                
                # 5. Make sure we're not part of a complex ring system
                ring_info = mol.GetRingInfo()
                if ring_info.IsAtomInRingOfSize(neighbor.GetIdx(), 6):
                    # If in a 6-membered ring, make sure it's not part of a complex/aromatic system
                    ring_atoms = ring_info.AtomRings()[0]
                    if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_atoms):
                        continue
                
                # If we passed all checks, this is an aliphatic nitrile
                return True, "Contains cyano group (C≡N) attached to an aliphatic carbon"
    
    return False, "No aliphatic nitrile groups found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:22677',
        'name': 'aliphatic nitrile',
        'definition': 'Any nitrile derived from an aliphatic compound.',
        'parents': ['CHEBI:33577']
    }
}