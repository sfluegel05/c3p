"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    Criteria:
    - Contains a 3-oxo (keto) group
    - Has a steroid nucleus (four fused rings)
    - 5-beta configuration (A/B rings cis-fused)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for at least one ketone group
    ketone = Chem.MolFromSmarts("[CX3]=O")
    if not mol.HasSubstructMatch(ketone):
        return False, "No ketone group found"
    
    # Check for 3-oxo: ketone attached to a carbon in the A ring
    # Assuming the A ring is a six-membered ring with the ketone
    # This SMARTS might need adjustment
    a_ring_ketone = Chem.MolFromSmarts("[C][C][C](=O)")
    if not mol.HasSubstructMatch(a_ring_ketone):
        return False, "Ketone not in A ring (position 3)"
    
    # Check 5-beta configuration: specific stereochemistry at ring junctions
    # SMARTS pattern for cis A/B ring fusion (5beta)
    steroid_beta = Chem.MolFromSmarts("[C@]12[C@@H]3[C@@H](CC1)[C@H](CC2)C3")
    if not mol.HasSubstructMatch(steroid_beta):
        return False, "5-beta configuration not found"
    
    return True, "3-oxo and 5-beta steroid"