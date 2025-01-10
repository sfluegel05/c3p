"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
"""
Classifies: CHEBI:35341 3-oxo-5beta-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    A 3-oxo-5beta-steroid has a ketone at position 3 and beta configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5beta-steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES with stereochemistry
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic steroid core pattern (four fused rings)
    steroid_core = Chem.MolFromSmarts("C1CC2CCC3C(CCC4CCCC43C)C2C1")
    if steroid_core is None:
        return False, "Invalid steroid core SMARTS pattern"
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Check for ketone at position 3
    # Pattern for 3-keto group in context of the A ring
    ketone_pattern = Chem.MolFromSmarts("C1CC(=O)CC2")
    if ketone_pattern is None:
        return False, "Invalid ketone SMARTS pattern"
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone group at position 3"

    # Check for 5-beta configuration
    # In 5-beta steroids, rings A/B are cis-fused
    # The hydrogen at C5 is alpha (pointing up), making the fusion beta
    # We use multiple patterns to catch different representations
    beta_patterns = [
        # Pattern looking for the A/B ring junction with explicit beta stereochemistry
        "[H][C@@]12CC[C@@H]",  # C5 with alpha H (meaning beta fusion)
        "[C@@H]1[C@@H]2CC",    # Alternative representation
        "[C@]1([H])[C@@H]2"    # Another alternative
    ]
    
    beta_matches = False
    for pattern in beta_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat is not None and mol.HasSubstructMatch(pat, useChirality=True):
            beta_matches = True
            break
    
    if not beta_matches:
        return False, "No 5-beta configuration found"

    # Additional checks for steroid-like properties
    
    # Count carbons (steroids typically have 19-27 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 19 or carbon_count > 30:
        return False, f"Carbon count ({carbon_count}) outside typical steroid range (19-30)"

    # Check molecular weight
    mol_weight = Chem.Descriptors.ExactMolWt(mol)
    if mol_weight < 250 or mol_weight > 500:
        return False, f"Molecular weight {mol_weight:.1f} outside typical steroid range (250-500)"

    # Check ring count (steroids should have 4 main rings)
    ring_count = len(Chem.GetSymmSSSR(mol))
    if ring_count < 4:
        return False, f"Too few rings ({ring_count}) for a steroid structure"

    # If all checks pass
    return True, "Molecule contains steroid core with 3-oxo group and 5-beta configuration"