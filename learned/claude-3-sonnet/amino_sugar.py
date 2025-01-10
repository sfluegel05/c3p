"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies: CHEBI:22044 amino sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is a sugar having one or more alcoholic hydroxy groups 
    replaced by substituted or unsubstituted amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for cyclic sugar pattern
    sugar_pattern = Chem.MolFromSmarts("[C]1[C][C]([OH0,NH0])[C]([OH0,NH0])[C]([OH0,NH0])[C]1")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar ring structure found"

    # Look for amino groups (both -NH2 and -NHR)
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")
    amide_pattern = Chem.MolFromSmarts("[NX3;H1;$(NC=O)]")
    
    has_amine = mol.HasSubstructMatch(amine_pattern)
    has_amide = mol.HasSubstructMatch(amide_pattern)
    
    if not (has_amine or has_amide):
        return False, "No amino or substituted amino groups found"

    # Count hydroxyl groups to confirm sugar nature
    hydroxyl_pattern = Chem.MolFromSmarts("[OH1]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    if hydroxyl_matches < 2:
        return False, "Too few hydroxyl groups for a sugar"

    # Verify cyclic nature and size
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No rings found"
    
    # Check ring sizes - sugars typically have 5 or 6-membered rings
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    if not any(size in [5,6] for size in ring_sizes):
        return False, "No suitable sugar ring size found"

    # Additional check for carbons and oxygens typical in sugars
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 5:
        return False, "Too few carbons for a sugar structure"
    if o_count < 4:
        return False, "Too few oxygens for a sugar structure"

    # Check for connectivity pattern typical in sugars
    # (each carbon should typically have at least one O or N)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            neighbors = [neighbor.GetAtomicNum() for neighbor in atom.GetNeighbors()]
            if 8 not in neighbors and 7 not in neighbors:
                connected_to_o_or_n = False
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() in [8, 7]:
                        connected_to_o_or_n = True
                        break
                if not connected_to_o_or_n:
                    continue  # Allow some carbons without O/N connections

    return True, "Contains sugar ring structure with amino group(s) replacing hydroxyl position(s)"