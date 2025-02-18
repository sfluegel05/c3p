"""
Classifies: CHEBI:46633 carbapenems
"""
"""
Classifies: carbapenems (CHEBI:60861)
Beta-lactam antibiotics with a carbapenem skeleton substituted at positions 3, 4, and 6.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    Carbapenems must contain a beta-lactam ring in a bicyclo[3.2.0]heptane system.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a carbapenem, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for beta-lactam ring (4-membered ring with N and C=O)
    beta_lactam_pattern = Chem.MolFromSmarts("[#7]1-&@[C](=O)-&@[C]-&@[C]-&@1")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring detected"

    # Check bicyclo[3.2.0] system using fused rings and bridge atoms
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Find rings that share bridge atoms (characteristic of bicyclo systems)
    fused_pairs = []
    for i, r1 in enumerate(atom_rings):
        for j, r2 in enumerate(atom_rings[i+1:], i+1):
            common = set(r1) & set(r2)
            if len(common) >= 2:  # Fused rings share at least two atoms
                fused_pairs.append((r1, r2))

    # Verify bicyclo[3.2.0] topology (at least two fused rings: 4-membered beta-lactam + 5-membered)
    for r1, r2 in fused_pairs:
        sizes = {len(r1), len(r2)}
        if {4,5} <= sizes or {4,6} <= sizes:  # Allow some variation
            # Verify beta-lactam is part of the fused system
            beta_atoms = set(mol.GetSubstructMatch(beta_lactam_pattern))
            if beta_atoms.intersection(r1) and beta_atoms.intersection(r2):
                return True, "Contains bicyclo[3.2.0] system with beta-lactam core"

    # Check characteristic carbapenem substitution pattern (positions 3,4,6)
    # Simplified check for common substituents (sulfur-containing side chains)
    sulfur_pattern = Chem.MolFromSmarts("[S]")
    if not mol.HasSubstructMatch(sulfur_pattern):
        return False, "Missing common sulfur-containing substituent"

    # Final check for molecular weight range (typical carbapenems: 300-500 Da)
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if not (300 < mw < 700):
        return False, f"Molecular weight {mw:.1f} outside typical range"

    return False, "Does not meet all carbapenem structural criteria"