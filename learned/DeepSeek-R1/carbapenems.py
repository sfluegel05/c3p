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

    # Check bicyclo[3.2.0] system using fused rings of sizes 4 and 5
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Find pairs of fused rings with sizes 4 and 5
    fused_4_5 = False
    for i, r1 in enumerate(atom_rings):
        for j, r2 in enumerate(atom_rings[i+1:], i+1):
            common = set(r1) & set(r2)
            if len(common) >= 2:  # Fused rings share at least two atoms
                sizes = {len(r1), len(r2)}
                if sizes == {4, 5}:
                    # Verify beta-lactam atoms are in both rings
                    beta_match = mol.GetSubstructMatch(beta_lactam_pattern)
                    beta_atoms = set(beta_match)
                    if beta_atoms.intersection(r1) and beta_atoms.intersection(r2):
                        fused_4_5 = True
                        break
        if fused_4_5:
            break

    if fused_4_5:
        return True, "Contains bicyclo[3.2.0] system with beta-lactam core"
    else:
        return False, "Missing fused 4- and 5-membered rings forming bicyclo[3.2.0] system"