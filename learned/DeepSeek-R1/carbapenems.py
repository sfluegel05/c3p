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
    Carbapenems must contain a beta-lactam ring in a bicyclo[3.2.0]hept-2-ene system.
    
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

    # Check bicyclo[3.2.0] system with 4- and 5-membered fused rings
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    fused_4_5 = False
    for i, r1 in enumerate(atom_rings):
        for j, r2 in enumerate(atom_rings[i+1:], i+1):
            common = set(r1) & set(r2)
            if len(common) >= 2:  # Fused rings share at least two atoms
                sizes = {len(r1), len(r2)}
                if sizes == {4, 5}:
                    # Identify which ring is 5-membered
                    five_membered = r1 if len(r1) == 5 else r2
                    four_membered = r2 if len(r2) == 4 else r1

                    # Check for double bond in 5-membered ring
                    has_double = False
                    for bond in mol.GetBonds():
                        if bond.GetBeginAtomIdx() in five_membered and bond.GetEndAtomIdx() in five_membered:
                            if bond.GetBondType() == Chem.BondType.DOUBLE:
                                has_double = True
                                break
                    if not has_double:
                        continue

                    # Check 5-membered ring contains no sulfur
                    has_sulfur = any(mol.GetAtomWithIdx(a).GetAtomicNum() == 16 for a in five_membered)
                    if has_sulfur:
                        continue

                    # Verify beta-lactam is in the 4-membered ring
                    beta_match = mol.GetSubstructMatch(beta_lactam_pattern)
                    beta_in_four = all(atom in four_membered for atom in beta_match)
                    if beta_in_four:
                        fused_4_5 = True
                        break
        if fused_4_5:
            break

    if fused_4_5:
        return True, "Contains bicyclo[3.2.0]hept-2-ene system with beta-lactam core"
    else:
        return False, "Missing required fused ring system with double bond in 5-membered ring"