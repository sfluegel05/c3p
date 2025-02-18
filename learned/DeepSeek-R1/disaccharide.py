"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: CHEBI:15953 disaccharide
"""
from rdkit import Chem
from rdkit.Chem import Mol

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide consists of two monosaccharide units joined by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
    except:
        return False, "Invalid SMILES"
    if not mol:
        return False, "Invalid SMILES"
    
    ring_info = mol.GetRingInfo()
    valid_bonds = 0

    # Iterate through potential bridging oxygens
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        
        # Oxygen must not be in any ring
        if any(atom.GetIdx() in ring for ring in ring_info.AtomRings()):
            continue

        # Must have exactly two neighbors
        if len(atom.GetNeighbors()) != 2:
            continue

        # Create modified molecule without the oxygen
        modified = Chem.RWMol(mol)
        modified.RemoveAtom(atom.GetIdx())
        parts = list(Chem.GetMolFrags(modified.GetMol(), asMols=True, sanitizeFrags=False))

        if len(parts) != 2:
            continue

        valid = True
        for part in parts:
            try:
                Chem.SanitizeMol(part)
            except:
                valid = False
                break

            # Check for at least one 5/6-membered ring with oxygen
            part_ring_info = part.GetRingInfo()
            found_valid_ring = False
            for ring in part_ring_info.AtomRings():
                if len(ring) not in (5,6):
                    continue
                if any(part.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring):
                    found_valid_ring = True
                    break
            if not found_valid_ring:
                valid = False

            # Check minimum oxygen count (2 for deoxy sugars)
            if sum(1 for a in part.GetAtoms() if a.GetAtomicNum() == 8) < 2:
                valid = False

        if valid:
            valid_bonds += 1

    if valid_bonds == 1:
        return True, "Two monosaccharide units connected by a glycosidic bond"
    else:
        return False, f"Found {valid_bonds} valid glycosidic bonds, need exactly 1"