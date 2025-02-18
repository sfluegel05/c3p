"""
Classifies: CHEBI:18154 polysaccharide
"""
"""
Classifies: CHEBI:18154 polysaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Mol

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    A polysaccharide is a biomacromolecule with >10 monosaccharide residues linked glycosidically.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find potential sugar units: 5/6-membered rings with oxygen and multiple hydroxyls
    # More flexible pattern: any 5/6-membered ring with at least two oxygen atoms (including ring O)
    sugar_units = 0
    rings = mol.GetRingInfo().AtomRings()
    for ring in rings:
        if len(ring) not in (5,6):
            continue
        has_ring_o = any(mol.GetAtomWithIdx(a).GetAtomicNum() == 8 for a in ring)
        if not has_ring_o:
            continue
        # Count oxygen-containing groups (OH, O-linked, etc.)
        o_count = 0
        for a in ring:
            atom = mol.GetAtomWithIdx(a)
            if atom.GetAtomicNum() == 8:  # Oxygen in ring
                o_count +=1
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetIdx() not in ring:  # Exclude ring oxygen
                    o_count +=1
        if o_count >= 3:  # At least two oxygen-containing groups besides ring O
            sugar_units +=1
    
    if sugar_units <= 10:
        return False, f"Only {sugar_units} sugar units found (needs >10)"
    
    # Find glycosidic bonds (oxygen connecting two sugar rings)
    glycosidic_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        if bond.GetBeginAtom().GetAtomicNum() != 8 or bond.GetEndAtom().GetAtomicNum() != 8:
            continue
        # Check if oxygen is connecting two different rings
        a1 = bond.GetBeginAtom().GetIdx()
        a2 = bond.GetEndAtom().GetIdx()
        in_rings = [any(a1 in r for r in rings), any(a2 in r for r in rings)]
        if all(in_rings) and not any(a1 in r and a2 in r for r in rings):
            glycosidic_bonds +=1
    
    if glycosidic_bonds < sugar_units -1:  # At least n-1 bonds for linear chain
        return False, f"Insufficient glycosidic bonds ({glycosidic_bonds} for {sugar_units} units)"
    
    # Check molecular weight (approximate for >10 units)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1500:  # Lowered threshold to account for different sugars
        return False, f"Molecular weight ({mol_wt:.1f} Da) too low for polysaccharide"
    
    return True, f"Contains {sugar_units} sugar units with {glycosidic_bonds} glycosidic bonds"