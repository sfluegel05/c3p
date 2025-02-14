"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    A glycosaminoglycan is a polysaccharide containing a substantial proportion of aminomonosaccharide residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosaminoglycan, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify sugar rings (5 or 6-membered rings with oxygen)
    ring_info = mol.GetRingInfo()
    ring_atoms = ring_info.AtomRings()
    sugar_rings = []
    for ring in ring_atoms:
        if len(ring) == 5 or len(ring) == 6:
            oxygens_in_ring = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetSymbol() == 'O':
                    oxygens_in_ring += 1
            if oxygens_in_ring >= 1:
                sugar_rings.append(set(ring))

    if len(sugar_rings) == 0:
        return False, "No sugar rings found"

    # Identify amino sugars (sugar rings with attached amino groups)
    amino_sugar_count = 0
    for ring in sugar_rings:
        has_amino_group = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() == 'C':
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'N':
                        # Check if nitrogen is part of an amino group (-NH2 or -NH-)
                        if neighbor.GetDegree() <= 2:
                            has_amino_group = True
                            break
                if has_amino_group:
                    break
        if has_amino_group:
            amino_sugar_count += 1

    if amino_sugar_count == 0:
        return False, "No amino sugars found"

    # Check for glycosidic linkages between sugar rings
    # A glycosidic linkage is an oxygen atom connecting two sugar rings
    glycosidic_bonds = 0
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            if atom1.GetSymbol() == 'O' and atom2.GetSymbol() == 'C':
                rings1 = [ring for ring in sugar_rings if atom1.GetIdx() in ring]
                rings2 = [ring for ring in sugar_rings if atom2.GetIdx() in ring]
                if rings1 and rings2 and rings1 != rings2:
                    glycosidic_bonds += 1
            elif atom2.GetSymbol() == 'O' and atom1.GetSymbol() == 'C':
                rings1 = [ring for ring in sugar_rings if atom2.GetIdx() in ring]
                rings2 = [ring for ring in sugar_rings if atom1.GetIdx() in ring]
                if rings1 and rings2 and rings1 != rings2:
                    glycosidic_bonds += 1

    if glycosidic_bonds == 0:
        return False, "No glycosidic linkages between sugar units found"

    # Check for polysaccharide chain (more than 3 sugar units)
    total_sugars = len(sugar_rings)
    if total_sugars < 3:
        return False, f"Only {total_sugars} sugar units found, not enough for a polysaccharide"

    # Check proportion of amino sugars
    proportion = amino_sugar_count / total_sugars
    if proportion < 0.5:
        return False, f"Only {amino_sugar_count} out of {total_sugars} sugars are amino sugars"

    return True, f"Contains polysaccharide chain with {amino_sugar_count} amino sugars out of {total_sugars} sugar units"