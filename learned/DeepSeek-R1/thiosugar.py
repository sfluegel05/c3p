"""
Classifies: CHEBI:73754 thiosugar
"""
"""
Classifies: CHEBI:79126 thiosugar
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar is a carbohydrate derivative with oxygen/hydroxyl groups replaced by sulfur.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Basic carbohydrate check: multiple hydroxyl groups and ring structure
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() 
                        if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1)
    if hydroxyl_count < 3:
        return False, "Insufficient hydroxyl groups for carbohydrate"

    # Check for sulfur presence
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    if not sulfur_atoms:
        return False, "No sulfur atoms found"

    # Improved carbohydrate backbone detection
    sugar_ring = False
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) not in (5, 6):
            continue
            
        # Count hydroxyls on ring carbons
        ring_oh = 0
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:
                continue
            for bond in atom.GetBonds():
                if bond.GetBondType() != Chem.BondType.SINGLE:
                    continue
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
                    ring_oh += 1
                    break  # Count one OH per carbon
                    
        if ring_oh >= 3:
            sugar_ring = True
            # Check for sulfur in ring structure
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() == 16:
                    return True, "Sulfur in carbohydrate ring"
            # Check sulfur attachments to ring carbons
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() != 6:
                    continue
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 16:
                        return True, "Sulfur attached to ring carbon"

    # Check for sulfur in glycosidic/thioether positions
    glycosidic_s = Chem.MolFromSmarts("[C][SX2][C]")
    if mol.HasSubstructMatch(glycosidic_s):
        return True, "Sulfur in glycosidic linkage"

    # Check for sulfur in sidechains of potential carbohydrate structure
    for sulfur in sulfur_atoms:
        # Look for S connected to carbons with multiple O neighbors
        for bond in sulfur.GetBonds():
            carbon = bond.GetOtherAtom(sulfur)
            if carbon.GetAtomicNum() != 6:
                continue
            # Check if carbon has carbohydrate-like features
            o_neighbors = sum(1 for n in carbon.GetNeighbors() 
                            if n.GetAtomicNum() == 8 and n.GetTotalNumHs() >= 1)
            if o_neighbors >= 2:
                return True, "Sulfur in carbohydrate sidechain"

    return False, "No sulfur substitution in carbohydrate structure"