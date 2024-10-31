from rdkit import Chem
from rdkit.Chem import AllChem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid derivative.
    Quinic acid is a cyclitol carboxylic acid with a cyclohexane ring substituted with 
    hydroxy groups and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinic acid derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of cyclohexane ring
    rings = mol.GetRingInfo()
    has_cyclohexane = False
    for ring in rings.AtomRings():
        if len(ring) == 6:
            ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
            # Check if ring is saturated (sp3) and made of carbons
            if all(atom.GetHybridization() == Chem.HybridizationType.SP3 and 
                  atom.GetSymbol() == 'C' for atom in ring_atoms):
                has_cyclohexane = True
                cyclohexane_atoms = set(ring)
                break
    
    if not has_cyclohexane:
        return False, "No cyclohexane ring found"

    # Check for carboxylic acid group attached to cyclohexane ring
    has_carboxyl = False
    for atom_idx in cyclohexane_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                o_count = 0
                oh_count = 0
                for n2 in neighbor.GetNeighbors():
                    if n2.GetSymbol() == 'O':
                        if n2.GetTotalNumHs() == 1:
                            oh_count += 1
                        else:
                            o_count += 1
                if o_count == 1 and oh_count == 1:
                    has_carboxyl = True
                    break
    
    if not has_carboxyl:
        return False, "No carboxylic acid group found attached to cyclohexane ring"

    # Check for hydroxyl groups on cyclohexane ring
    hydroxy_count = 0
    for atom_idx in cyclohexane_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == 'O' and neighbor.GetTotalNumHs() == 1:
                hydroxy_count += 1

    if hydroxy_count < 1:
        return False, "Insufficient hydroxyl groups on cyclohexane ring"

    # Additional check for ester derivatives
    ester_count = 0
    for atom_idx in cyclohexane_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == 'O':
                for n2 in neighbor.GetNeighbors():
                    if n2.GetSymbol() == 'C' and any(n3.GetSymbol() == 'O' for n3 in n2.GetNeighbors()):
                        ester_count += 1

    return True, f"Quinic acid derivative with {hydroxy_count} hydroxyl groups and {ester_count} ester groups"
# Pr=1.0
# Recall=0.8571428571428571