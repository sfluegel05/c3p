from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cyclohexenylalkanol(smiles: str):
    """
    Determines if a molecule is a cyclohexenylalkanol (an alkanol substituted with a cyclohexenyl group)

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclohexenylalkanol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for presence of OH group
    if not any(atom.GetSymbol() == 'O' and atom.GetTotalNumHs() == 1 for atom in mol.GetAtoms()):
        return False, "No hydroxyl (OH) group found"

    # Check for cyclohexene ring
    rings = mol.GetRingInfo()
    cyclohexene_found = False
    
    for ring in rings.AtomRings():
        if len(ring) == 6:  # 6-membered ring
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            # Check if all atoms in ring are carbon
            if not all(atom.GetSymbol() == 'C' for atom in atoms):
                continue
                
            # Check for one double bond in ring
            ring_bonds = []
            for i in range(len(ring)):
                for j in range(i+1, len(ring)):
                    bond = mol.GetBondBetweenAtoms(ring[i], ring[j])
                    if bond is not None:
                        ring_bonds.append(bond)
            
            double_bonds = sum(1 for bond in ring_bonds if bond.GetBondType() == Chem.BondType.DOUBLE)
            if double_bonds == 1:
                cyclohexene_found = True
                break

    if not cyclohexene_found:
        return False, "No cyclohexene ring found"

    # Check if OH group is connected to an alkyl chain
    oh_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetTotalNumHs() == 1:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    oh_carbons.append(neighbor)

    if not oh_carbons:
        return False, "Hydroxyl group not connected to carbon"

    # Check if OH is part of an alkanol group
    for c in oh_carbons:
        if c.IsInRing():
            return False, "Hydroxyl group directly attached to ring"

    return True, "Molecule contains cyclohexene ring and alkanol group"
# Pr=1.0
# Recall=1.0