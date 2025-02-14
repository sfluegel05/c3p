"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the core SMARTS pattern: Glycerol with alkyl ether at C1, acyl ester at C2 and phosphocholine at C3.
    # Note that [C@H] specifies the R stereochemistry
    # This pattern enforces proper stereochemistry and connectivity
    pattern = Chem.MolFromSmarts("[CX4]-[OX2]-[CH2][C@H](-[OX2]-[CX3](=[OX1])-[CX4])-[CH2]-[OX2]-P(=O)(O)(OCC[N+](C)(C)C)[O-]")

    
    match = mol.GetSubstructMatch(pattern)

    if match:
       #Check if the chain length of alkyl and acyl is greater than 3
        glycerol_c1 = match[2]
        glycerol_c2 = match[3]
        ether_oxygen = match[1]
        ester_oxygen = match[4]
        carbonyl_carbon = match[5]
        alkyl_carbon = match[0]
        acyl_carbon = match[6]


        alkyl_chain_atoms = []
        queue = [mol.GetAtomWithIdx(alkyl_carbon)]
        visited = {mol.GetAtomWithIdx(alkyl_carbon).GetIdx()}
        while len(queue) > 0:
            current_atom = queue.pop(0)
            alkyl_chain_atoms.append(current_atom)
            for neighbor in current_atom.GetNeighbors():
                if neighbor.GetIdx() not in visited and neighbor.GetAtomicNum()==6:
                    queue.append(neighbor)
                    visited.add(neighbor.GetIdx())
        
        alkyl_rotatable_bonds = 0
        for atom in alkyl_chain_atoms:
            for bond in atom.GetBonds():
                if bond.GetBeginAtom().GetIdx() in visited and bond.GetEndAtom().GetIdx() in visited and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    alkyl_rotatable_bonds += 1


        acyl_chain_atoms = []
        queue = [mol.GetAtomWithIdx(acyl_carbon)]
        visited = {mol.GetAtomWithIdx(acyl_carbon).GetIdx()}
        while len(queue) > 0:
            current_atom = queue.pop(0)
            acyl_chain_atoms.append(current_atom)
            for neighbor in current_atom.GetNeighbors():
              if neighbor.GetIdx() not in visited and neighbor.GetAtomicNum()==6:
                queue.append(neighbor)
                visited.add(neighbor.GetIdx())

        acyl_rotatable_bonds = 0
        for atom in acyl_chain_atoms:
            for bond in atom.GetBonds():
              if bond.GetBeginAtom().GetIdx() in visited and bond.GetEndAtom().GetIdx() in visited and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                  acyl_rotatable_bonds += 1
        
        if alkyl_rotatable_bonds < 3 or acyl_rotatable_bonds < 3:
            return False, "Alkyl and/or Acyl chains too short."


        return True, "Matches 2-acyl-1-alkyl-sn-glycero-3-phosphocholine criteria"
    else:
        return False, "Does not match 2-acyl-1-alkyl-sn-glycero-3-phosphocholine structure"