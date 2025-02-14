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
    
    # Define a relaxed SMARTS pattern for the glycerol backbone, phosphocholine and ether/ester groups.
    # The [C] allows any stereochemistry on the glycerol's C2
    pattern = Chem.MolFromSmarts("[CX4]-[OX2]-[CH2]-[C](-[OX2]-[CX3](=[OX1])-[CX4])-[CH2]-[OX2]-P(=O)(O)(OCC[N+](C)(C)C)[O-]")
    match = mol.GetSubstructMatch(pattern)

    if match:
        # Check for correct stereochemistry
        glycerol_c2 = match[3]
        glycerol_c2_atom = mol.GetAtomWithIdx(glycerol_c2)
        if not glycerol_c2_atom.HasProp('_Chirality') or glycerol_c2_atom.GetProp('_Chirality') != 'R':
            return False, "Incorrect stereochemistry at glycerol C2"
        
        # Obtain the atom indices for the alkyl and acyl groups for rotatable bond calculation.
        alkyl_carbon = match[0]
        acyl_carbon = match[6]

        # Function to recursively find the carbon chain attached to the specified carbon atom and return the number of rotatable bonds
        def get_carbon_chain_rotatable_bonds(mol, start_atom_idx):
            chain_atoms = []
            queue = [mol.GetAtomWithIdx(start_atom_idx)]
            visited = {mol.GetAtomWithIdx(start_atom_idx).GetIdx()}
            while len(queue) > 0:
                current_atom = queue.pop(0)
                chain_atoms.append(current_atom)
                for neighbor in current_atom.GetNeighbors():
                  if neighbor.GetIdx() not in visited and neighbor.GetAtomicNum()==6:
                    queue.append(neighbor)
                    visited.add(neighbor.GetIdx())

            rotatable_bonds = 0
            for atom in chain_atoms:
               for bond in atom.GetBonds():
                  if bond.GetBeginAtom().GetIdx() in visited and bond.GetEndAtom().GetIdx() in visited and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                      rotatable_bonds += 1

            return rotatable_bonds


        alkyl_rotatable_bonds = get_carbon_chain_rotatable_bonds(mol, alkyl_carbon)
        acyl_rotatable_bonds = get_carbon_chain_rotatable_bonds(mol, acyl_carbon)

        if alkyl_rotatable_bonds < 3 or acyl_rotatable_bonds < 3:
              return False, "Alkyl and/or Acyl chains too short."

        return True, "Matches 2-acyl-1-alkyl-sn-glycero-3-phosphocholine criteria"
    else:
        return False, "Does not match 2-acyl-1-alkyl-sn-glycero-3-phosphocholine structure"