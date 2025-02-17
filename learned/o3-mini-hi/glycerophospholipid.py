"""
Classifies: CHEBI:37739 glycerophospholipid
"""
#!/usr/bin/env python
"""
Classifies: glycerophospholipid
Definition:
  Any glycerolipid having a phosphate group ester-linked to a terminal carbon of the glycerol backbone.
Heuristic criteria used:
  1. The molecule must contain a phosphorus atom bonded via an oxygen to a primary (CH2) carbon.
  2. The molecule should contain at least one acyl ester group (OC(=O)) aside from the phosphate ester.
Note: Because of the chemical diversity in glycerophospholipid structures, this is a heuristic filter.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycerophospholipid(smiles: str):
    """
    Determines if a molecule is a glycerophospholipid based on its SMILES string.
    A glycerophospholipid is defined as any glycerolipid having a phosphate group 
    ester-linked to a terminal carbon (CH2) of the glycerol backbone.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a glycerophospholipid; False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # ------------ Criterion 1: Look for a phosphate group attached via an oxygen to a primary carbon ------------
    # Iterate through all bonds to find a bond where one atom is oxygen and the other is phosphorus.
    phosphate_linked = False
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Check if one atom is O and the other P
        if (a1.GetSymbol() == 'O' and a2.GetSymbol() == 'P') or (a1.GetSymbol() == 'P' and a2.GetSymbol() == 'O'):
            # Identify the oxygen (o_atom) and phosphorus (p_atom)
            o_atom = a1 if a1.GetSymbol() == 'O' else a2
            p_atom = a1 if a1.GetSymbol() == 'P' else a2
            
            # We also want to ensure that the phosphorus atom looks like part of a phosphate. 
            # A common pattern is P(=O)(O)(O), so at least one double bond to O should be present.
            p_neighbors = [n for n in p_atom.GetNeighbors() if n.GetSymbol() == 'O']
            double_bonded_count = 0
            for n in p_neighbors:
                bond_ptype = mol.GetBondBetweenAtoms(p_atom.GetIdx(), n.GetIdx())
                if bond_ptype.GetBondType() == Chem.BondType.DOUBLE:
                    double_bonded_count += 1
            if double_bonded_count < 1:
                continue  # not a typical phosphate environment; skip

            # Now check the oxygen atom.
            # Look for an oxygen that is also attached to a carbon.
            for nbr in o_atom.GetNeighbors():
                if nbr.GetSymbol() == 'C':
                    # Check if this carbon is primary (i.e. a CH2 group) 
                    # (RDKit might not explicitly have Hs -- use GetTotalNumHs)
                    if nbr.GetTotalNumHs() == 2:
                        phosphate_linked = True
                        break
            if phosphate_linked:
                break
    if not phosphate_linked:
        return False, "No phosphate group found ester-linked to a primary carbon (CH2) as expected for a glycerol backbone"
    
    # ------------ Criterion 2: Check for the presence of acyl ester groups (OC(=O)) ------------
    # The acyl ester group is represented by the SMARTS: oxygen bound to a carbonyl carbon.
    ester_smarts = "O[C;!$(OP(*))](=O)"  # O attached to C(=O); exclude those already part of the phosphate.
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No acyl ester (OC(=O)) group found that would indicate a fatty-acyl chain"
    
    # ------------ (Optional) Additional filters ------------
    # For example, one may wish to require a minimum molecular weight or count of rotatable bonds
    # to reinforce the expectation for a lipid. These criteria can be commented/uncommented as needed.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight is too low to be a glycerophospholipid"
    
    # If all criteria are met, we classify the molecule as a glycerophospholipid.
    return True, "Molecule contains a phosphate group estered to a primary (CH2) carbon and an acyl ester group, consistent with a glycerophospholipid structure"

# Example usage (for testing):
if __name__ == "__main__":
    # SMILES example from one of the provided glycerophospholipid examples
    test_smiles = "CCCCCCCCCCCCCCCC(=O)OCC(COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCC(=O)\\C=C\\C(O)=O"
    result, reason = is_glycerophospholipid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)