"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
#!/usr/bin/env python
"""
Classifies: 3-hydroxy fatty acyl-CoA
Definition: A hydroxy fatty acyl-CoA that results from the formal condensation of 
  the thiol group of coenzyme A with the carboxy group of any 3-hydroxy fatty acid.
Improved strategy:
  1. Parse the molecule and add explicit hydrogens.
  2. Identify the thioester group (a carbonyl directly bonded to a sulfur).
  3. For each thioester match, consider the carbonyl carbon (C1) and then select 
     the neighboring carbon (excluding sulfur and the carbonyl oxygen) as the α–carbon.
     Among the neighbors of this candidate (excluding C1), at least one must bear an –OH;
     this oxygen must be directly attached to a hydrogen (i.e. an -OH group) – this corresponds 
     to the 3–hydroxy (β–carbon) condition.
  4. Exclude any candidate that is part of the Coenzyme A fragment by pre‐identifying a characteristic 
     phosphate–ribose substructure.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES.
    
    Requirements:
      1. Must contain a thioester group, i.e., a carbonyl (C(=O)) directly bonded to a sulfur.
      2. The fatty acyl chain (the side attached to the carbonyl, not the CoA portion)
         must be long enough such that, counting the carbonyl carbon as C1, there is an α–carbon (C2)
         and a β–carbon (C3) that bears an –OH group.
      3. The molecule must contain a fragment diagnostic of Coenzyme A.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the heuristic for a 3-hydroxy fatty acyl-CoA, else False.
        str: Reason message.
    """
    # Parse the SMILES string and add explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Identify the Coenzyme A portion by its phosphate-ribose motif.
    coa_smarts = "COP(O)(=O)OP(O)(=O)OC"
    patt_coa = Chem.MolFromSmarts(coa_smarts)
    coa_matches = mol.GetSubstructMatches(patt_coa)
    coa_atoms = set()
    for match in coa_matches:
        coa_atoms.update(match)
    
    # Identify the thioester group using a SMARTS for a carbonyl bound to sulfur.
    # Here '[CX3](=O)[SX2]' finds a trigonal carbon (C(=O)) bound to a sulfur (S).
    thioester_smarts = "[CX3](=O)[SX2]"
    patt_thioester = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(patt_thioester)
    if not thioester_matches:
        return False, "No thioester (acyl-CoA) functional group found"
    
    # Helper function: Check if a given atom (candidate for β–carbon) bears an -OH group.
    def has_OH(atom):
        # Look for any oxygen neighbor that is connected to a hydrogen.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:  # Oxygen
                # Check if this oxygen directly bonds to a hydrogen.
                for o_nbr in nbr.GetNeighbors():
                    if o_nbr.GetAtomicNum() == 1:
                        return True
        return False
    
    found_acyl = False
    # Loop over each thioester match.
    for match in thioester_matches:
        # In "[CX3](=O)[SX2]", match[0] is the carbonyl carbon and match[1] is the sulfur.
        carbonyl = mol.GetAtomWithIdx(match[0])
        # Look for the acyl chain side. Of the neighbors of the carbonyl, skip:
        #   - The double-bonded oxygen and the sulfur atom.
        #   - Any atom that is part of the CoA fragment.
        for nbr in carbonyl.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(carbonyl.GetIdx(), nbr.GetIdx())
            if nbr.GetAtomicNum() == 16:  # sulfur, skip it.
                continue
            if nbr.GetAtomicNum() == 8 and bond.GetBondType() == rdchem.BondType.DOUBLE:
                continue
            if nbr.GetIdx() in coa_atoms:
                continue
            # Candidate for the α–carbon should be a carbon atom.
            if nbr.GetAtomicNum() != 6:
                continue
            alpha = nbr
            # Now, from the α–carbon, find a neighbor (other than the carbonyl) that could act as the β–carbon.
            for beta in alpha.GetNeighbors():
                if beta.GetIdx() == carbonyl.GetIdx():
                    continue
                if beta.GetAtomicNum() != 6:
                    continue
                # Check if the β–carbon directly bears an -OH group.
                if has_OH(beta):
                    found_acyl = True
                    break
            if found_acyl:
                break
        if found_acyl:
            break

    if not found_acyl:
        return False, "No acyl chain found with an -OH on the β–carbon (3-hydroxy condition)"
    
    # Finally, verify that the Coenzyme A fragment is present.
    if not mol.HasSubstructMatch(patt_coa):
        return False, "Coenzyme A moiety not detected (fragment pattern missing)"
    
    return True, ("Contains a thioester functional group with an acyl chain bearing an -OH on "
                   "the β–carbon and a CoA moiety")
    
if __name__ == "__main__":
    # Test examples from the provided list.
    test_smiles = [
        # (R)-3-hydroxyhexanoyl-CoA
        "CCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12",
        # (3R,9Z,12Z,15Z,18Z)-3-hydroxytetracosatetraenoyl-CoA
        "CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12",
        # (R)-carnitinyl-CoA betaine (which is not expected to fulfill the 3-hydroxy fatty acyl-CoA criteria)
        "CC(C)(COP([O-])(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)C[C@@H](O)C[N+](C)(C)C"
    ]
    
    for smi in test_smiles:
        result, reason = is_3_hydroxy_fatty_acyl_CoA(smi)
        print(result, "->", reason)