"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: 3-hydroxy fatty acyl-CoA
Definition: A hydroxy fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A 
with the carboxy group of any 3-hydroxy fatty acid.
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    Requirements:
      1. It must contain an acyl thioester bond (i.e. R-C(=O)-S-R').
      2. The fatty acyl group (R) must be at least three carbons long such that, 
         counting the carbonyl carbon as position 1, position 3 carries an OH.
      3. It must contain a fragment diagnostic of Coenzyme A.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a 3-hydroxy fatty acyl-CoA, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that -OH groups can be detected
    mol = Chem.AddHs(mol)
    
    # Look for a thioester group: a carbonyl (C(=O)) bound to a sulfur (S)
    thioester_smarts = "[CX3](=O)S"  
    patt_thioester = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(patt_thioester)
    if not thioester_matches:
        return False, "No thioester (acyl-CoA) functional group found"
    
    # For each thioester match, see whether the acyl chain has an OH at the 3-position.
    # In an acyl-CoA, the fatty acyl group is R-C(=O)-S-CoA.
    # We number the acyl from the carbonyl carbon as C1. To get to C3, we must
    # follow the bond from the carbonyl to the acyl chain (the neighbor that is a carbon),
    # then one more bond along that chain.
    found_hydroxy = False
    for match in thioester_matches:
        # match[0] is the carbonyl carbon (C) from the pattern [CX3](=O)S
        carbonyl_atom = mol.GetAtomWithIdx(match[0])
        # Identify the neighbor that belongs to the acyl chain (exclude S and the doubly bonded O)
        acyl_neighbors = []
        for nbr in carbonyl_atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(carbonyl_atom.GetIdx(), nbr.GetIdx())
            # Exclude the oxygen from the C=O as well as the sulfur (which goes to CoA)
            if nbr.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2.0:
                continue
            if nbr.GetAtomicNum() == 16:  # sulfur atom, part of the thioester bond to CoA
                continue
            if nbr.GetAtomicNum() == 6:  # carbon; assume acyl chain atoms are carbons
                acyl_neighbors.append(nbr)
        if not acyl_neighbors:
            continue  # go to next thioester match
        
        # For each candidate acyl chain start (should normally be one)
        for alpha in acyl_neighbors:
            # From the alpha carbon (position 2), look for the next carbon (position 3)
            beta_candidates = []
            for nbr in alpha.GetNeighbors():
                if nbr.GetIdx() == carbonyl_atom.GetIdx():
                    continue
                if nbr.GetAtomicNum() == 6:
                    beta_candidates.append(nbr)
            # Ideally there is one continuous chain; if multiple branches, check all
            for beta in beta_candidates:
                # Check if beta (position 3) carries an OH group.
                # Look at neighbors of beta for an oxygen. We require that at least one such oxygen
                # is bound via a single bond and (after adding Hs) appears as an -OH group.
                has_OH = False
                for nbr in beta.GetNeighbors():
                    if nbr.GetAtomicNum() == 8:
                        # Check bond type is single; and check if hydrogen(s) is attached to that O atom.
                        bond = mol.GetBondBetweenAtoms(beta.GetIdx(), nbr.GetIdx())
                        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                            continue
                        # Check if this oxygen has an explicit hydrogen.
                        # (After AddHs, -OH groups should have at least one hydrogen neighbor)
                        for o_nbr in nbr.GetNeighbors():
                            if o_nbr.GetAtomicNum() == 1:
                                has_OH = True
                                break
                    if has_OH:
                        break
                if has_OH:
                    found_hydroxy = True
                    break  # we found one valid acyl chain with 3-OH
            if found_hydroxy:
                break
        if found_hydroxy:
            break

    if not found_hydroxy:
        return False, "No acyl chain with an OH at the 3-position was found"
    
    # Now check for the coenzyme A moiety.
    # Most acyl-CoA molecules share a common substructure fragment from Coenzyme A.
    # We use a SMARTS pattern that matches a fragment typically found in CoA.
    # (This pattern looks for, e.g., the SCCNC(=O)CCNC(=O) unit present in many CoA derivatives.)
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    patt_coa = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(patt_coa):
        return False, "Coenzyme A moiety not detected"
    
    return True, "Contains a thioester group with a 3-hydroxy fatty acyl chain and a CoA moiety"

# Example (you can test with one or more SMILES strings from the examples provided):
if __name__ == "__main__":
    smiles_example = "CCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    result, reason = is_3_hydroxy_fatty_acyl_CoA(smiles_example)
    print(result, "->", reason)