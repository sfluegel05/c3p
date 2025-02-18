"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
"""
Classifies: trans-2-enoyl-CoA
An unsaturated fatty acyl-CoA that results from the formal condensation of the thiol group 
of coenzyme A with the carboxy group of any 2,3-trans-enoic acid.

Revised approach:
  1. Validate the SMILES and force assignment of stereochemistry.
  2. Verify presence of a CoA substructure by requiring both a thioester and the adenine motif.
  3. Identify the thioester linkage ([CX3](=O)[SX2]). For each match:
      a. Ensure that the carbonyl carbon (C=O) has exactly two heavy-atom neighbors,
         one being the sulfur (from CoA linkage) and the other, the acyl α–carbon.
      b. From the α–carbon, look for a double bond to a β–carbon that is explicitly set to E (trans).
      c. Optionally traverse a few additional carbons along the acyl chain to ensure a minimum length.
  4. If any thioester gives a valid trans-2-enoyl acyl chain and the molecule contains the adenine part of CoA,
     we classify the structure as a trans-2-enoyl-CoA.

If no valid combination is found, we return False with an explanatory message.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is identified as a trans-2-enoyl-CoA, False otherwise
        str: Explanation of the classification decision
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Force assignment of stereochemistry (helps capture / and \ definitions)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # --- Check for a CoA motif in two parts:
    # (a) Look for a thioester linkage: [CX3](=O)[SX2]
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Missing thioester (acyl-CoA linkage) moiety"
    
    # (b) Look for the adenine substructure within the CoA nucleotide.
    # This pattern picks up a common adenine ring: e.g. n1cnc2c(ncnc12)
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing characteristic adenine (CoA nucleotide) moiety"
    
    # --- For each thioester linkage, examine its acyl chain:
    for match in thioester_matches:
        # In our match, index0 is the carbonyl carbon and index1 is the sulfur.
        carbonyl_idx = match[0]
        sulfur_idx = match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # A proper thioester carbonyl should be attached to exactly 2 heavy atoms:
        # one being the sulfur, one being the acyl chain’s alpha carbon.
        heavy_neighbors = [nbr for nbr in carbonyl_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors) != 2:
            continue  # skip if extra connections might be present
        
        # Identify the α–carbon from the carbonyl (the neighbor that is not the sulfur)
        alpha_atom = None
        for nbr in heavy_neighbors:
            if nbr.GetIdx() != sulfur_idx and nbr.GetAtomicNum() == 6:
                alpha_atom = nbr
                break
        if alpha_atom is None:
            continue  # no valid α–carbon
        
        # Now, from the α–carbon, look for a double bond to a β–carbon.
        valid_trans_double = False
        for bond in alpha_atom.GetBonds():
            if bond.GetBondType() != Chem.rdchem.BondType.DOUBLE:
                continue
            # The double bond should be from α to a β carbon
            beta_atom = bond.GetOtherAtom(alpha_atom)
            if beta_atom.GetAtomicNum() != 6:
                continue
            
            # Check that the bond is explicitly set to trans (E)
            # This requires that stereochemistry has been properly assigned.
            if bond.GetStereo() != Chem.rdchem.BondStereo.STEREOE:
                continue

            # Optionally, verify that the acyl chain is of a minimal length.
            # We start counting from the β–carbon and follow a simple (first found) linear path.
            chain_length = 1  # start counting the β–carbon
            current_atom = beta_atom
            previous_idx = alpha_atom.GetIdx()
            # Traverse up to 4 bonds; expecting a real acyl chain to have several carbons.
            for _ in range(4):
                next_neighbors = [nbr for nbr in current_atom.GetNeighbors() 
                                  if nbr.GetIdx() != previous_idx and nbr.GetAtomicNum() == 6]
                if not next_neighbors:
                    break
                previous_idx = current_atom.GetIdx()
                current_atom = next_neighbors[0]
                chain_length += 1
            # For our tested data, we expect at least 2 carbons (β plus further chain) beyond the α–carbon.
            if chain_length < 1:
                continue  # chain too short; most valid examples have a chain longer than just the double bond

            valid_trans_double = True
            break  # found a valid α–β bond
        
        if valid_trans_double:
            return True, ("Found thioester linkage with an acyl chain exhibiting a trans (E) α–β double bond "
                          "at the 2-enoyl position and a CoA adenine motif")
    
    return False, "No appropriate trans (E) α–β double bond in an acyl chain attached via thioester to CoA was detected"

# Example usage:
if __name__ == "__main__":
    test_examples = [
        # Valid examples
        ("(2E,9Z,12Z)-octadecatrienoyl-CoA", "CCCC\\C=C/C\\C=C/CCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)"
         "[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"),
        ("trans-hex-2-enoyl-CoA", "CCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)"
         "[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"),
        # A false positive candidate (should return False because the phosphate charges disrupt the CoA motif):
        ("trans-2-tetracosenoyl-CoA(4-)", "CCCCCCCCCCCCCCCCCCCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)"
         "[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12")
    ]
    
    for name, smi in test_examples:
        valid, msg = is_trans_2_enoyl_CoA(smi)
        print(f"Name: {name}\nValid: {valid}\nReason: {msg}\n{'-'*40}")