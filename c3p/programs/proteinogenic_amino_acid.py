"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
#!/usr/bin/env python3
"""
Classifies: proteinogenic amino acid
Defined as one of the 23 canonical α–amino acids incorporated into proteins.
This implementation:
  - Parses the SMILES and removes explicit hydrogens and isotope labels.
  - Enforces a small size (≤25 heavy atoms).
  - Searches for one of several free α–amino acid motifs:
      * Canonical pattern: a chiral α–carbon with unsubstituted amino and carboxyl groups ([C@H](N)(C(=O)[O]) or its mirror)
      * Glycine pattern: an achiral motif (NCC(=O)[O])
      * Proline patterns: the cyclic motif where the α–carbon is part of a five–membered ring 
        (described here as [C@H]1CCCNC1C(=O)[O] or its mirror)
  - For chiral (non–glycine) cases, the α–carbon must have L configuration (CIP code “S”).
If any check fails the function returns False with an explanation.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a proteinogenic amino acid, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # Remove explicit hydrogens for consistency
    mol = Chem.RemoveHs(mol)
    # Remove isotopic labels (e.g. [2H] becomes H) to improve matching for deuterated compounds
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)
        
    # Check allowed elements (C, N, O, S, Se) and small size (free amino acids are small, ≤25 heavy atoms)
    allowed_atomic_nums = {6, 7, 8, 16, 34}
    heavy_atom_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() > 1:
            heavy_atom_count += 1
        if atom.GetAtomicNum() != 1 and atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Disallowed element found: atomic number {atom.GetAtomicNum()}"
    if heavy_atom_count > 25:
        return False, "Molecule is too large to be a single proteinogenic amino acid"
    
    # Ensure stereochemistry is assigned so CIP codes are computed
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Define a list of SMARTS patterns for the free α–amino acid motif.
    # We include patterns for:
    #  (1) canonical amino acids: the α–carbon (index 0) is chiral and then attached to an unsubstituted amino (index 1)
    #      and free carboxyl (index 2) in the motif "[C@H](N)(C(=O)[O])" (and its mirror).
    #  (2) glycine: pattern "NCC(=O)[O]" where the α–carbon is at index 1.
    #  (3) proline: we define patterns so that the cyclic α–carbon is first.
    patterns = [
        {"smarts": "[C@H](N)(C(=O)[O])", "type": "canonical", "alpha_idx": 0, "n_idx": 1, "c_idx": 2},
        {"smarts": "[C@@H](N)(C(=O)[O])", "type": "canonical", "alpha_idx": 0, "n_idx": 1, "c_idx": 2},
        {"smarts": "NCC(=O)[O]", "type": "glycine", "alpha_idx": 1, "n_idx": 0, "c_idx": 2},
        # Proline patterns – here we force the α–carbon to be the first atom in the SMARTS.
        # One way to depict proline is to write the cyclic system with the carboxyl group attached exocyclically.
        {"smarts": "[C@H]1CCCNC1C(=O)[O]", "type": "proline", "alpha_idx": 0, "c_idx": 5},
        {"smarts": "[C@@H]1CCCNC1C(=O)[O]", "type": "proline", "alpha_idx": 0, "c_idx": 5},
    ]
    
    # Iterate over each pattern to try and find a matching free α–amino acid motif.
    for entry in patterns:
        patt = Chem.MolFromSmarts(entry["smarts"])
        if patt is None:
            continue
        # useChirality flag enforces stereochemistry when matching.
        matches = mol.GetSubstructMatches(patt, useChirality=True)
        if not matches:
            continue
        
        # Evaluate each match.
        for match in matches:
            # For canonical and glycine motifs, we further check that:
            # (a) the amino nitrogen is unsubstituted (its only heavy neighbor should be the α–carbon).
            # (b) the carboxyl carbon is properly bonded to two oxygens (one double and one single bond).
            # (c) for non–glycine (canonical and proline) the α–carbon must have CIP (L configuration, "S").
            typ = entry["type"]
            
            # Identify atoms from the match according to our pattern mapping:
            if typ in ("canonical", "glycine"):
                alpha_atom = mol.GetAtomWithIdx(match[ entry["alpha_idx"] ])
                amino_atom = mol.GetAtomWithIdx(match[ entry["n_idx"] ])
                carbo_atom = mol.GetAtomWithIdx(match[ entry["c_idx"] ])
            elif typ == "proline":
                alpha_atom = mol.GetAtomWithIdx(match[ entry["alpha_idx"] ])
                # In proline the amino nitrogen is in the ring; we don’t know the index by our pattern.
                # Instead, find a nitrogen neighbor of the α–carbon that is in a ring.
                amino_atom = None
                for nbr in alpha_atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 7 and nbr.IsInRing():
                        amino_atom = nbr
                        break
                # The carboxyl carbon is taken from the match; for our proline patterns we recorded its index.
                carbo_atom = mol.GetAtomWithIdx(match[ entry["c_idx"] ])
            else:
                continue
            
            # (a) Check that the amino nitrogen is unsubstituted.
            # For a free amino group, the only heavy atom neighbor should be the α–carbon.
            heavy_nbrs = [n.GetIdx() for n in amino_atom.GetNeighbors() if n.GetAtomicNum() > 1]
            if heavy_nbrs != [alpha_atom.GetIdx()] and heavy_nbrs != [alpha_atom.GetIdx(),]:
                continue  # the nitrogen appears substituted
            
            # (b) Check the carboxyl carbon connectivity.
            # It should have exactly 3 heavy neighbors: the α–carbon and two oxygens.
            c_neighbors = [n for n in carbo_atom.GetNeighbors() if n.GetAtomicNum() > 1]
            if len(c_neighbors) != 3:
                continue
            oxygens = [n for n in c_neighbors if n.GetAtomicNum() == 8]
            if len(oxygens) != 2:
                continue
            # Verify one bond is double and the other is single.
            dbl_found = False
            sgl_found = False
            for o in oxygens:
                bond = mol.GetBondBetweenAtoms(carbo_atom.GetIdx(), o.GetIdx())
                if bond is None:
                    continue
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    dbl_found = True
                elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    sgl_found = True
            if not (dbl_found and sgl_found):
                continue
                
            # (c) For non–glycine amino acids, verify the α–carbon is chiral and has CIP code "S" (L configuration)
            if typ in ("canonical", "proline"):
                if not alpha_atom.HasProp('_CIPCode'):
                    continue
                if alpha_atom.GetProp('_CIPCode') != "S":
                    continue  # wrong configuration – likely D or undetermined
            
            # If we reach here, we have found a valid free α–amino acid motif.
            return True, f"Matches free α–amino acid criteria ({typ}): {entry['smarts']}"
        
    return False, "No appropriate free α–amino acid motif (with unsubstituted amino and carboxyl groups and L configuration) found"

# Example ad-hoc testing when running as main:
if __name__ == "__main__":
    test_smiles = {
        "L-glutamic acid": "N[C@@H](CCC(O)=O)C(O)=O",
        "L-histidine": "N[C@@H](Cc1c[nH]cn1)C(O)=O",
        "L-phenylalanine-d5": "OC(=O)[C@@H](N)CC1=C(C(=C(C(=C1[2H])[2H])[2H])[2H])[2H]",
        "L-arginine-d7": "OC(=O)[C@](N([2H])[2H])(C(C(CN=C(N)N)([2H])[2H])([2H])[2H])[2H]",
        "L-proline": "OC(=O)[C@@H]1CCCN1",
        "L-arginine": "N[C@@H](CCCNC(N)=N)C(O)=O",
        "L-leucine": "CC(C)C[C@H](N)C(O)=O",
        "L-isoleucine": "OC([C@H]([C@H](CC)C)N)=O",
        "L-valine": "CC(C)[C@H](N)C(O)=O",
        "L-alanine": "C[C@H](N)C(O)=O",
        "L-valine-d8": "OC(=O)[C@@](N)(C(C([2H])([2H])[2H])(C([2H])([2H])[2H])[2H])[2H]",
        "L-asparagine": "N[C@@H](CC(N)=O)C(O)=O",
        "L-methionine-d3": "S(CC[C@H](N)C(O)=O)C([2H])([2H])[2H]",
        "L-glutamine": "N[C@@H](CCC(N)=O)C(O)=O",
        "L-tyrosine-d4": "OC(=O)[C@@H](N)CC1=C(C(=C(O)C(=C1[2H])[2H])[2H])[2H]",
        "L-aspartic acid-d7": "O(C(=O)[C@@](N([2H])[2H])(C(C(O[2H])=O)([2H])[2H])[2H])[2H]",
        "L-lysine": "NCCCC[C@H](N)C(O)=O",
        "L-serine": "N[C@@H](CO)C(O)=O",
        "L-pyrrolysine": "C(=O)([C@@H](N)CCCCNC([C@H]1[C@@H](CC=N1)C)=O)O",
        "glycine-13C2,15N": "O[13C](=O)[13CH2][15NH2]",
        "L-threonine": "C[C@@H](O)[C@H](N)C(O)=O",
        "L-cysteine": "N[C@@H](CS)C(O)=O",
        "glycine-d5": "C(C(N([2H])[2H])([2H])[2H])(=O)O[2H]",
        "L-aspartic acid": "N[C@@H](CC(O)=O)C(O)=O",
        "L-leucine-d3": "OC(=O)[C@@H](N)CC(C([2H])([2H])[2H])C",
        "L-methionine": "CSCC[C@H](N)C(O)=O",
        # Some examples that previously were false positives:
        "gamma-Glu-Trp(1-)": "[NH3+][C@H](C([O-])=O)CCC(=O)N[C@H](C(=O)[O-])CC1=CNC2=C1C=CC=C2",
        "4-dimethylamino-L-phenylalanine zwitterion": "C1=CC(=CC=C1C[C@@H](C([O-])=O)[NH3+])N(C)C",
    }
    
    for name, smi in test_smiles.items():
        result, reason = is_proteinogenic_amino_acid(smi)
        print(f"{name}: {result} ({reason})")