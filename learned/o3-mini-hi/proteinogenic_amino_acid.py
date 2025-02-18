"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
#!/usr/bin/env python3
"""
Classifies: proteinogenic amino acid
Defined as one of the 23 canonical α–amino acids (including glycine, L–proline, selenocysteine,
pyrrolysine, and N–formylmethionine) incorporated into proteins by translation.
This implementation verifies that:
  - The molecule parses correctly.
  - It contains only allowed heavy atoms (C, N, O, S, and Se) and is small (≤25 heavy atoms).
  - It contains one of the free α–amino acid motifs (using SMARTS patterns for canonical, glycine, and proline).
  - For canonical motifs (and for proline) the free –NH2 and –COOH groups are unsubstituted and
    the α–carbon connectivity is as expected.
  - For chiral (non–glycine) amino acids, the α–carbon must have the L configuration (CIP code “S”).
If any check fails the function returns False with a reason.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a proteinogenic amino acid, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse SMILES and remove explicit hydrogens (for consistency)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.RemoveHs(mol)
    
    # Check allowed heavy atoms (C, N, O, S, Se)
    allowed_atomic_nums = {6, 7, 8, 16, 34}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Disallowed element found: atomic number {atom.GetAtomicNum()}"
    
    # Check molecule size: free amino acids are small (≤25 heavy atoms)
    heavy_atom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    if heavy_atom_count > 25:
        return False, "Molecule is too large to be a single proteinogenic amino acid"
    
    # Assign stereochemistry so that CIP codes (R/S) are computed.
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Define SMARTS patterns along with type tags.
    # For canonical amino acids (chiral, non glycine):
    #    Expected order: Atom0 = α–carbon, Atom1 = free amino nitrogen, Atom2 = carboxyl carbon.
    # For glycine, pattern "NCC(=O)[O]":
    #    Expected order: Atom0 = amino nitrogen, Atom1 = α–carbon, Atom2 = carboxyl carbon.
    # For proline, the pattern covers the cyclic structure; here we assume the chiral (α–carbon)
    #    is at index 2.
    patterns = [
        {"smarts": "[C@H](N)C(=O)[O]", "type": "canonical"},
        {"smarts": "[C@@H](N)C(=O)[O]", "type": "canonical"},
        {"smarts": "NCC(=O)[O]", "type": "glycine"},
        {"smarts": "OC(=O)[C@H]1CCCN1", "type": "proline"},
        {"smarts": "OC(=O)[C@@H]1CCCN1", "type": "proline"},
    ]
    
    # Iterate through each pattern and try to match.
    for entry in patterns:
        patt = Chem.MolFromSmarts(entry["smarts"])
        if patt is None:
            continue
        matches = mol.GetSubstructMatches(patt, useChirality=True)
        if not matches:
            continue
        
        # Evaluate each match.
        for match in matches:
            if entry["type"] in ("canonical", "glycine"):
                # Ensure the match contains at least three atoms.
                if len(match) < 3:
                    continue
                
                if entry["type"] == "canonical":
                    # For canonical patterns, the order is: [0]: α–carbon, [1]: amino nitrogen, [2]: carboxyl carbon.
                    alpha_idx = match[0]
                    amino_idx = match[1]
                    carboxyl_idx = match[2]
                    expected_alpha_degree = 3  # free amino acid: α–carbon has three heavy neighbors.
                else:  # glycine
                    # For glycine pattern "NCC(=O)[O]", order: [0]: amino nitrogen, [1]: α–carbon, [2]: carboxyl carbon.
                    amino_idx = match[0]
                    alpha_idx = match[1]
                    carboxyl_idx = match[2]
                    expected_alpha_degree = 2  # glycine: no side chain.
                
                # (a) Check amino nitrogen: it should be free (only attached to α–carbon).
                amino_atom = mol.GetAtomWithIdx(amino_idx)
                amino_neighbors = [nbr.GetIdx() for nbr in amino_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
                if set(amino_neighbors) != {alpha_idx}:
                    continue  # amino nitrogen appears substituted.
                
                # (b) Check carboxyl carbon connectivity: should be bonded to α–carbon and two oxygens,
                # with one O in a double bond and one in a single bond.
                carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
                nbrs = [nbr for nbr in carboxyl_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
                if len(nbrs) != 3:
                    continue
                oxygens = [nbr for nbr in nbrs if nbr.GetAtomicNum() == 8]
                if len(oxygens) != 2:
                    continue
                dbl_found = False
                sgl_found = False
                for o in oxygens:
                    bond = mol.GetBondBetweenAtoms(carboxyl_atom.GetIdx(), o.GetIdx())
                    if bond is None:
                        continue
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        dbl_found = True
                    elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                        sgl_found = True
                if not (dbl_found and sgl_found):
                    continue
                
                # (c) Check α–carbon connectivity: it should have the expected number of heavy atom neighbors.
                alpha_atom = mol.GetAtomWithIdx(alpha_idx)
                alpha_neighbors = [nbr.GetIdx() for nbr in alpha_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
                if len(alpha_neighbors) != expected_alpha_degree:
                    continue
                
                # (d) For canonical amino acids (non–glycine), check that the α–carbon has L configuration,
                # i.e. CIP code should be "S". (Glycine is achiral.)
                if entry["type"] == "canonical":
                    if not alpha_atom.HasProp('_CIPCode'):
                        continue
                    if alpha_atom.GetProp('_CIPCode') != "S":
                        continue  # wrong configuration (likely D amino acid)
                        
                return True, f"Matches canonical amino acid criteria: {entry['smarts']}"
                
            elif entry["type"] == "proline":
                # For proline, use the match and assume the cyclic backbone.
                # We assume that the chiral (α–carbon) appears at index 2 in the substructure match.
                if len(match) < 3:
                    continue
                alpha_idx = match[2]
                alpha_atom = mol.GetAtomWithIdx(alpha_idx)
                
                # For proline, we also check that the free carboxyl (or its ester equivalent) and amino groups are not substituted.
                # We will simply check that the amino nitrogen (which is in the ring) is not extra substituted.
                # Identify the amino nitrogen in the ring.
                ring_nitrogens = [nbr for nbr in alpha_atom.GetNeighbors() if nbr.GetAtomicNum() == 7]
                valid_ring_nitrogen = False
                for n in ring_nitrogens:
                    n_neighbors = [nbr.GetIdx() for nbr in n.GetNeighbors() if nbr.GetAtomicNum() > 1]
                    # If the nitrogen is connected only to atoms within the ring (typically 2 heavy neighbors) then accept.
                    if len(n_neighbors) <= 2:
                        valid_ring_nitrogen = True
                if not valid_ring_nitrogen:
                    continue
                    
                # Also enforce L configuration (α–carbon CIP code "S")
                if not alpha_atom.HasProp('_CIPCode'):
                    continue
                if alpha_atom.GetProp('_CIPCode') != "S":
                    continue
                    
                return True, f"Matches canonical amino acid criteria (proline): {entry['smarts']}"
    
    return False, "No appropriate free α–amino acid motif (with unsubstituted amino and carboxyl groups and L configuration) found"

# For ad-hoc testing when run as a script:
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
        # Examples noted as false positives previously:
        "gamma-Glu-Trp(1-)": "[NH3+][C@H](C([O-])=O)CCC(=O)N[C@H](C(=O)[O-])CC1=CNC2=C1C=CC=C2",
        "4-dimethylamino-L-phenylalanine zwitterion": "C1=CC(=CC=C1C[C@@H](C([O-])=O)[NH3+])N(C)C",
        "D-citrulline": "N[C@H](CCCNC(N)=O)C(O)=O",
    }
    
    for name, smi in test_smiles.items():
        result, reason = is_proteinogenic_amino_acid(smi)
        print(f"{name}: {result} ({reason})")