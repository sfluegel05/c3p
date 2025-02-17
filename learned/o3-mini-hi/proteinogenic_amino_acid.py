"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
#!/usr/bin/env python3
"""
Classifies: proteinogenic amino acid
Defined as one of the 23 canonical α‐amino acids (including glycine, L-proline, selenocysteine,
pyrrolysine, and N-formylmethionine) that are incorporated into proteins by translation.
This implementation checks that the molecule:
  - Parses correctly.
  - Contains only allowed elements (C, N, O, S, and Se).
  - Is small (≤25 heavy atoms).
  - Contains one of the free α–amino acid motifs.
For canonical (or glycine) motifs it further verifies that the free -NH2 and -COOH groups are not substituted.
For proline (cyclic) cases a separate SMARTS is used.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.

    The function first ensures that the molecule only contains allowed heavy atoms (C, N, O, S, Se)
    and that its heavy atom count is ≤25. It then uses a set of SMARTS patterns that capture the free
    α–amino acid motif. For canonical chiral amino acids and for glycine, further connectivity checks are made:
      - For chiral amino acids (canonical), the α–carbon (first atom in the pattern) must have exactly
        three heavy neighbours (amino N, carboxyl C, and the side chain) and the terminal amino nitrogen must
        be bonded only to that α–carbon.
      - For glycine, the α–carbon (which in the glycine pattern is the second atom) must have exactly two heavy neighbours.
      - For both, the carboxyl carbon must have exactly three heavy neighbours (the α–carbon and two O atoms),
        and one O bond should be double while the other single.
    For proline patterns the substructure is accepted as is.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a proteinogenic amino acid, False otherwise.
        str: Description of the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove explicit hydrogens for consistent connectivity (this also minimizes isotope label issues).
    mol = Chem.RemoveHs(mol)
    
    # Check for allowed elements (heavy atoms only: C, N, O, S, and optionally Se).
    allowed_atomic_nums = {6, 7, 8, 16, 34}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:  # Skip hydrogen.
            continue
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains disallowed element: atomic number {atom.GetAtomicNum()}"
    
    # Check overall size: free amino acids are small (≤25 heavy atoms).
    heavy_atom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    if heavy_atom_count > 25:
        return False, "Molecule is too large to be a single proteinogenic amino acid"
    
    # Define SMARTS patterns along with a type tag.
    # For canonical amino acids the expected ordering is:
    #    canonical patterns: [C@H](N)C(=O)[O] or [C@@H](N)C(=O)[O]
    #       match order: alpha-carbon, amino nitrogen, carboxyl carbon.
    #    For glycine, which is achiral: the pattern "NCC(=O)[O]" yields
    #       match order: amino nitrogen, alpha-carbon, carboxyl carbon.
    #    For proline we use patterns that account for the cyclic structure.
    patterns = [
        {"smarts": "[C@H](N)C(=O)[O]", "type": "canonical"},
        {"smarts": "[C@@H](N)C(=O)[O]", "type": "canonical"},
        {"smarts": "NCC(=O)[O]", "type": "glycine"},
        {"smarts": "OC(=O)[C@H]1CCCN1", "type": "proline"},
        {"smarts": "OC(=O)[C@@H]1CCCN1", "type": "proline"},
    ]
    
    # Try to match one of the defined patterns.
    for entry in patterns:
        patt = Chem.MolFromSmarts(entry["smarts"])
        if patt is None:
            continue  # Skip bad SMARTS.
        matches = mol.GetSubstructMatches(patt, useChirality=True)
        if not matches:
            continue  # No match for this pattern.
        # For each match, perform additional connectivity validations for canonical and glycine types.
        for match in matches:
            if entry["type"] in ("canonical", "glycine"):
                # Require at least three atoms in the match.
                if len(match) < 3:
                    continue
                # Use proper ordering of atoms:
                if entry["type"] == "canonical":
                    # Expected order: [0]: α–carbon, [1]: amino nitrogen, [2]: carboxyl carbon.
                    alpha_idx = match[0]
                    amino_idx = match[1]
                    carboxyl_idx = match[2]
                    expected_alpha_degree = 3  # chiral α–carbon has an extra side-chain group.
                else:  # glycine
                    # For glycine "NCC(=O)[O]", order: [0]: amino nitrogen, [1]: α–carbon, [2]: carboxyl carbon.
                    amino_idx = match[0]
                    alpha_idx = match[1]
                    carboxyl_idx = match[2]
                    expected_alpha_degree = 2  # glycine has no side-chain.
                
                # (a) Check that the amino nitrogen is “free” (bonded only to the α–carbon).
                amino_atom = mol.GetAtomWithIdx(amino_idx)
                amino_heavy_neighbors = [nbr.GetIdx() for nbr in amino_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
                if set(amino_heavy_neighbors) != {alpha_idx}:
                    continue  # The amino nitrogen appears substituted.
                
                # (b) Check the carboxyl carbon – it must bond only to the α–carbon and two oxygens.
                carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
                nbh = [nbr for nbr in carboxyl_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
                if len(nbh) != 3:
                    continue
                # Identify oxygen neighbors.
                oxygens = [nbr for nbr in nbh if nbr.GetAtomicNum() == 8]
                if len(oxygens) != 2:
                    continue
                # Check that one oxygen is double-bonded and the other is single-bonded.
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
                
                # (c) Check the α–carbon connectivity.
                alpha_atom = mol.GetAtomWithIdx(alpha_idx)
                alpha_heavy_neighbors = [nbr.GetIdx() for nbr in alpha_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
                if len(alpha_heavy_neighbors) != expected_alpha_degree:
                    continue
                
                # Passed all checks for a free α–amino acid motif.
                return True, f"Matches canonical amino acid criteria: {entry['smarts']}"
            
            elif entry["type"] == "proline":
                # For proline, we accept the SMARTS match without extra connectivity checks.
                return True, f"Matches canonical amino acid criteria (proline): {entry['smarts']}"
    
    # If no pattern passed all checks, report failure.
    return False, "No appropriate free α–amino acid motif (with terminal amino and carboxyl groups) found"

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
        # Some examples noted as false positives previously:
        "gamma-Glu-Trp(1-)": "[NH3+][C@H](C([O-])=O)CCC(=O)N[C@H](C(=O)[O-])CC1=CNC2=C1C=CC=C2",
        "4-dimethylamino-L-phenylalanine zwitterion": "C1=CC(=CC=C1C[C@@H](C([O-])=O)[NH3+])N(C)C",
        "D-citrulline": "N[C@H](CCCNC(N)=O)C(O)=O",
    }
    
    for name, smi in test_smiles.items():
        result, reason = is_proteinogenic_amino_acid(smi)
        print(f"{name}: {result} ({reason})")