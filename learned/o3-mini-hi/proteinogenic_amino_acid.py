"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
"""
Classifies: proteinogenic amino acid
Defined as one of the 23 canonical α‐amino acids (including glycine, L-proline, selenocysteine,
pyrrolysine, and N-formylmethionine) that are incorporated into proteins by translation.
In our heuristic we require that the SMILES has a free α–amino acid motif (with terminal amino and 
carboxyl groups), that the number of heavy atoms is small (to rule out peptides), and that only allowed 
elements are present.
"""

from rdkit import Chem

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.
    We check that:
      - The molecule parses correctly.
      - It only contains allowed elements (C, N, O, S, and optionally Se for selenocysteine).
      - Its heavy atom count is not too high (free amino acids are small; we use a threshold of 25).
      - It contains one of a set of canonical SMARTS patterns for a free α–amino acid:
            • conventional (chiral) α–amino acid: "[C@H](N)C(=O)[O]" or "[C@@H](N)C(=O)[O]"
            • glycine (achiral): "NCC(=O)[O]"
            • proline (cyclic): "[C@H]1CCCN1C(=O)[O]" or "[C@@H]1CCCN1C(=O)[O]"
      - In addition, for non-proline cases it checks that the amino nitrogen is only attached to the α–carbon,
        and that the α–carbon is attached exactly to three heavy groups (amino, carboxyl, and side chain).
      - It also verifies that the carboxyl carbon is attached to only three heavy atoms (the α–carbon and two oxygens,
        one of which is double-bonded).
    This helps to distinguish free amino acids from peptides or other carboxylated compounds.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a proteinogenic amino acid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Optionally add explicit hydrogens so that connectivity counts are more complete.
    mol = Chem.AddHs(mol)
    
    # Allowed atomic numbers: C (6), N (7), O (8), S (16) and Se (34) (for selenocysteine)
    allowed_atomic_nums = {6, 7, 8, 16, 34}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:  # hydrogen
            continue
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains disallowed element: atomic number {atom.GetAtomicNum()}"
    
    # Check overall size: free amino acids are small. (Most canonical ones have <=25 heavy atoms)
    heavy_atom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    if heavy_atom_count > 25:
        return False, "Molecule is too large to be a single proteinogenic amino acid"
    
    # Define candidate SMARTS patterns.
    # For typical amino acids, we expect a free carboxyl group attached to a chiral (or in glycine, achiral)
    # α–carbon which in turn is attached to a terminal (unsubstituted) amino group.
    # We allow separate patterns for glycine (which lacks a chiral tag) and for proline (cyclic).
    patterns = [
        {"smarts": "[C@H](N)C(=O)[O]", "type": "canonical"},
        {"smarts": "[C@@H](N)C(=O)[O]", "type": "canonical"},
        {"smarts": "NCC(=O)[O]", "type": "glycine"},
        {"smarts": "[C@H]1CCCN1C(=O)[O]", "type": "proline"},
        {"smarts": "[C@@H]1CCCN1C(=O)[O]", "type": "proline"},
    ]
    
    # Try each SMARTS pattern.
    for entry in patterns:
        patt = Chem.MolFromSmarts(entry["smarts"])
        if patt is None:
            continue
        matches = mol.GetSubstructMatches(patt)
        if not matches:
            continue  # no match of this pattern
        # For each match found, check additional connectivity constraints.
        for match in matches:
            # For the canonical and glycine patterns we assume:
            #   match[0]: α–carbon
            #   match[1]: amino nitrogen
            #   match[2]: carboxyl carbon
            # (For the glycine pattern "NCC(=O)[O]", the α–carbon is the middle C.)
            # For the proline patterns we accept the match as is without additional neighbor checks.
            if entry["type"] in ("canonical", "glycine"):
                if len(match) < 3:
                    continue  # something wrong with the match indexing
                
                alpha_idx = match[0]
                amino_idx = match[1]
                carboxyl_idx = match[2]
                
                alpha_atom = mol.GetAtomWithIdx(alpha_idx)
                amino_atom = mol.GetAtomWithIdx(amino_idx)
                carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
                
                # (a) Check that the amino nitrogen is terminal.
                # It should be connected to exactly one heavy atom – the α–carbon.
                amino_neighbors = [nbr.GetIdx() for nbr in amino_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
                if set(amino_neighbors) != {alpha_idx}:
                    continue  # amino nitrogen appears substituted or part of a peptide bond
                
                # (b) Check the carboxyl carbon: should be connected to exactly 3 heavy atoms
                # (the α–carbon and two oxygens). Also, one oxygen must be double-bonded.
                carboxyl_neighbors = [nbr for nbr in carboxyl_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
                if len(carboxyl_neighbors) != 3:
                    continue
                # Identify oxygens
                oxygens = [nbr for nbr in carboxyl_neighbors if nbr.GetAtomicNum() == 8]
                if len(oxygens) != 2:
                    continue
                # Check bond types – one double and one single between carboxyl carbon and oxygens.
                dbl_bond_found = False
                sgl_bond_found = False
                for o in oxygens:
                    bond = mol.GetBondBetweenAtoms(carboxyl_atom.GetIdx(), o.GetIdx())
                    if bond is None:
                        continue
                    # Using GetBondType() gives a BondType value; compare to Chem.rdchem.BondType.DOUBLE
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        dbl_bond_found = True
                    elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                        sgl_bond_found = True
                if not (dbl_bond_found and sgl_bond_found):
                    continue
                
                # (c) Check the α–carbon has exactly three heavy-atom neighbors.
                alpha_neighbors = [nbr.GetIdx() for nbr in alpha_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
                # We expect the α–carbon to be connected to:
                #   - the amino nitrogen,
                #   - the carboxyl carbon,
                #   - one side chain group.
                if len(alpha_neighbors) != 3:
                    continue

                # Passed all the extra checks for a free amino acid motif.
                return True, f"Matches canonical amino acid criteria: {entry['smarts']}"
            
            elif entry["type"] == "proline":
                # For proline we accept the SMARTS match without extra neighbor checks.
                # (Proline’s amino group is part of the cyclic structure so its nitrogen connectivity differs.)
                return True, f"Matches canonical amino acid criteria (proline): {entry['smarts']}"
    
    # If no pattern with extra checks was found, then the molecule does not appear to be a free (canonical) amino acid.
    return False, "No appropriate free α-amino acid motif (with terminal amino and carboxyl groups) found"

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