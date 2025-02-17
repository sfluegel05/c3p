"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
"""
Classifies: proteinogenic amino acid
Defined as one of the 23 canonical α‐amino acids that are incorporated into proteins.
Examples: L-glutamic acid, L-histidine, glycine (non-chiral), L-proline (cyclic) etc.
The program uses an improved heuristic that checks
  (a) for the presence of an α–carbon attached to an amino group and a carboxyl group,
  (b) that the amino group is terminal (i.e. not part of a peptide bond, unless it’s proline),
  (c) that the carboxyl group is free (only attached to the α–carbon and its two oxygens),
  (d) that the overall molecule is relatively small (free amino acids have few heavy atoms),
  (e) that no disallowed elements (e.g. halogens) are present.
"""

from rdkit import Chem

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a canonical (proteinogenic) amino acid based on its SMILES string.
    The heuristic checks for an α–carbon bearing an amino group and a carboxyl group in a free (non‐peptidic)
    setting. In peptides, the amino group or carboxyl group will be further substituted, and the overall
    heavy-atom count will be larger.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a proteinogenic amino acid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule only contains allowed elements (C, N, O, S).
    allowed_atomic_nums = {6, 7, 8, 16}  # Carbon, Nitrogen, Oxygen, Sulfur
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue  # Hydrogens (explicit or implicit)
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains disallowed element: atomic number {atom.GetAtomicNum()}"
    
    # Check overall molecule size: free amino acids are small.
    # Count heavy atoms (i.e. non-hydrogen atoms).
    heavy_atom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    if heavy_atom_count > 25:
        return False, "Molecule is too large to be a single proteinogenic amino acid"
    
    # Define backbone SMARTS patterns with a flag indicating if this pattern is used for proline.
    # We include both chiral and non-chiral (glycine) variants.
    patterns = [
        {"smarts": "[C@H](N)C(=O)[O]", "proline": False},
        {"smarts": "[C@@H](N)C(=O)[O]", "proline": False},
        {"smarts": "NC(C(=O)[O])", "proline": False},  # Glycine does not have a chiral tag.
        {"smarts": "[C@H]1CCCN1C(=O)[O]", "proline": True},
        {"smarts": "[C@@H]1CCCN1C(=O)[O]", "proline": True}
    ]
    
    # Try each pattern.
    for entry in patterns:
        patt = Chem.MolFromSmarts(entry["smarts"])
        if patt is None:
            continue
        matches = mol.GetSubstructMatches(patt)
        if not matches:
            continue
        
        # For each match found, perform additional checks.
        # We assume the SMARTS patterns are written so that:
        #   index 0: the α-carbon; index 1: the amino nitrogen; index 2: the carboxyl carbon.
        for match in matches:
            alpha_idx = match[0]
            amino_idx = match[1]
            carboxyl_idx = match[2]

            alpha_atom = mol.GetAtomWithIdx(alpha_idx)
            amino_atom = mol.GetAtomWithIdx(amino_idx)
            carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
            
            # For non-proline cases, the amino nitrogen should only be bonded to the α–carbon.
            # In a free amino acid the nitrogen is not further substituted (other than implicit hydrogens).
            if not entry["proline"]:
                # Get heavy-atom neighbors of the amino nitrogen (exclude implicit hydrogens).
                amino_neighbors = [nbr.GetIdx() for nbr in amino_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
                # It should be bonded only to the α–carbon.
                if set(amino_neighbors) != {alpha_idx}:
                    # This nitrogen is substituted (e.g. peptide bond) so not a free amino group.
                    continue
            
            # Check the carboxyl carbon.
            # It should be connected to the α-carbon plus exactly two oxygens in a free carboxyl group.
            carboxyl_neighbors = [nbr for nbr in carboxyl_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
            # Expecting one neighbor to be the α-carbon and two oxygens.
            if len(carboxyl_neighbors) != 3:
                continue
            # Identify the oxygen atoms:
            oxygen_neighbors = [nbr for nbr in carboxyl_neighbors if nbr.GetAtomicNum() == 8]
            if len(oxygen_neighbors) != 2:
                continue
            # Check that one of the oxygens is double-bonded (C=O) and the other is single-bonded.
            dbl_bond_found = False
            sgl_bond_found = False
            for nbr in oxygen_neighbors:
                bond = mol.GetBondBetweenAtoms(carboxyl_atom.GetIdx(), nbr.GetIdx())
                if bond is not None:
                    if bond.GetBondTypeAsDouble() == 1.0:  # double bond check
                        dbl_bond_found = True
                    elif bond.GetBondTypeAsDouble() != 1.0:
                        sgl_bond_found = True
            if not (dbl_bond_found and sgl_bond_found):
                continue
            
            # Also check that the α-carbon is connected to exactly three groups:
            # the amino group, the carboxyl group, and a side chain.
            alpha_neighbors = [nbr.GetIdx() for nbr in alpha_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
            # Remove known neighbors: amino nitrogen and carboxyl carbon.
            side_chain_neighbors = set(alpha_neighbors) - {amino_idx, carboxyl_idx}
            if len(side_chain_neighbors) != 1:
                # For example, if more than one substituent is attached, it could be a di- or tri-peptide fragment.
                continue
            
            # Passed all the extra checks. Return True.
            return True, f"Matches canonical amino acid criteria: {entry['smarts']}"
    
    # If no pattern passed the additional validations, we conclude that the molecule is not a proteinogenic amino acid.
    return False, "No appropriate free α-amino acid motif (with terminal amino and carboxyl groups) found"

# For ad-hoc testing when run as a script:
if __name__ == "__main__":
    # Test a few examples
    test_smiles = {
        "L-glutamic acid": "N[C@@H](CCC(O)=O)C(O)=O",
        "L-histidine": "N[C@@H](Cc1c[nH]cn1)C(O)=O",
        "Aminofluoropropionic acid (false positive expected)": "FC(N)(C)C(O)=O",
        "Thr-Leu-Trp (peptide, false positive expected)": "O=C(N[C@@H](CC=1C=2C(NC1)=CC=CC2)C(O)=O)[C@@H](NC(=O)[C@@H](N)[C@H](O)C)CC(C)C",
        "Invalid": "NotAValidSMILES"
    }
    for name, smi in test_smiles.items():
        result, reason = is_proteinogenic_amino_acid(smi)
        print(f"{name}: {result} ({reason})")