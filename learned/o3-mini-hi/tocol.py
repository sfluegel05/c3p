"""
Classifies: CHEBI:39437 tocol
"""
"""
Classifies: tocol
Definition: A chromanol with a chroman-6-ol skeleton that is substituted at position 2 
            by a saturated or triply-unsaturated hydrocarbon chain consisting of 
            three isoprenoid units.
            
This program uses rdkit to:
  1. Check for the chromanol (benzopyran) core with an aromatic –OH.
  2. Locate the substituent at position 2 using an atom-mapped SMARTS query.
  3. Extract the side chain (avoiding core atoms) and verify:
       • It is aliphatic, composed solely of carbon atoms.
       • It is long enough (heuristically, at least 10 carbons) to accommodate three isoprenoid units.
       • It is either saturated (0 C=C bonds) or has exactly three C=C double bonds.
       
If any of these checks fail, the function returns False with an explanation.
Note: Tocol chemistry is complex; this is one heuristic implementation.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdchem import BondType

def is_tocol(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) matches the tocol class
    (a chromanol with a specific 2‑position side chain as defined).

    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if the molecule meets the tocol definition.
       str: A reason explaining the decision.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Step 1. Look for the chromanol (benzopyran) core.
    # We use a SMARTS with an atom mapping to indicate the substitution at position 2.
    # The query "O1CC([*:1])c2ccccc2C1" attempts to capture a 2H-1-benzopyran system,
    # where the atom labelled [*:1] represents the expected attachment point at position 2.
    core_smarts = "O1CC([*:1])c2ccccc2C1"
    core_query = Chem.MolFromSmarts(core_smarts)
    # Correctly use GetSubstructMatch
    core_match = mol.GetSubstructMatch(core_query)
    if not core_match:
        return False, "Chromanol core (benzopyran) not found"

    # Record the atoms that are part of the core pattern.
    core_atom_indices = set(core_match)

    # Step 2. Check that the benzene (aromatic) part of the core bears a hydroxyl group.
    # Look for an aromatic atom in a phenolic OH.
    phenol_query = Chem.MolFromSmarts("c[OH]")
    if not mol.HasSubstructMatch(phenol_query):
        return False, "Aromatic hydroxyl group (phenolic OH) not found in the chromanol core"

    # Step 3. Locate the side chain attached to the 2-position.
    # The [*:1] atom on the query should appear in the molecule as an atom with a molAtomMapNumber property set to "1".
    substitution_root = None
    for atom in mol.GetAtoms():
        if atom.HasProp("molAtomMapNumber"):
            try:
                if int(atom.GetProp("molAtomMapNumber")) == 1:
                    substitution_root = atom
                    break
            except ValueError:
                continue
    if substitution_root is None:
        return False, "Could not locate the 2-position substituent from the core SMARTS match"

    # The substitution root may itself be embedded in the core.
    # So find its neighbor that is not part of the core to get the beginning of the side chain.
    chain_root = None
    for neighbor in substitution_root.GetNeighbors():
        if neighbor.GetIdx() not in core_atom_indices:
            chain_root = neighbor
            break
    if chain_root is None:
        return False, "No side chain found attached at the 2-position"

    # Step 4. Extract the entire side chain by traversing from the chain_root,
    # avoiding any atoms that are part of the chromanol core.
    side_chain_idxs = set()
    to_visit = [chain_root.GetIdx()]
    while to_visit:
        current = to_visit.pop()
        if current in side_chain_idxs:
            continue
        side_chain_idxs.add(current)
        current_atom = mol.GetAtomWithIdx(current)
        for nbr in current_atom.GetNeighbors():
            if nbr.GetIdx() in core_atom_indices:
                continue
            if nbr.GetIdx() not in side_chain_idxs:
                to_visit.append(nbr.GetIdx())

    # Step 5. Check that the side chain is aliphatic (only carbon atoms).
    for idx in side_chain_idxs:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6:
            return False, f"Side chain contains non-carbon atom: {atom.GetSymbol()} (atom idx: {idx})"
    n_chain_carbons = len(side_chain_idxs)
    if n_chain_carbons < 10:
        return False, f"Side chain is too short (found {n_chain_carbons} carbons; expected at least 10)"

    # Step 6. Count the number of double bonds in the side chain.
    double_bond_count = 0
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        # Count only bonds occurring entirely within the side chain.
        if a1 in side_chain_idxs and a2 in side_chain_idxs:
            if bond.GetBondType() == BondType.DOUBLE:
                double_bond_count += 1

    # Accept either a saturated side chain (0 double bonds) or one with exactly three double bonds.
    if double_bond_count not in (0, 3):
        return False, f"Side chain has {double_bond_count} double bond(s); expected 0 (saturated) or 3 (triply unsaturated)"

    return True, "Molecule contains chromanol core with appropriate aromatic -OH and a 2-substituent chain matching tocol definition"

# For testing (remove or modify these lines as needed):
if __name__ == "__main__":
    # Test with (S,R,S)-alpha-tocopherol SMILES as provided in the prompt.
    test_smiles = "CC(C)CCC[C@H](C)CCC[C@@H](C)CCC[C@@]1(C)CCC2=C(C)C(O)=C(C)C(C)=C2O1"
    flag, reason = is_tocol(test_smiles)
    print(flag, reason)